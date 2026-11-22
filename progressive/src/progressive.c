
/* progressive: tiled, resumable software renderer.
**
** The image is split into independent tiles (blocks). Pixel shading is a
** pure function of pixel coordinates: per-pixel hash-seeded sampling and
** no cross-pixel mutable state, so tile order, interruption and resume
** cannot change any pixel. A render stopped and resumed from a checkpoint
** is byte-for-byte identical to a one-shot row-major render.
**
** Cancellation is cooperative and polled at tile boundaries, so every
** tile that is marked done (and every checkpoint saved) is complete:
** a tile is either fully rendered or never started.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define JOB_MAGIC   0x4A54524Du            /* 'MRTJ' */
#define JOB_VERSION 1u
#define MAX_DIM     65536
#define SPP         4                       /* samples per pixel */

typedef struct { unsigned char r, g, b, a; } RGBA8;

typedef struct {
    int width, height;
    int tile_w, tile_h;
    int tiles_x, tiles_y, ntiles;
    int next_tile;                         /* scheduler cursor */
    int done_count;
    int abort_request;                     /* set by caller; polled per tile */
    unsigned char *done;                   /* ntiles bytes, 0/1 */
    RGBA8 *fb;                             /* width * height pixels */
} RenderJob;

typedef struct {
    int tiles_done, tiles_total, tiles_x, tiles_y;
    float fraction;
} Progress;

/* ---------- deterministic per-pixel mini ray tracer ---------- */

typedef struct { float x, y, z, r; float albedo[3]; } Sphere;

static const Sphere SPHERES[] = {
    {  0.00f,  0.00f, 3.0f, 1.00f, { 0.85f, 0.35f, 0.25f } },
    { -1.10f, -0.60f, 2.6f, 0.45f, { 0.30f, 0.75f, 0.45f } },
    {  1.05f, -0.50f, 2.3f, 0.40f, { 0.35f, 0.45f, 0.90f } },
};
#define N_SPHERES (int)(sizeof(SPHERES) / sizeof(SPHERES[0]))

static uint32_t hash32(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352dU;
    x ^= x >> 15; x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}

static float rand01(uint32_t seed) {
    return (float)(hash32(seed) >> 8) * (1.0f / 16777216.0f);
}

static int sphere_hit(const Sphere* s, const float o[3], const float d[3], float* t) {
    float lx = o[0] - s->x, ly = o[1] - s->y, lz = o[2] - s->z;
    float b = lx * d[0] + ly * d[1] + lz * d[2];
    float c = lx * lx + ly * ly + lz * lz - s->r * s->r;
    float disc = b * b - c;
    if (disc < 0.0f) return 0;
    float sq = sqrtf(disc);
    float th = -b - sq;
    if (th <= 1e-4f) th = -b + sq;
    if (th <= 1e-4f) return 0;
    *t = th;
    return 1;
}

/* pure function of screen position + image dimensions */
static void trace_scene(float sx, float sy, int w, int h, float out[3]) {
    float aspect = (float)w / (float)h;
    float d[3] = { sx * aspect * 0.9f, sy * 0.9f, 1.0f };
    float dl = sqrtf(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    d[0] /= dl; d[1] /= dl; d[2] /= dl;
    float o[3] = { 0.0f, 0.0f, 0.0f };

    /* vertical gradient background */
    float t = (d[1] + 1.0f) * 0.5f;
    out[0] = 0.02f + (0.25f - 0.02f) * t;
    out[1] = 0.03f + (0.45f - 0.03f) * t;
    out[2] = 0.06f + (0.75f - 0.06f) * t;

    int hit = -1;
    float best = 1e30f;
    for (int i = 0; i < N_SPHERES; i++) {
        float th;
        if (sphere_hit(&SPHERES[i], o, d, &th) && th < best) { best = th; hit = i; }
    }
    if (hit < 0) return;

    float px = d[0] * best, py = d[1] * best, pz = d[2] * best;
    float nx = (px - SPHERES[hit].x) / SPHERES[hit].r;
    float ny = (py - SPHERES[hit].y) / SPHERES[hit].r;
    float nz = (pz - SPHERES[hit].z) / SPHERES[hit].r;
    float lx = -0.5f, ly = 0.8f, lz = 0.6f;
    float ll = sqrtf(lx * lx + ly * ly + lz * lz);
    lx /= ll; ly /= ll; lz /= ll;
    float diff = nx * lx + ny * ly + nz * lz;
    if (diff < 0.0f) diff = 0.0f;
    float k = 0.15f + diff;
    out[0] = SPHERES[hit].albedo[0] * k;
    out[1] = SPHERES[hit].albedo[1] * k;
    out[2] = SPHERES[hit].albedo[2] * k;
}

/* Pure pixel shader: depends only on (x, y, w, h). Never touches global
** mutable state, so calling it in tile order, row order, after a resume,
** etc. always yields the same byte. */
static RGBA8 shade_pixel(int x, int y, int w, int h) {
    double acc[3] = { 0.0, 0.0, 0.0 };
    for (int s = 0; s < SPP; s++) {
        uint32_t seed = (uint32_t)x * 73856093u
                      ^ (uint32_t)y * 19349663u
                      ^ (uint32_t)(s + 1) * 83492791u;
        float jx = rand01(hash32(seed));
        float jy = rand01(hash32(seed ^ 0x9e3779b9u));
        float sx = 2.0f * ((float)x + jx) / (float)w - 1.0f;
        float sy = 1.0f - 2.0f * ((float)y + jy) / (float)h;
        float c[3];
        trace_scene(sx, sy, w, h, c);
        acc[0] += c[0]; acc[1] += c[1]; acc[2] += c[2];
    }
    RGBA8 px;
    unsigned char* ch = &px.r;
    for (int k = 0; k < 3; k++) {
        float v = (float)(acc[k] / SPP);
        if (v > 1.0f) v = 1.0f;
        ch[k] = (unsigned char)(sqrtf(v) * 255.0f + 0.5f);
    }
    px.a = 255;
    return px;
}

/* ---------- job lifecycle ---------- */

static void job_free(RenderJob* j) {
    if (!j) return;
    free(j->done);
    free(j->fb);
    memset(j, 0, sizeof(*j));
}

/* Rejects zero or negative block sizes (and image dimensions). */
static int job_init(RenderJob* j, int w, int h, int tile_w, int tile_h) {
    if (!j) return -1;
    if (w <= 0 || h <= 0 || tile_w <= 0 || tile_h <= 0) return -1;
    if (w > MAX_DIM || h > MAX_DIM) return -1;
    memset(j, 0, sizeof(*j));
    j->width = w; j->height = h;
    j->tile_w = tile_w; j->tile_h = tile_h;
    j->tiles_x = (w + tile_w - 1) / tile_w;
    j->tiles_y = (h + tile_h - 1) / tile_h;
    j->ntiles = j->tiles_x * j->tiles_y;
    j->done = calloc((size_t)j->ntiles, 1);
    j->fb = malloc((size_t)w * (size_t)h * 4u);
    if (!j->done || !j->fb) { job_free(j); return -1; }
    return 0;
}

static void tile_rect(const RenderJob* j, int t, int* x0, int* y0, int* x1, int* y1) {
    int tx = t % j->tiles_x;
    int ty = t / j->tiles_x;
    *x0 = tx * j->tile_w;
    *y0 = ty * j->tile_h;
    *x1 = *x0 + j->tile_w;
    *y1 = *y0 + j->tile_h;
    if (*x1 > j->width)  *x1 = j->width;
    if (*y1 > j->height) *y1 = j->height;
}

/* Renders one tile. Idempotent: a finished tile is skipped. Returns 1 if
** work was done, 0 if the tile was already complete. */
static int render_tile(RenderJob* j, int t) {
    if (!j || t < 0 || t >= j->ntiles || j->done[t]) return 0;
    int x0, y0, x1, y1;
    tile_rect(j, t, &x0, &y0, &x1, &y1);
    for (int y = y0; y < y1; y++)
        for (int x = x0; x < x1; x++)
            j->fb[y * j->width + x] = shade_pixel(x, y, j->width, j->height);
    j->done[t] = 1;
    j->done_count++;
    return 1;
}

/* Renders up to max_tiles pending tiles, in deterministic tile order.
** Aborts at a tile boundary when abort_request is set, so no partial
** tile is ever observable. Returns the number of tiles rendered. */
static int render_step(RenderJob* j, int max_tiles) {
    if (!j) return -1;
    int n = 0;
    while (max_tiles > 0 && n < max_tiles) {
        while (j->next_tile < j->ntiles && j->done[j->next_tile])
            j->next_tile++;
        if (j->next_tile >= j->ntiles) break;
        if (j->abort_request) break;
        if (!render_tile(j, j->next_tile)) break;
        j->next_tile++;
        n++;
    }
    return n;
}

static int job_is_complete(const RenderJob* j) {
    return j && j->done_count == j->ntiles;
}

static Progress job_progress(const RenderJob* j) {
    Progress p;
    p.tiles_done = j->done_count;
    p.tiles_total = j->ntiles;
    p.tiles_x = j->tiles_x;
    p.tiles_y = j->tiles_y;
    p.fraction = j->ntiles > 0 ? (float)j->done_count / (float)j->ntiles : 0.0f;
    return p;
}

/* ---------- checkpoint save / load ---------- */

/* File layout: 7 x uint32 header, done bitmap (ntiles bytes), RGBA8 fb. */
static int job_save(const RenderJob* j, const char* path) {
    if (!j || !path) return -1;
    FILE* f = fopen(path, "wb");
    if (!f) return -1;
    uint32_t hdr[7] = {
        JOB_MAGIC, JOB_VERSION,
        (uint32_t)j->width, (uint32_t)j->height,
        (uint32_t)j->tile_w, (uint32_t)j->tile_h,
        (uint32_t)j->ntiles
    };
    int ok = fwrite(hdr, sizeof(hdr), 1, f) == 1
          && fwrite(j->done, 1, (size_t)j->ntiles, f) == (size_t)j->ntiles
          && fwrite(j->fb, 1, (size_t)j->width * (size_t)j->height * 4u, f)
             == (size_t)j->width * (size_t)j->height * 4u
          && ferror(f) == 0;
    if (fclose(f) != 0) ok = 0;
    return ok ? 0 : -1;
}

static int job_load(RenderJob* j, const char* path) {
    if (!j || !path) return -1;
    memset(j, 0, sizeof(*j));
    FILE* f = fopen(path, "rb");
    if (!f) return -1;

    uint32_t hdr[7];
    if (fread(hdr, sizeof(hdr), 1, f) != 1) goto fail;
    if (hdr[0] != JOB_MAGIC || hdr[1] != JOB_VERSION) goto fail;

    int w = (int)hdr[2], h = (int)hdr[3], tw = (int)hdr[4], th = (int)hdr[5];
    int nt = (int)hdr[6];
    if (w <= 0 || h <= 0 || tw <= 0 || th <= 0 || w > MAX_DIM || h > MAX_DIM)
        goto fail;
    int tx = (w + tw - 1) / tw, ty = (h + th - 1) / th;
    if (nt != tx * ty || nt <= 0) goto fail;

    j->width = w; j->height = h; j->tile_w = tw; j->tile_h = th;
    j->tiles_x = tx; j->tiles_y = ty; j->ntiles = nt;
    j->done = calloc((size_t)nt, 1);
    j->fb = malloc((size_t)w * (size_t)h * 4u);
    if (!j->done || !j->fb) goto fail;

    if (fread(j->done, 1, (size_t)nt, f) != (size_t)nt) goto fail;
    size_t fb_bytes = (size_t)w * (size_t)h * 4u;
    if (fread(j->fb, 1, fb_bytes, f) != fb_bytes) goto fail;
    if (fgetc(f) != EOF || ferror(f)) goto fail;      /* truncated or trailing junk */
    fclose(f);

    for (int i = 0; i < nt; i++) {
        if (j->done[i] > 1) { job_free(j); return -1; }
        if (j->done[i]) j->done_count++;
    }
    j->next_tile = 0;
    while (j->next_tile < nt && j->done[j->next_tile]) j->next_tile++;
    return 0;

fail:
    fclose(f);
    job_free(j);
    return -1;
}

/* one-shot row-major render, the pre-progressive baseline */
static void render_reference(RenderJob* j) {
    for (int y = 0; y < j->height; y++)
        for (int x = 0; x < j->width; x++)
            j->fb[y * j->width + x] = shade_pixel(x, y, j->width, j->height);
}

/* ---------- self tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define CKPT_PATH  "/tmp/mini_rt_prog_checkpoint.bin"
#define JUNK_PATH  "/tmp/mini_rt_prog_junk.bin"
#define TRUNC_PATH "/tmp/mini_rt_prog_truncated.bin"

int main(void) {
    /* 1. zero / negative block sizes must be rejected */
    RenderJob bad;
    check("zero_tile_w_rejected", job_init(&bad, 16, 16, 0, 4) == -1);
    check("zero_tile_h_rejected", job_init(&bad, 16, 16, 4, 0) == -1);
    check("zero_image_rejected",   job_init(&bad, 0, 0, 4, 4) == -1);
    check("negative_tile_rejected", job_init(&bad, 16, 16, -2, 4) == -1);

    /* 30x20 image, 8x6 tiles -> 4x4 grid, edge tiles clipped, 16 tiles */
    const int W = 30, H = 20, TW = 8, TH = 6;
    RenderJob ref;
    check("init_ok", job_init(&ref, W, H, TW, TH) == 0);
    check("tile_grid", ref.tiles_x == 4 && ref.tiles_y == 4 && ref.ntiles == 16);

    /* 2. tile rectangles partition the image: every pixel covered once */
    unsigned char* mask = calloc((size_t)W * H, 1);
    int cover_ok = mask != NULL;
    for (int t = 0; t < ref.ntiles && cover_ok; t++) {
        int x0, y0, x1, y1;
        tile_rect(&ref, t, &x0, &y0, &x1, &y1);
        for (int y = y0; y < y1; y++)
            for (int x = x0; x < x1; x++) {
                int idx = y * W + x;
                if (mask[idx]) cover_ok = 0;
                mask[idx] = 1;
            }
    }
    for (int i = 0; i < W * H; i++)
        if (!mask[i]) cover_ok = 0;
    check("tiles_partition_image", cover_ok);
    free(mask);

    /* 3. one-shot row-major reference, every pixel equals the pure shader */
    render_reference(&ref);
    int ref_ok = 1;
    for (int y = 0; y < H && ref_ok; y++)
        for (int x = 0; x < W; x++) {
            RGBA8 want = shade_pixel(x, y, W, H);
            if (memcmp(&ref.fb[y * W + x], &want, sizeof(RGBA8)) != 0) {
                ref_ok = 0; break;
            }
        }
    check("reference_pixels", ref_ok);

    size_t fb_bytes = (size_t)W * H * 4u;

    /* 4. uninterrupted tiled render is pixel-identical to one-shot */
    RenderJob one;
    check("tiled_init", job_init(&one, W, H, TW, TH) == 0);
    check("tiled_renders_all", render_step(&one, 1000) == 16);
    check("tiled_matches_oneshot", memcmp(one.fb, ref.fb, fb_bytes) == 0);

    /* 5. progress queryable at start, mid-run and completion */
    RenderJob prog;
    job_init(&prog, W, H, TW, TH);
    Progress p0 = job_progress(&prog);
    check("progress_zero", p0.tiles_done == 0 && p0.tiles_total == 16 &&
                           p0.tiles_x == 4 && p0.tiles_y == 4 &&
                           p0.fraction == 0.0f);
    check("step_budget", render_step(&prog, 3) == 3);
    Progress p1 = job_progress(&prog);
    check("progress_mid", p1.tiles_done == 3 &&
                          fabsf(p1.fraction - 3.0f / 16.0f) < 1e-6f);
    check("zero_budget_noop", render_step(&prog, 0) == 0);
    check("progress_unchanged", job_progress(&prog).tiles_done == 3);

    /* 6. abort at tile boundary: nothing renders, then resume continues */
    prog.abort_request = 1;
    check("abort_holds", render_step(&prog, 5) == 0);
    check("abort_no_progress", job_progress(&prog).tiles_done == 3);
    prog.abort_request = 0;
    check("resume_after_abort", render_step(&prog, 1000) == 13);
    Progress p2 = job_progress(&prog);
    check("progress_complete", p2.tiles_done == 16 &&
                               p2.fraction == 1.0f &&
                               job_is_complete(&prog));
    check("resumed_matches_oneshot", memcmp(prog.fb, ref.fb, fb_bytes) == 0);

    /* 7. finished tiles are idempotent and bytes stay put */
    check("step_when_done_zero", render_tile(&prog, 7) == 0 &&
                                 render_step(&prog, 10) == 0);
    check("idempotent_bytes", memcmp(prog.fb, ref.fb, fb_bytes) == 0);

    /* 8. interrupt -> save -> load (new process-equivalent) -> resume */
    RenderJob run;
    job_init(&run, W, H, TW, TH);
    check("partial_run", render_step(&run, 5) == 5);
    run.abort_request = 1;
    check("abort_stops_run", render_step(&run, 100) == 0);
    run.abort_request = 0;
    check("checkpoint_save", job_save(&run, CKPT_PATH) == 0);

    RenderJob loaded;
    check("checkpoint_load", job_load(&loaded, CKPT_PATH) == 0);
    check("loaded_progress", job_progress(&loaded).tiles_done == 5);
    check("loaded_bitmap", memcmp(loaded.done, run.done,
                                  (size_t)loaded.ntiles) == 0);
    check("loaded_pixels", memcmp(loaded.fb, run.fb, fb_bytes) == 0);
    check("loaded_resume", render_step(&loaded, 1000) == 11);
    check("loaded_is_complete", job_is_complete(&loaded));
    check("resumed_pixel_identical", memcmp(loaded.fb, ref.fb, fb_bytes) == 0);

    /* 9. corrupt / truncated checkpoints are rejected, never crash */
    FILE* jf = fopen(JUNK_PATH, "wb");
    if (jf) { fputs("not a checkpoint", jf); fclose(jf); }
    check("load_bad_magic_rejected", job_load(&bad, JUNK_PATH) == -1);

    int trunc_ok = 0;
    FILE* src = fopen(CKPT_PATH, "rb");
    FILE* dst = NULL;
    if (src) {
        fseek(src, 0, SEEK_END);
        long sz = ftell(src);
        fseek(src, 0, SEEK_SET);
        dst = fopen(TRUNC_PATH, "wb");
        if (dst && sz > 32) {
            char buf[256];
            long left = sz - 17;                 /* drop the last 17 bytes */
            trunc_ok = 1;
            while (left > 0) {
                size_t want = (size_t)(left < (long)sizeof(buf) ? left : (long)sizeof(buf));
                size_t got = fread(buf, 1, want, src);
                if (got == 0 || fwrite(buf, 1, got, dst) != got) { trunc_ok = 0; break; }
                left -= (long)got;
            }
        }
        fclose(src);
    }
    if (dst) fclose(dst);
    check("truncated_made", trunc_ok);
    check("load_truncated_rejected", job_load(&bad, TRUNC_PATH) == -1);

    remove(CKPT_PATH);
    remove(JUNK_PATH);
    remove(TRUNC_PATH);

    job_free(&ref);
    job_free(&one);
    job_free(&prog);
    job_free(&run);
    job_free(&loaded);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
