
/* atlas: area-sorted shelf packing of many small textures into one atlas,
   so a scene draws in 1 batch instead of 1 bind per texture.
   Verifies the atlas render is pixel-identical to per-texture rendering
   and reports utilization / waste / batch counts. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { unsigned char r, g, b; } RGB;

typedef struct {
    int id;
    int w, h;
    RGB* px;    /* w*h texels */
    int ax, ay; /* atlas position, set by pack() */
} Tex;

typedef struct { int w, h; RGB* px; } Image;

static RGB* img_alloc(int w, int h) { return (RGB*)calloc((size_t)w * (size_t)h, sizeof(RGB)); }

/* deterministic procedural content: per-texture hue + checker + gradient */
static void gen_tex(Tex* t) {
    for (int y = 0; y < t->h; y++)
        for (int x = 0; x < t->w; x++) {
            int c = ((x >> 3) + (y >> 3)) & 1;
            RGB* p = &t->px[(size_t)y * t->w + x];
            p->r = (unsigned char)((t->id * 47 + x * 255 / (t->w > 1 ? t->w - 1 : 1) / 2 + (c ? 80 : 0)) & 0xFF);
            p->g = (unsigned char)((t->id * 91 + y * 255 / (t->h > 1 ? t->h - 1 : 1) / 2 + (c ? 40 : 0)) & 0xFF);
            p->b = (unsigned char)((t->id * 29 + (x + y) + (c ? 120 : 0)) & 0xFF);
        }
}

/* ---- shelf packer, input sorted by area descending ---- */

static int cmp_area_desc(const void* a, const void* b) {
    const Tex* ta = *(const Tex* const*)a;
    const Tex* tb = *(const Tex* const*)b;
    long da = (long)tb->w * tb->h - (long)ta->w * ta->h;
    if (da != 0) return da < 0 ? -1 : 1;
    return tb->h - ta->h; /* tie-break: taller first, deterministic */
}

/* shelf-pack n textures into an atlas of width atlas_w, gutters of `pad`.
   Returns 0 and the used height via out_h, or -1 if any texture is
   zero-sized (w <= 0 or h <= 0): zero-size input is rejected, not packed. */
static int pack(Tex** ts, int n, int pad, int atlas_w, int* out_h) {
    for (int i = 0; i < n; i++)
        if (ts[i]->w <= 0 || ts[i]->h <= 0) return -1;
    int x = pad, y = pad, shelf_h = 0;
    for (int i = 0; i < n; i++) {
        Tex* t = ts[i];
        if (x + t->w + pad > atlas_w) { /* close shelf, start a new one */
            x = pad;
            y += shelf_h + pad;
            shelf_h = 0;
        }
        t->ax = x;
        t->ay = y;
        x += t->w + pad;
        if (t->h > shelf_h) shelf_h = t->h;
    }
    *out_h = y + shelf_h + pad;
    return 0;
}

static int rects_overlap(const Tex* a, const Tex* b) {
    return a->ax < b->ax + b->w && b->ax < a->ax + a->w &&
           a->ay < b->ay + b->h && b->ay < a->ay + a->h;
}

/* ---- sampling: identical texel indexing on both paths ---- */

static int texel_i(float u, int n) {
    int i = (int)(u * (float)n);
    return i < n ? i : n - 1; /* u in [0,1), clamp the u==1 edge */
}

static RGB sample_individual(const Tex* t, float u, float v) {
    return t->px[(size_t)texel_i(v, t->h) * t->w + texel_i(u, t->w)];
}

static RGB sample_atlas(const Image* a, const Tex* t, float u, float v) {
    int x = t->ax + texel_i(u, t->w);
    int y = t->ay + texel_i(v, t->h);
    return a->px[(size_t)y * a->w + x];
}

/* ---- scene: grid of quads, one texture per quad ---- */

#define COLS 4
#define ROWS 3

static void render(Image* fb, Tex** ts, int n, const Image* atlas, int* batches) {
    int cw = fb->w / COLS, ch = fb->h / ROWS;
    *batches = 0;
    if (atlas) (*batches)++; /* one atlas bind covers every quad */
    for (int q = 0; q < n && q < COLS * ROWS; q++) {
        const Tex* t = ts[q];
        if (!atlas) (*batches)++; /* one bind per texture */
        int x0 = (q % COLS) * cw, y0 = (q / COLS) * ch;
        for (int y = 0; y < ch; y++)
            for (int x = 0; x < cw; x++) {
                float u = (float)x / (float)cw;
                float v = (float)y / (float)ch;
                fb->px[(size_t)(y0 + y) * fb->w + (x0 + x)] =
                    atlas ? sample_atlas(atlas, t, u, v)
                          : sample_individual(t, u, v);
            }
    }
}

static void write_ppm(const char* path, const Image* im) {
    FILE* f = fopen(path, "wb");
    if (!f) return;
    fprintf(f, "P6\n%d %d\n255\n", im->w, im->h);
    fwrite(im->px, sizeof(RGB), (size_t)im->w * im->h, f);
    fclose(f);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* 12 textures, mixed sizes, none power-of-two-aligned on purpose */
    int sizes[][2] = {
        { 64, 64 }, { 48, 48 }, { 80, 32 }, { 32, 80 },
        { 56, 40 }, { 24, 24 }, { 96, 16 }, { 40, 56 },
        { 16, 96 }, { 72, 72 }, { 36, 28 }, { 20, 44 },
    };
    const int n = (int)(sizeof(sizes) / sizeof(sizes[0]));
    const int pad = 1;

    Tex* tex = (Tex*)calloc(n, sizeof(Tex));
    Tex** sorted = (Tex**)malloc((size_t)n * sizeof(Tex*));
    long used_area = 0;
    for (int i = 0; i < n; i++) {
        tex[i].id = i;
        tex[i].w = sizes[i][0];
        tex[i].h = sizes[i][1];
        tex[i].px = img_alloc(tex[i].w, tex[i].h);
        tex[i].ax = tex[i].ay = -1;
        gen_tex(&tex[i]);
        sorted[i] = &tex[i];
        used_area += (long)tex[i].w * tex[i].h;
    }

    /* try candidate widths around sqrt(used area), keep the tightest atlas */
    qsort(sorted, n, sizeof(Tex*), cmp_area_desc); /* pack by area, largest first */
    int wmin = (int)ceil(sqrt((double)used_area));
    int atlas_w = 0, atlas_h = 0;
    long best_area = 0;
    for (int w = wmin; w <= 2 * wmin; w++) {
        int h = 0;
        if (pack(sorted, n, pad, w, &h) != 0) continue;
        long area = (long)w * h;
        if (!best_area || area < best_area) { best_area = area; atlas_w = w; atlas_h = h; }
    }
    check("pack_ok", best_area > 0);
    pack(sorted, n, pad, atlas_w, &atlas_h); /* repack at the chosen width */

    int all_placed = 1, any_overlap = 0;
    for (int i = 0; i < n; i++) {
        if (tex[i].ax < 0 || tex[i].ay < 0) all_placed = 0;
        for (int j = i + 1; j < n; j++)
            if (rects_overlap(&tex[i], &tex[j])) any_overlap = 1;
    }
    check("all_placed", all_placed);
    check("no_overlap", !any_overlap);

    /* zero-size textures must be rejected, not packed */
    Tex bad1 = { .id = 100, .w = 0, .h = 16, .px = NULL };
    Tex bad2 = { .id = 101, .w = 16, .h = 0, .px = NULL };
    Tex bad3 = { .id = 102, .w = 0, .h = 0, .px = NULL };
    Tex* with_bad[2];
    int h_dummy = 0;
    with_bad[0] = &tex[0]; with_bad[1] = &bad1;
    int rej1 = pack(with_bad, 2, pad, atlas_w, &h_dummy) == -1;
    with_bad[1] = &bad2;
    int rej2 = pack(with_bad, 2, pad, atlas_w, &h_dummy) == -1;
    with_bad[1] = &bad3;
    int rej3 = pack(with_bad, 2, pad, atlas_w, &h_dummy) == -1;
    check("zero_size_rejected", rej1 && rej2 && rej3);

    /* blit into the atlas image */
    Image atlas = { atlas_w, atlas_h, img_alloc(atlas_w, atlas_h) };
    for (int i = 0; i < n; i++) {
        Tex* t = &tex[i];
        for (int y = 0; y < t->h; y++)
            memcpy(&atlas.px[(size_t)(t->ay + y) * atlas.w + t->ax],
                   &t->px[(size_t)y * t->w], (size_t)t->w * sizeof(RGB));
    }

    /* render the same scene both ways */
    Image fb_ind = { COLS * 128, ROWS * 96, NULL };
    Image fb_atl = { COLS * 128, ROWS * 96, NULL };
    fb_ind.px = img_alloc(fb_ind.w, fb_ind.h);
    fb_atl.px = img_alloc(fb_atl.w, fb_atl.h);
    int batches_ind = 0, batches_atl = 0;
    render(&fb_ind, sorted, n, NULL, &batches_ind);
    render(&fb_atl, sorted, n, &atlas, &batches_atl);

    check("pixel_identical",
          memcmp(fb_ind.px, fb_atl.px, (size_t)fb_ind.w * fb_ind.h * sizeof(RGB)) == 0);
    check("batches_reduced", batches_atl < batches_ind);

    /* ---- data: utilization, waste, batches ---- */
    long atlas_area = (long)atlas.w * atlas.h;
    double util = (double)used_area / (double)atlas_area;
    printf("textures=%d\n", n);
    printf("atlas_w=%d\n", atlas.w);
    printf("atlas_h=%d\n", atlas.h);
    printf("used_area=%ld\n", used_area);
    printf("atlas_area=%ld\n", atlas_area);
    printf("utilization=%.4f\n", util);
    printf("waste=%.4f\n", 1.0 - util);
    printf("batches_per_texture=%d\n", batches_ind);
    printf("batches_atlas=%d\n", batches_atl);
    printf("batch_reduction=%.4f\n", 1.0 - (double)batches_atl / (double)batches_ind);
    check("utilization_above_half", util > 0.5);

    /* ---- visual comparison artifacts ---- */
    write_ppm("atlas.ppm", &atlas);
    write_ppm("render_individual.ppm", &fb_ind);
    write_ppm("render_atlas.ppm", &fb_atl);

    /* side-by-side: individual | atlas, with amplified diff strip below */
    Image cmp = { fb_ind.w * 2 + 8, fb_ind.h + 8 + 48, img_alloc(fb_ind.w * 2 + 8, fb_ind.h + 56) };
    for (int y = 0; y < fb_ind.h; y++) {
        memcpy(&cmp.px[(size_t)y * cmp.w], &fb_ind.px[(size_t)y * fb_ind.w],
               (size_t)fb_ind.w * sizeof(RGB));
        memcpy(&cmp.px[(size_t)y * cmp.w + fb_ind.w + 8], &fb_atl.px[(size_t)y * fb_atl.w],
               (size_t)fb_atl.w * sizeof(RGB));
    }
    for (int y = 0; y < 48; y++) /* diff x16: identical renders -> black strip */
        for (int x = 0; x < fb_ind.w; x++) {
            RGB a = fb_ind.px[(size_t)y * fb_ind.w + x];
            RGB b = fb_atl.px[(size_t)y * fb_atl.w + x];
            RGB d = { (unsigned char)(abs(a.r - b.r) * 16),
                      (unsigned char)(abs(a.g - b.g) * 16),
                      (unsigned char)(abs(a.b - b.b) * 16) };
            cmp.px[(size_t)(fb_ind.h + 8 + y) * cmp.w + x] = d;
        }
    write_ppm("compare.ppm", &cmp);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
