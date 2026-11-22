/* hybrid: baked lightmap + real-time direct light.
 *
 * A naive pipeline stacks baked direct + baked indirect + unshadowed
 * real-time direct.  Two defects follow:
 *   - baked direct is double counted (real-time adds it again),
 *   - real-time direct ignores occlusion, so light leaks into shadows.
 *
 * The hybrid path used here:
 *   1. drop baked direct, keep baked indirect (it is the only baked channel
 *      shipped in the final result),
 *   2. multiply real-time direct by the baked visibility,
 *   3. compensate the interpolation leak in proportion to the leak ratio
 *      lambda = 4 V (1 - V), which peaks where bilinear visibility is
 *      ambiguous (shadow boundary),
 *   4. compare against a pure real-time reference (per-pixel shadow rays,
 *      i.e. the accurate real-time baseline) and report an error upper
 *      bound from the triangle inequality,
 *   5. fall back to pure real-time where the bake is missing and annotate
 *      each such pixel.
 *
 * The indirect field is a stylized analytic model (lit term x exact AO);
 * the point of the module is the hybrid mechanics, not GI transport.
 * "Pure real-time" here means direct light with exact per-pixel visibility,
 * not an unshadowed pass: comparing against the naive path would be circular.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define GRID 64
#define TEXEL 8
#define NTEX (GRID / TEXEL)

#define NLIGHT 2
#define NBOX 3

#define AO_SAMPLES 32
#define AO_MAX_DIST 40.0f

#define IND_K 0.22f
#define LEAK_TAU 0.02f      /* direct light above this in a shadow counts as a leak */
static float g_c_leak = 1.0f;   /* leak-ratio compensation strength */
#define SHADOW_EPS 1e-3f

typedef struct { float x, y, z; } V3;
typedef struct { V3 lo, hi; } Box;
typedef struct { V3 pos; float intensity; } Light;

typedef struct {
    float direct[NLIGHT][NTEX][NTEX];   /* baked direct radiance at texel centers */
    float vis[NLIGHT][NTEX][NTEX];      /* baked visibility, 0 or 1 */
    float indirect[NTEX][NTEX];         /* baked indirect radiance */
    unsigned char valid[NTEX][NTEX];    /* 0 = bake missing */
} Lightmap;

typedef struct {
    Box boxes[NBOX];
    Light lights[NLIGHT];
    Lightmap lm;
    unsigned char fallback[GRID][GRID]; /* per-pixel bake-missing annotation */
} Scene;

typedef struct {
    int fallback_pixels;
    int invalid_texels;
} RenderStats;

typedef struct { int leak_pixels; int over_pixels; int shadow_pixels; double leaked_energy; } LeakStats;
typedef struct { double max_err; double mean_err; } ErrorStats;

/* ---------- math ---------- */

static float clampf(float v, float lo, float hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}
static V3 sub3(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

static V3 pixel_pos(int x, int z) {
    V3 p = { (float)x + 0.5f, 0.0f, (float)z + 0.5f };
    return p;
}
static V3 texel_center(int tx, int tz) {
    V3 p = { ((float)tx + 0.5f) * TEXEL, 0.0f, ((float)tz + 0.5f) * TEXEL };
    return p;
}

/* Slab test for a segment o + t*d, t in (eps, 1).  Returns 1 on hit. */
static int segment_hits_box(V3 o, V3 d, const Box* b) {
    float t0 = SHADOW_EPS, t1 = 1.0f - SHADOW_EPS;
    float ov[3] = { o.x, o.y, o.z }, dv[3] = { d.x, d.y, d.z };
    float lo[3] = { b->lo.x, b->lo.y, b->lo.z };
    float hi[3] = { b->hi.x, b->hi.y, b->hi.z };
    for (int i = 0; i < 3; i++) {
        if (fabsf(dv[i]) < 1e-8f) {
            if (ov[i] < lo[i] || ov[i] > hi[i]) return 0;
        } else {
            float inv = 1.0f / dv[i];
            float ta = (lo[i] - ov[i]) * inv;
            float tb = (hi[i] - ov[i]) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    return 1;
}

/* Ray test used by AO.  `d` is scaled to the AO range, so the unit-parameter
   segment test bounds the distance at AO_MAX_DIST. */
static int ray_hits_any_box(const Scene* s, V3 o, V3 d) {
    for (int i = 0; i < NBOX; i++)
        if (segment_hits_box(o, d, &s->boxes[i])) return 1;
    return 0;
}

static float halton(unsigned i, unsigned base) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= base; r += f * (i % base); i /= base; }
    return (float)r;
}

/* Cosine-weighted hemisphere occlusion, 0 = fully open, 1 = fully blocked. */
static float ao_at(const Scene* s, V3 p) {
    int occl = 0;
    V3 o = { p.x, p.y + 1e-4f, p.z };
    for (int i = 1; i <= AO_SAMPLES; i++) {
        float u1 = halton((unsigned)i, 2);
        float u2 = halton((unsigned)i, 3);
        float r = sqrtf(u1);
        float phi = 6.2831853f * u2;
        V3 dir = { r * cosf(phi) * AO_MAX_DIST,
                   sqrtf(clampf(1.0f - u1, 0.0f, 1.0f)) * AO_MAX_DIST,
                   r * sinf(phi) * AO_MAX_DIST };
        if (ray_hits_any_box(s, o, dir)) occl++;
    }
    return (float)occl / (float)AO_SAMPLES;
}

static void scene_init(Scene* s) {
    memset(s, 0, sizeof(*s));
    s->boxes[0].lo = (V3){ 18.0f, 0.0f, 14.0f };
    s->boxes[0].hi = (V3){ 30.0f, 12.0f, 22.0f };
    s->boxes[1].lo = (V3){ 42.0f, 0.0f, 34.0f };
    s->boxes[1].hi = (V3){ 50.0f, 12.0f, 46.0f };
    s->boxes[2].lo = (V3){ 8.0f, 0.0f, 52.0f };
    s->boxes[2].hi = (V3){ 56.0f, 10.0f, 56.0f };

    s->lights[0].pos = (V3){ 10.0f, 14.0f, 10.0f };
    s->lights[0].intensity = 500.0f;
    s->lights[1].pos = (V3){ 54.0f, 18.0f, 20.0f };
    s->lights[1].intensity = 800.0f;
}

/* Unshadowed direct geometry term for light l at p. */
static float light_geometry(const Scene* s, int l, V3 p) {
    V3 to = sub3(s->lights[l].pos, p);
    float d2 = dot3(to, to);
    float cosn = to.y / sqrtf(d2 + 1e-12f);
    if (cosn <= 0.0f) return 0.0f;
    return s->lights[l].intensity * cosn / (d2 + 1.0f);
}

static int light_shadowed(const Scene* s, int l, V3 p) {
    V3 d = sub3(s->lights[l].pos, p);
    for (int i = 0; i < NBOX; i++)
        if (segment_hits_box(p, d, &s->boxes[i])) return 1;
    return 0;
}

static float geometry_sum(const Scene* s, V3 p) {
    float g = 0.0f;
    for (int l = 0; l < NLIGHT; l++) g += light_geometry(s, l, p);
    return g;
}

static float direct_truth(const Scene* s, V3 p) {
    float d = 0.0f;
    for (int l = 0; l < NLIGHT; l++)
        if (!light_shadowed(s, l, p)) d += light_geometry(s, l, p);
    return d;
}

/* Stylized indirect: lit geometry term shaped by exact AO. */
static float indirect_truth(const Scene* s, V3 p) {
    return IND_K * (0.15f + 0.85f * (1.0f - ao_at(s, p))) * geometry_sum(s, p);
}

/* ---------- bake ---------- */

static void bake(Scene* s) {
    for (int tz = 0; tz < NTEX; tz++) {
        for (int tx = 0; tx < NTEX; tx++) {
            V3 c = texel_center(tx, tz);
            for (int l = 0; l < NLIGHT; l++) {
                int sh = light_shadowed(s, l, c);
                s->lm.vis[l][tx][tz] = sh ? 0.0f : 1.0f;
                s->lm.direct[l][tx][tz] = sh ? 0.0f : light_geometry(s, l, c);
            }
            s->lm.indirect[tx][tz] = indirect_truth(s, c);
            s->lm.valid[tx][tz] = 1;
        }
    }
}

static void invalidate_tile(Scene* s, int tx, int tz) {
    if (tx < 0 || tx >= NTEX || tz < 0 || tz >= NTEX) return;
    s->lm.valid[tx][tz] = 0;
}

/* Bilinear sample.  Returns 0 (and leaves *out untouched) if any of the four
   contributing texels is invalid, so the caller can fall back. */
static int sample_texel(const float g[NTEX][NTEX],
                        const unsigned char valid[NTEX][NTEX],
                        float fx, float fz, float* out) {
    int x0 = (int)floorf(fx), z0 = (int)floorf(fz);
    float ax = fx - (float)x0, az = fz - (float)z0;
    int ix0 = x0 < 0 ? 0 : (x0 >= NTEX ? NTEX - 1 : x0);
    int iz0 = z0 < 0 ? 0 : (z0 >= NTEX ? NTEX - 1 : z0);
    int ix1 = (x0 + 1) < 0 ? 0 : ((x0 + 1) >= NTEX ? NTEX - 1 : (x0 + 1));
    int iz1 = (z0 + 1) < 0 ? 0 : ((z0 + 1) >= NTEX ? NTEX - 1 : (z0 + 1));
    if (!valid[ix0][iz0] || !valid[ix1][iz0] || !valid[ix0][iz1] || !valid[ix1][iz1])
        return 0;
    float v00 = g[ix0][iz0], v10 = g[ix1][iz0];
    float v01 = g[ix0][iz1], v11 = g[ix1][iz1];
    float a = v00 * (1.0f - ax) + v10 * ax;
    float b = v01 * (1.0f - ax) + v11 * ax;
    *out = a * (1.0f - az) + b * az;
    return 1;
}

static void pixel_texel_coord(int x, int z, float* fx, float* fz) {
    *fx = ((float)x + 0.5f) / (float)TEXEL - 0.5f;
    *fz = ((float)z + 0.5f) / (float)TEXEL - 0.5f;
}

/* ---------- visibility correction and leak compensation ---------- */

/* Leak ratio: how much of a bilinear visibility sample came from blending
   texels with different values (0 at 0/1, 1 at 0.5). */
static float leak_ratio(float v) { return 4.0f * v * (1.0f - v); }

static float vis_compensated(float v) {
    float vc = v * (1.0f - g_c_leak * leak_ratio(v));
    return clampf(vc, 0.0f, 1.0f);
}

/* ---------- renderers ---------- */

static void render_reference(const Scene* s, float out[GRID][GRID]) {
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            V3 p = pixel_pos(x, z);
            out[x][z] = direct_truth(s, p) + indirect_truth(s, p);
        }
}

/* The bug: baked direct + baked indirect + unshadowed real-time direct. */
static void render_naive(const Scene* s, float out[GRID][GRID], float ind[GRID][GRID]) {
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            V3 p = pixel_pos(x, z);
            float fx, fz;
            pixel_texel_coord(x, z, &fx, &fz);
            float bd = 0.0f, bi = 0.0f, rt = 0.0f;
            for (int l = 0; l < NLIGHT; l++) {
                float v = 0.0f;
                (void)sample_texel(s->lm.direct[l], s->lm.valid, fx, fz, &v);
                bd += v;
                rt += light_geometry(s, l, p);
            }
            (void)sample_texel(s->lm.indirect, s->lm.valid, fx, fz, &bi);
            out[x][z] = bd + bi + rt;
            ind[x][z] = bi;
        }
}

static void render_hybrid(Scene* s, float out[GRID][GRID], float ind[GRID][GRID],
                          RenderStats* st) {
    memset(s->fallback, 0, sizeof(s->fallback));
    st->fallback_pixels = 0;
    st->invalid_texels = 0;
    for (int tz = 0; tz < NTEX; tz++)
        for (int tx = 0; tx < NTEX; tx++)
            if (!s->lm.valid[tx][tz]) st->invalid_texels++;

    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            V3 p = pixel_pos(x, z);
            float fx, fz;
            pixel_texel_coord(x, z, &fx, &fz);

            float bi = 0.0f;
            int ok = sample_texel(s->lm.indirect, s->lm.valid, fx, fz, &bi);
            float dir = 0.0f;
            if (ok) {
                for (int l = 0; l < NLIGHT; l++) {
                    float v = 0.0f;
                    if (!sample_texel(s->lm.vis[l], s->lm.valid, fx, fz, &v)) {
                        ok = 0;
                        break;
                    }
                    dir += light_geometry(s, l, p) * vis_compensated(clampf(v, 0.0f, 1.0f));
                }
            }
            if (!ok) {
                /* Bake missing: pure real-time for this pixel, annotated. */
                s->fallback[x][z] = 1;
                st->fallback_pixels++;
                ind[x][z] = indirect_truth(s, p);
                out[x][z] = direct_truth(s, p) + ind[x][z];
            } else {
                ind[x][z] = bi;
                out[x][z] = dir + bi;
            }
        }
}

/* ---------- metrics ---------- */

/* A lit pixel whose direct term lost more than OVER_TAU to the mask is
   over-shadowed (the opposite failure mode of a leak). */
#define OVER_TAU 0.02f

static LeakStats count_leaks(const Scene* s, const float out[GRID][GRID],
                             const float ind[GRID][GRID]) {
    LeakStats st;
    memset(&st, 0, sizeof(st));
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            V3 p = pixel_pos(x, z);
            float rd = direct_truth(s, p);
            float od = out[x][z] - ind[x][z];
            if (rd < SHADOW_EPS) {
                st.shadow_pixels++;
                if (od > LEAK_TAU) {
                    st.leak_pixels++;
                    st.leaked_energy += od;
                }
            } else if (od < rd - OVER_TAU) {
                st.over_pixels++;
            }
        }
    return st;
}

static ErrorStats compare(const float a[GRID][GRID], const float b[GRID][GRID]) {
    ErrorStats st;
    memset(&st, 0, sizeof(st));
    double sum = 0.0;
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            double e = fabs((double)a[x][z] - (double)b[x][z]);
            if (e > st.max_err) st.max_err = e;
            sum += e;
        }
    st.mean_err = sum / (double)(GRID * GRID);
    return st;
}

/* Per-pixel triangle-inequality upper bound on |hybrid - reference|.
 * Valid pixel:   |B_i - I|  + sum_l |D_l Vc_l - D_l s_l|
 * Fallback:     0 (hybrid reproduces reference exactly).
 */
static void hybrid_error_bound(const Scene* s, const float out_ind[GRID][GRID],
                               int with_fallback_mask, double bound[GRID * GRID]) {
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            if (with_fallback_mask && s->fallback[x][z]) {
                bound[x * GRID + z] = 0.0;
                continue;
            }
            V3 p = pixel_pos(x, z);
            float fx, fz;
            pixel_texel_coord(x, z, &fx, &fz);
            double bi = out_ind[x][z];
            double t = fabs(bi - indirect_truth(s, p));
            for (int l = 0; l < NLIGHT; l++) {
                float v = 0.0f;
                if (!sample_texel(s->lm.vis[l], s->lm.valid, fx, fz, &v)) continue;
                double d = light_geometry(s, l, p);
                double est = d * vis_compensated(clampf(v, 0.0f, 1.0f));
                double tru = light_shadowed(s, l, p) ? 0.0 : d;
                t += fabs(est - tru);
            }
            bound[x * GRID + z] = t;
        }
}

/* ---------- debug images ---------- */

static void write_ppm(const char* dir, const char* name, const float img[GRID][GRID],
                      float scale) {
    char path[512];
    snprintf(path, sizeof(path), "%s/%s.ppm", dir, name);
    FILE* f = fopen(path, "wb");
    if (!f) { printf("dump_failed=%s\n", path); return; }
    fprintf(f, "P6\n%d %d\n255\n", GRID, GRID);
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            float v = clampf(img[x][z] * scale, 0.0f, 1.0f);
            v = powf(v, 1.0f / 2.2f);
            unsigned char c = (unsigned char)(v * 255.0f + 0.5f);
            unsigned char px[3] = { c, c, c };
            fwrite(px, 1, 3, f);
        }
    fclose(f);
}

static void write_mask_ppm(const char* dir, const Scene* s) {
    char path[512];
    snprintf(path, sizeof(path), "%s/fallback_mask.ppm", dir);
    FILE* f = fopen(path, "wb");
    if (!f) { printf("dump_failed=%s\n", path); return; }
    fprintf(f, "P6\n%d %d\n255\n", GRID, GRID);
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            unsigned char px[3] = { 0, 0, 0 };
            if (s->fallback[x][z]) { px[0] = 255; px[1] = 0; px[2] = 255; } /* magenta */
            fwrite(px, 1, 3, f);
        }
    fclose(f);
}

/* ---------- self test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(int argc, char** argv) {
    const char* dump_dir = NULL;
    int sweep = 0;
    if (argc >= 3 && strcmp(argv[1], "--dump") == 0) dump_dir = argv[2];
    if (argc >= 2 && strcmp(argv[1], "--sweep") == 0) sweep = 1;

    Scene s;
    scene_init(&s);
    bake(&s);

    /* 1. bake sanity */
    int vis_binary = 1, vis_center_ok = 1;
    for (int tz = 0; tz < NTEX; tz++)
        for (int tx = 0; tx < NTEX; tx++)
            for (int l = 0; l < NLIGHT; l++) {
                float v = s.lm.vis[l][tx][tz];
                if (v != 0.0f && v != 1.0f) vis_binary = 0;
                V3 c = texel_center(tx, tz);
                if (fabsf(v - (light_shadowed(&s, l, c) ? 0.0f : 1.0f)) > 1e-6f)
                    vis_center_ok = 0;
            }
    check("bake_vis_binary", vis_binary);
    check("bake_vis_center_exact", vis_center_ok);

    int bilinear_ok = 1;
    for (int tz = 0; tz < NTEX; tz++)
        for (int tx = 0; tx < NTEX; tx++) {
            float fx = (float)tx, fz = (float)tz;
            float v = -1.0f;
            if (!sample_texel(s.lm.indirect, s.lm.valid, fx, fz, &v) ||
                fabsf(v - s.lm.indirect[tx][tz]) > 1e-5f) bilinear_ok = 0;
        }
    check("bilinear_texel_centers", bilinear_ok);

    /* leak ratio and compensation invariants */
    check("leak_ratio_endpoints",
          fabsf(leak_ratio(0.0f)) < 1e-6f && fabsf(leak_ratio(1.0f)) < 1e-6f);
    check("leak_ratio_peak", fabsf(leak_ratio(0.5f) - 1.0f) < 1e-6f);
    int comp_ok = 1;
    for (int i = 0; i <= 100; i++) {
        float v = (float)i / 100.0f;
        float vc = vis_compensated(v);
        if (vc < 0.0f || vc > 1.0f || vc > v + 1e-6f) comp_ok = 0;
    }
    /* compensation never brightens above the baked visibility */
    check("compensation_never_above_baked", comp_ok);
    check("compensation_endpoints",
          fabsf(vis_compensated(0.0f)) < 1e-6f &&
          fabsf(vis_compensated(1.0f) - 1.0f) < 1e-6f);

    /* 2. reference and naive */
    static float ref[GRID][GRID], naive[GRID][GRID], hyb[GRID][GRID];
    static float naive_ind[GRID][GRID], hyb_ind[GRID][GRID];

    render_reference(&s, ref);
    render_naive(&s, naive, naive_ind);

    int ref_shadow = 0;
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++)
            if (direct_truth(&s, pixel_pos(x, z)) < SHADOW_EPS) ref_shadow++;
    check("ref_has_shadow", ref_shadow > 50 && ref_shadow < GRID * GRID);


    if (sweep) {
        printf("c_leak,leak_pixels,over_pixels,leak_energy,mean_err,max_err\n");
        for (int i = 0; i <= 10; i++) {
            g_c_leak = (float)i / 10.0f;
            RenderStats ss;
            float so[GRID][GRID], si[GRID][GRID];
            render_hybrid(&s, so, si, &ss);
            LeakStats sl = count_leaks(&s, so, si);
            ErrorStats se = compare(so, ref);
            printf("%.1f,%d,%d,%.4f,%.4f,%.4f\n", (double)g_c_leak,
                   sl.leak_pixels, sl.over_pixels, sl.leaked_energy,
                   se.mean_err, se.max_err);
        }
        g_c_leak = 1.0f;
        return 0;
    }

    LeakStats nl = count_leaks(&s, naive, naive_ind);
    ErrorStats ne = compare(naive, ref);
    printf("shadow_pixels=%d\n", nl.shadow_pixels);
    printf("naive_leak_pixels=%d\n", nl.leak_pixels);
    printf("naive_leak_energy=%.4f\n", nl.leaked_energy);
    printf("naive_over_pixels=%d\n", nl.over_pixels);
    printf("naive_max_err=%.4f\n", ne.max_err);
    printf("naive_mean_err=%.4f\n", ne.mean_err);
    check("naive_leaks_pixels", nl.leak_pixels > 100);

    /* 3. hybrid, all bake valid */
    RenderStats hs;
    render_hybrid(&s, hyb, hyb_ind, &hs);
    LeakStats hl = count_leaks(&s, hyb, hyb_ind);
    ErrorStats he = compare(hyb, ref);
    double reduction = nl.leak_pixels > 0
        ? 100.0 * (1.0 - (double)hl.leak_pixels / (double)nl.leak_pixels) : 0.0;

    printf("hybrid_fallback_pixels=%d\n", hs.fallback_pixels);
    printf("hybrid_leak_pixels=%d\n", hl.leak_pixels);
    printf("hybrid_leak_energy=%.4f\n", hl.leaked_energy);
    printf("hybrid_over_pixels=%d\n", hl.over_pixels);
    printf("leak_reduction_pct=%.2f\n", reduction);
    printf("leak_energy_reduction_pct=%.2f\n",
           nl.leaked_energy > 0.0
               ? 100.0 * (1.0 - hl.leaked_energy / nl.leaked_energy) : 0.0);
    printf("hybrid_max_err=%.4f\n", he.max_err);
    printf("hybrid_mean_err=%.4f\n", he.mean_err);
    check("hybrid_no_fallback_when_valid", hs.fallback_pixels == 0);
    check("hybrid_leak_reduced", nl.leak_pixels > 0 &&
          hl.leak_pixels * 2 <= nl.leak_pixels);
    check("hybrid_mean_err_below_naive", he.mean_err < ne.mean_err);
    check("hybrid_max_err_below_naive", he.max_err < ne.max_err);

    /* 4. error upper bound vs pure real-time */
    static double bound[GRID * GRID];
    hybrid_error_bound(&s, hyb_ind, 0, bound);
    int bound_holds = 1;
    double bound_max = 0.0, err_max = 0.0;
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            double e = fabs((double)hyb[x][z] - (double)ref[x][z]);
            double b = bound[x * GRID + z];
            if (e > b + 1e-4) bound_holds = 0;
            if (b > bound_max) bound_max = b;
            if (e > err_max) err_max = e;
        }
    printf("hybrid_err_max_vs_realtime=%.4f\n", err_max);
    printf("hybrid_error_bound_max=%.4f\n", bound_max);
    check("error_bound_holds", bound_holds);
    check("error_bound_finite", bound_max < 2.0);

    /* 5. fallback when the bake is missing */
    static float valid_hyb[GRID][GRID];
    memcpy(valid_hyb, hyb, sizeof(hyb));
    invalidate_tile(&s, 3, 3);
    invalidate_tile(&s, 6, 1);
    RenderStats fs;
    float hyb_fb[GRID][GRID], hyb_fb_ind[GRID][GRID];
    render_hybrid(&s, hyb_fb, hyb_fb_ind, &fs);

    int annotated = 0, annotation_matches = 1, fallback_exact = 1, valid_kept = 1;
    double fb_err_max = 0.0;
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            if (s.fallback[x][z]) {
                annotated++;
                double e = fabs((double)hyb_fb[x][z] - (double)ref[x][z]);
                if (e > fb_err_max) fb_err_max = e;
                if (e > 1e-6) fallback_exact = 0;
            } else {
                if (fabs((double)hyb_fb[x][z] - (double)valid_hyb[x][z]) > 1e-6)
                    valid_kept = 0;
            }
        }
    if (annotated != fs.fallback_pixels) annotation_matches = 0;

    printf("invalid_texels=%d\n", fs.invalid_texels);
    printf("fallback_pixels=%d\n", fs.fallback_pixels);
    printf("fallback_annotated_pixels=%d\n", annotated);
    printf("fallback_max_err_vs_realtime=%.6f\n", fb_err_max);
    check("fallback_triggered", fs.fallback_pixels > 0 && fs.invalid_texels == 2);
    check("fallback_annotation_matches", annotation_matches);
    check("fallback_equals_realtime", fallback_exact);
    check("fallback_only_touched", valid_kept);

    /* fallback bound also holds */
    hybrid_error_bound(&s, hyb_fb_ind, 1, bound);
    int fb_bound_ok = 1;
    for (int z = 0; z < GRID; z++)
        for (int x = 0; x < GRID; x++) {
            double e = fabs((double)hyb_fb[x][z] - (double)ref[x][z]);
            if (e > bound[x * GRID + z] + 1e-4) fb_bound_ok = 0;
        }
    check("fallback_error_bound_holds", fb_bound_ok);

    if (dump_dir) {
        write_ppm(dump_dir, "reference", ref, 1.0f / 3.0f);
        write_ppm(dump_dir, "naive", naive, 1.0f / 3.0f);
        write_ppm(dump_dir, "hybrid", hyb_fb, 1.0f / 3.0f);
        write_mask_ppm(dump_dir, &s);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
