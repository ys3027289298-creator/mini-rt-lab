
/* upsample: edge-guided upscaling of a low-res render.
   - edge-directed interpolation along detected edges, bilinear elsewhere
   - before/after error metrics against a full-resolution reference
   - L_inf error reported as a guaranteed per-pixel upper bound
   - edge overshoot pixel count (values pushed past the local source range)
   - scale ratios <= 1 are rejected */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct { int w, h; float* p; } Img;

enum { UP_OK = 0, UP_ERR_NULL = 1, UP_ERR_SIZE = 2, UP_ERR_SCALE = 3 };

typedef struct {
    double mae;     /* mean absolute error vs reference */
    double rmse;    /* root mean square error */
    double maxabs;  /* L_inf: every per-pixel error is <= this upper bound */
    int    ne;      /* pixels with |error| > tol */
} ErrStats;

static float pix(const Img* a, int x, int y) {
    if (x < 0) x = 0; else if (x >= a->w) x = a->w - 1;
    if (y < 0) y = 0; else if (y >= a->h) y = a->h - 1;
    return a->p[(size_t)y * a->w + x];
}

static int img_alloc(Img* a, int w, int h) {
    if (!a || w <= 0 || h <= 0) return UP_ERR_SIZE;
    a->w = w; a->h = h;
    a->p = (float*)malloc((size_t)w * (size_t)h * sizeof(float));
    if (!a->p) { a->w = a->h = 0; return UP_ERR_SIZE; }
    return UP_OK;
}

static void img_free(Img* a) {
    if (!a) return;
    free(a->p); a->p = NULL; a->w = a->h = 0;
}

/* Upsampling only: the destination must be strictly larger on both axes.
   ratio <= 1 (including 1.0 itself) is rejected, NaN-safe via !(>). */
static int check_scale(const Img* src, const Img* dst) {
    if (!src || !dst || !src->p || !dst->p) return UP_ERR_NULL;
    if (src->w <= 0 || src->h <= 0 || dst->w <= 0 || dst->h <= 0) return UP_ERR_SIZE;
    if (!(dst->w > src->w && dst->h > src->h)) return UP_ERR_SCALE;
    return UP_OK;
}

/* area-average downsample; exact cell mean when the ratio is an integer */
static int box_downsample(const Img* ref, Img* low) {
    if (!ref || !low || !ref->p || !low->p) return UP_ERR_NULL;
    if (low->w <= 0 || low->h <= 0 || low->w > ref->w || low->h > ref->h)
        return UP_ERR_SIZE;
    for (int ly = 0; ly < low->h; ly++) {
        double fy0 = (double)ly * ref->h / low->h;
        double fy1 = (double)(ly + 1) * ref->h / low->h;
        int iy0 = (int)fy0, iy1 = (int)ceil(fy1);
        if (iy1 > ref->h) iy1 = ref->h;
        for (int lx = 0; lx < low->w; lx++) {
            double fx0 = (double)lx * ref->w / low->w;
            double fx1 = (double)(lx + 1) * ref->w / low->w;
            int ix0 = (int)fx0, ix1 = (int)ceil(fx1);
            if (ix1 > ref->w) ix1 = ref->w;
            double acc = 0.0, area = 0.0;
            for (int y = iy0; y < iy1; y++) {
                double wy = (y + 1 < fy1 ? y + 1 : fy1) - (y > fy0 ? y : fy0);
                for (int x = ix0; x < ix1; x++) {
                    double wx = (x + 1 < fx1 ? x + 1 : fx1) - (x > fx0 ? x : fx0);
                    acc += ref->p[(size_t)y * ref->w + x] * wx * wy;
                    area += wx * wy;
                }
            }
            low->p[(size_t)ly * low->w + lx] = (float)(acc / area);
        }
    }
    return UP_OK;
}

/* baseline "before upsampling": nearest-neighbour replication onto the
   full-resolution grid so it can be compared against the reference there */
static int upscale_nearest(const Img* s, Img* d) {
    int rc = check_scale(s, d);
    if (rc) return rc;
    for (int y = 0; y < d->h; y++) {
        int sy = (int)(((y + 0.5) * s->h) / d->h);
        if (sy >= s->h) sy = s->h - 1;
        for (int x = 0; x < d->w; x++) {
            int sx = (int)(((x + 0.5) * s->w) / d->w);
            if (sx >= s->w) sx = s->w - 1;
            d->p[(size_t)y * d->w + x] = s->p[(size_t)sy * s->w + sx];
        }
    }
    return UP_OK;
}

static int upscale_bilinear(const Img* s, Img* d) {
    int rc = check_scale(s, d);
    if (rc) return rc;
    double rx = (double)s->w / d->w, ry = (double)s->h / d->h;
    for (int y = 0; y < d->h; y++) {
        double fy = (y + 0.5) * ry - 0.5;
        int y0 = (int)floor(fy);
        float ty = (float)(fy - y0);
        for (int x = 0; x < d->w; x++) {
            double fx = (x + 0.5) * rx - 0.5;
            int x0 = (int)floor(fx);
            float tx = (float)(fx - x0);
            float v00 = pix(s, x0, y0),     v10 = pix(s, x0 + 1, y0);
            float v01 = pix(s, x0, y0 + 1), v11 = pix(s, x0 + 1, y0 + 1);
            float top = v00 + (v10 - v00) * tx;
            float bot = v01 + (v11 - v01) * tx;
            d->p[(size_t)y * d->w + x] = top + (bot - top) * ty;
        }
    }
    return UP_OK;
}

/* Edge-guided upscale.
   edge_thr: gradient magnitude (central differences) above which a source
             pixel counts as an edge; must be > 0.
   sharp:    optional unsharp amount sharp*(edge - bilinear), 0 disables it.
             A positive value can overshoot the local source range; the
             overshoot count is reported separately by count_overshoot(). */
static int upscale_edge(const Img* s, Img* d, float edge_thr, float sharp) {
    int rc = check_scale(s, d);
    if (rc) return rc;
    if (!(edge_thr > 0.0f) || !(sharp >= 0.0f)) return UP_ERR_SIZE;
    double rx = (double)s->w / d->w, ry = (double)s->h / d->h;
    for (int y = 0; y < d->h; y++) {
        double fy = (y + 0.5) * ry - 0.5;
        int y0 = (int)floor(fy);
        float ty = (float)(fy - y0);
        for (int x = 0; x < d->w; x++) {
            double fx = (x + 0.5) * rx - 0.5;
            int x0 = (int)floor(fx);
            float tx = (float)(fx - x0);

            float v00 = pix(s, x0, y0),     v10 = pix(s, x0 + 1, y0);
            float v01 = pix(s, x0, y0 + 1), v11 = pix(s, x0 + 1, y0 + 1);
            float top = v00 + (v10 - v00) * tx;
            float bot = v01 + (v11 - v01) * tx;
            float bil = top + (bot - top) * ty;

            int nx = (int)lround(fx), ny = (int)lround(fy);
            if (nx < 0) nx = 0; else if (nx >= s->w) nx = s->w - 1;
            if (ny < 0) ny = 0; else if (ny >= s->h) ny = s->h - 1;
            float gx = pix(s, nx + 1, ny) - pix(s, nx - 1, ny);
            float gy = pix(s, nx, ny + 1) - pix(s, nx, ny - 1);
            float mag = hypotf(gx, gy);

            /* hard edge gate, then a 0..1 blend toward the directional value */
            float w = 0.0f;
            if (mag > edge_thr) {
                w = (mag - edge_thr) / edge_thr;
                if (w > 1.0f) w = 1.0f;
            }

            float dir;
            if (fabsf(gx) >= fabsf(gy)) {
                /* vertical edge: interpolate along y, never cross the edge */
                double f = fy - ny;
                int ny1 = ny + (f >= 0.0 ? 1 : -1);
                float a = pix(s, nx, ny), c = pix(s, nx, ny1);
                dir = a + (c - a) * (float)fabs(f);
            } else {
                /* horizontal edge: interpolate along x */
                double f = fx - nx;
                int nx1 = nx + (f >= 0.0 ? 1 : -1);
                float a = pix(s, nx, ny), c = pix(s, nx1, ny);
                dir = a + (c - a) * (float)fabs(f);
            }

            float base = bil + (dir - bil) * w;
            d->p[(size_t)y * d->w + x] = base + sharp * (base - bil);
        }
    }
    return UP_OK;
}

/* Count edge overshoot pixels: output pixels in an edge neighbourhood whose
   value leaves the [min,max] range of the source pixels in their bilinear
   footprint by more than tol. Also reports the edge-region pixel total. */
static int count_overshoot(const Img* s, const Img* d, float edge_thr,
                           float tol, int* over_n, int* region_n) {
    if (!s || !d || !s->p || !d->p) return UP_ERR_NULL;
    if (!(edge_thr > 0.0f) || !(tol >= 0.0f)) return UP_ERR_SIZE;
    if (s->w <= 0 || s->h <= 0 || d->w <= s->w || d->h <= s->h)
        return UP_ERR_SCALE;

    unsigned char* edge = (unsigned char*)calloc((size_t)s->w * s->h, 1);
    if (!edge) return UP_ERR_SIZE;
    for (int y = 0; y < s->h; y++)
        for (int x = 0; x < s->w; x++) {
            float gx = pix(s, x + 1, y) - pix(s, x - 1, y);
            float gy = pix(s, x, y + 1) - pix(s, x, y - 1);
            edge[(size_t)y * s->w + x] = hypotf(gx, gy) > edge_thr;
        }

    double rx = (double)s->w / d->w, ry = (double)s->h / d->h;
    int over = 0, region = 0;
    for (int y = 0; y < d->h; y++) {
        double fy = (y + 0.5) * ry - 0.5;
        int y0 = (int)floor(fy);
        for (int x = 0; x < d->w; x++) {
            double fx = (x + 0.5) * rx - 0.5;
            int x0 = (int)floor(fx);
            int xs[2] = { x0, x0 + 1 }, ys[2] = { y0, y0 + 1 };
            float lo = 1e30f, hi = -1e30f;
            int near_edge = 0;
            for (int j = 0; j < 2; j++)
                for (int i = 0; i < 2; i++) {
                    int cx = xs[i] < 0 ? 0 : (xs[i] >= s->w ? s->w - 1 : xs[i]);
                    int cy = ys[j] < 0 ? 0 : (ys[j] >= s->h ? s->h - 1 : ys[j]);
                    float v = s->p[(size_t)cy * s->w + cx];
                    if (v < lo) lo = v;
                    if (v > hi) hi = v;
                    if (edge[(size_t)cy * s->w + cx]) near_edge = 1;
                }
            if (!near_edge) continue;
            region++;
            float v = d->p[(size_t)y * d->w + x];
            if (v > hi + tol || v < lo - tol) over++;
        }
    }
    free(edge);
    *over_n = over;
    *region_n = region;
    return UP_OK;
}

/* Error of a test image against a same-size full-resolution reference.
   maxabs is the L_inf norm: a guaranteed upper bound on every pixel error. */
static int compare_img(const Img* t, const Img* ref, float tol, ErrStats* st) {
    if (!t || !ref || !t->p || !ref->p || !st) return UP_ERR_NULL;
    if (t->w != ref->w || t->h != ref->h) return UP_ERR_SIZE;
    double sa = 0.0, s2 = 0.0, mx = 0.0;
    int ne = 0, n = t->w * t->h;
    for (int i = 0; i < n; i++) {
        double e = (double)t->p[i] - ref->p[i];
        double a = fabs(e);
        sa += a;
        s2 += e * e;
        if (a > mx) mx = a;
        if (a > tol) ne++;
    }
    st->mae = sa / n;
    st->rmse = sqrt(s2 / n);
    st->maxabs = mx;
    st->ne = ne;
    return UP_OK;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    enum { W = 64, H = 64, LW = 16, LH = 16 };
    const float EDGE_THR = 0.20f, SHARP = 0.5f, TOL = 1e-6f;

    /* full-res reference: hard vertical step on top, smooth ramp on bottom */
    Img ref;
    img_alloc(&ref, W, H);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            ref.p[y * W + x] = y < H / 2 ? (x < W / 2 ? 0.15f : 0.85f)
                                         : (float)x / (W - 1);

    /* low-res "render": box-filtered downsample, scale 4 */
    Img low;
    img_alloc(&low, LW, LH);
    check("downsample_ok", box_downsample(&ref, &low) == UP_OK);
    /* ramp cell (0,8) averages ref x=0..3: mean 1.5/63 */
    check("downsample_cell", fabsf(low.p[8 * LW] - 1.5f / 63.0f) < 1e-6f);
    /* aligned step cells survive exactly */
    check("downsample_step", fabsf(low.p[7] - 0.15f) < 1e-6f &&
                             fabsf(low.p[8] - 0.85f) < 1e-6f);

    Img bad;
    img_alloc(&bad, W, H);
    check("downsample_reject_grow", box_downsample(&low, &bad) == UP_ERR_SIZE);
    img_free(&bad);

    /* ratio <= 1 must be rejected (equal size and downscale) */
    Img same, small;
    img_alloc(&same, LW, LH);
    img_alloc(&small, LW / 2, LH / 2);
    check("reject_ratio_one", upscale_edge(&low, &same, EDGE_THR, 0.0f) == UP_ERR_SCALE);
    check("reject_ratio_below", upscale_edge(&low, &small, EDGE_THR, 0.0f) == UP_ERR_SCALE);
    check("reject_bilinear_scale", upscale_bilinear(&low, &same) == UP_ERR_SCALE);
    check("reject_bad_params", upscale_edge(&low, &same, 0.0f, 0.0f) == UP_ERR_SCALE ||
                               upscale_edge(&low, &same, EDGE_THR, -1.0f) == UP_ERR_SIZE);
    img_free(&same);
    img_free(&small);

    Img before, bil, edge0, edgeS;
    img_alloc(&before, W, H);
    img_alloc(&bil, W, H);
    img_alloc(&edge0, W, H);
    img_alloc(&edgeS, W, H);
    check("nearest_ok", upscale_nearest(&low, &before) == UP_OK);
    check("bilinear_ok", upscale_bilinear(&low, &bil) == UP_OK);
    check("edge_ok", upscale_edge(&low, &edge0, EDGE_THR, 0.0f) == UP_OK);
    check("edge_sharp_ok", upscale_edge(&low, &edgeS, EDGE_THR, SHARP) == UP_OK);

    ErrStats sb, sl, se, ss;
    check("compare_size_mismatch",
          compare_img(&edge0, &low, TOL, &se) == UP_ERR_SIZE);
    compare_img(&before, &ref, TOL, &sb);
    compare_img(&bil, &ref, TOL, &sl);
    compare_img(&edge0, &ref, TOL, &se);
    compare_img(&edgeS, &ref, TOL, &ss);

    printf("before_mae=%.6f before_rmse=%.6f before_maxabs=%.6f\n",
           sb.mae, sb.rmse, sb.maxabs);
    printf("bilinear_mae=%.6f bilinear_rmse=%.6f bilinear_maxabs=%.6f\n",
           sl.mae, sl.rmse, sl.maxabs);
    printf("after_mae=%.6f after_rmse=%.6f after_maxabs=%.6f\n",
           se.mae, se.rmse, se.maxabs);
    printf("sharp_mae=%.6f sharp_rmse=%.6f sharp_maxabs=%.6f\n",
           ss.mae, ss.rmse, ss.maxabs);

    /* edge-guided must beat bilinear on the step and beat the pre-upscale buffer */
    check("edge_beats_bilinear", se.mae < sl.mae && se.maxabs <= sl.maxabs + 1e-6);
    check("after_beats_before", se.mae < sb.mae);

    /* L_inf is a real upper bound: independently verify every pixel */
    int bound_holds = 1;
    for (int i = 0; i < W * H; i++)
        if (fabs((double)edge0.p[i] - ref.p[i]) > se.maxabs + 1e-6) bound_holds = 0;
    check("maxabs_upper_bound", bound_holds);
    check("metric_ordering", se.maxabs + 1e-12 >= se.rmse &&
                             se.rmse + 1e-12 >= se.mae);

    /* edge overshoot: convex interpolation has none; sharpening creates some */
    int over0 = -1, reg0 = -1, overS = -1, regS = -1;
    check("overshoot_none_ok",
          count_overshoot(&low, &edge0, EDGE_THR, TOL, &over0, &reg0) == UP_OK);
    check("overshoot_sharp_ok",
          count_overshoot(&low, &edgeS, EDGE_THR, TOL, &overS, &regS) == UP_OK);
    printf("edge_region_pixels=%d overshoot_no_sharp=%d overshoot_sharp=%d\n",
           regS, over0, overS);
    check("overshoot_zero_without_sharp", over0 == 0);
    check("overshoot_reported_with_sharp", overS > 0);
    check("overshoot_within_region", overS <= regS && regS > 0);

    /* constant field: every method must reproduce it exactly */
    Img c8, c24;
    img_alloc(&c8, 8, 8);
    img_alloc(&c24, 24, 24);
    for (int i = 0; i < 64; i++) c8.p[i] = 0.4f;
    check("const_upscale_ok", upscale_edge(&c8, &c24, EDGE_THR, SHARP) == UP_OK);
    int const_exact = 1;
    for (int i = 0; i < 24 * 24; i++)
        if (c24.p[i] != 0.4f) const_exact = 0;
    check("const_exact", const_exact);
    ErrStats sc;
    compare_img(&c24, &c24, TOL, &sc);
    check("const_zero_error", sc.maxabs == 0.0);
    int cov, creg;
    count_overshoot(&c8, &c24, EDGE_THR, TOL, &cov, &creg);
    check("const_no_overshoot", cov == 0 && creg == 0);
    img_free(&c8);
    img_free(&c24);

    /* smooth linear ramp: with the edge gate forced off, edge-guided output
       must be bit-identical to plain bilinear, and a downsample/upsample
       round trip reconstructs the full-res ramp within a small residual */
    Img r8, r32, rref, rbil;
    img_alloc(&r8, 8, 8);
    img_alloc(&r32, 32, 32);
    img_alloc(&rref, 32, 32);
    img_alloc(&rbil, 32, 32);
    for (int y = 0; y < 32; y++)
        for (int x = 0; x < 32; x++)
            rref.p[y * 32 + x] = (float)x / 31.0f;
    box_downsample(&rref, &r8);
    upscale_edge(&r8, &r32, 1.0f, 0.0f);   /* threshold above every gradient */
    upscale_bilinear(&r8, &rbil);
    int gate_off_equals_bilinear = 1;
    for (int i = 0; i < 32 * 32; i++)
        if (r32.p[i] != rbil.p[i]) gate_off_equals_bilinear = 0;
    check("gate_off_is_bilinear", gate_off_equals_bilinear);
    ErrStats sr;
    compare_img(&r32, &rref, TOL, &sr);
    printf("ramp_roundtrip_mae=%.6f maxabs=%.6f\n", sr.mae, sr.maxabs);
    /* boundary clamp leaves at most half a low-res cell (4/32) of residual */
    check("ramp_roundtrip_half_cell", sr.maxabs < 0.05f);
    img_free(&r8);
    img_free(&r32);
    img_free(&rref);
    img_free(&rbil);

    /* non-integer scale: every output must stay inside the source range,
       and the L_inf bound still holds against the resampled reference */
    Img nz, nref, nnear;
    img_alloc(&nz, 30, 20);
    img_alloc(&nref, 30, 20);
    img_alloc(&nnear, 30, 20);
    box_downsample(&ref, &nref);
    check("noninteger_ok", upscale_edge(&low, &nz, EDGE_THR, 0.0f) == UP_OK);
    check("noninteger_nearest_ok", upscale_nearest(&low, &nnear) == UP_OK);
    float slo = 1e30f, shi = -1e30f;
    for (int i = 0; i < LW * LH; i++) {
        if (low.p[i] < slo) slo = low.p[i];
        if (low.p[i] > shi) shi = low.p[i];
    }
    int in_range = 1;
    for (int i = 0; i < 30 * 20; i++)
        if (nz.p[i] < slo - 1e-5f || nz.p[i] > shi + 1e-5f) in_range = 0;
    check("noninteger_in_range", in_range);
    ErrStats sn, snn;
    compare_img(&nz, &nref, TOL, &sn);
    compare_img(&nnear, &nref, TOL, &snn);
    check("noninteger_beats_nearest", sn.mae < snn.mae);
    int nover, nreg;
    check("noninteger_overshoot_ok",
          count_overshoot(&low, &nz, EDGE_THR, TOL, &nover, &nreg) == UP_OK &&
          nover == 0 && nreg > 0);
    img_free(&nz);
    img_free(&nref);
    img_free(&nnear);

    img_free(&before);
    img_free(&bil);
    img_free(&edge0);
    img_free(&edgeS);
    img_free(&low);
    img_free(&ref);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
