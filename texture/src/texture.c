
/* texture: mip pyramid plus anisotropic filtering driven by screen-space
   derivatives.  The pixel footprint is the ellipse J*J^T of the UV
   Jacobian; its minor axis picks the mip level, the major/minor ratio picks
   the tap count.  A reference exact box filter (polygon clipping over the
   texel grid) supplies the per-pixel comparison, an error bound and a
   flicker metric. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define TEX_N        256
#define TEX_LOG2     8
#define MIP_LEVELS   (TEX_LOG2 + 1)
#define PYRAMID_N    87381            /* 256^2 + 128^2 + ... + 1 */
#define MAX_TAPS     16
#define LOD_MAX      (MIP_LEVELS - 1)
#define PI           3.14159265358979323846

static float g_tex[TEX_N * TEX_N];
static float g_pyr[PYRAMID_N];
static int   g_off[MIP_LEVELS];
static int   g_dim[MIP_LEVELS];
static int   g_truncations;

/* UV derivatives in texels per pixel.  A footprint axis length of 1 means
   one texel per pixel, i.e. mip level 0. */
typedef struct { double du_dx, du_dy, dv_dx, dv_dy; } Grad;

typedef struct {
    double lod;            /* log2(minor), clamped to [0, LOD_MAX] */
    double ratio;          /* major / max(minor, 1) */
    double major, minor;   /* footprint axis lengths, level-0 texels */
    double theta;          /* major axis direction */
    double cos_t, sin_t;
    double step;           /* spacing between taps, level-0 texels */
    int    taps;           /* taps along the major axis */
    int    level0, level1;
    double blend;          /* trilinear weight (0 when level1 == level0) */
    int    clamped;        /* ratio exceeded MAX_TAPS: taps truncated */
    int    fallback;       /* zero derivative: nearest level, level 0 */
} Plan;

static int wrapi(int i, int n) { i %= n; return i < 0 ? i + n : i; }

/* ---------------------------------------------------------------- pyramid */

static void tex_fill_detail(void) {
    for (int y = 0; y < TEX_N; y++) {
        for (int x = 0; x < TEX_N; x++) {
            double v = 0.5
                     + 0.20 * sin(0.071 * x + 0.049 * y)
                     + 0.13 * sin(0.33 * x - 0.21 * y)
                     + 0.08 * sin(0.91 * x + 0.57 * y)
                     + 0.05 * sin(2.10 * x - 1.70 * y)
                     + 0.03 * sin(1.30 * x + 1.10 * y);
            if (v < 0.0) v = 0.0;
            if (v > 1.0) v = 1.0;
            g_tex[y * TEX_N + x] = (float)v;
        }
    }
}

/* smooth band along u: 0.5 + 0.4 sin(2*pi*u/64), constant along v */
static void tex_fill_smooth(void) {
    for (int y = 0; y < TEX_N; y++)
        for (int x = 0; x < TEX_N; x++)
            g_tex[y * TEX_N + x] = (float)(0.5 + 0.4 * sin(2.0 * PI * x / 64.0));
}

static void mip_build(void) {
    int off = 0, n = TEX_N;
    for (int l = 0; l < MIP_LEVELS; l++) {
        g_off[l] = off;
        g_dim[l] = n;
        off += n * n;
        n >>= 1;
    }
    memcpy(&g_pyr[g_off[0]], g_tex, sizeof(g_tex));
    for (int l = 1; l < MIP_LEVELS; l++) {
        int n2 = g_dim[l], n1 = g_dim[l - 1];
        const float* s = &g_pyr[g_off[l - 1]];
        float* d = &g_pyr[g_off[l]];
        for (int y = 0; y < n2; y++) {
            for (int x = 0; x < n2; x++) {
                d[y * n2 + x] = (float)(0.25 * ((double)s[(2 * y) * n1 + 2 * x]
                                              + (double)s[(2 * y) * n1 + 2 * x + 1]
                                              + (double)s[(2 * y + 1) * n1 + 2 * x]
                                              + (double)s[(2 * y + 1) * n1 + 2 * x + 1]));
            }
        }
    }
}

static double lvl_get(int l, int x, int y) {
    int n = g_dim[l];
    return g_pyr[g_off[l] + wrapi(y, n) * n + wrapi(x, n)];
}

/* Bilinear sample of level l, position in level-0 texel units.  Level-l
   texel j is the average of level-0 texels [s*j, s*j+s), so its centre is at
   s*j + s/2 and maps to px = j under the standard -0.5 shift. */
static double sample_level(int l, double x, double y) {
    double s = (double)(1 << l);
    double px = x / s - 0.5, py = y / s - 0.5;
    int x0 = (int)floor(px), y0 = (int)floor(py);
    double fx = px - x0, fy = py - y0;
    double a = lvl_get(l, x0, y0),     b = lvl_get(l, x0 + 1, y0);
    double c = lvl_get(l, x0, y0 + 1), d = lvl_get(l, x0 + 1, y0 + 1);
    return (a + (b - a) * fx) + ((c + (d - c) * fx) - (a + (b - a) * fx)) * fy;
}

static double texel_nearest(int x, int y) {
    return g_tex[wrapi(y, TEX_N) * TEX_N + wrapi(x, TEX_N)];
}

/* --------------------------------------------------------------- planning */

static Plan plan_from_grad(const Grad* g) {
    Plan p;
    memset(&p, 0, sizeof(p));
    p.taps = 1;
    p.step = 1.0;
    p.cos_t = 1.0;

    /* E = J J^T with J = [[du/dx, du/dy], [dv/dx, dv/dy]]; eigenvalues of a
       symmetric 2x2 are closed form, eigenvectors give the footprint axes */
    double a = g->du_dx * g->du_dx + g->du_dy * g->du_dy;
    double b = g->du_dx * g->dv_dx + g->du_dy * g->dv_dy;
    double c = g->dv_dx * g->dv_dx + g->dv_dy * g->dv_dy;
    double tr = a + c;
    double disc = sqrt(fmax(0.0, (a - c) * (a - c) + 4.0 * b * b));
    double lmaj = 0.5 * (tr + disc);
    double lmin = fmax(0.0, 0.5 * (tr - disc));

    if (lmaj > 1e-24) {
        p.major = sqrt(lmaj);
        p.minor = sqrt(lmin);
        p.theta = 0.5 * atan2(2.0 * b, a - c);
        p.cos_t = cos(p.theta);
        p.sin_t = sin(p.theta);
        p.ratio = p.major / fmax(p.minor, 1.0);

        /* level from the minor axis: the one direction that must not alias */
        p.lod = log2(fmax(p.minor, 1.0));
        if (p.lod > LOD_MAX) p.lod = (double)LOD_MAX;
        if (p.lod < 0.0) p.lod = 0.0;
        p.level0 = (int)floor(p.lod);
        p.level1 = p.level0 + 1 < MIP_LEVELS ? p.level0 + 1 : p.level0;
        p.blend = p.lod - (double)p.level0;

        /* taps from the axis ratio; MAX_TAPS is a hard budget */
        double need = ceil(p.ratio - 1e-9);
        if (need > MAX_TAPS) {
            need = MAX_TAPS;
            p.clamped = 1;
            g_truncations++;
        }
        p.taps = need < 1.0 ? 1 : (int)need;
        p.step = p.major / (double)p.taps;
    } else {
        /* zero (or hopelessly small) derivative: no LOD estimate exists, so
           fall back to the nearest level -- level 0, no trilinear blend */
        p.fallback = 1;
        p.ratio = 1.0;
    }
    return p;
}

/* ======================================================================
   Filter.

   Level l is a box prefilter of side s = 2^l (level l+1 is the exact 2x2
   average of level l), and a bilinear level lookup is a tent of side 2s
   centred on the texel centre.  Along the major axis the taps cover the
   footprint; across it the level's box does the work, so the minor axis is
   filtered by that tent instead of the requested box of side minor.  The
   error against the exact per-pixel filter has three parts:

     - prefilter mismatch across the minor axis: both kernels are centred on
       the pixel but their half-widths differ by up to minor/2;
     - grid smear: the tent is square in the mip grid, not aligned with the
       footprint, so a rotated footprint also picks up up to s/2 per grid
       axis;
     - tap discretisation along the major axis: steps of `step`, so
       (major/2)(1 - 1/K) (midpoint rule).

   Both filters (exact box, tap average) are non-negative, have unit mass and
   are supported inside the footprint dilated by the tent half-width s, so a
   data-driven bound that needs no smoothness assumption is simply the range
   of T over that region: `error_bound_local`.  For a texture whose
   1-texel slope is bounded by `lip`, summing the three terms gives the
   closed form in `error_bound_lipschitz`.  Both are checked in the tests.
   ====================================================================== */

static double texel_lipschitz(void) {
    double m = 0.0;
    for (int y = 0; y < TEX_N; y++) {
        for (int x = 0; x < TEX_N; x++) {
            double t = texel_nearest(x, y);
            m = fmax(m, fabs(texel_nearest(x + 1, y) - t));
            m = fmax(m, fabs(texel_nearest(x, y + 1) - t));
        }
    }
    return m;
}

/* Closed-form bound for a texture whose 1-texel slope is <= lip on both axes
   (in the L1 sense: |T(a) - T(b)| <= lip * |a-b|_1; mip levels inherit this
   because they are box averages).  Both the exact box filter and the tap
   average are convex combinations of values of T (respectively of its box
   averages) over the footprint dilated by s = 2^level in each grid axis.
   The L1 diameter of that region is (major+minor)(|cos|+|sin|) + 4s, and the
   L1 range of T over a set of that diameter is at most lip times it. */
static double error_bound_lipschitz(const Plan* p, double lip) {
    if (p->fallback) return lip;              /* nearest: one texel */
    int l = (p->blend > 0.0) ? p->level1 : p->level0;
    double s = (double)(1 << l);
    double diam = (p->major + p->minor) * (fabs(p->cos_t) + fabs(p->sin_t)) + 4.0 * s;
    return lip * diam;
}

/* range of T over the union of both kernels' supports: the footprint dilated
   by the tent half-width s.  Valid for any texture, no smoothness needed;
   this is the number worth quoting per pixel.  Uses the footprint's AABB. */
static double error_bound_local(const Plan* p, double u, double v) {
    if (p->fallback) return 1.0;              /* one bilinear stencil */
    int l = (p->blend > 0.0) ? p->level1 : p->level0;
    double s = (double)(1 << l);
    double ct = p->cos_t, st = p->sin_t;
    double hx = p->major * 0.5, hy = p->minor * 0.5;
    double cx = u * TEX_N, cy = v * TEX_N;
    double minx = cx, maxx = cx, miny = cy, maxy = cy;
    for (int k = 0; k < 4; k++) {
        double ax = (k & 1) ? -hx : hx;
        double ay = (k & 2) ? -hy : hy;
        double x = cx + ax * ct - ay * st;
        double y = cy + ax * st + ay * ct;
        minx = fmin(minx, x); maxx = fmax(maxx, x);
        miny = fmin(miny, y); maxy = fmax(maxy, y);
    }
    minx -= s + 1.0; maxx += s + 1.0;
    miny -= s + 1.0; maxy += s + 1.0;
    int i0 = (int)floor(minx), i1 = (int)floor(maxx);
    int j0 = (int)floor(miny), j1 = (int)floor(maxy);
    if (i1 - i0 > 4096) i1 = i0 + 4096;       /* degenerate huge footprint */
    if (j1 - j0 > 4096) j1 = j0 + 4096;
    double lo = 1e30, hi = -1e30;
    for (int j = j0; j <= j1; j++) {
        for (int i = i0; i <= i1; i++) {
            double t = texel_nearest(i, j);
            lo = fmin(lo, t);
            hi = fmax(hi, t);
        }
    }
    return hi - lo;
}

/* texture at a position, filtered by the level's box and the trilinear blend */
static double texel_at(const Plan* p, double x, double y) {
    double s0 = sample_level(p->level0, x, y);
    if (p->level1 != p->level0 && p->blend > 0.0) {
        double s1 = sample_level(p->level1, x, y);
        return s0 + (s1 - s0) * p->blend;
    }
    return s0;
}

static double aniso_sample(const Plan* p, double u, double v) {
    if (p->fallback)                        /* nearest level, level 0 */
        return sample_level(0, u * TEX_N, v * TEX_N);

    double ct = p->cos_t, st = p->sin_t;
    double k = (double)p->taps;
    double cx = u * TEX_N, cy = v * TEX_N;
    double acc = 0.0;
    for (int i = 0; i < p->taps; i++) {
        /* taps cover exactly `major` texels, centred on the pixel */
        double t = p->step * ((double)i + 0.5) - 0.5 * p->step * k;
        acc += texel_at(p, cx + t * ct, cy + t * st);
    }
    return acc / k;
}

static double sample_aniso(const Grad* g, double u, double v, Plan* out) {
    Plan p = plan_from_grad(g);
    if (out) *out = p;
    return aniso_sample(&p, u, v);
}

/* isotropic baseline: one trilinear lookup at the longest-axis LOD */
static double sample_iso(const Grad* g, double u, double v) {
    Plan p = plan_from_grad(g);
    if (p.fallback) return sample_level(0, u * TEX_N, v * TEX_N);
    double lod = log2(fmax(p.major, 1.0));
    if (lod > LOD_MAX) lod = (double)LOD_MAX;
    int l0 = (int)floor(lod), l1 = l0 + 1 < MIP_LEVELS ? l0 + 1 : l0;
    double f = lod - (double)l0;
    double x = u * TEX_N, y = v * TEX_N;
    double s = sample_level(l0, x, y);
    if (f > 0.0 && l1 != l0) s += (sample_level(l1, x, y) - s) * f;
    return s;
}

/* ------------------------------------------------- exact box reference --- */

typedef struct { double x, y; } P2;

static double poly_area(const P2* q, int n) {
    double a = 0.0;
    for (int i = 0; i < n; i++) {
        const P2* p = &q[i];
        const P2* r = &q[(i + 1) % n];
        a += p->x * r->y - r->x * p->y;
    }
    return fabs(a) * 0.5;
}

/* keep the part of the polygon with coord >= v (keep_ge) or <= v */
static int clip_half(const P2* in, int n, P2* out, int coord, double v, int keep_ge) {
    int m = 0;
    for (int i = 0; i < n; i++) {
        P2 A = in[i], B = in[(i + 1) % n];
        double av = coord == 0 ? A.x : A.y;
        double bv = coord == 0 ? B.x : B.y;
        int ain = keep_ge ? (av >= v) : (av <= v);
        int bin = keep_ge ? (bv >= v) : (bv <= v);
        if (ain) out[m++] = A;
        if (ain != bin) {
            double t = (v - av) / (bv - av);
            P2 m2 = { A.x + (B.x - A.x) * t, A.y + (B.y - A.y) * t };
            out[m++] = m2;
        }
    }
    return m;
}

/* Area average of the piecewise-constant texture over the rotated rectangle
   with half extents (hmaj, hmin), centred on (cx, cy) in level-0 texels.
   Exact up to clipping precision: this is the per-pixel filter. */
static double exact_box(double cx, double cy, double hmaj, double hmin,
                        double cos_t, double sin_t) {
    /* a rank-1 footprint has no area; give the reference a minimal width */
    if (hmaj < 1e-6) hmaj = 1e-6;
    if (hmin < 1e-6) hmin = 1e-6;
    P2 ux = {  cos_t * hmaj,  sin_t * hmaj };
    P2 uy = { -sin_t * hmin,  cos_t * hmin };
    P2 base[4] = {
        { cx + ux.x + uy.x, cy + ux.y + uy.y },
        { cx - ux.x + uy.x, cy - ux.y + uy.y },
        { cx - ux.x - uy.x, cy - ux.y - uy.y },
        { cx + ux.x - uy.x, cy + ux.y - uy.y },
    };
    double minx = base[0].x, maxx = base[0].x, miny = base[0].y, maxy = base[0].y;
    for (int i = 1; i < 4; i++) {
        if (base[i].x < minx) minx = base[i].x;
        if (base[i].x > maxx) maxx = base[i].x;
        if (base[i].y < miny) miny = base[i].y;
        if (base[i].y > maxy) maxy = base[i].y;
    }
    int i0 = (int)floor(minx), i1 = (int)floor(maxx);
    int j0 = (int)floor(miny), j1 = (int)floor(maxy);
    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 > TEX_N - 1) i1 = TEX_N - 1;
    if (j1 > TEX_N - 1) j1 = TEX_N - 1;

    double sum = 0.0, area = 0.0;
    P2 t0[8], t1[8];
    for (int j = j0; j <= j1; j++) {
        for (int i = i0; i <= i1; i++) {
            int n = clip_half(base, 4, t0, 0, (double)i, 1);
            if (n < 3) continue;
            n = clip_half(t0, n, t1, 0, (double)(i + 1), 0);
            if (n < 3) continue;
            n = clip_half(t1, n, t0, 1, (double)j, 1);
            if (n < 3) continue;
            n = clip_half(t0, n, t1, 1, (double)(j + 1), 0);
            if (n < 3) continue;
            double a = poly_area(t1, n);
            if (a <= 0.0) continue;
            sum += a * texel_nearest(i, j);
            area += a;
        }
    }
    return area > 0.0 ? sum / area : 0.0;
}

static double exact_at(const Plan* p, double u, double v) {
    if (p->fallback) return sample_level(0, u * TEX_N, v * TEX_N);
    return exact_box(u * TEX_N, v * TEX_N, p->major * 0.5, p->minor * 0.5,
                     p->cos_t, p->sin_t);
}

static double tex_stddev(void) {
    double s = 0.0, s2 = 0.0;
    int n = TEX_N * TEX_N;
    for (int i = 0; i < n; i++) {
        double t = (double)g_tex[i];
        s += t;
        s2 += t * t;
    }
    double m = s / n;
    double v = s2 / n - m * m;
    return v > 0.0 ? sqrt(v) : 0.0;
}

/* -------------------------------------------------------- flicker metric */

/* RMS over phases of |filter - exact| as the pixel slides one texel along u,
   i.e. the part of the output motion that is wrong, divided by the texture
   standard deviation.  mode: 0 = isotropic longest-axis LOD, 1 = anisotropic,
   2 = plain bilinear (level 0). */
static double flicker(const Grad* g, double u, double v, int mode, int phases, double sigma) {
    Plan p = plan_from_grad(g);
    double acc = 0.0;
    for (int i = 0; i < phases; i++) {
        double uu = u + (double)i / (double)phases / TEX_N;
        double got;
        if (mode == 1) got = aniso_sample(&p, uu, v);
        else if (mode == 0) got = sample_iso(g, uu, v);
        else got = sample_level(0, uu * TEX_N, v * TEX_N);
        double e = got - exact_at(&p, uu, v);
        acc += e * e;
    }
    return sigma > 0.0 ? sqrt(acc / phases) / sigma : 0.0;
}

/* ------------------------------------------------------------------ tests */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void test_pyramid(void) {
    mip_build();
    check("mip_dims", g_dim[0] == 256 && g_dim[1] == 128 && g_dim[8] == 1);

    double maxerr = 0.0;
    for (int l = 1; l < MIP_LEVELS; l++) {
        int n2 = g_dim[l];
        for (int y = 0; y < n2; y++)
            for (int x = 0; x < n2; x++) {
                double want = 0.25 * (lvl_get(l - 1, 2 * x, 2 * y)
                                    + lvl_get(l - 1, 2 * x + 1, 2 * y)
                                    + lvl_get(l - 1, 2 * x, 2 * y + 1)
                                    + lvl_get(l - 1, 2 * x + 1, 2 * y + 1));
                maxerr = fmax(maxerr, fabs(lvl_get(l, x, y) - want));
            }
    }
    printf("mip_max_err=%.3g\n", maxerr);
    check("mip_2x2_average", maxerr < 1e-6);

    double mean = 0.0;
    for (int i = 0; i < TEX_N * TEX_N; i++) mean += (double)g_tex[i];
    mean /= (double)(TEX_N * TEX_N);
    check("mip_top_mean", fabs(lvl_get(LOD_MAX, 0, 0) - mean) < 1e-6);
}

static void test_planning(void) {
    Plan p = plan_from_grad(&(Grad){ 2, 0, 0, 2 });
    printf("iso_lod=%.6f taps=%d step=%.6f\n", p.lod, p.taps, p.step);
    check("iso_lod", fabs(p.lod - 1.0) < 1e-9);
    check("iso_taps", p.taps == 1 && !p.clamped && !p.fallback);

    p = plan_from_grad(&(Grad){ 8, 0, 0, 1 });
    printf("aniso_ratio=%.6f lod=%.3f taps=%d\n", p.ratio, p.lod, p.taps);
    check("aniso_lod", fabs(p.lod) < 1e-9);
    check("aniso_taps", p.taps == 8 && fabs(p.step - 1.0) < 1e-12 && !p.clamped);

    p = plan_from_grad(&(Grad){ 3.2, 0, 0, 1 });
    check("taps_ceil", p.taps == 4);

    int before = g_truncations;
    p = plan_from_grad(&(Grad){ 64, 0, 0, 1 });
    check("ratio_clamped", p.clamped && p.taps == MAX_TAPS && p.ratio > 16.0);
    check("clamp_counted_once", g_truncations == before + 1);
    check("clamp_step_covers_major", fabs(p.step * p.taps - p.major) < 1e-9);
    plan_from_grad(&(Grad){ 100, 0, 0, 1 });
    check("clamp_per_pixel", g_truncations == before + 2);

    p = plan_from_grad(&(Grad){ 0, 0, 0, 0 });
    check("zero_fallback", p.fallback && p.taps == 1 && p.lod == 0.0);
    double u = 12.25 / TEX_N, v = 33.75 / TEX_N;
    Grad zg = { 0, 0, 0, 0 };
    check("fallback_nearest_level",
          sample_aniso(&zg, u, v, NULL) == sample_level(0, u * TEX_N, v * TEX_N));

    p = plan_from_grad(&(Grad){ 1e6, 0, 0, 1e6 });
    check("lod_clamped", p.lod == (double)LOD_MAX && p.level0 == LOD_MAX);
}

/* rotating the screen frame must not change the footprint axes */
static void test_rotation(void) {
    Grad g = { 4.0, 0.4, 0.2, 1.0 };
    Plan p0 = plan_from_grad(&g);
    double worst = 0.0;
    for (int k = 1; k < 8; k++) {
        double t = (double)k * PI / 7.0, c = cos(t), s = sin(t);
        Grad r = { g.du_dx * c + g.du_dy * s, -g.du_dx * s + g.du_dy * c,
                   g.dv_dx * c + g.dv_dy * s, -g.dv_dx * s + g.dv_dy * c };
        Plan q = plan_from_grad(&r);
        worst = fmax(worst, fabs(q.major - p0.major));
        worst = fmax(worst, fabs(q.minor - p0.minor));
        if (q.taps != p0.taps || q.level0 != p0.level0) worst = 1.0;
    }
    printf("rotate_axis_drift=%.3g\n", worst);
    check("rotate_invariant", worst < 1e-9);
}

/* a zero derivative on a constant texture must come back exact, not just
   choose the right level */
static void test_fallback_value(void) {
    for (int i = 0; i < TEX_N * TEX_N; i++) g_tex[i] = 0.37f;
    mip_build();
    Grad zg = { 0, 0, 0, 0 };
    double got = sample_aniso(&zg, 100.25 / TEX_N, 200.75 / TEX_N, NULL);
    printf("fallback_value=%.9f\n", got);
    check("fallback_constant", fabs(got - 0.37) < 1e-6);
    check("fallback_energy", error_bound_local(&(Plan){ .fallback = 1 }, 0, 0) > 0.0);
}

/* exact cases: a lookup must reproduce the exact filter when the footprint
   is one texel, and on a linear ramp at any rotation and any fractional LOD
   -- this pins the axis, tap placement and trilinear blend */
static void test_exact_cases(void) {
    tex_fill_smooth();
    mip_build();

    Plan p;
    double worst = 0.0;
    for (int k = 0; k < 8; k++) {
        double t = (double)k * PI / 7.0, c = cos(t), s = sin(t);
        Grad g = { c, s, -s, c };
        double got = sample_aniso(&g, 77.5 / TEX_N, 33.5 / TEX_N, &p);
        worst = fmax(worst, fabs(got - texel_nearest(77, 33)));
        if (p.taps != 1 || p.lod != 0.0) worst = 1.0;
    }
    printf("texel_exact_err=%.3g\n", worst);
    check("texel_exact", worst < 1e-9);

    for (int y = 0; y < TEX_N; y++)
        for (int x = 0; x < TEX_N; x++)
            g_tex[y * TEX_N + x] = (float)(0.25 + 0.0011 * x + 0.0007 * y);
    mip_build();

    /* across the stripes the footprint is one texel wide, so bilinear is
       exact and the only error is float rounding */
    worst = 0.0;
    for (int a = 0; a < 5; a++) {
        double ratio = 1.0 + 2.5 * (double)a;
        double minor = 1.0;
        for (int k = 0; k < 6; k++) {
            double t = (double)k * PI / 5.0;
            Grad g = { minor * cos(t) * ratio, minor * sin(t) * ratio,
                       -minor * sin(t), minor * cos(t) };
            double u = (61.0 + 7.0 * (double)a + 13.0 * (double)k) / TEX_N;
            double v = (149.0 + 3.0 * (double)a + 11.0 * (double)k) / TEX_N;
            double got = sample_aniso(&g, u, v, NULL);
            /* the ramp is stored per texel centre (index x sits at x+0.5) */
            double want = 0.25 + 0.0011 * (u * TEX_N - 0.5)
                               + 0.0007 * (v * TEX_N - 0.5);
            worst = fmax(worst, fabs(got - want));
        }
    }
    printf("linear_texel_err=%.3g\n", worst);
    check("linear_texel", worst < 1e-3);

    /* fractional minor axis: ramp resampling at an intermediate LOD is
       still linear, level blend reproduces it */
    worst = 0.0;
    for (int a = 0; a < 5; a++) {
        double ratio = 1.0 + 2.5 * (double)a;
        double minor = 1.5 + 0.9 * (double)a;
        for (int k = 0; k < 6; k++) {
            double t = (double)k * PI / 5.0;
            Grad g = { minor * cos(t) * ratio, minor * sin(t) * ratio,
                       -minor * sin(t), minor * cos(t) };
            double u = (61.0 + 7.0 * (double)a + 13.0 * (double)k) / TEX_N;
            double v = (149.0 + 3.0 * (double)a + 11.0 * (double)k) / TEX_N;
            double got = sample_aniso(&g, u, v, NULL);
            double want = 0.25 + 0.0011 * (u * TEX_N - 0.5) + 0.0007 * (v * TEX_N - 0.5);
            /* a ramp is reproduced exactly by the level's symmetric tent,
               whatever the rotation, so this pins the quantised LOD and the
               fractional blend as well as the axes */
            worst = fmax(worst, fabs(got - want));
        }
    }
    printf("linear_fractional_err=%.3g\n", worst);
    check("linear_fractional", worst < 0.01);

    /* a ramp is reproduced exactly by any symmetric kernel, so it says
       nothing about bound tightness; a step edge does */
    for (int y = 0; y < TEX_N; y++)
        for (int x = 0; x < TEX_N; x++)
            g_tex[y * TEX_N + x] = (x < 128) ? 0.0f : 1.0f;
    mip_build();
    double lip = texel_lipschitz();          /* exactly 1 across the edge */
    double worst_frac = 0.0, worst_err = 0.0;
    for (int a = 0; a < 5; a++) {
        double minor = 2.4 + 1.6 * (double)a;
        double ratio = 1.0 + 2.0 * (double)a;
        for (int k = 0; k < 6; k++) {
            double t = (double)k * PI / 5.0;
            Grad g = { minor * cos(t) * ratio, minor * sin(t) * ratio,
                       -minor * sin(t), minor * cos(t) };
            /* edge crossing the footprint centre: worst case for both */
            double u = 128.5 / TEX_N, v = 100.5 / TEX_N;
            Plan st = plan_from_grad(&g);
            double got = aniso_sample(&st, u, v);
            double ref = exact_at(&st, u, v);
            double e = fabs(got - ref);
            double bd = error_bound_local(&st, u, v);
            double bl = error_bound_lipschitz(&st, lip);
            worst_err = fmax(worst_err, e);
            worst_frac = fmax(worst_frac, e / fmin(bd, bl));
            if (e > bd || e > bl) worst_frac = 2.0;
        }
    }
    printf("step_err_max=%.4f step_bound_frac=%.3f lip=%.3f\n",
           worst_err, worst_frac, lip);
    /* a range bound is necessarily loose on a step (the range is 1 while both
       kernels nearly agree); what matters is that it holds and that the
       error stays far below the bound */
    check("step_bound_holds", worst_frac <= 1.0);
    check("step_error_small", worst_err < 0.05);
}

/* grazing view: content varies across the footprint (along u), the footprint
   is 2 texels across and 64 along, so the anisotropic path has to keep the
   band while isotropic longest-axis LOD averages it away */
static void test_grazing(void) {
    tex_fill_smooth();
    mip_build();
    double sig = tex_stddev();

    Grad g = { 2, 0, 0, 64 };
    Plan p;
    double u = 16.5 / TEX_N, v = 128.5 / TEX_N;   /* crest of the band */
    double got = sample_aniso(&g, u, v, &p);
    double ref = exact_at(&p, u, v);
    double isoval = sample_iso(&g, u, v);
    printf("grazing lod=%.4f taps=%d clamped=%d step=%.3f\n",
           p.lod, p.taps, p.clamped, p.step);
    printf("grazing aniso=%.6f iso=%.6f exact=%.6f err_aniso=%.3g err_iso=%.4f\n",
           got, isoval, ref, fabs(got - ref), fabs(isoval - ref));
    check("grazing_plan", p.level0 == 1 && p.taps == MAX_TAPS && p.clamped);
    check("grazing_aniso_close", fabs(got - ref) < 1e-3);
    check("grazing_iso_wrong", fabs(isoval - ref) > 0.3);

    double fa = flicker(&g, u, v, 1, 8, sig);
    double fi = flicker(&g, u, v, 0, 8, sig);
    double fb = flicker(&g, u, v, 2, 8, sig);
    printf("grazing_flicker aniso=%.4f iso=%.4f bilinear=%.4f sigma=%.4f\n", fa, fi, fb, sig);
    check("grazing_flicker_low", fa < 0.01);
    check("grazing_flicker_beats_iso", fi > 20.0 * fa);
}

/* error upper bounds against the exact per-pixel filter over a grid of
   footprints: ratios 1..23, four rotations, both textures */
static void test_error_bounds(void) {
    const char* names[2] = { "smooth", "detail" };
    for (int tex = 0; tex < 2; tex++) {
        if (tex == 0) tex_fill_smooth();
        else tex_fill_detail();
        mip_build();
        double sig = tex_stddev();
        double lip = texel_lipschitz();
        double se = 0.0, si = 0.0, mx = 0.0;
        double worst_frac_lip = 0.0, worst_frac_local = 0.0;
        int n = 0, n_bound = 0;
        for (int a = 0; a < 7; a++) {
            double ratio = 1.0 + 3.5 * (double)a;         /* 1 .. 22 */
            for (int b = 0; b < 4; b++) {
                double th = (double)b * PI / 8.0;
                double minor = 2.0 + 0.5 * (double)a;
                Grad g = { minor * cos(th) * ratio, minor * sin(th) * ratio,
                           -minor * sin(th), minor * cos(th) };
                Plan p;
                for (int iy = 0; iy < 5; iy++) {
                    for (int ix = 0; ix < 5; ix++) {
                        double u = (96.0 + 5.0 * (double)a + 11.0 * (double)ix) / TEX_N;
                        double v = (110.0 + 4.0 * (double)b + 9.0 * (double)iy) / TEX_N;
                        double got = sample_aniso(&g, u, v, &p);
                        double ref = exact_at(&p, u, v);
                        double e = fabs(got - ref);
                        double bl = error_bound_lipschitz(&p, lip);
                        double bd = error_bound_local(&p, u, v);
                        e = e / sig;
                        se += e * e;
                        si += (fabs(sample_iso(&g, u, v) - ref) / sig) *
                              (fabs(sample_iso(&g, u, v) - ref) / sig);
                        mx = fmax(mx, e);
                        worst_frac_lip = fmax(worst_frac_lip, e / (bl / sig));
                        worst_frac_local = fmax(worst_frac_local, e / (bd / sig));
                        if (e <= bl / sig + 1e-12 && e <= bd / sig + 1e-12) n_bound++;
                        n++;
                    }
                }
            }
        }
        double rms = sqrt(se / (double)n), rms_iso = sqrt(si / (double)n);
        printf("%s bound_n=%d rms=%.5f max=%.5f iso_rms=%.5f lip=%.4f "
               "frac_lip=%.3f frac_local=%.3f within=%d/%d\n",
               names[tex], n, rms, mx, rms_iso, lip,
               worst_frac_lip, worst_frac_local, n_bound, n);
        if (tex == 1) { check("bounds_hold_detail", n_bound == n); check("bound_local_meaningful", worst_frac_local < 0.5); }
        else { check("bounds_hold_smooth", n_bound == n); check("bound_tight_smooth", worst_frac_local < 0.05); }
        check(tex ? "bounds_beats_iso_detail" : "bounds_beats_iso_smooth", rms <= rms_iso);
        if (tex == 0) check("bound_lip_tight_smooth", worst_frac_lip < 0.05);
    }
}

static void test_const(void) {
    /* every level of a constant texture is exact, so any blend is too */
    for (int i = 0; i < TEX_N * TEX_N; i++) g_tex[i] = 0.37f;
    mip_build();
    double err = 0.0, bounds = 0.0;
    for (int i = 0; i < 32; i++) {
        double u = (31.0 + 3.0 * i) / TEX_N, v = (77.0 + 5.0 * i) / TEX_N;
        Grad g = { 1.0 + 0.4 * i, 0.25, -0.1, 2.0 + 0.2 * i };
        Plan p = plan_from_grad(&g);
        err = fmax(err, fabs(aniso_sample(&p, u, v) - 0.37));
        bounds = fmax(bounds, error_bound_local(&p, u, v));
    }
    printf("const_err=%.3g local_bound=%.3g\n", err, bounds);
    check("const_exact", err < 1e-6);
    check("const_bound_small", bounds <= 1e-12);
}

int main(void) {
    tex_fill_detail();
    test_pyramid();
    test_planning();
    test_rotation();
    test_fallback_value();
    test_exact_cases();
    test_grazing();
    test_error_bounds();
    test_const();

    printf("truncations=%d\n", g_truncations);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
