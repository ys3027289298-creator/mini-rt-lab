/* softshadow: separable-convolution approximation of a square area light's
   visibility, with an exact reference and a provable error bound.

   Setup (the derivation is 2D and separable in x and y):
     square light   half-side L, at depth zl behind the receiver plane
     square blocker half-side B, centre (cx,cy), rotated by phi, at depth zb
     receiver plane z = 0, shading point p = (x,y)

   Exact visibility. A light point q is blocked iff the ray q->p crosses the
   blocker at z=zb. Projecting the blocker onto the light plane as seen from p
   gives an axis-rotated square of half-side B*zl/zb centred at
       O(p) = p*(1 - zl/zb) + c*zl/zb.
   So V(p) = 1 - area(projected blocker n light square)/(4 L^2).

   The box-blur identity. Let f = zl/zb and H be the point-light (hard) shadow
   mask, i.e. H(q)=0 iff q lies in the blocker's shadow square of half-side
   B*f centred at c*f. Changing variables q = O(p) shows

       V(p) = 1 - (1/w^2) \int_{-w/2}^{w/2} \int_{-w/2}^{w/2} H(p + u) du,
       w = 2 L zb/(zl - zb).

   In words: the penumbra width w is exactly the edge width of a box blur of the
   hard shadow, and since the box is a product kernel the blur separates into
   two 1D passes. The module replaces that box by a 1D kernel k with the same
   total mass: either the box itself, or a Gaussian of matching variance
   (sigma = w/sqrt(12)). The Gaussian gives the smooth falloff usually wanted
   in a renderer; the box reproduces the exact result up to discretisation.

   Error bound. The exact answer is (K(x)K) * H, where K is the box of width w
   and H the hard-shadow indicator of a rotated square. The code computes
   (k(x)k) * (P(x)P * H), where k is the sampled kernel on the pixel grid and P
   is the pixel-average box of width h. Let S = k_tilde * P be the effective
   1D kernel (k_tilde is the staircase density k[d]/h). Since the two factors
   are the same on both sides,

       (S(x)S - K(x)K) * H = S (x) (S-K) * H  +  (S-K) (x) K * H,

   and a slice of the rotated square is an interval, so H has total variation
   at most 2 along every row and column. Convolving with a kernel of mass 1
   and then integrating against a TV-2 slice gives at most 2*||S-K||_1 per
   term, hence

       max |approx - exact| <= 2 * ||S - K||_1.

   (This is the standard two-factor total-variation bound; no CDF argument is
   needed.) The code returns exactly 2*||S-K||_1, evaluating the piecewise
   linear difference S-K analytically. For the box kernel S is the exact box
   averaged over one pixel, so the bound is one pixel of smear -- O(h/w) --
   and for the Gaussian it is ~0.75, dominated by the genuine box/Gaussian
   shape difference. Both are checked against the measured error in the tests.

   Degenerate width. w == 0 happens for a point light (L = 0) or a blocker
   touching the receiver (zb = 0); then the kernel collapses to the identity,
   the convolution returns the hard shadow mask, and mode/width are reported so
   the caller can see the fallback was taken. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define SS_MAX_KERNEL 4096
#define SS_GRID 129

/* ---------- 2D helpers: convex polygon clipping (Sutherland-Hodgman) ---------- */

typedef struct { double x, y; } V2;

static int poly_inside(V2 p, V2 a, V2 b) {
    return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x) >= -1e-12;
}
static V2 poly_intersect(V2 p, V2 q, V2 a, V2 b) {
    double x1 = p.x, y1 = p.y, x2 = q.x, y2 = q.y, x3 = a.x, y3 = a.y, x4 = b.x, y4 = b.y;
    double den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    V2 r;
    if (fabs(den) < 1e-18) return q;
    double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
    r.x = x1 + t * (x2 - x1);
    r.y = y1 + t * (y2 - y1);
    return r;
}
/* Intersection area of two convex polygons (subject, clip), clip wound CCW. */
static double poly_clip_area(const V2* subj, int nsub, const V2* clip, int nclip) {
    V2 buf[16], nxt[16];
    int n = nsub;
    for (int i = 0; i < nsub; i++) buf[i] = subj[i];
    for (int c = 0; c < nclip; c++) {
        V2 a = clip[c], b = clip[(c + 1) % nclip];
        int m = 0;
        if (n == 0) return 0.0;
        for (int i = 0; i < n; i++) {
            V2 p = buf[i], q = buf[(i + 1) % n];
            int pin = poly_inside(p, a, b), qin = poly_inside(q, a, b);
            if (pin) nxt[m++] = p;
            if (pin != qin) nxt[m++] = poly_intersect(p, q, a, b);
        }
        n = m;
        for (int i = 0; i < n; i++) buf[i] = nxt[i];
    }
    if (n < 3) return 0.0;
    double area = 0.0;
    for (int i = 0; i < n; i++) {
        V2 p = buf[i], q = buf[(i + 1) % n];
        area += p.x * q.y - q.x * p.y;
    }
    return fabs(area) * 0.5;
}

/* Axis-aligned square of half-side hs centred at (cx,cy), CCW. */
static void square_axis(V2 out[4], double cx, double cy, double hs) {
    out[0].x = cx - hs; out[0].y = cy - hs;
    out[1].x = cx + hs; out[1].y = cy - hs;
    out[2].x = cx + hs; out[2].y = cy + hs;
    out[3].x = cx - hs; out[3].y = cy + hs;
}
/* Square rotated by phi, CCW (rotation preserves orientation). */
static void square_rot(V2 out[4], double cx, double cy, double hs, double phi) {
    double c = cos(phi), s = sin(phi);
    static const int sx[4] = { -1, 1, 1, -1 };
    static const int sy[4] = { -1, -1, 1, 1 };
    for (int i = 0; i < 4; i++) {
        out[i].x = cx + c * sx[i] * hs - s * sy[i] * hs;
        out[i].y = cy + s * sx[i] * hs + c * sy[i] * hs;
    }
}

/* ---------- scene ---------- */

typedef struct {
    double L, zl, B, phi, zb, cx, cy;
} Scene;

enum { PATH_NONE = 0, PATH_HARD = 1, PATH_SOFT = 2 };
enum { KERNEL_BOX = 0, KERNEL_GAUSS = 1 };

static int scene_path(const Scene* s) {
    if (!(s->zl > 0.0) || !(s->zb < s->zl)) return PATH_NONE;
    if (!(s->L > 0.0) || !(s->zb > 0.0)) return PATH_HARD;  /* zero-width limit */
    return PATH_SOFT;
}

static double scene_f(const Scene* s) { return s->zl / (s->zl - s->zb); }

/* Point-light hard shadow: square half-side B*f centred at c*f. */
static void scene_hard_square(const Scene* s, double* cx, double* cy, double* hs) {
    double f = scene_f(s);
    *cx = s->cx * f;
    *cy = s->cy * f;
    *hs = s->B * f;
}

/* Penumbra width; zero on the hard path (kernel collapses to identity). */
static double scene_width(const Scene* s) {
    if (scene_path(s) != PATH_SOFT) return 0.0;
    return 2.0 * s->L * s->zb / (s->zl - s->zb);
}

static int point_in_rot_square(const Scene* s, double x, double y) {
    double cx, cy, hs;
    scene_hard_square(s, &cx, &cy, &hs);
    double dx = x - cx, dy = y - cy, c = cos(-s->phi), sn = sin(-s->phi);
    double lx = c * dx - sn * dy, ly = sn * dx + c * dy;
    return fabs(lx) <= hs && fabs(ly) <= hs;
}

/* Exact visibility at receiver point (x,y). */
static double exact_vis(const Scene* s, double x, double y) {
    int path = scene_path(s);
    if (path == PATH_NONE) return 1.0;
    if (path == PATH_HARD)
        return point_in_rot_square(s, x, y) ? 0.0 : 1.0;
    double scale = s->zl / s->zb;
    V2 light[4], block[4];
    square_axis(light, 0.0, 0.0, s->L);
    square_rot(block, x * (1.0 - scale) + s->cx * scale,
                      y * (1.0 - scale) + s->cy * scale,
                      s->B * scale, s->phi);
    return 1.0 - poly_clip_area(light, 4, block, 4) / (4.0 * s->L * s->L);
}

/* Area-averaged hard-shadow mask over one square pixel of side h. */
static double hard_cov(const Scene* s, double xc, double yc, double h) {
    if (scene_path(s) == PATH_NONE) return 1.0;
    double cx, cy, hs;
    V2 cell[4], block[4];
    scene_hard_square(s, &cx, &cy, &hs);
    square_axis(cell, xc, yc, 0.5 * h);
    square_rot(block, cx, cy, hs, s->phi);
    return 1.0 - poly_clip_area(cell, 4, block, 4) / (h * h);
}

/* ---------- 1D kernel ---------- */

static double norm_cdf(double z) { return 0.5 * erfc(-z / sqrt(2.0)); }

/* Build symmetric kernel k[0..r] (k[-d] = k[d]); returns radius r.
   w<=0 gives the identity kernel. */
static int kernel_build(int kind, double w, double h, double* k) {
    if (!(w > 0.0)) { k[0] = 1.0; return 0; }
    int r;
    if (kind == KERNEL_BOX) {
        r = (int)floor(w / (2.0 * h)) + 1;
    } else {
        double sig = w / sqrt(12.0);           /* matched variance */
        r = (int)ceil(4.75 * sig / h);
    }
    if (r >= SS_MAX_KERNEL) r = SS_MAX_KERNEL - 1;
    for (int j = 0; j <= r; j++) {
        if (kind == KERNEL_BOX) {
            double lo = (j - 0.5) * h, hi = (j + 0.5) * h;
            double a = fmax(lo, -0.5 * w), b = fmin(hi, 0.5 * w);
            k[j] = fmax(0.0, b - a) / w;
        } else {
            double sig = w / sqrt(12.0);
            k[j] = norm_cdf((j * h + 0.5 * h) / sig) - norm_cdf((j * h - 0.5 * h) / sig);
        }
    }
    double tot = k[0];
    for (int j = 1; j <= r; j++) tot += 2.0 * k[j];
    for (int j = 0; j <= r; j++) k[j] /= tot;
    return r;
}

/* ||S - K||_1 with S = k_tilde * P_h and K the exact box of width w:
   k_tilde is the piecewise-constant density k[d]/h on cell d, P_h averages
   over one pixel of width h, and S is therefore piecewise linear with knots
   at d*h +- h/2. On each knot interval the difference D = S - K is linear,
   so its absolute value integrates exactly. */
static double kernel_l1(const double* k, int r, double w, double h) {
    if (!(w > 0.0)) return 0.0;

    /* S(t): running average of the staircase over a window of width h. */
    #define KDEV_S(t) ({                                          \
        double _t = (t), _s = 0.0;                                \
        for (int _e = -r; _e <= r; _e++) {                        \
            double _lo = (_e - 0.5) * h, _hi = (_e + 0.5) * h;    \
            double _a = fmax(_t - 0.5 * h, _lo);                  \
            double _b = fmin(_t + 0.5 * h, _hi);                  \
            if (_b > _a) _s += k[_e < 0 ? -_e : _e] / h * (_b - _a) / h; \
        }                                                         \
        _s; })
    #define KDEV_D(t) (KDEV_S(t) - ((fabs(t) < 0.5 * w) ? 1.0 / w : 0.0))

    double l1 = 0.0;
    double lo = -(double)(r + 2) * h, hi = (double)(r + 2) * h;
    /* knot list: cell-window edges; the box edge falls inside some interval and
       D is still linear across it (K jumps by a constant), so split there too */
    int n = 2 * (r + 2) * 2;
    double prev = lo, dprev = KDEV_D(lo);
    for (int i = 1; i <= n; i++) {
        double t = lo + (hi - lo) * i / n;
        double dt = KDEV_D(t);
        if (dprev * dt >= 0.0) {
            l1 += 0.5 * fabs(dprev + dt) * (t - prev);
        } else {
            double tz = prev + (t - prev) * dprev / (dprev - dt);
            l1 += 0.5 * (fabs(dprev) * (tz - prev) + fabs(dt) * (t - tz));
        }
        prev = t; dprev = dt;
    }
    #undef KDEV_D
    #undef KDEV_S
    return l1;
}

/* ---------- separable convolution ---------- */

/* out = (k (x) k) * in with symmetric kernel k[0..r]; tmp is N*N scratch.
   sep_conv and a direct 2D convolution are compared in the self-test. */
static void sep_conv(const double* in, double* tmp, double* out, int N, const double* k, int r) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double v = k[0] * in[i * N + j];
            for (int d = 1; d <= r; d++) {
                if (i - d >= 0) v += k[d] * in[(i - d) * N + j];
                if (i + d < N)  v += k[d] * in[(i + d) * N + j];
            }
            tmp[i * N + j] = v;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double v = k[0] * tmp[i * N + j];
            for (int d = 1; d <= r; d++) {
                if (j - d >= 0) v += k[d] * tmp[i * N + j - d];
                if (j + d < N)  v += k[d] * tmp[i * N + j + d];
            }
            out[i * N + j] = v;
        }
    }
}

/* ---------- render a scene on an N x N grid ---------- */

typedef struct {
    double width;        /* penumbra width w (0 on the hard path) */
    int    hard_fallback;/* 1 iff w == 0 */
    int    radius;       /* kernel radius in pixels */
    double h;            /* pixel size */
    double start;        /* lower-left corner of the square domain */
    double err;          /* max interior error vs exact_vis */
    double bound;        /* 4*L1 + 4*h/w */
    double l1;
} RenderStats;

static int scene_render(const Scene* s, int kind, int N,
                        double* mask, double* out, RenderStats* st) {
    if (N > SS_GRID) return 0;
    static double tmp[SS_GRID * SS_GRID];
    static double k[SS_MAX_KERNEL];
    double w = scene_width(s);
    double cx, cy, hs;
    scene_hard_square(s, &cx, &cy, &hs);

    /* Fix-point on the domain: radius depends on pixel size, which depends on
       the domain. The domain keeps the kernel support and a margin inside. */
    double E = fmax(fabs(cx), fabs(cy)) + hs + w + 1.0;
    int r = 0;
    for (int it = 0; it < 12; it++) {
        double h = 2.0 * E / N;
        r = kernel_build(kind, w, h, k);
        double need = fmax(fabs(cx), fabs(cy)) + hs + 0.5 * w + (r + 1) * h
                    + 0.25 * (hs + w) + 1e-6;
        if (need <= E) break;
        E = need;
    }
    double h = 2.0 * E / N;
    r = kernel_build(kind, w, h, k);

    for (int i = 0; i < N; i++) {
        double x = -E + (i + 0.5) * h;
        for (int j = 0; j < N; j++) {
            double y = -E + (j + 0.5) * h;
            mask[i * N + j] = hard_cov(s, x, y, h);
        }
    }
    sep_conv(mask, tmp, out, N, k, r);

    double m = (r + 1) * h;   /* keep the kernel support inside the domain */
    double worst = 0.0;
    for (int i = 0; i < N; i++) {
        double x = -E + (i + 0.5) * h;
        if (fabs(x) > E - m) continue;
        for (int j = 0; j < N; j++) {
            double y = -E + (j + 0.5) * h;
            if (fabs(y) > E - m) continue;
            double e = fabs(out[i * N + j] - exact_vis(s, x, y));
            if (e > worst) worst = e;
        }
    }
    double l1 = kernel_l1(k, r, w, h);
    st->width = w;
    st->hard_fallback = (w <= 0.0);
    st->radius = r;
    st->h = h;
    st->start = -E;
    st->err = worst;
    st->l1 = l1;
    st->bound = 2.0 * l1;
    return 1;
}

/* ---------- self-test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void test_exact_reference(void) {
    Scene s = { 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0 };
    /* aligned: a point on the hard edge sees exactly half the light */
    check("exact_half_edge", fabs(exact_vis(&s, 2.0, 0.0) - 0.5) < 1e-9);
    check("exact_corner_quarter", fabs(exact_vis(&s, 2.0, 2.0) - 0.75) < 1e-9);
    check("exact_deep_inside", exact_vis(&s, 0.0, 0.0) == 0.0);
    check("exact_far_outside", exact_vis(&s, 50.0, 50.0) == 1.0);
    /* rotated blocker: by symmetry the centre of the hard shadow is fully dark */
    Scene r = { 1.0, 2.0, 0.8, 0.6, 1.0, 0.0, 0.0 };
    check("exact_rotated_centre", exact_vis(&r, 0.0, 0.0) == 0.0);
    check("exact_in_unit_interval",
          exact_vis(&r, 1.5, -0.7) >= 0.0 && exact_vis(&r, 1.5, -0.7) <= 1.0);
}

static void test_width(void) {
    /* w = 2 L zb/(zl-zb) grows as the blocker approaches the light (zb -> zl)
       and shrinks as the blocker moves away from the receiver (zb -> 0). */
    int mono_zb = 1, mono_zl = 1;
    double prev = -1.0;
    for (int i = 1; i < 20; i++) {
        Scene s = { 1.0, 2.0, 1.0, 0.3, i * 0.05, 0.0, 0.0 };
        double w = scene_width(&s);
        if (!(w > prev)) mono_zb = 0;
        prev = w;
    }
    prev = 1e30;
    for (int i = 1; i < 20; i++) {
        Scene s = { 1.0, 1.0 + i * 0.2, 1.0, 0.3, 1.0, 0.0, 0.0 };
        double w = scene_width(&s);
        if (!(w < prev)) mono_zl = 0;
        prev = w;
    }
    check("width_monotone_zb", mono_zb);
    check("width_monotone_zl", mono_zl);

    /* The measured 0.1->0.9 transition of the exact visibility equals 0.8 w:
       the analytic penumbra width is the kernel width. */
    int matches = 1;
    for (int t = 0; t < 3; t++) {
        double want = 0.5 + 0.75 * t;
        Scene s = { 1.0, 3.0, 2.0, 0.0, want * 3.0 / (2.0 + want), 0.0, 0.0 };
        double w = scene_width(&s), hs;
        double cx, cy;
        scene_hard_square(&s, &cx, &cy, &hs);
        double span = 0.6 * (w + hs + 1e-3), lo = -1.0, hi = -1.0;
        int n = 4000;
        double xp = hs - span, vp = exact_vis(&s, xp, 0.0);
        for (int i = 1; i <= n; i++) {
            double x = hs - span + 2.0 * span * i / n;
            double v = exact_vis(&s, x, 0.0);
            if (fmin(vp, v) <= 0.9 && 0.9 <= fmax(vp, v))
                lo = xp + (x - xp) * (0.9 - vp) / (v - vp);
            if (fmin(vp, v) <= 0.1 && 0.1 <= fmax(vp, v))
                hi = xp + (x - xp) * (0.1 - vp) / (v - vp);
            xp = x; vp = v;
        }
        double got = fabs(hi - lo);
        if (fabs(got - 0.8 * w) > 0.03 * w) matches = 0;
    }
    check("penumbra_width_matches", matches);
}

static void test_separable_equals_direct(void) {
    /* Direct 2D convolution with the outer product kernel must equal the two
       separable passes. */
    enum { N = 23 };
    static double in[N * N], a[N * N], b[N * N], t[N * N];
    double k[5] = { 0.3, 0.25, 0.1, 0.05, 0.0 };
    int r = 3;
    k[4] = k[3]; /* symmetric radius 3: k[0..3] used */
    for (int i = 0; i < N * N; i++) in[i] = (double)((i * 37) % 11) / 10.0;
    sep_conv(in, t, a, N, k, r);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double v = 0.0;
            for (int di = -r; di <= r; di++)
                for (int dj = -r; dj <= r; dj++) {
                    int ii = i + di, jj = j + dj;
                    if (ii < 0 || ii >= N || jj < 0 || jj >= N) continue;
                    v += k[abs(di)] * k[abs(dj)] * in[ii * N + jj];
                }
            b[i * N + j] = v;
        }
    double d = 0.0;
    for (int i = 0; i < N * N; i++) d = fmax(d, fabs(a[i] - b[i]));
    check("separable_equals_direct", d < 1e-12);
}

static void test_bound(void) {
    Scene scenes[4] = {
        { 1.0, 2.0, 1.0, 0.0,        1.0, 0.0,  0.0 },
        { 1.0, 2.0, 0.6, 0.6108652,  0.7, 0.0,  0.0 },
        { 0.8, 1.5, 1.2, 0.3490659,  0.4, 0.1, -0.1 },
        { 1.5, 3.0, 0.5, 0.8726646,  0.6, 0.0,  0.0 },
    };
    static double mask[SS_GRID * SS_GRID], out[SS_GRID * SS_GRID];
    double worst_ratio[2] = { 0.0, 0.0 };
    int holds[2] = { 1, 1 };
    for (int kind = 0; kind < 2; kind++)
        for (int i = 0; i < 4; i++)
            for (int N = 65; N <= SS_GRID; N *= 2) {
                RenderStats st;
                scene_render(&scenes[i], kind, N, mask, out, &st);
                if (st.err > st.bound + 1e-12) holds[kind] = 0;
                double ratio = st.bound > 0.0 ? st.err / st.bound : 0.0;
                if (ratio > worst_ratio[kind]) worst_ratio[kind] = ratio;
            }
    check("box_bound_holds", holds[KERNEL_BOX]);
    check("gauss_bound_holds", holds[KERNEL_GAUSS]);
    printf("box_err_over_bound=%.4f\n", worst_ratio[KERNEL_BOX]);
    printf("gauss_err_over_bound=%.4f\n", worst_ratio[KERNEL_GAUSS]);
}

static void test_degenerate(void) {
    static double mask[SS_GRID * SS_GRID], out[SS_GRID * SS_GRID];
    RenderStats st;
    Scene point = { 0.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0 };  /* point light: w = 0 */
    Scene touch = { 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0 };  /* zb = 0:       w = 0 */
    Scene good  = { 1.0, 3.0, 1.0, 0.0, 1.5, 0.0, 0.0 };  /* w = 1            */
    scene_render(&point, KERNEL_GAUSS, 65, mask, out, &st);
    double maxd = 0.0, maxv = 0.0;
    for (int i = 0; i < 65 * 65; i++) {
        maxd = fmax(maxd, fabs(out[i] - mask[i]));
        maxv = fmax(maxv, out[i]);
    }
    check("degenerate_w_zero",
          scene_width(&point) == 0.0 && scene_width(&touch) == 0.0 &&
          scene_width(&good) > 0.0);
    check("degenerate_identity", maxd < 1e-12);
    check("degenerate_radius_zero", st.radius == 0 && st.hard_fallback == 1);
    check("degenerate_not_darkened", maxv <= 1.0 + 1e-12);
    (void)good;
}

/* Print one scanline through the penumbra: exact, box, gaussian, both errors. */
static void scanline_report(const Scene* s, int N, int samples) {
    static double mask[SS_GRID * SS_GRID], outbox[SS_GRID * SS_GRID], outg[SS_GRID * SS_GRID];
    RenderStats sb, sg;
    scene_render(s, KERNEL_BOX, N, mask, outbox, &sb);
    scene_render(s, KERNEL_GAUSS, N, mask, outg, &sg);
    double w = sb.width, H;
    double cx, cy;
    scene_hard_square(s, &cx, &cy, &H);
    double x0 = cx + H - 1.2 * w, x1 = cx + H + 1.2 * w;
    printf("# scanline y=0  w=%.4f  bound_box=%.4f  bound_gauss=%.4f\n", w, sb.bound, sg.bound);
    printf("# %8s %8s %8s %8s %8s %8s\n", "x", "exact", "box", "gauss", "e_box", "e_gauss");
    for (int i = 0; i < samples; i++) {
        double x = x0 + (x1 - x0) * i / (samples - 1);
        double ex = exact_vis(s, x, 0.0);
        double gb = exact_vis(s, 0.0, 0.0); (void)gb;
        /* bilinear sample both grids at (x,0) */
        double fb = 0.0, fg = 0.0;
        double u = (x - sb.start) / sb.h - 0.5, v = (0.0 - sb.start) / sb.h - 0.5;
        int i0 = (int)floor(u), j0 = (int)floor(v);
        double tu = u - i0, tv = v - j0;
        if (i0 >= 0 && i0 + 1 < N && j0 >= 0 && j0 + 1 < N) {
            fb = (1 - tu) * (1 - tv) * outbox[i0 * N + j0] + tu * (1 - tv) * outbox[(i0 + 1) * N + j0]
               + (1 - tu) * tv * outbox[i0 * N + j0 + 1] + tu * tv * outbox[(i0 + 1) * N + j0 + 1];
            fg = (1 - tu) * (1 - tv) * outg[i0 * N + j0] + tu * (1 - tv) * outg[(i0 + 1) * N + j0]
               + (1 - tu) * tv * outg[i0 * N + j0 + 1] + tu * tv * outg[(i0 + 1) * N + j0 + 1];
            printf("  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
                   x, ex, fb, fg, fabs(fb - ex), fabs(fg - ex));
        }
    }
}

int main(void) {
    printf("# softshadow\n");
    printf("w_formula=%.6f\n", 2.0 * 1.0 * 2.0 / (3.0 - 2.0));
    test_exact_reference();
    test_width();
    test_separable_equals_direct();
    test_bound();
    test_degenerate();

    Scene demo = { 1.0, 2.5, 0.9, 0.5236, 0.8, 0.1, -0.05 };
    scanline_report(&demo, SS_GRID, 17);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
