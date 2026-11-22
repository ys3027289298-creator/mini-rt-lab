/* btf: bidirectional texture parameterization by surface height and view angle.

   Why the highlight direction on cloth / velvet comes out wrong: the 2D
   lookups are being fed raw (u, v) that ignore how height and view actually
   stretch the appearance. A streak that should run along the pile is indexed
   with an isotropic lattice, and at the texture seam the first derivative
   jumps, which flips the streak on the far side of the wrap.

   This module replaces that with a separable, monotone warp of the two axes

       u = a(h)   h = normalized surface height  in [0,1]
       v = b(w)   w = normalized view elevation  in [0,1]

       a(h) = 6h^5 - 15h^4 + 10h^3            (smootherstep)
       b(w) = (1 - cos(pi w)) / 2             (sin^2(pi w / 2))

   Both warps are C2, strictly increasing on (0,1), and have zero first
   derivative at both ends, so the map is a bijection [0,1]^2 -> [0,1]^2 and
   C1 across the seam. Everything below is checked against the closed form,
   and the interpolation error is bounded from the actual second and third
   derivatives rather than asserted.

   Rejection rule: a zero-dimensional parameterization (dims == 0, or a range
   that has been pinned to a point so the map degenerates to a line/point) is
   refused instead of quietly producing a 1D lookup.

   Condition rule: the Jacobian det = a'(h) b'(w) vanishes exactly on the
   boundary of the unit square. A condition-number ceiling KAPPA_MAX is
   enforced; sample points whose kappa exceeds it are reported as
   near-singular instead of being used. */

#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI_F 3.14159265358979323846

/* condition-number ceiling for the 2x2 Jacobian (2-norm) */
#define KAPPA_MAX 10.0

/* global Lipschitz constants for the second derivatives, used to make the
   Taylor bounds rigorous:
     max_{[0,1]} |a'''(h)| = 60   (endpoints of 60(6h^2 - 6h + 1))
     max_{[0,1]} |b'''(w)| = pi^3 / 2  (max |sin| = 1)                    */
#define MAX_A3 60.0
#define MAX_B3 (PI_F * PI_F * PI_F * 0.5)

/* ---------- the 1D warps and their derivatives ---------- */

static double a_warp(double h) {
    if (h <= 0.0) return 0.0;
    if (h >= 1.0) return 1.0;
    return h * h * h * (10.0 + h * (-15.0 + 6.0 * h));
}
static double a_d1(double h) {
    if (h < 0.0 || h > 1.0) return 0.0;
    double t = h * (1.0 - h);
    return 30.0 * t * t;
}
static double a_d2(double h) {
    if (h < 0.0 || h > 1.0) return 0.0;
    return 60.0 * h * (1.0 - h) * (1.0 - 2.0 * h);
}

static double b_warp(double w) {
    if (w <= 0.0) return 0.0;
    if (w >= 1.0) return 1.0;
    return 0.5 * (1.0 - cos(PI_F * w));
}
static double b_d1(double w) {
    if (w < 0.0 || w > 1.0) return 0.0;
    return 0.5 * PI_F * sin(PI_F * w);
}
static double b_d2(double w) {
    if (w < 0.0 || w > 1.0) return 0.0;
    return 0.5 * PI_F * PI_F * cos(PI_F * w);
}

/* ---------- inverses (exist by strict monotonicity; bisection) ---------- */

static double inv_warp(double (*f)(double), double u) {
    if (u <= 0.0) return 0.0;
    if (u >= 1.0) return 1.0;
    double lo = 0.0, hi = 1.0;
    for (int i = 0; i < 80; i++) {
        double mid = 0.5 * (lo + hi);
        if (f(mid) < u) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}
static double a_inv(double u) { return inv_warp(a_warp, u); }
static double b_inv(double u) { return inv_warp(b_warp, u); }

/* ---------- the 2D map, Jacobian, condition number ---------- */

typedef struct { double u, v; } V2;

static V2 map_hw(double h, double w) {
    V2 p;
    p.u = a_warp(h);
    p.v = b_warp(w);
    return p;
}

/* J = [[a'(h), 0], [0, b'(w)]].  Returns the signed determinant. */
static double jacobian_hw(double h, double w) {
    return a_d1(h) * b_d1(w);
}

/* 2-norm condition number = max(singular) / min(singular); singular values of
   a diagonal matrix are |a'| and |b'|.  INFINITY when either vanishes. */
static double cond_hw(double h, double w) {
    double sa = fabs(a_d1(h));
    double sb = fabs(b_d1(w));
    double big = sa > sb ? sa : sb;
    double small = sa > sb ? sb : sa;
    if (small <= 0.0) return INFINITY;
    return big / small;
}

/* ---------- parameterization handle + zero-dim rejection ---------- */

typedef struct {
    int dims;          /* only 2 is accepted; 0 is rejected */
    double hmin, hmax; /* height range */
    double wmin, wmax; /* view range  */
} Param;

/* Returns 1 on success. Rejects dims <= 0 explicitly (zero-dimensional), any
   dims other than 2, and any range collapsed to a point (which would make the
   map effectively 1D or 0D). */
static int param_init(Param* p, int dims, double hmin, double hmax, double wmin, double wmax) {
    memset(p, 0, sizeof(*p));
    if (dims <= 0) return 0;
    if (dims != 2) return 0;
    if (!(hmax > hmin) || !(wmax > wmin)) return 0;
    p->dims = 2;
    p->hmin = hmin; p->hmax = hmax;
    p->wmin = wmin; p->wmax = wmax;
    return 1;
}

/* Map a raw (height, view) pair. Returns 0 if the handle was rejected. */
static int param_map(const Param* p, double h, double w, V2* out) {
    if (p->dims != 2) return 0;
    double zh = (h - p->hmin) / (p->hmax - p->hmin);
    double zw = (w - p->wmin) / (p->wmax - p->wmin);
    if (zh < 0.0) zh = 0.0; else if (zh > 1.0) zh = 1.0;
    if (zw < 0.0) zw = 0.0; else if (zw > 1.0) zw = 1.0;
    *out = map_hw(zh, zw);
    return 1;
}

/* ---------- error bound against pointwise sampling ---------- */

/* Pointwise reference sampling: take the exact map at the cell corners, form
   the linear interpolant (what a pointwise-sampled 2D table would reconstruct
   inside a cell), and measure its deviation from the exact map. Because the
   map is separable and each output depends on one axis only, the 2D bilinear
   error equals the 1D error of that axis, so max over axes is the answer.

   The upper bound is the Lagrange/Taylor bound for linear interpolation,
   |f - L| <= (step^2 / 8) * max |f''|, with max|f''| on the cell bounded
   rigorously via the global third-derivative constants:
       max|f''| <= |f''(cell start)| + step * max|f'''|.                  */
static int sample_error(int n, double* max_measured, double* max_bound) {
    if (n < 2 || n > 100000) return 0;
    double step = 1.0 / (double)n;
    double mm = 0.0, mb = 0.0;

    for (int k = 0; k < n; k++) {
        double x0 = (double)k * step;
        double x1 = x0 + step;
        /* height axis */
        double f0 = a_warp(x0), f1 = a_warp(x1);
        for (int q = 1; q <= 3; q++) {
            double t = (double)q / 4.0;
            double x = x0 + t * step;
            double lin = f0 + t * (f1 - f0);
            double e = fabs(a_warp(x) - lin);
            if (e > mm) mm = e;
        }
        double m2 = fabs(a_d2(x0)) + step * MAX_A3;
        double bound = (step * step / 8.0) * m2;
        if (bound > mb) mb = bound;

        /* view axis */
        f0 = b_warp(x0); f1 = b_warp(x1);
        for (int q = 1; q <= 3; q++) {
            double t = (double)q / 4.0;
            double x = x0 + t * step;
            double lin = f0 + t * (f1 - f0);
            double e = fabs(b_warp(x) - lin);
            if (e > mm) mm = e;
        }
        m2 = fabs(b_d2(x0)) + step * MAX_B3;
        bound = (step * step / 8.0) * m2;
        if (bound > mb) mb = bound;
    }
    *max_measured = mm;
    *max_bound = mb;
    return 1;
}

/* ---------- conditioning report ---------- */

typedef struct {
    int n_boundary_singular;   /* samples exactly on the singular set */
    int n_near_singular;       /* interior samples with kappa > KAPPA_MAX */
    int n_admissible;          /* interior samples with kappa <= KAPPA_MAX */
    double kappa_max_admissible;   /* largest kappa among admissible ones */
    double kappa_worst_near;       /* largest finite kappa above the ceiling */
    double h_worst, w_worst;       /* where kappa_worst_near occurs */
} CondReport;

/* Scan an n x n lattice. The boundary of [0,1]^2 is the exact singular set
   (det = 0); interior points are classified by the positive ceiling. */
static void cond_report(int n, CondReport* r) {
    memset(r, 0, sizeof(*r));
    r->kappa_worst_near = 0.0;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            double h = (double)i / (double)n;
            double w = (double)j / (double)n;
            int boundary = (i == 0 || i == n || j == 0 || j == n);
            double k = cond_hw(h, w);
            if (boundary) {
                if (!isfinite(k)) r->n_boundary_singular++;
                continue;
            }
            if (!isfinite(k) || k > KAPPA_MAX) {
                r->n_near_singular++;
                if (isfinite(k) && k > r->kappa_worst_near) {
                    r->kappa_worst_near = k;
                    r->h_worst = h; r->w_worst = w;
                }
            } else {
                r->n_admissible++;
                if (k > r->kappa_max_admissible) r->kappa_max_admissible = k;
            }
        }
    }
}

/* ---------- self checks ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* --- bijection: endpoints and strict monotonicity --- */
    check("a_endpoints", a_warp(0.0) == 0.0 && a_warp(1.0) == 1.0);
    check("b_endpoints", b_warp(0.0) == 0.0 && b_warp(1.0) == 1.0);

    int mono_a = 1, mono_b = 1;
    double pa = -1.0, pb = -1.0;
    for (int i = 0; i <= 2000; i++) {
        double x = (double)i / 2000.0;
        double va = a_warp(x), vb = b_warp(x);
        if (va <= pa && i > 0) mono_a = 0;
        if (vb <= pb && i > 0) mono_b = 0;
        pa = va; pb = vb;
    }
    check("a_strict_increasing", mono_a);
    check("b_strict_increasing", mono_b);

    /* --- derivatives agree with central differences --- */
    double hd = 1e-6;
    double fda = (a_warp(0.3 + hd) - a_warp(0.3 - hd)) / (2.0 * hd);
    double fdb = (b_warp(0.4 + hd) - b_warp(0.4 - hd)) / (2.0 * hd);
    check("a_deriv_matches_fd", fabs(fda - a_d1(0.3)) < 1e-5);
    check("b_deriv_matches_fd", fabs(fdb - b_d1(0.4)) < 1e-5);

    /* --- boundary C1: first derivative vanishes at both ends of both axes --- */
    check("a_boundary_c1", a_d1(0.0) == 0.0 && a_d1(1.0) == 0.0);
    check("b_boundary_c1", fabs(b_d1(0.0)) < 1e-15 && fabs(b_d1(1.0)) < 1e-15);
    /* second derivative stays finite at the seam (C1 needs only the first) */
    check("boundary_second_deriv_finite",
          isfinite(a_d2(0.0)) && isfinite(a_d2(1.0)) &&
          isfinite(b_d2(0.0)) && isfinite(b_d2(1.0)));

    /* --- Jacobian: positive in the interior, zero on the boundary --- */
    check("jacobian_interior_positive", jacobian_hw(0.5, 0.5) > 0.0);
    check("jacobian_boundary_zero",
          fabs(jacobian_hw(0.0, 0.5)) < 1e-15 && fabs(jacobian_hw(1.0, 0.5)) < 1e-15 &&
          fabs(jacobian_hw(0.5, 0.0)) < 1e-15 && fabs(jacobian_hw(0.5, 1.0)) < 1e-15);

    /* --- inverse: a(a_inv(u)) == u over the domain --- */
    double inv_err = 0.0;
    for (int i = 1; i < 200; i++) {
        double u = (double)i / 200.0;
        double e1 = fabs(a_warp(a_inv(u)) - u);
        double e2 = fabs(b_warp(b_inv(u)) - u);
        if (e1 > inv_err) inv_err = e1;
        if (e2 > inv_err) inv_err = e2;
    }
    check("inverse_roundtrip", inv_err < 1e-12);

    /* --- condition number: finite interior value, ceiling respected --- */
    double kc = cond_hw(0.5, 0.5);
    /* a'(1/2) = 30/16 = 1.875 ; b'(1/2) = pi/2 ; kappa = 1.875 / (pi/2) */
    check("cond_center_value", fabs(kc - 1.875 / (PI_F * 0.5)) < 1e-12);
    check("cond_center_below_ceiling", kc <= KAPPA_MAX);

    /* near-singular points are detected and reported */
    check("cond_boundary_infinite", isinf(cond_hw(0.0, 0.5)));
    check("cond_near_singular_detected", cond_hw(1e-3, 0.5) > KAPPA_MAX);

    CondReport rep;
    cond_report(64, &rep);
    check("report_boundary_singular", rep.n_boundary_singular > 0);
    check("report_has_admissible", rep.n_admissible > 0);
    check("report_ceiling_respected", rep.kappa_max_admissible <= KAPPA_MAX);
    check("report_near_singular_found", rep.n_near_singular > 0);
    check("report_worst_near_above_ceiling", rep.kappa_worst_near > KAPPA_MAX);

    /* --- error vs pointwise sampling: measured <= rigorous bound, order 2 --- */
    double m32, b32, m64, b64;
    check("sample_error_runs", sample_error(64, &m64, &b64) == 1);
    check("sample_error_bound_holds", m64 <= b64 + 1e-15);
    sample_error(32, &m32, &b32);
    check("sample_error_second_order",
          m64 > 0.0 && m32 > m64 && (m32 / m64) > 3.5);

    /* --- map matches closed form --- */
    V2 p = map_hw(0.5, 0.5);
    check("map_u_matches", fabs(p.u - 0.5) < 1e-12);
    check("map_v_matches", fabs(p.v - 0.5) < 1e-12);

    /* --- zero-dimensional parameterizations are rejected --- */
    Param par, bad;
    check("accept_two_dimensional", param_init(&par, 2, 0.0, 1.0, 0.0, 1.0) == 1);
    check("reject_dims_zero", param_init(&bad, 0, 0.0, 1.0, 0.0, 1.0) == 0);
    check("reject_dims_one", param_init(&bad, 1, 0.0, 1.0, 0.0, 1.0) == 0);
    check("reject_constant_height", param_init(&bad, 2, 0.5, 0.5, 0.0, 1.0) == 0);
    check("reject_constant_view", param_init(&bad, 2, 0.0, 1.0, 0.3, 0.3) == 0);

    /* a rejected handle refuses to map */
    V2 dummy;
    param_init(&bad, 0, 0.0, 1.0, 0.0, 1.0);
    check("rejected_handle_no_map", param_map(&bad, 0.5, 0.5, &dummy) == 0);
    check("accepted_handle_maps", param_map(&par, 0.25, 0.75, &dummy) == 1 &&
          fabs(dummy.u - a_warp(0.25)) < 1e-12 && fabs(dummy.v - b_warp(0.75)) < 1e-12);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
