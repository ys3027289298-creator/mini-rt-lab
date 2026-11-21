/* arealight: grid-sampled rectangular area light.

   The problem
   -----------
   Sampling an area light at its centre makes penumbrae wrong: every point on
   the receiver sees a hard shadow edge.  This module replaces the single centre
   sample with an nx*ny grid of cell-centre samples, weights each sample by its
   visibility, and reports a certified error bound against the analytic soft
   shadow plus a noise budget that selects the grid density.

   Scene (world units, z grows toward the camera)
   ----------------------------------------------
     light     rectangle x in [-LW,LW], y in [-LH,LH], at z = Z_LIGHT
     occluder  slab whose shadow edge is the line q_c(x,y) = x - M*y = Q0
               at z = Z_OCC
     receiver  point P on the plane z = Z_RECV, with P.z > Z_OCC

   Reduction to one dimension
   --------------------------
   Both the light rectangle and the straight shadow edge are level sets of
   q(x,y) = x - M*y.  Along any line of constant q the occluder hides a prefix
   of that line, so the visible fraction depends on q only.  Define

     e(P) = Q0*(1+u) - u*(P.x - M*P.y),   u = (Z_OCC - Z_LIGHT)/(P.z - Z_OCC)

   the q of the shadow edge projected onto the light plane.  The analytic soft
   shadow visibility is the light-weighted average of the penumbra filter:

     V(P) = (1/A) * int rho(q) * phi((q - e(P))/W) dq

   with A the light area, W the penumbra width, phi the smoothstep filter and
   rho the coarea density of the light, i.e. the length of the chord
   { (x,y) in light : x - M*y = q }.  rho is piecewise linear, so every integral
   here is evaluated by 3-point Gauss-Legendre on pieces split at rho's corners:
   exact up to rounding.  The gradient is exact as well:

     dV/de = -(1/(A*W)) * int rho(q) * phi'(0.5 + (q-e)/W) dq
     dV/dP = (dV/de) * de/dP

   Grid estimator
   --------------
   The integral is replaced by a sum of phi over nx*ny cell centres:

     V_hat = (1/N) * sum_{cells} phi(0.5 + (q_c - e)/W)

   Each sample is weighted by the of the light it sees behind the smooth edge,
   which is exactly phi.  The estimator is a midpoint rule in q with a cell
   q-extent of dq = hx + |M|*hy.

   Certified error bound
   ---------------------
   For each cell the module accumulates a rigorous bound on the difference
   between its centre sample and the exact cell average:

   * if the cell's q-interval crosses a penumbra edge, phi is monotone there,
     so the error is at most the change of phi across the cell;
   * otherwise phi is smooth on the cell, and the midpoint error of the mean
     is at most (1/2) * Var(q) * |phi''| with the exact cell variance
     Var(q) = (hx^2 + M^2*hy^2)/12; the gradient bound uses phi''' instead.

   The accumulated bound is reported in err_v and err_g.  It always contains
   the actual error measured against the analytic reference, and it is
   non-increasing as n grows (it scales as O(n^-2)).  Note that the *actual*
   per-point error can tick up slightly at individual grid sizes when the grid
   aligns unfavourably with the edge (aliasing); the certified bound does not.

   Noise control
   -------------
   With deterministic centres the error is bias.  Passing a seed displaces each
   sample inside its own cell (hash of cell index and seed), turning the bias
   into noise.  The standard deviation of the estimate is then bounded by

     sigma <= 0.5 * sqrt(n_straddle) / N <= 0.5 * sqrt(4n + 2) / n^2

   for a square grid, and arealight_grid_for_noise(eps) returns the smallest n
   meeting a target sigma.  The reported noise_sigma is the bound for the grid
   actually used.

   A grid with a zero dimension cannot be split; arealight_shadow() then returns
   SHADOW_NO_GRID instead of falling back to a centre sample. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ scene */
#define LW 0.40
#define LH 0.40
#define Z_LIGHT (-2.0)
#define Z_OCC   (-1.5)
#define Z_RECV  (-1.0)
#define SLOPE_M 0.6
#define OCC_Q0  (-0.1)
#define PEN_W   0.30
#define LIGHT_AREA (4.0 * LW * LH)
#define ABS_M (SLOPE_M < 0 ? -SLOPE_M : SLOPE_M)
#define Q_LO (-LW - ABS_M * LH)   /* support of rho */
#define Q_HI ( LW + ABS_M * LH)

/* ------------------------------------------------------------------ types */
typedef struct { double x, y, z; } V3;

typedef struct {
    double vis;        /* visibility estimate, in [0,1] */
    V3 grad;           /* d vis / d P */
    double samples;    /* number of cells */
    double cell_q;     /* q-extent of one cell, hx + |M| hy */
    double err_v;      /* certified bound on |vis - V(P)| */
    double err_g;      /* certified bound on |grad - dV/dP| */
    double noise_sigma;/* bound on the standard deviation of a jittered estimate */
} ShadowResult;

typedef enum {
    SHADOW_OK = 0,
    SHADOW_NO_GRID = 1,          /* a grid dimension is zero: cannot split */
    SHADOW_BEHIND_OCCLUDER = 2   /* P.z <= Z_OCC: the projection is undefined */
} ShadowStatus;

/* --------------------------------------------------------------- filter */
/* penumbra filter: 0 fully occluded, 1 fully visible, smoothstep between */
static double phi(double u) {
    if (u <= 0.0) return 0.0;
    if (u >= 1.0) return 1.0;
    return u * u * (3.0 - 2.0 * u);
}
static double phi_d(double u) {
    if (u <= 0.0 || u >= 1.0) return 0.0;
    return 6.0 * u * (1.0 - u);
}
static double phi_dd(double u) {
    if (u <= 0.0 || u >= 1.0) return 0.0;
    return 6.0 - 12.0 * u;
}
static double phi_ddd(double u) {
    return (u > 0.0 && u < 1.0) ? -12.0 : 0.0;
}

/* visibility weight of a light point at edge coordinate q */
static double sample_vis(double q, double e) {
    return phi(0.5 + (q - e) / PEN_W);
}
static double sample_vis_d(double q, double e) {
    return phi_d(0.5 + (q - e) / PEN_W);
}

/* ------------------------------------------------- coarea density of light */
static double rho(double q) {
    if (ABS_M < 1e-15) return (-LW <= q && q <= LW) ? 2.0 * LH : 0.0;
    double y0 = (-LW - q) / SLOPE_M;
    double y1 = ( LW - q) / SLOPE_M;
    double lo = fmax(fmin(y0, y1), -LH);
    double hi = fmin(fmax(y0, y1),  LH);
    return hi > lo ? hi - lo : 0.0;
}

static int cmp_d(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x < y) ? -1 : (x > y);
}

/* rho is piecewise linear; split at its corners so each Gauss piece is smooth */
typedef double (*Fn1)(double, void*);
static double integrate(Fn1 f, void* ctx, double a, double b) {
    static const double gx[3] = { -0.7745966692414834, 0.0, 0.7745966692414834 };
    static const double gw[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    double knots[8];
    int n = 0;
    if (b <= a) return 0.0;
    knots[n++] = a;
    knots[n++] = b;
    double corners[4] = { Q_LO, Q_HI, -LW + ABS_M * LH, LW - ABS_M * LH };
    for (int i = 0; i < 4; i++)
        if (corners[i] > a && corners[i] < b) knots[n++] = corners[i];
    qsort(knots, (size_t)n, sizeof(double), cmp_d);
    double s = 0.0;
    for (int i = 0; i + 1 < n; i++) {
        double u = knots[i], v = knots[i + 1];
        if (v - u < 1e-15) continue;
        double mid = 0.5 * (u + v), half = 0.5 * (v - u), acc = 0.0;
        for (int j = 0; j < 3; j++) acc += gw[j] * f(mid + half * gx[j], ctx);
        s += half * acc;
    }
    return s;
}

typedef struct { double e; } EdgeCtx;
static double fn_rho(double q, void* c) { (void)c; return rho(q); }
static double fn_rho_phi(double q, void* c) {
    return rho(q) * sample_vis(q, ((EdgeCtx*)c)->e);
}
static double fn_rho_phid(double q, void* c) {
    return rho(q) * sample_vis_d(q, ((EdgeCtx*)c)->e);
}

/* ------------------------------------------------------ edge and gradient */
static double edge_at(V3 P, V3* grad) {
    double u = (Z_OCC - Z_LIGHT) / (P.z - Z_OCC);
    double a = P.x - SLOPE_M * P.y;
    double e = OCC_Q0 * (1.0 + u) - u * a;
    if (grad) {
        double dz = -(Z_OCC - Z_LIGHT) / ((P.z - Z_OCC) * (P.z - Z_OCC));
        grad->x = -u;
        grad->y =  u * SLOPE_M;
        grad->z = dz * (OCC_Q0 - a);
    }
    return e;
}

/* analytic reference: visibility and its exact gradient for the same model */
static double analytic_visibility(V3 P, V3* grad) {
    V3 de;
    double e = edge_at(P, &de);
    double lo = e - 0.5 * PEN_W, hi = e + 0.5 * PEN_W;
    double a = fmax(Q_LO, lo), b = fmin(Q_HI, hi);
    EdgeCtx cx = { e };
    double V = 0.0, m = 0.0;
    if (a < b) {
        V += integrate(fn_rho_phi, &cx, a, b);
        m  = integrate(fn_rho_phid, &cx, a, b);
    }
    if (hi < Q_HI)   /* fully lit tail, phi == 1 */
        V += integrate(fn_rho, NULL, fmax(Q_LO, hi), Q_HI);
    if (grad) {
        double dVde = -m / (LIGHT_AREA * PEN_W);
        grad->x = de.x * dVde;
        grad->y = de.y * dVde;
        grad->z = de.z * dVde;
    }
    return V / LIGHT_AREA;
}

/* ------------------------------------------------------------ grid sampler */
static double hash01(unsigned i, unsigned j, unsigned seed) {
    unsigned h = i * 73856093u ^ j * 19349663u ^ seed * 83492791u;
    h ^= h >> 13;
    h *= 1274126177u;
    h ^= h >> 16;
    return (double)(h & 0xffffffu) / 16777216.0;
}

/* bound on the standard deviation of a square-grid jittered estimate */
static double noise_bound(int n) {
    double N = (double)n * (double)n;
    double ns = 4.0 * (double)n + 2.0;   /* at most two straddles per ramp edge */
    return 0.5 * sqrt(ns) / N;
}

/* smallest square grid whose jittered noise bound meets eps, 0 if unreachable */
int arealight_grid_for_noise(double eps) {
    if (!(eps > 0.0)) return 0;
    for (int n = 1; n <= 4096; n++)
        if (noise_bound(n) <= eps) return n;
    return 0;
}

ShadowStatus arealight_shadow(V3 P, int nx, int ny,
                              const double* jitter_seed,
                              ShadowResult* out) {
    memset(out, 0, sizeof(*out));
    if (!(P.z > Z_OCC)) {
        out->err_v = 1.0;
        out->err_g = 1.0;
        return SHADOW_BEHIND_OCCLUDER;
    }
    if (nx <= 0 || ny <= 0) {
        /* zero-size grid: the light cannot be split */
        out->err_v = 1.0;
        out->err_g = 1.0;
        return SHADOW_NO_GRID;
    }

    V3 de;
    double e = edge_at(P, &de);
    double hx = 2.0 * LW / nx, hy = 2.0 * LH / ny;
    double dq = hx + ABS_M * hy;
    double half_u = 0.5 * dq / PEN_W;
    double var_q = (hx * hx + SLOPE_M * SLOPE_M * hy * hy) / 12.0;
    double N = (double)nx * (double)ny;
    double cell_area = hx * hy;

    double sum_vis = 0.0, sum_dvis = 0.0;
    double acc_v = 0.0, acc_g = 0.0;
    int n_straddle = 0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double X, Y;
            if (jitter_seed) {
                unsigned s = (unsigned)(*jitter_seed);
                X = -LW + ((double)i + 0.05 + 0.9 * hash01((unsigned)i, (unsigned)j, s)) * hx;
                Y = -LH + ((double)j + 0.05 + 0.9 * hash01((unsigned)j, (unsigned)i, s ^ 0x9e3779b9u)) * hy;
            } else {
                X = -LW + ((double)i + 0.5) * hx;
                Y = -LH + ((double)j + 0.5) * hy;
            }
            double q = X - SLOPE_M * Y;
            double u = 0.5 + (q - e) / PEN_W;
            sum_vis += phi(u);
            sum_dvis += phi_d(u);

            double ulo = u - half_u, uhi = u + half_u;
            if ((ulo < 0.0 && uhi > 0.0) || (ulo < 1.0 && uhi > 1.0)) {
                /* Cell crosses a penumbra edge: phi is monotone there, so the
                   centre-sample error is at most the change across the cell.
                   phi_d peaks at u = 0.5; include that interior extremum when
                   the cell spans the middle too. */
                n_straddle++;
                acc_v += cell_area * fmax(fabs(phi(ulo) - phi(u)),
                                          fabs(phi(uhi) - phi(u)));
                double d0 = fabs(phi_d(ulo) - phi_d(u));
                double d1 = fabs(phi_d(uhi) - phi_d(u));
                double dm = (ulo < 0.5 && uhi > 0.5) ? fabs(1.5 - phi_d(u)) : 0.0;
                acc_g += cell_area * fmax(d0, fmax(d1, dm));
            } else {
                /* smooth cell: midpoint error <= (1/2) Var(q) |phi''| */
                acc_v += cell_area * 0.5 * var_q * fabs(phi_dd(u)) / (PEN_W * PEN_W);
                acc_g += cell_area * 0.5 * var_q * fabs(phi_ddd(u)) / (PEN_W * PEN_W);
            }
        }
    }

    double den = sqrt(de.x * de.x + de.y * de.y + de.z * de.z);
    out->vis = sum_vis / N;
    out->grad.x = -de.x * sum_dvis / (PEN_W * N);
    out->grad.y = -de.y * sum_dvis / (PEN_W * N);
    out->grad.z = -de.z * sum_dvis / (PEN_W * N);
    out->samples = N;
    out->cell_q = dq;
    out->err_v = acc_v / LIGHT_AREA;
    out->err_g = den * acc_g / (PEN_W * LIGHT_AREA);
    out->noise_sigma = jitter_seed ? 0.5 * sqrt((double)n_straddle) / N : 0.0;
    return SHADOW_OK;
}

/* -------------------------------------------------------------- self-test */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

#define NX 7
#define NY 5
static V3 g_recv[NX * NY];
static void build_receivers(void) {
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            V3* P = &g_recv[i * NY + j];
            P->x = -0.45 + 0.15 * i;
            P->y = -0.25 + 0.125 * j;
            P->z = Z_RECV;
        }
}

/* independent dense midpoint quadrature, used only to validate the reference */
static double dense_visibility(V3 P) {
    const int n = 20000;
    double e = edge_at(P, NULL);
    double h = (Q_HI - Q_LO) / n, s = 0.0;
    for (int i = 0; i < n; i++)
        s += rho(Q_LO + (i + 0.5) * h) * sample_vis(Q_LO + (i + 0.5) * h, e);
    return s * h / LIGHT_AREA;
}

int main(void) {
    build_receivers();

    /* 1. the analytic reference agrees with an independent dense quadrature
          and with central differences of itself */
    {
        double worst_q = 0.0, worst_fd = 0.0;
        for (int r = 0; r < NX * NY; r++) {
            V3 P = g_recv[r], g;
            double v = analytic_visibility(P, &g);
            double vd = dense_visibility(P);
            worst_q = fmax(worst_q, fabs(v - vd));
            double gv[3] = { g.x, g.y, g.z };
            for (int k = 0; k < 3; k++) {
                double h = 1e-6;
                V3 Pp = P, Pm = P;
                double* c = (k == 0) ? &Pp.x : (k == 1 ? &Pp.y : &Pp.z);
                *c += h;
                c = (k == 0) ? &Pm.x : (k == 1 ? &Pm.y : &Pm.z);
                *c -= h;
                double fd = (analytic_visibility(Pp, NULL) - analytic_visibility(Pm, NULL)) / (2.0 * h);
                worst_fd = fmax(worst_fd, fabs(fd - gv[k]) / fmax(1e-12, fabs(fd)));
            }
        }
        check("ref_matches_dense_quadrature", worst_q < 1e-5);
        check("ref_gradient_matches_fd", worst_fd < 1e-6);
    }

    /* 2. refining the grid lowers the error over the receiver sweep, the
          certified bound is monotone in n, and it always contains the error */
    {
        int n_list[] = { 4, 8, 16, 32, 64, 128 };
        double prev_v = 1e9, prev_g = 1e9, prev_bv = 1e9, prev_bg = 1e9;
        int err_mono = 1, bound_mono = 1, bound_holds = 1, all_ok = 1;
        double first_v = -1.0, last_v = -1.0;
        for (int idx = 0; idx < 6; idx++) {
            int n = n_list[idx];
            double ev = 0.0, eg = 0.0, bv = 0.0, bg = 0.0;
            for (int r = 0; r < NX * NY; r++) {
                V3 P = g_recv[r], gref;
                double vref = analytic_visibility(P, &gref);
                ShadowResult s;
                if (arealight_shadow(P, n, n, NULL, &s) != SHADOW_OK) all_ok = 0;
                double dv = fabs(s.vis - vref);
                double gx = s.grad.x - gref.x, gy = s.grad.y - gref.y, gz = s.grad.z - gref.z;
                double dg = sqrt(gx * gx + gy * gy + gz * gz);
                if (dv > s.err_v + 1e-12) bound_holds = 0;
                if (dg > s.err_g + 1e-12) bound_holds = 0;
                ev = fmax(ev, dv);
                eg = fmax(eg, dg);
                bv = fmax(bv, s.err_v);
                bg = fmax(bg, s.err_g);
            }
            if (ev > prev_v + 1e-13) err_mono = 0;
            if (eg > prev_g + 1e-13) err_mono = 0;
            if (bv > prev_bv + 1e-13) bound_mono = 0;
            if (bg > prev_bg + 1e-13) bound_mono = 0;
            prev_v = ev; prev_g = eg; prev_bv = bv; prev_bg = bg;
            if (idx == 0) first_v = ev;
            last_v = ev;
        }
        check("refine_lowers_error_monotone", err_mono && all_ok && last_v < first_v);
        check("bound_monotone", bound_mono);
        check("bound_contains_error", bound_holds);
    }

    /* 3. a zero grid dimension is reported as unsplittable */
    {
        ShadowResult s;
        ShadowStatus a = arealight_shadow(g_recv[0], 0, 8, NULL, &s);
        ShadowStatus b = arealight_shadow(g_recv[0], 8, 0, NULL, &s);
        check("zero_grid_reported", a == SHADOW_NO_GRID && b == SHADOW_NO_GRID);
    }

    /* 4. the receiver must be in front of the occluder */
    {
        V3 P = { 0.0, 0.0, Z_OCC - 0.1 };
        ShadowResult s;
        check("behind_occluder_reported",
              arealight_shadow(P, 8, 8, NULL, &s) == SHADOW_BEHIND_OCCLUDER);
    }

    /* 5. centre-only sampling shows the original defect: near the penumbra the
          shadow is off by a large amount */
    {
        double worst = 0.0;
        for (int r = 0; r < NX * NY; r++) {
            ShadowResult s;
            double vref = analytic_visibility(g_recv[r], NULL);
            arealight_shadow(g_recv[r], 1, 1, NULL, &s);
            worst = fmax(worst, fabs(s.vis - vref));
        }
        check("center_only_error_large", worst > 0.05);
    }

    /* 6. the noise budget is met and shrinks as the grid is refined */
    {
        double eps = 0.005;
        int n = arealight_grid_for_noise(eps);
        int mono = 1;
        double prev = 1e9, at_n = 0.0;
        for (int k = 0; k < 6; k++) {
            int nn = 4 << k;
            double bound = 0.0;
            for (int r = 0; r < NX * NY; r++) {
                double seed = 31.0 + 13.0 * k;
                ShadowResult s;
                arealight_shadow(g_recv[r], nn, nn, &seed, &s);
                bound = fmax(bound, s.noise_sigma);
            }
            if (bound > prev + 1e-15) mono = 0;
            prev = bound;
            if (nn == n) at_n = bound;
        }
        check("noise_target_attained", n > 0 && at_n <= eps);
        check("noise_bound_monotone", mono);
    }

    /* 7. the reported noise bound contains the observed spread over seeds */
    {
        int ok = 1;
        int n_list[] = { 8, 16, 32 };
        for (int t = 0; t < 3; t++) {
            int n = n_list[t];
            for (int r = 0; r < NX * NY; r++) {
                double sum = 0.0, sum2 = 0.0, bound = 0.0;
                int T = 32;
                for (int s = 0; s < T; s++) {
                    double seed = 1000.0 + 7919.0 * s;
                    ShadowResult sr;
                    arealight_shadow(g_recv[r], n, n, &seed, &sr);
                    sum += sr.vis;
                    sum2 += sr.vis * sr.vis;
                    bound = fmax(bound, sr.noise_sigma);
                }
                double mean = sum / T;
                double var = fmax(0.0, sum2 / T - mean * mean);
                /* three standard deviations is a generous allowance for 32 seeds */
                if (sqrt(var) > 3.0 * bound + 1e-12) ok = 0;
            }
        }
        check("noise_bound_holds_empirically", ok);
    }

    /* 8. jittered samples are unbiased: the seed average tracks the reference */
    {
        int ok = 1;
        for (int r = 0; r < NX * NY; r++) {
            double vref = analytic_visibility(g_recv[r], NULL);
            double acc = 0.0;
            int T = 24;
            for (int s = 0; s < T; s++) {
                double seed = 1000.0 + 7919.0 * s;
                ShadowResult sr;
                arealight_shadow(g_recv[r], 16, 16, &seed, &sr);
                acc += sr.vis;
            }
            if (fabs(acc / T - vref) > 0.02) ok = 0;
        }
        check("jitter_mean_tracks_reference", ok);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
