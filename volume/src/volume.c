/* volume: participating-medium transmittance.
 *
 * Optical depth along a ray segment is tau = \int_{t0}^{t1} sigma(t) dt and
 * transmittance is T = exp(-tau).  The exponential is applied once to the
 * integrated density instead of multiplying per-step transmittances (or using
 * 1 - sigma*dt), so a coarse step cannot accumulate a far-field bias.
 *
 * Two integrators:
 *   fixed    - composite Simpson on n, 2n and 4n intervals.  The successive
 *              differences d1 = |S(2n) - S(n)| and d2 = |S(4n) - S(2n)| fall
 *              by 16 in Simpson's asymptotic regime, and d2 is 15x the error
 *              of the returned S(4n): a conservative a-posteriori bound.
 *              `converged` reports whether that 16-fold ratio was observed.
 *   adaptive - recursive Simpson driven by |S_left + S_right - S| / 15, with a
 *              forced minimum refinement so a narrow feature cannot hide
 *              between the three initial samples.
 *
 * Because |exp(-a) - exp(-b)| <= |a - b| for a, b >= 0, the tau error bound
 * carries straight over to T: err_T == err_tau.
 *
 * A medium whose known upper bound is <= 0 (empty) returns T = 1 without a
 * single density evaluation.  We do not guess emptiness from a few probes:
 * the same aliasing that makes a cheap probe unsafe is why adaptive Simpson
 * has min_depth.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define VOL_PI 3.14159265358979323846

typedef struct {
    double (*sigma)(double t, const void *ctx);
    const void *ctx;
    double max_sigma;   /* known upper bound on sigma; <= 0 means empty */
} Medium;

typedef struct {
    double tau;         /* estimated optical depth */
    double T;           /* exp(-tau) */
    double err_tau;     /* a-posteriori absolute error estimate on tau */
    double err_T;       /* |T - T_true| <= err_T (Lipschitz in tau) */
    long evals;         /* density evaluations performed */
    int zero;           /* density vanished at every sample / empty shortcut */
    int converged;      /* error estimate sits in the asymptotic regime */
} Transmittance;

typedef struct {
    double tol_tau;     /* target absolute error on optical depth */
    int min_depth;      /* forced uniform refinement: 2^depth cells */
    int max_depth;
} AdaptiveOpts;

static AdaptiveOpts adaptive_opts(double tol_tau, int min_depth, int max_depth) {
    AdaptiveOpts o;
    o.tol_tau = tol_tau;
    o.min_depth = min_depth;
    o.max_depth = max_depth;
    return o;
}

typedef struct { long evals; int positive; } Acc;

static double sample_sigma(const Medium *m, double t, Acc *acc) {
    double s = m->sigma(t, m->ctx);
    acc->evals++;
    if (s > 0.0) { acc->positive = 1; return s; }
    return 0.0;   /* densities are non-negative; clamp roundoff / NaN to zero */
}

static Transmittance empty_result(void) {
    Transmittance r;
    memset(&r, 0, sizeof(r));
    r.T = 1.0;
    r.zero = 1;
    r.converged = 1;
    return r;
}

/* Composite Simpson on nested grids, one sampling pass.
 * n is even; evaluate sigma on a 4n grid with step h and derive S(n), S(2n),
 * S(4n) from the same values: node i sits at t0 + i*h, and the coarser grids
 * use i divisible by 4 and 2 respectively. */
static void simpson_triplet(const Medium *m, double t0, double t1, int n,
                            double *sn, double *s2n, double *s4n, Acc *acc) {
    int N = 4 * n;
    double h = (t1 - t0) / (double)N;
    double f0 = sample_sigma(m, t0, acc);
    double fN = sample_sigma(m, t1, acc);
    double odd4 = 0.0, even4 = 0.0;
    double odd2 = 0.0, even2 = 0.0;
    double oddn = 0.0, evenn = 0.0;

    for (int i = 1; i < N; i++) {
        double f = sample_sigma(m, t0 + (double)i * h, acc);
        int r = i & 7;
        if (i & 1) {
            odd4 += f;
        } else {
            even4 += f;
            if (r == 2 || r == 6) odd2 += f;
            else                  even2 += f;
            if (r == 4) oddn += f;
            if (r == 0) evenn += f;
        }
    }
    double ends = f0 + fN;
    *s4n = h / 3.0 * (ends + 4.0 * odd4 + 2.0 * even4);
    *s2n = 2.0 * h / 3.0 * (ends + 4.0 * odd2 + 2.0 * even2);
    *sn  = 4.0 * h / 3.0 * (ends + 4.0 * oddn + 2.0 * evenn);
}

Transmittance transmittance_fixed(const Medium *m, double t0, double t1, int n) {
    if (!(t1 > t0) || !(m->max_sigma > 0.0)) return empty_result();

    if (n < 2) n = 2;
    if (n & 1) n++;                 /* Simpson needs an even interval count */
    if (n > (1 << 24)) n = 1 << 24; /* keep the 4n node walk finite */

    Acc acc;
    memset(&acc, 0, sizeof(acc));
    double s1, s2, s4;
    simpson_triplet(m, t0, t1, n, &s1, &s2, &s4, &acc);

    Transmittance r;
    memset(&r, 0, sizeof(r));
    r.tau = s4;
    r.T = exp(-s4);
    r.evals = acc.evals;
    r.zero = acc.positive ? 0 : 1;

    /* In Simpson's asymptotic regime s1-s2 ~ 15*C*h^4 and s2-s4 ~ 15*C*(h/2)^4,
       so the successive differences fall by ~16 and d2 = |s4 - s2| ~ 15x the
       error of s4: d2 alone is a conservative bound.  Across a jump the
       differences stop decaying, the ratio collapses, and d2 can understate
       the error; r.converged reports whether the asymptotic ratio was seen.
       When it was not, fall back to the larger of the two differences, which
       is still only heuristic -- hence the flag. */
    double d1 = fabs(s2 - s1);
    double d2 = fabs(s4 - s2);
    r.converged = (d1 == 0.0 && d2 == 0.0) || (d2 > 0.0 && d1 > 8.0 * d2);
    r.err_tau = r.converged ? d2 : fmax(d1, d2);
    r.err_T = r.err_tau;
    return r;
}

static double adaptive_rec(const Medium *m, double a, double b,
                           double fa, double fm, double fb, double whole,
                           double tol, int depth, const AdaptiveOpts *o,
                           Acc *acc, double *err) {
    double mid = 0.5 * (a + b);
    double lm = 0.5 * (a + mid);
    double rm = 0.5 * (mid + b);
    double flm = sample_sigma(m, lm, acc);
    double frm = sample_sigma(m, rm, acc);

    double left  = (mid - a) / 6.0 * (fa + 4.0 * flm + fm);
    double right = (b - mid) / 6.0 * (fm + 4.0 * frm + fb);
    double delta = left + right - whole;
    double local = fabs(delta) / 15.0;

    if (depth >= o->max_depth || (depth >= o->min_depth && local <= tol)) {
        *err += local;
        return left + right + delta / 15.0;   /* Richardson correction */
    }
    /* Passing tol/2 to both halves bounds the sum of the children; but the
       Simpson difference is also ~0 on a flat interval, so a feature lying
       between the three samples looks converged and is never seen.  The
       forced min_depth cells are what protect against that. */
    return adaptive_rec(m, a, mid, fa, flm, fm, left, tol * 0.5,
                        depth + 1, o, acc, err)
         + adaptive_rec(m, mid, b, fm, frm, fb, right, tol * 0.5,
                        depth + 1, o, acc, err);
}

Transmittance transmittance_adaptive(const Medium *m, double t0, double t1,
                                     AdaptiveOpts o) {
    if (!(t1 > t0) || !(m->max_sigma > 0.0)) return empty_result();

    if (o.min_depth < 0) o.min_depth = 0;
    if (o.max_depth < o.min_depth) o.max_depth = o.min_depth;
    if (!(o.tol_tau > 0.0)) o.tol_tau = 1e-12;

    Acc acc;
    memset(&acc, 0, sizeof(acc));

    double fa = sample_sigma(m, t0, &acc);
    double fb = sample_sigma(m, t1, &acc);
    double mid = 0.5 * (t0 + t1);
    double fm = sample_sigma(m, mid, &acc);
    double whole = (t1 - t0) / 6.0 * (fa + 4.0 * fm + fb);

    double err = 0.0;
    double tau = adaptive_rec(m, t0, t1, fa, fm, fb, whole, o.tol_tau,
                              0, &o, &acc, &err);

    Transmittance r;
    memset(&r, 0, sizeof(r));
    r.tau = tau;
    r.T = exp(-tau);
    r.err_tau = err;
    r.err_T = err;
    r.evals = acc.evals;
    r.zero = acc.positive ? 0 : 1;
    /* Adaptive tol_tau is a per-pass budget; the sampled estimate is the only
       a-posteriori measure available, so trust it only if it met the target. */
    r.converged = (r.err_tau <= o.tol_tau);
    return r;
}

/* ------------------------------------------------------------------ tests */

typedef struct { double s0; } ConstCtx;
typedef struct { double s0, k; } ExpCtx;
typedef struct { double s0, amp, L; } SinCtx;
typedef struct { double A, tc, w; } GaussCtx;
typedef struct { double s0, td; } StepCtx;

static double sig_const(double t, const void *ctx) {
    (void)t;
    return ((const ConstCtx *)ctx)->s0;
}
static double sig_exp_decay(double t, const void *ctx) {
    const ExpCtx *c = (const ExpCtx *)ctx;
    return c->s0 * exp(-c->k * t);
}
static double sig_exp_grow(double t, const void *ctx) {
    const ExpCtx *c = (const ExpCtx *)ctx;
    return c->s0 * exp(c->k * t);
}
static double sig_sin(double t, const void *ctx) {
    const SinCtx *c = (const SinCtx *)ctx;
    return c->s0 * (1.0 + c->amp * sin(VOL_PI * t / c->L));
}
static double sig_gauss(double t, const void *ctx) {
    const GaussCtx *c = (const GaussCtx *)ctx;
    double x = (t - c->tc) / c->w;
    return c->A * exp(-x * x);
}
static double sig_zero(double t, const void *ctx) { (void)t; (void)ctx; return 0.0; }
static double sig_step(double t, const void *ctx) {
    const StepCtx *c = (const StepCtx *)ctx;
    return (t < c->td) ? c->s0 : 0.0;
}

static int g_calls = 0;
static double sig_counted(double t, const void *ctx) {
    (void)t; (void)ctx;
    g_calls++;
    return 1.0;   /* must never run on the empty shortcut */
}

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    char name[80];
    int i;

    /* --- empty density: T = 1 with no integration at all --- */
    g_calls = 0;
    Medium empty = { sig_counted, NULL, 0.0 };
    Transmittance e1 = transmittance_fixed(&empty, 0.0, 10.0, 256);
    check("empty_fixed_T_one", e1.T == 1.0);
    check("empty_fixed_tau_zero", e1.tau == 0.0);
    check("empty_fixed_flag", e1.zero == 1);
    check("empty_fixed_no_evals", e1.evals == 0 && g_calls == 0);

    Transmittance e2 = transmittance_adaptive(&empty, 0.0, 10.0,
                                              adaptive_opts(1e-9, 8, 24));
    check("empty_adaptive_T_one", e2.T == 1.0 && e2.tau == 0.0);
    check("empty_adaptive_flag", e2.zero == 1);
    check("empty_adaptive_no_evals", e2.evals == 0 && g_calls == 0);

    Medium neg = { sig_counted, NULL, -2.0 };
    Transmittance e3 = transmittance_fixed(&neg, -1.0, 1.0, 64);
    check("negative_bound_is_empty", e3.T == 1.0 && e3.evals == 0 && g_calls == 0);

    /* Zero density at every sample, but a positive declared bound: the
       integral runs and exp(-0) = 1 exactly. */
    Medium ns = { sig_zero, NULL, 1.0 };
    Transmittance e4 = transmittance_fixed(&ns, 0.0, 1.0, 32);
    check("zero_field_T_one", e4.T == 1.0 && e4.tau == 0.0);
    check("zero_field_flag", e4.zero == 1);
    check("zero_field_integrated", e4.evals > 0);

    /* Zero-length segment: nothing to integrate. */
    ConstCtx cc = { 0.5 };
    Medium cm = { sig_const, &cc, 0.5 };
    Transmittance zl = transmittance_fixed(&cm, 2.0, 2.0, 32);
    check("zero_length_T_one", zl.T == 1.0 && zl.evals == 0);

    /* --- constant density: Simpson is exact, T = exp(-tau) applied once --- */
    ConstCtx c1 = { 0.7 };
    Medium m1 = { sig_const, &c1, 0.7 };
    double L1 = 3.5;
    double tau_c = c1.s0 * L1;
    Transmittance cf = transmittance_fixed(&m1, 0.0, L1, 8);
    check("const_tau_exact", fabs(cf.tau - tau_c) < 1e-12);
    check("const_err_zero", cf.err_tau < 1e-12);
    check("const_exp_form", fabs(cf.T - exp(-cf.tau)) <= 1e-15);
    check("const_T_analytic", fabs(cf.T - exp(-tau_c)) < 1e-12);
    check("const_T_in_unit", cf.T > 0.0 && cf.T <= 1.0);

    /* --- smooth exponential decay: bound vs analytic, 4th-order convergence --- */
    ExpCtx ed = { 1.5, 0.8 };
    Medium me = { sig_exp_decay, &ed, ed.s0 };
    double Le = 4.0;
    double tau_e = ed.s0 / ed.k * (1.0 - exp(-ed.k * Le));
    int steps[4] = { 16, 32, 64, 128 };
    double act[4], est[4];

    for (i = 0; i < 4; i++) {
        Transmittance r = transmittance_fixed(&me, 0.0, Le, steps[i]);
        act[i] = fabs(r.tau - tau_e);
        est[i] = r.err_tau;

        snprintf(name, sizeof(name), "exp_decay_tau_bound_n%d", steps[i]);
        check(name, act[i] <= est[i] + 1e-14);
        snprintf(name, sizeof(name), "exp_decay_T_bound_n%d", steps[i]);
        check(name, fabs(r.T - exp(-tau_e)) <= r.err_T + 1e-14);
        snprintf(name, sizeof(name), "exp_decay_exp_form_n%d", steps[i]);
        check(name, fabs(r.T - exp(-r.tau)) <= 1e-15);
    }
    check("exp_decay_error_decreases",
          act[1] < act[0] && act[2] < act[1] && act[3] < act[2]);
    check("exp_decay_order4", act[0] > 1e-8 && act[0] / act[1] > 12.0);
    check("exp_decay_est_converges",
          est[0] > est[1] && est[1] > est[2] && est[0] / est[1] > 8.0);

    /* --- exponentially growing density: far end dominates --- */
    ExpCtx eg = { 1.0, 0.5 };
    double Lg = 6.0;
    Medium mg = { sig_exp_grow, &eg, exp(0.5 * Lg) };
    double tau_g = eg.s0 / eg.k * (exp(eg.k * Lg) - 1.0);
    Transmittance fg = transmittance_fixed(&mg, 0.0, Lg, 32);
    check("grow_tau_bound", fabs(fg.tau - tau_g) <= fg.err_tau + 1e-13);
    check("grow_T_bound", fabs(fg.T - exp(-tau_g)) <= fg.err_T + 1e-13);

    /* --- smooth sinusoid: bound vs analytic --- */
    SinCtx sc = { 0.9, 0.6, 3.0 };
    Medium ms = { sig_sin, &sc, sc.s0 * 1.6 };
    double tau_s = sc.s0 * (sc.L + sc.amp * 2.0 * sc.L / VOL_PI);
    Transmittance fs = transmittance_fixed(&ms, 0.0, sc.L, 64);
    check("sin_tau_bound", fabs(fs.tau - tau_s) <= fs.err_tau + 1e-13);
    check("sin_T_bound", fabs(fs.T - exp(-tau_s)) <= fs.err_T + 1e-13);

    /* --- adaptive on a smooth density: accuracy, bound, cost --- */
    Transmittance ag = transmittance_adaptive(&mg, 0.0, Lg,
                                              adaptive_opts(1e-6, 4, 24));
    double ag_err = fabs(ag.tau - tau_g);
    check("adaptive_exp_accurate", ag_err < 1e-9);
    check("adaptive_exp_bound", ag_err <= ag.err_tau + 1e-12);
    Transmittance ag16 = transmittance_fixed(&mg, 0.0, Lg, 16);
    check("adaptive_beats_fixed16",
          ag_err < fabs(ag16.tau - tau_g) && ag.evals < 4 * ag16.evals);

    /* --- narrow feature at the far end: the fixed grid's bias --- *
     * Gaussian pulse, 1/sigma half-width 2.5e-3 of the segment; grid spacing
     * at n=32 is 3.1e-2, so the pulse lives between two nodes.  min_depth=2
     * (the default) is far too coarse for it, but the error estimate goes to
     * zero instead of exposing the miss -- the flat-interval failure. */
    GaussCtx gp = { 1500.0, 0.97, 0.0025 };
    Medium mp = { sig_gauss, &gp, gp.A };
    double tau_p = gp.A * gp.w * sqrt(VOL_PI) * 0.5 *
                   (erf((1.0 - gp.tc) / gp.w) - erf((0.0 - gp.tc) / gp.w));
    Transmittance fp = transmittance_fixed(&mp, 0.0, 1.0, 32);
    check("far_pulse_fixed_biased", fabs(fp.T - exp(-tau_p)) / exp(-tau_p) > 1e-2);

    Transmittance aps = transmittance_adaptive(&mp, 0.0, 1.0,
                                              adaptive_opts(1e-6, 2, 24));
    check("far_pulse_aliased_bound_zero", aps.err_tau < 1e-9);
    check("far_pulse_aliased_wrong", fabs(aps.tau - tau_p) / tau_p > 0.1);
    /* The estimator is honest about what it saw, not about what it missed:
       with the pulse between all samples it reports convergence.  Only the
       forced min_depth and a density bound (max_sigma) can catch this. */
    check("far_pulse_aliased_false_converged", aps.converged == 1);

    Transmittance ap = transmittance_adaptive(&mp, 0.0, 1.0,
                                              adaptive_opts(1e-6, 8, 24));
    double ap_err = fabs(ap.T - exp(-tau_p));
    check("far_pulse_adaptive_accurate", ap_err < 1e-6);
    check("far_pulse_adaptive_bound", ap_err <= ap.err_T + 1e-9);

    /* A discontinuous density puts the whole optical depth into one step.
       Simpson degrades to first order across the jump, so the reported bound
       no longer shrinks as h^4; `converged` is clear and the estimate is
       still >= the true error. */
    StepCtx st = { 2.0, 0.4137 };
    Medium mst = { sig_step, &st, st.s0 };
    double tau_st = st.s0 * st.td;
    Transmittance s64 = transmittance_fixed(&mst, 0.0, 1.0, 64);
    Transmittance s128 = transmittance_fixed(&mst, 0.0, 1.0, 128);
    check("step_bound_holds",
          fabs(s64.tau - tau_st) <= s64.err_tau + 1e-13 &&
          fabs(s128.tau - tau_st) <= s128.err_tau + 1e-13);
    check("step_not_converged", !s64.converged && !s128.converged);
    check("step_bound_slower_than_order4",
          s64.err_tau / s128.err_tau < 10.0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
