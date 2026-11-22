
/* density: optimal adaptive sampling density from per-cell error contribution.
 *
 * A region is split into k cells.  Cell i has local standard deviation sigma_i,
 * measure weight w_i > 0, and receives N_i samples with N = sum N_i.  Its error
 * contribution is modelled as
 *
 *     e_i(N_i) = w_i * sigma_i^a / N_i^b .
 *
 * Two objectives:
 *     variance  (a = 2, b = 1  ): e_i is the estimator variance
 *     std       (a = 1, b = 1/2): e_i is the standard error
 *
 * Minimising sum_i e_i under sum_i N_i = N gives N_i proportional to
 * c_i^(1/(b+1)) with c_i = w_i sigma_i^a, i.e. density
 *
 *     rho_i proportional to sigma_i^(a/(b+1)) = sigma_i^p
 *
 * for equal weights (p = 1 for variance, p = 2/3 for std).  The continuous
 * optimum has error E* = (sum_i c_i^(1/(b+1)))^(b+1) / N^b; uniform density
 * rho_i = 1/k gives E_u = (k/N)^b sum_i c_i, hence
 *
 *     E_u / E* = k^b sum c_i / (sum c_i^(1/(b+1)))^(b+1) >= 1
 *
 * by Holder, with equality iff all c_i are equal.  With x_i = N*rho_i the
 * continuous shares and n_i the integer allocation, writing delta_i = n_i/x_i-1,
 *
 *     E_int / E* = sum_i rho_i * (1 + delta_i)^(-b) ,
 *
 * so if every active cell satisfies n_i >= x_i - 1 then delta_i >= -1/x_i and
 *
 *     E_int / E* <= (1 - 1/min_i x_i)^(-b)     for min_i x_i > 1.
 *
 * The largest-remainder allocation guarantees |n_i - x_i| < 1 and (once the
 * budget covers every active cell) n_i >= 1, so the bound above is certified
 * whenever min_share = min_i x_i > 1.  The exhaustive integer optimum lies
 * between E* and E_int; the self-test checks all three against enumeration.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define K_MAX 32

typedef enum { OBJ_VARIANCE = 0, OBJ_STD = 1 } Objective;
typedef enum { DENS_OK = 0, DENS_EMPTY = -1, DENS_BAD_INPUT = -2, DENS_DEGENERATE = -3,
               DENS_UNDERFUNDED = -4 } DensStatus;

typedef struct {
    int k;
    double sigma[K_MAX];   /* local standard deviation, >= 0 */
    double weight[K_MAX];  /* cell measure, > 0 */
} Region;

typedef struct {
    DensStatus status;
    int k;
    Objective obj;
    int budget;
    double a, b, p;        /* objective exponents and density power */
    double coef[K_MAX];    /* c_i = w_i * sigma_i^a */
    double rho[K_MAX];     /* optimal continuous density, sums to 1 */
    int n_int[K_MAX];      /* largest-remainder integer allocation */
    double err_cont;       /* continuous optimal error E* */
    double err_uniform;    /* uniform-density error E_u */
    double err_int;        /* rounded integer allocation error */
    double ratio_uniform;  /* E_u / E* */
    double ratio_int;      /* E_int / E_int-continuous */
    double min_share;      /* min_i N * rho_i */
    double round_bound;    /* upper bound on ratio_int, INF if not certified */
} Dens;

static void obj_exponents(Objective obj, double* a, double* b) {
    if (obj == OBJ_VARIANCE) { *a = 2.0; *b = 1.0; }
    else                     { *a = 1.0; *b = 0.5; }
}

static DensStatus region_check(const Region* r) {
    if (r->k <= 0) return DENS_EMPTY;
    double wsum = 0.0;
    for (int i = 0; i < r->k; i++) {
        if (r->weight[i] < 0.0 || r->sigma[i] < 0.0) return DENS_BAD_INPUT;
        wsum += r->weight[i];
    }
    if (wsum <= 0.0) return DENS_EMPTY;   /* zero measure: empty region */
    return DENS_OK;
}

static double err_int_alloc(const Dens* d) {
    double e = 0.0;
    for (int i = 0; i < d->k; i++) {
        if (d->coef[i] == 0.0) continue;
        if (d->n_int[i] <= 0) return INFINITY;
        e += d->coef[i] / pow((double)d->n_int[i], d->b);
    }
    return e;
}

static void dens_solve(const Region* r, Objective obj, int budget, Dens* d) {
    memset(d, 0, sizeof(*d));
    d->k = r->k;
    d->obj = obj;
    d->budget = budget;
    obj_exponents(obj, &d->a, &d->b);
    d->p = d->a / (d->b + 1.0);
    d->round_bound = INFINITY;

    DensStatus st = region_check(r);
    if (st != DENS_OK) { d->status = st; return; }
    if (budget <= 0) { d->status = DENS_BAD_INPUT; return; }

    /* c_i and the optimal share numerator t_i = c_i^(1/(b+1)). */
    double inv = 1.0 / (d->b + 1.0);
    double t[K_MAX], T = 0.0, csum = 0.0;
    for (int i = 0; i < d->k; i++) {
        d->coef[i] = r->weight[i] * pow(r->sigma[i], d->a);
        t[i] = pow(d->coef[i], inv);
        T += t[i];
        csum += d->coef[i];
    }
    if (T <= 0.0) { d->status = DENS_DEGENERATE; return; }  /* nothing to allocate */

    d->status = DENS_OK;
    for (int i = 0; i < d->k; i++) d->rho[i] = t[i] / T;

    d->err_cont    = pow(T, d->b + 1.0) / pow((double)budget, d->b);
    d->err_uniform = pow((double)d->k / (double)budget, d->b) * csum;
    d->ratio_uniform = d->err_cont > 0.0 ? d->err_uniform / d->err_cont : 0.0;

    d->min_share = INFINITY;
    for (int i = 0; i < d->k; i++) {
        if (d->coef[i] == 0.0) continue;
        double share = (double)budget * d->rho[i];
        if (share < d->min_share) d->min_share = share;
    }

    /* Largest-remainder integer allocation: floor every share, then hand the
       residual to the largest fractional parts.  A cell whose continuous share
       is below one sample must still receive one, so seed those first; if the
       budget cannot cover every active cell the allocation is infeasible. */
    int active = 0;
    for (int i = 0; i < d->k; i++) if (d->coef[i] > 0.0) active++;
    if (budget < active) { d->status = DENS_UNDERFUNDED; return; }

    int rem = budget;
    double frac[K_MAX];
    for (int i = 0; i < d->k; i++) {
        double x = (double)budget * d->rho[i];
        d->n_int[i] = (int)floor(x);
        if (d->coef[i] > 0.0 && d->n_int[i] < 1) d->n_int[i] = 1;
        rem -= d->n_int[i];
        frac[i] = x - floor(x);
    }
    for (int left = rem; left > 0; left--) {
        int best = -1;
        double bf = -1.0;
        for (int i = 0; i < d->k; i++) if (frac[i] > bf) { bf = frac[i]; best = i; }
        if (best < 0) break;
        d->n_int[best]++;
        frac[best] = -1.0;
    }

    d->err_int = err_int_alloc(d);
    d->ratio_int = d->err_cont > 0.0 ? d->err_int / d->err_cont : 0.0;

    if (d->min_share > 1.0)
        d->round_bound = pow(1.0 / (1.0 - 1.0 / d->min_share), d->b);
}

/* ---- brute-force integer optimum, for the discrete-bound test ---- */

typedef struct { const Dens* d; int n[K_MAX]; double best; } Enum;

static void enum_rec(Enum* e, int idx, int remaining) {
    const Dens* d = e->d;
    if (idx == d->k) {
        if (remaining != 0) return;
        double err = 0.0;
        for (int i = 0; i < d->k; i++) {
            if (d->coef[i] == 0.0) continue;
            if (e->n[i] <= 0) return;
            err += d->coef[i] / pow((double)e->n[i], d->b);
        }
        if (err < e->best) e->best = err;
        return;
    }
    int rest = d->k - idx - 1;
    for (int v = 1; v <= remaining - rest; v++) {
        e->n[idx] = v;
        enum_rec(e, idx + 1, remaining - v);
    }
}

static double exhaustive_best(const Dens* d, int budget) {
    if (budget < d->k) return INFINITY;
    Enum e;
    memset(&e, 0, sizeof(e));
    e.d = d;
    e.best = INFINITY;
    enum_rec(&e, 0, budget);
    return e.best;
}

/* ---- self-test ---- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static Region make_region(const double* sig, const double* w, int k) {
    Region r;
    memset(&r, 0, sizeof(r));
    r.k = k;
    for (int i = 0; i < k; i++) { r.sigma[i] = sig[i]; r.weight[i] = w ? w[i] : 1.0; }
    return r;
}

static int near(double x, double y, double tol) { return fabs(x - y) <= tol; }

int main(void) {
    /* empty / invalid regions are rejected */
    Region empty;
    memset(&empty, 0, sizeof(empty));
    Dens d0;
    dens_solve(&empty, OBJ_VARIANCE, 16, &d0);
    check("reject_empty_k0", d0.status == DENS_EMPTY);

    double zero_w[3] = { 0.0, 0.0, 0.0 };
    double sig3[3] = { 1.0, 2.0, 3.0 };
    Region rzw = make_region(sig3, zero_w, 3);
    dens_solve(&rzw, OBJ_VARIANCE, 16, &d0);
    check("reject_zero_measure", d0.status == DENS_EMPTY);

    double neg[2] = { 1.0, -1.0 };
    Region rneg = make_region(neg, NULL, 2);
    dens_solve(&rneg, OBJ_VARIANCE, 16, &d0);
    check("reject_negative_sigma", d0.status == DENS_BAD_INPUT);

    double zero_sig[2] = { 0.0, 0.0 };
    Region rzs = make_region(zero_sig, NULL, 2);
    dens_solve(&rzs, OBJ_VARIANCE, 16, &d0);
    check("reject_degenerate", d0.status == DENS_DEGENERATE);

    /* density is proportional to sigma^p */
    double s124[3] = { 1.0, 2.0, 4.0 };
    Region r124 = make_region(s124, NULL, 3);
    Dens dv, ds;
    dens_solve(&r124, OBJ_VARIANCE, 24, &dv);
    dens_solve(&r124, OBJ_STD, 24, &ds);

    printf("info_p_variance=%.6f\n", dv.p);
    printf("info_p_std=%.6f\n", ds.p);
    check("power_variance_p1", near(dv.p, 1.0, 1e-12));
    check("power_std_p2_3", near(ds.p, 2.0 / 3.0, 1e-12));

    check("density_prop_sigma_var",
          near(dv.rho[1] / dv.rho[0], 2.0, 1e-9) &&
          near(dv.rho[2] / dv.rho[0], 4.0, 1e-9));
    check("density_prop_sigma_23",
          near(ds.rho[1] / ds.rho[0], pow(2.0, 2.0 / 3.0), 1e-9) &&
          near(ds.rho[2] / ds.rho[0], pow(4.0, 2.0 / 3.0), 1e-9));

    /* first-order optimality: N_i^(b+1) / c_i is equal across cells */
    for (int pass = 0; pass < 2; pass++) {
        Dens* d = pass ? &ds : &dv;
        double qmin = INFINITY, qmax = 0.0;
        for (int i = 0; i < d->k; i++) {
            if (d->coef[i] == 0.0) continue;
            double N = (double)d->budget * d->rho[i];
            double q = pow(N, d->b + 1.0) / d->coef[i];
            if (q < qmin) qmin = q;
            if (q > qmax) qmax = q;
        }
        check(pass ? "optimality_std" : "optimality_variance", qmax / qmin - 1.0 < 1e-12);
    }

    /* optimal density beats uniform; ratio matches the closed form */
    double s113[4] = { 1.0, 1.0, 1.0, 3.0 };
    Region r113 = make_region(s113, NULL, 4);
    Dens du, du2;
    dens_solve(&r113, OBJ_VARIANCE, 60, &du);
    dens_solve(&r113, OBJ_STD, 60, &du2);

    double closed_var = 4.0 * (1 + 1 + 1 + 9) / (1 + 1 + 1 + 3) / (1 + 1 + 1 + 3);
    double a_sum = 1 + 1 + 1 + 3;
    double t_sum = 3 * pow(1.0, 2.0 / 3.0) + pow(3.0, 2.0 / 3.0);
    double closed_std = sqrt(4.0) * a_sum / pow(t_sum, 1.5);

    printf("info_ratio_uniform_var=%.9f\n", du.ratio_uniform);
    printf("info_ratio_uniform_std=%.9f\n", du2.ratio_uniform);
    check("uniform_ratio_matches_var", near(du.ratio_uniform, closed_var, 1e-9));
    check("uniform_ratio_matches_std", near(du2.ratio_uniform, closed_std, 1e-9));
    check("density_beats_uniform_var", du.ratio_uniform > 1.0 + 1e-9);
    check("density_beats_uniform_std", du2.ratio_uniform > 1.0 + 1e-9);
    printf("ratio_uniform_gain_var=%.6f\n", du.ratio_uniform);
    printf("ratio_uniform_gain_std=%.6f\n", du2.ratio_uniform);

    /* homogeneous cells: uniform density is already optimal (ratio == 1) */
    double sh[4] = { 2.0, 2.0, 2.0, 2.0 };
    Region rh = make_region(sh, NULL, 4);
    Dens dh, dh2;
    dens_solve(&rh, OBJ_VARIANCE, 40, &dh);
    dens_solve(&rh, OBJ_STD, 40, &dh2);
    check("uniform_optimal_when_equal_var", near(dh.ratio_uniform, 1.0, 1e-12));
    check("uniform_optimal_when_equal_std", near(dh2.ratio_uniform, 1.0, 1e-12));

    /* discrete allocation: continuous E* is a lower bound, the rounded
       allocation and the exhaustive optimum both respect the bound. */
    double s1248[4] = { 1.0, 2.0, 4.0, 8.0 };
    Region r1248 = make_region(s1248, NULL, 4);
    Dens dtest[2];
    dens_solve(&r1248, OBJ_VARIANCE, 32, &dtest[0]);
    dens_solve(&r1248, OBJ_STD, 32, &dtest[1]);
    for (int pass = 0; pass < 2; pass++) {
        Dens* d = &dtest[pass];
        double ex = exhaustive_best(d, d->budget);
        const char* tag = pass ? "std" : "var";
        printf("info_%s_min_share=%.6f\n", tag, d->min_share);
        printf("info_%s_round_bound=%.6f\n", tag, d->round_bound);
        printf("info_%s_ratio_int=%.9f\n", tag, d->ratio_int);
        printf("info_%s_exhaustive_ratio=%.9f\n", tag, ex / d->err_cont);
        check(pass ? "discrete_bound_std" : "discrete_bound_var",
              d->round_bound > 1.0 && d->ratio_int <= d->round_bound + 1e-12);
        check(pass ? "discrete_cont_lower_std" : "discrete_cont_lower_var",
              d->err_cont <= ex + 1e-12);
        check(pass ? "discrete_exhaustive_le_round_std" : "discrete_exhaustive_le_round_var",
              ex <= d->err_int + 1e-12);
        check(pass ? "discrete_ratio_ge_one_std" : "discrete_ratio_ge_one_var",
              d->ratio_int >= 1.0 - 1e-12);
    }

    /* more budget strictly lowers the continuous error */
    Dens lo, hi;
    dens_solve(&r1248, OBJ_VARIANCE, 16, &lo);
    dens_solve(&r1248, OBJ_VARIANCE, 64, &hi);
    check("error_drops_with_budget", hi.err_cont < lo.err_cont);

    /* a budget smaller than the number of active cells cannot be allocated */
    Dens under;
    dens_solve(&r1248, OBJ_VARIANCE, 3, &under);
    check("reject_underfunded", under.status == DENS_UNDERFUNDED);

    /* Case where largest-remainder rounding is strictly worse than the
       exhaustive integer optimum, yet still inside the certified bound:
       the gap is real, so the exhaustive < rounded check has teeth. */
    double sg[4] = { 1.0, 1.5, 2.0, 3.2 };
    Region rgap = make_region(sg, NULL, 4);
    Dens dg;
    dens_solve(&rgap, OBJ_VARIANCE, 32, &dg);
    double exg = exhaustive_best(&dg, dg.budget);
    printf("info_gap_min_share=%.6f\n", dg.min_share);
    printf("info_gap_round_bound=%.6f\n", dg.round_bound);
    printf("info_gap_ratio_int=%.9f\n", dg.ratio_int);
    printf("info_gap_exhaustive_ratio=%.9f\n", exg / dg.err_cont);
    check("gap_exhaustive_strictly_better", exg < dg.err_int - 1e-12);
    check("gap_bound_certified", dg.round_bound > 1.0 && dg.ratio_int <= dg.round_bound + 1e-12);
    check("gap_exhaustive_within_bound", exg / dg.err_cont <= dg.round_bound + 1e-12);

    /* the integer allocation covers every active cell and sums to the budget */
    int sum = 0, covered = 1;
    for (int i = 0; i < dv.k; i++) {
        sum += dv.n_int[i];
        if (dv.coef[i] > 0.0 && dv.n_int[i] < 1) covered = 0;
    }
    check("int_allocation_sums", sum == dv.budget);
    check("int_allocation_covers", covered);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
