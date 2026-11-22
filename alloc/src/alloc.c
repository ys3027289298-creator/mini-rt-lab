
/* alloc: adaptive per-pixel dimension (sample) budget allocation.
 *
 * Giving every pixel the same number of dimensions is wasteful: a
 * low-complexity pixel (small pilot standard deviation sigma) converges
 * faster than a high-complexity one.  With n_i dimensions the variance of
 * pixel i's mean estimate is the error-vs-dimension relation
 *
 *     v_i(n_i) = sigma_i^2 / n_i,
 *
 * so under a global budget B (sum n_i <= B, every n_i >= 1) the total
 * MSE  E = sum sigma_i^2 / n_i  is minimized when n_i is proportional to
 * sigma_i.  We find that water-filling level, floor to integers, then hand
 * out the remaining dimensions one at a time by largest marginal error
 * reduction sigma_i^2 / (n_i (n_i + 1)).  Zero dimensions are never
 * produced, and requests that cannot give every pixel at least one
 * dimension are rejected.
 */
#include <stdio.h>
#include <math.h>

#define ALLOC_MAX_PIX 64
#define ALLOC_MAX_BUDGET 1000000

enum { ALLOC_OK = 0, ALLOC_EBAD = -1 };

/* unbiased sample standard deviation from a pilot pass (0 for k < 2) */
static double pilot_std(const double* x, int k) {
    if (k < 2) return 0.0;
    double mean = 0.0;
    for (int i = 0; i < k; i++) mean += x[i];
    mean /= (double)k;
    double acc = 0.0;
    for (int i = 0; i < k; i++) { double d = x[i] - mean; acc += d * d; }
    double var = acc / (double)(k - 1);
    return sqrt(var > 0.0 ? var : 0.0);
}

/* predicted per-pixel MSE: the error-vs-dimension curve v(n) = sigma^2 / n */
static double pixel_mse(double sigma, int n) {
    return sigma * sigma / (double)n;
}

static double total_mse(const double* sigma, const int* dims, int npix) {
    double e = 0.0;
    for (int i = 0; i < npix; i++) e += pixel_mse(sigma[i], dims[i]);
    return e;
}

/* uniform reference: every pixel receives budget/npix dimensions */
static double uniform_mse(const double* sigma, int npix, int budget) {
    double per = (double)budget / (double)npix;
    double e = 0.0;
    for (int i = 0; i < npix; i++) e += sigma[i] * sigma[i] / per;
    return e;
}

static int valid_inputs(const double* sigma, int npix, int budget, int n_min) {
    if (npix <= 0 || npix > ALLOC_MAX_PIX) return 0;
    if (n_min < 1) return 0;
    if (budget < npix * n_min || budget > ALLOC_MAX_BUDGET) return 0;
    for (int i = 0; i < npix; i++)
        if (!isfinite(sigma[i]) || sigma[i] < 0.0) return 0;
    return 1;
}

static double cont_sum(const double* sigma, int npix, int n_min, double lambda) {
    double s = 0.0;
    for (int i = 0; i < npix; i++) {
        double n = sigma[i] / lambda;
        s += n > (double)n_min ? n : (double)n_min;
    }
    return s;
}

/* allocate `budget` dimensions across `npix` pixels.
 * returns ALLOC_EBAD on invalid input; on success out[i] >= n_min and
 * sum(out) == budget: never over budget, never a zero-dimensional pixel. */
static int alloc_dims(const double* sigma, int npix, int budget, int n_min, int* out) {
    if (!valid_inputs(sigma, npix, budget, n_min)) return ALLOC_EBAD;

    /* water-filling level lambda: sum max(n_min, sigma_i / lambda) = budget */
    double lo = 1e-300, hi = 1.0;
    while (cont_sum(sigma, npix, n_min, hi) > (double)budget) hi *= 2.0;
    for (int it = 0; it < 200; it++) {
        double mid = sqrt(lo * hi);
        if (cont_sum(sigma, npix, n_min, mid) > (double)budget) lo = mid;
        else hi = mid;
    }
    double lambda = hi;

    long used = 0;
    for (int i = 0; i < npix; i++) {
        double target = sigma[i] / lambda;
        if (target < (double)n_min) target = (double)n_min;
        out[i] = (int)floor(target);
        if (out[i] < n_min) out[i] = n_min;
        used += out[i];
    }

    /* repair any float overshoot, then fill the remainder greedily:
     * one more dimension on pixel i cuts MSE by sigma_i^2/(n_i(n_i+1)) */
    while (used > (long)budget) {
        int pick = -1;
        double smallest_loss = 1e300;
        for (int i = 0; i < npix; i++) {
            if (out[i] <= n_min) continue;
            double loss = sigma[i] * sigma[i] / (double)(out[i] * (out[i] - 1));
            if (loss < smallest_loss) { smallest_loss = loss; pick = i; }
        }
        if (pick < 0) return ALLOC_EBAD;
        out[pick]--;
        used--;
    }
    while (used < (long)budget) {
        int pick = 0;
        double best_gain = -1.0;
        for (int i = 0; i < npix; i++) {
            double gain = sigma[i] * sigma[i] / (double)(out[i] * (out[i] + 1));
            if (gain > best_gain) { best_gain = gain; pick = i; }
        }
        out[pick]++;
        used++;
    }
    return ALLOC_OK;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    double s4[4] = {1.0, 1.0, 1.0, 9.0};
    int dims[ALLOC_MAX_PIX];

    /* ---- rejection: zero / impossible dimensions are refused ---- */
    check("reject_zero_budget", alloc_dims(s4, 4, 0, 1, dims) == ALLOC_EBAD);
    check("reject_budget_below_pixels", alloc_dims(s4, 4, 3, 1, dims) == ALLOC_EBAD);
    check("reject_zero_nmin", alloc_dims(s4, 4, 8, 0, dims) == ALLOC_EBAD);
    check("reject_bad_npix", alloc_dims(s4, 0, 8, 1, dims) == ALLOC_EBAD);
    double neg[2] = {1.0, -2.0};
    check("reject_neg_sigma", alloc_dims(neg, 2, 8, 1, dims) == ALLOC_EBAD);
    double nanv[2] = {1.0, NAN};
    check("reject_nan_sigma", alloc_dims(nanv, 2, 8, 1, dims) == ALLOC_EBAD);

    /* ---- pilot complexity estimate ---- */
    double xs[3] = {2.0, 4.0, 6.0};
    check("pilot_std_known", fabs(pilot_std(xs, 3) - 2.0) < 1e-12);
    double cst[3] = {5.0, 5.0, 5.0};
    check("pilot_std_const", pilot_std(cst, 3) == 0.0);
    check("pilot_std_single", pilot_std(xs, 1) == 0.0);

    /* ---- error vs dimension: v(n) = sigma^2 / n ---- */
    check("mse_formula", fabs(pixel_mse(3.0, 4) - 2.25) < 1e-12);
    check("mse_halves_when_dims_double", fabs(pixel_mse(3.0, 8) - 1.125) < 1e-12);
    double gain = pixel_mse(3.0, 4) - pixel_mse(3.0, 5);
    check("mse_marginal_gain", fabs(gain - 9.0 / (4.0 * 5.0)) < 1e-12);

    /* ---- heterogeneous scene: one complex pixel among three simple ones ---- */
    const int N = 4, B = 40;
    int rc = alloc_dims(s4, N, B, 1, dims);
    long used = 0;
    int no_zero = 1;
    for (int i = 0; i < N; i++) {
        used += dims[i];
        if (dims[i] < 1) no_zero = 0;
    }
    check("alloc_ok", rc == ALLOC_OK);
    check("within_budget", used <= B);
    check("budget_fully_used", used == B);
    check("no_zero_dims", no_zero);
    check("complex_gets_most", dims[3] > dims[0] && dims[0] >= 1);

    double e_adp = total_mse(s4, dims, N);
    double e_uni = uniform_mse(s4, N, B);
    double ratio = e_adp / e_uni;
    double s_sum = 12.0, q_sum = 84.0;            /* sum sigma, sum sigma^2 */
    double e_opt = s_sum * s_sum / (double)B;      /* continuous optimum */
    double ratio_opt = s_sum * s_sum / (double)N / q_sum;
    check("adaptive_beats_uniform", e_adp < e_uni);
    check("ratio_band", ratio > 0.42 && ratio < 0.45);
    check("above_continuous_opt", e_adp >= e_opt - 1e-12);
    check("opt_ratio_theory", fabs(ratio_opt - (3.6 / 8.4)) < 1e-12);

    printf("\nscene sigma={1,1,1,9} budget=%d n_min=1\n", B);
    printf(" pix  sigma  dims  mse(sigma^2/n)\n");
    for (int i = 0; i < N; i++)
        printf("%4d %6.1f %5d %14.6f\n", i, s4[i], dims[i], pixel_mse(s4[i], dims[i]));
    printf("mse_uniform=%.6f mse_adaptive=%.6f ratio=%.4f reduction=%.1f%%\n",
           e_uni, e_adp, ratio, (1.0 - ratio) * 100.0);

    /* ---- uniform scene: nothing to gain, allocation stays even ---- */
    double su[4] = {2.0, 2.0, 2.0, 2.0};
    alloc_dims(su, 4, 16, 1, dims);
    check("uniform_scene_even", dims[0] == 4 && dims[1] == 4 &&
                                dims[2] == 4 && dims[3] == 4);
    double ru = total_mse(su, dims, 4) / uniform_mse(su, 4, 16);
    check("uniform_scene_ratio_one", fabs(ru - 1.0) < 1e-12);
    alloc_dims(su, 4, 17, 1, dims);
    int evenish = 1;
    for (int i = 0; i < 4; i++)
        if (dims[i] < 4 || dims[i] > 5) evenish = 0;
    check("uniform_scene_remainder", evenish);

    /* ---- zero-complexity pixel: allowed, still gets n_min, never 0 ---- */
    double sz[2] = {0.0, 4.0};
    rc = alloc_dims(sz, 2, 10, 1, dims);
    check("zero_sigma_ok", rc == ALLOC_OK && dims[0] == 1 && dims[1] == 9);
    check("zero_sigma_mse", pixel_mse(0.0, dims[0]) == 0.0);

    /* ---- budget sweep: error never rises, ratio always <= 1, sum == B ---- */
    double prev = 1e300;
    int mono = 1, sweep_ok = 1;
    printf("\nbudget sweep on sigma={1,1,1,9}:\n");
    printf(" budget  mse_uniform  mse_adaptive   ratio   sum_dims\n");
    for (int b = 8; b <= 64; b *= 2) {
        alloc_dims(s4, 4, b, 1, dims);
        long u = 0;
        for (int i = 0; i < 4; i++) u += dims[i];
        double ea = total_mse(s4, dims, 4);
        double eu = uniform_mse(s4, 4, b);
        if (ea > prev + 1e-12) mono = 0;
        prev = ea;
        double r = ea / eu;
        if (r > 1.0 + 1e-12 || u != (long)b) sweep_ok = 0;
        printf("%7d %12.5f %13.5f %8.4f %9ld\n", b, eu, ea, r, u);
    }
    check("mse_monotone_in_budget", mono);
    check("sweep_valid", sweep_ok);
    /* integer allocation ratio converges to the continuous optimum 0.4286 */
    alloc_dims(s4, 4, 64, 1, dims);
    double r64 = total_mse(s4, dims, 4) / uniform_mse(s4, 4, 64);
    check("ratio_converges", fabs(r64 - ratio_opt) < 0.02);

    /* ---- error-vs-dimension curve for one complex pixel ---- */
    printf("\nerror vs dimensions (sigma=9):\n");
    for (int n = 1; n <= 8; n++)
        printf(" n=%d mse=%.4f\n", n, pixel_mse(9.0, n));

    printf("\nall_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
