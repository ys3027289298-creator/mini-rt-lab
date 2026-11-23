/* site_stats.c -- accuracy of one site over many work items.
 *
 * Prints mean, analytic target, signed bias, standard error and RMSE, plus a
 * verdict against a 5-standard-error bound.  Run with the legacy build and the
 * unified build to compare before/after accuracy:
 *
 *   site_stats 0 SITE_SIGMA   (default 1.7)
 *   site_stats 1 SITE_N       (default 8)
 *   site_stats 2 SITE_STRATA  (default 16)
 *
 * Usage: site_stats <site> <param> [items]
 *
 * Always single-threaded: the legacy path has one ordered stream, so that is
 * the only setting in which it has a defined value at all.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "rt_kernel.h"

int main(int argc, char** argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <site 0|1|2> <param> [items]\n", argv[0]);
        return 2;
    }
    int site = atoi(argv[1]);
    double param = atof(argv[2]);
    int n = argc > 3 ? atoi(argv[3]) : 20000;
    if (n < 2) return 2;

    rt_rng_configure(0x5EEDu);

    double target;
    if (site == 0)      target = 0.5;   /* E[exp(-sigma t)], t~Exp(sigma) */
    else if (site == 1) target = 2.0 / 3.0; /* E[cos] under cosine sampling */
    else                target = 1.0 / 3.0; /* integral u^2 over [0,1) */

    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; i++) {
        double v;
        if (site == 0)      v = rt_atten_sample(param, (uint64_t)i).weight;
        else if (site == 1) v = rt_material_mean_cos((int)param, (uint64_t)i);
        else                v = rt_strat_quad((int)param, (uint64_t)i);
        if (!isfinite(v)) { printf("nonfinite=1\n"); return 1; }
        sum += v; sum2 += v * v;
    }
    double mean = sum / n;
    double var = sum2 / n - mean * mean;
    double se = var > 0 ? sqrt(var / n) : 0.0;
    double bias = mean - target;
    double rmse = sqrt(bias * bias + var / n);

    printf("site=%d param=%.17g items=%d\n", site, param, n);
    printf("mean=%.17g target=%.17g bias=%.6g stderr=%.6g rmse=%.6g\n",
           mean, target, bias, se, rmse);
    printf("verdict=%s\n", fabs(bias) <= 5.0 * se + 1e-12 ? "PASS" : "FAIL");
    return fabs(bias) <= 5.0 * se + 1e-12 ? 0 : 1;
}
