
/* rr: Russian roulette path termination vs fixed-depth truncation.
   Toy path integral: each bounce scales throughput by r and adds e*throughput,
   so the analytic value is r*e/(1-r). Fixed depth drops the tail (biased);
   Russian roulette survives with probability p and compensates weight by 1/p,
   keeping the estimator unbiased while lowering the average bounce count. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_DEPTH 64   /* safety cap only; r^64 is far below test tolerances */
#define N_PATHS 20000

typedef struct { unsigned state; } Rng;
static void rng_seed(Rng* g, unsigned s) { g->state = s ? s : 1u; }
static double rng_next(Rng* g) {
    g->state = g->state * 1664525u + 1013904223u;
    return (double)(g->state >> 8) / 16777216.0;
}

typedef struct { double sum, sum2, bounces; int n; } Est;
static void est_add(Est* s, double v, int b) {
    s->sum += v; s->sum2 += v * v; s->bounces += b; s->n++;
}
static double est_mean(Est* s) { return s->sum / s->n; }
static double est_bounces(Est* s) { return s->bounces / s->n; }

/* survival probability must lie in (0, 1]; reject zero, negative, or > 1 */
static int rr_valid(double p) { return (p > 0.0 && p <= 1.0) ? 0 : -1; }

static double analytic(double r, double e) { return r * e / (1.0 - r); }

/* fixed depth: exactly `depth` bounces, tail dropped -> biased low */
static double trace_fixed(double r, double e, int depth, int* bounces) {
    double w = 1.0, L = 0.0;
    for (int k = 0; k < depth; k++) {
        w *= r;
        L += w * e;
    }
    *bounces = depth;
    return L;
}

/* russian roulette: after each bounce, survive with prob p, weight *= 1/p */
static double trace_rr(Rng* g, double r, double e, double p, int* bounces) {
    double w = 1.0, L = 0.0;
    int k;
    for (k = 1; k <= MAX_DEPTH; k++) {
        w *= r;
        L += w * e;
        if (k == MAX_DEPTH || rng_next(g) >= p) break;
        w /= p;
    }
    *bounces = k;
    return L;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    check("rr_reject_zero", rr_valid(0.0) == -1);
    check("rr_reject_negative", rr_valid(-0.5) == -1);
    check("rr_reject_over_one", rr_valid(1.0 + 1e-9) == -1);
    check("rr_accept_valid", rr_valid(0.9) == 0 && rr_valid(1.0) == 0);

    const double r = 0.8, e = 1.0;
    const double truth = analytic(r, e); /* 4.0 */

    /* p = 1: never killed, must match the truncated series at MAX_DEPTH */
    Rng g1; rng_seed(&g1, 7);
    int b1 = 0;
    double L1 = trace_rr(&g1, r, e, 1.0, &b1);
    check("rr_p1_full_depth", b1 == MAX_DEPTH);
    check("rr_p1_matches_series", fabs(L1 - truth) < 1e-4);

    /* fixed depth 5: deterministic, truncates tail -> visibly biased */
    Est fixed; memset(&fixed, 0, sizeof(fixed));
    for (int i = 0; i < N_PATHS; i++) {
        int b; double L = trace_fixed(r, e, 5, &b);
        est_add(&fixed, L, b);
    }
    check("fixed_depth_biased", fabs(est_mean(&fixed) - truth) > 0.5);

    /* russian roulette, p = 0.9: unbiased, fewer bounces on average */
    Est rr; memset(&rr, 0, sizeof(rr));
    Rng g2; rng_seed(&g2, 12345);
    for (int i = 0; i < N_PATHS; i++) {
        int b; double L = trace_rr(&g2, r, e, 0.9, &b);
        est_add(&rr, L, b);
    }
    check("rr_unbiased", fabs(est_mean(&rr) - truth) < 0.06);

    /* fixed depth 20 as the deep-but-costly baseline for bounce counts */
    Est deep; memset(&deep, 0, sizeof(deep));
    for (int i = 0; i < N_PATHS; i++) {
        int b; double L = trace_fixed(r, e, 20, &b);
        est_add(&deep, L, b);
    }
    check("rr_fewer_bounces", est_bounces(&rr) < est_bounces(&deep));

    printf("data truth=%.4f\n", truth);
    printf("data fixed_d5 mean=%.4f avg_bounces=%.2f bias=%.4f\n",
           est_mean(&fixed), est_bounces(&fixed), est_mean(&fixed) - truth);
    printf("data fixed_d20 mean=%.4f avg_bounces=%.2f bias=%.4f\n",
           est_mean(&deep), est_bounces(&deep), est_mean(&deep) - truth);
    printf("data rr_p0.9 mean=%.4f avg_bounces=%.2f bias=%.4f\n",
           est_mean(&rr), est_bounces(&rr), est_mean(&rr) - truth);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
