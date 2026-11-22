/* alloc_test: self-checks for cost-optimal multi-light sample allocation. */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "../src/alloc.h"

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static int near(double a, double b, double rel) {
    double scale = fabs(b) > 1e-12 ? fabs(b) : 1.0;
    return fabs(a - b) <= rel * scale;
}

/* Deterministic RNG for the empirical-variance check. */
static unsigned long long rng_state = 0x9E3779B97F4A7C15ULL;
static double rng_uniform(void) {
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    unsigned long long x = rng_state * 0x2545F4914F6CDD1DULL;
    double mantissa = (double)(x >> 11);
    return (mantissa + 0.5) * (1.0 / 9007199254740992.0);
}
static double rng_normal(void) {
    double u1 = rng_uniform(), u2 = rng_uniform();
    return sqrt(-2.0 * log(u1)) * cos(6.283185307179586 * u2);
}

static void test_continuous_proportional(void) {
    /* Equal costs, m=1: n_i* proportional to sigma_i and costs sum to B. */
    LightInput in[3] = {
        {4.0, 1, 1},
        {1.0, 1, 1},
        {0.25, 1, 1},
    };
    AllocResult r;
    int ok = alloc_build(in, 3, 1000, 1, 0.0, &r) == 0;
    double s = r.continuous[0] + r.continuous[1] + r.continuous[2];
    check("continuous_sum_matches_budget", ok && near(s, 1000.0, 1e-9));
    check("continuous_ratio_sigma",
          near(r.continuous[0] / r.continuous[1], 2.0, 1e-2) &&
              near(r.continuous[1] / r.continuous[2], 2.0, 1e-2));
    check("continuous_matches_theory_value",
          near(r.theoretical_var,
               (2.0 + 1.0 + 0.5) * (2.0 + 1.0 + 0.5) / 1000.0, 1e-6));
}

static void test_cost_weighting(void) {
    /* Cost 4 light gets half the samples of the same-variance cost-1 light. */
    LightInput in[2] = {
        {1.0, 1, 1},
        {1.0, 4, 1},
    };
    AllocResult r;
    int ok = alloc_build(in, 2, 600, 1, 0.0, &r) == 0;
    check("cost_weighting_inverse_sqrt",
          ok && near(r.continuous[0] / r.continuous[1], 2.0, 1e-2));
    check("cost_budget_respected", ok && r.cost_used <= 600);
}

static void test_integer_vs_bruteforce(void) {
    LightInput in[3] = {
        {4.0, 1, 1},
        {1.0, 2, 1},
        {0.25, 3, 1},
    };
    const long long budget = 48;
    AllocResult r;
    int best[ALLOC_MAX_LIGHTS];
    double bf_var = 0.0;
    long long examined =
        alloc_bruteforce(in, 3, budget, 2, 0.0, best, &bf_var, 20000000);
    int ok = alloc_build(in, 3, budget, 2, 0.0, &r) == 0;

    check("bruteforce_ran", examined > 0 && isfinite(bf_var));
    check("integer_respects_min", ok && r.samples[0] >= 2 && r.samples[1] >= 2 &&
                                      r.samples[2] >= 2);
    check("integer_respects_budget", ok && r.cost_used <= budget);
    check("theoretical_is_lower_bound", ok && bf_var + 1e-12 >= r.theoretical_var);
    check("allocated_near_bruteforce", ok && near(r.allocated_var, bf_var, 0.05));

    double bound = alloc_deviation_bound(&r);
    double actual_gap = (r.allocated_var - bf_var) / bf_var;
    check("deviation_bound_covers_bruteforce",
          ok && bound + 1e-12 >= actual_gap && bound >= 0.0);
    check("deviation_bound_tight", ok && bound < 0.10);

    /* The motivating failure: splitting the budget evenly by cost is
     * measurably worse than optimizing by variance. */
    int even[3];
    for (int i = 0; i < 3; i++)
        even[i] = (int)(budget / (3 * in[i].cost));
    double even_var = alloc_variance(r.effective_var, even, 3);
    check("optimal_beats_even_split",
          ok && r.allocated_var < even_var * 0.95);
}

static void test_allocation_shape(void) {
    /* Heterogeneous variances: the noisy light must get the most samples. */
    LightInput in[3] = {
        {16.0, 1, 1},
        {1.0, 1, 1},
        {0.01, 1, 1},
    };
    AllocResult r;
    int ok = alloc_build(in, 3, 900, 4, 0.0, &r) == 0;
    check("noisy_light_dominates",
          ok && r.samples[0] > r.samples[1] && r.samples[1] > r.samples[2]);
    check("quiet_light_small_share", ok && r.samples[2] < r.samples[1] / 5);
    check("samples_sum_within_budget", ok && r.cost_used <= 900);
}

static void test_zero_variance(void) {
    LightInput in[2] = {
        {0.0, 1, 1},
        {2.0, 1, 1},
    };
    AllocResult r;
    int ok = alloc_build(in, 2, 100, 3, 0.0, &r) == 0;
    check("zero_variance_gets_min", ok && r.samples[0] == 3);
    check("zero_variance_allocation_finite",
          ok && isfinite(r.allocated_var) && near(r.allocated_var, 2.0 / 97.0, 1e-9));
}

static void test_fallback_missing_variance(void) {
    LightInput in[3] = {
        {4.0, 1, 1},
        {-1.0, 1, 0},
        {1.0, 1, 0},
    };
    AllocResult r;
    int ok = alloc_build(in, 3, 600, 1, 0.0, &r) == 0;
    check("fallback_uses_max_known", ok && near(r.effective_var[1], 4.0, 1e-12) &&
                                         near(r.effective_var[2], 4.0, 1e-12));
    check("fallback_recorded",
          ok && r.fallback_count == 2 && (r.fallback_mask & 0x6u) == 0x6u &&
              near(r.fallback_var, 4.0, 1e-12));
    check("fallback_not_conservative_for_known",
          ok && (r.fallback_mask & 0x1u) == 0);
    check("fallback_allocation_feasible", ok && r.cost_used <= 600);

    /* Explicit conservative value wins over the heuristic. */
    LightInput in2[2] = {
        {1.0, 1, 1},
        {0.0, 1, 0},
    };
    AllocResult r2;
    int ok2 = alloc_build(in2, 2, 100, 1, 25.0, &r2) == 0;
    check("fallback_requested_value", ok2 && near(r2.effective_var[1], 25.0, 1e-12));

    /* NaN is treated as missing too. */
    LightInput in3[1] = {{NAN, 1, 1}};
    AllocResult r3;
    int ok3 = alloc_build(in3, 1, 10, 1, 0.0, &r3) == 0;
    check("fallback_nan", ok3 && r3.fallback_count == 1 &&
                               near(r3.effective_var[0], 1.0, 1e-12));

    /* No missing estimates -> no fallback recorded. */
    LightInput all_known[2] = {
        {1.0, 1, 1},
        {2.0, 1, 1},
    };
    AllocResult r4;
    int ok4 = alloc_build(all_known, 2, 100, 1, 0.0, &r4) >= 0;
    check("fallback_absent_when_known", ok4 && r4.fallback_count == 0 &&
                                            r4.fallback_mask == 0 &&
                                            r4.fallback_var == 0.0);

    double fv = alloc_fallback_var(in, 3, 0.0);
    check("fallback_helper", near(fv, 4.0, 1e-12));
}

static void test_infeasible(void) {
    LightInput in[2] = {
        {1.0, 10, 1},
        {1.0, 10, 1},
    };
    AllocResult r;
    int rc = alloc_build(in, 2, 15, 1, 0.0, &r);
    check("infeasible_rejected", rc == -1 && r.feasible == 0);
    check("infeasible_reports_cost", r.cost_used == 20);

    int bad = alloc_build(in, 2, 100, 0, 0.0, &r);
    check("zero_min_rejected", bad == -1);
}

static void test_empirical_variance(void) {
    LightInput in[3] = {
        {4.0, 1, 1},
        {1.0, 1, 1},
        {0.25, 4, 1},
    };
    const int reps = 20000;
    AllocResult r;
    int ok = alloc_build(in, 3, 1200, 8, 0.0, &r) == 0;
    if (!ok) { check("empirical_setup", 0); return; }

    /* Each replication draws n_i normal samples per light and averages them;
     * the estimator variance should match sum var_i / n_i. */
    double sum = 0.0, sum2 = 0.0;
    for (int k = 0; k < reps; k++) {
        double est = 0.0;
        for (int i = 0; i < 3; i++) {
            double mean = 0.0;
            for (int s = 0; s < r.samples[i]; s++)
                mean += sqrt(r.effective_var[i]) * rng_normal();
            est += mean / (double)r.samples[i];
        }
        sum += est;
        sum2 += est * est;
    }
    double mean = sum / reps;
    double emp = (sum2 / reps - mean * mean) * reps / (reps - 1);
    check("empirical_matches_theory", near(emp, r.allocated_var, 0.08));
    check("empirical_above_continuous_bound",
          emp >= r.theoretical_var * 0.97);
    check("empirical_mean_near_zero", fabs(mean) < 0.05);
}

int main(void) {
    test_continuous_proportional();
    test_cost_weighting();
    test_integer_vs_bruteforce();
    test_allocation_shape();
    test_zero_variance();
    test_fallback_missing_variance();
    test_infeasible();
    test_empirical_variance();

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
