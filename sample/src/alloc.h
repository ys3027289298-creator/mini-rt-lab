#ifndef ALLOC_H
#define ALLOC_H

#define ALLOC_MAX_LIGHTS 64

/*
 * Per-light sampling statistics.
 *
 * var is the per-sample variance of the light's contribution estimator
 * (sigma^2).  cost is the relative cost of one sample in budget units.
 * Missing/invalid estimates (var_known == 0, negative, NaN or infinite)
 * are replaced by a conservative fallback and reported in AllocResult.
 */
typedef struct {
    double var;
    long long cost;
    int var_known;
} LightInput;

typedef struct {
    int n_lights;
    int samples[ALLOC_MAX_LIGHTS];      /* integer allocation */
    double continuous[ALLOC_MAX_LIGHTS];/* continuous optimum, n_i* */
    double effective_var[ALLOC_MAX_LIGHTS]; /* variance actually used */
    double theoretical_var;  /* variance of the continuous optimum (lower bound) */
    double allocated_var;    /* predicted variance of the integer allocation */
    double fallback_var;     /* conservative variance substituted (0 if none) */
    unsigned long long fallback_mask; /* bit i set when light i used fallback */
    int fallback_count;
    long long cost_used;
    int feasible;            /* 0 if budget cannot cover min samples */
} AllocResult;

/* Choose the conservative variance for missing estimates.  requested > 0 is
 * used directly; otherwise the largest known variance is used, or 1.0 when
 * nothing is known.  Returns 0 when no light needs a fallback. */
double alloc_fallback_var(const LightInput *lights, int n, double requested);

/* Allocate samples with fixed cost budget:
 *   minimize  sum_i var_i / n_i
 *   subject to sum_i cost_i * n_i <= budget,  n_i >= min_samples.
 * The continuous optimum is n_i* proportional to sigma_i / sqrt(cost_i).
 * Missing variances use the fallback rule above and are recorded.
 * Returns 0 on success, -1 on invalid arguments or infeasible budget. */
int alloc_build(const LightInput *lights, int n, long long budget,
                int min_samples, double fallback_var, AllocResult *out);

/* sum_i var_i / samples_i (infinite if sample_i == 0 and var_i > 0). */
double alloc_variance(const double *var, const int *samples, int n);

/* Exhaustively enumerate integer allocations for small budgets.  Returns the
 * number of allocations examined, or -1 if infeasible / max_enum exceeded. */
long long alloc_bruteforce(const LightInput *lights, int n, long long budget,
                           int min_samples, double fallback_var,
                           int *best_samples, double *best_var,
                           long long max_enum);

/* Rigorous upper bound on the relative variance gap to the integer optimum:
 * (allocated_var - theoretical_var) / theoretical_var.  Since the continuous
 * optimum is a lower bound and the exhaustive optimum is feasible, the bound
 * also covers (allocated_var - brute_force_var) / brute_force_var. */
double alloc_deviation_bound(const AllocResult *r);

#endif
