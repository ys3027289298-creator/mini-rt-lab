/*
 * alloc.c: cost-optimal multi-light sample allocation.
 *
 * With independent lights the estimator variance is
 *
 *     V(n) = sum_i var_i / n_i
 *
 * and the shadow-ray budget is C(n) = sum_i cost_i * n_i.  Minimizing V
 * subject to C <= B and n_i >= m is separable and convex.  The continuous
 * optimum is found by water-filling: with a_i = sigma_i / sqrt(cost_i),
 *
 *     n_i(lambda) = max(m, a_i * lambda),
 *
 * and lambda is the root of sum_i cost_i * n_i(lambda) = B.  Its variance
 * V* is a lower bound on every feasible integer allocation, so it is the
 * reference for the optimality-gap bound.
 *
 * The integer allocation floors the continuous solution and then spends the
 * remaining budget greedily by marginal variance reduction per unit cost.
 */
#include "alloc.h"

#include <limits.h>
#include <math.h>
#include <string.h>

#define N_SAMPLE_CAP 1073741824 /* 2^30, keeps counts inside int */

/* A known variance may legitimately be zero (e.g. a light fully occluded from
 * the current shading point).  Only negative, NaN or infinite values are
 * unusable and fall back to the conservative estimate. */
static int valid_var(double v) { return v >= 0.0 && isfinite(v); }

double alloc_fallback_var(const LightInput *lights, int n, double requested) {
    int need = 0;
    double max_known = 0.0;
    for (int i = 0; i < n; i++)
        if (!lights[i].var_known || !valid_var(lights[i].var)) {
            need = 1;
        } else if (lights[i].var > max_known) {
            max_known = lights[i].var;
        }
    if (!need) return 0.0;
    if (requested > 0.0 && isfinite(requested)) return requested;
    return max_known > 0.0 ? max_known : 1.0;
}

/* Cost of the continuous water-filling solution with a given lambda. */
static double water_cost(const double *a, const long long *cost, int n,
                         int m, double lambda) {
    double total = 0.0;
    for (int i = 0; i < n; i++) {
        double ni = a[i] * lambda;
        if (ni < (double)m) ni = (double)m;
        total += (double)cost[i] * ni;
    }
    return total;
}

/* Continuous optimum variance subject to n_i >= m and the cost budget.
 * Returns V*; fills n_i* into cont (may be NULL). */
static double continuous_opt(const double *var, const long long *cost, int n,
                             int m, long long budget, double *cont) {
    double a[ALLOC_MAX_LIGHTS];
    int any_var = 0;
    for (int i = 0; i < n; i++) {
        a[i] = sqrt(var[i]) / sqrt((double)cost[i]);
        if (var[i] > 0.0) any_var = 1;
    }
    if (!any_var) {
        for (int i = 0; i < n; i++) if (cont) cont[i] = (double)m;
        return 0.0;
    }

    /* Bracket lambda: f(0) = cost(m); grow hi until it reaches budget. */
    double lo = 0.0, hi = 1.0;
    while (water_cost(a, cost, n, m, hi) < (double)budget) {
        lo = hi;
        hi *= 2.0;
        if (hi > 1e300) break;
    }
    for (int it = 0; it < 200; it++) {
        double mid = 0.5 * (lo + hi);
        if (water_cost(a, cost, n, m, mid) < (double)budget)
            lo = mid;
        else
            hi = mid;
    }
    double lambda = 0.5 * (lo + hi);

    double v = 0.0;
    for (int i = 0; i < n; i++) {
        double ni = a[i] * lambda;
        if (ni < (double)m) ni = (double)m;
        if (cont) cont[i] = ni;
        if (var[i] > 0.0) v += (ni > 0.0) ? var[i] / ni : INFINITY;
    }
    return v;
}

int alloc_build(const LightInput *lights, int n, long long budget,
                int min_samples, double fallback_var, AllocResult *out) {
    if (!lights || !out || n <= 0 || n > ALLOC_MAX_LIGHTS) return -1;
    if (budget < 0 || min_samples < 1) return -1;

    memset(out, 0, sizeof(*out));
    out->n_lights = n;

    /* Resolve variances first so fallbacks are recorded even on failure. */
    double fb = alloc_fallback_var(lights, n, fallback_var);
    out->fallback_var = fb;
    for (int i = 0; i < n; i++) {
        const LightInput *l = &lights[i];
        if (l->cost <= 0) return -1;
        if (l->var_known && valid_var(l->var)) {
            out->effective_var[i] = l->var;
        } else {
            out->effective_var[i] = fb;
            out->fallback_mask |= 1ULL << i;
            out->fallback_count++;
        }
    }

    /* Feasibility: every light needs min_samples. */
    long long base_cost = 0;
    for (int i = 0; i < n; i++) base_cost += lights[i].cost * (long long)min_samples;
    if (base_cost > budget) {
        out->feasible = 0;
        out->cost_used = base_cost;
        return -1;
    }
    out->feasible = 1;

    long long cost[ALLOC_MAX_LIGHTS];
    for (int i = 0; i < n; i++) cost[i] = lights[i].cost;

    out->theoretical_var = continuous_opt(out->effective_var, cost, n,
                                          min_samples, budget, out->continuous);

    /* Integer stage: floor the continuous solution, respecting feasibility. */
    long long spent = 0;
    for (int i = 0; i < n; i++) {
        long long extra = (long long)floor(out->continuous[i] - (double)min_samples);
        if (extra < 0) extra = 0;
        if (extra > (long long)N_SAMPLE_CAP - min_samples)
            extra = (long long)N_SAMPLE_CAP - min_samples;
        out->samples[i] = min_samples + (int)extra;
        spent += cost[i] * (long long)out->samples[i];
    }
    /* Sequential flooring can overshoot the budget; trim from the lights
     * whose reduction costs the least variance until it fits. */
    while (spent > budget) {
        int best = -1;
        double best_loss = INFINITY;
        for (int i = 0; i < n; i++) {
            if (out->samples[i] <= min_samples) continue;
            double ni = (double)out->samples[i];
            double loss = out->effective_var[i] / ((double)cost[i] * ni * (ni - 1.0));
            if (loss < best_loss) { best_loss = loss; best = i; }
        }
        if (best < 0) break;
        out->samples[best]--;
        spent -= cost[best];
    }

    /* Spend the remaining budget by marginal variance reduction per cost:
     * V(n) - V(n+1) = var_i / (n_i (n_i + 1)). */
    while (spent < budget) {
        int best = -1;
        double best_gain = -1.0;
        for (int i = 0; i < n; i++) {
            if (cost[i] > budget - spent) continue;
            if (out->samples[i] >= N_SAMPLE_CAP) continue;
            double ni = (double)out->samples[i];
            double gain = out->effective_var[i] / ((double)cost[i] * ni * (ni + 1.0));
            if (gain > best_gain) { best_gain = gain; best = i; }
        }
        if (best < 0 || best_gain <= 0.0) break;
        out->samples[best]++;
        spent += cost[best];
    }
    out->cost_used = spent;

    out->allocated_var = alloc_variance(out->effective_var, out->samples, n);
    return 0;
}

double alloc_variance(const double *var, const int *samples, int n) {
    double v = 0.0;
    for (int i = 0; i < n; i++) {
        if (var[i] <= 0.0) continue; /* zero-variance light adds nothing */
        if (samples[i] <= 0) return INFINITY;
        v += var[i] / (double)samples[i];
    }
    return v;
}

/* Exhaustive search over integer allocations, pruning budgets that cannot
 * still cover the remaining lights' minimum samples. */
static int bf_rec(const long long *cost, const double *var, int n,
                  long long rem, int min_samples, int idx,
                  int *cur, int *best, long long *count,
                  double *best_var, long long max_enum) {
    if (idx == n) {
        if (*count >= max_enum) return -1;
        double v = alloc_variance(var, cur, n);
        (*count)++;
        if (v < *best_var) {
            *best_var = v;
            memcpy(best, cur, (size_t)n * sizeof(int));
        }
        return 0;
    }
    long long tail_min = 0;
    for (int j = idx + 1; j < n; j++) tail_min += cost[j] * (long long)min_samples;
    if (cost[idx] * (long long)min_samples + tail_min > rem) return 0;
    long long max_n = (rem - tail_min) / cost[idx];
    if (max_n > N_SAMPLE_CAP) max_n = N_SAMPLE_CAP;
    for (long long k = min_samples; k <= max_n; k++) {
        cur[idx] = (int)k;
        long long before = *count;
        if (bf_rec(cost, var, n, rem - cost[idx] * k, min_samples, idx + 1,
                   cur, best, count, best_var, max_enum) != 0)
            return -1;
        if (*count == before && k > min_samples) break; /* tail infeasible */
    }
    return 0;
}

long long alloc_bruteforce(const LightInput *lights, int n, long long budget,
                           int min_samples, double fallback_var,
                           int *best_samples, double *best_var,
                           long long max_enum) {
    if (!lights || !best_samples || !best_var || n <= 0 || n > ALLOC_MAX_LIGHTS)
        return -1;
    if (budget < 0 || min_samples < 1 || max_enum <= 0) return -1;

    double fb = alloc_fallback_var(lights, n, fallback_var);
    long long cost[ALLOC_MAX_LIGHTS];
    double var[ALLOC_MAX_LIGHTS];
    long long min_cost = 0;
    for (int i = 0; i < n; i++) {
        if (lights[i].cost <= 0) return -1;
        cost[i] = lights[i].cost;
        var[i] = (lights[i].var_known && valid_var(lights[i].var))
                     ? lights[i].var
                     : fb;
        min_cost += cost[i] * (long long)min_samples;
    }
    if (min_cost > budget) return -1;

    int cur[ALLOC_MAX_LIGHTS] = {0};
    *best_var = INFINITY;
    long long count = 0;
    if (bf_rec(cost, var, n, budget, min_samples, 0, cur, best_samples,
               &count, best_var, max_enum) != 0)
        return -1;
    return count;
}

double alloc_deviation_bound(const AllocResult *r) {
    if (!r || r->theoretical_var <= 0.0) return 0.0;
    double d = (r->allocated_var - r->theoretical_var) / r->theoretical_var;
    return d > 0.0 ? d : 0.0;
}
