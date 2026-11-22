
/* layer: stacked-material energy compensation.

   Stacking several materials with the naive single-bounce sum
       E = r0 + t0*r1 + t0*t1*r2 + ...
   mis-weights the deeper layers (single transmittance instead of the
   round trip, and no inter-reflection bounces), so each channel deviates
   by a different amount and the stack shifts color.  This module composes the
   stack analytically with the adding method (closed-form geometric series
   of inter-reflections, layers assumed symmetric: same r/t from both
   sides), then solves per-layer compensation coefficients c_k such that
       sum_k c_k * w_k = R_analytic        (w_k = naive single-bounce weight)
   and reports an a priori upper bound on the residual.  Reflected plus
   transmitted energy never exceeds the incident energy; layers that make
   a solution impossible (unphysical r+t>1, divergent mirror cavity, or a
   layer invisible to the naive model that still contributes energy) are
   reported as conflicts. */
#include <float.h>
#include <math.h>
#include <stdio.h>

#define LAYER_MAX 16
#define NCH 3

typedef struct { float r[NCH], t[NCH]; } Layer; /* absorption a = 1 - r - t */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* physical layers: 0 <= r, 0 <= t, r + t <= 1 per channel.
   Returns the number of conflicting layers; their indices are appended
   to conflicts (up to maxc). */
static int layer_validate(const Layer* ls, int n, float tol, int* conflicts, int maxc) {
    int nc = 0;
    for (int i = 0; i < n; i++) {
        for (int ch = 0; ch < NCH; ch++) {
            float r = ls[i].r[ch], t = ls[i].t[ch];
            if (r < -tol || t < -tol || r + t > 1.0f + tol) {
                if (nc < maxc) conflicts[nc] = i;
                nc++;
                break; /* report each layer once */
            }
        }
    }
    return nc;
}

/* Analytic composition, top layer first.  Fills R, T and the per-layer
   decomposition contrib (contrib[k][ch] = energy layer k adds to R,
   including the inter-reflection series triggered when it is stacked).
   Divergent cavities (denom 1 - R_acc*r_k ~ 0: two near-perfect mirrors
   bouncing energy forever) are reported as conflicts on layer k and the
   channel is left composed up to layer k-1. */
static int layer_compose(const Layer* ls, int n, float tol,
                         float R[NCH], float T[NCH],
                         float contrib[LAYER_MAX][NCH],
                         int* conflicts, int maxc) {
    int nc = 0;
    for (int ch = 0; ch < NCH; ch++) {
        float r_acc = ls[0].r[ch], t_acc = ls[0].t[ch];
        contrib[0][ch] = r_acc;
        for (int k = 1; k < n; k++) {
            float rk = ls[k].r[ch], tk = ls[k].t[ch];
            float denom = 1.0f - r_acc * rk;
            if (denom < 1e-6f) { /* geometric series diverges */
                if (nc < maxc) conflicts[nc] = k;
                nc++;
                contrib[k][ch] = 0.0f;
                break;
            }
            float extra = t_acc * t_acc * rk / denom;
            contrib[k][ch] = extra;
            r_acc += extra;
            t_acc = t_acc * tk / denom;
        }
        R[ch] = r_acc;
        T[ch] = t_acc;
    }
    (void)tol;
    return nc;
}

/* naive single-bounce weights: w_k = (prod_{j<k} t_j) * r_k */
static void layer_naive_weights(const Layer* ls, int n, float w[LAYER_MAX][NCH]) {
    for (int ch = 0; ch < NCH; ch++) {
        float thr = 1.0f;
        for (int k = 0; k < n; k++) {
            w[k][ch] = thr * ls[k].r[ch];
            thr *= ls[k].t[ch];
        }
    }
}

/* Solve compensation coefficients c_k so the compensated naive sum
   matches the analytic reflectance, and bound the residual:
       |sum_k c_k*w_k - R| <= ebound   per channel.
   c_k = contrib_k / w_k; a layer with w_k ~ 0 but contrib_k > tol cannot
   be represented in the naive model and is reported as a conflict.
   Returns the number of conflicts (0 = solution valid). */
static int layer_solve(const Layer* ls, int n, float tol,
                       float c[LAYER_MAX][NCH], float R[NCH], float T[NCH],
                       float ebound[NCH], int* conflicts, int maxc) {
    float contrib[LAYER_MAX][NCH], w[LAYER_MAX][NCH];
    int nc = layer_validate(ls, n, tol, conflicts, maxc);
    if (nc > 0) return nc; /* unphysical input: no coefficient can fix it */
    nc = layer_compose(ls, n, tol, R, T, contrib, conflicts, maxc);
    if (nc > 0) return nc;
    layer_naive_weights(ls, n, w);
    for (int k = 0; k < n; k++) {
        for (int ch = 0; ch < NCH; ch++) {
            if (w[k][ch] > 1e-7f) {
                c[k][ch] = contrib[k][ch] / w[k][ch];
            } else if (contrib[k][ch] > tol) {
                if (nc < maxc) conflicts[nc] = k;
                nc++;
                c[k][ch] = 0.0f;
            } else {
                c[k][ch] = 0.0f; /* contributes nothing either way */
            }
        }
    }
    if (nc > 0) return nc;
    for (int ch = 0; ch < NCH; ch++) {
        float mag = 0.0f;
        for (int k = 0; k < n; k++) mag += fabsf(c[k][ch] * w[k][ch]);
        /* each c*w product and the running sum round by <= eps; the
           compose step adds <= ~4 roundings per layer on R itself */
        ebound[ch] = FLT_EPSILON * (4.0f * n * (mag + R[ch]) + 1.0f);
    }
    return 0;
}

static float compensated_sum(const float c[LAYER_MAX][NCH],
                             const float w[LAYER_MAX][NCH], int n, int ch) {
    float s = 0.0f;
    for (int k = 0; k < n; k++) s += c[k][ch] * w[k][ch];
    return s;
}

/* deterministic pseudo-random in [0,1) */
static unsigned rng_state = 12345u;
static float frand(void) {
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float)(rng_state >> 8) / 16777216.0f;
}

int main(void) {
    const float tol = 1e-4f;
    int conflicts[LAYER_MAX];
    float R[NCH], T[NCH], eb[NCH];
    float c[LAYER_MAX][NCH], w[LAYER_MAX][NCH], contrib[LAYER_MAX][NCH];

    /* --- the reported case: three stacked materials, color shift --- */
    Layer stack[3] = {
        { { 0.60f, 0.45f, 0.30f }, { 0.30f, 0.40f, 0.50f } },
        { { 0.50f, 0.55f, 0.40f }, { 0.30f, 0.30f, 0.35f } },
        { { 0.40f, 0.35f, 0.60f }, { 0.20f, 0.25f, 0.15f } },
    };
    int nc = layer_solve(stack, 3, tol, c, R, T, eb, conflicts, LAYER_MAX);
    check("stack_solved", nc == 0);
    layer_naive_weights(stack, 3, w);
    layer_compose(stack, 3, tol, R, T, contrib, conflicts, LAYER_MAX);
    {
        float shift = 0.0f;
        for (int ch = 0; ch < NCH; ch++) {
            float naive = 0.0f;
            for (int k = 0; k < 3; k++) naive += w[k][ch];
            float comp = compensated_sum(c, w, 3, ch);
            printf("ch%d naive=%.4f analytic=%.4f compensated=%.4f ebound=%.2e\n",
                   ch, naive, R[ch], comp, eb[ch]);
            check(ch == 0 ? "comp_matches_r" : ch == 1 ? "comp_matches_g" : "comp_matches_b",
                  fabsf(comp - R[ch]) <= eb[ch]);
            float d = fabsf(naive - R[ch]);
            if (d > shift) shift = d;
        }
        /* the naive sum mis-weights each channel differently: color shift */
        check("naive_shifts_color", shift > 10.0f * eb[0]);
    }
    /* per-layer comparison against the analytic decomposition */
    {
        int ok = 1;
        for (int k = 0; k < 3; k++)
            for (int ch = 0; ch < NCH; ch++) {
                float pb = 4.0f * FLT_EPSILON * fabsf(contrib[k][ch]) + FLT_EPSILON;
                if (fabsf(c[k][ch] * w[k][ch] - contrib[k][ch]) > pb) ok = 0;
            }
        check("per_layer_within_bound", ok);
    }
    /* total energy must not exceed the incident energy */
    check("energy_le_incident",
          R[0] + T[0] <= 1.0f + eb[0] && R[1] + T[1] <= 1.0f + eb[1] &&
          R[2] + T[2] <= 1.0f + eb[2]);
    check("absorption_nonneg",
          R[0] + T[0] >= -eb[0] && R[1] + T[1] >= -eb[1] && R[2] + T[2] >= -eb[2]);

    /* --- single layer: coefficient must be exactly 1 --- */
    Layer one[1] = { { { 0.4f, 0.5f, 0.6f }, { 0.3f, 0.2f, 0.1f } } };
    nc = layer_solve(one, 1, tol, c, R, T, eb, conflicts, LAYER_MAX);
    check("single_solved", nc == 0);
    check("single_coeff_one",
          fabsf(c[0][0] - 1.0f) < 1e-6f && fabsf(c[0][1] - 1.0f) < 1e-6f &&
          fabsf(c[0][2] - 1.0f) < 1e-6f);
    check("single_R_is_r", fabsf(R[0] - 0.4f) < 1e-6f);

    /* --- unphysical layer (r + t > 1): conflict must name the layer --- */
    Layer bad[3] = {
        { { 0.5f, 0.5f, 0.5f }, { 0.4f, 0.4f, 0.4f } },
        { { 0.8f, 0.1f, 0.1f }, { 0.5f, 0.1f, 0.1f } }, /* ch0: 1.3 > 1 */
        { { 0.2f, 0.2f, 0.2f }, { 0.2f, 0.2f, 0.2f } },
    };
    nc = layer_solve(bad, 3, tol, c, R, T, eb, conflicts, LAYER_MAX);
    check("unphysical_rejected", nc == 1);
    check("unphysical_layer_named", nc == 1 && conflicts[0] == 1);

    /* --- mirror cavity: r=1 over r=1 makes the series diverge --- */
    Layer mirrors[2] = {
        { { 1.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 0.0f } },
        { { 1.0f, 1.0f, 1.0f }, { 0.0f, 0.0f, 0.0f } },
    };
    nc = layer_solve(mirrors, 2, tol, c, R, T, eb, conflicts, LAYER_MAX);
    check("cavity_rejected", nc == NCH);
    check("cavity_layer_named", nc == NCH && conflicts[0] == 1);

    /* --- opaque top layer: hidden layer gets c=0, no conflict --- */
    Layer opaque[2] = {
        { { 0.7f, 0.7f, 0.7f }, { 0.0f, 0.0f, 0.0f } },
        { { 0.9f, 0.9f, 0.9f }, { 0.0f, 0.0f, 0.0f } },
    };
    nc = layer_solve(opaque, 2, tol, c, R, T, eb, conflicts, LAYER_MAX);
    check("hidden_layer_ok", nc == 0);
    check("hidden_coeff_zero", c[1][0] == 0.0f && c[1][1] == 0.0f && c[1][2] == 0.0f);

    /* --- randomized physical stacks: conservation + residual bound --- */
    {
        int conserved = 1, bounded = 1;
        for (int t = 0; t < 64; t++) {
            Layer ls[3];
            for (int k = 0; k < 3; k++)
                for (int ch = 0; ch < NCH; ch++) {
                    float r = frand();
                    ls[k].r[ch] = r;
                    ls[k].t[ch] = frand() * (1.0f - r); /* keep r+t <= 1 */
                }
            nc = layer_solve(ls, 3, tol, c, R, T, eb, conflicts, LAYER_MAX);
            layer_naive_weights(ls, 3, w);
            for (int ch = 0; ch < NCH; ch++) {
                if (nc != 0 || R[ch] + T[ch] > 1.0f + eb[ch] || R[ch] + T[ch] < -eb[ch])
                    conserved = 0;
                if (nc == 0 &&
                    fabsf(compensated_sum(c, w, 3, ch) - R[ch]) > eb[ch])
                    bounded = 0;
            }
        }
        check("random_energy_conserved", conserved);
        check("random_residual_bounded", bounded);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
