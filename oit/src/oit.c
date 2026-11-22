
/* oit: order-independent transparency by depth-weighted accumulation.

   A stack of transparent layers is composited exactly ("over", sorted
   front to back) and compared against weighted blending.  Layers are put
   into a canonical order before any arithmetic, so permuting the input
   leaves every output bit unchanged.

   Exact sorted composite:
       T_i     = prod_{j<i} (1 - a_j)          front transmittance
       f_i     = a_i * T_i                     exact front weight
       T_total = prod_i (1 - a_i)
       C_exact = sum_i f_i c_i + T_total c_bg
   with sum_i f_i + T_total = 1.

   Weighted accumulation at extinction kappa > 0:
       w_i  = a_i * exp(-kappa t_i)
       W_bg = exp(-kappa * max_i t_i)
       C_oit = (sum_i w_i c_i + W_bg c_bg) / (sum_i w_i + W_bg)

   Weight / occlusion relation.  Sort the layers by depth.  Since
   exp(-kappa t_1) <= 1 = T_1, and by induction
       exp(-kappa t_i) = exp(-kappa t_{i-1}) exp(-kappa (t_i - t_{i-1}))
                       <= T_{i-1} (1 - a_{i-1}) = T_i
   whenever kappa (t_i - t_{i-1}) >= -ln(1 - a_{i-1}), we get
       w_i <= f_i            for every layer, and
       sum_i w_i <= 1 - T_total.
   So the exponential weight never overstates how much of the front stack
   a layer actually occludes; equality holds exactly when the depth is the
   layer's optical depth, exp(-kappa t_i) = T_i.  oit_kappa_for() returns
   the smallest kappa that certifies the inequality for a given stack.

   Error bound.  Write W_i = w_i / (sum w + W_bg) and B = W_bg / (...).
   Both C_oit and C_exact are convex combinations over the same colours
   (layers plus background), so with Cmax the largest colour component,
       |C_oit - C_exact| <= Cmax/2 * (sum_i |W_i - f_i| + |B - T_total|).
   oit_error_bound() evaluates the right-hand side; the tests check that
   the measured error never exceeds it.

   Zero layers are rejected: an empty stack is a caller error, not a
   background-only image. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct { double depth, alpha, color[3]; } OitLayer;

typedef enum {
    OIT_OK = 0,
    OIT_E_LAYERS = -1,  /* n <= 0: empty stack rejected */
    OIT_E_DEPTH = -2,
    OIT_E_ALPHA = -3,
    OIT_E_KAPPA = -4,
    OIT_E_COLOR = -5,
    OIT_E_ALLOC = -6
} OitErr;

static int cmp_double(double a, double b) {
    return a < b ? -1 : (a > b ? 1 : 0);
}

static int cmp_layer(const void* pa, const void* pb) {
    const OitLayer* a = (const OitLayer*)pa;
    const OitLayer* b = (const OitLayer*)pb;
    int c = cmp_double(a->depth, b->depth);
    if (c) return c;
    c = cmp_double(a->alpha, b->alpha);
    if (c) return c;
    for (int k = 0; k < 3; k++) {
        c = cmp_double(a->color[k], b->color[k]);
        if (c) return c;
    }
    return 0;
}

/* Fold -0.0 to +0.0 so that bitwise-identical layers compare identical. */
static void canonicalize(OitLayer* l) {
    if (l->depth == 0.0) l->depth = 0.0;
    if (l->alpha == 0.0) l->alpha = 0.0;
    for (int c = 0; c < 3; c++) if (l->color[c] == 0.0) l->color[c] = 0.0;
}

/* Validate and return a canonically sorted copy. */
static OitErr oit_prepare(const OitLayer* in, int n, OitLayer** sorted) {
    *sorted = NULL;
    if (in == NULL || n <= 0) return OIT_E_LAYERS;
    for (int i = 0; i < n; i++) {
        const OitLayer* l = &in[i];
        if (!isfinite(l->depth) || l->depth < 0.0) return OIT_E_DEPTH;
        if (!isfinite(l->alpha) || l->alpha < 0.0 || l->alpha > 1.0) return OIT_E_ALPHA;
        for (int c = 0; c < 3; c++)
            if (!isfinite(l->color[c])) return OIT_E_COLOR;
    }
    OitLayer* s = (OitLayer*)malloc((size_t)n * sizeof *s);
    if (s == NULL) return OIT_E_ALLOC;
    for (int i = 0; i < n; i++) { s[i] = in[i]; canonicalize(&s[i]); }
    qsort(s, (size_t)n, sizeof *s, cmp_layer);
    *sorted = s;
    return OIT_OK;
}

/* Exact front weights f_i (canonical order) and total transmittance.
   sum_i f_i + *t_total == 1. */
static OitErr oit_exact_weights(const OitLayer* layers, int n,
                                double* f, double* t_total) {
    OitLayer* s;
    OitErr e = oit_prepare(layers, n, &s);
    if (e != OIT_OK) return e;
    double T = 1.0;
    for (int i = 0; i < n; i++) {
        f[i] = s[i].alpha * T;
        T *= 1.0 - s[i].alpha;
    }
    *t_total = T;
    free(s);
    return OIT_OK;
}

/* Exact "over" composite, layers sorted front to back. */
static OitErr oit_exact_sorted(const OitLayer* layers, int n,
                               const double bg[3], double out[3]) {
    OitLayer* s;
    OitErr e = oit_prepare(layers, n, &s);
    if (e != OIT_OK) return e;
    double T = 1.0, acc[3] = { 0.0, 0.0, 0.0 };
    for (int i = 0; i < n; i++) {
        double f = s[i].alpha * T;
        for (int c = 0; c < 3; c++) acc[c] += f * s[i].color[c];
        T *= 1.0 - s[i].alpha;
    }
    for (int c = 0; c < 3; c++) out[c] = acc[c] + T * bg[c];
    free(s);
    return OIT_OK;
}

/* Depth-weighted accumulation; identical for every permutation of layers. */
static OitErr oit_weighted(const OitLayer* layers, int n, double kappa,
                           const double bg[3], double out[3]) {
    if (!isfinite(kappa) || !(kappa > 0.0)) return OIT_E_KAPPA;
    OitLayer* s;
    OitErr e = oit_prepare(layers, n, &s);
    if (e != OIT_OK) return e;
    double A = 0.0, acc[3] = { 0.0, 0.0, 0.0 };
    double t_max = s[0].depth;
    for (int i = 0; i < n; i++) {
        double w = s[i].alpha * exp(-kappa * s[i].depth);
        A += w;
        for (int c = 0; c < 3; c++) acc[c] += w * s[i].color[c];
        if (s[i].depth > t_max) t_max = s[i].depth;
    }
    double w_bg = exp(-kappa * t_max);
    double denom = A + w_bg;
    if (denom > 0.0) {
        for (int c = 0; c < 3; c++) out[c] = (acc[c] + w_bg * bg[c]) / denom;
    } else {
        for (int c = 0; c < 3; c++) out[c] = bg[c];  /* all weights underflowed */
    }
    free(s);
    return OIT_OK;
}

/* Smallest kappa that certifies w_i <= f_i for every layer:
   kappa (t_i - t_{i-1}) >= -ln(1 - a_{i-1}) for adjacent sorted layers. */
static OitErr oit_kappa_for(const OitLayer* layers, int n, double* kappa_out) {
    OitLayer* s;
    OitErr e = oit_prepare(layers, n, &s);
    if (e != OIT_OK) return e;
    double k = 0.0;
    int have = 0;
    for (int i = 1; i < n; i++) {
        double a = s[i - 1].alpha;
        if (a <= 0.0) continue;                  /* factor 1, no constraint */
        if (a >= 1.0) { free(s); return OIT_E_ALPHA; } /* no finite kappa */
        double gap = s[i].depth - s[i - 1].depth;
        if (!(gap > 0.0)) { free(s); return OIT_E_DEPTH; }
        double req = -log1p(-a) / gap;
        if (!have || req > k) { k = req; have = 1; }
    }
    free(s);
    *kappa_out = have ? k : 1.0;
    return OIT_OK;
}

/* Certified upper bound on the per-channel difference between
   oit_weighted() and oit_exact_sorted(). */
static OitErr oit_error_bound(const OitLayer* layers, int n, double kappa,
                              const double bg[3], double* bound) {
    if (!isfinite(kappa) || !(kappa > 0.0)) return OIT_E_KAPPA;
    OitLayer* s;
    OitErr e = oit_prepare(layers, n, &s);
    if (e != OIT_OK) return e;

    double cmax = 0.0, A = 0.0, t_max = s[0].depth;
    for (int c = 0; c < 3; c++) {
        double v = fabs(bg[c]);
        if (v > cmax) cmax = v;
    }
    for (int i = 0; i < n; i++) {
        for (int c = 0; c < 3; c++) {
            double v = fabs(s[i].color[c]);
            if (v > cmax) cmax = v;
        }
        A += s[i].alpha * exp(-kappa * s[i].depth);
        if (s[i].depth > t_max) t_max = s[i].depth;
    }

    double w_bg = exp(-kappa * t_max);
    double denom = A + w_bg;
    if (!(denom > 0.0)) {
        /* weighted() degrades to the background; the exact composite stays
           inside [-cmax, cmax], so the gap is at most 2*cmax. */
        *bound = 2.0 * cmax;
        free(s);
        return OIT_OK;
    }

    double tv = 0.0, T = 1.0;
    for (int i = 0; i < n; i++) {
        double f = s[i].alpha * T;
        double w = s[i].alpha * exp(-kappa * s[i].depth);
        tv += fabs(w / denom - f);
        T *= 1.0 - s[i].alpha;
    }
    tv += fabs(w_bg / denom - T);
    free(s);
    *bound = 0.5 * cmax * tv;
    return OIT_OK;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static OitLayer layer(double depth, double alpha, double r, double g, double b) {
    OitLayer l;
    l.depth = depth; l.alpha = alpha;
    l.color[0] = r; l.color[1] = g; l.color[2] = b;
    return l;
}

static int same3(const double a[3], const double b[3]) {
    return memcmp(a, b, 3 * sizeof(double)) == 0;
}

static double maxdiff3(const double a[3], const double b[3]) {
    double d = 0.0;
    for (int c = 0; c < 3; c++) {
        double v = fabs(a[c] - b[c]);
        if (v > d) d = v;
    }
    return d;
}

int main(void) {
    const double bg[3] = { 0.1, 0.2, 0.3 };

    /* --- layer count and domain rejection ---------------------------- */
    double out[3] = { 42.0, 42.0, 42.0 };
    double bound = -1.0, T = -1.0, f[8];
    check("reject_zero_exact", oit_exact_sorted(NULL, 0, bg, out) == OIT_E_LAYERS);
    check("reject_zero_weighted", oit_weighted(NULL, 0, 1.0, bg, out) == OIT_E_LAYERS);
    check("reject_zero_front", oit_exact_weights(NULL, 0, f, &T) == OIT_E_LAYERS);
    check("reject_zero_bound", oit_error_bound(NULL, 0, 1.0, bg, &bound) == OIT_E_LAYERS);
    check("reject_zero_kappa_for", oit_kappa_for(NULL, 0, &bound) == OIT_E_LAYERS);
    check("reject_negative_n", oit_weighted(NULL, -4, 1.0, bg, out) == OIT_E_LAYERS);
    check("reject_no_write", out[0] == 42.0 && out[1] == 42.0 && out[2] == 42.0);

    OitLayer good[1] = { layer(1.0, 0.5, 0.4, 0.6, 0.8) };
    check("single_layer_ok", oit_exact_sorted(good, 1, bg, out) == OIT_OK &&
                             oit_weighted(good, 1, 1.0, bg, out) == OIT_OK);

    OitLayer bad_alpha[1] = { layer(1.0, 1.5, 0.0, 0.0, 0.0) };
    check("reject_alpha", oit_exact_sorted(bad_alpha, 1, bg, out) == OIT_E_ALPHA &&
                          oit_weighted(bad_alpha, 1, 1.0, bg, out) == OIT_E_ALPHA);
    OitLayer bad_depth[1] = { layer(-1.0, 0.5, 0.0, 0.0, 0.0) };
    check("reject_depth", oit_exact_sorted(bad_depth, 1, bg, out) == OIT_E_DEPTH);
    OitLayer nan_depth[1] = { layer(NAN, 0.5, 0.0, 0.0, 0.0) };
    check("reject_nan_depth", oit_exact_sorted(nan_depth, 1, bg, out) == OIT_E_DEPTH);
    check("reject_kappa_zero", oit_weighted(good, 1, 0.0, bg, out) == OIT_E_KAPPA);
    check("reject_kappa_neg", oit_weighted(good, 1, -1.0, bg, out) == OIT_E_KAPPA);

    /* --- exact partition of the sorted composite -------------------- */
    OitLayer stack[6] = {
        layer(0.5, 0.30, 0.9, 0.1, 0.1),
        layer(1.5, 0.55, 0.1, 0.9, 0.1),
        layer(2.5, 0.20, 0.1, 0.1, 0.9),
        layer(4.0, 0.70, 0.8, 0.8, 0.2),
        layer(6.0, 0.15, 0.2, 0.5, 0.7),
        layer(9.0, 0.40, 0.6, 0.3, 0.9)
    };
    check("front_partition", oit_exact_weights(stack, 6, f, &T) == OIT_OK &&
                             fabs(f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + T - 1.0) < 1e-12);

    double kappa = 0.0;
    check("kappa_for_ok", oit_kappa_for(stack, 6, &kappa) == OIT_OK && kappa > 0.0);

    /* --- order independence, bit for bit ---------------------------- */
    double w_ref[3], w_perm[3], e_ref[3], e_perm[3];
    oit_weighted(stack, 6, kappa, bg, w_ref);
    oit_exact_sorted(stack, 6, bg, e_ref);

    OitLayer rev[6];
    for (int i = 0; i < 6; i++) rev[i] = stack[5 - i];
    oit_weighted(rev, 6, kappa, bg, w_perm);
    oit_exact_sorted(rev, 6, bg, e_perm);
    check("order_invariant_weighted", same3(w_ref, w_perm));
    check("order_invariant_exact", same3(e_ref, e_perm));

    int idx[6] = { 3, 0, 5, 2, 1, 4 };
    OitLayer shuffled[6];
    for (int i = 0; i < 6; i++) shuffled[i] = stack[idx[i]];
    oit_weighted(shuffled, 6, kappa, bg, w_perm);
    oit_exact_sorted(shuffled, 6, bg, e_perm);
    check("shuffle_invariant_weighted", same3(w_ref, w_perm));
    check("shuffle_invariant_exact", same3(e_ref, e_perm));
    printf("sorted_vs_weighted_err=%.3e\n", maxdiff3(e_ref, w_ref));

    /* --- certified error bound -------------------------------------- */
    double err = maxdiff3(e_ref, w_ref);
    check("bound_ok", oit_error_bound(stack, 6, kappa, bg, &bound) == OIT_OK);
    printf("bound=%.3e\n", bound);
    check("error_within_bound", err <= bound + 1e-12);

    /* An under-estimated kappa is the worst case; the bound must still hold. */
    double bad_kappa = 0.02;
    oit_weighted(stack, 6, bad_kappa, bg, w_perm);
    oit_exact_sorted(stack, 6, bg, e_perm);
    double err_bad = maxdiff3(e_perm, w_perm);
    double bound_bad = 0.0;
    check("bound_ok_bad_kappa", oit_error_bound(stack, 6, bad_kappa, bg, &bound_bad) == OIT_OK);
    printf("bad_kappa_err=%.3e bad_kappa_bound=%.3e\n", err_bad, bound_bad);
    check("bad_kappa_error_nontrivial", err_bad > 1e-3);
    check("bad_kappa_within_bound", err_bad <= bound_bad + 1e-12);

    /* --- weight vs occlusion ---------------------------------------- */
    int all_le = 1;
    for (int i = 0; i < 6; i++) {
        double w = stack[i].alpha * exp(-kappa * stack[i].depth);
        if (w > f[i] + 1e-12) all_le = 0;
    }
    check("weight_le_front", all_le);

    double wsum = 0.0;
    for (int i = 0; i < 6; i++) wsum += stack[i].alpha * exp(-kappa * stack[i].depth);
    check("weight_sum_le_occlusion", wsum <= (1.0 - T) + 1e-12);

    OitLayer even[4] = {
        layer(1.0, 0.3, 1.0, 1.0, 1.0),
        layer(2.0, 0.3, 1.0, 1.0, 1.0),
        layer(3.0, 0.3, 1.0, 1.0, 1.0),
        layer(4.0, 0.3, 1.0, 1.0, 1.0)
    };
    double we[4];
    for (int i = 0; i < 4; i++) we[i] = even[i].alpha * exp(-kappa * even[i].depth);
    check("weight_monotone_with_depth", we[0] > we[1] && we[1] > we[2] && we[2] > we[3]);

    /* --- optical-depth stack: weights reproduce the exact composite -- */
    const double k0 = 0.8;
    const double alphas[4] = { 0.20, 0.35, 0.10, 0.45 };
    OitLayer beer[5];
    double tb = 1.0;
    for (int i = 0; i < 4; i++) {
        beer[i] = layer(-log(tb) / k0, alphas[i], 0.1 * (i + 1), 0.9 - 0.1 * i, 0.2 * i);
        tb *= 1.0 - alphas[i];
    }
    beer[4] = layer(-log(tb) / k0, 0.0, 0.0, 0.0, 0.0);  /* background marker */

    double fb[5], Tb = 0.0;
    check("beer_weights_ok", oit_exact_weights(beer, 5, fb, &Tb) == OIT_OK);
    double beer_w_err = 0.0, beer_wsum = 0.0;
    for (int i = 0; i < 4; i++) {
        double w = alphas[i] * exp(-k0 * beer[i].depth);
        beer_wsum += w;
        double d = fabs(w - fb[i]);
        if (d > beer_w_err) beer_w_err = d;
    }
    check("beer_weights_match_front", beer_w_err < 1e-12);
    check("beer_weights_sum_occlusion", fabs(beer_wsum - (1.0 - Tb)) < 1e-12);

    double kfit = 0.0;
    check("beer_kappa_for", oit_kappa_for(beer, 5, &kfit) == OIT_OK &&
                            kfit >= k0 - 1e-9 && kfit <= k0 + 1e-9);

    double beer_e[3], beer_w[3], beer_bound = -1.0;
    oit_exact_sorted(beer, 5, bg, beer_e);
    oit_weighted(beer, 5, k0, bg, beer_w);
    double beer_err = maxdiff3(beer_e, beer_w);
    oit_error_bound(beer, 5, k0, bg, &beer_bound);
    printf("beer_err=%.3e beer_bound=%.3e\n", beer_err, beer_bound);
    check("beer_exact_match", beer_err < 1e-9);
    check("beer_bound_tight", beer_bound < 1e-9);
    check("beer_within_bound", beer_err <= beer_bound + 1e-12);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
