
/* temporal: variance-optimal history weighting for temporal accumulation.

   Blend: out = (1 - a) * current + a * history, a = history weight.
   With current variance vc and history variance vh,

       E[(out - truth)^2] = (1 - a)^2 * vc + a^2 * vh,

   whose minimizer and minimum are

       a*     = vc / (vc + vh)
       E_min  = vc * vh / (vc + vh).

   If vc + vh == 0 the estimate is exact; fall back to the current frame
   (a = 0) instead of evaluating 0 / 0.

   Fixed-weight recursion with iid current variance s2:

       V_{k+1} = (1 - a)^2 * s2 + a^2 * V_k
       V_inf   = s2 * (1 - a) / (1 + a).

   So a capped weight a <= a_max has the steady-state floor
   s2 * (1 - a_max) / (1 + a_max).  Without a cap the per-frame optimal
   weight telescopes to V_n = s2 / n, i.e. no positive floor: the bound
   only exists because motion handling caps the history weight.

   Ghosting is measured against the current frame as truth:
   |out - current| > tau marks a ghosted pixel.  Under camera motion the
   stale history is biased, so the count rises with motion; an alarm is
   raised once the ghosted fraction exceeds a threshold. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_PIX 96
#define BRUTE_STEPS 2048
#define BRUTE_H (1.0 / (double)BRUTE_STEPS)
#define PI 3.14159265358979323846

typedef struct {
    double cur;       /* current frame sample, treated as truth */
    double hist;      /* reprojected history sample */
    double var_cur;   /* variance of the current sample */
    double var_hist;  /* variance of the history sample */
} TemporalPixel;

typedef struct {
    TemporalPixel px[MAX_PIX];
    int n;
} TemporalFrame;

typedef struct {
    double value;
    double weight;  /* history weight actually used */
    int fallback;   /* 1 when both variances are zero and cur is returned */
} BlendResult;

typedef struct {
    int count;      /* pixels with |out - cur| > tau */
    double fraction;
    int alarm;
} GhostReport;

static double clamp_var(double v) { return v > 0.0 ? v : 0.0; }

/* a* = vc / (vc + vh); zero total variance -> current frame only. */
static double optimal_history_weight(double vc, double vh) {
    vc = clamp_var(vc);
    vh = clamp_var(vh);
    double denom = vc + vh;
    if (!(denom > 0.0)) return 0.0;
    return vc / denom;
}

static double blend_variance(double vc, double vh, double a) {
    vc = clamp_var(vc);
    vh = clamp_var(vh);
    double b = 1.0 - a;
    return b * b * vc + a * a * vh;
}

static double min_variance(double vc, double vh) {
    vc = clamp_var(vc);
    vh = clamp_var(vh);
    double denom = vc + vh;
    if (!(denom > 0.0)) return 0.0;
    return vc * vh / denom;
}

static BlendResult blend_pixel(const TemporalPixel* p) {
    BlendResult r;
    double vc = clamp_var(p->var_cur);
    double vh = clamp_var(p->var_hist);
    if (!(vc + vh > 0.0)) {
        r.value = p->cur;
        r.weight = 0.0;
        r.fallback = 1;
        return r;
    }
    double a = vc / (vc + vh);
    r.value = (1.0 - a) * p->cur + a * p->hist;
    r.weight = a;
    r.fallback = 0;
    return r;
}

/* Exhaustive grid search over the history weight. */
static double brute_min_variance(double vc, double vh, double* best_a) {
    double best = HUGE_VAL;
    double ba = 0.0;
    for (int i = 0; i <= BRUTE_STEPS; i++) {
        double a = (double)i / (double)BRUTE_STEPS;
        double v = blend_variance(vc, vh, a);
        if (v < best) {
            best = v;
            ba = a;
        }
    }
    if (best_a) *best_a = ba;
    return best;
}

/* A grid point lies within h / 2 of a*, and the quadratic has second
   derivative 2 * (vc + vh), so f(grid) - f(a*) <= (vc + vh) * h^2 / 4. */
static double brute_deviation_bound(double vc, double vh) {
    return (clamp_var(vc) + clamp_var(vh)) * BRUTE_H * BRUTE_H / 4.0;
}

static double steady_state_variance(double s2, double a) {
    if (a < 0.0) a = 0.0;
    if (a > 1.0) a = 1.0;
    return s2 * (1.0 - a) / (1.0 + a);
}

static double steady_state_floor(double s2, double a_max) {
    return steady_state_variance(s2, a_max);
}

/* Deterministic band-limited scene: neighboring pixels see different
   gradients, so displacement pushes pixels over tau at different times. */
static double texture(double x) {
    return 0.5
         + 0.18 * sin(2.0 * PI * 2.0 * x)
         + 0.11 * sin(2.0 * PI * 5.0 * x + 0.7)
         + 0.06 * sin(2.0 * PI * 11.0 * x + 2.1);
}

static void make_motion_frame(TemporalFrame* f, int n, double motion) {
    f->n = n;
    for (int i = 0; i < n; i++) {
        double x = ((double)i + 0.5) / (double)n;
        double vc = 0.04 + 0.02 * (0.5 + 0.5 * sin(2.0 * PI * 3.0 * x + 0.3));
        double vh = 0.004 + 0.004 * (0.5 + 0.5 * cos(2.0 * PI * 2.0 * x));
        f->px[i].cur = texture(x);
        f->px[i].hist = texture(x - motion);
        f->px[i].var_cur = vc;
        f->px[i].var_hist = vh;
    }
}

static GhostReport ghost_report(const TemporalFrame* f, double tau, double alarm_frac) {
    GhostReport r;
    r.count = 0;
    r.fraction = 0.0;
    r.alarm = 0;
    for (int i = 0; i < f->n; i++) {
        BlendResult b = blend_pixel(&f->px[i]);
        double err = fabs(b.value - f->px[i].cur);
        if (err > tau) r.count++;
    }
    if (f->n > 0) r.fraction = (double)r.count / (double)f->n;
    r.alarm = (r.fraction > alarm_frac) ? 1 : 0;
    return r;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* Optimal weight: equal variances split, noisy current trusts history,
       an exact history gets all the weight. */
    check("weight_equal_var", fabs(optimal_history_weight(1.0, 1.0) - 0.5) < 1e-12);
    check("weight_noisy_current", fabs(optimal_history_weight(4.0, 1.0) - 0.8) < 1e-12);
    check("weight_exact_history", fabs(optimal_history_weight(1.0, 0.0) - 1.0) < 1e-12);
    check("weight_exact_current", fabs(optimal_history_weight(0.0, 1.0) - 0.0) < 1e-12);
    check("weight_zero_both", fabs(optimal_history_weight(0.0, 0.0)) < 1e-12);

    /* Minimum variance is the harmonic-style combination, never above
       either input variance. */
    check("min_var_harmonic",
          fabs(min_variance(2.0, 3.0) - 2.0 * 3.0 / 5.0) < 1e-12);
    check("min_var_le_inputs",
          min_variance(2.0, 3.0) <= 2.0 + 1e-12 &&
          min_variance(2.0, 3.0) <= 3.0 + 1e-12);
    check("min_var_zero", min_variance(0.0, 0.0) == 0.0);

    /* Brute-force cross-check over a grid of weights: analytic <= brute,
       and the gap is below the rigourous grid bound. */
    static const double pairs[][2] = {
        { 1.0, 1.0 }, { 0.25, 0.75 }, { 2.5, 0.1 }, { 1e-3, 3.0 },
        { 0.0, 1.0 }, { 1.0, 0.0 }, { 0.0, 0.0 }
    };
    int brute_ok = 1, bound_ok = 1, arg_ok = 1;
    double max_dev = 0.0, max_bound = 0.0;
    int npairs = (int)(sizeof(pairs) / sizeof(pairs[0]));
    for (int i = 0; i < npairs; i++) {
        double vc = pairs[i][0], vh = pairs[i][1];
        double a = optimal_history_weight(vc, vh);
        double mn = min_variance(vc, vh);
        double ba = 0.0;
        double bm = brute_min_variance(vc, vh, &ba);
        double dev = bm - mn;
        double bound = brute_deviation_bound(vc, vh);
        if (dev < -1e-12) brute_ok = 0;
        if (dev > bound + 1e-12) bound_ok = 0;
        if (fabs(ba - a) > BRUTE_H / 2.0 + 1e-12) arg_ok = 0;
        if (dev > max_dev) max_dev = dev;
        if (bound > max_bound) max_bound = bound;
    }
    check("brute_ge_analytic", brute_ok);
    check("brute_dev_bound", bound_ok);
    check("brute_argmin_close", arg_ok);
    printf("brute_max_dev=%.3e bound=%.3e\n", max_dev, max_bound);

    /* Fixed-weight recursion converges to its analytic steady state. */
    double s2 = 1.0, a = 0.8;
    double v = 0.0;
    for (int i = 0; i < 2000; i++)
        v = (1.0 - a) * (1.0 - a) * s2 + a * a * v;
    check("steady_state_limit",
          fabs(v - steady_state_variance(s2, a)) < 1e-9);

    /* With a <= a_max, the floor is attained at a_max. */
    double a_max = 0.9;
    double floor = steady_state_floor(s2, a_max);
    static const double alphas[] = { 0.0, 0.2, 0.5, 0.8, 0.9 };
    int floor_ok = 1;
    for (int i = 0; i < 5; i++)
        if (steady_state_variance(s2, alphas[i]) + 1e-12 < floor) floor_ok = 0;
    check("steady_state_floor", floor_ok);

    /* Uncapped per-frame optimal weighting telescopes to s2 / n. */
    int tele_ok = 1;
    v = s2;
    for (int n = 2; n <= 32; n++) {
        double w = optimal_history_weight(s2, v);
        v = blend_variance(s2, v, w);
        if (fabs(v - s2 / (double)n) > 1e-12) tele_ok = 0;
    }
    check("optimal_recursion_telescopes", tele_ok);

    /* Zero variance: fall back to the current frame, no 0/0. */
    TemporalFrame z;
    memset(&z, 0, sizeof(z));
    z.n = 4;
    int zero_ok = 1;
    for (int i = 0; i < z.n; i++) {
        z.px[i].cur = 10.0 + (double)i;
        z.px[i].hist = 100.0 + (double)i;
        z.px[i].var_cur = 0.0;
        z.px[i].var_hist = 0.0;
        BlendResult b = blend_pixel(&z.px[i]);
        if (!b.fallback || b.weight != 0.0 || b.value != z.px[i].cur)
            zero_ok = 0;
    }
    check("zero_var_fallback", zero_ok);

    TemporalPixel exact_hist = { 0.0, 7.0, 1.0, 0.0 };
    BlendResult eh = blend_pixel(&exact_hist);
    check("exact_history_used", !eh.fallback && eh.weight == 1.0 && eh.value == 7.0);

    TemporalPixel exact_cur = { 3.0, 7.0, 0.0, 1.0 };
    BlendResult ec = blend_pixel(&exact_cur);
    check("exact_current_used", !ec.fallback && ec.weight == 0.0 && ec.value == 3.0);

    /* Ghost count must be monotone in motion, and high motion must alarm. */
    static const double motions[] = { 0.0, 0.001, 0.002, 0.004, 0.008, 0.016 };
    int nm = (int)(sizeof(motions) / sizeof(motions[0]));
    int counts[6], alarms[6];
    for (int i = 0; i < nm; i++) {
        TemporalFrame f;
        make_motion_frame(&f, 64, motions[i]);
        GhostReport r = ghost_report(&f, 0.005, 0.5);
        counts[i] = r.count;
        alarms[i] = r.alarm;
    }
    int mono_ok = 1;
    for (int i = 0; i + 1 < nm; i++)
        if (counts[i] > counts[i + 1]) mono_ok = 0;
    check("ghost_zero_motion", counts[0] == 0 && alarms[0] == 0);
    check("ghost_monotone", mono_ok && counts[nm - 1] > counts[0]);
    check("ghost_alarm_rises", alarms[0] == 0 && alarms[nm - 1] == 1);
    printf("ghost_counts=%d,%d,%d,%d,%d,%d\n",
           counts[0], counts[1], counts[2], counts[3], counts[4], counts[5]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
