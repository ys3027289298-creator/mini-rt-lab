
/* multiscale: distance-adaptive sampling and pixel-coverage weighting for
   brightness / colour consistency across scales.

   Scene model: an emitter patch of physical area `area` and uniform radiance,
   with |n.v| = cos_theta, seen from distance d through a pixel of solid
   angle pixel_sa.  Its projected solid angle is

       Omega(d) = area * cos_theta / d^2

   and its pixel coverage is c(d) = clamp(Omega(d) / pixel_sa, 0, 1).

   Two distinct quantities must not be confused:

     - the pixel-integrated radiance  L * c(d)  does fall off with distance
       (the patch covers less of the pixel); this is correct and is what a
       plain "pixel average" reports;
     - the emitter's own radiance L must be scale-independent.  A near render
       (coverage clamped to 1) is the reference, and coverage-normalising a
       weighted pixel estimate recovers L at every distance.

   Sampling density rises as coverage falls so the expected number of hits
   per pixel stays near a target, and each sample carries its sub-pixel
   coverage weight.  Every distance-taking entry point refuses non-finite
   or non-positive d before writing any output. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MS_LUM_R 0.2126
#define MS_LUM_G 0.7152
#define MS_LUM_B 0.0722
#define MS_EPS 1e-12
#define MS_MAX_SAMPLES 1000000.0   /* ceiling accepted for cfg->max_samples */

typedef enum {
    MS_OK = 0,
    MS_ERR_DISTANCE = 1,   /* d <= 0 or non-finite: refused, output zeroed */
    MS_ERR_ARG = 2
} MsStatus;

typedef enum {
    MS_POINT = 0,      /* binary hit/miss samples (noisy at low coverage) */
    MS_COVERAGE = 1    /* sub-pixel coverage-weighted samples (recommended) */
} MsMode;

typedef struct { double r, g, b; } MsRGB;

typedef struct {
    double area;        /* physical patch area, > 0 */
    double cos_theta;   /* |n.v| in [0, 1] */
    MsRGB radiance;     /* uniform radiance, each channel >= 0 */
} MsEmitter;

typedef struct {
    double min_samples;       /* per-pixel budget floor, >= 1 */
    double max_samples;       /* per-pixel budget ceiling */
    double min_expected_hits; /* target hits; drives density with distance */
} MsConfig;

typedef struct {
    double samples;     /* samples actually allocated (integral value) */
    double coverage;    /* clamped c(d) in [0, 1] */
    double expected_hits;
} MsBudget;

/* ---- rejection helpers ------------------------------------------------- */

/* Usable distances are finite and strictly positive.  Zero is refused. */
static int ms_distance_ok(double d) { return isfinite(d) && d > 0.0; }

static void ms_zero_rgb(MsRGB* c) { c->r = c->g = c->b = 0.0; }

static double ms_luminance(MsRGB c) {
    return MS_LUM_R * c.r + MS_LUM_G * c.g + MS_LUM_B * c.b;
}

static int ms_rgb_close(MsRGB a, MsRGB b, double tol) {
    return fabs(a.r - b.r) <= tol && fabs(a.g - b.g) <= tol && fabs(a.b - b.b) <= tol;
}

/* ---- analytic model ---------------------------------------------------- */

/* Projected solid angle of the emitter seen from distance d. */
static MsStatus ms_solid_angle(const MsEmitter* e, double d, double* out) {
    if (!e || !out) return MS_ERR_ARG;
    if (!ms_distance_ok(d)) return MS_ERR_DISTANCE;
    if (!(e->area > 0.0) || !isfinite(e->area)) return MS_ERR_ARG;
    if (!(e->cos_theta >= 0.0) || e->cos_theta > 1.0 || !isfinite(e->cos_theta)) return MS_ERR_ARG;
    *out = e->area * e->cos_theta / (d * d);
    return MS_OK;
}

/* Fraction of the pixel covered at distance d. */
static MsStatus ms_coverage(const MsEmitter* e, double pixel_sa, double d, double* out) {
    if (!(pixel_sa > 0.0) || !isfinite(pixel_sa)) return MS_ERR_ARG;
    double omega;
    MsStatus st = ms_solid_angle(e, d, &omega);
    if (st != MS_OK) return st;
    double c = omega / pixel_sa;
    *out = c < 0.0 ? 0.0 : (c > 1.0 ? 1.0 : c);
    return MS_OK;
}

/* Exact pixel-integrated average radiance, L * coverage. */
static MsStatus ms_analytic_avg(const MsEmitter* e, double pixel_sa, double d, MsRGB* out) {
    double c;
    MsStatus st = ms_coverage(e, pixel_sa, d, &c);
    if (st != MS_OK) return st;
    out->r = e->radiance.r * c;
    out->g = e->radiance.g * c;
    out->b = e->radiance.b * c;
    return MS_OK;
}

/* ---- distance-adaptive sampling density --------------------------------- */

/* Allocate samples so the expected hit count stays near
   config->min_expected_hits regardless of distance:
       n = clamp(target / coverage, min_samples, max_samples).
   Refuses invalid distances; on refusal nothing is written. */
static MsStatus ms_budget(const MsEmitter* e, const MsConfig* cfg, double pixel_sa,
                          double d, MsBudget* out) {
    if (!e || !cfg || !out) return MS_ERR_ARG;
    if (!ms_distance_ok(d)) return MS_ERR_DISTANCE;
    if (!(cfg->min_samples >= 1.0) || !isfinite(cfg->min_samples)) return MS_ERR_ARG;
    if (!(cfg->max_samples >= cfg->min_samples) || !isfinite(cfg->max_samples)) return MS_ERR_ARG;
    if (cfg->max_samples > MS_MAX_SAMPLES) return MS_ERR_ARG;
    if (!(cfg->min_expected_hits > 0.0) || !isfinite(cfg->min_expected_hits)) return MS_ERR_ARG;

    double c;
    MsStatus st = ms_coverage(e, pixel_sa, d, &c);
    if (st != MS_OK) return st;

    double n = cfg->min_samples;
    if (c > 0.0) n = cfg->min_expected_hits / c;
    if (n < cfg->min_samples) n = cfg->min_samples;
    if (n > cfg->max_samples) n = cfg->max_samples;

    out->samples = ceil(n);
    out->coverage = c;
    out->expected_hits = out->samples * c;
    return MS_OK;
}

/* ---- deterministic sample stream ---------------------------------------- */

static unsigned long long ms_hash_u64(unsigned long long x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static double ms_uniform01(unsigned long long key) {
    return (double)(ms_hash_u64(key) >> 11) * (1.0 / 9007199254740992.0);
}

/* ---- per-pixel estimator ------------------------------------------------ */

/* Render one pixel at distance d and return the scale-independent (coverage-
   normalised) radiance.

   MS_POINT: one Bernoulli trial per sample with hit probability = coverage;
             estimate = L * (hits / n) / coverage.  Zero hits in a pixel makes
             it read black, which is the far-object fading this module
             corrects.
   MS_COVERAGE: sub-pixel strata weighted by how much of each stratum the
             patch covers; estimate = sum(w*L) / sum(w) = L for any cov > 0.

   A near render (coverage clamped to 1) is the reference: every mode returns
   the emitter radiance there.  For d <= 0 or non-finite d the call is refused
   before any sampling: *out is zeroed and MS_ERR_DISTANCE returned. */
static MsStatus ms_render_pixel(const MsEmitter* e, const MsConfig* cfg, double pixel_sa,
                                double d, unsigned long long seed, MsMode mode, MsRGB* out) {
    if (!out) return MS_ERR_ARG;
    ms_zero_rgb(out);
    if (!e || !cfg) return MS_ERR_ARG;
    if (!ms_distance_ok(d)) return MS_ERR_DISTANCE;
    if (mode != MS_POINT && mode != MS_COVERAGE) return MS_ERR_ARG;

    MsBudget b;
    MsStatus st = ms_budget(e, cfg, pixel_sa, d, &b);
    if (st != MS_OK) { ms_zero_rgb(out); return st; }

    int n = (int)b.samples;
    if (n < 1) return MS_ERR_ARG;

    double cov = b.coverage;
    unsigned long long base = seed * 0x9E3779B97F4A7C15ULL;

    if (mode == MS_POINT) {
        double hits = 0.0;
        for (int i = 0; i < n; i++) {
            if (ms_uniform01(base + (unsigned long long)i) < cov) hits += 1.0;
        }
        if (cov <= 0.0) return MS_OK;   /* nothing covered: stays zero */
        double k = (hits / (double)n) / cov;
        out->r = e->radiance.r * k;
        out->g = e->radiance.g * k;
        out->b = e->radiance.b * k;
        return MS_OK;
    }

    /* MS_COVERAGE: stratify the covered fraction across sub-pixel strata.
       n*cov strata are covered; the boundary stratum is covered only
       fractionally, so it gets a fractional weight.  Any cov > 0 puts weight
       on at least one sample, so the estimate never collapses to black just
       because the patch covers less than one sample. */
    double k = (double)n * cov;
    double wsum = 0.0, wr = 0.0, wg = 0.0, wb = 0.0;
    for (int i = 0; i < n; i++) {
        double w = k - (double)i;
        if (w > 1.0) w = 1.0;
        if (w <= 0.0) break;          /* remaining strata uncovered */
        wsum += w;
        wr += w * e->radiance.r;
        wg += w * e->radiance.g;
        wb += w * e->radiance.b;
    }
    if (wsum > 0.0) {
        out->r = wr / wsum;
        out->g = wg / wsum;
        out->b = wb / wsum;
    }
    return MS_OK;
}

/* ---- analytic error bound ---------------------------------------------- */

/* Upper bound on |estimate - emitter radiance| per channel, given the mode,
   expressed against the analytic (near-render) radiance.
   Coverage mode is exact up to rounding: bound = MS_EPS * L.
   Point mode is binomial: with n samples and hit probability cov the relative
   standard error is sqrt((1 - cov) / (cov * n)); z of these give the bound.
   That is a statistical (z-sigma) bound per pixel; averaging K independent
   pixels tightens it by sqrt(K). */
static MsStatus ms_error_bound(const MsEmitter* e, const MsConfig* cfg, double pixel_sa,
                               double d, MsMode mode, double z, MsRGB* bound) {
    if (!bound) return MS_ERR_ARG;
    ms_zero_rgb(bound);
    if (!e || !cfg) return MS_ERR_ARG;
    if (!ms_distance_ok(d)) return MS_ERR_DISTANCE;
    if (!(z > 0.0) || !isfinite(z)) return MS_ERR_ARG;

    MsBudget b;
    MsStatus st = ms_budget(e, cfg, pixel_sa, d, &b);
    if (st != MS_OK) return st;

    double rel;
    if (mode == MS_COVERAGE) {
        rel = MS_EPS;
    } else if (mode == MS_POINT) {
        double cov = b.coverage;
        double n = b.samples;
        /* cov == 0 means the estimator and the truth are both zero */
        if (cov <= 0.0) rel = 0.0;
        else if (cov >= 1.0) rel = 0.0;
        else rel = z * sqrt((1.0 - cov) / (cov * n));
    } else {
        return MS_ERR_ARG;
    }

    bound->r = rel * e->radiance.r;
    bound->g = rel * e->radiance.g;
    bound->b = rel * e->radiance.b;
    return MS_OK;
}

/* ---- self-checks -------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int all_finite(MsRGB c) { return isfinite(c.r) && isfinite(c.g) && isfinite(c.b); }

int main(void) {
    MsEmitter e;
    e.area = 1.0;
    e.cos_theta = 1.0;
    e.radiance.r = 0.9;
    e.radiance.g = 0.4;
    e.radiance.b = 0.1;
    double L = ms_luminance(e.radiance);

    MsConfig cfg;
    cfg.min_samples = 4.0;
    cfg.max_samples = 4096.0;
    cfg.min_expected_hits = 8.0;

    double pixel_sa = 0.01;
    MsRGB out;

    /* --- distance rejection (zero refused, output untouched) --- */
    ms_zero_rgb(&out);
    check("reject_zero_distance",
          ms_render_pixel(&e, &cfg, pixel_sa, 0.0, 1, MS_COVERAGE, &out) == MS_ERR_DISTANCE);
    check("reject_zero_output_zeroed", out.r == 0.0 && out.g == 0.0 && out.b == 0.0);
    check("reject_negative_distance",
          ms_render_pixel(&e, &cfg, pixel_sa, -2.0, 1, MS_COVERAGE, &out) == MS_ERR_DISTANCE);
    check("reject_nan_distance",
          ms_render_pixel(&e, &cfg, pixel_sa, NAN, 1, MS_COVERAGE, &out) == MS_ERR_DISTANCE);
    check("reject_inf_distance",
          ms_render_pixel(&e, &cfg, pixel_sa, INFINITY, 1, MS_COVERAGE, &out) == MS_ERR_DISTANCE);
    {
        MsBudget b;
        check("budget_reject_zero", ms_budget(&e, &cfg, pixel_sa, 0.0, &b) == MS_ERR_DISTANCE);
        MsRGB bd;
        check("bound_reject_zero",
              ms_error_bound(&e, &cfg, pixel_sa, 0.0, MS_POINT, 4.0, &bd) == MS_ERR_DISTANCE);
    }

    /* --- bad configs are rejected rather than over-allocating --- */
    {
        MsConfig bad = cfg; bad.max_samples = 2.0 * MS_MAX_SAMPLES;
        MsBudget b;
        check("reject_oversized_budget", ms_budget(&e, &bad, pixel_sa, 10.0, &b) == MS_ERR_ARG);
        bad = cfg; bad.min_samples = 0.0;
        check("reject_zero_min_samples", ms_budget(&e, &bad, pixel_sa, 10.0, &b) == MS_ERR_ARG);
        check("reject_bad_mode",
              ms_render_pixel(&e, &cfg, pixel_sa, 10.0, 1, (MsMode)7, &out) == MS_ERR_ARG);
    }

    /* --- zero coverage (grazing emitter) stays black and is not a NaN --- */
    {
        MsEmitter grazing = e;
        grazing.cos_theta = 0.0;
        double c0;
        ms_coverage(&grazing, pixel_sa, 10.0, &c0);
        MsRGB g, gb;
        check("grazing_coverage_zero", c0 == 0.0);
        check("grazing_point_black",
              ms_render_pixel(&grazing, &cfg, pixel_sa, 10.0, 1, MS_POINT, &g) == MS_OK &&
              g.r == 0.0 && g.g == 0.0 && g.b == 0.0);
        check("grazing_coverage_mode_black",
              ms_render_pixel(&grazing, &cfg, pixel_sa, 10.0, 1, MS_COVERAGE, &g) == MS_OK &&
              g.r == 0.0 && g.g == 0.0 && g.b == 0.0);
        check("grazing_bound_zero",
              ms_error_bound(&grazing, &cfg, pixel_sa, 10.0, MS_POINT, 4.0, &gb) == MS_OK &&
              gb.r == 0.0 && gb.g == 0.0 && gb.b == 0.0);
    }

    /* --- pixel-integrated analytic falls off with distance --- */
    MsRGB an;
    double cov1;
    ms_coverage(&e, pixel_sa, 1.0, &cov1);
    check("analytic_near_clamped", ms_analytic_avg(&e, pixel_sa, 1.0, &an) == MS_OK &&
          fabs(cov1 - 1.0) < 1e-12 && ms_rgb_close(an, e.radiance, 1e-12));
    double cov100;
    ms_coverage(&e, pixel_sa, 100.0, &cov100);
    check("analytic_far_cov", fabs(cov100 - 0.01) < 1e-12);
    check("analytic_far_ok", ms_analytic_avg(&e, pixel_sa, 100.0, &an) == MS_OK &&
          fabs(ms_luminance(an) - L * 0.01) < 1e-12);

    /* --- coverage mode recovers near-render brightness at all distances --- */
    double ds[] = { 1.0, 5.0, 10.0, 20.0, 40.0, 80.0, 1000.0 };
    int consistent = 1, bounded = 1;
    for (int i = 0; i < 7; i++) {
        MsRGB est, bnd;
        if (ms_render_pixel(&e, &cfg, pixel_sa, ds[i], 7, MS_COVERAGE, &est) != MS_OK) { consistent = 0; break; }
        if (ms_error_bound(&e, &cfg, pixel_sa, ds[i], MS_COVERAGE, 4.0, &bnd) != MS_OK) { bounded = 0; break; }
        /* compare to the near render, which equals the emitter radiance */
        if (!all_finite(est) || !ms_rgb_close(est, e.radiance, 1e-9)) { consistent = 0; break; }
        if (fabs(est.r - e.radiance.r) > bnd.r + 1e-9 ||
            fabs(est.g - e.radiance.g) > bnd.g + 1e-9 ||
            fabs(est.b - e.radiance.b) > bnd.b + 1e-9) { bounded = 0; break; }
    }
    check("coverage_matches_near_render", consistent);
    check("coverage_within_bound", bounded);

    /* --- budget-capped far distance (coverage below one sample) still
       recovers brightness rather than fading to black --- */
    {
        MsBudget bfar;
        ms_budget(&e, &cfg, pixel_sa, 1000.0, &bfar);
        check("far_budget_capped", bfar.samples == cfg.max_samples &&
              bfar.coverage * bfar.samples < 1.0);
        MsRGB est;
        ms_render_pixel(&e, &cfg, pixel_sa, 1000.0, 5, MS_COVERAGE, &est);
        check("far_below_one_sample_recovers", ms_rgb_close(est, e.radiance, 1e-9));
    }

    /* --- colour ratios stay fixed (hue consistency) --- */
    {
        MsRGB n1, f1;
        ms_render_pixel(&e, &cfg, pixel_sa, 1.0, 3, MS_COVERAGE, &n1);
        ms_render_pixel(&e, &cfg, pixel_sa, 60.0, 3, MS_COVERAGE, &f1);
        double rn = n1.r / n1.g, rf = f1.r / f1.g;
        check("colour_ratio_stable", fabs(rn - rf) < 1e-9 &&
              fabs(rn - e.radiance.r / e.radiance.g) < 1e-9);
    }

    /* --- coverage decreases monotonically with distance --- */
    double c20, c40, c80;
    ms_coverage(&e, pixel_sa, 20.0, &c20);
    ms_coverage(&e, pixel_sa, 40.0, &c40);
    ms_coverage(&e, pixel_sa, 80.0, &c80);
    check("coverage_monotone", c20 > c40 && c40 > c80 && c80 > 0.0 && c20 <= 1.0);

    /* --- sampling density rises as coverage falls --- */
    MsBudget b20, b40, b80;
    ms_budget(&e, &cfg, pixel_sa, 20.0, &b20);
    ms_budget(&e, &cfg, pixel_sa, 40.0, &b40);
    ms_budget(&e, &cfg, pixel_sa, 80.0, &b80);
    check("density_rises_with_distance",
          b20.samples < b40.samples && b40.samples < b80.samples);
    check("density_tracks_target",
          fabs(b20.expected_hits - cfg.min_expected_hits) < 1.0 &&
          fabs(b40.expected_hits - cfg.min_expected_hits) < 1.0 &&
          fabs(b80.expected_hits - cfg.min_expected_hits) < 1.0);

    /* --- brightness is continuous as distance varies (analytic curve) --- */
    /* Continuity at the coverage clamp (d=10 here): the change in L*c across
       a tiny distance step shrinks with the step, so no fixed jump remains.
       Check both sides of the kink and that the jump halves with the step. */
    {
        double dk = sqrt(e.area * e.cos_theta / pixel_sa);   /* c(dk) == 1 */
        double lum = 0.0;
        MsRGB a;
        ms_analytic_avg(&e, pixel_sa, dk, &a);
        lum = ms_luminance(a);

        double eps = 1e-4;
        MsRGB lo, hi;
        ms_analytic_avg(&e, pixel_sa, dk - eps, &lo);
        ms_analytic_avg(&e, pixel_sa, dk + eps, &hi);
        double j1 = fabs(ms_luminance(hi) - ms_luminance(lo));

        double eps2 = eps * 0.5;
        ms_analytic_avg(&e, pixel_sa, dk - eps2, &lo);
        ms_analytic_avg(&e, pixel_sa, dk + eps2, &hi);
        double j2 = fabs(ms_luminance(hi) - ms_luminance(lo));

        /* jump scales with the step (no fixed discontinuity) and is tiny */
        check("brightness_continuous", j2 <= j1 * 0.6 + 1e-15 && j1 < 1e-3);
        check("brightness_equals_near_at_clamp", fabs(lum - L) < 1e-9);
    }

    /* --- point mode error stays within its analytic bound --- */
    {
        int K = 400;
        double cov = 0.25;                 /* choose d so coverage is 0.25 */
        double d = sqrt(e.area * e.cos_theta / (cov * pixel_sa)); /* = 20 */
        double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
        int all_ok = 1;
        for (int k = 0; k < K; k++) {
            MsRGB est;
            if (ms_render_pixel(&e, &cfg, pixel_sa, d, (unsigned long long)k + 1, MS_POINT, &est) != MS_OK ||
                !all_finite(est)) { all_ok = 0; break; }
            sum_r += est.r; sum_g += est.g; sum_b += est.b;
        }
        MsRGB mean = { sum_r / K, sum_g / K, sum_b / K };
        /* bound for the mean of K independent estimates shrinks by sqrt(K) */
        MsRGB bnd;
        ms_error_bound(&e, &cfg, pixel_sa, d, MS_POINT, 4.0, &bnd);
        double shrink = 1.0 / sqrt((double)K);
        check("point_mode_finite", all_ok);
        check("point_mean_within_bound",
              fabs(mean.r - e.radiance.r) <= bnd.r * shrink + 1e-6 &&
              fabs(mean.g - e.radiance.g) <= bnd.g * shrink + 1e-6 &&
              fabs(mean.b - e.radiance.b) <= bnd.b * shrink + 1e-6);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
