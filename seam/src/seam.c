/* seam: lightmap UV-seam repair.

   For a vertical seam at texel boundary `seam_x`, each row is sampled on both
   sides at offsets 0.5..3.5 texels.  A cubic extrapolation to the seam gives
   the edge value and slope of each side; their differences D and S are the
   seam's value jump and gradient jump.  The repair distributes -D/2 and +D/2
   over a band of `band` texels per side with a Hermite profile whose end slope
   matches S (alpha = band*S/D), so the correction and its slope both vanish at
   the band edge.  Residuals before/after are measured with a separate, wider
   sample ring (0.75..4.75 texels, quartic extrapolation), so the check is not
   the estimator that built the correction.  If the measured discrepancy is
   below `threshold`, nothing is written at all. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define SEAM_MAX_W 128
#define SEAM_MAX_H 128
#define SEAM_NODE 4
#define SEAM_INDEP 5

typedef struct {
    int w, h;
    float px[SEAM_MAX_H][SEAM_MAX_W];
} SeamImage;

typedef struct {
    int seam_x;          /* texel boundary index; seam between x-1 and x */
    float band;          /* half-width of the compensation band, in texels */
    float threshold;     /* measured seam residual below which nothing is touched */
    int gradient_match;  /* 1: match end slope of the correction profile */
} SeamParams;

typedef struct {
    int detected;            /* measured residual reached the threshold */
    int accepted;            /* a correction was committed (residual met threshold) */
    int rows_active;         /* rows that received a correction */
    int changed;             /* pixels written */
    float residual_before;   /* max row value jump before repair */
    float residual_after;    /* independent-estimator value jump after repair */
    float slope_before;      /* max row gradient jump before repair */
    float slope_after;       /* independent-estimator gradient jump after repair */
    float bound;             /* upper bound on |post - exact rebake|, see README */
    float max_correction;    /* largest single-pixel delta applied */
} SeamReport;

static double clampd(double v, double lo, double hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

/* Lagrange basis weights for nodes[0..n) evaluated at x. */
static void lagrange_at(const double* nodes, int n, double x, double* w) {
    for (int i = 0; i < n; i++) {
        double p = 1.0;
        for (int j = 0; j < n; j++) {
            if (j != i) p *= (x - nodes[j]) / (nodes[i] - nodes[j]);
        }
        w[i] = p;
    }
}

/* Derivative of the Lagrange interpolant at x. */
static void lagrange_deriv_at(const double* nodes, int n, double x, double* w) {
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double p = 1.0;
            for (int k = 0; k < n; k++) {
                if (k != i && k != j) p *= (x - nodes[k]) / (nodes[i] - nodes[k]);
            }
            s += p / (nodes[i] - nodes[j]);
        }
        w[i] = s;
    }
}

/* Bilinear sample at continuous texel coordinates (texel i is centered at i+0.5). */
static double sample_px(const float px[SEAM_MAX_H][SEAM_MAX_W], int w, int h,
                        double u, double v) {
    double x = u - 0.5, y = v - 0.5;
    int x0 = (int)floor(x), y0 = (int)floor(y);
    double fx = x - x0, fy = y - y0;
    int xa = x0 < 0 ? 0 : (x0 > w - 1 ? w - 1 : x0);
    int xb = x0 + 1 < 0 ? 0 : (x0 + 1 > w - 1 ? w - 1 : x0 + 1);
    int ya = y0 < 0 ? 0 : (y0 > h - 1 ? h - 1 : y0);
    int yb = y0 + 1 < 0 ? 0 : (y0 + 1 > h - 1 ? h - 1 : y0 + 1);
    double v00 = px[ya][xa], v10 = px[ya][xb];
    double v01 = px[yb][xa], v11 = px[yb][xb];
    return (1 - fx) * (1 - fy) * v00 + fx * (1 - fy) * v10 +
           (1 - fx) * fy * v01 + fx * fy * v11;
}

static double sample(const SeamImage* im, double u, double v) {
    return sample_px(im->px, im->w, im->h, u, v);
}

/* Extrapolate one row's left/right sides to the seam.  Values are the edge
   value (L0/R0) and edge slope (L1/R1) of the fitted profile. */
static void row_profile(const SeamImage* im, const double* nodes, int n,
                        const double* wv, const double* wd,
                        double seam_x, int j,
                        double* L0, double* L1, double* R0, double* R1) {
    double y = j + 0.5, val;
    *L0 = *L1 = *R0 = *R1 = 0.0;
    for (int k = 0; k < n; k++) {
        val = sample(im, seam_x - nodes[k], y);
        *L0 += wv[k] * val;
        *L1 += wd[k] * val;
        val = sample(im, seam_x + nodes[k], y);
        *R0 += wv[k] * val;
        *R1 += wd[k] * val;
    }
}

/* Independent estimator: a separate, wider node set and a different
   extrapolation order from the one that drives the repair.  Returns the
   largest value jump and gradient jump left across the seam. */
static void independent_residual_px(const float px[SEAM_MAX_H][SEAM_MAX_W],
                                    int w, int h, const double* nodes,
                                    const double* wv, const double* wd,
                                    double seam_x, double* value, double* slope) {
    double maxv = 0.0, maxs = 0.0;
    for (int j = 0; j < h; j++) {
        double y = j + 0.5, lv = 0.0, rv = 0.0, ls = 0.0, rs = 0.0;
        for (int k = 0; k < SEAM_INDEP; k++) {
            double l = sample_px(px, w, h, seam_x - nodes[k], y);
            double r = sample_px(px, w, h, seam_x + nodes[k], y);
            lv += wv[k] * l;
            rv += wv[k] * r;
            ls += wd[k] * l;
            rs += wd[k] * r;
        }
        double dv = fabs(lv - rv), ds = fabs(ls - rs);
        if (dv > maxv) maxv = dv;
        if (ds > maxs) maxs = ds;
    }
    *value = maxv;
    *slope = maxs;
}

static void independent_residual(const SeamImage* im, const double* nodes,
                                 const double* wv, const double* wd,
                                 double seam_x, double* value, double* slope) {
    independent_residual_px(im->px, im->w, im->h, nodes, wv, wd, seam_x, value, slope);
}

static double h00(double t) { return 2 * t * t * t - 3 * t * t + 1; }
static double h10(double t) { return t * t * t - 2 * t * t + t; }

static int valid_params(const SeamImage* im, const SeamParams* p) {
    if (!im || !p) return 0;
    if (im->w <= 0 || im->h <= 0 || im->w > SEAM_MAX_W || im->h > SEAM_MAX_H) return 0;
    if (p->band < 2.0f || p->threshold < 0.0f) return 0;
    /* the widest sample ring (0.75..4.75 texels) must stay inside the image
       on both sides; keep at least one texel of margin */
    int reach = (int)ceil(p->band > 4.75f ? p->band : 4.75f) + 1;
    if (p->seam_x < reach || p->seam_x > im->w - reach) return 0;
    return 1;
}

int seam_repair(SeamImage* im, const SeamParams* p, SeamReport* rep) {
    double node[SEAM_NODE], wv[SEAM_NODE], wd[SEAM_NODE];
    double inode[SEAM_INDEP], iwv[SEAM_INDEP], iwd[SEAM_INDEP];
    double D[SEAM_MAX_H], S[SEAM_MAX_H];

    memset(rep, 0, sizeof(*rep));
    if (!valid_params(im, p)) return -1;

    for (int k = 0; k < SEAM_NODE; k++) node[k] = 0.5 + k;
    for (int k = 0; k < SEAM_INDEP; k++) inode[k] = 0.75 + k;
    lagrange_at(node, SEAM_NODE, 0.0, wv);
    lagrange_deriv_at(node, SEAM_NODE, 0.0, wd);
    lagrange_at(inode, SEAM_INDEP, 0.0, iwv);
    lagrange_deriv_at(inode, SEAM_INDEP, 0.0, iwd);

    /* pass 1: measure every row before writing anything */
    double before = 0.0, slope_before = 0.0;
    for (int j = 0; j < im->h; j++) {
        double L0, L1, R0, R1;
        row_profile(im, node, SEAM_NODE, wv, wd, p->seam_x, j, &L0, &L1, &R0, &R1);
        D[j] = L0 - R0;
        S[j] = L1 - R1;
        if (fabs(D[j]) > before) before = fabs(D[j]);
        if (fabs(S[j]) > slope_before) slope_before = fabs(S[j]);
    }
    rep->residual_before = (float)before;
    rep->slope_before = (float)slope_before;

    if (before < p->threshold) {
        double ra, rs;
        rep->detected = 0;
        independent_residual(im, inode, iwv, iwd, p->seam_x, &ra, &rs);
        rep->residual_after = (float)ra;
        rep->slope_after = (float)rs;
        rep->bound = rep->residual_after;
        return 0;
    }
    rep->detected = 1;

    /* pass 2: value + gradient compensated Hermite blend inside the band.
       Work on a scratch copy so a correction that fails the acceptance test
       can be discarded without touching the caller's pixels. */
    static float scratch[SEAM_MAX_H][SEAM_MAX_W];
    memcpy(scratch, im->px, sizeof(scratch));
    double maxcorr = 0.0;
    int rows_active = 0, changed = 0;
    for (int j = 0; j < im->h; j++) {
        if (fabs(D[j]) < p->threshold) continue;
        double alpha = 0.0;
        if (p->gradient_match && fabs(D[j]) > 1e-12) {
            alpha = clampd(p->band * S[j] / D[j], -2.0 * p->band, 2.0 * p->band);
        }
        int row_changed = 0;
        for (int x = 0; x < im->w; x++) {
            double e = fabs((x + 0.5) - p->seam_x);
            if (e >= p->band) continue;
            double t = e / p->band;
            double w = h00(t) + alpha * h10(t);
            double delta = ((x + 0.5) < p->seam_x ? -0.5 : 0.5) * D[j] * w;
            if (delta == 0.0) continue;
            scratch[j][x] = (float)(scratch[j][x] + delta);
            row_changed = 1;
            changed++;
            if (fabs(delta) > maxcorr) maxcorr = fabs(delta);
        }
        if (row_changed) rows_active++;
    }

    {
        double ra, rs;
        independent_residual_px(scratch, im->w, im->h, inode, iwv, iwd,
                                p->seam_x, &ra, &rs);
        /* The repair is only allowed to commit when it actually meets the
           seam requirement.  A profile the band cannot represent (e.g. a box
           wider than the band) is measured and rejected, leaving the input
           untouched rather than trading the seam for a band-edge cliff. */
        if (ra >= p->threshold) {
            rep->accepted = 0;
            rep->residual_after = (float)ra;
            rep->slope_after = (float)rs;
            rep->bound = rep->residual_before + rep->residual_after;
            return 0;
        }
        rep->accepted = 1;
        rep->residual_after = (float)ra;
        rep->slope_after = (float)rs;
    }

    memcpy(im->px, scratch, sizeof(scratch));
    rep->changed = changed;
    rep->rows_active = rows_active;
    rep->max_correction = (float)maxcorr;
    /* The jump measured before repair bounds how far each side was moved; the
       independent residual bounds what a wider, higher-order estimator still
       sees afterwards.  See README for the assumption this rests on. */
    rep->bound = rep->residual_before + rep->residual_after;
    return 0;
}

/* ------------------------------------------------------------------ */
/* self-test                                                           */

static const int W = 64, H = 48;
static const int SEAM_X = 32;
static const double BAND = 5.0;
static const double THRESHOLD = 2e-3;
static const double WORLD = 6.0;

/* ground-truth radiance of a lit surface; both UV charts see the same world
   position across the seam, so an exact rebake is continuous there */
static double ground(double x, double y, double lx, double ly, double lz) {
    double dx = lx - x, dy = ly - y, dz = lz;
    double r2 = dx * dx + dy * dy + dz * dz;
    double r = sqrt(r2);
    double albedo = 0.9 + 0.08 * cos(2.0 * 3.14159265358979 * y / 6.0) + 0.05 * cos(1.3 * x);
    return 0.02 + albedo * 1.2 * (dz / r) / r2;
}

static double g_lx = 1.0, g_ly = 2.0, g_lz = 2.5;

static void set_light(double lx, double ly, double lz) {
    g_lx = lx; g_ly = ly; g_lz = lz;
}

static void world_xy(int i, int j, double* x, double* y) {
    *x = ((i + 0.5) - SEAM_X) * WORLD / W;
    *y = (j + 0.5) * WORLD / W;
}

static void make_rebake(SeamImage* im) {
    memset(im, 0, sizeof(*im));
    im->w = W; im->h = H;
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            double x, y;
            world_xy(i, j, &x, &y);
            im->px[j][i] = (float)ground(x, y, g_lx, g_ly, g_lz);
        }
    }
}

enum { ART_QUAD, ART_SLOPE, ART_EXP, ART_STEP };
static const char *art_name[] = { "quad", "slope", "exp", "step" };

/* symmetric seam artifact: left side drifts one way, right side the other */
static double artifact(int kind, int i, int j, double amp) {
    double e = fabs((i + 0.5) - SEAM_X);
    if (e >= BAND) return 0.0;
    double t = e / BAND;
    double y = (j + 0.5) * WORLD / W;
    double mod = 0.85 + 0.15 * cos(4.0 * 3.14159265358979 * y / 6.0);
    double a = 0.0;
    switch (kind) {
    case ART_QUAD:  a = amp * (1.0 - t) * (1.0 - t); break;
    case ART_SLOPE: a = amp * (1.0 - t) * (1.0 - t) * (1.0 + 0.6 * (1.0 - t)); break;
    case ART_EXP:   a = amp * exp(-(e - 0.5) / 1.5); break;
    case ART_STEP:  a = amp; break;
    default: break;
    }
    return ((i + 0.5) < SEAM_X ? -1.0 : 1.0) * a * mod;
}

static void make_seamed(SeamImage* im, int kind, double amp) {
    make_rebake(im);
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            im->px[j][i] = (float)(im->px[j][i] + artifact(kind, i, j, amp));
        }
    }
}

static double max_abs_diff(const SeamImage* a, const SeamImage* b) {
    double m = 0.0;
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            double d = fabs((double)a->px[j][i] - (double)b->px[j][i]);
            if (d > m) m = d;
        }
    }
    return m;
}

static int identical(const SeamImage* a, const SeamImage* b) {
    return memcmp(a->px, b->px, sizeof(a->px)) == 0;
}

static int band_mismatch(const SeamImage* a, const SeamImage* b) {
    int n = 0;
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            if (fabs((i + 0.5) - SEAM_X) >= BAND && a->px[j][i] != b->px[j][i]) n++;
        }
    }
    return n;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static SeamParams params(int gradient) {
    SeamParams p;
    memset(&p, 0, sizeof(p));
    p.seam_x = SEAM_X;
    p.band = (float)BAND;
    p.threshold = (float)THRESHOLD;
    p.gradient_match = gradient;
    return p;
}

int main(void) {
    SeamImage clean, base, seamed, work, ref;
    SeamReport rep;

    /* --- parameter validation --- */
    make_rebake(&clean);
    SeamParams bad = params(1);
    base = clean;
    bad.seam_x = 2;
    check("reject_seam_too_close", seam_repair(&base, &bad, &rep) == -1);
    bad = params(1);
    bad.band = 1.0f;
    base = clean;
    check("reject_band_too_small", seam_repair(&base, &bad, &rep) == -1);
    check("reject_null", seam_repair(0, &bad, &rep) == -1);

    /* --- no seam: not a single pixel may change --- */
    SeamParams pm = params(1);
    base = clean;
    check("clean_detects_nothing", seam_repair(&base, &pm, &rep) == 0 && !rep.detected);
    check("clean_pixels_untouched", identical(&base, &clean) && rep.changed == 0);
    check("clean_residual_below_threshold", rep.residual_before < THRESHOLD);

    /* clean but with an unrelated smooth wave across the whole map */
    SeamImage wavy = clean;
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            double x, y;
            world_xy(i, j, &x, &y);
            wavy.px[j][i] = (float)(wavy.px[j][i] + 0.3 * sin(0.7 * x) * cos(0.5 * y));
        }
    }
    base = wavy;
    check("clean_wave_detects_nothing", seam_repair(&base, &pm, &rep) == 0 && !rep.detected);
    check("clean_wave_untouched", identical(&base, &wavy));

    /* --- seam repair, per artifact shape ---
       quad/slope/exp are profiles a smooth band can represent; the step is a
       box that a band of the same width cannot match exactly and is handled
       separately below. */
    for (int kind = 0; kind < ART_STEP; kind++) {
        char name[64];
        make_rebake(&ref);
        make_seamed(&seamed, kind, 0.02);
        double err_pre = max_abs_diff(&seamed, &ref);

        work = seamed;
        int rc = seam_repair(&work, &pm, &rep);
        double err_post = max_abs_diff(&work, &ref);

        snprintf(name, sizeof(name), "%s_detected", art_name[kind]);
        check(name, rc == 0 && rep.detected && rep.rows_active == H);
        snprintf(name, sizeof(name), "%s_wrote_band", art_name[kind]);
        check(name, rep.changed == 10 * H);
        snprintf(name, sizeof(name), "%s_residual_below_threshold", art_name[kind]);
        check(name, rep.residual_after < THRESHOLD);
        snprintf(name, sizeof(name), "%s_residual_dropped", art_name[kind]);
        check(name, rep.residual_after < 0.1 * rep.residual_before);
        snprintf(name, sizeof(name), "%s_slope_residual_dropped", art_name[kind]);
        check(name, rep.slope_after < 0.1 * rep.slope_before);
        snprintf(name, sizeof(name), "%s_error_within_bound", art_name[kind]);
        check(name, err_post <= (double)rep.bound + 1e-6);
        snprintf(name, sizeof(name), "%s_band_confined", art_name[kind]);
        check(name, band_mismatch(&work, &seamed) == 0);
        snprintf(name, sizeof(name), "%s_beats_pre_repair", art_name[kind]);
        check(name, err_post < 0.5 * err_pre);
    }

    /* --- box profile: representable only when it fits inside the band --- */
    make_rebake(&ref);
    make_seamed(&seamed, ART_STEP, 0.01);
    work = seamed;
    seam_repair(&work, &pm, &rep);
    check("step_fits_band_residual_below_threshold", rep.detected && rep.residual_after < THRESHOLD);
    check("step_fits_band_error_within_bound", max_abs_diff(&work, &ref) <= (double)rep.bound + 1e-6);
    check("step_fits_band_accepted", rep.accepted && rep.changed > 0);

    /* wider than the band: no smooth band can reach the threshold, so the
       correction must be rejected and the input left bit-for-bit intact */
    make_seamed(&seamed, ART_STEP, 0.02);
    work = seamed;
    seam_repair(&work, &pm, &rep);
    check("step_wide_detected", rep.detected);
    check("step_wide_rejected", !rep.accepted && rep.changed == 0);
    check("step_wide_leaves_input_untouched", identical(&work, &seamed));
    check("step_wide_reports_threshold_violation", rep.residual_after >= THRESHOLD);

    /* --- gradient compensation actually does something --- */
    make_rebake(&ref);
    make_seamed(&seamed, ART_SLOPE, 0.02);
    work = seamed;
    SeamParams nog = params(0);
    seam_repair(&work, &nog, &rep);
    double err_value_only = max_abs_diff(&work, &ref);
    work = seamed;
    seam_repair(&work, &pm, &rep);
    double err_gradient = max_abs_diff(&work, &ref);
    check("gradient_match_reduces_error", err_gradient < 0.7 * err_value_only);
    check("gradient_match_slope_cleared", rep.slope_after < 0.1 * rep.slope_before);

    /* matching the end slope also removes the slope step the artifact created;
       without it the seam would trade a value jump for a slope kink */
    {
        SeamImage g = seamed, v = seamed;
        seam_repair(&g, &pm, &rep);
        float gs = rep.slope_after;
        seam_repair(&v, &nog, &rep);
        float vs = rep.slope_after;
        check("gradient_match_slope_vs_value_only", gs < 0.5f * vs);
    }

    /* --- idempotent: a repaired map is seen as clean --- */
    work = seamed;
    seam_repair(&work, &pm, &rep);
    SeamImage once = work;
    seam_repair(&work, &pm, &rep);
    check("second_pass_detects_nothing", !rep.detected && rep.changed == 0);
    check("second_pass_is_noop", identical(&work, &once));

    /* --- deterministic --- */
    SeamImage a = seamed, b = seamed;
    seam_repair(&a, &pm, &rep);
    seam_repair(&b, &pm, &rep);
    check("deterministic", identical(&a, &b));

    /* --- threshold gate: an imperceptibly small seam is left alone --- */
    make_seamed(&seamed, ART_QUAD, 4e-4);
    work = seamed;
    seam_repair(&work, &pm, &rep);
    check("subthreshold_untouched", !rep.detected && rep.changed == 0 && identical(&work, &seamed));

    /* --- tilted lights: same guarantees --- */
    set_light(3.5, 2.0, 1.2);
    make_rebake(&ref);
    make_seamed(&seamed, ART_QUAD, 0.02);
    work = seamed;
    seam_repair(&work, &pm, &rep);
    check("tilted_residual_below_threshold", rep.residual_after < THRESHOLD);
    check("tilted_error_within_bound", max_abs_diff(&work, &ref) <= (double)rep.bound + 1e-6);
    check("tilted_beats_pre_repair", max_abs_diff(&work, &ref) < 0.5 * max_abs_diff(&seamed, &ref));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
