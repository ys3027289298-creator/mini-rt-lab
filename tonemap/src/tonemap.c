
/* tonemap: exposure-consistent tone mapping.
   Pipeline per pixel:
     1. linear exposure scaling: c' = c * exposure (exposure is a linear gain,
        e.g. 2^ev for ev stops; must be > 0, zero/negative/NaN/Inf is rejected)
     2. luminance-only compression on Rec.709 luma L, via a single per-pixel
        gain g = curve(L) / L applied to all three channels, so the output is
        collinear with the exposed color and hue is preserved exactly
     3. analytic reference curve: Reinhard L/(1+L), evaluated in double, used
        to bound the error of the float pipeline; hue shift is reported as the
   angle (radians) between the pre/post-compression normalized RGB directions. */
#include <stdio.h>
#include <math.h>

typedef struct { float r, g, b; } Rgb;

/* Rec.709 luma weights */
static const float LUMA_R = 0.2126f;
static const float LUMA_G = 0.7152f;
static const float LUMA_B = 0.0722f;

/* explicit upper bounds on the float-pipeline error vs the double analytic
   Reinhard curve, verified over the self-test sweep */
#define ERR_ABS_BOUND 1e-5f
#define ERR_REL_BOUND 1e-5f
/* hue shift is structural (one shared gain), only float round-off remains */
#define HUE_SHIFT_BOUND 1e-5f

static float luma(Rgb c) {
    return LUMA_R * c.r + LUMA_G * c.g + LUMA_B * c.b;
}

/* analytic compression curve, double precision reference: Reinhard */
static double reinhard_ref(double L) { return L / (1.0 + L); }

/* same curve in the pipeline's float precision, stable reciprocal form */
static float reinhard_f(float L) { return L / (1.0f + L); }

/* returns 0 on success, -1 when exposure is invalid (out is zeroed). */
static int tonemap(Rgb in, float exposure, Rgb* out) {
    out->r = out->g = out->b = 0.0f;
    if (!(exposure > 0.0f) || isinf(exposure)) return -1;

    Rgb e = { in.r * exposure, in.g * exposure, in.b * exposure };
    float L = luma(e);
    float gain = (L > 0.0f) ? reinhard_f(L) / L : 0.0f;
    out->r = e.r * gain;
    out->g = e.g * gain;
    out->b = e.b * gain;
    return 0;
}

/* hue shift metric: angle between normalized RGB directions before and after
   compression, in radians; 0 for black or collinear colors */
static double hue_angle(Rgb a, Rgb b) {
    double na = sqrt((double)a.r * a.r + (double)a.g * a.g + (double)a.b * a.b);
    double nb = sqrt((double)b.r * b.r + (double)b.g * b.g + (double)b.b * b.b);
    if (na <= 0.0 || nb <= 0.0) return 0.0;
    double cosang = ((double)a.r * b.r + (double)a.g * b.g + (double)a.b * b.b) / (na * nb);
    if (cosang > 1.0) cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    return acos(cosang);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- exposure zero / invalid must be rejected ---- */
    Rgb in = { 0.5f, 1.0f, 2.0f }, out;
    check("reject_zero", tonemap(in, 0.0f, &out) == -1 &&
          out.r == 0.0f && out.g == 0.0f && out.b == 0.0f);
    check("reject_negative", tonemap(in, -2.0f, &out) == -1);
    check("reject_nan", tonemap(in, NAN, &out) == -1);
    check("reject_inf", tonemap(in, INFINITY, &out) == -1);
    check("accept_positive", tonemap(in, 1.0f, &out) == 0);

    /* ---- black stays black and never moves the hue metric ---- */
    Rgb black = { 0, 0, 0 };
    check("black_stays_black", tonemap(black, 4.0f, &out) == 0 &&
          out.r == 0.0f && out.g == 0.0f && out.b == 0.0f);

    /* ---- analytic curve spot checks ---- */
    check("curve_half", fabs(reinhard_f(1.0f) - 0.5f) < 1e-7f);
    check("curve_zero", reinhard_f(0.0f) == 0.0f);
    check("curve_saturates", reinhard_f(1e6f) < 1.0f);

    /* ---- main sweep: curve error bounds + hue shift metric ---- */
    Rgb colors[] = {
        { 0.0f, 0.0f, 0.0f },
        { 0.01f, 0.01f, 0.01f },
        { 0.2f, 0.5f, 0.1f },
        { 1.0f, 0.2f, 0.05f },
        { 0.1f, 0.8f, 2.5f },
        { 4.0f, 4.0f, 4.0f },
        { 8.0f, 1.5f, 0.2f },
        { 16.0f, 32.0f, 64.0f }
    };
    float exposures[] = { 0.125f, 0.5f, 1.0f, 2.0f, 8.0f, 64.0f };
    int ncol = (int)(sizeof(colors) / sizeof(colors[0]));
    int nexp = (int)(sizeof(exposures) / sizeof(exposures[0]));

    double max_abs = 0.0, max_rel = 0.0, max_hue = 0.0;
    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nexp; j++) {
            Rgb c = colors[i];
            Rgb exposed = { c.r * exposures[j], c.g * exposures[j], c.b * exposures[j] };
            Rgb mapped;
            int rc = tonemap(c, exposures[j], &mapped);
            if (rc != 0) { check("sweep_valid", 0); continue; }

            /* analytic reference: same exposure, then Reinhard on luma,
               rebuilt through the same shared gain so the per-channel
               reference stays exactly hue preserving */
            double Ld = LUMA_R * exposed.r + LUMA_G * exposed.g + LUMA_B * exposed.b;
            double gd = (Ld > 0.0) ? reinhard_ref(Ld) / Ld : 0.0;
            double ref[3] = { exposed.r * gd, exposed.g * gd, exposed.b * gd };

            double eabs = fabs((double)mapped.r - ref[0]);
            double e2 = fabs((double)mapped.g - ref[1]);
            double e3 = fabs((double)mapped.b - ref[2]);
            if (e2 > eabs) eabs = e2;
            if (e3 > eabs) eabs = e3;
            if (eabs > max_abs) max_abs = eabs;

            double mv[3] = { mapped.r, mapped.g, mapped.b };
            for (int k = 0; k < 3; k++) {
                double rv = ref[k];
                if (rv > 1e-3) {
                    double er = fabs((mv[k] - rv) / rv);
                    if (er > max_rel) max_rel = er;
                }
            }

            double dh = hue_angle(exposed, mapped);
            if (dh > max_hue) max_hue = dh;
        }
    }
    printf("curve_err_abs_max=%.3e bound=%.0e\n", max_abs, (double)ERR_ABS_BOUND);
    printf("curve_err_rel_max=%.3e bound=%.0e\n", max_rel, (double)ERR_REL_BOUND);
    printf("hue_shift_rad_max=%.3e bound=%.0e\n", max_hue, (double)HUE_SHIFT_BOUND);
    check("curve_err_abs_bounded", max_abs < (double)ERR_ABS_BOUND);
    check("curve_err_rel_bounded", max_rel < (double)ERR_REL_BOUND);
    check("hue_preserved", max_hue < (double)HUE_SHIFT_BOUND);

    /* ---- luminance-only: dark values pass through, all values compress ---- */
    Rgb dark = { 0.02f, 0.02f, 0.02f }, dm;
    tonemap(dark, 1.0f, &dm);
    check("dark_passthrough", fabsf(dm.r - reinhard_f(0.02f)) < 1e-6f);

    Rgb bright = { 10.0f, 10.0f, 10.0f }, bm;
    tonemap(bright, 1.0f, &bm);
    check("bright_compressed", bm.r < 1.0f && bm.r > 0.0f);

    /* ---- linear exposure scaling: pipeline output scales with exposure ---- */
    Rgb o1, o2;
    tonemap(in, 1.0f, &o1);
    tonemap(in, 2.0f, &o2);
    check("exposure_changes_output",
          fabsf(o2.r - o1.r) > 1e-4f || fabsf(o2.g - o1.g) > 1e-4f);

    /* exposure is a pure linear pre-scale: doubling input equals doubling gain */
    Rgb doubled = { in.r * 2.0f, in.g * 2.0f, in.b * 2.0f }, od;
    tonemap(doubled, 1.0f, &od);
    check("exposure_linear_prescale",
          fabsf(od.r - o2.r) < 1e-5f && fabsf(od.g - o2.g) < 1e-5f &&
          fabsf(od.b - o2.b) < 1e-5f);

    /* ---- monotone in exposure ---- */
    Rgb prev = { 0, 0, 0 }, cur;
    int mono = 1;
    for (int j = 0; j < nexp; j++) {
        tonemap(in, exposures[j], &cur);
        if (j > 0 && (cur.r < prev.r - 1e-6f || cur.g < prev.g - 1e-6f ||
                      cur.b < prev.b - 1e-6f))
            mono = 0;
        prev = cur;
    }
    check("exposure_monotone", mono);

    /* ---- channel ordering preserved (no cross-channel remix) ---- */
    Rgb ordered = { 0.5f, 1.0f, 2.0f }, om;
    tonemap(ordered, 4.0f, &om);
    check("channel_order_preserved", om.r <= om.g && om.g <= om.b);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
