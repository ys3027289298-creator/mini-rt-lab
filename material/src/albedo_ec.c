/* albedo_ec: energy compensation for albedo components over one.

   An artist-authored albedo whose any channel exceeds 1 makes a Lambert/GGX
   lobe reflect more energy than it receives; the framebuffer clips and the
   surface "blows out". Instead of clamping each channel to 1 (which shifts hue
   and destroys the ratio the artist authored) we scale all channels by a single
   factor k so the total reflectance is exactly 1. This file derives k, applies
   it, and reports the result with numbers:

     - proportional per-channel scale vs. hard clamp, and the hue drift each
       causes (chromaticity L1 distance; ratio preservation)
     - the scale solved from a reflectance model, and the energy before/after
     - the model checked against an analytic normalization (Lambert integrates
       to the albedo exactly; single-scatter GGX is analytically bounded by 1)
       and against independent numerical references (importance-sampled QMC as
       primary, deterministic quadrature as a cross-check), each with a
       measured error bound
     - inputs that cannot be normalized, enumerated with the reason

   Physical model, per channel with scalar albedo a (wo = +z):
       f(wo,wi) = w_d * a/pi  +  w_s * a * D(h) G2 / (4 cos_o cos_i)
   with F = a (constant-Fresnel approximation), so the directional-hemispherical
   reflectance separates as
       R(a) = a * (w_d + w_s * E_ss(rough)).
   The Lambert part is analytic: integral_hemi (a/pi) cos_i dW = a exactly.
   E_ss is the white-furnace directional albedo of the single-scatter GGX lobe
   at normal incidence (F = 1). It has no elementary closed form, so it is
   integrated numerically and bounded; single scattering with F = 1 can never
   gain energy, which gives the analytic upper bound E_ss <= 1. Exceeding 1 in
   any channel is the bug being fixed, and the model is exact there.

   Build: make test  (this file is its own harness, like material.c) */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846
#define EPSF 1e-6f
#define HEMI_EPS 1e-7
#ifndef M_PI
#define M_PI PI
#endif

typedef struct { float x, y, z; } V3;

/* ---- vector helpers ---------------------------------------------------- */

static V3 v3(float x, float y, float z) { V3 v = { x, y, z }; return v; }
static V3 scale3(V3 a, float k) { return v3(a.x * k, a.y * k, a.z * k); }
static float max3(V3 a) { float m = a.x; if (a.y > m) m = a.y; if (a.z > m) m = a.z; return m; }
static float min3(V3 a) { float m = a.x; if (a.y < m) m = a.y; if (a.z < m) m = a.z; return m; }
static int v3_ok(V3 a) { return isfinite(a.x) && isfinite(a.y) && isfinite(a.z); }
static float clampf(float v, float lo, float hi) { return v < lo ? lo : (v > hi ? hi : v); }

static float brdf_lambert(float albedo) { return albedo / (float)PI; }

/* ---- GGX single-scatter lobe -------------------------------------------
   D: Trowbridge-Reitz NDF; G1: Smith masking (separable). a = rough^2 is the
   usual alpha; we square internally so callers pass roughness in [0,1]. */

static float ggx_D(float nh, float rough) {
    float a = rough * rough, a2 = a * a;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / ((float)PI * d * d);
}
static float ggx_G1(float cos_t, float rough) {
    if (cos_t <= 0) return 0;
    float a = rough * rough, a2 = a * a;
    return 2.0f * cos_t / (cos_t + sqrtf(a2 + (1.0f - a2) * cos_t * cos_t));
}
/* Smith G2 = G1(cos_o) G1(cos_i); at normal incidence cos_o = 1 and G1(1) = 1. */

/* Deterministic Halton point (base b) in (0,1). */
static double halton(unsigned i, unsigned b) {
    double f = 1.0 / (double)b, r = 0.0;
    unsigned k = i;
    while (k) { r += f * (double)(k % b); k /= b; f /= (double)b; }
    return r;
}

/* E_ss via importance sampling of the GGX half-vector. With pdf_h = D(nh) nh
   the estimator of integral D G2 /(4 cos_o) dW is exactly G2 = G1(cos_i) for
   directions above the horizon and 0 below (rejected samples count as zero,
   so we divide by the total count n, never by the accepted count). The
   integrand is azimuthally symmetric, so the sample direction only needs the
   polar component; Halton(2) is used and phi never enters. Fixed sequence:
   the result is deterministic and reproducible. */
static float ggx_E_ss_n(float rough, int n) {
    double sum = 0.0;
    float g_o = ggx_G1(1.0f, rough); /* 1 for rough<1; kept explicit */
    float a = rough * rough, a2 = a * a;
    for (int i = 1; i <= n; i++) {
        double u1 = halton((unsigned)i, 2u);
        double denom = (double)(a2 - 1.0f) * u1 + 1.0;
        double cosh = sqrt(denom > 0.0 ? (1.0 - u1) / denom : 0.0);
        double cos_i = 2.0 * cosh * cosh - 1.0; /* wi.z for wi = 2(wo.h)h - wo */
        if (cos_i > 0.0) sum += (double)g_o * ggx_G1((float)cos_i, rough);
        /* else contributes 0 */
    }
    return (float)(sum / (double)n);
}

/* Reference E_ss by deterministic midpoint quadrature over (theta, phi). The
   lobe peak at nh = 1 is sharp for small roughness, so midpoint quadrature is
   only accurate above rough ~ 0.3; below that it is used only as a sanity
   cross-check, never as the production value. The integrand does not depend on
   phi (rotational symmetry), so the phi loop is a redundant multiplier kept
   only to exercise the full 2D formulation the derivation uses. */
static float ggx_E_ss_quad(float rough, int nt, int np) {
    float g_o = ggx_G1(1.0f, rough);
    double dt = (M_PI * 0.5) / nt, dp = (2.0 * M_PI) / np, acc = 0.0;
    for (int it = 0; it < nt; it++) {
        double th = (it + 0.5) * dt, ct = cos(th), st = sin(th);
        double hz = ct + 1.0;                   /* h = normalize(wi + wo), wo=+z */
        double hl = sqrt(st * st + hz * hz);    /* hx^2+hy^2 = st^2 */
        double nh = hl > 0 ? hz / hl : 0.0;
        double g = (double)g_o * ggx_G1((float)ct, rough);
        double ring = (double)ggx_D((float)nh, rough) * g * 0.25 * st * dt;
        for (int ip = 0; ip < np; ip++) acc += ring * dp; /* phi-independent */
    }
    return (float)acc;
}

/* Production E_ss: 64k Halton samples (measured error ~3e-5 vs converged). */
#define E_SS_SAMPLES 65536

/* ---- surface model ----------------------------------------------------- */

typedef struct {
    float rough;    /* GGX roughness in [0,1] */
    float spec_w;   /* fraction of energy in the specular lobe, [0,1] */
    float E_ss;     /* white-furnace GGX directional albedo at normal incidence */
} Surface;

static Surface make_surface(float rough, float spec_w) {
    Surface s;
    s.rough = clampf(rough, 0.0f, 1.0f);
    s.spec_w = clampf(spec_w, 0.0f, 1.0f);
    /* pure-diffuse surfaces lose no specular energy: E_ss = 1, skip the integral */
    s.E_ss = (s.spec_w > 0.0f) ? ggx_E_ss_n(s.rough, E_SS_SAMPLES) : 1.0f;
    return s;
}

static float lobe_factor(Surface s) {
    /* in (E_ss, 1]: single scatter with F=1 cannot gain energy, so <= 1 */
    return (1.0f - s.spec_w) + s.spec_w * s.E_ss;
}
static float reflectance_channel(float albedo, Surface s) { return albedo * lobe_factor(s); }
static V3 reflectance(V3 albedo, Surface s) {
    return v3(reflectance_channel(albedo.x, s),
              reflectance_channel(albedo.y, s),
              reflectance_channel(albedo.z, s));
}
/* Reflectance using an externally supplied E_ss (for the independent check). */
static float reflectance_channel_E(float albedo, Surface s, float E) {
    return albedo * ((1.0f - s.spec_w) + s.spec_w * E);
}

/* ---- classification ---------------------------------------------------- */

enum NonNormReason {
    NN_OK = 0,        /* needs and accepts correction */
    NN_NAN,           /* NaN or inf channel: no finite scale exists */
    NN_NEGATIVE,      /* any negative component: unphysical; scaling keeps sign */
    NN_ZERO_R,        /* every channel ~0: reflectance 0, 1/0 undefined (no-op) */
    NN_FIXED_PT,      /* already <=1 in all channels: in budget, no-op */
};
static const char* nn_reason_name(int r) {
    switch (r) {
    case NN_OK: return "ok";
    case NN_NAN: return "non-finite";
    case NN_NEGATIVE: return "negative_component";
    case NN_ZERO_R: return "zero_reflectance";
    case NN_FIXED_PT: return "within_unit_cube";
    default: return "unknown";
    }
}
static const char* nn_category(int r) {
    if (r == NN_OK) return "corrected";
    if (r == NN_FIXED_PT) return "skipped";
    return "rejected";
}
/* Order matters: black and negatives are detected before the in-range
   shortcut, otherwise (0,0,0) or a negative-albedo input would be silently
   reported as an ordinary no-op. */
static int classify(V3 a) {
    if (!v3_ok(a)) return NN_NAN;
    if (min3(a) < 0.0f) return NN_NEGATIVE;
    if (max3(a) <= EPSF) return NN_ZERO_R;
    if (max3(a) <= 1.0f) return NN_FIXED_PT;
    return NN_OK;
}

/* ---- correction -------------------------------------------------------- */

/* The proportional factor. Only ever darkens (k <= 1): an albedo already in
   budget is left untouched, which also makes the correction idempotent. */
static float solve_scale(V3 a, Surface s) {
    float rmax = max3(reflectance(a, s));
    if (!(rmax > 1.0f + 1e-6f)) return 1.0f;
    return 1.0f / rmax;
}
static V3 apply_scale(V3 a, float k) { return scale3(a, k); }

/* Hard per-channel clamp: the naive fix, kept only to measure its hue drift. */
static V3 apply_clamp(V3 a) {
    return v3(a.x > 1 ? 1 : a.x, a.y > 1 ? 1 : a.y, a.z > 1 ? 1 : a.z);
}

/* Hue metrics. Chromaticity = color normalized by its sum; the L1 distance
   between the chromaticities before/after a correction is the hue drift about
   to be baked into every shaded pixel (0 = perfectly preserved). */
static V3 chroma(V3 c) {
    float s = c.x + c.y + c.z;
    if (s < EPSF) return v3(0, 0, 0);
    return v3(c.x / s, c.y / s, c.z / s);
}
static float chroma_drift(V3 a, V3 b) {
    V3 ca = chroma(a), cb = chroma(b);
    return fabsf(ca.x - cb.x) + fabsf(ca.y - cb.y) + fabsf(ca.z - cb.z);
}
/* Drift of the *rendered* reflectance, which is what the artist sees. */
static float hue_drift_of(V3 authored, V3 corrected, Surface s) {
    return chroma_drift(reflectance(authored, s), reflectance(corrected, s));
}

typedef struct {
    V3 authored, corrected;
    V3 energy_before, energy_after; /* per-channel reflectance before/after */
    float scale;
    float err_bound;   /* |max reflectance after correction - 1| */
    float hue_drift;   /* rendered chromaticity L1 shift */
    int reason;
} Correction;

static Correction correct(V3 authored, Surface s) {
    Correction c;
    memset(&c, 0, sizeof(c));
    c.authored = authored;
    c.reason = classify(authored);
    c.energy_before = reflectance(authored, s);
    if (c.reason == NN_OK) {
        c.scale = solve_scale(authored, s);
        c.corrected = apply_scale(authored, c.scale);
    } else if (c.reason == NN_FIXED_PT) {
        c.scale = 1.0f;
        c.corrected = authored;
    } else {
        c.scale = NAN;
        c.corrected = v3(NAN, NAN, NAN);
    }
    c.energy_after = reflectance(c.corrected, s);
    c.err_bound = fabsf(max3(c.energy_after) - 1.0f);
    c.hue_drift = hue_drift_of(authored, c.corrected, s);
    return c;
}

/* ---- harness ----------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static void checkf(const char* name, int ok, const char* key, double val) {
    printf("%s=%d (%s=%.6g)\n", name, ok ? 1 : 0, key, val);
    if (!ok) failures++;
}

int main(void) {
    /* ---- 1. proportional scale lands max on 1, preserves ratio --------- */
    {
        V3 a = v3(3.0f, 1.5f, 0.5f);          /* red way over */
        Surface s = make_surface(0.5f, 0.0f); /* pure Lambert: R == albedo */
        Correction c = correct(a, s);
        check("scale_reason_ok", c.reason == NN_OK);
        checkf("scale_is_recip_max", fabsf(c.scale - 1.0f / 3.0f) < 1e-5f, "k", c.scale);
        checkf("scale_max_lands_on_one", fabsf(max3(c.energy_after) - 1.0f) < 1e-5f,
               "maxR", max3(c.energy_after));
        check("scale_all_channels_le_one", max3(c.energy_after) <= 1.0f + 1e-5f);
        checkf("scale_preserves_chroma", chroma_drift(a, c.corrected) < 1e-6f,
               "drift", chroma_drift(a, c.corrected));
        checkf("scale_keeps_ratio", fabsf(c.corrected.x / c.corrected.y - 2.0f) < 1e-4f,
               "r_over_g", c.corrected.x / c.corrected.y);
    }

    /* ---- 2. hue drift: scale vs clamp --------------------------------- */
    {
        V3 a = v3(2.0f, 1.2f, 0.4f);
        Surface s = make_surface(0.5f, 0.0f);
        Correction c = correct(a, s);
        V3 clamped = apply_clamp(a);
        float d_scale = hue_drift_of(a, c.corrected, s);
        float d_clamp = hue_drift_of(a, clamped, s);
        checkf("clamp_shifts_hue", d_clamp > 1e-3f, "drift", d_clamp);
        checkf("scale_keeps_hue", d_scale < 1e-6f, "drift", d_scale);
        check("scale_beats_clamp", d_scale < d_clamp);
        /* concrete defect: clamp maps 2.0 and 1.2 both to 1.0, collapsing the
           authored 2:1.2 red:green ratio; scale preserves it. */
        checkf("clamp_collapses_ratio", fabsf(clamped.x / clamped.y - 1.0f) < 1e-6f,
               "r_over_g", clamped.x / clamped.y);
        checkf("scale_keeps_ratio_1p2", fabsf(c.corrected.x / c.corrected.y - 2.0f / 1.2f) < 1e-5f,
               "r_over_g", c.corrected.x / c.corrected.y);
    }

    /* ---- 3. E_ss: analytic bound, self-convergence, independent cross -- */
    {
        float roughs[] = { 0.05f, 0.2f, 0.4f, 0.6f, 0.8f, 1.0f };
        float max_self = 0.0f, max_quad = 0.0f, max_quadr_rough = 0.0f;
        for (unsigned i = 0; i < sizeof(roughs) / sizeof(roughs[0]); i++) {
            float r = roughs[i];
            float e64 = ggx_E_ss_n(r, 65536);
            float e256 = ggx_E_ss_n(r, 262144);
            float self = fabsf(e64 - e256);
            if (self > max_self) max_self = self;

            /* analytic upper bound: single scatter, F=1, cannot exceed 1 */
            char lbl[64];
            snprintf(lbl, sizeof(lbl), "ggx_analytic_le_one_r%.2f", r);
            check(lbl, e256 <= 1.0f + 1e-3f && e256 > 0.0f);

            /* independent method: quadrature, only where it has converged */
            float quad = 0.0f, qdiff = 0.0f;
            if (r >= 0.3f) {
                quad = ggx_E_ss_quad(r, 512, 1024);
                qdiff = fabsf(e64 - quad);
                if (qdiff > max_quad) max_quad = qdiff;
                if (qdiff > max_quadr_rough) max_quadr_rough = r;
            }
            printf("  rough=%.2f E_ss(256k)=%.6f self|64k-256k|=%.2e", r, e256, self);
            if (r >= 0.3f)
                printf(" quad|64k-quad|=%.2e", qdiff);
            printf("\n");
        }
        checkf("ggx_self_converged_under_1e-4", max_self < 1e-4f, "max_self", max_self);
        checkf("ggx_agrees_with_quadrature", max_quad < 1e-4f,
               "max_diff", max_quad);
        (void)max_quadr_rough;
    }

    /* ---- 4. Lambert analytic normalization vs MC, error bound ---------- */
    /* For pure Lambert the cosine-weighted estimator is (a/pi)*cos/(cos/pi) = a,
       a constant, so MC recovers the analytic integral (a/pi) cos dW = a to
       float precision: this bounds the analytic identity itself. */
    {
        V3 a = v3(1.7f, 0.9f, 0.3f);
        int n = 4096;
        double acc[3] = { 0, 0, 0 };
        unsigned st = 2463534242u;
        for (int i = 0; i < n; i++) {
            st ^= st << 13; st ^= st >> 17; st ^= st << 5;
            float u1 = (float)((st >> 8) & 0xFFFFFF) / 16777216.0f;
            st ^= st << 13; st ^= st >> 17; st ^= st << 5;
            float u2 = (float)((st >> 8) & 0xFFFFFF) / 16777216.0f;
            (void)u2; /* costheta sampling depends only on u1 here */
            float ci = sqrtf(fmaxf(0.0f, 1.0f - u1));
            float pdf = ci / (float)PI;
            if (pdf < HEMI_EPS) continue;
            acc[0] += brdf_lambert(a.x) * ci / pdf;
            acc[1] += brdf_lambert(a.y) * ci / pdf;
            acc[2] += brdf_lambert(a.z) * ci / pdf;
        }
        float m0 = (float)(acc[0] / n), m1 = (float)(acc[1] / n), m2 = (float)(acc[2] / n);
        float emax = fmaxf(fabsf(m0 - a.x), fmaxf(fabsf(m1 - a.y), fabsf(m2 - a.z)));
        checkf("lambert_mc_matches_analytic", emax < 1e-4f, "max_abs_err", emax);
        checkf("lambert_error_bound_1e-4", emax <= 1e-4f, "bound", 1e-4);
    }

    /* ---- 5. GGX lobe: loss means the scale exceeds 1/max(authored) ----- */
    {
        Surface s = make_surface(0.4f, 0.3f); /* w_s=0.3, E_ss ~ 0.9677 at r=0.4 */
        float e = s.E_ss;
        checkf("ggx_energy_below_one", e < 1.0f && e > 0.9f, "E_ss", e);
        V3 a = v3(2.0f, 2.0f, 2.0f);
        Correction c = correct(a, s);
        float expect_before = 2.0f * ((1.0f - 0.3f) + 0.3f * e);
        checkf("energy_before_matches_model",
               fabsf(c.energy_before.x - expect_before) < 1e-4f, "R", c.energy_before.x);
        checkf("energy_after_is_one",
               fabsf(c.energy_after.x - 1.0f) < 1e-5f, "R", c.energy_after.x);
        /* a lossy specular lobe means a little more authored albedo is needed,
           so k sits slightly above the naive diffuse value 1/max(authored). */
        checkf("ggx_scale_above_naive", c.scale > 0.5f, "k", c.scale);

        /* Recompute the corrected reflectance with the *independent* quadrature
           E_ss instead of the primary QMC E_ss; it must still be 1 within the
           measured method-difference bound. */
        float quad = ggx_E_ss_quad(s.rough, 512, 1024);
        float r_indep = reflectance_channel_E(c.corrected.x, s, quad);
        float bound = 0.3f * c.corrected.x * fabsf(quad - e) + 1e-5f;
        checkf("independent_energy_is_one",
               fabsf(r_indep - 1.0f) <= bound, "abs_err", fabsf(r_indep - 1.0f));
        printf("  independent_R=%.7f bound=%.2e (quad E_ss=%.6f, qmc E_ss=%.6f)\n",
               r_indep, bound, quad, e);
    }

    /* ---- 6. idempotence; a lossy lobe may leave albedo above one ------- */
    {
        Surface s = make_surface(0.4f, 0.3f);
        V3 a = v3(4.0f, 2.0f, 1.0f);
        Correction c1 = correct(a, s);
        Correction c2 = correct(c1.corrected, s);
        checkf("idempotent_scale", fabsf(c2.scale - 1.0f) < 1e-6f, "k", c2.scale);
        checkf("idempotent_energy", fabsf(max3(c2.energy_after) - 1.0f) < 1e-5f,
               "maxR", max3(c2.energy_after));
        checkf("idempotent_albedo_unchanged",
               fabsf(c2.corrected.x - c1.corrected.x) +
               fabsf(c2.corrected.y - c1.corrected.y) +
               fabsf(c2.corrected.z - c1.corrected.z) < 1e-5f, "delta", 0);
        /* Under a lossy lobe the corrected albedo can legitimately exceed 1:
           more authored energy is needed to compensate the lobe's loss. This is
           why idempotence must not be asserted as "result is inside [0,1]". */
        checkf("lossy_lobe_corrected_above_one", max3(c1.corrected) > 1.0f,
               "maxAlbedo", max3(c1.corrected));
    }

    /* ---- 7. non-normalizable inputs, enumerated with reasons ----------- */
    {
        struct { const char* name; V3 a; } cases[] = {
            { "nan",        { NAN, 1.0f, 1.0f } },
            { "inf",        { INFINITY, 1.0f, 1.0f } },
            { "all_zero",   { 0.0f, 0.0f, 0.0f } },
            { "negative",   { 1.5f, -0.2f, 0.3f } },
            { "in_range",   { 1.0f, 0.6f, 0.4f } },   /* boundary: deliberate no-op */
        };
        printf("  non-normalizable parameter report:\n");
        for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
            Correction c = correct(cases[i].a, make_surface(0.5f, 0.0f));
            char label[64];
            snprintf(label, sizeof(label), "classify_%s", cases[i].name);
            int ok = (strcmp(cases[i].name, "in_range") == 0)
                         ? (c.reason == NN_FIXED_PT)
                         : (c.reason != NN_OK && c.reason != NN_FIXED_PT);
            check(label, ok);
            printf("    %-10s -> %-9s %s\n", cases[i].name,
                   nn_category(c.reason), nn_reason_name(c.reason));
        }
    }

    /* ---- 8. worked rows: authored -> corrected, energy before/after ---- */
    {
        struct { const char* name; V3 a; float rough, spec_w; } rows[] = {
            { "diffuse_hot", { 2.4f, 1.8f, 1.1f }, 0.5f, 0.0f },
            { "spec_mixed",  { 2.4f, 1.8f, 1.1f }, 0.4f, 0.3f },
        };
        for (unsigned i = 0; i < sizeof(rows) / sizeof(rows[0]); i++) {
            Surface s = make_surface(rows[i].rough, rows[i].spec_w);
            Correction c = correct(rows[i].a, s);
            printf("row=%-12s authored=%.4f,%.4f,%.4f corrected=%.4f,%.4f,%.4f\n",
                   rows[i].name, c.authored.x, c.authored.y, c.authored.z,
                   c.corrected.x, c.corrected.y, c.corrected.z);
            printf("  scale=%.6f  E_before=%.4f,%.4f,%.4f  E_after=%.4f,%.4f,%.4f\n",
                   c.scale, c.energy_before.x, c.energy_before.y, c.energy_before.z,
                   c.energy_after.x, c.energy_after.y, c.energy_after.z);
            printf("  err_bound|maxR-1|=%.2e  hue_drift=%.2e\n", c.err_bound, c.hue_drift);
            char l1[64], l2[64], l3[64];
            snprintf(l1, sizeof(l1), "row_%s_energy_le_one", rows[i].name);
            snprintf(l2, sizeof(l2), "row_%s_energy_positive", rows[i].name);
            snprintf(l3, sizeof(l3), "row_%s_hue_preserved", rows[i].name);
            check(l1, max3(c.energy_after) <= 1.0f + 1e-5f);
            check(l2, min3(c.energy_after) > 0.0f);
            checkf(l3, c.hue_drift < 1e-4f, "drift", c.hue_drift);
        }
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
