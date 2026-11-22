
/* material: Lambert + GGX BRDF evaluation and cosine-weighted sampling,
   with a PDF check against the analytic cosine distribution, plus a
   subsurface diffusion approximation validated against Monte Carlo. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

typedef struct { float x, y, z; } V3;

static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}

/* Lambert: albedo / pi */
static float brdf_lambert(float albedo) { return albedo / 3.14159265f; }

/* GGX normal distribution */
static float ggx_ndf(float nh, float rough) {
    float a = rough * rough;
    float a2 = a * a;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / (3.14159265f * d * d);
}

/* cosine-weighted hemisphere: pdf = cos/pi */
static float cosine_pdf(float cos_theta) { return cos_theta / 3.14159265f; }

static V3 cosine_sample(float u1, float u2) {
    float r = sqrtf(u1);
    float phi = 2.0f * 3.14159265f * u2;
    V3 v = { r * cosf(phi), r * sinf(phi), sqrtf(1.0f - u1 > 0 ? 1.0f - u1 : 0) };
    return norm3(v);
}

/* ---- subsurface diffusion approximation ----
   Homogeneous, isotropically scattering half-space with a matched boundary
   (no refractive-index mismatch, so no internal reflection).  With
   sigma_t = 1/mfp, D = 1/(3 sigma_t) and sigma_tr = sqrt(3 (1-a) sigma_t^2),
   the energy from a point source at depth z_r = 1/sigma_t diffusing back out
   of the surface (dipole with image source at z_v = z_r + 4D) follows

       p(r) = [ z_r (1+sigma_tr d_r) e^{-sigma_tr d_r} / d_r^3
              + z_v (1+sigma_tr d_v) e^{-sigma_tr d_v} / d_v^3 ] / (4 pi),

   with d^2 = r^2 + z^2.  Its radial integral over the disk r <= R is the
   closed form in sss_capture(); over the whole plane it is
   Z = (e^{-sigma_tr z_r} + e^{-sigma_tr z_v}) / 2, so the model lets
   a*Z of the incident energy leave the medium and 1 - a*Z is lost to
   absorption.  Because the profile depends on R only through R/mfp, a
   capture radius chosen from mfp scales linearly with mfp. */

#define PI 3.14159265358979323846

typedef struct {
    float albedo;       /* single-scattering albedo in [0,1] */
    float mfp;          /* mean free path, world units */
    double sigma_t;     /* 1 / mfp */
    double sigma_tr;    /* effective attenuation */
    double D;           /* diffusion coefficient */
    double z_r, z_v;    /* source and image depth */
    double Z;           /* integral of the profile over the whole plane */
    double escaped;     /* albedo * Z: energy that can leave the medium */
} Sss;

/* Fill a medium description.  Returns NULL on success, otherwise the name of
   the offending parameter (albedo outside [0,1], mfp not finite and > 0). */
static const char* sss_init(Sss* s, float albedo, float mfp) {
    if (!isfinite(albedo) || albedo < 0.0f || albedo > 1.0f) return "albedo";
    if (!isfinite(mfp) || mfp <= 0.0f) return "mfp";
    s->albedo = albedo;
    s->mfp = mfp;
    s->sigma_t = 1.0 / mfp;
    s->D = 1.0 / (3.0 * s->sigma_t);
    s->sigma_tr = sqrt(3.0 * (1.0 - (double)albedo) * s->sigma_t * s->sigma_t);
    s->z_r = 1.0 / s->sigma_t;
    s->z_v = s->z_r + 4.0 * s->D;
    s->Z = 0.5 * (exp(-s->sigma_tr * s->z_r) + exp(-s->sigma_tr * s->z_v));
    s->escaped = albedo * s->Z;
    return NULL;
}

/* Outgoing energy per unit area at surface distance r from the source. */
static double sss_profile(const Sss* s, double r) {
    double dr = hypot(r, s->z_r), dv = hypot(r, s->z_v);
    double er = exp(-s->sigma_tr * dr), ev = exp(-s->sigma_tr * dv);
    return (s->z_r * (1.0 + s->sigma_tr * dr) * er / (dr * dr * dr)
          + s->z_v * (1.0 + s->sigma_tr * dv) * ev / (dv * dv * dv)) / (4.0 * PI);
}

/* Fraction of the escaping energy that exits within radius R: the radial
   integral of sss_profile(), in closed form.  Monotone in R, 0 at R = 0 and
   tending to 1 as R grows. */
static double sss_capture(const Sss* s, double R) {
    if (!(R > 0.0)) return 0.0;
    double dr = hypot(R, s->z_r), dv = hypot(R, s->z_v);
    double er = exp(-s->sigma_tr * s->z_r), ev = exp(-s->sigma_tr * s->z_v);
    double fr = er - s->z_r * exp(-s->sigma_tr * dr) / dr;
    double fv = ev - s->z_v * exp(-s->sigma_tr * dv) / dv;
    double c = (fr + fv) / (er + ev);
    /* near R = 0 the two terms cancel and roundoff can nudge c out of range */
    if (c < 0.0) return 0.0;
    if (c > 1.0) return 1.0;
    return c;
}

/* Same integral by Simpson quadrature, used to check the closed form. */
static double sss_capture_numeric(const Sss* s, double R) {
    if (!(R > 0.0)) return 0.0;
    int n = 20000;
    double h = R / n, sum = 0.0;
    for (int i = 0; i <= n; i++) {
        double r = i * h;
        double w = (i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
        sum += w * 2.0 * PI * r * sss_profile(s, r);
    }
    return (sum * h / 3.0) / s->Z;
}

/* Smallest radius whose disk captures `target` of the escaping energy.
   Returns 0 for target <= 0 and -1 when no finite radius reaches target. */
static double sss_radius_for_capture(const Sss* s, double target) {
    if (!(target > 0.0)) return 0.0;
    if (target >= 1.0) return -1.0;
    double lo = 0.0, hi = s->mfp;
    while (sss_capture(s, hi) < target && hi < 1e12) hi *= 2.0;
    for (int i = 0; i < 80; i++) {
        double mid = 0.5 * (lo + hi);
        if (sss_capture(s, mid) < target) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

/* Energy bookkeeping for a finite capture radius.  The three terms always sum
   to 1: `profile` is the escaped energy resolved inside R, `diffuse` is the
   escaped energy emitted as a Lambertian lobe (no lateral spread resolved),
   and `absorbed` is the loss.  At R = 0 nothing is resolved laterally and the
   material is purely diffuse; as R grows the diffuse term vanishes. */
typedef struct { double diffuse, profile, absorbed; } SssRefl;

static SssRefl sss_reflectance(const Sss* s, double R) {
    SssRefl o;
    double c = sss_capture(s, R);
    o.profile = s->escaped * c;
    o.diffuse = s->escaped * (1.0 - c);
    o.absorbed = 1.0 - s->escaped;
    return o;
}

/* BRDF of the degenerate case: no diffusion radius means no lateral
   transport, so the whole escaped energy leaves as a Lambertian lobe. */
static double sss_diffuse_brdf(const Sss* s) { return s->escaped / PI; }

/* Deterministic xorshift64 so the Monte Carlo check is reproducible. */
static uint64_t sss_rng = 0x9E3779B97F4A7C15ULL;
static double sss_uniform(void) {
    sss_rng ^= sss_rng << 13;
    sss_rng ^= sss_rng >> 7;
    sss_rng ^= sss_rng << 17;
    return (double)(sss_rng >> 11) * (1.0 / 9007199254740992.0);
}

/* Monte Carlo reference for the same transport: a point source at depth z_r
   emitting isotropically, isotropic scattering, analog absorption.  Returns
   the fraction of escaping paths with exit radius <= radii[i] and folds the
   four-sigma noise of each fraction into noise4[i].  Supports up to
   SSS_MC_MAX_RADII radii. */
#define SSS_MC_MAX_RADII 8

static void sss_mc_reference(const Sss* s, const double* radii, double* frac,
                             int n, long paths, double* escape, double* noise4) {
    double cnt[SSS_MC_MAX_RADII] = { 0 };
    long escapes = 0;
    if (n > SSS_MC_MAX_RADII) n = SSS_MC_MAX_RADII;
    sss_rng = 0x9E3779B97F4A7C15ULL;
    for (long p = 0; p < paths; p++) {
        double x = 0.0, y = 0.0, z = s->z_r;
        double u = 2.0 * sss_uniform() - 1.0;
        double phi = 2.0 * PI * sss_uniform();
        double sq = sqrt(1.0 - u * u);
        double dx = sq * cos(phi), dy = sq * sin(phi), dz = u;
        for (;;) {
            double L = -log(1.0 - sss_uniform()) * s->mfp;
            double nx = x + L * dx, ny = y + L * dy, nz = z + L * dz;
            if (nz < 0.0) {
                double t = z / (z - nz);          /* pierce the surface exactly */
                double ex = x + t * (nx - x), ey = y + t * (ny - y);
                double r = hypot(ex, ey);
                escapes++;
                for (int i = 0; i < n; i++) if (r <= radii[i]) cnt[i] += 1.0;
                break;
            }
            x = nx; y = ny; z = nz;
            if (sss_uniform() > s->albedo) break; /* absorbed */
            u = 2.0 * sss_uniform() - 1.0;
            phi = 2.0 * PI * sss_uniform();
            sq = sqrt(1.0 - u * u);
            dx = sq * cos(phi); dy = sq * sin(phi); dz = u;
        }
    }
    *escape = (double)escapes / (double)paths;
    for (int i = 0; i < n; i++) {
        frac[i] = escapes ? cnt[i] / (double)escapes : 0.0;
        noise4[i] = 4.0 * sqrt(frac[i] * (1.0 - frac[i]) / (double)(escapes ? escapes : 1));
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    check("lambert_unit", fabsf(brdf_lambert(3.14159265f) - 1.0f) < 1e-4f);
    check("lambert_zero", brdf_lambert(0.0f) == 0.0f);

    check("ggx_peak", ggx_ndf(1.0f, 0.5f) > ggx_ndf(0.2f, 0.5f));
    check("ggx_positive", ggx_ndf(0.5f, 0.3f) > 0.0f);

    check("pdf_half", fabsf(cosine_pdf(0.5f) - 0.5f / 3.14159265f) < 1e-6f);

    V3 s = cosine_sample(0.5f, 0.25f);
    check("sample_unit", fabsf(sqrtf(dot3(s, s)) - 1.0f) < 1e-4f);
    check("sample_z_nonneg", s.z >= 0.0f);

    V3 up = { 0, 0, 1 };
    check("sample_hemisphere", dot3(s, up) >= -1e-6f);

    /* mean of cos over samples ~ 2/3 for cosine sampling */
    double acc = 0.0; int n = 512;
    for (int i = 0; i < n; i++) {
        float u1 = (float)((i % 16) + 0.5f) / 16.0f;
        float u2 = (float)((i / 16) + 0.5f) / 32.0f;
        V3 v = cosine_sample(u1, u2);
        acc += v.z;
    }
    check("mean_cos_close", fabs(acc / n - 2.0 / 3.0) < 0.02);

    /* --- subsurface diffusion --- */

    /* closed form agrees with direct radial integration */
    Sss sss;
    check("sss_init_ok", sss_init(&sss, 0.9f, 1.0f) == NULL);
    double max_quad = 0.0;
    double qr[4] = { 0.5, 1.0, 2.0, 4.0 };
    for (int i = 0; i < 4; i++) {
        double a = sss_capture(&sss, qr[i]);
        double b = sss_capture_numeric(&sss, qr[i]);
        if (fabs(a - b) > max_quad) max_quad = fabs(a - b);
    }
    check("sss_capture_quad", max_quad < 1e-6);

    check("sss_capture_zero", sss_capture(&sss, 0.0) == 0.0);
    check("sss_capture_monotone", sss_capture(&sss, 1.0) < sss_capture(&sss, 2.0) &&
          sss_capture(&sss, 2.0) < sss_capture(&sss, 4.0));
    check("sss_capture_total", sss_capture(&sss, 1e6) > 0.9999 && sss_capture(&sss, 1e6) <= 1.0);

    /* the profile depends on R only through R/mfp */
    Sss s2;
    sss_init(&s2, 0.9f, 2.5f);
    check("sss_scale_with_mfp", fabs(sss_capture(&s2, 5.0) - sss_capture(&sss, 2.0)) < 1e-9);

    /* energy bookkeeping always accounts for everything exactly once */
    for (int i = 0; i < 4; i++) {
        SssRefl r = sss_reflectance(&sss, qr[i] * 3.0);
        check("sss_energy_conserved", fabs(r.diffuse + r.profile + r.absorbed - 1.0) < 1e-9);
    }
    check("sss_loss_positive", sss_reflectance(&sss, 100.0).absorbed > 0.0);

    /* R = 0 is exactly the Lambertian limit */
    SssRefl r0 = sss_reflectance(&sss, 0.0);
    check("sss_zero_radius_diffuse",
          fabs(r0.profile) < 1e-12 && fabs(r0.diffuse - (double)sss.escaped) < 1e-12 &&
          fabs(sss_diffuse_brdf(&sss) - (double)sss.escaped / 3.14159265) < 1e-9 &&
          fabs(sss_radius_for_capture(&sss, 0.0)) < 1e-12);

    /* radius is monotone in the capture target and vanishes with the target */
    double ra = sss_radius_for_capture(&sss, 0.5), rb = sss_radius_for_capture(&sss, 0.8);
    check("sss_radius_monotone", ra > 0.0 && rb > ra);
    check("sss_radius_roundtrip", fabs(sss_capture(&sss, ra) - 0.5) < 1e-9);

    /* invalid parameters name the parameter */
    Sss bad;
    check("sss_bad_mfp_zero", strcmp(sss_init(&bad, 0.9f, 0.0f), "mfp") == 0);
    check("sss_bad_mfp_neg", strcmp(sss_init(&bad, 0.9f, -1.0f), "mfp") == 0);
    check("sss_bad_mfp_nan", strcmp(sss_init(&bad, 0.9f, NAN), "mfp") == 0);
    check("sss_bad_albedo", strcmp(sss_init(&bad, 1.5f, 1.0f), "albedo") == 0);

    /* error bound against Monte Carlo.  The diffusion model's cumulative exit
       distribution must stay within its calibration error plus 4-sigma MC
       noise at every radius.  Measured for this half-space at a=0.9: the
       cumulative fraction is within 0.009 (max, around R = 1.2 mfp) and never
       more than 0.006 outside MC noise; the 0.02 bound is a safe upper bound,
       not a tight match. */
    double radii[6] = { 0.5, 1.0, 2.0, 3.0, 4.0, 6.0 };
    double frac[6], noise4[6], escape;
    sss_mc_reference(&sss, radii, frac, 6, 200000, &escape, noise4);
    double max_err = 0.0, max_excess = 0.0;
    for (int i = 0; i < 6; i++) {
        double model = sss_capture(&sss, radii[i]);
        double err = fabs(model - frac[i]);
        if (err > max_err) max_err = err;
        double excess = err - noise4[i];
        if (excess > max_excess) max_excess = excess;
    }
    check("sss_mc_error_within_bound", max_excess < 0.02);

    /* total escape is the model's weakest quantity: it under-predicts by
       ~1.8% (absolute) at a=0.9 and ~0.9% at a=0.95, and by ~5.7% absolute
       at a=0.5, where the diffusion limit itself breaks down.  Bound it
       tightly where the model is valid and only loosely elsewhere rather
       than hiding the breakdown behind one loose tolerance. */
    check("sss_mc_total_valid", fabs(escape - sss.escaped) < 0.02);
    Sss slow;
    sss_init(&slow, 0.5f, 1.0f);
    double slow_frac[6], slow_noise[6], slow_escape;
    sss_mc_reference(&slow, radii, slow_frac, 6, 200000, &slow_escape, slow_noise);
    check("sss_mc_total_diffuse_limit", fabs(slow_escape - slow.escaped) < 0.06);

    printf("sss_mc_max_abs_err=%.4f sss_mc_max_excess=%.4f\n", max_err, max_excess);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
