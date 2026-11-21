
/* material: Lambert + GGX BRDF evaluation and cosine-weighted sampling,
   with a PDF check against the analytic cosine distribution. */
#include <stdio.h>
#include <math.h>
#include <string.h>

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

/* ---- energy conservation validation ---- */

#define MAT_PI 3.14159265f
#define ENERGY_GAIN_TOL 0.02   /* reflectance above 1 + tol raises an alarm */
#define ENERGY_RATIO_TOL 0.05  /* MC integral must sit within 5% of reference */

typedef enum { MAT_LAMBERT, MAT_GGX } MatType;

typedef struct {
    MatType type;
    float albedo;    /* reflectance scale, valid range [0,1] */
    float roughness; /* GGX only, in (0,1] */
} Material;

typedef struct {
    int rejected;       /* albedo out of [0,1] or negative BRDF: refuse to use */
    int energy_alarm;   /* reflectance above 1 + ENERGY_GAIN_TOL: creates energy */
    double reflectance; /* MC estimate of integral of f*cos over the hemisphere */
    double reference;   /* analytic (Lambert) or dense-quadrature (GGX) value */
    double ratio;       /* reflectance / reference, should be ~1 */
} EnergyReport;

/* Smith masking-shadowing G1 (alpha = roughness^2, matching ggx_ndf) */
static float ggx_g1(float cos_v, float rough) {
    float a = rough * rough;
    float a2 = a * a;
    float c = cos_v > 1e-6f ? cos_v : 1e-6f;
    return 2.0f * c / (c + sqrtf(a2 + (1.0f - a2) * c * c));
}

/* full BRDF f(wi, wo) about the +z normal; GGX uses constant F = albedo */
static float brdf_eval(const Material* m, V3 wi, V3 wo) {
    if (m->type == MAT_LAMBERT) return brdf_lambert(m->albedo);
    V3 h = { wi.x + wo.x, wi.y + wo.y, wi.z + wo.z };
    h = norm3(h);
    float ci = wi.z > 1e-6f ? wi.z : 1e-6f;
    float co = wo.z > 1e-6f ? wo.z : 1e-6f;
    float d = ggx_ndf(h.z, m->roughness);
    float g = ggx_g1(ci, m->roughness) * ggx_g1(co, m->roughness);
    return m->albedo * d * g / (4.0f * ci * co);
}

/* R = integral of f*cos over the hemisphere. Cosine-weighted sampling makes
   the estimator f*cos/pdf = f*pi. Also flags any negative BRDF value. */
static double reflectance_mc(const Material* m, V3 wo, int* saw_negative) {
    double acc = 0.0;
    int n = 512;
    for (int i = 0; i < n; i++) {
        float u1 = (float)((i % 16) + 0.5f) / 16.0f;
        float u2 = (float)((i / 16) + 0.5f) / 32.0f;
        V3 wi = cosine_sample(u1, u2);
        float f = brdf_eval(m, wi, wo);
        if (f < -1e-6f) *saw_negative = 1;
        acc += (double)f * MAT_PI;
    }
    return acc / n;
}

/* reference integral: exact for Lambert, dense (theta, phi) quadrature for GGX */
static double reflectance_ref(const Material* m, V3 wo) {
    if (m->type == MAT_LAMBERT) return m->albedo; /* albedo/pi integrated against cos */
    double acc = 0.0;
    int nt = 128, np = 256;
    for (int it = 0; it < nt; it++) {
        double theta = (it + 0.5) / nt * (MAT_PI / 2.0);
        double st = sin(theta), ct = cos(theta);
        for (int ip = 0; ip < np; ip++) {
            double phi = (ip + 0.5) / np * (2.0 * MAT_PI);
            V3 wi = { (float)(st * cos(phi)), (float)(st * sin(phi)), (float)ct };
            acc += brdf_eval(m, wi, wo) * ct * st;
        }
    }
    return acc * (MAT_PI / 2.0 / nt) * (2.0 * MAT_PI / np);
}

static EnergyReport energy_check(const Material* m, V3 wo) {
    EnergyReport r;
    memset(&r, 0, sizeof(r));
    int neg = 0;
    r.reflectance = reflectance_mc(m, wo, &neg);
    if (m->albedo < 0.0f || m->albedo > 1.0f || neg) {
        r.rejected = 1;
        return r;
    }
    r.reference = reflectance_ref(m, wo);
    r.ratio = r.reference > 0.0 ? r.reflectance / r.reference : 0.0;
    if (r.reflectance > 1.0 + ENERGY_GAIN_TOL) r.energy_alarm = 1;
    return r;
}

/* mixed material: sum of weighted reflectances must not exceed 1 */
static int mix_energy_alarm(const Material* mats, const float* w, int n, V3 wo,
                            double* total_out) {
    double total = 0.0;
    for (int i = 0; i < n; i++) {
        EnergyReport r;
        if (w[i] < 0.0f) return 1;
        r = energy_check(&mats[i], wo);
        if (r.rejected) return 1;
        total += (double)w[i] * r.reflectance;
    }
    *total_out = total;
    return total > 1.0 + ENERGY_GAIN_TOL;
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

    /* energy conservation: per-material ratio of MC integral vs reference */
    V3 wo = { 0, 0, 1 };
    Material lam = { MAT_LAMBERT, 0.8f, 0.0f };
    Material ggx = { MAT_GGX, 1.0f, 0.5f };
    EnergyReport rl = energy_check(&lam, wo);
    EnergyReport rg = energy_check(&ggx, wo);
    printf("energy_ratio[lambert]=%.4f\n", rl.ratio);
    printf("energy_ratio[ggx]=%.4f\n", rg.ratio);
    printf("reflectance[lambert]=%.4f\n", rl.reflectance);
    printf("reflectance[ggx]=%.4f\n", rg.reflectance);
    check("lambert_energy_close", !rl.rejected && fabs(rl.ratio - 1.0) < ENERGY_RATIO_TOL);
    check("ggx_energy_close", !rg.rejected && fabs(rg.ratio - 1.0) < ENERGY_RATIO_TOL);
    check("lambert_no_gain", !rl.energy_alarm);
    check("ggx_no_gain", !rg.energy_alarm);

    /* rejection: albedo out of [0,1] or negative BRDF */
    Material bad_hi = { MAT_LAMBERT, 1.5f, 0.0f };
    Material bad_neg = { MAT_GGX, -0.2f, 0.5f };
    check("reject_albedo_high", energy_check(&bad_hi, wo).rejected);
    check("reject_albedo_neg", energy_check(&bad_neg, wo).rejected);

    /* mixing: weighted sum of reflectances must stay <= 1 */
    Material pair[2] = { { MAT_LAMBERT, 0.8f, 0.0f }, { MAT_GGX, 1.0f, 0.5f } };
    float w_ok[2] = { 0.5f, 0.5f };
    float w_bad[2] = { 1.0f, 1.0f };
    double mix_total = 0.0;
    check("mix_energy_ok", mix_energy_alarm(pair, w_ok, 2, wo, &mix_total) == 0);
    printf("mix_reflectance[ok]=%.4f\n", mix_total);
    check("mix_energy_alarm", mix_energy_alarm(pair, w_bad, 2, wo, &mix_total) == 1);
    printf("mix_reflectance[bad]=%.4f\n", mix_total);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
