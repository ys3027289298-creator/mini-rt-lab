
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
static float brdf_lambert(float albedo) { return albedo * 3.14159265f; }

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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
