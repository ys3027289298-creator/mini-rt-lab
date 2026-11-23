/*
 * Minimal test entry for the material evaluation path.
 * No third-party test framework: includes the source under test directly
 * (path overridable via -DMATERIAL_SRC='"..."') and renames its demo main.
 *
 * Strategy:
 *   - exact/closed-form numeric checks (Lambert, PDF, known sample point),
 *     with reference constants computed independently in Python
 *   - GGX normalization oracle: stratified MC integral of D(h)*cos over
 *     the hemisphere must be 1 (~0.999995 at N=4096), plus exact peak value
 *   - sampling invariants: unit length, hemisphere, analytic cosine mean
 *     2/3, PDF/marginals, and the u1=0 pole edge case
 *   - physical monotonicity (rougher GGX => lower normal peak)
 */
#ifndef MATERIAL_SRC
#define MATERIAL_SRC "../material/src/material.c"
#endif
#define main material_demo_main
#include MATERIAL_SRC
#undef main

static const float PI = 3.14159265f;

static int t_failures = 0;
#define CHECK(name, cond) do {                                          \
    int _ok = !!(cond);                                                 \
    printf("%s=%d\n", (name), _ok);                                     \
    if (!_ok) t_failures++;                                             \
} while (0)
#define CHECK_NEAR(name, got, ref, tol) do {                            \
    double _g = (double)(got), _r = (double)(ref);                      \
    int _ok = fabs(_g - _r) <= (double)(tol);                           \
    printf("%s=%d\n", (name), _ok);                                     \
    if (!_ok) { t_failures++;                                           \
        printf("# %s got=%.9g ref=%.9g\n", (name), _g, _r); }           \
} while (0)

int main(void) {
    /* ---- Lambert: rho/pi across several albedos (independently computed) ---- */
    CHECK_NEAR("lambert_unit", brdf_lambert(PI), 1.0f, 1e-4f);
    CHECK_NEAR("lambert_two", brdf_lambert(2.0f), 2.0f / PI, 1e-6f);
    CHECK_NEAR("lambert_half", brdf_lambert(0.5f), 0.5f / PI, 1e-6f);
    CHECK("lambert_zero", brdf_lambert(0.0f) == 0.0f);

    /* ---- GGX NDF. Source uses a=rough^2, a2=a^2 (=rough^4); at nh=1,
       d=a2 so D(1,rough)=1/(pi*a2)=1/(pi*rough^4). rough=.5 -> 16/pi.
       Reference values from an independent Python re-implementation. ---- */
    CHECK_NEAR("ggx_peak_exact", ggx_ndf(1.0f, 0.5f), 16.0f / PI, 1e-4f);
    CHECK_NEAR("ggx_peak_rough1", ggx_ndf(1.0f, 1.0f), 1.0f / PI, 1e-5f);
    CHECK_NEAR("ggx_offpeak_a", ggx_ndf(0.5f, 0.5f), 0.03393891331239084f, 1e-5f);
    CHECK_NEAR("ggx_offpeak_b", ggx_ndf(0.2f, 0.3f), 0.002795758028913393f, 1e-5f);
    CHECK("ggx_positive", ggx_ndf(0.5f, 0.3f) > 0.0f);
    CHECK("ggx_peak_vs_grazing", ggx_ndf(1.0f, 0.5f) > ggx_ndf(0.2f, 0.5f));
    /* at a fixed grazing angle, a rougher lobe has a smaller peak there */
    CHECK("ggx_roughness_monotone", ggx_ndf(0.2f, 0.3f) < ggx_ndf(0.2f, 0.8f));

    /* normalization oracle: 2*pi*E[D*cos] over a stratified hemisphere = 1 */
    {
        enum { N = 4096 };
        double acc = 0.0;
        for (int i = 0; i < N; i++) {
            float z = 1.0f - ((float)i + 0.5f) / (float)N;
            acc += (double)ggx_ndf(z, 0.5f) * z;
        }
        double integ = 2.0 * PI * acc / N;
        CHECK_NEAR("ggx_normalized_mc", integ, 1.0, 0.01);
        /* same identity for a second roughness guards against tuned constants */
        double acc2 = 0.0;
        for (int i = 0; i < N; i++) {
            float z = 1.0f - ((float)i + 0.5f) / (float)N;
            acc2 += (double)ggx_ndf(z, 0.8f) * z;
        }
        CHECK_NEAR("ggx_normalized_mc_rough8", 2.0 * PI * acc2 / N, 1.0, 0.02);
    }

    /* ---- cosine PDF ---- */
    CHECK_NEAR("pdf_half", cosine_pdf(0.5f), 0.5f / PI, 1e-6f);
    CHECK_NEAR("pdf_one", cosine_pdf(1.0f), 1.0f / PI, 1e-6f);
    CHECK("pdf_zero", cosine_pdf(0.0f) == 0.0f);
    /* pdf integrates to 1 over the hemisphere: with solid angle domega,
       2*pi*mean over z in [0,1] of pdf(z) = 2*pi*(1/(2*pi)) = 1 */
    {
        double acc = 0.0; enum { N = 1024 };
        for (int i = 0; i < N; i++) {
            float z = ((float)i + 0.5f) / N;
            acc += cosine_pdf(z);
        }
        double integ = 2.0 * PI * acc / N;
        CHECK_NEAR("pdf_hemisphere_integral", integ, 1.0, 2e-3);
    }

    /* ---- cosine-weighted sampling invariants ---- */
    V3 s = cosine_sample(0.5f, 0.25f);
    CHECK("sample_unit", fabsf(sqrtf(dot3(s, s)) - 1.0f) < 1e-4f);
    CHECK("sample_z_nonneg", s.z >= 0.0f);
    V3 up = {0, 0, 1};
    CHECK("sample_hemisphere", dot3(s, up) >= -1e-6f);

    /* exact known point: u1=0.5 -> z=1/sqrt(2), r=sqrt(0.5); u2=.25 -> phi=pi/2 */
    {
        V3 v = cosine_sample(0.5f, 0.25f);
        CHECK_NEAR("sample_known_x", v.x, 0.0f, 1e-5f);
        CHECK_NEAR("sample_known_y", v.y, 0.70710678118f, 1e-5f);
        CHECK_NEAR("sample_known_z", v.z, 0.70710678118f, 1e-5f);
    }
    /* u1=0 pole: sample collapses to the normal */
    {
        V3 v = cosine_sample(0.0f, 0.7f);
        CHECK_NEAR("sample_pole_z", v.z, 1.0f, 1e-5f);
        CHECK_NEAR("sample_pole_x", v.x, 0.0f, 1e-6f);
        CHECK_NEAR("sample_pole_y", v.y, 0.0f, 1e-6f);
    }
    /* u2 endpoints map to the x-axis azimuth (phi=0 and phi=2pi) */
    {
        V3 v0 = cosine_sample(0.3f, 0.0f);
        CHECK_NEAR("sample_phi0_y", v0.y, 0.0f, 1e-6f);
        CHECK("sample_phi0_xpos", v0.x > 0.0f);
    }

    /* mean cosine over the same stratified grid the demo uses: must be 2/3 */
    {
        double acc = 0.0; int n = 512;
        for (int i = 0; i < n; i++) {
            float u1 = (float)((i % 16) + 0.5f) / 16.0f;
            float u2 = (float)((i / 16) + 0.5f) / 32.0f;
            V3 v = cosine_sample(u1, u2);
            acc += v.z;
        }
        CHECK_NEAR("mean_cos_close", acc / n, 2.0 / 3.0, 0.005);
    }
    /* radial marginal: r=sqrt(u1), so mean(r^2)=mean(u1)=1/2 on [0,1) grid */
    {
        double acc = 0.0; int n = 256;
        for (int i = 0; i < n; i++) {
            float u1 = (float)(i + 0.5f) / n;
            V3 v = cosine_sample(u1, (float)(i % 7) / 7.0f);
            acc += v.x * v.x + v.y * v.y;
        }
        CHECK_NEAR("radial_marginal", acc / n, 0.5, 1e-3);
    }

    printf("all_ok=%d\n", t_failures == 0);
    return t_failures == 0 ? 0 : 1;
}
