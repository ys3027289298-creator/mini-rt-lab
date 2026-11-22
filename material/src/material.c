
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

/* ---- multi-layer compositing -------------------------------------------
   Each layer carries an albedo and a coverage (alpha).  The composite is the
   coverage-weighted mean of the layer albedos.  When the coverages sum above
   one they are scaled down proportionally and the scale is recorded, so the
   reported reflectance can never exceed one.  A separate per-layer "over"
   reference and an error bound let the caller check the approximation. */

#define LAYER_MAX 8

typedef struct {
    float albedo[3];
    float coverage;
} Layer;

typedef struct {
    float out[3];              /* composited albedo, each channel in [0,1] */
    float coverage_sum;        /* raw sum of coverages before any scaling */
    float norm_scale;          /* 1/sum when normalized, else 1 */
    int   normalized;          /* 1 if coverage_sum > 1 and was rescaled */
    int   layer_count;         /* number of layers considered */
    float layer_cov[LAYER_MAX];/* per-layer coverage actually used */
} LayerReport;

/* total reflectance proxy: sum of the absolute channel values */
static float energy(const float v[3]) {
    return fabsf(v[0]) + fabsf(v[1]) + fabsf(v[2]);
}

/* Composite `count` layers.  Returns 0 on invalid input (bad coverage,
   out-of-range albedo, NaN, empty or overflowing layer list), else 1. */
static int layers_composite(const Layer* layers, int count, LayerReport* r) {
    for (int k = 0; k < 3; k++) r->out[k] = 0.0f;
    r->coverage_sum = 0.0f;
    r->norm_scale = 1.0f;
    r->normalized = 0;
    r->layer_count = count;
    for (int i = 0; i < LAYER_MAX; i++) r->layer_cov[i] = 0.0f;

    if (count < 1 || count > LAYER_MAX) return 0;

    Layer tmp[LAYER_MAX];
    float covsum = 0.0f;
    for (int i = 0; i < count; i++) {
        float cov = layers[i].coverage;
        /* written so NaN fails the range test; large coverages are allowed
           and handled by proportional normalization below */
        if (!(cov >= 0.0f) || cov > 1.0e30f) return 0;
        for (int k = 0; k < 3; k++) {
            float a = layers[i].albedo[k];
            if (!(a >= 0.0f) || a > 1.0f) return 0;
        }
        tmp[i] = layers[i];
        covsum += cov;
    }

    r->coverage_sum = covsum;
    if (covsum > 1.0f) {
        /* proportional normalization; record the scale we applied */
        r->normalized = 1;
        r->norm_scale = 1.0f / covsum;
        for (int i = 0; i < count; i++) tmp[i].coverage *= r->norm_scale;
    }
    for (int i = 0; i < count; i++) r->layer_cov[i] = tmp[i].coverage;

    /* accumulate in double; coverage-weighted mean of albedo */
    double acc[3] = { 0.0, 0.0, 0.0 };
    double total = 0.0;
    for (int i = 0; i < count; i++) {
        double cov = (double)tmp[i].coverage;
        for (int k = 0; k < 3; k++) acc[k] += (double)tmp[i].albedo[k] * cov;
        total += cov;
    }

    for (int k = 0; k < 3; k++) {
        double v = (total >= 1.0) ? acc[k] / total : acc[k];
        float f = (float)v;
        if (f < 0.0f) f = 0.0f;
        if (f > 1.0f) f = 1.0f;
        r->out[k] = f;
    }
    return 1;
}

/* Analytic per-layer "over" composite: layer i contributes
   coverage_i * prod_{j<i}(1 - coverage_j) of its albedo.  No normalization. */
static void layers_reference(const Layer* layers, int count, float out[3]) {
    double w[LAYER_MAX];
    double rem = 1.0;
    for (int i = 0; i < count; i++) {
        double cov = (double)layers[i].coverage;
        w[i] = rem * cov;
        rem *= (1.0 - cov);
    }
    for (int k = 0; k < 3; k++) {
        double s = 0.0;
        for (int i = 0; i < count; i++) s += w[i] * (double)layers[i].albedo[k];
        out[k] = (float)s;
    }
}

/* Upper bound on |composite - reference| per channel.  Both are coverage
   weighted means of the same albedos, so the difference is bounded by the
   extreme albedos times the accumulated positive / negative coverage
   discrepancy. */
static void layers_error_bound(const Layer* layers, int count, float bound[3]) {
    const double eps = 3.0 / 1048576.0; /* absorbs float rounding */
    double cov[LAYER_MAX];
    double sum = 0.0;
    for (int i = 0; i < count; i++) {
        cov[i] = (double)layers[i].coverage;
        sum += cov[i];
    }
    if (sum >= 1.0) {
        for (int i = 0; i < count; i++) cov[i] /= sum;
    }

    double s_pos = 0.0, s_neg = 0.0, rem = 1.0;
    for (int i = 0; i < count; i++) {
        double raw = (double)layers[i].coverage;
        double d = cov[i] - rem * raw; /* normalized - over weight */
        if (d >= 0.0) s_pos += d;
        else          s_neg += d;
        rem *= (1.0 - raw);
    }

    for (int k = 0; k < 3; k++) {
        double lo = (double)layers[0].albedo[k];
        double hi = lo;
        for (int i = 1; i < count; i++) {
            double a = (double)layers[i].albedo[k];
            if (a < lo) lo = a;
            if (a > hi) hi = a;
        }
        double t1 = hi * s_pos + lo * s_neg;
        double t2 = lo * s_pos + hi * s_neg;
        double e = fmax(fabs(t1), fabs(t2)) + eps;
        bound[k] = (float)e;
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

    /* ---- multi-layer compositing ---- */
    LayerReport r;

    Layer single = { { 0.2f, 0.5f, 0.9f }, 1.0f };
    check("layer_single_ok", layers_composite(&single, 1, &r) == 1);
    check("layer_single_exact",
          r.out[0] == 0.2f && r.out[1] == 0.5f && r.out[2] == 0.9f);
    check("layer_single_no_norm", !r.normalized && r.norm_scale == 1.0f);

    Layer part = { { 0.4f, 0.4f, 0.4f }, 0.25f };
    layers_composite(&part, 1, &r);
    check("layer_partial_blend", fabsf(r.out[0] - 0.1f) < 1e-6f);

    Layer three[3] = {
        { { 0.90f, 0.90f, 0.90f }, 0.80f },
        { { 0.85f, 0.85f, 0.85f }, 0.80f },
        { { 0.80f, 0.80f, 0.80f }, 0.80f },
    };
    check("layer_three_ok", layers_composite(three, 3, &r) == 1);
    check("layer_three_normalized", r.normalized == 1);
    check("layer_three_scale", fabsf(r.norm_scale - 0.416667f) < 1e-4f);
    check("layer_three_sum_recorded", fabsf(r.coverage_sum - 2.4f) < 1e-5f);
    check("layer_three_no_white", r.out[0] < 1.0f && r.out[1] < 1.0f && r.out[2] < 1.0f);
    check("layer_naive_overshot",
          (three[0].albedo[0] * three[0].coverage +
           three[1].albedo[0] * three[1].coverage +
           three[2].albedo[0] * three[2].coverage) > 1.0f);
    check("layer_energy_nonincrease", energy((float[3]){ 1, 1, 1 }) - energy(r.out) >= -1e-6f);
    check("layer_channels_clamped", r.out[0] <= 1.0f && r.out[1] <= 1.0f && r.out[2] <= 1.0f);
    check("layer_report_cov",
          fabsf(r.layer_cov[0] + r.layer_cov[1] + r.layer_cov[2] - 1.0f) < 1e-5f);

    float ref[3], bound[3];
    layers_reference(three, 3, ref);
    layers_error_bound(three, 3, bound);
    int within = 1, differs = 0;
    for (int k = 0; k < 3; k++) {
        if (fabsf(r.out[k] - ref[k]) > bound[k]) within = 0;
        if (fabsf(r.out[k] - ref[k]) > 1e-3f) differs = 1;
    }
    check("layer_ref_within_bound", within);
    check("layer_ref_differs", differs);
    check("layer_bound_positive", bound[0] > 0.0f);
    check("layer_bound_sane", bound[0] < 3.40282347e38f);

    Layer low[2] = {
        { { 1.0f, 1.0f, 1.0f }, 0.3f },
        { { 0.0f, 0.0f, 0.0f }, 0.2f },
    };
    check("layer_low_ok", layers_composite(low, 2, &r) == 1);
    check("layer_low_not_normalized", !r.normalized);
    check("layer_low_keeps_remainder", fabsf(r.out[0] - 0.3f) < 1e-6f);

    Layer bad_cov[1] = { { { 0.5f, 0.5f, 0.5f }, -0.1f } };
    Layer bad_alb[1] = { { { 1.5f, 0.5f, 0.5f }, 0.5f } };
    Layer nan_cov[1] = { { { 0.5f, 0.5f, 0.5f }, NAN } };
    check("layer_reject_neg_cov", layers_composite(bad_cov, 1, &r) == 0);
    check("layer_reject_high_albedo", layers_composite(bad_alb, 1, &r) == 0);
    check("layer_reject_nan", layers_composite(nan_cov, 1, &r) == 0);
    check("layer_reject_empty", layers_composite(bad_cov, 0, &r) == 0);
    check("layer_reject_overflow", layers_composite(bad_cov, LAYER_MAX + 1, &r) == 0);

    Layer reordered[3] = { three[2], three[0], three[1] };
    LayerReport r2, r3;
    layers_composite(three, 3, &r3);
    layers_composite(reordered, 3, &r2);
    check("layer_blend_order_free",
          fabsf(r3.out[0] - r2.out[0]) < 1e-6f &&
          fabsf(r3.out[1] - r2.out[1]) < 1e-6f &&
          fabsf(r3.out[2] - r2.out[2]) < 1e-6f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
