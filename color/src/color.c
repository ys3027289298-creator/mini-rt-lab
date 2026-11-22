
/* color: sRGB color decoding consistency for texture sampling.
   - sRGB EOTF/OETF per IEC 61966-2-1 (analytic float + double reference)
   - 8-bit LUT sampling path, checked per pixel against the analytic decoder
   - bilinear texture blending in linear space vs. the dark gamma-space blend
   - generic gamma transfer that rejects illegal gamma and reports the value
   - round-trip errors with explicit upper bounds */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define LUT_N 256

/* ---------------- analytic sRGB transfer ---------------- */

/* EOTF: nonlinear [0,1] -> linear [0,1] */
static float srgb_decode(float c) {
    if (c <= 0.04045f) return c / 12.92f;
    return powf((c + 0.055f) / 1.055f, 2.4f);
}

/* OETF: linear [0,1] -> nonlinear [0,1], exact inverse of srgb_decode */
static float srgb_encode(float l) {
    if (l <= 0.0031308f) return 12.92f * l;
    return 1.055f * powf(l, 1.0f / 2.4f) - 0.055f;
}

/* same spec in double precision, independent reference for the float path */
static double srgb_decode_ref(double c) {
    if (c <= 0.04045) return c / 12.92;
    return pow((c + 0.055) / 1.055, 2.4);
}

/* ---------------- 8-bit texture LUT sampling path ---------------- */

static float srgb_lut[LUT_N];

static void srgb_lut_build(void) {
    for (int i = 0; i < LUT_N; i++)
        srgb_lut[i] = srgb_decode((float)i / 255.0f);
}

/* texture sampling goes through the table, never the raw byte value */
static float srgb_decode_u8(unsigned char b) { return srgb_lut[b]; }

/* ---------------- generic gamma transfer, validated ---------------- */

typedef struct { int count; float last_gamma; } GammaErr;

static int gamma_valid(float gamma) { return isfinite(gamma) && gamma > 0.0f; }

/* returns 0 on success, 1 on illegal gamma; the offending value is
   printed and handed back via the error record so callers can surface it */
static int gamma_decode(float c, float gamma, float* out, GammaErr* e) {
    if (!gamma_valid(gamma)) {
        if (e) { e->count++; e->last_gamma = gamma; }
        printf("gamma_invalid op=decode val=%g\n", (double)gamma);
        return 1;
    }
    *out = powf(c, gamma);
    return 0;
}

static int gamma_encode(float l, float gamma, float* out, GammaErr* e) {
    if (!gamma_valid(gamma)) {
        if (e) { e->count++; e->last_gamma = gamma; }
        printf("gamma_invalid op=encode val=%g\n", (double)gamma);
        return 1;
    }
    *out = powf(l, 1.0f / gamma);
    return 0;
}

/* ---------------- bilinear texture blending ---------------- */

typedef struct { unsigned char c[2][2]; } Tex2x2;

static void bilerp_weights(float u, float v, float w[4]) {
    w[0] = (1.0f - u) * (1.0f - v); /* c00 */
    w[1] = u * (1.0f - v);          /* c10 */
    w[2] = (1.0f - u) * v;          /* c01 */
    w[3] = u * v;                   /* c11 */
}

/* correct pipeline: decode every texel, then blend in linear space */
static float bilerp_linear(const Tex2x2* t, float u, float v) {
    float px[4] = {
        srgb_decode_u8(t->c[0][0]), srgb_decode_u8(t->c[1][0]),
        srgb_decode_u8(t->c[0][1]), srgb_decode_u8(t->c[1][1])
    };
    float w[4]; bilerp_weights(u, v, w);
    float r = 0.0f;
    for (int i = 0; i < 4; i++) r += px[i] * w[i];
    return r;
}

/* buggy pipeline kept for comparison: blend encoded bytes, decode once.
   the sRGB EOTF is convex, so decode(E[x]) <= E[decode(x)]: this comes
   out darker than the linear blend and is why textures look muddy */
static float bilerp_srgb_space(const Tex2x2* t, float u, float v) {
    float sx[4] = {
        t->c[0][0] / 255.0f, t->c[1][0] / 255.0f,
        t->c[0][1] / 255.0f, t->c[1][1] / 255.0f
    };
    float w[4]; bilerp_weights(u, v, w);
    float s = 0.0f;
    for (int i = 0; i < 4; i++) s += sx[i] * w[i];
    return srgb_decode(s);
}

/* independent blend: analytic decoder instead of the LUT path */
static float bilerp_analytic(const Tex2x2* t, float u, float v) {
    float px[4] = {
        srgb_decode(t->c[0][0] / 255.0f), srgb_decode(t->c[1][0] / 255.0f),
        srgb_decode(t->c[0][1] / 255.0f), srgb_decode(t->c[1][1] / 255.0f)
    };
    float w[4]; bilerp_weights(u, v, w);
    float r = 0.0f;
    for (int i = 0; i < 4; i++) r += px[i] * w[i];
    return r;
}

/* ---------------- test scaffolding ---------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int bits_eq(float a, float b) { return memcmp(&a, &b, sizeof(float)) == 0; }
static int same_gamma(float a, float b) {
    if (isnan(a)) return isnan(b);
    return a == b;
}
static int clamp_byte(int b) {
    if (b < 0) return 0;
    if (b > 255) return 255;
    return b;
}
static int idiff(int a, int b) { int d = a - b; return d < 0 ? -d : d; }

/* observed round-trip error must stay under these declared upper bounds */
#define RT_SRGB_FLOAT_BOUND 1e-6f   /* float encode/decode, dense grid */
#define RT_GAMMA_FLOAT_BOUND 1e-6f  /* gamma 2.2 encode/decode, dense grid */
#define RT_BYTE_BOUND 1             /* 8-bit round trip, in LSBs */

int main(void) {
    /* ---- sRGB decode: endpoints, known value, branch continuity ---- */
    check("decode_endpoints", srgb_decode(0.0f) == 0.0f &&
                              fabsf(srgb_decode(1.0f) - 1.0f) < 1e-7f);
    /* published reference: sRGB 0.5 -> linear 0.21404114 */
    check("decode_known_half", fabsf(srgb_decode(0.5f) - 0.21404114f) < 1e-6f);
    float lo = 0.04045f / 12.92f;
    float hi = powf((0.04045f + 0.055f) / 1.055f, 2.4f);
    check("decode_branch_continuous", fabsf(lo - hi) < 1e-6f);

    int mono = 1;
    for (int i = 1; i < LUT_N; i++)
        if (!(srgb_decode((float)i / 255.0f) > srgb_decode((float)(i - 1) / 255.0f)))
            mono = 0;
    check("decode_monotonic", mono);

    /* ---- float analytic decoder vs double-precision reference ---- */
    double ref_max = 0.0;
    for (int i = 0; i < LUT_N; i++) {
        double c = i / 255.0;
        ref_max = fmax(ref_max, fabs((double)srgb_decode((float)c) - srgb_decode_ref(c)));
    }
    printf("decode_ref_max_diff=%.3g\n", ref_max);
    check("decode_vs_double_ref", ref_max < 3e-7);

    /* ---- LUT sampling path: per-pixel equality with analytic decode ---- */
    srgb_lut_build();
    int exact = 1;
    float lut_max = 0.0f;
    for (int i = 0; i < LUT_N; i++) {
        float a = srgb_decode((float)i / 255.0f);
        float l = srgb_decode_u8((unsigned char)i);
        if (!bits_eq(a, l)) exact = 0;
        lut_max = fmaxf(lut_max, fabsf(a - l));
    }
    printf("lut_pixel_max_diff=%.3g\n", (double)lut_max);
    check("lut_pixel_exact", exact);

    /* ---- sRGB round trip: nonlinear -> linear -> nonlinear ---- */
    float rt_max = 0.0f;
    for (int k = 0; k <= 10000; k++) {
        float c = (float)k / 10000.0f;
        rt_max = fmaxf(rt_max, fabsf(srgb_encode(srgb_decode(c)) - c));
    }
    printf("srgb_rt_float_max=%.3g bound=%.0e\n", (double)rt_max, (double)RT_SRGB_FLOAT_BOUND);
    check("srgb_rt_float_bound", rt_max <= RT_SRGB_FLOAT_BOUND);

    int byte_worst = 0;
    for (int i = 0; i < LUT_N; i++) {
        float back = 255.0f * srgb_encode(srgb_decode((float)i / 255.0f));
        int b = clamp_byte((int)lroundf(back));
        byte_worst = (int)fmaxf((float)byte_worst, (float)idiff(b, i));
    }
    printf("srgb_rt_byte_worst_lsb=%d bound=%d\n", byte_worst, RT_BYTE_BOUND);
    check("srgb_rt_byte_bound", byte_worst <= RT_BYTE_BOUND);

    /* ---- linear-space blending matches analytic blend per pixel ---- */
    Tex2x2 tex = { .c = { { 64, 192 }, { 32, 240 } } };
    float blend_max = 0.0f;
    for (int iu = 0; iu <= 10; iu++) {
        for (int iv = 0; iv <= 10; iv++) {
            float u = (float)iu / 10.0f, v = (float)iv / 10.0f;
            blend_max = fmaxf(blend_max,
                fabsf(bilerp_linear(&tex, u, v) - bilerp_analytic(&tex, u, v)));
        }
    }
    printf("blend_analytic_max_diff=%.3g\n", (double)blend_max);
    check("blend_matches_analytic", blend_max < 1e-7f);

    /* at texel centers both pipelines must agree with the texel itself */
    check("blend_center_exact",
        bits_eq(bilerp_linear(&tex, 0.0f, 0.0f), srgb_decode_u8(tex.c[0][0])) &&
        bits_eq(bilerp_srgb_space(&tex, 0.0f, 0.0f), srgb_decode_u8(tex.c[0][0])));

    /* the gamma-space blend is strictly darker in the interior (convex EOTF) */
    float lin = bilerp_linear(&tex, 0.5f, 0.5f);
    float dark = bilerp_srgb_space(&tex, 0.5f, 0.5f);
    printf("blend_mid_linear=%.6f blend_mid_gammaspace=%.6f ratio=%.4f\n",
           (double)lin, (double)dark, (double)(dark / lin));
    check("gammaspace_blend_darker", dark < lin);
    int convex_holds = 1;
    for (int iu = 1; iu < 10; iu++) {
        for (int iv = 1; iv < 10; iv++) {
            float u = (float)iu / 10.0f, v = (float)iv / 10.0f;
            if (bilerp_srgb_space(&tex, u, v) > bilerp_linear(&tex, u, v))
                convex_holds = 0;
        }
    }
    check("gammaspace_never_brighter", convex_holds);

    /* ---- generic gamma: legal path and round-trip bound ---- */
    float g_out;
    GammaErr ge; memset(&ge, 0, sizeof(ge));
    check("gamma_valid_path", gamma_decode(1.0f, 2.2f, &g_out, &ge) == 0 &&
                              fabsf(g_out - 1.0f) < 1e-7f && ge.count == 0);

    float grt_max = 0.0f;
    for (int k = 0; k <= 10000; k++) {
        float c = (float)k / 10000.0f, dl, el;
        gamma_decode(c, 2.2f, &dl, NULL);
        gamma_encode(dl, 2.2f, &el, NULL);
        grt_max = fmaxf(grt_max, fabsf(el - c));
    }
    printf("gamma22_rt_float_max=%.3g bound=%.0e\n", (double)grt_max, (double)RT_GAMMA_FLOAT_BOUND);
    check("gamma_rt_float_bound", grt_max <= RT_GAMMA_FLOAT_BOUND);

    int gbyte_worst = 0;
    for (int i = 0; i < LUT_N; i++) {
        float c = (float)i / 255.0f, dl, el;
        gamma_decode(c, 2.2f, &dl, NULL);
        gamma_encode(dl, 2.2f, &el, NULL);
        int b = clamp_byte((int)lroundf(255.0f * el));
        gbyte_worst = (int)fmaxf((float)gbyte_worst, (float)idiff(b, i));
    }
    printf("gamma22_rt_byte_worst_lsb=%d bound=%d\n", gbyte_worst, RT_BYTE_BOUND);
    check("gamma_rt_byte_bound", gbyte_worst <= RT_BYTE_BOUND);

    /* ---- illegal gamma must be rejected and its value reported ---- */
    float bad_vals[4] = { 0.0f, -2.2f, NAN, INFINITY };
    int report_ok = 1;
    for (int i = 0; i < 4; i++) {
        float tmp;
        int rd = gamma_decode(0.5f, bad_vals[i], &tmp, &ge);
        int re = gamma_encode(0.5f, bad_vals[i], &tmp, &ge);
        if (rd != 1 || re != 1 || !same_gamma(ge.last_gamma, bad_vals[i]))
            report_ok = 0;
    }
    printf("gamma_invalid_reports=%d values_tried=%d\n", ge.count, 8);
    check("gamma_invalid_reports_value", report_ok && ge.count == 8);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
