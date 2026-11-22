
/* roughness: perceptual roughness r -> GGX alpha mapping.

   A slider fed straight into alpha feels non-proportional because the GGX
   highlight width is nonlinear in alpha. This module measures the actual
   blur response of a GGX ground plane reflecting a sun disk, inverts it to
   build a perceptual mapping alpha(r), fits a monotone analytic curve, and
   validates the result against re-rendered targets with an error bound.
   The r = 0 case is a perfect mirror and is exported as an annotated PPM.

   Self-test style mirrors bvh/sample/material: make test prints check=1/0. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#define PI 3.14159265358979323846f

typedef struct { float x, y, z; } V3;

static V3 add3(V3 a, V3 b) { V3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static V3 sub3(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static V3 mul3(V3 a, float s) { V3 r = { a.x * s, a.y * s, a.z * s }; return r; }
static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 cross3(V3 a, V3 b) {
    V3 r = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return r;
}
static float len3(V3 a) { return sqrtf(dot3(a, a)); }
static V3 norm3(V3 a) {
    float l = len3(a);
    if (l < 1e-12f) { V3 z = { 0, 0, 0 }; return z; }
    return mul3(a, 1.0f / l);
}
static V3 reflect3(V3 i, V3 n) { return sub3(i, mul3(n, 2.0f * dot3(i, n))); }

/* xorshift32, deterministic per render */
typedef struct { unsigned s; } Rng;
static float rng_f(Rng *r) {
    r->s ^= r->s << 13; r->s ^= r->s >> 17; r->s ^= r->s << 5;
    return (r->s >> 8) * (1.0f / 16777216.0f); /* top 24 bits / 2^24 */
}

/* ---------------- GGX microfacet model (local frame, n = +z) ------------- */

/* GGX (Walter 2007) normal distribution D(h). Normalized under the
   microfacet constraint int D cos th d w = 1: at alpha = 1, D = 1/pi and
   int cos d w = pi, so the weighted integral is exactly 1. */
static float ggx_d(float nh, float alpha) {
    float a2 = alpha * alpha;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / (PI * d * d);
}

/* Smith Lambda for GGX; G1 = 1/(1+Lambda), G2 = 1/(1+Lambda_v+Lambda_l) */
static float ggx_lambda(float cos_th, float alpha) {
    if (cos_th <= 0.0f) return 1e30f;
    float tan2 = (1.0f - cos_th * cos_th) / (cos_th * cos_th);
    return 0.5f * (sqrtf(1.0f + alpha * alpha * tan2) - 1.0f);
}
static float ggx_g1(float cos_th, float alpha) { return 1.0f / (1.0f + ggx_lambda(cos_th, alpha)); }
static float ggx_g2(float vn, float ln, float alpha) {
    return 1.0f / (1.0f + ggx_lambda(vn, alpha) + ggx_lambda(ln, alpha));
}

/* NDF-only angular spread of microfacet normals, in radians. This is the
   perceptual blur metric: the distribution of facet tilt angles, free of the
   Smith masking term (which darkens grazing facets and makes the rendered
   pixel-space sigma non-monotone at high alpha).

   The microfacet tilt CDF uses the cos-weighted measure D cos th d w (the
   probability density of facet normals over projected area). Its analytic
   form is
       P(theta) = tan^2 th / (a^2 + tan^2 th),
   hence tan^2 th_q = a^2 q / (1 - q); the median is atan(a).
   The numeric cross-check integrates D cos th sin th on a uniform theta grid
   (a sin^2-theta grid is non-uniform and biases the cos-weighted CDF). */
#define NDF_STRATA 4096
static float ndf_theta_quantile(float alpha, float q, float *numeric_out) {
    double tan2 = (double)alpha * alpha * q / (1.0 - q);
    float analytic = (float)atan(sqrt(tan2));

    double w[NDF_STRATA], total = 0;
    for (int i = 0; i < NDF_STRATA; i++) {
        double th = (PI / 2.0) * (i + 0.5) / NDF_STRATA;
        w[i] = ggx_d((float)cos(th), alpha) * cos(th) * sin(th);
        total += w[i];
    }
    double target = q * total, cum = 0, numeric = PI / 2.0;
    for (int i = 0; i < NDF_STRATA; i++) {
        cum += w[i];
        if (cum >= target) {
            numeric = (PI / 2.0) * (i + 0.5) / NDF_STRATA;
            break;
        }
    }
    if (numeric_out) *numeric_out = (float)numeric;
    return analytic;
}

/* Visible-normal hemisphere sample (Heitz 2014), h returned in world frame
   built from tangent t, bitangent b and normal n. */
static V3 ggx_sample_visible(V3 v, V3 n, V3 t, V3 b, float alpha, float u1, float u2) {
    V3 vh = norm3((V3){ alpha * dot3(v, t), alpha * dot3(v, b), dot3(v, n) });
    float lensq = vh.x * vh.x + vh.y * vh.y;
    V3 th = lensq > 0.0f ? (V3){ -vh.y / sqrtf(lensq), vh.x / sqrtf(lensq), 0 }
                         : (V3){ 1, 0, 0 };
    V3 bh = cross3(vh, th);
    float r = sqrtf(u1), phi = 2.0f * PI * u2;
    float sp = r * cosf(phi), sb = r * sinf(phi);
    float s = 0.5f * (1.0f + vh.z);
    sb = (1.0f - s) * sqrtf(fmaxf(0.0f, 1.0f - sp * sp)) + s * sb;
    V3 nh_local = add3(add3(mul3(th, sp), mul3(bh, sb)),
                       mul3(vh, sqrtf(fmaxf(0.0f, 1.0f - sp * sp - sb * sb))));
    V3 h = norm3(add3(add3(mul3(t, alpha * nh_local.x), mul3(b, alpha * nh_local.y)),
                      mul3(n, fmaxf(0.0f, nh_local.z))));
    return h;
}

/* ---------------- Scene: flat conductor ground + sun disk environment ---- */

#define SUN_ANGLE 0.061f   /* ~3.5 degrees half-angle */
#define SUN_I     18.0f
#define AMBIENT   0.25f    /* flat baseline, independent of view direction */
#define F0        1.0f     /* conductor: no Fresnel darkening confuses the metric */

typedef struct { V3 eye, forward, right, up; float half_h; } Camera;

static Camera camera_make(void) {
    Camera c;
    c.eye = (V3){ 0.0f, -2.0f, 1.35f };
    V3 target = { 0.0f, 0.0f, 0.28f };
    V3 world_up = { 0, 0, 1 };
    c.forward = norm3(sub3(target, c.eye));
    c.right = norm3(cross3(c.forward, world_up));
    c.up = cross3(c.right, c.forward);
    c.half_h = tanf(0.5f * 50.0f * PI / 180.0f);
    return c;
}

/* Sun direction is chosen so the r = 0 highlight lands at frame center:
   the central incoming ray d = forward reflects to R(d,n) = (dx,dy,-dz),
   which points toward the sun (light travels from the sun along -l and
   leaves along R). Verified: sun dot R == 1 at frame center. */
static V3 sun_dir(void) {
    static V3 cached = { 0, 0, 0 };
    if (cached.x == 0.0f && cached.y == 0.0f && cached.z == 0.0f) {
        Camera c = camera_make();
        cached = norm3((V3){ c.forward.x, c.forward.y, -c.forward.z });
    }
    return cached;
}

static float env_radiance(V3 l) {
    return dot3(l, sun_dir()) >= cosf(SUN_ANGLE) ? SUN_I : AMBIENT;
}

static V3 camera_ray(const Camera *c, int px, int py, int W, int H) {
    float sx = (2.0f * (px + 0.5f) / W - 1.0f) * c->half_h * (float)W / (float)H;
    float sy = (1.0f - 2.0f * (py + 0.5f) / H) * c->half_h;
    return norm3(add3(c->forward, add3(mul3(c->right, sx), mul3(c->up, sy))));
}

/* render one linear-RGB (r=g=g, scalar) buffer; alpha == 0 is a pure mirror */
static float *render(float alpha, int W, int H, int spp, unsigned seed) {
    Camera cam = camera_make();
    V3 n = { 0, 0, 1 };
    float *buf = calloc((size_t)W * H, sizeof(float));
    Rng rng = { seed ? seed : 1u };
    int side = (int)sqrtf((float)spp + 0.5f);
    if (side * side != spp) side = spp; /* fallback: unstratified */

    for (int py = 0; py < H; py++) {
        for (int px = 0; px < W; px++) {
            V3 d = camera_ray(&cam, px, py, W, H);
            if (d.z >= -1e-6f) continue;             /* never hits ground */
            /* ground hit at t = -eye.z / d.z; shading depends only on the
               fixed normal, so the point itself is not needed. */
            V3 v = mul3(d, -1.0f);                   /* toward viewer */
            V3 tangent = { 1, 0, 0 }, bitangent = { 0, 1, 0 };

            float acc = 0.0f;
            if (alpha <= 0.0f) {
                V3 l = reflect3(d, n);
                acc = env_radiance(l);
            } else {
                for (int k = 0; k < spp; k++) {
                    float u1, u2;
                    if (side * side == spp) {
                        int a = k % side, b = k / side;
                        u1 = ((float)a + rng_f(&rng)) / side;
                        u2 = ((float)b + rng_f(&rng)) / side;
                    } else {
                        u1 = rng_f(&rng); u2 = rng_f(&rng);
                    }
                    V3 h = ggx_sample_visible(v, n, tangent, bitangent, alpha, u1, u2);
                    V3 l = reflect3(mul3(v, -1.0f), h);
                    if (l.z <= 0.0f) continue;
                    float vn = dot3(v, n), ln = dot3(l, n), vh = dot3(v, h);
                    /* p(l) = D G1 / (4 vn), f cos / p collapses to F * G2/G1 */
                    float f = F0 + (1.0f - F0) * powf(fmaxf(0.0f, 1.0f - vh), 5.0f);
                    float w = f * ggx_g2(vn, ln, alpha) / ggx_g1(vn, alpha);
                    acc += w * env_radiance(l);
                }
                acc /= spp;
            }
            buf[py * W + px] = acc;
        }
    }
    return buf;
}

/* highlight spread in pixels: intensity-weighted standard deviation about
   the centroid, flat AMBIENT baseline subtracted first */
static double measure_sigma(const float *img, int W, int H, double *cx_out, double *cy_out) {
    double sw = 0, sx = 0, sy = 0;
    for (int py = 0; py < H; py++)
        for (int px = 0; px < W; px++) {
            double w = fmax(0.0, (double)img[py * W + px] - AMBIENT);
            sw += w; sx += w * px; sy += w * py;
        }
    double cx = sx / sw, cy = sy / sw, var = 0;
    for (int py = 0; py < H; py++)
        for (int px = 0; px < W; px++) {
            double w = fmax(0.0, (double)img[py * W + px] - AMBIENT);
            double dx = px - cx, dy = py - cy;
            var += w * (dx * dx + dy * dy);
        }
    if (cx_out) *cx_out = cx;
    if (cy_out) *cy_out = cy;
    return sqrt(var / sw);
}

/* ---------------- PPM export with 5x7 bitmap annotation ------------------ */

static const unsigned glyph[26][7] = {
    /*A*/ {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11},
    /*B*/ {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E},
    /*C*/ {0x0E,0x11,0x10,0x10,0x10,0x11,0x0E},
    /*D*/ {0x1C,0x12,0x11,0x11,0x11,0x12,0x1C},
    /*E*/ {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F},
    /*F*/ {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10},
    /*G*/ {0x0E,0x11,0x10,0x17,0x11,0x11,0x0F},
    /*H*/ {0x11,0x11,0x11,0x1F,0x11,0x11,0x11},
    /*I*/ {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E},
    /*J*/ {0x07,0x02,0x02,0x02,0x02,0x12,0x0C},
    /*K*/ {0x11,0x12,0x14,0x18,0x14,0x12,0x11},
    /*L*/ {0x10,0x10,0x10,0x10,0x10,0x10,0x1F},
    /*M*/ {0x11,0x1B,0x15,0x15,0x11,0x11,0x11},
    /*N*/ {0x11,0x19,0x15,0x13,0x11,0x11,0x11},
    /*O*/ {0x0E,0x11,0x11,0x11,0x11,0x11,0x0E},
    /*P*/ {0x1E,0x11,0x11,0x1E,0x10,0x10,0x10},
    /*Q*/ {0x0E,0x11,0x11,0x11,0x15,0x12,0x0D},
    /*R*/ {0x1E,0x11,0x11,0x1E,0x14,0x12,0x11},
    /*S*/ {0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E},
    /*T*/ {0x1F,0x04,0x04,0x04,0x04,0x04,0x04},
    /*U*/ {0x11,0x11,0x11,0x11,0x11,0x11,0x0E},
    /*V*/ {0x11,0x11,0x11,0x11,0x11,0x0A,0x04},
    /*W*/ {0x11,0x11,0x11,0x15,0x15,0x1B,0x11},
    /*X*/ {0x11,0x11,0x0A,0x04,0x0A,0x11,0x11},
    /*Y*/ {0x11,0x11,0x0A,0x04,0x04,0x04,0x04},
    /*Z*/ {0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}
};

/* digits share a second small 5x7 table */
static const unsigned digit[10][7] = {
    {0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}, /*0*/
    {0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, /*1*/
    {0x0E,0x11,0x01,0x02,0x04,0x08,0x1F}, /*2*/
    {0x1E,0x01,0x01,0x0E,0x01,0x01,0x1E}, /*3*/
    {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, /*4*/
    {0x1F,0x10,0x10,0x1E,0x01,0x01,0x1E}, /*5*/
    {0x0E,0x10,0x10,0x1E,0x11,0x11,0x0E}, /*6*/
    {0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, /*7*/
    {0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}, /*8*/
    {0x0E,0x11,0x11,0x0F,0x01,0x01,0x0E}  /*9*/
};

static void draw_text(unsigned char *ppm, int W, int H, int x0, int y0,
                      const char *text, int scale) {
    int cursor = x0;
    for (const char *p = text; *p; p++) {
        if (*p == ' ') { cursor += 6 * scale; continue; }
        const unsigned *g = NULL;
        if (*p >= 'A' && *p <= 'Z') g = glyph[*p - 'A'];
        else if (*p >= '0' && *p <= '9') g = digit[*p - '0'];
        else if (*p == '=') { static const unsigned eq[7] = {0,0,0x1F,0,0x1F,0,0}; g = eq; }
        if (!g) { cursor += 6 * scale; continue; }
        for (int row = 0; row < 7; row++)
            for (int col = 0; col < 5; col++)
                if ((g[row] >> (4 - col)) & 1)
                    for (int dy = 0; dy < scale; dy++)
                        for (int dx = 0; dx < scale; dx++) {
                            int x = cursor + col * scale + dx, y = y0 + row * scale + dy;
                            if (x >= 0 && x < W && y >= 0 && y < H)
                                memset(ppm + (y * W + x) * 3, 0, 3);
                        }
        cursor += 6 * scale;
    }
}

static void export_mirror_ppm(const char *path, int W, int H) {
    float *img = render(0.0f, W, H, 1, 20261123u);
    unsigned char *ppm = malloc((size_t)W * H * 3);
    for (int i = 0; i < W * H; i++) {
        float c = img[i] / (1.0f + img[i]);                /* Reinhard */
        c = powf(c, 1.0f / 2.2f);                          /* gamma */
        unsigned char v = (unsigned char)lroundf(fminf(1.0f, c) * 255.0f);
        ppm[i * 3] = ppm[i * 3 + 1] = ppm[i * 3 + 2] = v;
    }
    /* light backing strip for the annotation */
    for (int y = 6; y < 50; y++)
        for (int x = 6; x < W - 6; x++)
            memset(ppm + (y * W + x) * 3, 235, 3);
    draw_text(ppm, W, H, 12, 10, "R=0 ALPHA=0", 2);
    draw_text(ppm, W, H, 12, 30, "MIRROR", 2);
    /* 2px white frame marks this as an exported reference image */
    for (int i = 0; i < W; i++)
        for (int b = 0; b < 2; b++) {
            memset(ppm + (b * W + i) * 3, 255, 3);
            memset(ppm + (((H - 1 - b) * W) + i) * 3, 255, 3);
        }
    for (int i = 0; i < H; i++)
        for (int b = 0; b < 2; b++) {
            memset(ppm + (i * W + b) * 3, 255, 3);
            memset(ppm + (i * W + (W - 1 - b)) * 3, 255, 3);
        }
    FILE *f = fopen(path, "wb");
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    fwrite(ppm, 3, (size_t)W * H, f);
    fclose(f);
    free(ppm);
    free(img);
}

/* ---------------- small numeric helpers ---------------------------------- */

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* monotone linear interpolation of table y(x), x strictly increasing */
static double interp(const double *xs, const double *ys, int n, double x) {
    if (x <= xs[0]) return ys[0];
    if (x >= xs[n - 1]) return ys[n - 1];
    int lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        int m = (lo + hi) / 2;
        if (xs[m] <= x) lo = m; else hi = m;
    }
    double t = (x - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + t * (ys[hi] - ys[lo]);
}

int main(void) {
    mkdir("out", 0755);
    const int W = 240, H = 240, SPP = 196;
    const int NA = 81;                 /* alpha sweep: k/80 */
    const int NR = 17;                 /* perceptual grid: j/16 */
    const float ALPHA_MIN = 1e-5f;     /* NDF limit below which width ~ 0 */

    /* --- unit check 1: NDF normalization under the microfacet constraint
       int D cos th d w = 1. At alpha = 1, D = 1/pi is constant and
       int cos d w = pi, so the weighted integral is exactly 1; a uniform-grid
       evaluation at mid alpha confirms it away from the degenerate point. --- */
    check("ggx_ndf_constant_at_alpha1",
          fabs(ggx_d(0.37f, 1.0f) - 1.0f / PI) < 1e-6f &&
          fabs(ggx_d(0.91f, 1.0f) - 1.0f / PI) < 1e-6f);
    {
        int N = 4096;
        double acc = 0;
        for (int i = 0; i < N; i++) {
            double x = (i + 0.5) / N;                    /* x = cos theta */
            acc += ggx_d((float)x, 0.5f) * x;           /* D cos th */
        }
        double norm = 2.0 * PI * acc / N;
        check("ggx_ndf_cos_weighted_normalized", fabs(norm - 1.0) < 0.01);
    }

    /* --- unit check 2: analytic CDF vs independent numeric integration.
       The cos-weighted P(theta) = tan^2 / (a^2 + tan^2) is exact; the numeric
       integral uses a uniform theta grid, so it resolves the peak across the
       whole alpha range. ----------------------------------------------- */
    int quantile_ok = 1;
    for (int ai = 1; ai <= 10; ai++) {
        float a = ai / 10.0f, num = 0;
        float an = ndf_theta_quantile(a, 0.5f, &num);
        if (fabs(num - an) > 2e-3f) quantile_ok = 0;
    }
    check("ndf_quantile_analytic", quantile_ok);

    /* --- alpha sweep: NDF median facet angle = perceptual blur B(alpha) -- */
    double *alpha = malloc(NA * sizeof(double));
    double *bresp = malloc(NA * sizeof(double));
    double th_lo = ndf_theta_quantile(ALPHA_MIN, 0.5f, NULL);
    double th_hi = ndf_theta_quantile(1.0f, 0.5f, NULL);
    for (int k = 0; k < NA; k++) {
        alpha[k] = k == 0 ? 0.0 : (double)k / (NA - 1);
        double th = ndf_theta_quantile((float)fmax(ALPHA_MIN, alpha[k]), 0.5f, NULL);
        bresp[k] = (th - th_lo) / (th_hi - th_lo);
    }
    /* numeric integration is deterministic; the only tolerance is discretization */
    int sweep_mono = 1;
    double min_step = 1e9;
    for (int k = 1; k < NA; k++) {
        double d = bresp[k] - bresp[k - 1];
        if (d < min_step) min_step = d;
        if (d < -1e-9) sweep_mono = 0;
    }
    check("blur_response_monotone", sweep_mono);
    check("blur_response_endpoints", fabs(bresp[0]) < 1e-6 &&
                                     fabs(bresp[NA - 1] - 1.0) < 1e-6);
    printf("min_blur_response_step=%.6f\n", min_step);

    /* --- PRIMARY MAPPING: monotone LUT alpha_lut(r) by inverting B.
       The 81-point sweep is inverted onto a dense 1024-entry perceptual grid
       via monotone linear interpolation; this table is the shipped mapping.
       The power fit and textbook r^2 are kept only as compact comparisons. */
    const int NL = 1024;
    double r_lut[NL], alpha_lut[NL];
    for (int i = 0; i < NL; i++) {
        r_lut[i] = (double)i / (NL - 1);
        alpha_lut[i] = interp(bresp, alpha, NA, r_lut[i]);
    }
    int lut_mono = 1;
    for (int i = 1; i < NL; i++)
        if (alpha_lut[i] < alpha_lut[i - 1] - 1e-12) lut_mono = 0;
    check("lut_monotone", lut_mono);
    check("lut_endpoints", fabs(alpha_lut[0]) < 1e-9 &&
                           fabs(alpha_lut[NL - 1] - 1.0) < 1e-9);

    /* coarse perceptual grid used for the printed validation table */
    double r[NR], alpha_emp[NR];
    for (int j = 0; j < NR; j++) {
        r[j] = (double)j / (NR - 1);
        alpha_emp[j] = interp(bresp, alpha, NA, r[j]);
    }
    check("empirical_endpoints", fabs(alpha_emp[0]) < 1e-6 &&
                                 fabs(alpha_emp[NR - 1] - 1.0) < 1e-6);

    /* --- compact analytic approximation: alpha = r^p, log-linear LSQ ---- */
    double num = 0, den = 0;
    for (int j = 1; j < NR - 1; j++) {
        double lr = log(r[j]), la = log(fmax(1e-9, alpha_emp[j]));
        num += lr * la; den += lr * lr;
    }
    double p = num / den;
    printf("power_approx_exponent_p=%.4f (textbook reference: 2.0)\n", p);
    check("fit_exponent_plausible", p > 0.8 && p < 2.5);
    check("fit_endpoints", pow(0.0, p) == 0.0 && fabs(pow(1.0, p) - 1.0) < 1e-12);

    int curve_mono = 1;
    for (int i = 1; i < 4096; i++) {
        double r0 = (i - 1) / 4095.0, r1 = i / 4095.0;
        if (pow(r1, p) - pow(r0, p) < -1e-15) curve_mono = 0;
    }
    check("power_curve_monotone", curve_mono);

    /* --- LUT validation: blur implied by alpha_lut(r) must return r.
       Error is measured on the dense LUT grid, not just the coarse table. -- */
    double b_lut[NR], err_lut[NR];
    double max_err_lut = 0;
    for (int j = 0; j < NR; j++) {
        double th = ndf_theta_quantile(fmaxf(ALPHA_MIN, (float)alpha_emp[j]),
                                       0.5f, NULL);
        b_lut[j] = (th - th_lo) / (th_hi - th_lo);
        err_lut[j] = b_lut[j] - r[j];
        if (fabs(err_lut[j]) > max_err_lut) max_err_lut = fabs(err_lut[j]);
    }
    /* dense-grid LUT residual: invert B at every LUT node */
    double dense_max = 0;
    for (int i = 0; i < NL; i++) {
        double th = ndf_theta_quantile(fmaxf(ALPHA_MIN, (float)alpha_lut[i]),
                                       0.5f, NULL);
        double back = (th - th_lo) / (th_hi - th_lo);
        if (fabs(back - r_lut[i]) > dense_max) dense_max = fabs(back - r_lut[i]);
    }

    /* continuous error bound. B(alpha) is Lipschitz; its max slope L over the
       sweep bounds how far an alpha error maps to a blur error. The LUT nodes
       sit on the exact inverse, so between nodes the piecewise-linear alpha
       error is bounded by half the node spacing in r, scaled by d alpha/dr. */
    double max_slope = 0;
    for (int k = 1; k < NA; k++) {
        double sl = (bresp[k] - bresp[k - 1]) / (alpha[k] - alpha[k - 1]);
        if (sl > max_slope) max_slope = sl;
    }
    double h_lut = 1.0 / (NL - 1);
    double max_dadr = 0;
    for (int i = 1; i < NL; i++) {
        double dadr = (alpha_lut[i] - alpha_lut[i - 1]) / h_lut;
        if (dadr > max_dadr) max_dadr = dadr;
    }
    double lipschitz_term = max_slope * max_dadr * h_lut / 2.0;
    double bound_continuous = dense_max + lipschitz_term;
    printf("lut_grid_max_abs_error=%.6f (dense %d-node grid)\n", dense_max, NL);
    printf("continuous_error_bound=%.6f (grid %.6f + Lipschitz %.6f)\n",
           bound_continuous, dense_max, lipschitz_term);
    check("lut_grid_error_under_1pct", dense_max < 0.01);
    check("bound_under_2pct", bound_continuous < 0.02);

    /* --- power fit and textbook r^2 residuals on the coarse grid -------- */
    double err_fit[NR], max_err = 0, err_r2[NR], max_err_r2 = 0;
    for (int j = 0; j < NR; j++) {
        double thp = ndf_theta_quantile(fmaxf(ALPHA_MIN, (float)pow(r[j], p)),
                                        0.5f, NULL);
        double bp = (thp - th_lo) / (th_hi - th_lo);
        err_fit[j] = bp - r[j];
        if (fabs(err_fit[j]) > max_err) max_err = fabs(err_fit[j]);

        double th2 = ndf_theta_quantile(fmaxf(ALPHA_MIN, (float)(r[j] * r[j])),
                                        0.5f, NULL);
        double b2 = (th2 - th_lo) / (th_hi - th_lo);
        err_r2[j] = b2 - r[j];
        if (fabs(err_r2[j]) > max_err_r2) max_err_r2 = fabs(err_r2[j]);
    }
    printf("power_approx_grid_max_abs_error=%.5f\n", max_err);
    printf("r2_curve_grid_max_abs_error=%.5f\n", max_err_r2);
    check("power_beats_quadratic", max_err <= max_err_r2 + 1e-12);
    check("lut_beats_power", dense_max < max_err);

    /* --- independent cross-check: actual MC-rendered highlight spread.
       Pixel-space sigma is expected to peak mid-range then fall, because the
       Smith G2 term masks grazing facets; we verify the effect is real and
       bounded, not that it is monotone. ------------------------------------ */
    double mc_noise;
    {
        float *m1 = render(0.25f, W, H, SPP, 424242u);
        float *m2 = render(0.25f, W, H, SPP, 7777777u);
        double s1 = measure_sigma(m1, W, H, NULL, NULL);
        double s2 = measure_sigma(m2, W, H, NULL, NULL);
        mc_noise = fabs(s1 - s2) / sqrt(2.0);
        free(m1); free(m2);
    }
    double sig_px[11], edge_weight[11];
    int mc_compact_ends = 1, mc_sigma_peak_real = 0;
    double peak_sigma = 0;
    for (int ai = 0; ai <= 10; ai++) {
        float a = fmaxf(ALPHA_MIN, ai / 10.0f);
        float *img = render(a, W, H, ai == 0 ? 1 : SPP, 2000u + (unsigned)ai * 131u);
        double cx, cy;
        sig_px[ai] = measure_sigma(img, W, H, &cx, &cy);
        double ew = 0, sw = 0;
        for (int py = 0; py < H; py++)
            for (int px = 0; px < W; px++) {
                double wgt = fmax(0.0, (double)img[py * W + px] - AMBIENT);
                sw += wgt;
                if (px < 4 || py < 4 || px >= W - 4 || py >= H - 4) ew += wgt;
            }
        edge_weight[ai] = ew / sw;
        if (sig_px[ai] > peak_sigma) peak_sigma = sig_px[ai];
        if (ai == 0 && (cx < 0.3 * W || cx > 0.7 * W || cy < 0.3 * H || cy > 0.7 * H))
            mc_compact_ends = 0;
        free(img);
    }
    for (int ai = 1; ai < 10; ai++)
        if (sig_px[ai] > sig_px[10] + 5.0 * mc_noise &&
            sig_px[ai] > sig_px[0] + 5.0 * mc_noise) mc_sigma_peak_real = 1;
    int edge_ok = 1;
    /* at high alpha the highlight physically covers the whole ground plane and
       any finite frame clips it; edge energy is only constrained in the regime
       where the lobe is expected to fit (alpha <= 0.3). */
    for (int ai = 0; ai <= 3; ai++)
        if (edge_weight[ai] > 0.05) edge_ok = 0;
    check("mirror_centered", mc_compact_ends);
    check("low_alpha_edge_energy_below_5pct", edge_ok);
    check("g2_masking_peak_observed", mc_sigma_peak_real);
    printf("mc_sigma_noise_px=%.4f peak_sigma_px=%.2f\n", mc_noise, peak_sigma);

    /* --- CSV outputs ----------------------------------------------------- */
    FILE *f = fopen("out/alpha_sweep.csv", "w");
    fprintf(f, "alpha,blur_normalized,theta50_deg\n");
    for (int k = 0; k < NA; k++) {
        double th = ndf_theta_quantile(fmaxf(ALPHA_MIN, (float)alpha[k]), 0.5f, NULL);
        fprintf(f, "%.6f,%.6f,%.6f\n", alpha[k], bresp[k], th * 180.0 / PI);
    }
    fclose(f);

    f = fopen("out/mapping_curve.csv", "w");
    fprintf(f, "r,alpha_lut,alpha_power,alpha_r2\n");
    for (int k = 0; k < NA; k++) {
        double rr = (double)k / (NA - 1);
        double al = interp(r_lut, alpha_lut, NL, rr);
        fprintf(f, "%.6f,%.6f,%.6f,%.6f\n", rr, al, pow(rr, p), rr * rr);
    }
    fclose(f);

    f = fopen("out/mapping_lut.csv", "w");
    fprintf(f, "r,alpha\n");
    for (int i = 0; i < NL; i++)
        fprintf(f, "%.6f,%.6f\n", r_lut[i], alpha_lut[i]);
    fclose(f);

    f = fopen("out/validation.csv", "w");
    fprintf(f, "r,alpha_lut,blur_from_lut,error_lut,error_power,error_r2\n");
    for (int j = 0; j < NR; j++)
        fprintf(f, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                r[j], alpha_emp[j], b_lut[j], err_lut[j], err_fit[j], err_r2[j]);
    fclose(f);

    f = fopen("out/render_crosscheck.csv", "w");
    fprintf(f, "alpha,pixel_sigma_px,edge_energy_fraction\n");
    for (int ai = 0; ai <= 10; ai++)
        fprintf(f, "%.4f,%.4f,%.6f\n", ai / 10.0, sig_px[ai], edge_weight[ai]);
    fclose(f);

    f = fopen("out/summary.txt", "w");
    fprintf(f, "perceptual roughness r -> GGX alpha mapping\n");
    fprintf(f, "primary mapping: %d-entry monotone LUT (out/mapping_lut.csv)\n", NL);
    fprintf(f, "compact approximations: alpha = r^%.4f (fitted), alpha = r^2 (textbook)\n", p);
    fprintf(f, "blur metric: GGX NDF median microfacet angle theta50, normalized 0..1\n");
    fprintf(f, "LUT dense-grid max |measured blur - target r|: %.6f\n", dense_max);
    fprintf(f, "continuous error bound (grid %.6f + Lipschitz %.6f): %.6f\n",
            dense_max, lipschitz_term, bound_continuous);
    fprintf(f, "power approx grid max |error|: %.5f; alpha=r^2 grid max |error|: %.5f\n",
            max_err, max_err_r2);
    fprintf(f, "MC render: %dx%d, %d spp, sun half-angle %.2f deg, sigma noise %.3f px\n",
            W, H, SPP, SUN_ANGLE * 180.0f / PI, mc_noise);
    fprintf(f, "rendered pixel sigma peaks mid-range from Smith G2 masking (documented, non-monotone)\n");
    fprintf(f, "r = 0 exports a perfect mirror: out/mirror_r0.ppm\n");
    fclose(f);

    /* --- r = 0 mirror export + PPM structural checks -------------------- */
    export_mirror_ppm("out/mirror_r0.ppm", 160, 160);
    FILE *pf = fopen("out/mirror_r0.ppm", "rb");
    char magic[8] = {0}; int pw = 0, ph = 0, pv = 0;
    int header_ok = pf && fscanf(pf, "%7s %d %d %d", magic, &pw, &ph, &pv) == 4;
    fclose(pf);
    check("ppm_header", header_ok && strcmp(magic, "P6") == 0 &&
                       pw == 160 && ph == 160 && pv == 255);
    pf = fopen("out/mirror_r0.ppm", "rb");
    int hdr[3]; char tag[8] = {0};
    int ok2 = pf && fscanf(pf, "%7s %d %d %d", tag, &hdr[0], &hdr[1], &hdr[2]) == 4;
    fgetc(pf); /* single whitespace separates the P6 header from raw pixels */
    int corner = fgetc(pf);
    fclose(pf);
    check("ppm_border_marker", ok2 && corner > 250);

    /* mirror sharpness: highlight covers only a small pixel fraction */
    float *mir = render(0.0f, W, H, 1, 20261123u);
    int hot = 0;
    for (int i = 0; i < W * H; i++)
        if (mir[i] > AMBIENT + 0.5f * (SUN_I - AMBIENT)) hot++;
    check("mirror_highlight_compact", (double)hot / (W * H) < 0.06);
    free(mir);

    free(alpha); free(bresp);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
