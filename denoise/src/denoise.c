
/* denoise: edge-preserving guided filter for low-spp images.
   Weights are guided by normal and depth discontinuities; the filter
   radius scales with per-pixel variance. Quality is measured against a
   clean reference: PSNR gain and edge over-blur pixel count. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MAX_W 128
#define MAX_H 128
#define MIN_RADIUS 1
#define MAX_RADIUS 6
#define VAR_REF 0.0025f /* variance at which the radius saturates */

typedef struct { float r, g, b; } RGB;

typedef struct {
    int w, h;
    RGB color[MAX_W * MAX_H];   /* noisy input */
    RGB normal[MAX_W * MAX_H];  /* guide: shading normal, xyz in r/g/b */
    float depth[MAX_W * MAX_H]; /* guide: view depth */
    float var[MAX_W * MAX_H];   /* per-pixel luminance variance */
    int has_normal, has_depth, has_var;
} Frame;

typedef struct {
    int missing_normal, missing_depth, missing_var;
} GuideReport;

/* radius from variance: 0 var -> MIN_RADIUS, VAR_REF -> MAX_RADIUS */
static int radius_for_var(float var) {
    if (var < 0.0f) var = 0.0f;
    float t = sqrtf(var / VAR_REF);
    int r = MIN_RADIUS + (int)(t * (MAX_RADIUS - MIN_RADIUS) + 0.5f);
    if (r < MIN_RADIUS) r = MIN_RADIUS;
    if (r > MAX_RADIUS) r = MAX_RADIUS;
    return r;
}

static float lum(RGB c) { return 0.2126f * c.r + 0.7152f * c.g + 0.0722f * c.b; }

/* guided bilateral filter; guides that are absent get neutral weight 1.
   Returns 0 on success, -1 if a required guide buffer is missing. */
static int denoise(const Frame* f, RGB* out, GuideReport* rep) {
    rep->missing_normal = !f->has_normal;
    rep->missing_depth = !f->has_depth;
    rep->missing_var = !f->has_var;
    if (rep->missing_normal || rep->missing_depth || rep->missing_var) return -1;

    const float sigma_c = 0.25f;  /* color distance */
    const float sigma_n = 0.15f;  /* normal angle */
    const float sigma_z = 0.05f;  /* relative depth */

    for (int y = 0; y < f->h; y++) {
        for (int x = 0; x < f->w; x++) {
            int i = y * f->w + x;
            int r = radius_for_var(f->var[i]);
            float sigma_s = 0.5f * (float)r;
            RGB cp = f->color[i];
            RGB np = f->normal[i];
            float zp = f->depth[i];
            float lp = lum(cp);
            float wsum = 0.0f, ar = 0.0f, ag = 0.0f, ab = 0.0f;
            for (int dy = -r; dy <= r; dy++) {
                for (int dx = -r; dx <= r; dx++) {
                    int qx = x + dx, qy = y + dy;
                    if (qx < 0 || qy < 0 || qx >= f->w || qy >= f->h) continue;
                    int j = qy * f->w + qx;
                    RGB cq = f->color[j];
                    float dl = lum(cq) - lp;
                    float w_s = expf(-(float)(dx * dx + dy * dy) / (2.0f * sigma_s * sigma_s));
                    float w_c = expf(-(dl * dl) / (2.0f * sigma_c * sigma_c));
                    RGB nq = f->normal[j];
                    float dn = np.r * nq.r + np.g * nq.g + np.b * nq.b;
                    if (dn < 0.0f) dn = 0.0f;
                    float w_n = expf(-(1.0f - dn) / (sigma_n * sigma_n));
                    float dz = (f->depth[j] - zp) / (zp > 1e-6f ? zp : 1e-6f);
                    float w_z = expf(-(dz * dz) / (2.0f * sigma_z * sigma_z));
                    float w = w_s * w_c * w_n * w_z;
                    ar += w * cq.r; ag += w * cq.g; ab += w * cq.b;
                    wsum += w;
                }
            }
            if (wsum > 0.0f) {
                out[i].r = ar / wsum;
                out[i].g = ag / wsum;
                out[i].b = ab / wsum;
            } else {
                out[i] = cp;
            }
        }
    }
    return 0;
}

static double psnr(const RGB* a, const RGB* b, int n) {
    double mse = 0.0;
    for (int i = 0; i < n; i++) {
        double dr = a[i].r - b[i].r, dg = a[i].g - b[i].g, db = a[i].b - b[i].b;
        mse += dr * dr + dg * dg + db * db;
    }
    mse /= (double)(n * 3);
    if (mse < 1e-12) return 99.0;
    return 10.0 * log10(1.0 / mse);
}

/* count reference-edge pixels where the output drifted too far from the
   reference: the filter blurred across the edge. */
static int edge_overblur(const RGB* out, const RGB* ref, int w, int h, float grad_th, float tol) {
    int bad = 0;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            int i = y * w + x;
            float gx = 0.0f, gy = 0.0f;
            if (x + 1 < w) gx = lum(ref[i + 1]) - lum(ref[i]);
            if (y + 1 < h) gy = lum(ref[i + w]) - lum(ref[i]);
            if (fabsf(gx) < grad_th && fabsf(gy) < grad_th) continue;
            float d = fabsf(lum(out[i]) - lum(ref[i]));
            if (d > tol) bad++;
        }
    }
    return bad;
}

/* deterministic rng (splitmix-ish) so tests are reproducible */
static unsigned rng_state = 12345u;
static float frand(void) {
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float)(rng_state >> 8) / 16777216.0f;
}
static float gauss(void) { /* Box-Muller, one sample */
    float u1 = frand(); if (u1 < 1e-7f) u1 = 1e-7f;
    return sqrtf(-2.0f * logf(u1)) * cosf(6.2831853f * frand());
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* radius follows variance */
    check("radius_min_at_zero_var", radius_for_var(0.0f) == MIN_RADIUS);
    check("radius_max_at_ref_var", radius_for_var(VAR_REF) == MAX_RADIUS);
    check("radius_clamps_high_var", radius_for_var(10.0f) == MAX_RADIUS);
    check("radius_grows", radius_for_var(VAR_REF * 0.25f) > radius_for_var(0.0f));

    /* scene: two flat patches split by a vertical edge; the right patch
       faces away (normal z 1 -> 0.2) and is farther (depth 1 -> 2). */
    const int W = 32, H = 32, N = W * H;
    Frame f;
    memset(&f, 0, sizeof(f));
    f.w = W; f.h = H;
    static RGB ref[MAX_W * MAX_H];
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            int i = y * W + x;
            int right = x >= W / 2;
            float v = right ? 0.8f : 0.2f;
            ref[i].r = ref[i].g = ref[i].b = v;
            f.normal[i].r = 0.0f; f.normal[i].g = 0.0f;
            f.normal[i].b = right ? 0.2f : 1.0f;
            f.depth[i] = right ? 2.0f : 1.0f;
            f.var[i] = 0.01f; /* low-spp noise level */
            float n = 0.1f * gauss();
            f.color[i].r = v + n; f.color[i].g = v + n; f.color[i].b = v + n;
        }
    }
    f.has_normal = f.has_depth = f.has_var = 1;

    static RGB out[MAX_W * MAX_H];
    GuideReport rep;
    check("denoise_ok", denoise(&f, out, &rep) == 0);

    double p_noisy = psnr(f.color, ref, N);
    double p_denoised = psnr(out, ref, N);
    printf("psnr_noisy=%.2f\n", p_noisy);
    printf("psnr_denoised=%.2f\n", p_denoised);
    printf("psnr_gain=%.2f\n", p_denoised - p_noisy);
    check("psnr_improves", p_denoised > p_noisy + 3.0);

    int bad = edge_overblur(out, ref, W, H, 0.3f, 0.05f);
    printf("edge_overblur_pixels=%d\n", bad);
    check("edge_overblur_zero", bad == 0);

    /* missing guide buffers must be reported, not silently ignored */
    Frame g = f;
    g.has_normal = 0;
    check("report_missing_normal", denoise(&g, out, &rep) == -1 && rep.missing_normal == 1);
    g = f; g.has_depth = 0;
    check("report_missing_depth", denoise(&g, out, &rep) == -1 && rep.missing_depth == 1);
    g = f; g.has_var = 0;
    check("report_missing_var", denoise(&g, out, &rep) == -1 && rep.missing_var == 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
