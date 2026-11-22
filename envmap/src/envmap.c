
/* envmap: equirectangular environment map sampling.
   Bilinear filtering with longitude wraparound (no seam line) and
   latitude clamping at the poles, verified per-pixel against an
   analytic reference with a derived interpolation error bound. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

/* accepted resolution ranges; longitude is the u (width) axis */
#define ENV_MIN_W 2
#define ENV_MAX_W 4096
#define ENV_MIN_H 2
#define ENV_MAX_H 4096

typedef struct { float r, g, b; } RGB;
typedef struct { float x, y, z; } V3;

typedef struct {
    RGB* tex;
    int w, h; /* w: longitude resolution, h: latitude resolution */
} Env;

/* reject bad resolutions and say which values are accepted */
static int env_check_res(int w, int h) {
    if (w < ENV_MIN_W || w > ENV_MAX_W) {
        fprintf(stderr, "envmap: longitude resolution %d invalid, valid range [%d, %d]\n",
                w, ENV_MIN_W, ENV_MAX_W);
        return 0;
    }
    if (h < ENV_MIN_H || h > ENV_MAX_H) {
        fprintf(stderr, "envmap: latitude resolution %d invalid, valid range [%d, %d]\n",
                h, ENV_MIN_H, ENV_MAX_H);
        return 0;
    }
    return 1;
}

static int env_init(Env* e, int w, int h) {
    if (!env_check_res(w, h)) return -1;
    e->tex = (RGB*)malloc(sizeof(RGB) * (size_t)w * (size_t)h);
    if (!e->tex) return -1;
    e->w = w;
    e->h = h;
    return 0;
}

static void env_free(Env* e) { free(e->tex); e->tex = 0; }

static RGB texel(const Env* e, int i, int j) { return e->tex[j * e->w + i]; }

static float lerp1(float a, float b, double t) { return (float)(a + t * (b - a)); }

/* bilinear sample at (u, v); longitude wraps, latitude clamps.
   u and u + k (integer k) select identical texels and weights, so the
   sampled color is exactly periodic across the seam. */
static RGB env_sample_uv(const Env* e, double u, double v) {
    u -= floor(u); /* wrap longitude into [0, 1) */
    if (v < 0.0) v = 0.0;
    if (v > 1.0) v = 1.0;

    double x = u * e->w - 0.5; /* texel-center coordinates */
    double y = v * e->h - 0.5;

    int i0 = (int)floor(x);
    double fx = x - i0;
    i0 %= e->w; if (i0 < 0) i0 += e->w;        /* column -1 is column w-1 */
    int i1 = i0 + 1; if (i1 >= e->w) i1 -= e->w; /* column w is column 0 */

    int j0 = (int)floor(y);
    double fy = y - j0;
    if (j0 < 0) { j0 = 0; fy = 0.0; }                /* north pole: clamp */
    else if (j0 >= e->h - 1) { j0 = e->h - 2; fy = 1.0; } /* south pole */
    int j1 = j0 + 1;

    RGB c00 = texel(e, i0, j0), c10 = texel(e, i1, j0);
    RGB c01 = texel(e, i0, j1), c11 = texel(e, i1, j1);
    RGB c;
    c.r = lerp1(lerp1(c00.r, c10.r, fx), lerp1(c01.r, c11.r, fx), fy);
    c.g = lerp1(lerp1(c00.g, c10.g, fx), lerp1(c01.g, c11.g, fx), fy);
    c.b = lerp1(lerp1(c00.b, c10.b, fx), lerp1(c01.b, c11.b, fx), fy);
    return c;
}

static void dir_to_uv(V3 d, double* u, double* v) {
    double l = sqrt((double)d.x * d.x + (double)d.y * d.y + (double)d.z * d.z);
    double y = l > 0.0 ? d.y / l : 0.0;
    if (y > 1.0) y = 1.0;
    if (y < -1.0) y = -1.0;
    double phi = l > 0.0 ? atan2((double)d.z, (double)d.x) : 0.0;
    *u = phi / (2.0 * PI) + 0.5;
    *v = acos(y) / PI;
}

static RGB env_sample(const Env* e, V3 d) {
    double u, v;
    dir_to_uv(d, &u, &v);
    return env_sample_uv(e, u, v);
}

/* analytic reference: smooth and 1-periodic in u, so the error bound
   below also holds across the wrap interval.
   every channel has |f_uu| <= 2*pi^2 and |f_vv| <= pi^2/2. */
static RGB analytic(double u, double v) {
    double s = sin(PI * v);
    RGB c = {
        (float)(0.5 + 0.5 * cos(2.0 * PI * u) * s),
        (float)(0.5 + 0.5 * sin(2.0 * PI * u) * s),
        (float)(0.5 + 0.5 * cos(PI * v))
    };
    return c;
}

/* error bound for sampling the analytic reference from a w x h grid
   with spacing (hu, hv):
   - interior: bilinear interpolation error
       |f - Bf| <= (hu^2 * max|f_uu| + hv^2 * max|f_vv|) / 8
     (1D linear interpolation error (1/8)h^2 M2 applied per axis;
     interpolation is a convex combination, so the per-axis bounds add)
   - pole rows: clamp-to-edge over the boundary half-interval adds a
     first-order term (hv/2) * max|f_v|
   longitude wraps, so there is no u-boundary term. */
static double analytic_bound(int w, int h) {
    double hu = 1.0 / w, hv = 1.0 / h;
    double m_uu = 2.0 * PI * PI, m_vv = 0.5 * PI * PI, m_v = 0.5 * PI;
    return (hu * hu * m_uu + hv * hv * m_vv) / 8.0 + 0.5 * hv * m_v;
}

static float rgb_maxdiff(RGB a, RGB b) {
    float d = fabsf(a.r - b.r);
    if (fabsf(a.g - b.g) > d) d = fabsf(a.g - b.g);
    if (fabsf(a.b - b.b) > d) d = fabsf(a.b - b.b);
    return d;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Env e;

    /* illegal longitude/latitude resolutions are rejected and the
       valid range is reported on stderr */
    check("res_w_zero_rejected", env_init(&e, 0, 32) != 0);
    check("res_w_one_rejected", env_init(&e, 1, 32) != 0);
    check("res_w_neg_rejected", env_init(&e, -16, 32) != 0);
    check("res_w_huge_rejected", env_init(&e, ENV_MAX_W + 1, 32) != 0);
    check("res_h_one_rejected", env_init(&e, 32, 1) != 0);
    check("res_min_ok", env_init(&e, ENV_MIN_W, ENV_MIN_H) == 0);
    env_free(&e);

    const int W = 64, H = 32;
    check("res_ok", env_init(&e, W, H) == 0);
    for (int j = 0; j < H; j++)
        for (int i = 0; i < W; i++)
            e.tex[j * W + i] = analytic((i + 0.5) / W, (j + 0.5) / H);

    /* seam: sampling must be exactly periodic in longitude —
       f(u), f(u+1), f(u-1) are bit-identical */
    int seam_ok = 1;
    for (int k = 0; k <= 256 && seam_ok; k++) {
        double u = k / 256.0;
        for (int j = 0; j < H; j++) {
            double v = (j + 0.5) / H;
            RGB a = env_sample_uv(&e, u, v);
            RGB b = env_sample_uv(&e, u + 1.0, v);
            RGB c = env_sample_uv(&e, u - 1.0, v);
            if (memcmp(&a, &b, sizeof(RGB)) != 0 || memcmp(&a, &c, sizeof(RGB)) != 0) {
                seam_ok = 0;
                break;
            }
        }
    }
    check("seam_periodic_exact", seam_ok);

    /* the two atan2 branch-cut sides of the seam direction (+0/-0 in z)
       land on u = 1.0 and u = 0.0 and must give the identical texel */
    V3 dpos = { -1.0f, 0.25f, 0.0f };
    V3 dneg = { -1.0f, 0.25f, -0.0f };
    RGB sp = env_sample(&e, dpos);
    RGB sn = env_sample(&e, dneg);
    check("seam_dir_exact", memcmp(&sp, &sn, sizeof(RGB)) == 0);

    /* across the seam the filter blends column w-1 with column 0 */
    {
        int j = 5;
        RGB got = env_sample_uv(&e, -0.25 / W, (j + 0.5) / H);
        RGB c0 = texel(&e, W - 1, j), c1 = texel(&e, 0, j);
        RGB want = { lerp1(c0.r, c1.r, 0.25), lerp1(c0.g, c1.g, 0.25), lerp1(c0.b, c1.b, 0.25) };
        check("seam_wrap_blend", memcmp(&got, &want, sizeof(RGB)) == 0);
    }

    /* per-pixel comparison against the analytic reference, through the
       full direction -> uv -> bilinear path */
    double maxerr = 0.0;
    const int PW = 200, PH = 100;
    for (int py = 0; py < PH; py++) {
        double v = (py + 0.5) / PH;
        double theta = v * PI;
        for (int px = 0; px < PW; px++) {
            double u = (px + 0.5) / PW;
            double phi = (u - 0.5) * 2.0 * PI;
            V3 d = { (float)(sin(theta) * cos(phi)), (float)cos(theta),
                     (float)(sin(theta) * sin(phi)) };
            RGB s = env_sample(&e, d);
            double err = rgb_maxdiff(s, analytic(u, v));
            if (err > maxerr) maxerr = err;
        }
    }
    double bound = analytic_bound(W, H);
    printf("analytic_max_err=%.6g err_bound=%.6g\n", maxerr, bound);
    check("analytic_within_bound", maxerr <= bound);

    /* poles: clamped latitude keeps the sample near-constant while
       longitude sweeps around the pole */
    double spread = 0.0;
    for (int pole = 0; pole < 2; pole++) {
        double v = pole == 0 ? 0.0 : 1.0;
        RGB mn = { 1e30f, 1e30f, 1e30f }, mx = { -1e30f, -1e30f, -1e30f };
        for (int k = 0; k < 64; k++) {
            RGB c = env_sample_uv(&e, k / 64.0, v);
            if (c.r < mn.r) mn.r = c.r;
            if (c.r > mx.r) mx.r = c.r;
            if (c.g < mn.g) mn.g = c.g;
            if (c.g > mx.g) mx.g = c.g;
            if (c.b < mn.b) mn.b = c.b;
            if (c.b > mx.b) mx.b = c.b;
        }
        float s = rgb_maxdiff(mx, mn);
        if (s > spread) spread = s;
    }
    /* pole row sits at v = 0.5/H, so spread <= sin(pi/(2H)) <= pi/(2H) */
    double pole_tol = PI / (2.0 * H);
    printf("pole_spread=%.6g pole_tol=%.6g\n", spread, pole_tol);
    check("pole_constant", spread <= pole_tol + 1e-6);

    V3 npole = { 0.0f, 1.0f, 0.0f }, spole = { 0.0f, -1.0f, 0.0f };
    RGB cn = env_sample(&e, npole), cs = env_sample(&e, spole);
    check("pole_finite", isfinite(cn.r) && isfinite(cn.g) && isfinite(cn.b) &&
                         isfinite(cs.r) && isfinite(cs.g) && isfinite(cs.b));

    env_free(&e);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
