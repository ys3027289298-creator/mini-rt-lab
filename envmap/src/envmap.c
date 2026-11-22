
/* envmap: equirectangular environment-map direction lookup.
   Builds a per-texel 3D direction LUT under a yaw rotation about +Y,
   with a pole/seam-safe inverse query and analytic error bounds.

   Parameterization (y-up):
     u = (i + 0.5) / W,  v = (j + 0.5) / H        texel centers, never
     phi = 2*pi*u,        theta = pi*v               land on a pole
     d0 = (sin(theta)cos(phi), cos(theta), sin(theta)sin(phi))
   Stored directions are d0 rotated by yaw around +Y; the inverse query
   un-rotates first, so rotating the map and rotating queries stay
   consistent. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define EM_MAX_W 256
#define EM_MAX_H 128
#define PI_F 3.14159265358979323846f
#define PI_D 3.14159265358979323846

typedef struct { float x, y, z; } V3;

typedef struct {
    int w, h;
    float yaw, cos_yaw, sin_yaw;
    V3 lut[EM_MAX_H][EM_MAX_W];
} EnvMap;

static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0.0f) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}

static float angle3(V3 a, V3 b) {
    /* atan2(|cross|, dot) in double: acos(dot) loses all precision when
       the directions nearly coincide (dot rounds to exactly 1 in float). */
    double ax = a.x, ay = a.y, az = a.z;
    double bx = b.x, by = b.y, bz = b.z;
    double cx = ay * bz - az * by;
    double cy = az * bx - ax * bz;
    double cz = ax * by - ay * bx;
    double s = sqrt(cx * cx + cy * cy + cz * cz);
    double d = ax * bx + ay * by + az * bz;
    return (float)atan2(s, d);
}

/* R_y(yaw): x' = c*x + s*z, y' = y, z' = -s*x + c*z */
static V3 rotate_y(V3 d, float c, float s) {
    V3 r = { c * d.x + s * d.z, d.y, -s * d.x + c * d.z };
    return r;
}

/* Build the direction LUT. Returns 0 on success, -1 on illegal
   resolution: each axis needs >= 2 texels (one row cannot reach both
   poles, one column cannot represent the azimuthal seam) and must fit
   the fixed-capacity table. */
static int env_build(EnvMap* em, int w, int h, float yaw) {
    if (w < 2 || h < 2 || w > EM_MAX_W || h > EM_MAX_H) return -1;
    em->w = w;
    em->h = h;
    em->yaw = yaw;
    em->cos_yaw = cosf(yaw);
    em->sin_yaw = sinf(yaw);
    for (int j = 0; j < h; j++) {
        float v = ((float)j + 0.5f) / (float)h;
        float theta = PI_F * v;
        float st = sinf(theta), ct = cosf(theta);
        for (int i = 0; i < w; i++) {
            float u = ((float)i + 0.5f) / (float)w;
            float phi = 2.0f * PI_F * u;
            V3 d = { st * cosf(phi), ct, st * sinf(phi) };
            d = rotate_y(d, em->cos_yaw, em->sin_yaw);
            /* renormalize: kills float drift, keeps pole-row texels unit */
            em->lut[j][i] = norm3(d);
        }
    }
    return 0;
}

/* Texel direction. Column wraps across the seam; row clamps at the
   poles so out-of-range neighbor lookups never read garbage. */
static V3 env_dir(const EnvMap* em, int i, int j) {
    if (j < 0) j = 0;
    if (j >= em->h) j = em->h - 1;
    i %= em->w;
    if (i < 0) i += em->w;
    return em->lut[j][i];
}

/* Inverse query: unit direction -> texel index.
   Returns -1 on a zero/near-zero input; otherwise always yields a
   valid texel. Exact poles take a dedicated branch because atan2(0,0)
   is undefined there. */
static int env_pixel(const EnvMap* em, V3 d, int* pi, int* pj) {
    float n = sqrtf(dot3(d, d));
    if (!(n > 1e-20f)) return -1;
    float x = d.x / n, y = d.y / n, z = d.z / n;

    /* undo yaw: R_y^T */
    float c = em->cos_yaw, s = em->sin_yaw;
    float x0 = c * x - s * z;
    float y0 = y;
    float z0 = s * x + c * z;

    float horiz = sqrtf(fmaxf(0.0f, 1.0f - y0 * y0));
    if (horiz < 1e-7f) {
        /* exact pole: longitude is undefined, pin column 0 */
        *pi = 0;
        *pj = (y0 > 0.0f) ? 0 : em->h - 1;
        return 0;
    }

    float theta = acosf(fmaxf(-1.0f, fminf(1.0f, y0)));
    float phi = atan2f(z0, x0);
    if (phi < 0.0f) phi += 2.0f * PI_F;          /* seam: [-pi,pi] -> [0,2pi) */

    int i = (int)floorf(phi / (2.0f * PI_F) * em->w);
    int j = (int)floorf(theta / PI_F * em->h);
    if (i >= em->w) i = 0;                        /* rounding at the seam */
    if (i < 0) i = 0;
    if (j >= em->h) j = em->h - 1;
    if (j < 0) j = 0;
    *pi = i;
    *pj = j;
    return 0;
}

/* Double-precision analytic reference, same convention as env_build. */
static void analytic_dir(int w, int h, int i, int j, double yaw, double out[3]) {
    double u = ((double)i + 0.5) / (double)w;
    double v = ((double)j + 0.5) / (double)h;
    double phi = 2.0 * PI_D * u;
    double theta = PI_D * v;
    double x0 = sin(theta) * cos(phi);
    double y0 = cos(theta);
    double z0 = sin(theta) * sin(phi);
    double c = cos(yaw), s = sin(yaw);
    out[0] = c * x0 + s * z0;
    out[1] = y0;
    out[2] = -s * x0 + c * z0;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- illegal resolution must be reported, not silently built ---- */
    int dims[][2] = {
        { 0, 8 }, { 8, 0 }, { -1, 8 }, { 8, -1 },
        { 1, 8 }, { 8, 1 },
        { EM_MAX_W + 1, 8 }, { 8, EM_MAX_H + 1 }
    };
    int rejected = 1;
    for (unsigned k = 0; k < sizeof(dims) / sizeof(dims[0]); k++) {
        EnvMap tmp;
        if (env_build(&tmp, dims[k][0], dims[k][1], 0.0f) != -1) rejected = 0;
    }
    check("res_reject_illegal", rejected);

    const int W = 64, H = 32;
    const float YAW = 0.7f;
    EnvMap em;
    memset(&em, 0, sizeof(em));
    check("build_valid", env_build(&em, W, H, YAW) == 0);

    EnvMap em0;
    memset(&em0, 0, sizeof(em0));
    check("build_zero_yaw", env_build(&em0, W, H, 0.0f) == 0);

    /* ---- LUT vs analytic reference: whole-sphere forward error ---- */
    int nondeg = 1;
    double max_fwd = 0.0;
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            V3 d = em.lut[j][i];
            if (!isfinite(d.x) || !isfinite(d.y) || !isfinite(d.z)) nondeg = 0;
            if (fabsf(sqrtf(dot3(d, d)) - 1.0f) > 1e-6f) nondeg = 0;
            double a[3];
            analytic_dir(W, H, i, j, (double)YAW, a);
            V3 ref = { (float)a[0], (float)a[1], (float)a[2] };
            double ang = (double)angle3(d, ref);
            if (ang > max_fwd) max_fwd = ang;
        }
    }
    check("lut_unit_nondeg", nondeg);
    check("lut_vs_analytic", max_fwd < 1e-6);

    /* ---- pole rows: texel centers keep a nonzero horizontal part ---- */
    int pole_bad = 0;
    int pole_rows[2] = { 0, H - 1 };
    for (unsigned r = 0; r < 2; r++) {
        int j = pole_rows[r];
        float expect_horiz = sinf(PI_F * ((float)j + 0.5f) / (float)H);
        for (int i = 0; i < W; i++) {
            V3 d = em.lut[j][i];
            float horiz = sqrtf(d.x * d.x + d.z * d.z);
            if (!(horiz > 0.0f)) pole_bad++;
            if (fabsf(horiz - expect_horiz) > 1e-6f) pole_bad++;
            if (j == 0 && d.y <= 0.0f) pole_bad++;
            if (j == H - 1 && d.y >= 0.0f) pole_bad++;
        }
    }
    check("pole_rows_nondeg", pole_bad == 0);

    /* ---- exact-pole queries are safe; zero vector is rejected ---- */
    V3 north = { 0, 1, 0 }, south = { 0, -1, 0 }, zero = { 0, 0, 0 };
    int pi, pj;
    int pole_q = env_pixel(&em, north, &pi, &pj) == 0 && pj == 0 &&
                 pi >= 0 && pi < W;
    pole_q &= env_pixel(&em, south, &pi, &pj) == 0 && pj == H - 1;
    pole_q &= env_pixel(&em, zero, &pi, &pj) == -1;
    check("pole_query_safe", pole_q);

    /* a yaw leaves the pole direction at the pole row */
    V3 rn = rotate_y(north, em.cos_yaw, em.sin_yaw);
    check("pole_yaw_invariant", env_pixel(&em, rn, &pi, &pj) == 0 && pj == 0);

    /* ---- seam: wrapped neighbor access and angular continuity ---- */
    V3 lm = env_dir(&em, -1, 5), l0 = env_dir(&em, W - 1, 5);
    V3 rm = env_dir(&em, W, 5), r0 = env_dir(&em, 0, 5);
    check("seam_wrap_access",
          lm.x == l0.x && lm.y == l0.y && lm.z == l0.z &&
          rm.x == r0.x && rm.y == r0.y && rm.z == r0.z);

    float seam_ang = angle3(em.lut[5][W - 1], em.lut[5][0]);
    float az_texel = 2.0f * PI_F / (float)W;
    check("seam_neighbor_angle", seam_ang <= az_texel + 1e-6f);

    /* inverse query on either side of the seam lands in the right column */
    V3 d_phi0 = { 1, 0, 0 };
    float phi_end = 2.0f * PI_F - 0.01f;
    V3 d_phie = { cosf(phi_end), 0, sinf(phi_end) };
    int i0, j0, ie, je;
    env_pixel(&em0, d_phi0, &i0, &j0);
    env_pixel(&em0, d_phie, &ie, &je);
    check("seam_inverse", i0 == 0 && ie == W - 1 && j0 == je);

    /* ---- texel centers round-trip exactly, yawed map included ---- */
    int rt_bad = 0;
    for (int j = 0; j < H; j++) {
        for (int i = 0; i < W; i++) {
            int qi, qj;
            if (env_pixel(&em, em.lut[j][i], &qi, &qj) != 0 ||
                qi != i || qj != j) rt_bad++;
        }
    }
    check("center_roundtrip", rt_bad == 0);

    /* ---- analytic upper bound on inverse-query angular error.
       A texel cell spans dtheta = pi/H and dphi = 2*pi/W; on the unit
       sphere ds^2 = dtheta^2 + sin(theta)^2 dphi^2 <= dtheta^2 + dphi^2,
       so any point in the cell is within half the cell diagonal of
       its center: B = sqrt((pi/2H)^2 + (pi/W)^2). ---- */
    float half_dtheta = PI_F / (2.0f * (float)H);
    float half_dphi = PI_F / (float)W;
    float cell_bound = sqrtf(half_dtheta * half_dtheta +
                             half_dphi * half_dphi);
    unsigned seed = 0x12345678u;
    float worst = 0.0f;
    int over = 0;
    for (int k = 0; k < 3000; k++) {
        seed = seed * 1664525u + 1013904223u;
        float r1 = (float)(seed >> 8) / 16777216.0f;
        seed = seed * 1664525u + 1013904223u;
        float r2 = (float)(seed >> 8) / 16777216.0f;
        float y = 1.0f - 2.0f * r1;
        float ph = 2.0f * PI_F * r2;
        float rh = sqrtf(fmaxf(0.0f, 1.0f - y * y));
        V3 d = { rh * cosf(ph), y, rh * sinf(ph) };
        int qi, qj;
        if (env_pixel(&em, d, &qi, &qj) != 0) { over++; continue; }
        float ang = angle3(d, env_dir(&em, qi, qj));
        if (ang > worst) worst = ang;
        if (ang > cell_bound + 1e-6f) over++;
    }
    check("cell_error_bound", over == 0);

    /* ---- rotation consistency: rotating the map equals rotating the
       query the other way (the regression behind sampled-position
       misalignment) ---- */
    int rot_bad = 0;
    for (int k = 0; k < 500; k++) {
        seed = seed * 1664525u + 1013904223u;
        float r1 = (float)(seed >> 8) / 16777216.0f;
        seed = seed * 1664525u + 1013904223u;
        float r2 = (float)(seed >> 8) / 16777216.0f;
        float y = 1.0f - 2.0f * r1;
        float ph = 2.0f * PI_F * r2;
        float rh = sqrtf(fmaxf(0.0f, 1.0f - y * y));
        V3 d = { rh * cosf(ph), y, rh * sinf(ph) };
        int i_a, j_a, i_b, j_b;
        env_pixel(&em0, d, &i_a, &j_a);
        V3 dr = rotate_y(d, em.cos_yaw, em.sin_yaw);
        env_pixel(&em, dr, &i_b, &j_b);
        if (i_a != i_b || j_a != j_b) rot_bad++;
    }
    check("rotation_consistency", rot_bad == 0);

    printf("max_forward_err=%.3e (bound 1e-6)\n", max_fwd);
    printf("cell_bound=%.6f max_cell_err=%.6f\n", (double)cell_bound, (double)worst);
    printf("seam_angle=%.6f azimuth_texel=%.6f\n", (double)seam_ang, (double)az_texel);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
