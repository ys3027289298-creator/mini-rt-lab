
/* light: rectangular area-light sampling with area->solid-angle PDF
   conversion. The Jacobian p_w = p_A * r^2 / |cos_light| supplies the
   cosine factor that pure area sampling misses; backface samples carry
   zero contribution; r -> 0 is reported as singular with a clamped
   fallback. Self-test compares a Monte-Carlo solid-angle estimate with
   the analytic rectangle solid angle and bounds the error. */
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct { float x, y, z; } V3;

static V3 sub3(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

typedef struct {
    V3 c;         /* center; light lies in the plane z = c.z */
    V3 n;         /* unit normal on the emission side, +z here */
    float hx, hy; /* half extents along local x/y */
} RectLight;

typedef enum { LIGHT_OK = 0, LIGHT_BACKFACE, LIGHT_SINGULAR } LightStatus;

typedef struct {
    V3 p;        /* sampled point on the light */
    float pdf_a; /* area pdf = 1/A */
} AreaSample;

typedef struct {
    float pdf_w;     /* solid-angle pdf, 0 on the backface */
    float cos_light; /* cosine at the light, <= 0 on the backface */
    float r2;        /* squared distance, clamped when singular */
    LightStatus status;
} SolidAnglePdf;

#define R2_MIN 1e-12f /* squared-distance floor for the singular fallback */

static float rect_area(RectLight l) { return 4.0f * l.hx * l.hy; }

/* uniform area sample on the rectangle */
static AreaSample rect_sample_area(RectLight l, float u1, float u2) {
    AreaSample s;
    s.p.x = l.c.x + (2.0f * u1 - 1.0f) * l.hx;
    s.p.y = l.c.y + (2.0f * u2 - 1.0f) * l.hy;
    s.p.z = l.c.z;
    s.pdf_a = 1.0f / rect_area(l);
    return s;
}

/* area pdf -> solid-angle pdf: p_w = p_A * r^2 / cos_light.
   Backface (cos_light <= 0) gives pdf 0 and zero contribution;
   r^2 below R2_MIN is reported singular and falls back to the clamp. */
static SolidAnglePdf pdf_area_to_solid_angle(RectLight l, V3 from, AreaSample s) {
    SolidAnglePdf out;
    V3 d = sub3(s.p, from);
    float r2 = dot3(d, d);
    out.status = LIGHT_OK;
    if (r2 < R2_MIN) {
        /* singular: shading point on the light, clamp and report */
        out.status = LIGHT_SINGULAR;
        out.r2 = R2_MIN;
        out.cos_light = 1.0f; /* head-on fallback keeps the pdf finite */
        out.pdf_w = s.pdf_a * R2_MIN;
        return out;
    }
    float r = sqrtf(r2);
    /* cosine at the light between n and the direction back to `from` */
    float cos_light = -dot3(l.n, d) / r;
    out.cos_light = cos_light;
    out.r2 = r2;
    if (cos_light <= 0.0f) {
        out.status = LIGHT_BACKFACE;
        out.pdf_w = 0.0f;
        return out;
    }
    out.pdf_w = s.pdf_a * r2 / cos_light;
    return out;
}

/* one-sample estimate of Le * cos_shade / p_w; exactly zero on the
   light's backface or below the shading horizon */
static float light_contrib(RectLight l, float le, V3 from, V3 ns, AreaSample s) {
    SolidAnglePdf c = pdf_area_to_solid_angle(l, from, s);
    if (c.pdf_w <= 0.0f) return 0.0f;
    V3 d = sub3(s.p, from);
    float r = sqrtf(dot3(d, d));
    float cos_shade = r > 0.0f ? dot3(ns, d) / r : 0.0f;
    if (cos_shade <= 0.0f) return 0.0f;
    return le * cos_shade / c.pdf_w;
}

/* analytic solid angle of the rectangle seen from `from`
   (van Oosterom & Strackee); 0 when `from` is behind the light plane */
static double rect_solid_angle_analytic(RectLight l, V3 from) {
    double d = (double)from.z - l.c.z;
    if (d <= 0.0) return 0.0;
    double x0 = -l.hx - ((double)from.x - l.c.x), x1 = l.hx - ((double)from.x - l.c.x);
    double y0 = -l.hy - ((double)from.y - l.c.y), y1 = l.hy - ((double)from.y - l.c.y);
    double r00 = sqrt(d * d + x0 * x0 + y0 * y0);
    double r10 = sqrt(d * d + x1 * x1 + y0 * y0);
    double r01 = sqrt(d * d + x0 * x0 + y1 * y1);
    double r11 = sqrt(d * d + x1 * x1 + y1 * y1);
    return atan(x0 * y0 / (d * r00)) - atan(x1 * y0 / (d * r10))
         - atan(x0 * y1 / (d * r01)) + atan(x1 * y1 / (d * r11));
}

static unsigned rng_state = 12345u;
static float frand(void) {
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float)(rng_state >> 8) / 16777216.0f;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    RectLight l = { { 0, 0, 0 }, { 0, 0, 1 }, 0.6f, 0.4f };
    V3 ns = { 0, 0, 1 };

    /* Jacobian: p_w == p_A * r^2 / cos_light exactly, and the cosine
       factor genuinely matters (dropping it is the classic bug) */
    V3 from = { 0.1f, -0.2f, 1.0f };
    AreaSample s = rect_sample_area(l, 0.3f, 0.7f);
    SolidAnglePdf c = pdf_area_to_solid_angle(l, from, s);
    float expect = s.pdf_a * c.r2 / c.cos_light;
    check("jacobian_exact", fabsf(c.pdf_w - expect) < 1e-5f * expect);
    check("cosine_factor_matters", fabsf(c.pdf_w - s.pdf_a * c.r2) > 0.05f * c.pdf_w);

    /* backface: pdf and contribution are exactly zero */
    V3 behind = { 0.0f, 0.0f, -1.0f };
    AreaSample sc = rect_sample_area(l, 0.5f, 0.5f);
    SolidAnglePdf cb = pdf_area_to_solid_angle(l, behind, sc);
    check("backface_pdf_zero", cb.status == LIGHT_BACKFACE && cb.pdf_w == 0.0f);
    check("backface_contrib_zero", light_contrib(l, 1.0f, behind, ns, sc) == 0.0f);

    /* singular: r -> 0 is reported and falls back to a finite pdf */
    V3 on_light = { 0.05f, 0.05f, 1e-9f };
    AreaSample sn = rect_sample_area(l, (0.05f / l.hx + 1.0f) * 0.5f,
                                        (0.05f / l.hy + 1.0f) * 0.5f);
    SolidAnglePdf cs = pdf_area_to_solid_angle(l, on_light, sn);
    check("singular_reported", cs.status == LIGHT_SINGULAR);
    check("singular_fallback_finite", isfinite(cs.pdf_w) && cs.pdf_w > 0.0f);
    check("singular_contrib_finite", isfinite(light_contrib(l, 1.0f, on_light, ns, sn)));

    /* analytic formula sanity: centered point, Omega = 4 atan(hx*hy/(d*R)) */
    RectLight lc = { { 0, 0, 0 }, { 0, 0, 1 }, 0.6f, 0.4f };
    V3 center = { 0, 0, 1.0f };
    double R = sqrt(1.0 + 0.36 + 0.16);
    double omega_c = 4.0 * atan(0.24 / R);
    check("analytic_centered", fabs(rect_solid_angle_analytic(lc, center) - omega_c) < 1e-6);

    /* MC solid angle vs analytic: Omega = E[1/p_w] over area samples,
       error bounded by 5 standard errors */
    V3 obs = { 0.3f, 0.1f, 0.8f };
    double omega_ref = rect_solid_angle_analytic(l, obs);
    int n = 200000;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; i++) {
        AreaSample si = rect_sample_area(l, frand(), frand());
        SolidAnglePdf ci = pdf_area_to_solid_angle(l, obs, si);
        double g = ci.pdf_w > 0.0f ? 1.0 / ci.pdf_w : 0.0; /* backface -> 0 */
        sum += g; sum2 += g * g;
    }
    double mean = sum / n;
    double var = (sum2 / n - mean * mean) * n / (n - 1);
    double se = sqrt(var / n);
    double err = fabs(mean - omega_ref);
    double bound = 5.0 * se;
    printf("omega_mc=%.8f omega_ref=%.8f err=%.3g bound=%.3g\n",
           mean, omega_ref, err, bound);
    check("omega_positive", omega_ref > 0.0 && mean > 0.0);
    check("omega_error_bounded", err <= bound);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
