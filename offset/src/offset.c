
/* offset: scene-scale aware normal offset for ray origins. Kills self-intersection
   at large world coordinates, bounds hit error against the analytic solution, and
   scans the offset size where nearby geometry starts being missed. */
#include <stdio.h>
#include <math.h>
#include <float.h>

#define IMG_W 48
#define IMG_H 48
#define OFFSET_K 16.0f   /* working offset multiplier, in ulps of scene scale */
#define SELF_TMAX 0.05f  /* secondary hit below this t counts as self-intersection */

typedef struct { float c[3]; float r; } Sphere;

static float dot3f(const float a[3], const float b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/* nearest positive root of ray/sphere, float; -1 on miss */
static float sphere_hit_f(const Sphere* s, const float ro[3], const float rd[3]) {
    float oc[3] = { ro[0]-s->c[0], ro[1]-s->c[1], ro[2]-s->c[2] };
    float b = dot3f(oc, rd);
    float c = dot3f(oc, oc) - s->r*s->r;
    float disc = b*b - c;
    if (disc < 0.0f) return -1.0f;
    float sq = sqrtf(disc);
    float t = -b - sq;
    if (t < 0.0f) t = -b + sq;
    return t >= 0.0f ? t : -1.0f;
}

/* same ray solved in double: analytic ground truth for the float hit */
static double sphere_hit_d(const Sphere* s, const float rof[3], const float rdf[3]) {
    double ro[3] = { rof[0], rof[1], rof[2] };
    double rd[3] = { rdf[0], rdf[1], rdf[2] };
    double oc[3] = { ro[0]-s->c[0], ro[1]-s->c[1], ro[2]-s->c[2] };
    double b = oc[0]*rd[0] + oc[1]*rd[1] + oc[2]*rd[2];
    double c = oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2] - (double)s->r*s->r;
    double disc = b*b - c;
    if (disc < 0.0) return -1.0;
    double sq = sqrt(disc);
    double t = -b - sq;
    if (t < 0.0) t = -b + sq;
    return t >= 0.0 ? t : -1.0;
}

/* scene scale: largest coordinate magnitude any geometry can reach */
static float scene_scale(const Sphere* s, int n) {
    float m = 1.0f;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 3; j++) {
            float a = fabsf(s[i].c[j]) + s[i].r;
            if (a > m) m = a;
        }
    return m;
}

/* epsilon scaled by scene magnitude: k ulps of the largest coordinate in play */
static float offset_eps(const float p[3], float scale, float k) {
    float m = scale;
    for (int i = 0; i < 3; i++) {
        float a = fabsf(p[i]);
        if (a > m) m = a;
    }
    return k * FLT_EPSILON * m;
}

static void offset_origin(float p[3], const float n[3], float eps) {
    for (int i = 0; i < 3; i++) p[i] += n[i] * eps;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* large-coordinate scene: unit-ish geometry placed far from the origin */
    Sphere s = { { 10000.5f, 9999.25f, 0.0f }, 1.5f };
    float scale = scene_scale(&s, 1);
    float ro[3] = { s.c[0], s.c[1], -6.0f };

    int self_plain = 0, self_offset = 0, primary_hits = 0;
    double max_err = 0.0;
    float eps_work = 0.0f;

    for (int y = 0; y < IMG_H; y++) {
        for (int x = 0; x < IMG_W; x++) {
            float u = ((x + 0.5f) / IMG_W - 0.5f);
            float v = ((y + 0.5f) / IMG_H - 0.5f);
            float rd[3] = { u, v, 1.0f };
            float il = 1.0f / sqrtf(dot3f(rd, rd));
            for (int i = 0; i < 3; i++) rd[i] *= il;

            float t = sphere_hit_f(&s, ro, rd);
            if (t < 0.0f) continue;
            primary_hits++;

            double td = sphere_hit_d(&s, ro, rd);
            double err = fabs((double)t - td);
            if (err > max_err) max_err = err;

            float p[3] = { ro[0] + t*rd[0], ro[1] + t*rd[1], ro[2] + t*rd[2] };
            float n[3] = { p[0]-s.c[0], p[1]-s.c[1], p[2]-s.c[2] };
            float in = 1.0f / sqrtf(dot3f(n, n));
            for (int i = 0; i < 3; i++) n[i] *= in;

            /* secondary ray along the outward normal must not re-hit the sphere */
            float t2 = sphere_hit_f(&s, p, n);
            if (t2 >= 0.0f && t2 < SELF_TMAX) self_plain++;

            float po[3] = { p[0], p[1], p[2] };
            float eps = offset_eps(p, scale, OFFSET_K);
            eps_work = eps;
            offset_origin(po, n, eps);
            float t3 = sphere_hit_f(&s, po, n);
            if (t3 >= 0.0f && t3 < SELF_TMAX) self_offset++;
        }
    }

    /* error upper bound: a few dozen ulps of the largest magnitude in play */
    double err_bound = 32.0 * FLT_EPSILON * ((double)scale + 6.0);
    printf("primary_hits=%d\n", primary_hits);
    printf("self_nooffset=%d\n", self_plain);
    printf("self_offset=%d\n", self_offset);
    printf("max_err=%.3e err_bound=%.3e\n", max_err, err_bound);
    printf("eps_working=%.6f\n", eps_work);

    check("primary_hits", primary_hits > IMG_W * IMG_H / 8);
    check("selfhit_present_without_offset", self_plain > 0);
    check("selfhit_zero_with_offset", self_offset == 0);
    check("analytic_err_bounded", max_err < err_bound);
    check("eps_below_feature_size", eps_work < 0.1f * s.r);

    /* missed-hit threshold: detail sphere floating a small gap above the surface;
       once the offset jumps past it, the secondary ray misses it */
    float gap = 0.1f, rdet = 0.05f;
    float n[3] = { 1.0f, 0.0f, 0.0f };
    float ps[3] = { s.c[0] + s.r, s.c[1], s.c[2] };
    Sphere det = { { ps[0] + gap + rdet, ps[1], ps[2] }, rdet };

    float pw[3] = { ps[0], ps[1], ps[2] };
    offset_origin(pw, n, offset_eps(ps, scale, OFFSET_K));
    check("detail_hit_at_working_offset", sphere_hit_f(&det, pw, n) >= 0.0f);

    float eps_thr = -1.0f;
    int k_thr = 0;
    for (int k = 1; k <= 65536; k *= 2) {
        float eps = offset_eps(ps, scale, (float)k);
        float po[3] = { ps[0], ps[1], ps[2] };
        offset_origin(po, n, eps);
        if (sphere_hit_f(&det, po, n) < 0.0f) { eps_thr = eps; k_thr = k; break; }
    }
    printf("miss_threshold_eps=%.6f k=%d\n", eps_thr, k_thr);
    /* true threshold is gap + 2*rdet (origin lands past the detail sphere);
       the doubling scan must report it within one octave */
    check("miss_threshold_found", eps_thr > 0.0f);
    check("miss_threshold_above_gap", eps_thr > gap + rdet);
    check("miss_threshold_within_octave", eps_thr <= 2.0f * (gap + 2.0f * rdet));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
