/*
 * golden.c -- a known-correct reference implementation of the renderer edge.
 *
 * This is NOT the thing under audit.  It exists so the calibration half of the
 * exercise is meaningful: if the invariant checker cannot pass this, no result
 * it reports about the shipped modules means anything.  It is deliberately
 * written independently of the module code (double precision, exact slab test,
 * correct radical inverse, normalized GGX) rather than by copying it.
 */
#include <math.h>
#include <string.h>

#include "inv_impl.h"

/* ---- geometry ---------------------------------------------------------- */

int golden_box_hit(const double lo[3], const double hi[3],
                   const double o[3], const double d[3],
                   double tmax, double* t) {
    double t0 = 0.0, t1 = tmax;
    for (int i = 0; i < 3; i++) {
        if (d[i] == 0.0) {
            if (o[i] < lo[i] || o[i] > hi[i]) return 0;
        } else {
            double inv = 1.0 / d[i];
            double ta = (lo[i] - o[i]) * inv;
            double tb = (hi[i] - o[i]) * inv;
            if (ta > tb) { double q = ta; ta = tb; tb = q; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    if (t) *t = t0;
    return 1;
}

/* ---- sampling ---------------------------------------------------------- */

static double vdc(unsigned i, unsigned base) {
    double f = 1.0 / (double)base, r = 0.0;
    while (i > 0) { r += f * (double)(i % base); i /= base; f /= (double)base; }
    return r;
}

double golden_halton(unsigned i, unsigned base) { return vdc(i, base); }

/* ---- material ---------------------------------------------------------- */

static double sat(double x) { return x < 0.0 ? 0.0 : x; }

double golden_ggx_D(double nh, double rough) {
    double a = rough * rough;
    double a2 = a * a;
    double nh2 = nh * nh;
    double denom = nh2 * (a2 - 1.0) + 1.0;
    if (denom <= 0.0) return 0.0;
    return a2 / (INV_PI * denom * denom);
}

double golden_lambert(double albedo) { return albedo / INV_PI; }

double golden_cosine_pdf(double cos_theta) { return sat(cos_theta) / INV_PI; }

void golden_cosine_sample(double u1, double u2, double out[3]) {
    /* Malley: concentric-free polar form, valid for u1,u2 in [0,1). */
    double r = sqrt(u1 < 0.0 ? 0.0 : u1);
    double phi = 2.0 * INV_PI * u2;
    double z = sqrt(1.0 - (u1 < 0.0 ? 0.0 : (u1 > 1.0 ? 1.0 : u1)));
    out[0] = r * cos(phi);
    out[1] = r * sin(phi);
    out[2] = z;
}

/* ---- tree build -------------------------------------------------------- */

int golden_build(InvTree* out) {
    /* This reference just mirrors the scene into a flat two-level tree.  The
     * geometry class checks *point-set and containment* properties, not a
     * particular split strategy, so any legal tree is accepted. */
    InvBox scene[INV_MAX_PRIM];
    int n = inv_scene_boxes(scene);
    out->nprims = n;
    for (int i = 0; i < n; i++) out->prim[i] = scene[i];

    InvBox whole;
    for (int a = 0; a < 3; a++) { whole.lo[a] = 1e30; whole.hi[a] = -1e30; }
    for (int i = 0; i < n; i++)
        for (int a = 0; a < 3; a++) {
            if (scene[i].lo[a] < whole.lo[a]) whole.lo[a] = scene[i].lo[a];
            if (scene[i].hi[a] > whole.hi[a]) whole.hi[a] = scene[i].hi[a];
        }
    out->node[0].box = whole;
    out->node[0].left = 1; out->node[0].right = 2;
    out->node[0].first = -1; out->node[0].count = 0;

    int half = n / 2;
    out->node[1].first = 0;    out->node[1].count = half;
    out->node[2].first = half; out->node[2].count = n - half;
    for (int k = 1; k <= 2; k++) {
        out->node[k].left = out->node[k].right = -1;
        for (int a = 0; a < 3; a++) { out->node[k].box.lo[a] = 1e30; out->node[k].box.hi[a] = -1e30; }
        for (int i = out->node[k].first; i < out->node[k].first + out->node[k].count; i++)
            for (int a = 0; a < 3; a++) {
                if (scene[i].lo[a] < out->node[k].box.lo[a]) out->node[k].box.lo[a] = scene[i].lo[a];
                if (scene[i].hi[a] > out->node[k].box.hi[a]) out->node[k].box.hi[a] = scene[i].hi[a];
            }
    }
    out->nnodes = 3;
    return 0;
}

/* ------------------------------------------------------------------------
 * Float control: the same physics done entirely in float, i.e. at the same
 * precision as the modules under audit, but implemented correctly.  Each
 * function is the canonical textbook form; none of it is copied from the
 * modules, so agreement is evidence, not tautology.
 * ---------------------------------------------------------------------- */

static int f_box_hit(const double lo[3], const double hi[3],
                     const double o[3], const double d[3],
                     double tmax, double* t) {
    float t0 = 0.0f, t1 = (float)tmax;
    for (int i = 0; i < 3; i++) {
        float lo_i = (float)lo[i], hi_i = (float)hi[i];
        float o_i = (float)o[i], d_i = (float)d[i];
        /* No special "parallel" epsilon: dividing by a tiny but nonzero d_i is
           exactly what the slab test is defined to do, and the resulting
           interval test is symmetric.  Only a true zero needs the special
           case. */
        if (d_i == 0.0f) {
            if (o_i < lo_i || o_i > hi_i) return 0;
        } else {
            float inv = 1.0f / d_i;
            float ta = (lo_i - o_i) * inv;
            float tb = (hi_i - o_i) * inv;
            if (ta > tb) { float q = ta; ta = tb; tb = q; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    if (t) *t = (double)t0;
    return 1;
}

static double f_halton(unsigned i, unsigned base) {
    float f = 1.0f / (float)base, r = 0.0f;
    while (i > 0) { r += f * (float)(i % base); i /= base; f /= (float)base; }
    return (double)r;
}

static double f_ggx_D(double nh, double rough) {
    float mu = (float)nh, r = (float)rough;
    float a = r * r, a2 = a * a;
    float denom = mu * mu * (a2 - 1.0f) + 1.0f;
    if (denom <= 0.0f) return 0.0;              /* delta limit -> finite 0 */
    return (double)(a2 / (3.14159265f * denom * denom));
}

static double f_lambert(double a) { return (double)((float)a / 3.14159265f); }
static double f_cos_pdf(double c) { float x = (float)c; return (double)((x < 0 ? 0.0f : x) / 3.14159265f); }

static void f_cos_sample(double u1, double u2, double out[3]) {
    float r = sqrtf((float)u1);
    float phi = 2.0f * 3.14159265f * (float)u2;
    float z = sqrtf(1.0f - (float)u1 > 0 ? 1.0f - (float)u1 : 0.0f);
    out[0] = (double)(r * cosf(phi));
    out[1] = (double)(r * sinf(phi));
    out[2] = (double)z;
}

const InvModel inv_float_control = {
    "float-control", f_box_hit, f_halton, f_ggx_D,
    f_lambert, f_cos_pdf, f_cos_sample, golden_build
};

const InvModel inv_golden_reference = {
    "golden-reference", golden_box_hit, golden_halton, golden_ggx_D,
    golden_lambert, golden_cosine_pdf, golden_cosine_sample, golden_build
};

/* The supervisor's reference answer for the geometry class.  It is the same
 * exact slab test; exposed under its own name so a reader can see that the
 * checker's notion of "correct" is not supplied by the code under audit. */
int inv_exact_box_hit(const double lo[3], const double hi[3],
                      const double o[3], const double d[3],
                      double tmax, double* t) {
    return golden_box_hit(lo, hi, o, d, tmax, t);
}
