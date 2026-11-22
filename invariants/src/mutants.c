/*
 * mutants.c -- calibration fixtures.
 *
 * Each mutant is the golden reference with exactly one behaviour deliberately
 * broken.  The invariant checker must flag the mutated class and leave the
 * others alone (a miss is a false negative; flagging an untouched class is a
 * false positive).  Mutants reuse the (correct) golden functions except for the
 * single overridden entry, so a failure can only come from that override.
 */
#include <math.h>
#include <string.h>

#include "inv_impl.h"

const int INV_MUTANT_CAP = 12;

/* ---- energy mutants ---------------------------------------------------- */

static double mut_lambert_gain(double a) { return golden_lambert(a) * 1.5; }
static double mut_pdf_double(double c) { return golden_cosine_pdf(c) * 2.0; }
static double mut_ggx_scale(double nh, double r) { return golden_ggx_D(nh, r) * 1.3; }
static double mut_ggx_nan(double nh, double r) {
    if (r == 0.0) return NAN;
    return golden_ggx_D(nh, r);
}

/* ---- geometry mutant: the classic near-parallel-cutoff tunnel ---------- */

static int mut_box_parallel_eps(const double lo[3], const double hi[3],
                                const double o[3], const double d[3],
                                double tmax, double* t) {
    /* Same math as the renderer, but treats |d_i|<1e-8 as parallel and then
       only does an origin-inside rejection.  For a grazing ray whose origin is
       outside the slab this silently returns "miss" even though it hits. */
    double t0 = 0.0, t1 = tmax;
    for (int i = 0; i < 3; i++) {
        if (fabs(d[i]) < 1e-8) {
            if (o[i] < lo[i] || o[i] > hi[i]) return 0;
        } else {
            double inv = 1.0 / d[i];
            double ta = (lo[i] - o[i]) * inv, tb = (hi[i] - o[i]) * inv;
            if (ta > tb) { double q = ta; ta = tb; tb = q; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    if (t) *t = t0;
    return 1;
}

/* ---- sampling mutants -------------------------------------------------- */

static double mut_halton_dup(unsigned i, unsigned b) { return golden_halton(i / 2u, b); }
static double mut_halton_bias(unsigned i, unsigned b) { return golden_halton(i, b) * 0.5; }
static double mut_halton_stride2(unsigned i, unsigned b) { return golden_halton(2u * i, b); }

/* ---- symmetry mutants -------------------------------------------------- */

static void mut_cosine_sample_skew(double u1, double u2, double out[3]) {
    /* φ = π·u2 maps the hemisphere onto a half-disc: the sampler is no longer
       azimuthally symmetric, so the mirror identity and E[x]=E[y]=0 break. */
    double r = sqrt(u1 < 0.0 ? 0.0 : u1);
    double phi = INV_PI * u2;
    double z = sqrt(1.0 - (u1 < 0.0 ? 0.0 : (u1 > 1.0 ? 1.0 : u1)));
    out[0] = r * cos(phi);
    out[1] = r * sin(phi);
    out[2] = z;
}

static double mut_ggx_odd(double nh, double r) {
    double d = golden_ggx_D(nh, r);
    return nh < 0.0 ? d * 1.05 : d;   /* D(h) != D(-h) */
}

/* ------------------------------------------------------------------------ */

/* Bind a mutant: start from a copy of the golden reference, then override the
 * one function named by the caller.  Overriding explicitly (rather than
 * conditionally) keeps the compiler from warning that a function address is
 * always true. */
static InvMutant mk(const char* name, const char* cls) {
    InvMutant m;
    m.name = name;
    m.target_cls = cls;
    m.m = inv_golden_reference;
    m.m.name = name;
    return m;
}

int inv_build_mutants(InvMutant* out, int cap) {
    int k = 0;
    if (cap < 10) return 0;

    out[k] = mk("energy/lambert_1.5x", "energy");
    out[k++].m.lambert_f = mut_lambert_gain;

    out[k] = mk("energy/pdf_2x", "energy");
    out[k++].m.cosine_pdf = mut_pdf_double;

    out[k] = mk("energy/ggx_scale_1.3", "energy");
    out[k++].m.ggx_D = mut_ggx_scale;

    out[k] = mk("energy/ggx_nan_at_rough0", "energy");
    out[k++].m.ggx_D = mut_ggx_nan;

    out[k] = mk("geometry/parallel_eps", "geometry");
    out[k++].m.box_hit = mut_box_parallel_eps;

    out[k] = mk("sampling/duplicate_pairs", "sampling");
    out[k++].m.halton_1d = mut_halton_dup;

    out[k] = mk("sampling/bias_half", "sampling");
    out[k++].m.halton_1d = mut_halton_bias;

    out[k] = mk("sampling/stride2", "sampling");
    out[k++].m.halton_1d = mut_halton_stride2;

    out[k] = mk("symmetry/skew_phi", "symmetry");
    out[k++].m.cosine_sample = mut_cosine_sample_skew;

    out[k] = mk("symmetry/ggx_odd", "symmetry");
    out[k++].m.ggx_D = mut_ggx_odd;

    return k;
}
