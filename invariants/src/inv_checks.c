/*
 * inv_checks.c -- the four invariant classes, written against InvModel.
 *
 * Conventions used everywhere below:
 *   - "measured" is always an error/distance where 0 is perfect, so a check
 *     passes iff measured <= threshold.  Recall is reported as 1-recall.
 *   - each check names its class; a model "fails a class" if any check in it
 *     fails.  Calibration asserts that a broken model fails its target class
 *     and that the known-correct reference fails none.
 *   - the checker never calls the renderer's own check()/main(); it only calls
 *     pure functions through InvModel.
 */
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "invariants.h"

/* ------------------------------------------------------------------ report */

static void addf(InvReport* rep, const char* cls, const char* name, int pass,
                 double measured, double threshold, const char* fmt, ...) {
    if (rep->n >= (int)(sizeof rep->r / sizeof rep->r[0])) return;
    InvResult* r = &rep->r[rep->n++];
    r->cls = cls; r->name = name; r->pass = pass;
    r->measured = measured; r->threshold = threshold;
    va_list ap; va_start(ap, fmt);
    vsnprintf(r->detail, sizeof r->detail, fmt, ap);
    va_end(ap);
}

int inv_class_total(const InvReport* rep, const char* cls) {
    int n = 0;
    for (int i = 0; i < rep->n; i++) if (strcmp(rep->r[i].cls, cls) == 0) n++;
    return n;
}

int inv_class_failed(const InvReport* rep, const char* cls) {
    for (int i = 0; i < rep->n; i++)
        if (strcmp(rep->r[i].cls, cls) == 0 && !rep->r[i].pass) return 1;
    return 0;
}

void inv_run_model(const InvModel* m, InvReport* rep) {
    inv_energy(m, rep);
    inv_geometry(m, rep);
    inv_sampling(m, rep);
    inv_symmetry(m, rep);
}

/* ------------------------------------------------------------------ scene */

/* Deterministic scene shared by every tree builder, so "same scene" is not a
 * matter of trust.  LCG constants are the classic Numerical Recipes ones. */
int inv_scene_boxes(InvBox out[INV_MAX_PRIM]) {
    unsigned long s = 0x2545F491u;
    int n = 24;
    for (int i = 0; i < n; i++) {
        s = s * 1664525u + 1013904223u; double a = (s >> 8) / 16777216.0;
        s = s * 1664525u + 1013904223u; double b = (s >> 8) / 16777216.0;
        s = s * 1664525u + 1013904223u; double c = (s >> 8) / 16777216.0;
        s = s * 1664525u + 1013904223u; double d = (s >> 8) / 16777216.0;
        double x = a * 16.0 - 8.0, y = b * 16.0 - 8.0, z = c * 16.0 - 8.0;
        double sx = 0.05 + d * 1.5, sy = 0.05 + a * 1.2, sz = 0.05 + b * 1.0;
        out[i].lo[0] = x;       out[i].lo[1] = y;       out[i].lo[2] = z;
        out[i].hi[0] = x + sx;  out[i].hi[1] = y + sy;  out[i].hi[2] = z + sz;
    }
    return n;
}

/* ================================================================== ENERGY */

/* Reflectance of a Lambert BRDF over the hemisphere is exactly its albedo:
 *   R = integral f cos dω = (albedo/π) · π = albedo.
 * Energy conservation requires R <= 1 for albedo <= 1.  We also build an MC
 * estimator from the model's *own* sampler and pdf; that estimator must equal
 * the analytic value, which catches a sampler/pdf mismatch as well. */
static void energy_lambert(const InvModel* m, InvReport* rep) {
    const double albedo[] = { 0.0, 0.2, 0.5, 1.0 };
    double worst_le1 = 0.0, worst_consist = 0.0, worst_neg = 0.0;
    int K = 64;

    for (int a = 0; a < 4; a++) {
        double f = m->lambert_f(albedo[a]);
        double R = f * INV_PI;                 /* analytic reflectance */
        if (R - 1.0 > worst_le1) worst_le1 = R - 1.0;
        if (-f > worst_neg) worst_neg = -f;

        double acc = 0.0; int n = 0;
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                double u1 = (i + 0.5) / K, u2 = (j + 0.5) / K, s[3];
                m->cosine_sample(u1, u2, s);
                double mu = s[2];
                double pdf = m->cosine_pdf(mu);
                /* f*cos/pdf, with cos = mu; skip the measure-zero μ==0 plane */
                if (pdf > 0.0) { acc += m->lambert_f(albedo[a]) * mu / pdf; n++; }
            }
        }
        if (n > 0) {
            double est = acc / n;
            double err = fabs(est - R);
            if (err > worst_consist) worst_consist = err;
        }
    }

    addf(rep, "energy", "lambert_reflectance_le_one",
         worst_le1 <= TH_ENERGY_INCREASE, worst_le1, TH_ENERGY_INCREASE,
         "max(R-1)=%.3e over albedo<=1", worst_le1);
    addf(rep, "energy", "lambert_mc_consistency",
         worst_consist <= TH_ENERGY_CONSIST, worst_consist, TH_ENERGY_CONSIST,
         "|MC reflectance - analytic| = %.3e", worst_consist);
    addf(rep, "energy", "brdf_nonnegative",
         worst_neg <= TH_ENERGY_NEG, worst_neg, TH_ENERGY_NEG,
         "max(-f)=%.3e", worst_neg);
}

/* ∫_hemisphere (cos/π) dω = 1, evaluated as 2π ∫_0^1 pdf(μ) dμ. */
static void energy_pdf_norm(const InvModel* m, InvReport* rep) {
    int n = 200000;
    double h = 1.0 / n, s = 0.0;
    int bad = 0;
    for (int i = 0; i <= n; i++) {
        double mu = i * h;
        double p = m->cosine_pdf(mu);
        if (!isfinite(p) || p < -TH_ENERGY_NEG) bad++;
        double w = (i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
        s += w * p;
    }
    double integral = 2.0 * INV_PI * s * h / 3.0;
    double err = fabs(integral - 1.0);
    if (!isfinite(err)) err = 1e9;

    addf(rep, "energy", "cosine_pdf_normalized",
         err <= TH_PDF_NORM && bad == 0, err, TH_PDF_NORM,
         "|2pi*int pdf dmu - 1|=%.3e integral=%.9f neg/illformed=%d",
         err, integral, bad);
}

/* GGX NDF normalization: 2π ∫_0^1 D(μ) μ dμ = 1 (the μ factor is n·h). */
static void energy_ndf_norm(const InvModel* m, InvReport* rep) {
    const double rough[] = { 0.05, 0.1, 0.2, 0.4, 0.7, 1.0 };
    double worst = 0.0, worst_rough = 0.0;
    int n = TH_NDF_SIMPSON_N;
    double h = 1.0 / n;

    for (int r = 0; r < 6; r++) {
        double s = 0.0; int na = 0;
        for (int i = 0; i <= n; i++) {
            double mu = i * h;
            double D = m->ggx_D(mu, rough[r]);
            if (!isfinite(D)) { na = 1; break; }
            double w = (i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0);
            s += w * D * mu;
        }
        double integral = na ? NAN : 2.0 * INV_PI * s * h / 3.0;
        double err = na ? 1e9 : fabs(integral - 1.0);
        if (err > worst) { worst = err; worst_rough = rough[r]; }
    }
    addf(rep, "energy", "ggx_ndf_normalized",
         worst <= TH_NDF_NORM, worst, TH_NDF_NORM,
         "max|integral-1|=%.3e (worst rough=%.2f)", worst, worst_rough);
}

/* BRDF/NDF values must be finite and non-negative on the whole legal domain.
 * A renderer multiplies these by finite weights and sums, so a single NaN or
 * Inf poisons the accumulated radiance -- an energy invariant cannot even be
 * stated with non-finite terms.  The rough->0 delta is represented by the
 * finite convention D=0 at the singular point (as a correct implementation
 * does); we require that convention here rather than allowing Inf. */
static void energy_finite(const InvModel* m, InvReport* rep) {
    const double rough[] = { 0.0, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0 };
    const double nh[] = { 0.0, 0.25, 0.5, 0.75, 0.99, 1.0 };
    int nan_bad = 0, inf_bad = 0, neg_bad = 0;
    double fr = 0.0, fn = 0.0; int have_first = 0;

    for (int r = 0; r < 7; r++)
        for (int i = 0; i < 6; i++) {
            double D = m->ggx_D(nh[i], rough[r]);
            int bad = 0;
            if (isnan(D)) { nan_bad++; bad = 1; }
            else if (!isfinite(D)) { inf_bad++; bad = 1; }
            else if (D < -TH_ENERGY_NEG) { neg_bad++; bad = 1; }
            if (bad && !have_first) { fr = rough[r]; fn = nh[i]; have_first = 1; }
        }
    int bad = nan_bad + inf_bad + neg_bad;
    addf(rep, "energy", "ggx_finite_nonnegative",
         bad == 0, (double)bad, 0.0,
         "nan=%d inf=%d neg=%d first at rough=%.3g nh=%.3g",
         nan_bad, inf_bad, neg_bad, fr, fn);
}

void inv_energy(const InvModel* m, InvReport* rep) {
    energy_lambert(m, rep);
    energy_pdf_norm(m, rep);
    energy_ndf_norm(m, rep);
    energy_finite(m, rep);
}

/* ================================================================ GEOMETRY */

/* ---- float rounding helpers: the renderer receives float inputs --------- */

static double fsnap(double x) { return (double)(float)x; }
static double fnext(double x, double dir) {
    return (double)nextafterf((float)x, dir > 0 ? INFINITY : -INFINITY);
}

/*
 * A ray/box case is "ambiguous" if the exact hit/miss answer is not stable
 * under a one-float-ulp change of any input component.  This is the decision
 * boundary for geometry: outside it, a correctly implemented float slab test
 * has a well-defined answer and must match exact arithmetic.  Inside it, no
 * implementation can be blamed.
 */
static int geom_ambiguous(const double lo[3], const double hi[3],
                          const double o[3], const double d[3], double tmax) {
    double base_t;
    int base = inv_exact_box_hit(lo, hi, o, d, tmax, &base_t);
    double op[3], dp[3];
    for (int c = 0; c < 3; c++) {
        for (int sgn = 0; sgn < 2; sgn++) {
            double dir = sgn ? -1.0 : 1.0;
            memcpy(op, o, sizeof op); memcpy(dp, d, sizeof dp);
            op[c] = fnext(o[c], dir);
            double t;
            if (inv_exact_box_hit(lo, hi, op, dp, tmax, &t) != base) return 1;
            memcpy(op, o, sizeof op); memcpy(dp, d, sizeof dp);
            dp[c] = fnext(d[c], dir);
            if (inv_exact_box_hit(lo, hi, op, dp, tmax, &t) != base) return 1;
        }
    }
    return 0;
}

typedef struct { long clear_hits, clear_misses, ambig, false_miss, false_hit;
                 double worst_t_err; } GeomAcc;

static void geom_case(GeomAcc* a, const InvModel* m,
                      const double lo[3], const double hi[3],
                      double ox, double oy, double oz,
                      double dx, double dy, double dz) {
    /* snap to float-representable doubles so the adapter's cast is lossless */
    double o[3] = { fsnap(ox), fsnap(oy), fsnap(oz) };
    double d[3] = { fsnap(dx), fsnap(dy), fsnap(dz) };
    const double tmax = 1e30;
    double te;
    int he = inv_exact_box_hit(lo, hi, o, d, tmax, &te);
    if (geom_ambiguous(lo, hi, o, d, tmax)) { a->ambig++; return; }

    double tf = 0.0;
    int hf = m->box_hit(lo, hi, o, d, tmax, &tf);
    if (he) {
        a->clear_hits++;
        if (!hf) { a->false_miss++; return; }
        double err = fabs(tf - te) / fmax(1.0, fabs(te));
        if (err > a->worst_t_err) a->worst_t_err = err;
    } else {
        a->clear_misses++;
        if (hf) a->false_hit++;
    }
}

void inv_geometry(const InvModel* m, InvReport* rep) {
    GeomAcc a; memset(&a, 0, sizeof a);
    const double ulo[3] = { -1, -1, -1 }, uhi[3] = { 1, 1, 1 };

    /* (a) random rays, unit box -- off the ambiguity band this is a pure
       differential test against exact double arithmetic */
    unsigned long s = 0x9E3779B97F4A7C15ull;
    for (int k = 0; k < TH_GEOM_RAYS; k++) {
        double o[3], d[3];
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        o[0] = ((s >> 11) / 9007199254740992.0) * 6.0 - 3.0;
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        o[1] = ((s >> 11) / 9007199254740992.0) * 6.0 - 3.0;
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        o[2] = ((s >> 11) / 9007199254740992.0) * 6.0 - 3.0;
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        d[0] = ((s >> 11) / 9007199254740992.0) * 2.0 - 1.0;
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        d[1] = ((s >> 11) / 9007199254740992.0) * 2.0 - 1.0;
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        d[2] = ((s >> 11) / 9007199254740992.0) * 2.0 - 1.0;
        double L = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        if (L < 1e-3) continue;
        for (int i = 0; i < 3; i++) d[i] /= L;
        geom_case(&a, m, ulo, uhi, o[0], o[1], o[2], d[0], d[1], d[2]);
    }

    /* (b) boundary/grazing suite: rays that pass through the box centre from
       very far away, so they genuinely hit with a large interior margin, yet
       one direction component is tiny (< the 1e-8 class of parallel cutoffs).
       These are exactly the "boundary sampling hit rate" cases. */
    const double dist[] = { 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9 };
    const double y0[]   = { 0.2, 0.6, 1.2, 1.5, 1.9 };
    const double z0[]   = { -0.5, 0.0, 0.5, 1.5 };
    for (int i = 0; i < 7; i++)
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 4; k++) {
                double D = dist[i];
                double dx = 1.0, dy = -y0[j] / D, dz = 0.0; /* hits (0,0,z0) */
                geom_case(&a, m, ulo, uhi, -D, y0[j], z0[k], dx, dy, dz);
            }

    /* (c) tangency suite: rays that graze a face/edge/corner exactly.  These
       land inside the ambiguity band on purpose -- the point is to prove the
       band is exercised and that neither the float control nor the golden
       reference is *penalised* for disagreeing on them.  A check whose band
       is never populated cannot claim to police the boundary. */
    {
        double e = 1e-3;   /* exact float-representable face coordinates */
        struct { double o[3], d[3]; } tg[] = {
            /* graze the +x face, direction entirely in y: never enters */
            {{ -2,  1 + e, 0 }, { 1, 0, 0 }},
            {{ -2,  1 - e, 0 }, { 1, 0, 0 }},
            /* graze the top edge running along z */
            {{ -2,  1 + e, 0 }, { 1, 0, 0 }},
            {{ -2, -1 - e, 0 }, { 1, 0, 0 }},
            /* exactly along a face plane */
            {{ -2,  1.0, 0 }, { 1, 0, 0 }},
            /* corner-grazing in 3D */
            {{ -2,  1 + e, 1 + e }, { 1, 0, 0 }},
            {{ -2,  1 - e, 1 - e }, { 1, 0, 0 }},
            /* tiny direction into the box, origin just outside a slab */
            {{ -2,  1 + e, 0 }, { 1, -e, 0 }},
            {{ -2,  1 + e, 0 }, { 1, -2 * e, 0 }},
            /* rays starting on the face and pointing outwards */
            {{  1,  0.5, 0.5 }, { 1, 0, 0 }},
            {{  1,  0.5, 0.5 }, { -1, 0, 0 }},
            {{ -1,  0.5, 0.5 }, { -1, 0, 0 }},
        };
        int nt = (int)(sizeof tg / sizeof tg[0]);
        for (int i = 0; i < nt; i++)
            geom_case(&a, m, ulo, uhi, tg[i].o[0], tg[i].o[1], tg[i].o[2],
                      tg[i].d[0], tg[i].d[1], tg[i].d[2]);
        addf(rep, "geometry", "ambiguity_band_exercised",
             a.ambig > 0, (double)a.ambig, 1.0,
             "tangency cases landing in the ambiguity band: %ld", a.ambig);
    }

    long decided = a.clear_hits + a.clear_misses;
    double recall_err = a.clear_hits ? (double)a.false_miss / a.clear_hits : 0.0;
    double false_hit = a.clear_misses ? (double)a.false_hit / a.clear_misses : 0.0;
    long mismatch = a.false_miss + a.false_hit;

    addf(rep, "geometry", "box_hit_vs_exact",
         mismatch <= TH_GEOM_MISMATCH, (double)mismatch, TH_GEOM_MISMATCH,
         "clear mismatches=%ld/%ld (ambig skipped=%ld)", mismatch, decided, a.ambig);
    addf(rep, "geometry", "boundary_hit_recall",
         recall_err <= (1.0 - TH_GEOM_RECALL), recall_err, 1.0 - TH_GEOM_RECALL,
         "false misses=%ld of %ld clear hits (recall=%.6f)",
         a.false_miss, a.clear_hits, 1.0 - recall_err);
    addf(rep, "geometry", "boundary_false_hit",
         false_hit <= TH_GEOM_FALSEHIT, false_hit, TH_GEOM_FALSEHIT,
         "false hits=%ld of %ld clear misses", a.false_hit, a.clear_misses);
    addf(rep, "geometry", "hit_t_accuracy",
         a.worst_t_err <= TH_GEOM_T_REL, a.worst_t_err, TH_GEOM_T_REL,
         "max relative |t-t_exact| = %.3e", a.worst_t_err);

    /* (c) acceleration-structure containment: every node box must contain its
       children / primitives (a gap is a place a ray can tunnel through). */
    InvTree t; memset(&t, 0, sizeof t);
    if (m->build_tree(&t) != 0 || t.nnodes <= 0 || t.nprims <= 0) {
        addf(rep, "geometry", "tree_build",
             0, 1.0, 0.0, "build_tree failed or produced an empty tree");
        return;
    }
    InvBox scene[INV_MAX_PRIM];
    int nprim_scene = inv_scene_boxes(scene);
    InvBox whole; for (int c = 0; c < 3; c++) { whole.lo[c] = 1e30; whole.hi[c] = -1e30; }
    for (int i = 0; i < nprim_scene; i++)
        for (int c = 0; c < 3; c++) {
            if (scene[i].lo[c] < whole.lo[c]) whole.lo[c] = scene[i].lo[c];
            if (scene[i].hi[c] > whole.hi[c]) whole.hi[c] = scene[i].hi[c];
        }
    double root_err = 0.0;
    for (int c = 0; c < 3; c++) {
        double ext = fmax(1e-6, whole.hi[c] - whole.lo[c]);
        root_err = fmax(root_err, fabs(t.node[0].box.lo[c] - whole.lo[c]) / ext);
        root_err = fmax(root_err, fabs(t.node[0].box.hi[c] - whole.hi[c]) / ext);
    }
    addf(rep, "geometry", "tree_root_covers_scene",
         root_err <= TH_TREE_TOL_REL, root_err, TH_TREE_TOL_REL,
         "relative root-box error=%.3e", root_err);

    double overhang = 0.0; long bad_leaf = 0, bad_ref = 0;
    int covered[INV_MAX_PRIM]; memset(covered, 0, sizeof covered);
    for (int i = 0; i < t.nnodes; i++) {
        InvNode* nd = &t.node[i];
        if (nd->left >= 0) {
            int kids[2] = { nd->left, nd->right };
            for (int q = 0; q < 2; q++) {
                if (kids[q] < 0 || kids[q] >= t.nnodes) { bad_ref++; continue; }
                InvBox kb = t.node[kids[q]].box;
                for (int c = 0; c < 3; c++) {
                    double ext = fmax(1e-6, nd->box.hi[c] - nd->box.lo[c]);
                    overhang = fmax(overhang, (nd->box.lo[c] - kb.lo[c]) / ext);
                    overhang = fmax(overhang, (kb.hi[c] - nd->box.hi[c]) / ext);
                }
            }
        } else {
            if (nd->first < 0 || nd->count < 0 || nd->first + nd->count > t.nprims) { bad_leaf++; continue; }
            for (int k = nd->first; k < nd->first + nd->count; k++) {
                if (k >= 0 && k < INV_MAX_PRIM) covered[k]++;
                InvBox pb = t.prim[k];
                for (int c = 0; c < 3; c++) {
                    double ext = fmax(1e-6, nd->box.hi[c] - nd->box.lo[c]);
                    overhang = fmax(overhang, (nd->box.lo[c] - pb.lo[c]) / ext);
                    overhang = fmax(overhang, (pb.hi[c] - nd->box.hi[c]) / ext);
                }
            }
        }
    }
    long part_bad = 0;
    for (int i = 0; i < t.nprims && i < INV_MAX_PRIM; i++) if (covered[i] != 1) part_bad++;

    addf(rep, "geometry", "tree_child_containment",
         overhang <= TH_TREE_TOL_REL && bad_ref == 0 && bad_leaf == 0,
         fmax(overhang, (double)(bad_ref + bad_leaf)), TH_TREE_TOL_REL,
         "max overhang=%.3e bad_refs=%ld bad_leaves=%ld", overhang, bad_ref, bad_leaf);
    addf(rep, "geometry", "tree_partition_exact",
         part_bad == 0, (double)part_bad, 0.0,
         "primitives not referenced exactly once: %ld/%d", part_bad, t.nprims);
}

/* ================================================================ SAMPLING */

static int cmp_d(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

void inv_sampling(const InvModel* m, InvReport* rep) {
    const unsigned bases[] = { 2, 3, 5 };
    const unsigned high[] = { 1000000u, 10000000u, 100000000u, 1000000000u,
                              2000000000u, 3000000000u, 4000000000u, 4294967295u };
    long out_of_domain = 0, domain_hard = 0;
    unsigned N = TH_SAMPLE_N_UNIF;

    /* Domain: every sample must be a valid probability coordinate.  A float
     * implementation cannot represent 1 - 2^-32 (the true value at the last
     * base-2 index), so exactly 1.0 may appear as a rounding of "almost 1".
     * That is a representation limit, not a defect, and 1.0 is a harmless
     * coordinate in practice.  Anything *above* 1.0 -- which is what an index
     * overflow produces (e.g. 2.09) -- has no legitimate explanation, so it is
     * counted separately and is the hard failure. */
    for (int b = 0; b < 3; b++) {
        for (unsigned i = 0; i < N; i++) {
            double v = m->halton_1d(i, bases[b]);
            if (!(v >= 0.0)) { out_of_domain++; domain_hard++; }
            else if (v > 1.0)  { out_of_domain++; domain_hard++; }
            else if (v == 1.0) { out_of_domain++; }   /* only as float rounding */
        }
        for (int i = 0; i < 8; i++) {
            double v = m->halton_1d(high[i], bases[b]);
            if (!(v >= 0.0) || v > 1.0) { out_of_domain++; domain_hard++; }
            else if (v == 1.0) { out_of_domain++; }
        }
    }
    addf(rep, "sampling", "sample_domain",
         domain_hard <= TH_SAMPLE_DOMAIN, (double)domain_hard, (double)TH_SAMPLE_DOMAIN,
         "samples >1 or <0: %ld ; samples rounding to exactly 1.0: %ld",
         domain_hard, out_of_domain - domain_hard);

    /* uniformity: mean 1/2, variance 1/12 over a full block */
    double worst_mean = 0.0, worst_var = 0.0;
    for (int b = 0; b < 2; b++) {
        double s = 0.0, s2 = 0.0;
        for (unsigned i = 0; i < N; i++) { double v = m->halton_1d(i, bases[b]); s += v; s2 += v * v; }
        double mean = s / N, var = s2 / N - mean * mean;
        worst_mean = fmax(worst_mean, fabs(mean - 0.5));
        worst_var  = fmax(worst_var,  fabs(var - 1.0 / 12.0));
    }
    addf(rep, "sampling", "uniform_mean",
         worst_mean <= TH_SAMPLE_MEAN, worst_mean, TH_SAMPLE_MEAN,
         "max |mean-1/2| = %.3e", worst_mean);
    addf(rep, "sampling", "uniform_variance",
         worst_var <= TH_SAMPLE_VAR, worst_var, TH_SAMPLE_VAR,
         "max |var-1/12| = %.3e", worst_var);

    /* 1-D star discrepancy D* over a block */
    unsigned M = TH_SAMPLE_N_DISC;
    double* xs = (double*)malloc(sizeof(double) * M);
    double worst_d = 0.0;
    if (xs) {
        for (int b = 0; b < 2; b++) {
            for (unsigned i = 0; i < M; i++) xs[i] = m->halton_1d(i, bases[b]);
            qsort(xs, M, sizeof(double), cmp_d);
            /* For sorted 0-based x_i, D* = max_i max( x_i - i/M , (i+1)/M - x_i ).
               Both sides are needed: a one-sided test would let an upward drift
               through as easily as a downward one. */
            double D = 0.0;
            for (unsigned i = 0; i < M; i++) {
                double d1 = fabs(xs[i] - (double)i / M);
                double d2 = fabs((double)(i + 1) / M - xs[i]);
                double e = d1 > d2 ? d1 : d2;
                if (e > D) D = e;
            }
            if (D > worst_d) worst_d = D;
        }
        free(xs);
    }
    addf(rep, "sampling", "discrepancy",
         worst_d <= TH_SAMPLE_DISC, worst_d, TH_SAMPLE_DISC,
         "max 1-D star discrepancy D* = %.4e", worst_d);

    /* distinctness: a block of M samples must contain ~M distinct values */
    unsigned U = 4096;
    double* ys = (double*)malloc(sizeof(double) * U);
    double uniq_err = 1.0;
    if (ys) {
        for (unsigned i = 0; i < U; i++) ys[i] = m->halton_1d(i, 2);
        qsort(ys, U, sizeof(double), cmp_d);
        unsigned uniq = 1;
        for (unsigned i = 1; i < U; i++) if (ys[i] != ys[i-1]) uniq++;
        uniq_err = 1.0 - (double)uniq / U;
        free(ys);
    }
    addf(rep, "sampling", "distinctness",
         uniq_err <= (1.0 - TH_SAMPLE_DISTINCT), uniq_err, 1.0 - TH_SAMPLE_DISTINCT,
         "1 - unique/N = %.4e", uniq_err);

    /* consistency: MC estimate of ∫_0^1 u² du must be 1/3 */
    double worst_int = 0.0;
    for (int b = 0; b < 2; b++) {
        double acc = 0.0;
        for (unsigned i = 0; i < N; i++) { double v = m->halton_1d(i, bases[b]); acc += v * v; }
        worst_int = fmax(worst_int, fabs(acc / N - 1.0 / 3.0));
    }
    addf(rep, "sampling", "integral_consistency",
         worst_int <= TH_SAMPLE_CONSIST, worst_int, TH_SAMPLE_CONSIST,
         "max |MC(u^2)-1/3| = %.3e", worst_int);
}

/* ================================================================ SYMMETRY */

/*
 * Reflection symmetry available in these modules:
 *   - the cosine sampler is azimuthally symmetric: with φ = 2πu2, mapping
 *     u2 -> 1-u2 mirrors the direction across the xz-plane
 *       s(u1,u2) and reflect_xz(s(u1,1-u2)) must coincide;
 *     hence E[x]=E[y]=0 and Var[x]=Var[y].
 *   - GGX D(h) depends on |n·h| only, so D(h)=D(-h).
 * (BRDF reciprocity is vacuous here: the only BRDF is constant Lambert.)
 */
void inv_symmetry(const InvModel* m, InvReport* rep) {
    int K = 181;
    double sx = 0, sy = 0, sxx = 0, syy = 0; long n = 0;
    double worst_mirror = 0.0;

    for (int i = 0; i < K; i++) {
        double u1 = (i + 0.5) / K;
        for (int j = 0; j < K; j++) {
            double u2 = (j + 0.5) / K;
            double s1[3], s2[3];
            m->cosine_sample(u1, u2, s1);
            m->cosine_sample(u1, 1.0 - u2, s2);
            sx += s1[0]; sy += s1[1];
            sxx += s1[0] * s1[0]; syy += s1[1] * s1[1];
            n++;
            /* reflection about the xz-plane: x,z same, y negated */
            double e = fabs(s1[0] - s2[0]);
            if (fabs(s1[2] - s2[2]) > e) e = fabs(s1[2] - s2[2]);
            if (fabs(s1[1] + s2[1]) > e) e = fabs(s1[1] + s2[1]);
            if (e > worst_mirror) worst_mirror = e;
        }
    }
    double mx = sx / n, my = sy / n;
    double vx = fmax(0.0, sxx / n - mx * mx), vy = fmax(0.0, syy / n - my * my);
    double varrel = fabs(vx - vy) / fmax(1e-12, fmax(vx, vy));
    double mean_err = fmax(fabs(mx), fabs(my));

    addf(rep, "symmetry", "sampler_azimuthal_mean",
         mean_err <= TH_SYM_MEAN, mean_err, TH_SYM_MEAN,
         "max(|E[x]|,|E[y]|)=%.3e", mean_err);
    addf(rep, "symmetry", "sampler_var_equal",
         varrel <= TH_SYM_VARREL, varrel, TH_SYM_VARREL,
         "|Var[x]-Var[y]|/max = %.3e", varrel);
    addf(rep, "symmetry", "sampler_mirror_pointwise",
         worst_mirror <= TH_SYM_POINT, worst_mirror, TH_SYM_POINT,
         "max |s(u)-reflect s(1-u)| = %.3e", worst_mirror);

    const double rough[] = { 0.1, 0.3, 0.5, 0.8, 1.0 };
    const double nh[] = { 0.05, 0.2, 0.45, 0.7, 0.95, 1.0 };
    double worst_even = 0.0;
    for (int r = 0; r < 5; r++)
        for (int i = 0; i < 6; i++) {
            double a = m->ggx_D(nh[i], rough[r]);
            double b = m->ggx_D(-nh[i], rough[r]);
            double e = fabs(a - b) / (1.0 + fabs(a));
            if (!isfinite(e)) e = 1e9;
            if (e > worst_even) worst_even = e;
        }
    addf(rep, "symmetry", "ggx_even_in_nh",
         worst_even <= TH_SYM_GGX_EVEN, worst_even, TH_SYM_GGX_EVEN,
         "max |D(h)-D(-h)|/(1+|D|) = %.3e", worst_even);
}
