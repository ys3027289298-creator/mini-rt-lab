/*
 * invariants.h -- invariant-checker API for the mini-rt-lab renderer modules.
 *
 * The checker is a *pure observer*.  It never edits the three module sources,
 * never calls their main()/check(), and never writes module state.  It is built
 * by compiling the unmodified sources into this translation unit with their
 * entry points renamed (renderer_edge.c), then calling only their pure maths.
 *
 * Every invariant is written against an InvModel so that the same code can be
 * pointed at (a) the shipped renderer, (b) a known-correct reference, and
 * (c) hand-broken variants.  Calibration then has a meaning independent of the
 * audit: broken variants must be flagged, the correct reference must not be.
 */
#ifndef INVARIANTS_H
#define INVARIANTS_H

#include <stddef.h>
#include <stdint.h>

#define INV_PI 3.14159265358979323846

/* -------------------------------------------------------------------------
 * Decision thresholds, one block per invariant class.
 *
 * These are *decision* thresholds, not fitted constants.  Each is chosen at
 * least two orders of magnitude above the numerical noise floor of its check
 * (float rounding, MC/quadrature error) and at least one order below the
 * smallest violation the calibration mutants inject, so the same threshold
 * separates "correct" from "broken" on both sides.
 * ---------------------------------------------------------------------- */

/* energy ---------------------------------------------------------------- */
#define TH_ENERGY_INCREASE   1e-3   /* reflectance may exceed 1 by this much   */
#define TH_ENERGY_CONSIST    1e-2   /* |MC reflectance - analytic albedo|       */
#define TH_ENERGY_NEG        1e-6   /* negative reflectance / pdf allowance     */
#define TH_PDF_NORM          2e-3   /* |integral of sampling pdf over sphere-1|  */
#define TH_NDF_NORM          5e-3   /* |integral of GGX D*cos over hemisphere-1| */
#define TH_NDF_SIMPSON_N     2000000/* even; resolves the rough=0.05 lobe       */

/* geometry -------------------------------------------------------------- */
#define TH_GEOM_MISMATCH     0      /* hit/miss disagreements on clear cases    */
#define TH_GEOM_RECALL       1.0    /* fraction of clear true hits detected     */
#define TH_GEOM_FALSEHIT     0.0    /* fraction of clear misses reported as hit */
#define TH_GEOM_T_REL        2e-3   /* relative error of the reported t         */
#define TH_GEOM_MARGIN_REL   1e-2   /* ambiguity band as fraction of extent     */
#define TH_GEOM_RAYS         20000
#define TH_TREE_TOL_REL      1e-4   /* containment slack for float boxes        */

/* sampling -------------------------------------------------------------- */
#define TH_SAMPLE_DOMAIN     0      /* samples outside [0,1) allowed            */
#define TH_SAMPLE_MEAN       1e-2   /* |sample mean - 0.5|                       */
#define TH_SAMPLE_VAR        5e-3   /* |sample variance - 1/12|                  */
#define TH_SAMPLE_DISC       1e-2   /* 1-D star discrepancy D*                   */
#define TH_SAMPLE_DISTINCT   0.999  /* unique/N                                 */
#define TH_SAMPLE_CONSIST    2e-2   /* abs err of MC integral vs analytic        */
#define TH_SAMPLE_N_UNIF     65536
#define TH_SAMPLE_N_DISC     4096

/* symmetry -------------------------------------------------------------- */
#define TH_SYM_MEAN          1e-2   /* |mean x|, |mean y| of the sampler         */
#define TH_SYM_VARREL        2e-2   /* |var x - var y| / max(var x, var y)       */
#define TH_SYM_POINT         1e-5   /* pointwise mirror identity                 */
#define TH_SYM_GGX_EVEN      1e-6   /* |D(h)-D(-h)| / (1+|D|)                     */

#define INV_MAX_PRIM 32
#define INV_MAX_NODE 64

typedef struct { double lo[3], hi[3]; } InvBox;
typedef struct { InvBox box; int left, right; int first, count; } InvNode;
typedef struct {
    int nnodes, nprims;
    InvNode node[INV_MAX_NODE];
    InvBox  prim[INV_MAX_PRIM];
} InvTree;

typedef struct {
    const char* cls;      /* "energy" | "geometry" | "sampling" | "symmetry" */
    const char* name;
    int    pass;
    double measured;
    double threshold;
    char   detail[160];
} InvResult;

typedef struct { InvResult r[512]; int n; } InvReport;

/*
 * The surface the checker is allowed to observe.  Each pointer is a pure
 * function of its arguments.  box_hit uses the same three-axis slab contract
 * as the renderer (entry t >= 0, t <= tmax, ray origin may start inside).
 * build_tree fills a canonical tree for a fixed deterministic scene and
 * returns 0 on success.
 */
typedef struct {
    const char* name;
    int    (*box_hit)(const double lo[3], const double hi[3],
                      const double o[3], const double d[3],
                      double tmax, double* t);
    double (*halton_1d)(unsigned i, unsigned base);
    double (*ggx_D)(double nh, double rough);
    double (*lambert_f)(double albedo);
    double (*cosine_pdf)(double cos_theta);
    void   (*cosine_sample)(double u1, double u2, double out[3]);
    int    (*build_tree)(InvTree* out);
} InvModel;

/* exact double-precision slab test: the reference for the geometry class */
int inv_exact_box_hit(const double lo[3], const double hi[3],
                      const double o[3], const double d[3],
                      double tmax, double* t);

/* fixed scene shared by every build_tree implementation */
int  inv_scene_boxes(InvBox out[INV_MAX_PRIM]);

void inv_energy  (const InvModel* m, InvReport* rep);
void inv_geometry(const InvModel* m, InvReport* rep);
void inv_sampling(const InvModel* m, InvReport* rep);
void inv_symmetry(const InvModel* m, InvReport* rep);
void inv_run_model(const InvModel* m, InvReport* rep);

/* classification helpers over a filled report */
int inv_class_failed(const InvReport* rep, const char* cls);
int inv_class_total (const InvReport* rep, const char* cls);

/* known-correct reference (double precision; all invariants must pass) */
extern const InvModel inv_golden_reference;
/* known-correct reference using ordinary float arithmetic, i.e. the same
   precision class the shipped modules use.  It exists to prove the thresholds
   admit a *legitimate* float implementation: if this control failed, every
   failure attributed to the shipped renderer would be suspect. */
extern const InvModel inv_float_control;
/* shipped modules, using the sampler the module's own tests exercise (halton2) */
extern const InvModel inv_real_renderer;
/* shipped modules, using the *other* shipped generator halton() */
extern const InvModel inv_real_halton;

/* Calibration mutants: copies of the golden reference with exactly one
   behaviour broken.  target_cls is the invariant class that must catch it. */
typedef struct { const char* name; const char* target_cls; InvModel m; } InvMutant;
int inv_build_mutants(InvMutant* out, int cap);
extern const int INV_MUTANT_CAP;

/* purity witnesses ------------------------------------------------------ */
/* Hash of a deterministic render workload (BVH build + ray casts + samples +
   BRDF evals).  Called before and after the checks; must be bit-identical.  */
uint64_t inv_render_fingerprint(void);
/* True iff none of the modules' own self-test counters were touched.        */
int inv_module_state_clean(void);

#endif
