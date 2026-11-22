
/* material: Lambert + GGX BRDF evaluation and cosine-weighted sampling,
   with a PDF check against the analytic cosine distribution.

   Also runs a reciprocity audit: sample the BRDF in two directions, swap
   them, and compare f(wi, wo) with f(wo, wi). Reciprocal models are
   evaluated in a swap-symmetric form, so their required deviation is an
   exact zero, not a tolerance. Models that are non-reciprocal by
   construction are reported under their own class and never counted as
   failures. A direction argument that is not unit length is reported with
   its role, vector and length before any evaluation happens. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI_F 3.14159265f
#define PI_D 3.14159265358979323846
#define UNIT_TOL 1e-4f
#define RECIP_MAX_PAIRS 4096
#define RECIP_BINS 5

typedef struct { float x, y, z; } V3;

static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}

static V3 add3(V3 a, V3 b) { V3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static V3 scale3(V3 a, float s) { V3 r = { a.x * s, a.y * s, a.z * s }; return r; }

/* dot product in double: IEEE multiply/add is commutative, so dot(a,b) and
   dot(b,a) are bit-identical, which the reciprocity audit relies on. */
static double dot3d(V3 a, V3 b) {
    return (double)a.x * (double)b.x + (double)a.y * (double)b.y + (double)a.z * (double)b.z;
}
static double len3(V3 a) { return sqrt(dot3d(a, a)); }

/* Lambert: albedo / pi */
static float brdf_lambert(float albedo) { return albedo / PI_F; }

/* GGX normal distribution */
static float ggx_ndf(float nh, float rough) {
    float a = rough * rough;
    float a2 = a * a;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / (PI_F * d * d);
}

/* cosine-weighted hemisphere: pdf = cos/pi */
static float cosine_pdf(float cos_theta) { return cos_theta / PI_F; }

static V3 cosine_sample(float u1, float u2) {
    float r = sqrtf(u1);
    float phi = 2.0f * PI_F * u2;
    V3 v = { r * cosf(phi), r * sinf(phi), sqrtf(1.0f - u1 > 0 ? 1.0f - u1 : 0) };
    return norm3(v);
}

/* ------------------------------------------------------------------ */
/* Materials and BRDF evaluation                                      */
/* ------------------------------------------------------------------ */

typedef enum { RECIP_ANALYTIC = 0, RECIP_NONRECIPROCAL = 1 } RecipClass;

typedef struct {
    const char* name;
    RecipClass cls;
    int ggx;
    float albedo;
    float rough;
    float f0;
    float sheen;   /* forward-scatter gain; breaks reciprocity when != 0 */
    float defect;  /* self-test only: deliberate directional asymmetry */
} Material;

static const Material MAT_LAMBERT = { "lambert_analytic",  RECIP_ANALYTIC,    0, 0.80f, 0.00f, 0.00f, 0.00f, 0.00f };
static const Material MAT_GGX_010 = { "ggx_r0.10_analytic", RECIP_ANALYTIC,    1, 0.00f, 0.10f, 0.04f, 0.00f, 0.00f };
static const Material MAT_GGX_035 = { "ggx_r0.35_analytic", RECIP_ANALYTIC,    1, 0.00f, 0.35f, 0.04f, 0.00f, 0.00f };
static const Material MAT_GGX_070 = { "ggx_r0.70_analytic", RECIP_ANALYTIC,    1, 0.00f, 0.70f, 0.04f, 0.00f, 0.00f };
static const Material MAT_GGX_100 = { "ggx_r1.00_analytic", RECIP_ANALYTIC,    1, 0.00f, 1.00f, 0.04f, 0.00f, 0.00f };
static const Material MAT_SHEEN   = { "sheen_forward_nr",  RECIP_NONRECIPROCAL, 1, 0.00f, 0.45f, 0.04f, 0.45f, 0.00f };

/* Sensitivity probe: claims analytic reciprocity but carries a directional
   term, so the audit MUST report FAIL. Its line in the output is expected to
   say FAIL_nonzero_deviation; the self-test asserts exactly that, which proves
   the audit is not vacuous. */
static const Material MAT_SELFTEST_DEFECT = { "selftest_directional_defect", RECIP_ANALYTIC, 1, 0.00f, 0.35f, 0.04f, 0.00f, 0.10f };

/* Exact Smith masking, alpha2 = (rough^2)^2. */
static double smith_g1(double cos_t, double alpha2) {
    double c2 = cos_t * cos_t;
    return 2.0 * cos_t / (cos_t + sqrt(alpha2 + (1.0 - alpha2) * c2));
}

/* Cook-Torrance GGX. Everything that is mathematically symmetric under
   wi <-> wo is computed from one shared value, so both argument orders go
   through the exact same floating-point operations. */
static double ggx_cook_torrance(const Material* m, V3 wi, V3 wo) {
    double ci = (double)wi.z;
    double co = (double)wo.z;
    double alpha = (double)m->rough * (double)m->rough;
    double alpha2 = alpha * alpha;

    V3 h = norm3(add3(wi, wo));   /* wi + wo is exactly wo + wi */
    double nh = (double)h.z;
    if (nh <= 0.0) return 0.0;

    /* Analytically wi.h == wo.h; averaging the two numerically evaluated
       values keeps the model analytic and the swapped result bit-identical. */
    double ch = 0.5 * (dot3d(wi, h) + dot3d(wo, h));
    if (ch < 0.0) ch = 0.0;
    if (ch > 1.0) ch = 1.0;

    double dd = nh * nh * (alpha2 - 1.0) + 1.0;
    double d = alpha2 / (PI_D * dd * dd);
    double g = smith_g1(ci, alpha2) * smith_g1(co, alpha2);   /* commutes exactly */
    double omc = 1.0 - ch;
    double f = (double)m->f0 + (1.0 - (double)m->f0) * omc * omc * omc * omc * omc;
    double denom = 4.0 * (ci * co);                           /* commutes exactly */
    if (denom <= 0.0) return 0.0;
    return d * g * f / denom;
}

/* f(wi, wo) in the local shading frame (+z = normal). Directions are
   expected to be unit length; use brdf_eval_checked() for untrusted input. */
static double brdf_eval(const Material* m, V3 wi, V3 wo) {
    if (wi.z <= 0.0f || wo.z <= 0.0f) return 0.0;
    if (!m->ggx) return (double)m->albedo / PI_D;
    double f = ggx_cook_torrance(m, wi, wo);
    if (m->sheen != 0.0f) f *= 1.0 + (double)m->sheen * (double)wo.z;
    if (m->defect != 0.0f) f *= 1.0 + (double)m->defect * ((double)wo.z - (double)wi.z);
    return f;
}

/* ------------------------------------------------------------------ */
/* Direction validation                                               */
/* ------------------------------------------------------------------ */

typedef enum {
    DIR_OK = 0,
    DIR_WI_NOT_UNIT,
    DIR_WO_NOT_UNIT,
    DIR_BOTH_NOT_UNIT
} DirStatus;

typedef struct {
    DirStatus status;
    V3 wi, wo;
    double len_i, len_o;
} DirReport;

static const char* dir_status_name(DirStatus s) {
    switch (s) {
        case DIR_OK: return "ok";
        case DIR_WI_NOT_UNIT: return "wi_not_unit";
        case DIR_WO_NOT_UNIT: return "wo_not_unit";
        case DIR_BOTH_NOT_UNIT: return "both_not_unit";
    }
    return "unknown";
}

static int nearly_unit(double len) { return fabs(len - 1.0) <= (double)UNIT_TOL; }

static DirReport check_directions(V3 wi, V3 wo) {
    DirReport r;
    memset(&r, 0, sizeof(r));
    r.wi = wi;
    r.wo = wo;
    r.len_i = len3(wi);
    r.len_o = len3(wo);
    int bad_i = !nearly_unit(r.len_i);
    int bad_o = !nearly_unit(r.len_o);
    r.status = (bad_i && bad_o) ? DIR_BOTH_NOT_UNIT
             : bad_i ? DIR_WI_NOT_UNIT
             : bad_o ? DIR_WO_NOT_UNIT
             : DIR_OK;
    return r;
}

/* Report the offending direction(s): role, vector and length. */
static void report_directions(const char* material, const DirReport* r) {
    if (r->status == DIR_OK) return;
    if (r->status == DIR_WI_NOT_UNIT || r->status == DIR_BOTH_NOT_UNIT)
        printf("direction_not_unit material=%s status=%s direction=wi vector=(%.6f, %.6f, %.6f) length=%.6f tolerance=%.0e\n",
               material, dir_status_name(r->status),
               (double)r->wi.x, (double)r->wi.y, (double)r->wi.z, r->len_i, (double)UNIT_TOL);
    if (r->status == DIR_WO_NOT_UNIT || r->status == DIR_BOTH_NOT_UNIT)
        printf("direction_not_unit material=%s status=%s direction=wo vector=(%.6f, %.6f, %.6f) length=%.6f tolerance=%.0e\n",
               material, dir_status_name(r->status),
               (double)r->wo.x, (double)r->wo.y, (double)r->wo.z, r->len_o, (double)UNIT_TOL);
}

/* Validate, then evaluate. Invalid input is refused (returns 0) and the
   report carries the offending direction for the caller to print. */
static double brdf_eval_checked(const Material* m, V3 wi, V3 wo, DirReport* rep) {
    *rep = check_directions(wi, wo);
    if (rep->status != DIR_OK) return 0.0;
    return brdf_eval(m, wi, wo);
}

/* ------------------------------------------------------------------ */
/* BRDF sampling                                                      */
/* ------------------------------------------------------------------ */

static V3 sample_ggx_half(float u1, float u2, float rough) {
    float alpha = rough * rough;
    float alpha2 = alpha * alpha;
    float phi = 2.0f * PI_F * u1;
    float denom = 1.0f + (alpha2 - 1.0f) * u2;
    float ct = sqrtf((1.0f - u2) / denom);
    float st = sqrtf(fmaxf(0.0f, 1.0f - ct * ct));
    V3 h = { st * cosf(phi), st * sinf(phi), ct };
    return norm3(h);
}

static V3 reflect_about(V3 wi, V3 h) {
    double d = dot3d(wi, h);
    V3 r = { (float)(2.0 * d * h.x - wi.x),
             (float)(2.0 * d * h.y - wi.y),
             (float)(2.0 * d * h.z - wi.z) };
    return norm3(r);
}

/* Sample an outgoing direction from the material's BRDF given wi. */
static V3 brdf_sample(const Material* m, V3 wi, float u1, float u2) {
    if (!m->ggx) {
        (void)wi;                 /* isotropic diffuse: around the normal */
        return cosine_sample(u1, u2);
    }
    V3 h = sample_ggx_half(u1, u2, m->rough);
    return reflect_about(wi, h);
}

/* ------------------------------------------------------------------ */
/* Reciprocity audit                                                  */
/* ------------------------------------------------------------------ */

typedef struct {
    int attempts, evaluated, skipped;
    int nonzero;         /* pairs whose two evaluations differ in any bit */
    int positive;        /* pairs where at least one evaluation is non-zero */
    double max_abs_dev, max_rel_dev;
    double p50_rel, p90_rel, p99_rel;
    int bins[RECIP_BINS];
} RecipStats;

static double rel_devs[RECIP_MAX_PAIRS];

static float frand(unsigned* state) {
    *state = *state * 1664525u + 1013904223u;
    return (float)((*state >> 8) * (1.0 / 16777216.0));
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a;
    double y = *(const double*)b;
    return (x > y) - (x < y);
}

static double percentile(const double* v, int n, double p) {
    if (n <= 0) return 0.0;
    int idx = (int)((n - 1) * p + 0.5);
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;
    return v[idx];
}

static int classify_bin(double abs_dev, double rel_dev) {
    if (abs_dev == 0.0) return 0;
    if (rel_dev <= 1e-6) return 1;
    if (rel_dev <= 1e-3) return 2;
    if (rel_dev <= 1e-1) return 3;
    return 4;
}

/* Generate direction pairs by sampling the BRDF itself: a first direction
   around the normal, a second direction from the BRDF given the first.
   Then evaluate both argument orders and collect the deviations. */
static RecipStats reciprocity_scan(const Material* m, int attempts, unsigned seed) {
    RecipStats s;
    memset(&s, 0, sizeof(s));
    if (attempts > RECIP_MAX_PAIRS) attempts = RECIP_MAX_PAIRS;
    s.attempts = attempts;

    unsigned state = seed;
    V3 n = { 0, 0, 1 };

    for (int i = 0; i < attempts; i++) {
        float u1 = frand(&state), u2 = frand(&state);
        float u3 = frand(&state), u4 = frand(&state);

        V3 wi = brdf_sample(m, n, u1, u2);
        V3 wo = brdf_sample(m, wi, u3, u4);
        if (wi.z <= 0.0f || wo.z <= 0.0f) { s.skipped++; continue; }
        if (check_directions(wi, wo).status != DIR_OK) { s.skipped++; continue; }

        double fab = brdf_eval(m, wi, wo);
        double fba = brdf_eval(m, wo, wi);
        if (!isfinite(fab) || !isfinite(fba)) { s.skipped++; continue; }

        double abs_dev = fabs(fab - fba);
        double den = fabs(fab) + fabs(fba);
        double rel_dev = den > 0.0 ? 2.0 * abs_dev / den : 0.0;

        if (fab > 0.0 || fba > 0.0) s.positive++;
        if (abs_dev > 0.0) s.nonzero++;
        if (abs_dev > s.max_abs_dev) s.max_abs_dev = abs_dev;
        if (rel_dev > s.max_rel_dev) s.max_rel_dev = rel_dev;
        s.bins[classify_bin(abs_dev, rel_dev)]++;

        rel_devs[s.evaluated++] = rel_dev;
    }

    qsort(rel_devs, (size_t)s.evaluated, sizeof(double), cmp_double);
    s.p50_rel = percentile(rel_devs, s.evaluated, 0.50);
    s.p90_rel = percentile(rel_devs, s.evaluated, 0.90);
    s.p99_rel = percentile(rel_devs, s.evaluated, 0.99);
    return s;
}

/* A material declared non-reciprocal cannot be judged reciprocal, but it is
   not an error either: the audit's job is to label it. It is "detected" when
   the scan found deviations on pairs that carry signal (one side non-zero),
   rather than underflow noise from pairs that landed at the horizon. */
static int nonreciprocal_detected(const RecipStats* s) {
    return s->positive > 0 && s->nonzero > 0;
}

/* Print the deviation distribution and max deviation, and return whether the
   audit got what it expects: exactly zero for a reciprocal analytic model, a
   detected deviation for a declared non-reciprocal one. */
static int reciprocity_report(const Material* m, const RecipStats* s) {
    int analytic = (m->cls == RECIP_ANALYTIC);
    /* An empty scan must never read as "exactly zero": require real pairs. */
    int exact_zero = (s->evaluated > 0 && s->nonzero == 0 && s->max_abs_dev == 0.0);
    const char* verdict;
    if (analytic) verdict = exact_zero ? "exact_zero" : "FAIL_nonzero_deviation";
    else if (nonreciprocal_detected(s)) verdict = "labeled_nonreciprocal";
    else verdict = "labeled_nonreciprocal_undetected";

    printf("reciprocity material=%s class=%s pairs=%d evaluated=%d skipped=%d "
           "max_abs_dev=%.6e max_rel_dev=%.6e "
           "p50_rel=%.6e p90_rel=%.6e p99_rel=%.6e "
           "distribution=exact:%d,le1e-6:%d,le1e-3:%d,le1e-1:%d,gt1e-1:%d verdict=%s\n",
           m->name,
           analytic ? "reciprocal_analytic" : "nonreciprocal",
           s->attempts, s->evaluated, s->skipped,
           s->max_abs_dev, s->max_rel_dev,
           s->p50_rel, s->p90_rel, s->p99_rel,
           s->bins[0], s->bins[1], s->bins[2], s->bins[3], s->bins[4],
           verdict);

    return analytic ? exact_zero : nonreciprocal_detected(s);
}

static int run_reciprocity(const Material* m, unsigned seed) {
    RecipStats s = reciprocity_scan(m, RECIP_MAX_PAIRS, seed);
    return reciprocity_report(m, &s);
}

/* ------------------------------------------------------------------ */
/* Self test                                                          */
/* ------------------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    check("lambert_unit", fabsf(brdf_lambert(PI_F) - 1.0f) < 1e-4f);
    check("lambert_zero", brdf_lambert(0.0f) == 0.0f);

    check("ggx_peak", ggx_ndf(1.0f, 0.5f) > ggx_ndf(0.2f, 0.5f));
    check("ggx_positive", ggx_ndf(0.5f, 0.3f) > 0.0f);

    check("pdf_half", fabsf(cosine_pdf(0.5f) - 0.5f / PI_F) < 1e-6f);

    V3 s = cosine_sample(0.5f, 0.25f);
    check("sample_unit", fabsf(sqrtf(dot3(s, s)) - 1.0f) < 1e-4f);
    check("sample_z_nonneg", s.z >= 0.0f);

    V3 up = { 0, 0, 1 };
    check("sample_hemisphere", dot3(s, up) >= -1e-6f);

    /* mean of cos over samples ~ 2/3 for cosine sampling */
    double acc = 0.0; int n = 512;
    for (int i = 0; i < n; i++) {
        float u1 = ((float)(i % 16) + 0.5f) / 16.0f;
        float u2 = ((float)(i / 16) + 0.5f) / 32.0f;
        V3 v = cosine_sample(u1, u2);
        acc += v.z;
    }
    check("mean_cos_close", fabs(acc / n - 2.0 / 3.0) < 0.02);

    /* ---- direct swap checks ---- */
    V3 a = norm3((V3){ 0.4f, 0.1f, 1.0f });
    V3 b = norm3((V3){ -0.2f, 0.5f, 0.7f });
    check("lambert_swap_bit_identical",
          brdf_eval(&MAT_LAMBERT, a, b) == brdf_eval(&MAT_LAMBERT, b, a));
    check("ggx_swap_bit_identical",
          brdf_eval(&MAT_GGX_035, a, b) == brdf_eval(&MAT_GGX_035, b, a) &&
          brdf_eval(&MAT_GGX_035, a, b) > 0.0);
    check("sheen_swap_differs",
          brdf_eval(&MAT_SHEEN, a, b) != brdf_eval(&MAT_SHEEN, b, a));

    /* ---- reciprocity audit: distribution and max deviation ---- */
    check("reciprocity_lambert_exact_zero", run_reciprocity(&MAT_LAMBERT, 0x12345678u));
    check("reciprocity_ggx_r010_exact_zero", run_reciprocity(&MAT_GGX_010, 0x9e3779b9u));
    check("reciprocity_ggx_r035_exact_zero", run_reciprocity(&MAT_GGX_035, 0x85ebca6bu));
    check("reciprocity_ggx_r070_exact_zero", run_reciprocity(&MAT_GGX_070, 0xc2b2ae35u));
    check("reciprocity_ggx_r100_exact_zero", run_reciprocity(&MAT_GGX_100, 0x27d4eb2fu));
    check("nonreciprocal_labeled_not_error", run_reciprocity(&MAT_SHEEN, 0x51ed2701u));

    /* The audit must actually fail a broken model; the probe's printed line is
       expected to read verdict=FAIL_nonzero_deviation. */
    check("audit_detects_directional_defect",
          run_reciprocity(&MAT_SELFTEST_DEFECT, 0xa5a5a5a5u) == 0);

    /* ---- un-normalized direction reporting ---- */
    DirReport rep;
    V3 wi_bad = scale3(a, 2.0f);
    double ev = brdf_eval_checked(&MAT_GGX_035, wi_bad, b, &rep);
    report_directions(MAT_GGX_035.name, &rep);
    check("unnormalized_wi_detected",
          rep.status == DIR_WI_NOT_UNIT && fabs(rep.len_i - 2.0) < 1e-5);
    check("unnormalized_wi_vector_kept",
          rep.wi.x == wi_bad.x && rep.wi.y == wi_bad.y && rep.wi.z == wi_bad.z &&
          ev == 0.0);

    V3 wo_bad = scale3(b, 0.5f);
    ev = brdf_eval_checked(&MAT_GGX_035, a, wo_bad, &rep);
    report_directions(MAT_GGX_035.name, &rep);
    check("unnormalized_wo_detected",
          rep.status == DIR_WO_NOT_UNIT && fabs(rep.len_o - 0.5) < 1e-5);
    check("unnormalized_eval_blocked", ev == 0.0);

    rep = check_directions(wi_bad, wo_bad);
    report_directions(MAT_GGX_035.name, &rep);
    check("unnormalized_both_detected", rep.status == DIR_BOTH_NOT_UNIT);

    V3 zero = { 0, 0, 0 };
    rep = check_directions(zero, b);
    report_directions(MAT_GGX_035.name, &rep);
    check("zero_direction_reported", rep.status == DIR_WI_NOT_UNIT && rep.len_i == 0.0);

    check("unit_directions_ok", check_directions(a, b).status == DIR_OK);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
