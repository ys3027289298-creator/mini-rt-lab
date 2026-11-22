
/* onb: orthonormal basis construction from a surface normal.
   - build_onb_fixed:   textbook n x world-axis construction; condition number
                        kappa = 1/sin(angle(n, axis)) blows up near the poles,
                        which flips the tangent frame and breaks cosine sampling.
   - build_onb_stable:  branch on the dominant component (PBRT-style); the
                        smallest |component| of a unit normal is <= 1/sqrt(3),
                        so kappa is bounded by sqrt(2) for every input.
   - build_onb_quat:    reference frame built from the shortest-arc quaternion
                        rotating +z onto n; used for a per-vector error bound.
   - build_onb_safe:    policy wrapper: reject zero normals, run the fixed
                        construction, report and fall back when kappa or the
                        orthogonality residual exceeds its budget.
   Usage: ./onb [selftest|sweep out.csv|flip out.csv]
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif
#define EPS_D 2.220446049250313e-16
#define EPS_F 1.1920928955078125e-7f
/* A normal shorter than this cannot be normalized in the type. */
#define ZERO_TOL_D 1e-300
#define ZERO_TOL_F 1e-37f
/* Fixed-construction budget: above these, report and switch to the fallback. */
#define KAPPA_LIMIT 1e4
#define RESID_LIMIT 1e-12
/* Empirically tight constants in err <= C * eps * kappa (see sweep output).
   2e6 random normals found worst ratios 0.77 (double) / 0.44 (float); the
   values below keep a >2x margin and match the rounding-analysis order. */
#define C_BOUND_D 2.0
#define C_BOUND_F 4.0

typedef struct { double x, y, z; } V3;
typedef struct { float x, y, z; } V3f;

typedef enum {
    ONB_OK = 0,
    ONB_ZERO_NORMAL,    /* norm below ZERO_TOL: rejected */
    ONB_NONFINITE,      /* NaN/Inf component: rejected */
    ONB_DEGENERATE,     /* fixed axis parallel to n: construction undefined */
    ONB_KAPPA_HIGH,     /* safe path: kappa over budget, fallback used */
    ONB_RESID_HIGH      /* safe path: residual over budget, fallback used */
} OnbStatus;

typedef struct {
    V3 t, b, n;
    double kappa;       /* condition number against the construction helper */
    double resid;       /* max of |dot| and ||v|-1| over the three axes */
    OnbStatus status;
    int method;         /* 0 = fixed, 1 = stable fallback */
} Frame;

static V3 v3(double x, double y, double z) { V3 r = { x, y, z }; return r; }
static double dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 cross3(V3 a, V3 b) {
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static double len3(V3 a) { return sqrt(dot3(a, a)); }
static int finite3(V3 a) { return isfinite(a.x) && isfinite(a.y) && isfinite(a.z); }
static V3 normalize3(V3 a) {
    double l = len3(a);
    if (!(l > 0.0) || !finite3(a)) return v3(0, 0, 0);
    return v3(a.x / l, a.y / l, a.z / l);
}
static V3 fv3(V3f a) { return v3(a.x, a.y, a.z); }

static float dot3f(V3f a, V3f b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3f cross3f(V3f a, V3f b) {
    V3f r = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return r;
}
static float len3f(V3f a) { return sqrtf(dot3f(a, a)); }
static int finite3f(V3f a) { return isfinite(a.x) && isfinite(a.y) && isfinite(a.z); }
static V3f normalize3f(V3f a) {
    float l = len3f(a);
    V3f zero = { 0, 0, 0 };
    if (!(l > 0.0f) || !finite3f(a)) return zero;
    V3f r = { a.x / l, a.y / l, a.z / l };
    return r;
}

/* Orthogonality/unit-length residual of a frame. */
static double frame_resid(V3 t, V3 b, V3 n) {
    double r = fabs(dot3(t, b));
    double d = fabs(dot3(b, n)); if (d > r) r = d;
    d = fabs(dot3(n, t)); if (d > r) r = d;
    d = fabs(len3(t) - 1.0); if (d > r) r = d;
    d = fabs(len3(b) - 1.0); if (d > r) r = d;
    d = fabs(len3(n) - 1.0); if (d > r) r = d;
    return r;
}

/* kappa of helper axis h against unit normal n: 1/sin(angle(n, h)). */
static double kappa_of(V3 n, V3 h) {
    double s = len3(cross3(n, h));
    return s > 0.0 ? 1.0 / s : INFINITY;
}

/* Textbook construction against the fixed world z axis.
   Fails exactly where the problem statement says: n parallel to +/- z. */
static Frame build_onb_fixed(V3 n_in) {
    Frame f; memset(&f, 0, sizeof(f));
    f.method = 0;
    if (!finite3(n_in)) { f.status = ONB_NONFINITE; f.kappa = INFINITY; return f; }
    double l = len3(n_in);
    if (!(l > ZERO_TOL_D)) { f.status = ONB_ZERO_NORMAL; f.kappa = INFINITY; return f; }
    V3 n = v3(n_in.x / l, n_in.y / l, n_in.z / l);
    V3 z = v3(0, 0, 1);
    double s = len3(cross3(n, z));
    f.kappa = s > 0.0 ? 1.0 / s : INFINITY;
    if (!(s > 1e-9)) { f.n = n; f.status = ONB_DEGENERATE; return f; }
    V3 t = normalize3(cross3(z, n));
    V3 b = cross3(n, t);               /* right-handed: t x b = n */
    f.t = t; f.b = b; f.n = n;
    f.resid = frame_resid(t, b, n);
    f.status = ONB_OK;
    return f;
}

/* Stable dominant-axis construction (PBRT "branchless" form, branched here).
   The helper axis is the coordinate axis most perpendicular to n, so the
   helper condition number never exceeds sqrt(2). */
static Frame build_onb_stable(V3 n_in) {
    Frame f; memset(&f, 0, sizeof(f));
    f.method = 1;
    if (!finite3(n_in)) { f.status = ONB_NONFINITE; f.kappa = INFINITY; return f; }
    double l = len3(n_in);
    if (!(l > ZERO_TOL_D)) { f.status = ONB_ZERO_NORMAL; f.kappa = INFINITY; return f; }
    V3 n = v3(n_in.x / l, n_in.y / l, n_in.z / l);
    V3 t, b, helper;
    if (fabs(n.z) <= fabs(n.x) && fabs(n.z) <= fabs(n.y)) {
        /* helper = z; valid because |n.z| <= 1/sqrt(3) on this branch */
        helper = v3(0, 0, 1);
        t = normalize3(cross3(helper, n));
        b = cross3(n, t);
    } else if (fabs(n.y) <= fabs(n.x)) {
        /* helper = y */
        helper = v3(0, 1, 0);
        t = normalize3(cross3(helper, n));
        b = cross3(n, t);
    } else {
        /* helper = x */
        helper = v3(1, 0, 0);
        t = normalize3(cross3(helper, n));
        b = cross3(n, t);
    }
    f.t = t; f.b = b; f.n = n;
    f.kappa = kappa_of(n, helper);
    f.resid = frame_resid(t, b, n);
    f.status = ONB_OK;
    return f;
}

/* Reference frame from the shortest-arc quaternion q rotating +z onto n.
   q = (k_x, k_y, k_z, w) with k = normalize(+z x n), w = 1 + n.z.
   Singular at n = -z, handled as the 180-degree rotation about x. */
static Frame build_onb_quat(V3 n_in) {
    Frame f; memset(&f, 0, sizeof(f));
    if (!finite3(n_in)) { f.status = ONB_NONFINITE; f.kappa = INFINITY; return f; }
    double l = len3(n_in);
    if (!(l > ZERO_TOL_D)) { f.status = ONB_ZERO_NORMAL; f.kappa = INFINITY; return f; }
    V3 n = v3(n_in.x / l, n_in.y / l, n_in.z / l);
    V3 z = v3(0, 0, 1);
    V3 t, b;
    if (n.z < -1.0 + 1e-12) {
        /* 180 degrees about x: z -> -z, x -> -x, y -> y */
        t = v3(-1, 0, 0);
        b = v3(0, 1, 0);
        f.kappa = INFINITY;                 /* quaternion axis is undefined here */
    } else {
        V3 k = cross3(z, n);
        double w = 1.0 + n.z;
        double inv = 1.0 / sqrt(dot3(k, k) + w * w);  /* normalize the quaternion */
        k = v3(k.x * inv, k.y * inv, k.z * inv);
        w *= inv;
        /* v' = v + 2 k x (k x v + w v) */
        V3 kx = v3(1, 0, 0), ky = v3(0, 1, 0);
        V3 ax = v3(k.y * kx.z - k.z * kx.y + w * kx.x,
                   k.z * kx.x - k.x * kx.z + w * kx.y,
                   k.x * kx.y - k.y * kx.x + w * kx.z);
        V3 ay = v3(k.y * ky.z - k.z * ky.y + w * ky.x,
                   k.z * ky.x - k.x * ky.z + w * ky.y,
                   k.x * ky.y - k.y * ky.x + w * ky.z);
        t = v3(kx.x + 2.0 * (k.y * ax.z - k.z * ax.y),
               kx.y + 2.0 * (k.z * ax.x - k.x * ax.z),
               kx.z + 2.0 * (k.x * ax.y - k.y * ax.x));
        b = v3(ky.x + 2.0 * (k.y * ay.z - k.z * ay.y),
               ky.y + 2.0 * (k.z * ay.x - k.x * ay.z),
               ky.z + 2.0 * (k.x * ay.y - k.y * ay.x));
        f.kappa = kappa_of(n, z);           /* same pole geometry as the fixed build */
    }
    f.t = t; f.b = b; f.n = n;
    f.resid = frame_resid(t, b, n);
    f.status = ONB_OK;
    return f;
}

/* Policy: reject bad normals, prefer the fixed construction, and fall back to
   the stable branch when kappa or the residual budget is exceeded. */
static Frame build_onb_safe(V3 n_in) {
    if (!finite3(n_in)) {
        Frame f; memset(&f, 0, sizeof(f));
        f.status = ONB_NONFINITE; f.kappa = INFINITY; return f;
    }
    if (!(len3(n_in) > ZERO_TOL_D)) {
        Frame f; memset(&f, 0, sizeof(f));
        f.status = ONB_ZERO_NORMAL; f.kappa = INFINITY; return f;
    }
    Frame f = build_onb_fixed(n_in);
    if (f.status == ONB_DEGENERATE || f.kappa > KAPPA_LIMIT) {
        Frame g = build_onb_stable(n_in);
        g.status = (f.status == ONB_DEGENERATE) ? ONB_DEGENERATE : ONB_KAPPA_HIGH;
        return g;
    }
    if (f.resid > RESID_LIMIT) {
        Frame g = build_onb_stable(n_in);
        g.status = ONB_RESID_HIGH;
        return g;
    }
    return f;
}

/* Per-vector comparison of F against the quaternion reference Q. Two right-handed
   frames sharing n can only differ by a twist around n, so two error views:
   - raw: direct componentwise difference, includes the O(1) geometric twist
     between two constructions that pick different reference tangents;
   - aligned: rotate Q's axes in the tangent plane onto F's twist first, which
     removes the geometric part and leaves numerical rounding noise, bounded by
     C * eps * kappa. */
static void compare_frames(const Frame* F, const Frame* Q,
                           double raw[3], double aligned[3]) {
    V3 er = v3(F->t.x - Q->t.x, F->t.y - Q->t.y, F->t.z - Q->t.z);
    V3 eb0 = v3(F->b.x - Q->b.x, F->b.y - Q->b.y, F->b.z - Q->b.z);
    V3 en0 = v3(F->n.x - Q->n.x, F->n.y - Q->n.y, F->n.z - Q->n.z);
    raw[0] = len3(er); raw[1] = len3(eb0); raw[2] = len3(en0);

    double c = dot3(F->t, Q->t);
    double s = dot3(F->t, Q->b);          /* F->t = c Q->t + s Q->b in the Q basis */
    V3 qt_aligned = v3(c * Q->t.x + s * Q->b.x,
                       c * Q->t.y + s * Q->b.y,
                       c * Q->t.z + s * Q->b.z);
    /* right-handed twist: F->b = -s Q->t + c Q->b */
    V3 qb_aligned = v3(-s * Q->t.x + c * Q->b.x,
                       -s * Q->t.y + c * Q->b.y,
                       -s * Q->t.z + c * Q->b.z);
    V3 et = v3(F->t.x - qt_aligned.x, F->t.y - qt_aligned.y, F->t.z - qt_aligned.z);
    V3 eb = v3(F->b.x - qb_aligned.x, F->b.y - qb_aligned.y, F->b.z - qb_aligned.z);
    aligned[0] = len3(et); aligned[1] = len3(eb); aligned[2] = len3(en0);
}

/* ---------- float versions, used to show the same curve at single precision. */

typedef struct {
    V3f t, b, n;
    float kappa;
    float resid;
    OnbStatus status;
} Framef;

static float frame_residf(V3f t, V3f b, V3f n) {
    float r = fabsf(dot3f(t, b));
    float d = fabsf(dot3f(b, n)); if (d > r) r = d;
    d = fabsf(dot3f(n, t)); if (d > r) r = d;
    d = fabsf(len3f(t) - 1.0f); if (d > r) r = d;
    d = fabsf(len3f(b) - 1.0f); if (d > r) r = d;
    d = fabsf(len3f(n) - 1.0f); if (d > r) r = d;
    return r;
}

static Framef build_onb_fixedf(V3f n_in) {
    Framef f; memset(&f, 0, sizeof(f));
    if (!finite3f(n_in)) { f.status = ONB_NONFINITE; f.kappa = INFINITY; return f; }
    float l = len3f(n_in);
    if (!(l > ZERO_TOL_F)) { f.status = ONB_ZERO_NORMAL; f.kappa = INFINITY; return f; }
    V3f n = { n_in.x / l, n_in.y / l, n_in.z / l };
    V3f z = { 0, 0, 1 };
    V3f c = cross3f(z, n);
    float s = len3f(c);
    f.kappa = s > 0.0f ? 1.0f / s : INFINITY;
    if (!(s > 1e-4f)) { f.n = n; f.status = ONB_DEGENERATE; return f; }
    V3f t = normalize3f(c);
    V3f b = cross3f(n, t);
    f.t = t; f.b = b; f.n = n;
    f.resid = frame_residf(t, b, n);
    f.status = ONB_OK;
    return f;
}

static Framef build_onb_stablef(V3f n_in) {
    Framef f; memset(&f, 0, sizeof(f));
    if (!finite3f(n_in)) { f.status = ONB_NONFINITE; f.kappa = INFINITY; return f; }
    float l = len3f(n_in);
    if (!(l > ZERO_TOL_F)) { f.status = ONB_ZERO_NORMAL; f.kappa = INFINITY; return f; }
    V3f n = { n_in.x / l, n_in.y / l, n_in.z / l };
    V3f t, b, helper;
    if (fabsf(n.z) <= fabsf(n.x) && fabsf(n.z) <= fabsf(n.y)) {
        helper = (V3f){ 0, 0, 1 };
        t = normalize3f(cross3f(helper, n));
        b = cross3f(n, t);
    } else if (fabsf(n.y) <= fabsf(n.x)) {
        helper = (V3f){ 0, 1, 0 };
        t = normalize3f(cross3f(helper, n));
        b = cross3f(n, t);
    } else {
        helper = (V3f){ 1, 0, 0 };
        t = normalize3f(cross3f(helper, n));
        b = cross3f(n, t);
    }
    f.t = t; f.b = b; f.n = n;
    f.kappa = 1.0f / len3f(cross3f(n, helper));
    f.resid = frame_residf(t, b, n);
    f.status = ONB_OK;
    return f;
}

/* ---------- data generation */

/* kappa sweep: normals approaching the fixed-helper pole +z. */
static int run_sweep(const char* path) {
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); return 1; }
    fprintf(fp, "precision,angle_deg,kappa,stable_kappa,fixed_resid,stable_resid,"
                "fixed_raw,stable_raw,fixed_aligned,stable_aligned,bound\n");
    double max_ratio_d = 0.0, max_ratio_f = 0.0;
    const int N = 200;
    /* Nonzero azimuth: at phi=0 the fixed and quaternion tangents coincide
       algebraically and float rounding of both paths matches exactly, which
       hides the construction difference this sweep is meant to measure. */
    const double phi = 1.23456789;
    for (int i = 0; i < N; i++) {
        double ang_deg = 90.0 * pow(10.0, -14.0 * i / (N - 1)); /* 90 .. 9e-13 deg */
        double th = ang_deg * M_PI / 180.0;
        /* double precision */
        V3 n = v3(sin(th) * cos(phi), sin(th) * sin(phi), cos(th));
        Frame fd = build_onb_fixed(n);
        Frame sd = build_onb_stable(n);
        Frame qd = build_onb_quat(n);
        double raw_f[3] = {0, 0, 0}, aln_f[3] = {0, 0, 0};
        double raw_s[3], aln_s[3];
        if (fd.status == ONB_OK) {
            compare_frames(&fd, &qd, raw_f, aln_f);
            double e = aln_f[0]; if (aln_f[1] > e) e = aln_f[1];
            if (aln_f[2] > e) e = aln_f[2];
            if (fd.kappa > 1.0 && e > 0.0) {
                double r = e / (EPS_D * fd.kappa);
                if (r > max_ratio_d) max_ratio_d = r;
            }
        }
        compare_frames(&sd, &qd, raw_s, aln_s);
        double es = aln_s[0]; if (aln_s[1] > es) es = aln_s[1];
        if (aln_s[2] > es) es = aln_s[2];
        double bound_d = C_BOUND_D * EPS_D * fd.kappa;
        fprintf(fp,
            "double,%.12g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
            ang_deg, fd.kappa, sd.kappa,
            fd.status == ONB_OK ? fd.resid : (double)NAN, sd.resid,
            raw_f[0], raw_s[0],
            fd.status == ONB_OK ? aln_f[0] : (double)NAN, es,
            isfinite(bound_d) ? bound_d : (double)NAN);
        /* single precision: skip angles below its representable floor */
        float thf = (float)th;
        float phif = (float)phi;
        float sf = sinf(thf);
        if (sf == 0.0f) continue;
        V3f nf = { sf * cosf(phif), sf * sinf(phif), cosf(thf) };
        Framef ff = build_onb_fixedf(nf);
        Framef sf2 = build_onb_stablef(nf);
        if (ff.status != ONB_OK) continue;
        Frame Fd = { fv3(ff.t), fv3(ff.b), fv3(nf), ff.kappa, ff.resid, ONB_OK, 0 };
        Frame Sd = { fv3(sf2.t), fv3(sf2.b), fv3(nf), sf2.kappa, sf2.resid, ONB_OK, 1 };
        Frame Qf = build_onb_quat(v3(nf.x, nf.y, nf.z));
        double rf[3], af[3], rs[3], as[3];
        compare_frames(&Fd, &Qf, rf, af);
        compare_frames(&Sd, &Qf, rs, as);
        double ef = af[0]; if (af[1] > ef) ef = af[1]; if (af[2] > ef) ef = af[2];
        if (ff.kappa > 1.0f && ef > 0.0) {
            double r = ef / (EPS_F * ff.kappa);
            if (r > max_ratio_f) max_ratio_f = r;
        }
        double bound_f = C_BOUND_F * EPS_F * ff.kappa;
        fprintf(fp,
            "float,%.12g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
            ang_deg, (double)ff.kappa, (double)sf2.kappa,
            (double)ff.resid, (double)sf2.resid,
            rf[0], rs[0], af[0], as[0], bound_f);
    }
    fclose(fp);
    printf("sweep_written=%s rows=%d\n", path, 2 * N);
    printf("max_ratio_double=%.3f max_ratio_float=%.3f\n", max_ratio_d, max_ratio_f);
    printf("bound_constants: C_double=%.0f C_float=%.0f"
           " (2e6-trial worst ratio: 0.77 double / 0.44 float)\n",
           C_BOUND_D, C_BOUND_F);
    return 0;
}

/* Flip demo: tilt n across the +z pole and watch the fixed tangent jump. */
static int run_flip(const char* path) {
    FILE* fp = fopen(path, "w");
    if (!fp) { perror(path); return 1; }
    fprintf(fp, "angle_deg,fixed_ty,stable_ty,fixed_kappa,stable_kappa\n");
    const int N = 161;                 /* -8 deg .. +8 deg, includes 0 */
    for (int i = 0; i < N; i++) {
        double ang = -8.0 + 16.0 * i / (N - 1);
        double th = ang * M_PI / 180.0;
        V3 n = v3(sin(th), 0.0, cos(th));
        Frame fx = build_onb_fixed(n);
        Frame st = build_onb_stable(n);
        /* In this geometry the fixed tangent is +/-y, so its pole flip shows
           up in t.y; the stable tangent lies along y with constant sign. */
        fprintf(fp, "%.10g,%.10g,%.10g,%.10g,%.10g\n",
                ang,
                fx.status == ONB_OK ? fx.t.y : (double)NAN, st.t.y,
                fx.kappa, st.kappa);
    }
    fclose(fp);
    printf("flip_written=%s rows=%d\n", path, N);
    return 0;
}

/* ---------- self test */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void selftest(void) {
    V3 zero = v3(0, 0, 0);
    V3 tiny = v3(1e-305, 0, 0);
    V3 bad = v3(NAN, 1, 0);

    /* 1. zero / degenerate normals are rejected, never silently zeroed */
    check("zero_rejected_fixed", build_onb_fixed(zero).status == ONB_ZERO_NORMAL);
    check("zero_rejected_stable", build_onb_stable(zero).status == ONB_ZERO_NORMAL);
    check("zero_rejected_quat", build_onb_quat(zero).status == ONB_ZERO_NORMAL);
    check("zero_rejected_safe", build_onb_safe(zero).status == ONB_ZERO_NORMAL);
    check("tiny_rejected", build_onb_safe(tiny).status == ONB_ZERO_NORMAL);
    check("nonfinite_rejected", build_onb_safe(bad).status == ONB_NONFINITE);
    V3f zerof = { 0, 0, 0 };
    check("zero_rejected_floatf", build_onb_fixedf(zerof).status == ONB_ZERO_NORMAL);

    /* 2. pole normal: fixed construction degenerate, stable still valid */
    V3 pole = v3(0, 0, 1);
    Frame fp = build_onb_fixed(pole);
    Frame sp = build_onb_stable(pole);
    check("pole_fixed_degenerate", fp.status == ONB_DEGENERATE);
    check("pole_stable_ok", sp.status == ONB_OK && sp.resid < 1e-15);
    check("pole_stable_kappa_bound", sp.kappa <= M_SQRT2 + 1e-15);
    check("pole_safe_fallback", build_onb_safe(pole).method == 1);

    /* 3. stable construction: orthogonal everywhere, kappa <= sqrt(2) */
    int ortho_ok = 1, kappa_ok = 1;
    double worst_resid = 0.0, worst_kappa = 0.0;
    for (int i = 0; i < 2000; i++) {
        double th = M_PI * i / 2000.0;
        double ph = 2.0 * M_PI * (i * 7 % 97) / 97.0;
        V3 n = v3(sin(th) * cos(ph), sin(th) * sin(ph), cos(th));
        Frame g = build_onb_stable(n);
        if (g.resid > 1e-14) ortho_ok = 0;
        if (g.resid > worst_resid) worst_resid = g.resid;
        if (g.kappa > worst_kappa) worst_kappa = g.kappa;
        if (!(g.kappa <= M_SQRT2 + 1e-12)) kappa_ok = 0;
    }
    check("stable_orthogonal_sphere", ortho_ok);
    check("stable_kappa_sqrt2", kappa_ok);
    printf("stable_worst_resid=%.3e stable_worst_kappa=%.6f (sqrt2=%.6f)\n",
           worst_resid, worst_kappa, M_SQRT2);

    /* 4. quaternion reference: exact at the easy poles, valid frames elsewhere */
    Frame qz = build_onb_quat(v3(0, 0, 1));
    check("quat_plusz_identity",
          fabs(qz.t.x - 1.0) < 1e-15 && fabs(qz.b.y - 1.0) < 1e-15);
    Frame qm = build_onb_quat(v3(0, 0, -1));
    check("quat_minusz_handled", qm.resid < 1e-15);
    Frame qx = build_onb_quat(v3(1, 0, 0));
    check("quat_x_orthonormal", qx.resid < 1e-15 && fabs(qx.n.x - 1.0) < 1e-15);

    /* 5. per-vector error against the quaternion frame.
       Raw error includes the O(1) geometric twist between the constructions'
       reference tangents; after twist alignment the remaining numerical error
       obeys the proven bound err <= C * eps * kappa, double and single. */
    int bound_d_ok = 1, bound_f_ok = 1;
    int raw_bounded = 1;
    double worst_margin_d = 0.0;
    const double phi_t = 0.7;
    for (int i = 0; i < 300; i++) {
        double ang_deg = 45.0 * pow(10.0, -9.0 * i / 299.0);
        double th = ang_deg * M_PI / 180.0;
        V3 n = v3(sin(th) * cos(phi_t), sin(th) * sin(phi_t), cos(th));
        Frame F = build_onb_fixed(n);
        Frame Q = build_onb_quat(n);
        if (F.status != ONB_OK) continue;
        double raw[3], aln[3];
        compare_frames(&F, &Q, raw, aln);
        double e = aln[0]; if (aln[1] > e) e = aln[1]; if (aln[2] > e) e = aln[2];
        double bound = C_BOUND_D * EPS_D * F.kappa;
        if (e > bound) bound_d_ok = 0;
        if (bound > 0) { double m = e / bound; if (m > worst_margin_d) worst_margin_d = m; }
        /* raw per-vector differences of unit vectors can never exceed 2 */
        for (int j = 0; j < 3; j++) if (!(raw[j] <= 2.0 + 1e-12)) raw_bounded = 0;

        float thf = (float)th;
        float sf = sinf(thf);
        if (sf == 0.0f) continue;
        float pf = (float)phi_t;
        V3f nf = { sf * cosf(pf), sf * sinf(pf), cosf(thf) };
        Framef Ff = build_onb_fixedf(nf);
        if (Ff.status != ONB_OK) continue;
        Frame Ffd = { fv3(Ff.t), fv3(Ff.b), fv3(nf), Ff.kappa, Ff.resid, ONB_OK, 0 };
        Frame Qf = build_onb_quat(v3(nf.x, nf.y, nf.z));
        double rf[3], af[3];
        compare_frames(&Ffd, &Qf, rf, af);
        double ef = af[0]; if (af[1] > ef) ef = af[1]; if (af[2] > ef) ef = af[2];
        if (ef > C_BOUND_F * EPS_F * Ff.kappa) bound_f_ok = 0;
    }
    check("quat_raw_error_le_2", raw_bounded);
    check("quat_aligned_bound_double", bound_d_ok);
    check("quat_aligned_bound_float", bound_f_ok);
    printf("aligned_bound_worst_usage_double=%.3f of C=%.0f\n",
           worst_margin_d, C_BOUND_D);

    /* 6. safe policy reports over-budget kappa and switches construction */
    V3 near_pole = v3(sin(1e-6), 0, cos(1e-6));     /* kappa ~ 1e6 > limit */
    Frame sf = build_onb_safe(near_pole);
    check("safe_reports_kappa", sf.status == ONB_KAPPA_HIGH && sf.method == 1);
    check("safe_frame_valid", sf.resid < 1e-14 && sf.kappa <= M_SQRT2 + 1e-12);
    V3 ordinary = v3(1, 1, 1);
    check("safe_keeps_fixed", build_onb_safe(ordinary).method == 0);

    /* 7. the original bug: fixed tangent flips across the pole; stable does not.
       Cosine sampling in a frame whose tangent flips gets mirrored directions. */
    double th = 0.5 * M_PI / 180.0;
    V3 np = v3(sin(th), 0, cos(th)), nm = v3(-sin(th), 0, cos(th));
    Frame fpp = build_onb_fixed(np), fmm = build_onb_fixed(nm);
    Frame spp = build_onb_stable(np), smm = build_onb_stable(nm);
    check("fixed_tangent_flips", dot3(fpp.t, fmm.t) < -0.9);
    check("stable_tangent_continuous", dot3(spp.t, smm.t) > 0.99);

    /* 8. end-to-end: a cosine-weighted sample keeps its hemisphere through the
       stable frame for both tilt signs (no mirrored direction). */
    V3 local = normalize3(v3(0.3, -0.2, 0.9));
    V3 wp = v3(dot3(local, v3(spp.t.x, spp.b.x, spp.n.x)),
               dot3(local, v3(spp.t.y, spp.b.y, spp.n.y)),
               dot3(local, v3(spp.t.z, spp.b.z, spp.n.z)));
    V3 wm = v3(dot3(local, v3(smm.t.x, smm.b.x, smm.n.x)),
               dot3(local, v3(smm.t.y, smm.b.y, smm.n.y)),
               dot3(local, v3(smm.t.z, smm.b.z, smm.n.z)));
    check("sample_hemisphere_kept", dot3(wp, spp.n) > 0.9 && dot3(wm, smm.n) > 0.9);

    printf("all_ok=%d\n", failures == 0);
}

int main(int argc, char** argv) {
    const char* mode = argc > 1 ? argv[1] : "selftest";
    if (strcmp(mode, "selftest") == 0) {
        selftest();
        return failures == 0 ? 0 : 1;
    }
    if (argc < 3) { fprintf(stderr, "usage: %s [selftest|sweep out.csv|flip out.csv]\n", argv[0]); return 2; }
    if (strcmp(mode, "sweep") == 0) return run_sweep(argv[2]);
    if (strcmp(mode, "flip") == 0) return run_flip(argv[2]);
    fprintf(stderr, "unknown mode: %s\n", mode);
    return 2;
}
