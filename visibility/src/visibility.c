/* visibility: piecewise-linear visibility for a spherical occluder in front of a
   segment area light.

   Geometry (all in world space, any dimension-free units):

       P          receiver (shading point)
       L0, L1     endpoints of a 1-D area light, parameterised by s in [0,1]
       C, r       occluder centre and radius

   For a light sample O(s) = L0 + s*(L1-L0) the shadow ray is P + t*(O(s)-P),
   t in (0,1).  It is blocked iff the ray meets the sphere.  Writing
   c = C-P, v(s) = O(s)-P, the blocked set is exactly

       phi(s) = |c x v(s)|^2 - r^2 |v(s)|^2 <= 0

   and phi is a quadratic in s, so the *exact* blocked measure of any interval
   is computable in closed form.  The module approximates the blocked set per
   segment by the chord through the two endpoint margins (a piecewise linear
   approximation), and certifies the result with a rigorous bound derived from
   the curvature of phi.

   Per segment, with h = 1/K and e = M2*h^2/8 (M2 an upper bound on |phi''|),
   the chord error at each band edge is at most 2*e / |slope|.  The certificate
   for a segment is the length of the sub-interval where the chord lies within
   [-e, e] -- the only place the piecewise-linear blocked set can be wrong.

   Properties the self-test pins down:
     * the certificate is an upper bound on |piecewise - exact|,
     * the certificate is an upper bound on |piecewise - per-point shadow rays|,
     * piecewise-accuracy is second order (e * K^2 stays bounded),
     * any segment whose visible share is exactly one is never integrated,
     * K <= 0 is rejected (returns 0 and an error code). */

#include <stdio.h>
#include <string.h>
#include <math.h>

/* ------------------------------------------------------------------ errors */

enum {
    VIS_OK = 0,
    VIS_ERR_K = -1,        /* segment count <= 0 */
    VIS_ERR_DEGENERATE = -2, /* zero-length light, non-positive radius, P == C */
    VIS_ERR_GEOMETRY = -3  /* occluder not between P and the light samples */
};

/* ---------------------------------------------------------------- vectors */

typedef struct { double x, y, z; } V3;

static V3 v3(double x, double y, double z) { V3 v = { x, y, z }; return v; }
static V3 vsub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static V3 vadd(V3 a, V3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static V3 vmul(V3 a, double s) { return v3(a.x * s, a.y * s, a.z * s); }
static double vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vcross(V3 a, V3 b) {
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

/* ---------------------------------------------------------------- problem */

typedef struct {
    V3 P;          /* receiver */
    V3 L0, L1;     /* light endpoints */
    V3 C;          /* occluder centre */
    double r;      /* occluder radius */

    double U2;     /* |L1-L0|^2 */
    double r2;     /* r^2 */
    double q2, q1, q0; /* phi(s) = q2*s^2 + q1*s + q0 */
    double M2;     /* bound on |phi''| = 2*max(|q2|, ...) -- see setup() */

    int error;
} Vis;

/* evaluate phi at s */
static double vis_phi(const Vis* p, double s) {
    return (p->q2 * s + p->q1) * s + p->q0;
}

/* measure of chord x0 + (x1-x0)*t/h <= 0 for t in [0,h] */
static double vis_under(const double x0, const double x1, const double h) {
    double m = (x1 - x0) / h;
    if (m == 0.0) return x0 <= 0.0 ? h : 0.0;
    double t = -x0 / m; /* crossing */
    if (m > 0.0) { /* chord increasing: below zero on [0, t] */
        if (t <= 0.0) return 0.0;
        if (t >= h) return h;
        return t;
    }
    /* chord decreasing: below zero on [t, h] */
    if (t <= 0.0) return h;
    if (t >= h) return 0.0;
    return h - t;
}

/* measure of {t in [0,h] : |chord(t)| <= e} */
static double vis_band(const double x0, const double x1, const double h, const double e) {
    double m = (x1 - x0) / h;
    if (m == 0.0) return fabs(x0) <= e ? h : 0.0;
    double a = (-e - x0) / m;
    double b = (e - x0) / m;
    if (a > b) { double tmp = a; a = b; b = tmp; }
    if (a < 0.0) a = 0.0;
    if (b > h) b = h;
    return b > a ? b - a : 0.0;
}

/* ---------------------------------------------------------------- setup */

/* Initialise and validate.  Returns VIS_OK or a negative VIS_ERR_* code and
   leaves *out->error set.  The quadratic coefficients are always filled in
   when the input is non-degenerate, so callers can still inspect phi. */
static int vis_setup(Vis* out, V3 P, V3 L0, V3 L1, V3 C, double r) {
    memset(out, 0, sizeof(*out));
    out->P = P; out->L0 = L0; out->L1 = L1; out->C = C; out->r = r;

    V3 u = vsub(L1, L0);
    V3 c = vsub(C, P);
    out->U2 = vdot(u, u);
    out->r2 = r * r;

    if (!(out->U2 > 1e-24) || !(out->r2 > 0.0)) {
        out->error = VIS_ERR_DEGENERATE;
        return out->error;
    }
    if (!(vdot(c, c) > out->r2)) { /* P inside or on the occluder */
        out->error = VIS_ERR_DEGENERATE;
        return out->error;
    }

    /* phi(s) = |c x v|^2 - r^2 |v|^2 with v = v0 + s*u, v0 = L0-P */
    V3 v0 = vsub(L0, P);
    V3 cu = vcross(c, u);
    V3 cv0 = vcross(c, v0);
    out->q2 = vdot(cu, cu) - out->r2 * out->U2;
    out->q1 = 2.0 * (vdot(cv0, cu) - out->r2 * vdot(v0, u));
    out->q0 = vdot(cv0, cv0) - out->r2 * vdot(v0, v0);

    /* |phi''| = 2*|q2|; add slack for the linear term's round-off. */
    out->M2 = 2.0 * fabs(out->q2) + 1e-12 * (fabs(out->q1) + fabs(out->q0) + 1.0);

    /* The ray parameter at closest approach to C is t*(s) = (c.v)/(v.v).
       The ray must reach the occluder strictly before the light sample and
       the occluder must sit in front of P, for every s in [0,1]; otherwise a
       single occluder can shadow only part of a sample's segment and the
       quadratic cone test stops matching a per-point ray test.  Also the
       light segment itself must stay outside the sphere. */

    /* (1) t*(s) < 1  <=>  c.v - |v|^2 < 0 ; a downward parabola (leading
       coefficient -U2), so check both endpoints. */
    double g_end0 = vdot(c, v0) - vdot(v0, v0);
    V3 v1 = vsub(L1, P);
    double g_end1 = vdot(c, v1) - vdot(v1, v1);
    if (g_end0 >= 0.0 || g_end1 >= 0.0) {
        out->error = VIS_ERR_GEOMETRY;
        return out->error;
    }

    /* (2) t*(s) > 0  <=>  c.v > 0 ; linear in s, so check both endpoints. */
    if (vdot(c, v0) <= 0.0 || vdot(c, v1) <= 0.0) {
        out->error = VIS_ERR_GEOMETRY;
        return out->error;
    }

    /* (3) min |O(s)-C|^2 > r^2.  |v0-c + s*u|^2 is a convex parabola. */
    V3 w = vsub(v0, c);
    double h0 = vdot(w, w) - out->r2;
    double h1 = vdot(vadd(w, u), vadd(w, u)) - out->r2;
    double hmin = h0 < h1 ? h0 : h1;
    if (out->U2 > 0.0) {
        double sv = -vdot(w, u) / out->U2;
        if (sv > 0.0 && sv < 1.0) {
            double hv = vdot(vadd(w, vmul(u, sv)), vadd(w, vmul(u, sv))) - out->r2;
            if (hv < hmin) hmin = hv;
        }
    }
    if (!(hmin > 0.0)) {
        out->error = VIS_ERR_GEOMETRY;
        return out->error;
    }

    out->error = VIS_OK;
    return VIS_OK;
}

/* --------------------------------------------------------------- result */

typedef struct {
    double visible;   /* piecewise-linear estimate of the visible fraction */
    double bound;     /* certified |estimate - exact| upper bound */
    int    n_full;    /* segments with visible share exactly 1 (skipped) */
    int    n_blocked; /* segments with visible share exactly 0 */
    int    n_partial; /* segments with a strictly fractional share */
    int    n_integrated; /* n_blocked + n_partial: segments that entered the sum */
} VisResult;

/* Write the per-segment visible share into shares[0..K-1] (band-relative:
   every segment has the same length h = 1/K).  Segments reported as full are
   written as exactly 1.0 -- they are skipped by the integrator.  Returns
   VIS_OK or a negative VIS_ERR_* code. */
static int vis_shares(const Vis* p, int K, double* shares) {
    if (p->error != VIS_OK) return p->error;
    if (K <= 0) return VIS_ERR_K;
    const double h = 1.0 / (double)K;
    double prev = vis_phi(p, 0.0);
    for (int i = 0; i < K; i++) {
        double next = vis_phi(p, (double)(i + 1) * h);
        double blocked = vis_under(prev, next, h);
        double share = 1.0;
        if (blocked >= h) share = 0.0;
        else if (blocked > 0.0) share = 1.0 - blocked / h;
        shares[i] = share;
        prev = next;
    }
    return VIS_OK;
}

/* Evaluate with K segments.  Returns VIS_OK, or VIS_ERR_K when K <= 0 (out is
   zeroed in that case).  Segments whose visible share is exactly one are not
   integrated: they keep share 1 and contribute nothing to the occlusion sum. */
static int vis_eval(const Vis* p, int K, VisResult* out) {
    memset(out, 0, sizeof(*out));
    if (p->error != VIS_OK) return p->error;
    if (K <= 0) return VIS_ERR_K;

    const double h = 1.0 / (double)K;
    const double e = p->M2 * h * h / 8.0;

    double occluded = 0.0;
    double bound = 0.0;
    double prev = vis_phi(p, 0.0);

    for (int i = 0; i < K; i++) {
        double next = vis_phi(p, (double)(i + 1) * h);

        double blocked = vis_under(prev, next, h);
        if (blocked <= 0.0) {
            /* fully visible: share one, never integrated */
            out->n_full++;
        } else if (blocked >= h) {
            out->n_blocked++;
            out->n_integrated++;
            occluded += h;
        } else {
            out->n_partial++;
            out->n_integrated++;
            occluded += blocked;
        }

        /* certificate: every point of this segment whose chord lies within
           e of zero may be misclassified.  Computed for every segment, even
           the skipped ones, because the bound must cover the whole domain. */
        bound += vis_band(prev, next, h, e);

        prev = next;
    }

    out->visible = 1.0 - occluded;
    out->bound = bound;
    return VIS_OK;
}

/* Exact visible fraction via the closed-form roots of phi (reference only). */
static double vis_exact(const Vis* p) {
    double q2 = p->q2, q1 = p->q1, q0 = p->q0;
    double roots[2];
    int nr = 0;
    double eps = 1e-14 * (fabs(q2) + fabs(q1) + fabs(q0) + 1.0);
    if (fabs(q2) <= eps) {
        if (fabs(q1) > eps) roots[nr++] = -q0 / q1;
    } else {
        double disc = q1 * q1 - 4.0 * q2 * q0;
        if (disc >= 0.0) {
            double sq = sqrt(disc);
            roots[nr++] = (-q1 - sq) / (2.0 * q2);
            roots[nr++] = (-q1 + sq) / (2.0 * q2);
        }
    }
    for (int i = 0; i < nr; i++) {
        if (roots[i] < 0.0) roots[i] = 0.0;
        if (roots[i] > 1.0) roots[i] = 1.0;
    }
    if (nr == 2 && roots[0] > roots[1]) { double t = roots[0]; roots[0] = roots[1]; roots[1] = t; }

    double pts[4];
    int np = 0;
    pts[np++] = 0.0;
    if (nr > 0) pts[np++] = roots[0];
    if (nr > 1) pts[np++] = roots[1];
    pts[np++] = 1.0;

    double blocked = 0.0;
    for (int i = 0; i + 1 < np; i++) {
        double a = pts[i], b = pts[i + 1];
        if (b - a < 1e-15) continue;
        if (vis_phi(p, 0.5 * (a + b)) <= 0.0) blocked += b - a;
    }
    return 1.0 - blocked;
}

/* Per-point shadow rays: M evenly spaced samples, each tested against the
   sphere with a true ray/sphere intersection.  Reference only. */
static double vis_pointwise(const Vis* p, int M) {
    if (M <= 0) return 0.0;
    int blocked = 0;
    for (int i = 0; i < M; i++) {
        double s = ((double)i + 0.5) / (double)M;
        V3 O = vadd(p->L0, vmul(vsub(p->L1, p->L0), s));
        V3 v = vsub(O, p->P);
        V3 f = vsub(p->P, p->C);
        double A = vdot(v, v);
        double B = 2.0 * vdot(f, v);
        double Cc = vdot(f, f) - p->r2;
        double disc = B * B - 4.0 * A * Cc;
        if (disc < 0.0) continue;
        double sq = sqrt(disc);
        double t0 = (-B - sq) / (2.0 * A);
        double t1 = (-B + sq) / (2.0 * A);
        if (t1 > 0.0 && t0 < 1.0) blocked++;
    }
    return 1.0 - (double)blocked / (double)M;
}

/* ----------------------------------------------------------------- tests */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* slope of log2(err) versus log2(K) over the given K values */
static double log2_slope(const int* Ks, const double* es, int n) {
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < n; i++) {
        double x = log2((double)Ks[i]);
        double y = log2(es[i]);
        sx += x; sy += y; sxx += x * x; sxy += x * y;
    }
    double den = (double)n * sxx - sx * sx;
    if (fabs(den) < 1e-30) return 0.0;
    return ((double)n * sxy - sx * sy) / den;
}

static Vis make(V3 P, V3 L0, V3 L1, V3 C, double r, int* err) {
    Vis v;
    *err = vis_setup(&v, P, L0, L1, C, r);
    return v;
}

/* one complete per-scene battery: certificate validity, pointwise agreement,
   shares, and second-order accuracy */
static void scene_checks(const char* tag, V3 P, V3 L0, V3 L1, V3 C, double r) {
    int err;
    Vis v = make(P, L0, L1, C, r, &err);
    char name[128];
    snprintf(name, sizeof(name), "%s_setup", tag);
    check(name, err == VIS_OK);
    if (err != VIS_OK) {
        printf("%s_error=%d\n", tag, err);
        return;
    }

    double exact = vis_exact(&v);
    double point = vis_pointwise(&v, 1 << 17);

    static const int Ks[7] = { 16, 32, 64, 128, 256, 512, 1024 };
    double err_exact[7], err_point[7];

    int cert_ok = 1, point_ok = 1, share_ok = 1, class_ok = 1;
    for (int i = 0; i < 7; i++) {
        VisResult res;
        int rc = vis_eval(&v, Ks[i], &res);
        if (rc != VIS_OK) { cert_ok = 0; break; }

        err_exact[i] = fabs(res.visible - exact);
        err_point[i] = fabs(res.visible - point);

        if (!(err_exact[i] <= res.bound + 1e-12)) cert_ok = 0;
        if (!(err_point[i] <= res.bound + 2.0 / (double)(1 << 17) + 1e-12)) point_ok = 0;
        if (res.n_full + res.n_blocked + res.n_partial != Ks[i]) class_ok = 0;
        if (res.n_integrated != res.n_blocked + res.n_partial) class_ok = 0;

        /* "share one" segments must not enter the integral: the visible
           fraction recomputed from the per-segment shares must match, and all
           fully visible segments must carry share exactly 1. */
        double shares[1024];
        if (vis_shares(&v, Ks[i], shares) != VIS_OK) { share_ok = 0; continue; }
        double sum_share = 0.0;
        double h = 1.0 / (double)Ks[i];
        for (int j = 0; j < Ks[i]; j++) {
            sum_share += h * shares[j];
            if (shares[j] >= 1.0 && shares[j] != 1.0) share_ok = 0;
        }
        if (!(fabs(sum_share - res.visible) < 1e-12)) share_ok = 0;
        if (res.n_full + res.n_integrated != Ks[i]) share_ok = 0;
    }

    snprintf(name, sizeof(name), "%s_cert", tag);   check(name, cert_ok);
    snprintf(name, sizeof(name), "%s_point", tag);  check(name, point_ok);
    snprintf(name, sizeof(name), "%s_class", tag);  check(name, class_ok);
    snprintf(name, sizeof(name), "%s_share", tag);  check(name, share_ok);

    /* second order: e*K^2 must stay bounded and roughly flat across octaves.
       The error at a fixed K is phase dependent (each chord straddles a band
       edge at a varying offset), so one fixed K can sit several times above
       or below its K^-2 envelope.  Average e*K^2 over each octave [K0,2K0):
       if the scheme were first order these octave means would double every
       octave, so requiring their spread to be modest pins down order two. */
    double oct[5];
    for (int j = 0; j < 5; j++) {
        int lo = 64 << j;
        double acc = 0.0;
        for (int K = lo; K < 2 * lo; K++) {
            VisResult res;
            vis_eval(&v, K, &res);
            acc += fabs(res.visible - exact) * (double)K * (double)K;
        }
        oct[j] = acc / (double)lo;
    }
    double omin = oct[0], omax = oct[0];
    for (int j = 1; j < 5; j++) {
        if (oct[j] < omin) omin = oct[j];
        if (oct[j] > omax) omax = oct[j];
    }
    snprintf(name, sizeof(name), "%s_order2", tag);
    check(name, omax < 2.5 && omax / omin < 2.0);
    printf("%s_envelope=%.4f..%.4f spread=%.2f\n", tag, omin, omax, omax / omin);

    /* a single-knot fit must still be at least first order */
    double slope = log2_slope(Ks + 2, err_exact + 2, 5);
    snprintf(name, sizeof(name), "%s_fit", tag);
    check(name, fabs(slope) > 1.0);
    printf("%s_order_value=%.3f\n", tag, fabs(slope));

    /* the bound must also shrink at second order */
    double b64, b1024;
    {
        VisResult a, b;
        vis_eval(&v, 64, &a);
        vis_eval(&v, 1024, &b);
        b64 = a.bound; b1024 = b.bound;
    }
    double border = (b64 > 0.0) ? log2(b64 / b1024) / log2(1024.0 / 64.0) : 2.0;
    snprintf(name, sizeof(name), "%s_bound_order", tag);
    check(name, border > 1.9);
}

int main(void) {
    /* --- rejection paths ------------------------------------------------ */
    {
        VisResult res;
        int err;
        Vis v = make(v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(0, 0, -1.5), 0.45, &err);
        check("setup_basic", err == VIS_OK);

        int rc = vis_eval(&v, 0, &res);
        check("zero_segments_rejected", rc == VIS_ERR_K);
        check("zero_segments_result_zeroed", res.visible == 0.0 && res.bound == 0.0);
        rc = vis_eval(&v, -4, &res);
        check("negative_segments_rejected", rc == VIS_ERR_K);
    }
    {
        int err;
        make(v3(0, 0, -3), v3(0, 0, 0), v3(0, 0, 0), v3(0, 0, -1), 0.4, &err);
        check("degenerate_light", err == VIS_ERR_DEGENERATE);
        make(v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(0, 0, -1), 0.0, &err);
        check("degenerate_radius", err == VIS_ERR_DEGENERATE);
        /* receiver inside the occluder */
        make(v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(0, 0, -3), 0.4, &err);
        check("receiver_inside", err == VIS_ERR_DEGENERATE);
    }
    {
        int err;
        /* occluder behind the receiver */
        make(v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(0, 0, -5), 0.4, &err);
        check("occluder_behind", err == VIS_ERR_GEOMETRY);
        /* light segment pokes into the sphere */
        make(v3(0, 0, -3), v3(-1, 0, -1.5), v3(1, 0, -1.5), v3(0, 0, -1.5), 0.6, &err);
        check("light_inside_occluder", err == VIS_ERR_GEOMETRY);
    }

    /* --- fully lit: every share is one, nothing is integrated ----------- */
    {
        int err;
        Vis v = make(v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(1.9, 0, -1.5), 0.25, &err);
        check("lit_setup", err == VIS_OK);
        VisResult res;
        int ok = 1;
        for (int K = 1; K <= 512; K *= 2) {
            if (vis_eval(&v, K, &res) != VIS_OK) { ok = 0; break; }
            if (res.n_full != K || res.n_integrated != 0) ok = 0;
            if (res.visible != 1.0 || res.bound != 0.0) ok = 0;
        }
        check("lit_skips_every_segment", ok);
    }

    /* --- fully shadowed ------------------------------------------------ */
    {
        int err;
        Vis v = make(v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(0, 0, -1.5), 1.0, &err);
        check("blocked_setup", err == VIS_OK);
        VisResult res;
        int ok = 1;
        for (int K = 1; K <= 512; K *= 2) {
            if (vis_eval(&v, K, &res) != VIS_OK) { ok = 0; break; }
            if (res.n_full != 0 || res.n_partial != 0 || res.n_blocked != K) ok = 0;
            if (fabs(res.visible) > 1e-15) ok = 0;
        }
        check("blocked_all_integrated", ok);
    }

    /* --- half-covered, offset light: the motivating case ---------------- */
    scene_checks("half", v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(0.05, 0, -1.5), 0.45);

    /* --- twist the light and occluder off the axis ---------------------- */
    scene_checks("skew", v3(0.2, 0.1, -3.2), v3(-1.2, -0.1, 0.1), v3(1.1, 0.2, -0.1),
                 v3(-0.1, 0.05, -1.6), 0.35);

    /* --- light tilted in depth, occluder grazing the shadow boundary ---- */
    scene_checks("grazing", v3(-0.3, 0, -2.6), v3(-1.5, 0.2, 0.2), v3(1.4, -0.2, 0.0),
                 v3(0.25, 0.1, -1.4), 0.5);

    /* --- taller receiver, light behind a raised occluder ---------------- */
    scene_checks("raised", v3(0, 0.3, -3.5), v3(-1.4, 0, 0), v3(1.4, 0, 0),
                 v3(0.15, 0.15, -1.5), 0.6);

    /* --- smoothness in the occluder position (no flicker) --------------- */
    {
        double prev_pw = -1.0, prev_pt = -1.0;
        double worst_pw = 0.0, worst_pt = 0.0, worst_gap = 0.0;
        int ok = 1;
        for (int j = 0; j <= 120; j++) {
            double cx = -1.2 + 2.4 * (double)j / 120.0;
            int err;
            Vis v = make(v3(0, 0, -3), v3(-1, 0, 0), v3(1, 0, 0), v3(cx, 0, -1.5), 0.45, &err);
            if (err != VIS_OK) { ok = 0; break; }
            VisResult res;
            if (vis_eval(&v, 256, &res) != VIS_OK) { ok = 0; break; }
            double pt = vis_pointwise(&v, 4096);
            if (prev_pw >= 0.0) {
                double a = fabs(res.visible - prev_pw);
                double b = fabs(pt - prev_pt);
                if (a > worst_pw) worst_pw = a;
                if (b > worst_pt) worst_pt = b;
            }
            /* pointwise error is bounded by 1/M per sample */
            if (!(fabs(res.visible - pt) <= res.bound + 1.0 / 4096.0 + 1e-12)) ok = 0;
            if (fabs(res.visible - pt) > worst_gap) worst_gap = fabs(res.visible - pt);
            prev_pw = res.visible;
            prev_pt = pt;
        }
        /* translate speed is 90% of the light width per full sweep, so a
           per-step jump of the piecewise estimate must not exceed the exact
           change by more than the per-point sampling noise. */
        check("sweep_smooth", ok && worst_pw <= worst_pt + 2.0 / 4096.0);
        printf("sweep_max_step=%.4f pointwise_max_step=%.4f worst_gap=%.4f\n",
               worst_pw, worst_pt, worst_gap);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
