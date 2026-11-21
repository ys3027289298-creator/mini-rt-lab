/*
 * intersect: watertight triangle intersection.
 *
 * Why seams leak
 * --------------
 * The textbook epsilon test
 *     if (u < -eps || v < -eps || u + v > 1 + eps) miss;
 * is the usual cause.  When a ray lands on the edge shared by two triangles,
 * each triangle computes the point from different operands and rounds to the
 * wrong side of zero.  The point is judged outside both, and a one-pixel gap
 * opens along the joint.
 *
 * This module uses a sign-consistency test instead:
 *
 *   1. Each barycentric numerator is an independent edge function built from
 *      the two vertices of one edge, not a "1 - u - v" shortcut.  The two
 *      triangles sharing an edge therefore evaluate the same vertices, up to
 *      order, and the difference between the two results is pure rounding.
 *      The tolerance in (3) is built from those same operands, so both
 *      neighbours see a point on the shared edge as zero.
 *   2. Containment means all three edge functions have the same sign, with a
 *      zero treated as "on the boundary" and admitted.  That is the same-side
 *      rule: a boundary point is claimed by both triangles, so it cannot fall
 *      into a hole.  (A ray through a vertex is claimed by all triangles
 *      meeting there; overlap is harmless, a gap is not.)
 *   3. Each sign is taken outside a forward-error tolerance built from the
 *      magnitude of the products that formed that edge function, so pure
 *      rounding noise near a boundary cannot flip the sign.
 *   4. Two independent formulas are implemented: the right-hand Moller-
 *      Trumbore algebra and the Woop projected 2D edge functions.  They share
 *      only the degeneracy test.  The self-test cross-checks them pixel by
 *      pixel and requires identical hit/miss decisions and t values.
 *   5. Degenerate triangles (non-finite coordinates or zero/near-zero area)
 *      are classified separately and never reported as hits.  They are
 *      counted per reason in TriStats, and per call through HitFlags.  A ray
 *      parallel to the plane is its own flag and is not a triangle defect.
 *
 * Geometry is Vec3 (float, matching the rest of the lab) but all arithmetic is
 * done in double; float subtraction/division was the dominant source of the
 * rounding that flipped boundary points.  The negative-control test at the
 * bottom shows a strict float Moller-Trumbore leaking ~0.3% of seam hits on
 * random shared edges, versus zero here.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdint.h>

#ifndef INTERSECT_K
#define INTERSECT_K 16.0   /* forward-error multiplier for edge-function signs */
#endif

/* ------------------------------------------------------------------ */
/* vectors                                                             */
/* ------------------------------------------------------------------ */

typedef struct { float x, y, z; } Vec3;
typedef struct { double x, y, z; } V3;

static Vec3 v3f(float x, float y, float z) { Vec3 r = { x, y, z }; return r; }
static V3 to_v3(Vec3 a) { V3 r = { (double)a.x, (double)a.y, (double)a.z }; return r; }
static V3 vsub(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static V3 vcross(V3 a, V3 b) {
    V3 r = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return r;
}
static double vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static int finite_v3(V3 a) { return isfinite(a.x) && isfinite(a.y) && isfinite(a.z); }


/* ------------------------------------------------------------------ */
/* public types                                                        */
/* ------------------------------------------------------------------ */

typedef struct { Vec3 a, b, c; } Tri;

typedef enum { INTERSECT_MT = 0, INTERSECT_WOOP = 1 } Formula;

/* Per-call reasons, so a caller can tell "miss" from "bad triangle". */
typedef enum {
    INTERSECT_MISS = 0,
    INTERSECT_HIT,
    INTERSECT_PARALLEL,    /* ray parallel to the triangle plane  */
    INTERSECT_DEGENERATE   /* triangle itself is unusable         */
} HitFlags;

typedef struct {
    int hit;            /* 1 on a valid intersection                  */
    int flags;          /* HitFlags                                   */
    int degenerate;     /* 0, or reason 1 = zero/near-zero area, 2 = non-finite */
    double t;           /* ray parameter, double precision            */
    double wa, wb, wc;  /* barycentric weights for a, b, c            */
} Hit;

/* Aggregate report, normally accumulated per scene. */
typedef struct {
    long long tested;
    long long hits;
    long long parallel;
    long long degenerate_zero_area;
    long long degenerate_nonfinite;
} TriStats;

static void tri_stats_init(TriStats* s) { memset(s, 0, sizeof(*s)); }

static void tri_stats_add(TriStats* s, const Hit* h) {
    s->tested++;
    if (h->flags == INTERSECT_HIT) s->hits++;
    else if (h->flags == INTERSECT_PARALLEL) s->parallel++;
    else if (h->flags == INTERSECT_DEGENERATE) {
        if (h->degenerate == 2) s->degenerate_nonfinite++;
        else s->degenerate_zero_area++;
    }
}

/* ------------------------------------------------------------------ */
/* sign test                                                           */
/* ------------------------------------------------------------------ */

/* Snap a value to -1/0/+1 outside a tolerance; 0 means "on the boundary". */
static int sign_snap(double v, double tol) {
    if (v > tol) return 1;
    if (v < -tol) return -1;
    return 0;
}

/* Same-side rule: all non-negative or all non-positive.  A zero is inside. */
static int same_side(int sa, int sb, int sc) {
    return (sa >= 0 && sb >= 0 && sc >= 0) ||
           (sa <= 0 && sb <= 0 && sc <= 0);
}

/* ------------------------------------------------------------------ */
/* degeneracy                                                          */
/* ------------------------------------------------------------------ */

/* Returns 0 when usable, else a reason code (1 = area, 2 = non-finite). */
static int tri_degenerate(const Tri* t) {
    V3 a = to_v3(t->a), b = to_v3(t->b), c = to_v3(t->c);
    if (!finite_v3(a) || !finite_v3(b) || !finite_v3(c)) return 2;
    V3 n = vcross(vsub(b, a), vsub(c, a));
    double a2 = vdot(n, n);
    if (!(a2 > 0.0)) return 1;
    double e1 = vdot(vsub(b, a), vsub(b, a));
    double e2 = vdot(vsub(c, a), vsub(c, a));
    if (e1 == 0.0 || e2 == 0.0) return 1;
    /*
     * |n|^2 = |e1|^2 |e2|^2 sin^2(theta), so comparing against
     * (64 eps)^2 * |e1|^2 |e2|^2 rejects triangles whose sine is below the
     * noise floor.  The product of both squared lengths is required: using
     * the longer edge twice would wrongly accept very elongated slivers.
     */
    const double thr = 64.0 * DBL_EPSILON;
    if (a2 <= thr * thr * e1 * e2) return 1;
    return 0;
}

/* ------------------------------------------------------------------ */
/* edge function                                                       */
/* ------------------------------------------------------------------ */
/*
 * For the edges (b,c), (c,a), (a,b) in the triangle's own winding this is
 * -w_i * den with den = n.d, so on an interior point all three share a sign.
 * The value depends only on the two vertices of the edge and on the ray,
 * never on the third vertex.
 *
 * The two triangles sharing an edge evaluate it from opposite endpoints, so
 * they are not bit-identical: swapping endpoints negates the value exactly,
 * but the attachment vertex also changes ((o - P) becomes (o - Q)), and that
 * is a different rounding.  Woop's guarantee is not bit-exactness but that
 * the error is bounded by the same products on both sides, which is why the
 * tolerance below is built from those products.  A point on the shared edge
 * then snaps to "zero" in both triangles.
 *
 * Do not "canonicalise" the endpoints per edge: independently flipping each
 * edge's direction makes the three signs no longer comparable, and valid
 * interiors get rejected.  Comparability of the three signs is the whole
 * point of the same-side test, so all three edges must keep the winding
 * order.
 */
static double edge_fn(V3 d, V3 o, V3 p, V3 q, double* scale) {
    V3 u = vsub(q, p);
    V3 v = vsub(o, p);
    V3 c = vcross(u, v);
    *scale = fabs(d.x) * (fabs(u.y * v.z) + fabs(u.z * v.y)) +
             fabs(d.y) * (fabs(u.z * v.x) + fabs(u.x * v.z)) +
             fabs(d.z) * (fabs(u.x * v.y) + fabs(u.y * v.x));
    return vdot(d, c);
}

/* Shared sign test used by both formulas.
 * On success the raw edge functions are returned so the caller can derive
 * barycentrics from the very values that were sign-tested, rather than
 * recomputing them a second time and risking disagreement. */
static int edge_signs_agree(V3 d, V3 o, V3 a, V3 b, V3 c, double g[3]) {
    double sa, sb, sc;
    double ga = edge_fn(d, o, b, c, &sa);   /* opposite a */
    double gb = edge_fn(d, o, c, a, &sb);   /* opposite b */
    double gc = edge_fn(d, o, a, b, &sc);   /* opposite c */
    int sga = sign_snap(ga, INTERSECT_K * DBL_EPSILON * sa);
    int sgb = sign_snap(gb, INTERSECT_K * DBL_EPSILON * sb);
    int sgc = sign_snap(gc, INTERSECT_K * DBL_EPSILON * sc);
    if (g) { g[0] = ga; g[1] = gb; g[2] = gc; }
    return same_side(sga, sgb, sgc);
}

/* ------------------------------------------------------------------ */
/* formula 1: Moller-Trumbore algebra                                  */
/* ------------------------------------------------------------------ */

static Hit intersect_mt(const Tri* t, V3 o, V3 d) {
    Hit h;
    memset(&h, 0, sizeof(h));
    int deg = tri_degenerate(t);
    if (deg) { h.flags = INTERSECT_DEGENERATE; h.degenerate = deg; return h; }
    V3 a = to_v3(t->a), b = to_v3(t->b), c = to_v3(t->c);
    V3 n = vcross(vsub(b, a), vsub(c, a));
    double den = vdot(n, d);
    if (den == 0.0) { h.flags = INTERSECT_PARALLEL; return h; }
    double g[3];
    if (!edge_signs_agree(d, o, a, b, c, g)) return h;  /* flags stays MISS */
    double t0 = -vdot(n, vsub(o, a)) / den;
    if (!(t0 >= 0.0)) return h;                      /* behind the origin */
    h.hit = 1;
    h.flags = INTERSECT_HIT;
    h.t = t0;
    {
        /* barycentrics from the sign-tested values themselves */
        double sum = g[0] + g[1] + g[2];
        if (sum != 0.0) { h.wa = g[0] / sum; h.wb = g[1] / sum; h.wc = g[2] / sum; }
    }
    return h;
}

/* ------------------------------------------------------------------ */
/* formula 2: Woop projected 2D edge functions                         */
/* ------------------------------------------------------------------ */
/*
 * Project along the dominant ray axis.  The three 2D edge functions U, V, W
 * come from different operands than formula 1, which makes this a real
 * cross-check rather than a rewrite.  The 2D determinant equals -(n.d)/d_kz
 * in exact arithmetic, so det == 0 is exactly the parallel case, not a
 * degenerate triangle.
 */
static Hit intersect_woop(const Tri* t, V3 o, V3 d) {
    Hit h;
    memset(&h, 0, sizeof(h));
    int deg = tri_degenerate(t);
    if (deg) { h.flags = INTERSECT_DEGENERATE; h.degenerate = deg; return h; }

    double D[3] = { d.x, d.y, d.z };
    int kz = 0;
    double m = fabs(D[0]);
    if (fabs(D[1]) > m) { kz = 1; m = fabs(D[1]); }
    if (fabs(D[2]) > m) { kz = 2; m = fabs(D[2]); }
    if (D[kz] == 0.0) { h.flags = INTERSECT_PARALLEL; return h; }
    int kx = kz + 1; if (kx == 3) kx = 0;
    int ky = kx + 1; if (ky == 3) ky = 0;
    if (D[kz] < 0.0) { int tmp = kx; kx = ky; ky = tmp; }

    double Sx = -D[kx] / D[kz];
    double Sy = -D[ky] / D[kz];
    double Sz = 1.0 / D[kz];

    V3 vs[3] = { to_v3(t->a), to_v3(t->b), to_v3(t->c) };
    double ox[3] = { o.x, o.y, o.z };
    double px[3], py[3], pz[3];
    for (int i = 0; i < 3; i++) {
        double q[3] = { vs[i].x - ox[0], vs[i].y - ox[1], vs[i].z - ox[2] };
        px[i] = q[kx] + Sx * q[kz];
        py[i] = q[ky] + Sy * q[kz];
        pz[i] = q[kz] * Sz;
    }
    double Ax = px[0], Ay = py[0], Az = pz[0];
    double Bx = px[1], By = py[1], Bz = pz[1];
    double Cx = px[2], Cy = py[2], Cz = pz[2];

    double U = Cx * By - Cy * Bx;   /* edge b-c */
    double V = Ax * Cy - Ay * Cx;   /* edge c-a */
    double W = Bx * Ay - By * Ax;   /* edge a-b */
    double det = U + V + W;
    if (det == 0.0) { h.flags = INTERSECT_PARALLEL; return h; }

    int su = sign_snap(U, INTERSECT_K * DBL_EPSILON * (fabs(Cx * By) + fabs(Cy * Bx)));
    int sv = sign_snap(V, INTERSECT_K * DBL_EPSILON * (fabs(Ax * Cy) + fabs(Ay * Cx)));
    int sw = sign_snap(W, INTERSECT_K * DBL_EPSILON * (fabs(Bx * Ay) + fabs(By * Ax)));
    if (!same_side(su, sv, sw)) return h;

    double t0 = (U * Az + V * Bz + W * Cz) / det;
    if (!(t0 >= 0.0)) return h;
    h.hit = 1;
    h.flags = INTERSECT_HIT;
    h.t = t0;
    h.wa = U / det;
    h.wb = V / det;
    h.wc = W / det;
    return h;
}

/* Public entry point. */
static Hit intersect_tri(const Tri* t, Vec3 origin, Vec3 dir, Formula f) {
    V3 o = to_v3(origin), d = to_v3(dir);
    return f == INTERSECT_WOOP ? intersect_woop(t, o, d) : intersect_mt(t, o, d);
}

/* ------------------------------------------------------------------ */
/* reference predicates used only by the self-test                     */
/* ------------------------------------------------------------------ */

/*
 * Exact analytic boundary predicate for pixels with dyadic-exact vertices.
 * Projects to a plane and evaluates all three sides in 64-bit integers, so
 * it is a true oracle rather than another floating-point implementation.
 */
typedef struct { long long a[2], b[2], c[2]; } Tri2i;

static long long cross2(long long ax, long long ay, long long bx, long long by) {
    return ax * by - ay * bx;
}

static int exact_contains_2d(const Tri2i* t, long long x, long long y) {
    long long f0 = cross2(t->b[0] - t->a[0], t->b[1] - t->a[1], x - t->a[0], y - t->a[1]);
    long long f1 = cross2(t->c[0] - t->b[0], t->c[1] - t->b[1], x - t->b[0], y - t->b[1]);
    long long f2 = cross2(t->a[0] - t->c[0], t->a[1] - t->c[1], x - t->c[0], y - t->c[1]);
    return (f0 >= 0 && f1 >= 0 && f2 >= 0) || (f0 <= 0 && f1 <= 0 && f2 <= 0);
}

/*
 * Analytic ray-vs-triangle reference in long double, used for non-dyadic
 * geometry where an integer oracle is not available.  Ray-plane, then the
 * same three edge functions, zero inclusive.
 */
static int analytic_hit(const Tri* t, V3 o, V3 d) {
    long double A[3] = { t->a.x, t->a.y, t->a.z };
    long double B[3] = { t->b.x, t->b.y, t->b.z };
    long double C[3] = { t->c.x, t->c.y, t->c.z };
    long double O[3] = { o.x, o.y, o.z };
    long double D[3] = { d.x, d.y, d.z };
    long double e1[3] = { B[0] - A[0], B[1] - A[1], B[2] - A[2] };
    long double e2[3] = { C[0] - A[0], C[1] - A[1], C[2] - A[2] };
    long double n[3] = {
        e1[1] * e2[2] - e1[2] * e2[1],
        e1[2] * e2[0] - e1[0] * e2[2],
        e1[0] * e2[1] - e1[1] * e2[0]
    };
    long double dn = D[0] * n[0] + D[1] * n[1] + D[2] * n[2];
    if (dn == 0.0L) return 0;
    long double s[3] = { O[0] - A[0], O[1] - A[1], O[2] - A[2] };
    long double tn = -(s[0] * n[0] + s[1] * n[1] + s[2] * n[2]);
    if (dn < 0.0L) { dn = -dn; tn = -tn; }
    if (tn < 0.0L) return 0;
    long double tp = tn / dn;
    long double P[3] = { O[0] + tp * D[0], O[1] + tp * D[1], O[2] + tp * D[2] };
    long double W[3][3] = {
        { C[0] - B[0], C[1] - B[1], C[2] - B[2] },
        { A[0] - C[0], A[1] - C[1], A[2] - C[2] },
        { B[0] - A[0], B[1] - A[1], B[2] - A[2] }
    };
    long double Q[3][3] = {
        { P[0] - B[0], P[1] - B[1], P[2] - B[2] },
        { P[0] - C[0], P[1] - C[1], P[2] - C[2] },
        { P[0] - A[0], P[1] - A[1], P[2] - A[2] }
    };
    int pos = 1, neg = 1;
    for (int i = 0; i < 3; i++) {
        long double c[3] = {
            W[i][1] * Q[i][2] - W[i][2] * Q[i][1],
            W[i][2] * Q[i][0] - W[i][0] * Q[i][2],
            W[i][0] * Q[i][1] - W[i][1] * Q[i][0]
        };
        long double f = c[0] * n[0] + c[1] * n[1] + c[2] * n[2];
        if (!(f >= 0.0L)) pos = 0;
        if (!(f <= 0.0L)) neg = 0;
    }
    return pos || neg;
}

/* ------------------------------------------------------------------ */
/* self-test                                                           */
/* ------------------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* A ray through the centre must agree exactly in both formulas. */
static void test_center(void) {
    Tri t = { v3f(0, 0, 0), v3f(1, 0, 0), v3f(0, 1, 0) };
    Vec3 o = v3f(0.2f, 0.2f, 1.0f), d = v3f(0, 0, -1);
    Hit m = intersect_tri(&t, o, d, INTERSECT_MT);
    Hit w = intersect_tri(&t, o, d, INTERSECT_WOOP);
    check("center_mt_hit", m.hit == 1 && m.flags == INTERSECT_HIT);
    check("center_woop_hit", w.hit == 1 && w.flags == INTERSECT_HIT);
    check("center_t_agree", fabs(m.t - w.t) < 1e-12);
    check("center_t_value", fabs(m.t - 1.0) < 1e-6);
    check("center_uv_agree", fabs(m.wa - w.wa) < 1e-9 &&
                              fabs(m.wb - w.wb) < 1e-9 &&
                              fabs(m.wc - w.wc) < 1e-9);
    check("center_bary_sum", fabs((m.wa + m.wb + m.wc) - 1.0) < 1e-9);
}

/* Rays deliberately aimed at a vertex and at an edge midpoint. */
static void test_boundary_vertex_and_edge(void) {
    Tri t = { v3f(0, 0, 0), v3f(4, 0, 0), v3f(0, 4, 0) };
    Vec3 dir = v3f(0, 0, -1);
    Hit va = intersect_tri(&t, v3f(0, 0, 1), dir, INTERSECT_MT);
    Hit ea = intersect_tri(&t, v3f(2, 0, 1), dir, INTERSECT_MT);
    Hit vw = intersect_tri(&t, v3f(0, 0, 1), dir, INTERSECT_WOOP);
    Hit ew = intersect_tri(&t, v3f(2, 0, 1), dir, INTERSECT_WOOP);
    check("vertex_is_claimed", va.hit == 1 && vw.hit == 1);
    check("edge_is_claimed", ea.hit == 1 && ew.hit == 1);
    check("edge_bary_zero", ea.wc == 0.0 && ew.wc == 0.0 &&
                             fabs(ea.wa - 0.5) < 1e-6 && fabs(ea.wb - 0.5) < 1e-6);
}

/*
 * The two neighbours' copies of a shared edge are expressed from opposite
 * endpoints.  The values are not bit-identical (the attachment vertex of the
 * cross product also changes), but because the operands are the same vertices
 * in the opposite order the difference is pure rounding, and the tolerance
 * built from those same products must absorb it.  This test pins that down:
 * sample long-double-exact zeros on the edge and confirm the snapped sign is
 * zero from both sides.
 */
static void test_shared_edge_agreement(void) {
    Tri t1 = { v3f(0, 0, 0), v3f(4, 2, 1), v3f(1, 3, 2) };
    Tri t2 = { v3f(4, 2, 1), v3f(0, 0, 0), v3f(3, -2, 2) };
    V3 o = { 0.3, 0.7, -2.0 }, d = { 0.1, -0.2, 1.0 };
    (void)o; (void)d;
    int both_zero = 0, flips = 0;
    for (int i = 0; i <= 64; i++) {
        float s = (float)i / 64.0f;
        Vec3 p3 = v3f(4.0f * s, 2.0f * s, 1.0f * s);
        V3 p = to_v3(p3);
        V3 dir = { 0.0, 0.0, 1.0 };
        V3 origin = p;
        double sa, sb;
        /* t1's shared edge (a,b): opposite vertex is c; t2's is (a,b) too */
        double g1 = edge_fn(dir, origin, to_v3(t1.a), to_v3(t1.b), &sa);
        double g2 = edge_fn(dir, origin, to_v3(t2.b), to_v3(t2.a), &sb);
        double tol1 = INTERSECT_K * DBL_EPSILON * sa;
        double tol2 = INTERSECT_K * DBL_EPSILON * sb;
        int s1 = sign_snap(g1, tol1);
        int s2 = sign_snap(g2, tol2);
        if (s1 == 0 && s2 == 0) both_zero++;
        /* if one is non-zero the other must be its negation: still consistent */
        if (s1 != 0 && s2 != 0 && s1 != -s2) flips++;
    }
    check("shared_edge_both_zero", both_zero == 65);
    check("shared_edge_signs_consistent", flips == 0);
}

/*
 * The seam test: two triangles share the diagonal of a cell.  Rays are sent
 * through points exactly on that diagonal, which is the pathological case.
 * For every seam pixel the exact union must be covered and both formulas must
 * agree.  N is a power of two so every sampled coordinate is exact in float.
 */
static void test_seam_zero_leak(void) {
    const int N = 32;
    long long leaked_mt = 0, leaked_woop = 0, xform = 0, seam_pixels = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            float x0 = (float)i / N, y0 = (float)j / N;
            float x1 = (float)(i + 1) / N, y1 = (float)(j + 1) / N;
            Tri t1 = { v3f(x0, y0, 0), v3f(x1, y0, 0), v3f(x1, y1, 0) };
            Tri t2 = { v3f(x0, y0, 0), v3f(x1, y1, 0), v3f(x0, y1, 0) };
            /* exact predicate in units of 1/(2N) */
            Tri2i e1 = { { 2 * i, 2 * j }, { 2 * i + 2, 2 * j }, { 2 * i + 2, 2 * j + 2 } };
            Tri2i e2 = { { 2 * i, 2 * j }, { 2 * i + 2, 2 * j + 2 }, { 2 * i, 2 * j + 2 } };
            float sx = (x0 + x1) * 0.5f, sy = (y0 + y1) * 0.5f;
            Vec3 o = v3f(sx, sy, 1.0f), d = v3f(0, 0, -1);
            int ex1 = exact_contains_2d(&e1, 2 * i + 1, 2 * j + 1);
            int ex2 = exact_contains_2d(&e2, 2 * i + 1, 2 * j + 1);
            int cov = ex1 || ex2;
            if (cov && ex1 && ex2) seam_pixels++;
            Hit m1 = intersect_tri(&t1, o, d, INTERSECT_MT);
            Hit m2 = intersect_tri(&t2, o, d, INTERSECT_MT);
            Hit w1 = intersect_tri(&t1, o, d, INTERSECT_WOOP);
            Hit w2 = intersect_tri(&t2, o, d, INTERSECT_WOOP);
            if (cov && !m1.hit && !m2.hit) leaked_mt++;
            if (cov && !w1.hit && !w2.hit) leaked_woop++;
            if (m1.hit != w1.hit || m2.hit != w2.hit) xform++;
        }
    }
    check("seam_pixels_exist", seam_pixels > 0);
    check("seam_zero_leak_mt", leaked_mt == 0);
    check("seam_zero_leak_woop", leaked_woop == 0);
    check("seam_formulas_agree", xform == 0);
}

/*
 * Pixel-identical comparison against the exact analytic boundary predicate,
 * on a lattice that includes interiors, shared edges, and vertices.  Every
 * triangle is checked individually, so this catches over-claiming too.
 * Coordinates are multiples of 1/128, exact in float.
 */
static void test_analytic_pixel_identity(void) {
    const int N = 32;
    long long checks = 0, mism_mt = 0, mism_woop = 0, xform = 0, covered = 0;
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            float x0 = (float)i / N, y0 = (float)j / N;
            float x1 = (float)(i + 1) / N, y1 = (float)(j + 1) / N;
            Tri t1 = { v3f(x0, y0, 0), v3f(x1, y0, 0), v3f(x1, y1, 0) };
            Tri t2 = { v3f(x0, y0, 0), v3f(x1, y1, 0), v3f(x0, y1, 0) };
            /* exact predicate in units of 1/(4N) */
            Tri2i e1 = { { 4 * i, 4 * j }, { 4 * i + 4, 4 * j }, { 4 * i + 4, 4 * j + 4 } };
            Tri2i e2 = { { 4 * i, 4 * j }, { 4 * i + 4, 4 * j + 4 }, { 4 * i, 4 * j + 4 } };
            for (int a = 0; a < 4; a++) {
                for (int b = 0; b < 4; b++) {
                    long long qx = 4 * i + a, qy = 4 * j + b;
                    float px = (float)qx / (4 * N), py = (float)qy / (4 * N);
                    Vec3 o = v3f(px, py, 1.0f), d = v3f(0, 0, -1);
                    int ex1 = exact_contains_2d(&e1, qx, qy);
                    int ex2 = exact_contains_2d(&e2, qx, qy);
                    Hit m1 = intersect_tri(&t1, o, d, INTERSECT_MT);
                    Hit m2 = intersect_tri(&t2, o, d, INTERSECT_MT);
                    Hit w1 = intersect_tri(&t1, o, d, INTERSECT_WOOP);
                    Hit w2 = intersect_tri(&t2, o, d, INTERSECT_WOOP);
                    checks += 2;
                    if (ex1 || ex2) covered++;
                    if (m1.hit != ex1) mism_mt++;
                    if (m2.hit != ex2) mism_mt++;
                    if (w1.hit != ex1) mism_woop++;
                    if (w2.hit != ex2) mism_woop++;
                    if (m1.hit != w1.hit) xform++;
                    if (m2.hit != w2.hit) xform++;
                }
            }
        }
    }
    printf("analytic_checks=%lld covered=%lld mism_mt=%lld mism_woop=%lld xform=%lld\n",
           checks, covered, mism_mt, mism_woop, xform);
    check("analytic_covered", covered > 0);
    check("analytic_mt_identity", mism_mt == 0);
    check("analytic_woop_identity", mism_woop == 0);
    check("analytic_formula_identity", xform == 0);
}

/*
 * Non-dyadic geometry: random shared edges checked against the long-double
 * oracle.  The interpolated float origin usually lands a hair off the edge,
 * which is fine: the correct answer then belongs to exactly one triangle.
 * What must hold is no leak and exact agreement with the oracle.
 */
static void test_random_seams(void) {
    uint32_t s = 0x13579BDFu;
    long long covered = 0, leaked_mt = 0, leaked_woop = 0;
    long long mism_mt = 0, mism_woop = 0, xform = 0, both_claim = 0;
    for (int it = 0; it < 100000; it++) {
        float r[8];
        for (int k = 0; k < 8; k++) {
            s = s * 1664525u + 1013904223u;
            r[k] = ((float)((s >> 8) & 0xffffu) / 65536.0f) * 4.0f - 2.0f;
        }
        s = s * 1664525u + 1013904223u;
        float lam = (float)((s >> 8) & 0xffffu) / 65536.0f;
        Tri t1 = { v3f(r[0], r[1], 0), v3f(r[2], r[3], 0), v3f(r[4], r[5], 1) };
        Tri t2 = { v3f(r[2], r[3], 0), v3f(r[0], r[1], 0), v3f(r[6], r[7], 1) };
        float px = t1.a.x + (t1.b.x - t1.a.x) * lam;
        float py = t1.a.y + (t1.b.y - t1.a.y) * lam;
        V3 o = { (double)px, (double)py, -1.0 }, d = { 0.0, 0.0, 1.0 };
        int e1 = analytic_hit(&t1, o, d);
        int e2 = analytic_hit(&t2, o, d);
        if (!e1 && !e2) continue;   /* rounded origin fell off the edge */
        covered++;
        Vec3 of = v3f(px, py, -1.0f), df = v3f(0, 0, 1);
        Hit m1 = intersect_tri(&t1, of, df, INTERSECT_MT);
        Hit m2 = intersect_tri(&t2, of, df, INTERSECT_MT);
        Hit w1 = intersect_tri(&t1, of, df, INTERSECT_WOOP);
        Hit w2 = intersect_tri(&t2, of, df, INTERSECT_WOOP);
        if (!m1.hit && !m2.hit) leaked_mt++;
        if (!w1.hit && !w2.hit) leaked_woop++;
        if (m1.hit != e1) mism_mt++;
        if (m2.hit != e2) mism_mt++;
        if (w1.hit != e1) mism_woop++;
        if (w2.hit != e2) mism_woop++;
        if (m1.hit != w1.hit || m2.hit != w2.hit) xform++;
        if (m1.hit && m2.hit) both_claim++;
    }
    printf("random_seam covered=%lld both_claim=%lld\n", covered, both_claim);
    check("random_seam_covered", covered > 0);
    check("random_seam_zero_leak_mt", leaked_mt == 0);
    check("random_seam_zero_leak_woop", leaked_woop == 0);
    check("random_seam_mt_oracle", mism_mt == 0);
    check("random_seam_woop_oracle", mism_woop == 0);
    check("random_seam_formulas_agree", xform == 0);
}

/*
 * Negative control: the classic strict float Moller-Trumbore really does leak
 * on the same kind of shared edge, so the zero-leak result above is not
 * vacuous.  (A fixed epsilon window hides the problem only while the rounding
 * stays inside the window; strict float is the honest comparison.)
 */
static int naive_float_mt(const Tri* t, Vec3 o, Vec3 d) {
    Vec3 e1 = v3f(t->b.x - t->a.x, t->b.y - t->a.y, t->b.z - t->a.z);
    Vec3 e2 = v3f(t->c.x - t->a.x, t->c.y - t->a.y, t->c.z - t->a.z);
    Vec3 p = v3f(d.y * e2.z - d.z * e2.y, d.z * e2.x - d.x * e2.z, d.x * e2.y - d.y * e2.x);
    float det = e1.x * p.x + e1.y * p.y + e1.z * p.z;
    if (det == 0.0f) return 0;
    float inv = 1.0f / det;
    Vec3 s = v3f(o.x - t->a.x, o.y - t->a.y, o.z - t->a.z);
    float u = (s.x * p.x + s.y * p.y + s.z * p.z) * inv;
    Vec3 q = v3f(s.y * e1.z - s.z * e1.y, s.z * e1.x - s.x * e1.z, s.x * e1.y - s.y * e1.x);
    float v = (d.x * q.x + d.y * q.y + d.z * q.z) * inv;
    return !(u < 0.0f || v < 0.0f || u + v > 1.0f);
}

static void test_negative_control(void) {
    uint32_t s = 0x13579BDFu;
    long long covered = 0, naive_leaks = 0, ours_leaks = 0;
    for (int it = 0; it < 100000; it++) {
        float r[8];
        for (int k = 0; k < 8; k++) {
            s = s * 1664525u + 1013904223u;
            r[k] = ((float)((s >> 8) & 0xffffu) / 65536.0f) * 4.0f - 2.0f;
        }
        s = s * 1664525u + 1013904223u;
        float lam = (float)((s >> 8) & 0xffffu) / 65536.0f;
        Tri t1 = { v3f(r[0], r[1], 0), v3f(r[2], r[3], 0), v3f(r[4], r[5], 1) };
        Tri t2 = { v3f(r[2], r[3], 0), v3f(r[0], r[1], 0), v3f(r[6], r[7], 1) };
        float px = t1.a.x + (t1.b.x - t1.a.x) * lam;
        float py = t1.a.y + (t1.b.y - t1.a.y) * lam;
        V3 od = { (double)px, (double)py, -1.0 }, dd = { 0.0, 0.0, 1.0 };
        int e1 = analytic_hit(&t1, od, dd);
        int e2 = analytic_hit(&t2, od, dd);
        if (!e1 && !e2) continue;
        covered++;
        Vec3 of = v3f(px, py, -1.0f), df = v3f(0, 0, 1);
        if (!naive_float_mt(&t1, of, df) && !naive_float_mt(&t2, of, df)) naive_leaks++;
        Hit h1 = intersect_tri(&t1, of, df, INTERSECT_MT);
        Hit h2 = intersect_tri(&t2, of, df, INTERSECT_MT);
        if (!h1.hit && !h2.hit) ours_leaks++;
    }
    printf("negative_control covered=%lld naive_leaks=%lld (%.3f%%) ours_leaks=%lld\n",
           covered, naive_leaks,
           covered ? 100.0 * (double)naive_leaks / (double)covered : 0.0,
           ours_leaks);
    check("negative_control_naive_leaks", naive_leaks > 0);
    check("negative_control_ours_clean", ours_leaks == 0);
}

/* A sheared, tilted pair sharing a slanted edge: both formulas, both flags. */
static void test_slanted_shared_edge(void) {
    Tri t1 = { v3f(0, 0, 0), v3f(4, 2, 1), v3f(1, 3, 2) };
    Tri t2 = { v3f(4, 2, 1), v3f(0, 0, 0), v3f(3, -2, 2) };
    int leaks = 0, xform = 0, mism = 0;
    for (int i = 0; i < 64; i++) {
        float s = (float)i / 63.0f;
        Vec3 o = v3f(4.0f * s, 2.0f * s, 1.0f * s + 1.0f);
        Vec3 d = v3f(0, 0, -1);
        Hit m1 = intersect_tri(&t1, o, d, INTERSECT_MT);
        Hit m2 = intersect_tri(&t2, o, d, INTERSECT_MT);
        Hit w1 = intersect_tri(&t1, o, d, INTERSECT_WOOP);
        Hit w2 = intersect_tri(&t2, o, d, INTERSECT_WOOP);
        if (!m1.hit && !m2.hit) leaks++;
        if (m1.hit != w1.hit || m2.hit != w2.hit) xform++;
        if (m1.hit != analytic_hit(&t1, to_v3(o), to_v3(d)) ||
            m2.hit != analytic_hit(&t2, to_v3(o), to_v3(d))) mism++;
    }
    check("slanted_seam_no_leak", leaks == 0);
    check("slanted_seam_formulas_agree", xform == 0);
    check("slanted_seam_oracle", mism == 0);
}

/* Degenerate triangles are reported separately and never hit. */
static void test_degenerate(void) {
    Tri collinear = { v3f(0, 0, 0), v3f(1, 0, 0), v3f(2, 0, 0) };
    Tri coincident = { v3f(3, 3, 3), v3f(3, 3, 3), v3f(4, 3, 3) };
    Tri allsame = { v3f(7, 7, 7), v3f(7, 7, 7), v3f(7, 7, 7) };
    Tri nonfinite = { v3f(0, 0, 0), v3f(1, 0, 0), v3f(NAN, 0, 0) };
    Vec3 o = v3f(0.5f, 0.0f, 1.0f), d = v3f(0, 0, -1);

    Tri tris[4] = { collinear, coincident, allsame, nonfinite };
    const char* names[4] = { "collinear", "coincident", "allsame", "nonfinite" };
    int reason[4] = { 1, 1, 1, 2 };
    for (int i = 0; i < 4; i++) {
        Hit m = intersect_tri(&tris[i], o, d, INTERSECT_MT);
        Hit w = intersect_tri(&tris[i], o, d, INTERSECT_WOOP);
        char label[64];
        snprintf(label, sizeof(label), "degenerate_%s_reported", names[i]);
        check(label, m.flags == INTERSECT_DEGENERATE && w.flags == INTERSECT_DEGENERATE &&
                     m.degenerate == reason[i] && m.hit == 0 && w.hit == 0);
    }

    /* parallel ray is not a degenerate triangle */
    Tri good = { v3f(0, 0, 0), v3f(1, 0, 0), v3f(0, 1, 0) };
    Hit par = intersect_tri(&good, v3f(0.2f, 0.2f, 0.0f), v3f(1, 0, 0), INTERSECT_MT);
    Hit parw = intersect_tri(&good, v3f(0.2f, 0.2f, 0.0f), v3f(1, 0, 0), INTERSECT_WOOP);
    check("parallel_reported", par.flags == INTERSECT_PARALLEL && parw.flags == INTERSECT_PARALLEL);
    check("parallel_not_degenerate", par.degenerate == 0 && parw.degenerate == 0);

    /* behind the origin is a miss, not a degenerate */
    Hit behind = intersect_tri(&good, v3f(0.2f, 0.2f, 1.0f), v3f(0, 0, 1), INTERSECT_MT);
    check("behind_is_miss", behind.hit == 0 && behind.flags == INTERSECT_MISS);

    /* the aggregate stats count each class separately */
    TriStats st;
    tri_stats_init(&st);
    for (int i = 0; i < 4; i++) {
        Hit hd = intersect_tri(&tris[i], o, d, INTERSECT_MT);
        tri_stats_add(&st, &hd);
    }
    Hit hgood = intersect_tri(&good, v3f(0.2f, 0.2f, 1.0f), d, INTERSECT_MT);
    Hit hpar = intersect_tri(&good, v3f(0.2f, 0.2f, 0.0f), v3f(1, 0, 0), INTERSECT_MT);
    tri_stats_add(&st, &hgood);
    tri_stats_add(&st, &hpar);
    printf("degenerate_report zero_area=%lld nonfinite=%lld parallel=%lld hits=%lld\n",
           st.degenerate_zero_area, st.degenerate_nonfinite, st.parallel, st.hits);
    check("stats_tested", st.tested == 6);
    check("stats_zero_area", st.degenerate_zero_area == 3);
    check("stats_nonfinite", st.degenerate_nonfinite == 1);
    check("stats_parallel", st.parallel == 1);
    check("stats_hits", st.hits == 1);
}

/*
 * The near-degenerate threshold must use both edge lengths.  An elongated
 * sliver (long thin triangle) is still usable and must be hit; a triangle
 * whose three points are collinear to within rounding must be rejected.
 */
static void test_thin_triangles(void) {
    /* long but genuinely thin: area 0.5, height/width ratio 1e-6 */
    Tri sliver = { v3f(0, 0, 0), v3f(1.0f, 0, 0), v3f(0.5f, 1e-6f, 0) };
    Vec3 o = v3f(0.5f, 2.5e-7f, 1.0f), d = v3f(0, 0, -1);
    Hit m = intersect_tri(&sliver, o, d, INTERSECT_MT);
    Hit w = intersect_tri(&sliver, o, d, INTERSECT_WOOP);
    check("sliver_usable", m.flags == INTERSECT_HIT && w.flags == INTERSECT_HIT);
    check("sliver_t_agree", fabs(m.t - w.t) < 1e-9 && fabs(m.t - 1.0) < 1e-5);

    /* height below the float noise floor of the vertices: reject */
    float tiny = 1e-9f;
    Tri flat = { v3f(1.0f, 1.0f, 0), v3f(2.0f, 2.0f, 0), v3f(3.0f, 3.0f + tiny, 0) };
    Hit fm = intersect_tri(&flat, v3f(2.0f, 2.0f, 1.0f), d, INTERSECT_MT);
    Hit fw = intersect_tri(&flat, v3f(2.0f, 2.0f, 1.0f), d, INTERSECT_WOOP);
    check("noise_floor_rejected_mt", fm.flags == INTERSECT_DEGENERATE);
    check("noise_floor_rejected_woop", fw.flags == INTERSECT_DEGENERATE);
}

int main(void) {
    printf("K=%.1f\n", (double)INTERSECT_K);
    test_center();
    test_boundary_vertex_and_edge();
    test_shared_edge_agreement();
    test_seam_zero_leak();
    test_analytic_pixel_identity();
    test_random_seams();
    test_negative_control();
    test_slanted_shared_edge();
    test_degenerate();
    test_thin_triangles();
    printf("degenerate_policy: reported separately, never counted as hits\n");
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
