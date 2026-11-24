
/* light: rectangular area light sampled on a uniform N x N grid.
   Each cell's radiance is weighted by the cosine at its center, and the
   shadow mask is the cosine-weighted visible fraction of the light.
   This replaces edge-midpoint sampling, which never sampled the light
   interior and left a dark spot at the mask center.
   Zero-area lights degrade to an explicit point-light fallback instead
   of dividing by zero. */
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct { double x, y, z; } V3;
static V3 sub3(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static double dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 mad3(V3 c, V3 u, double su, V3 v, double sv) {
    V3 r = { c.x + su * u.x + sv * v.x,
             c.y + su * u.y + sv * v.y,
             c.z + su * u.z + sv * v.z };
    return r;
}

typedef struct {
    V3 c;          /* center */
    V3 n;          /* unit normal, points toward the lit side */
    V3 u, v;       /* in-plane unit axes */
    double hx, hy; /* half extents */
} AreaLight;

typedef int (*OccluderFn)(V3 from, V3 to, void* ctx);

/* center of cell (i,j) in an n x n grid; local fractions stay in (-1,1),
   so every sample lies strictly inside the light */
static V3 cell_center(const AreaLight* L, int n, int i, int j) {
    double fu = 2.0 * ((double)i + 0.5) / n - 1.0;
    double fv = 2.0 * ((double)j + 0.5) / n - 1.0;
    return mad3(L->c, L->u, fu * L->hx, L->v, fv * L->hy);
}

/* cosine-weighted visibility of the light as seen from p, on an n x n
   grid. Degenerate cases are explicit, never a division by zero:
   - zero area or n < 1: point-light fallback, one shadow ray, hard 0/1;
   - no contributing cell (light fully back-facing): defined as 0. */
static double visibility(const AreaLight* L, V3 p, int n, OccluderFn occ, void* ctx) {
    double area = 4.0 * L->hx * L->hy;
    if (area <= 0.0 || n < 1) {
        V3 d = sub3(p, L->c);
        double len = sqrt(dot3(d, d));
        if (len <= 0.0) return 0.0;
        if (dot3(L->n, d) / len <= 0.0) return 0.0;
        return (occ && occ(p, L->c, ctx)) ? 0.0 : 1.0;
    }
    double acell = area / ((double)n * n);
    double wsum = 0.0, wvis = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            V3 q = cell_center(L, n, i, j);
            V3 d = sub3(p, q);
            double len = sqrt(dot3(d, d));
            double c = dot3(L->n, d) / len;
            if (c <= 0.0) continue;
            double w = c * acell;
            wsum += w;
            if (!occ || !occ(p, q, ctx)) wvis += w;
        }
    }
    if (wsum <= 0.0) return 0.0;
    return wvis / wsum;
}

/* straight-edge occluder: blocks cells whose center has x > x0 */
typedef struct { double x0; } EdgeCtx;
static int edge_occ(V3 from, V3 to, void* ctx_) {
    EdgeCtx* e = ctx_;
    (void)from;
    return to.x > e->x0;
}

static int full_occ(V3 from, V3 to, void* ctx) {
    (void)from; (void)to; (void)ctx;
    return 1;
}

/* wall occluder: half-plane z = z_wall, x < 0, between light and receiver */
typedef struct { double z; } WallCtx;
static int wall_occ(V3 from, V3 to, void* ctx_) {
    WallCtx* w = ctx_;
    double dz = to.z - from.z;
    if (fabs(dz) < 1e-12) return 0;
    double t = (w->z - from.z) / dz;
    if (t <= 0.0 || t >= 1.0) return 0;
    double x = from.x + t * (to.x - from.x);
    return x < 0.0;
}

/* analytic reference for the straight-edge mask. With the light at
   height h, half extent a, weight cos = h/sqrt(h^2+x^2+y^2), the inner
   integral is closed form:
     I(x0) = int_{-a}^{x0} 2h asinh(a/sqrt(h^2+x^2)) dx
   V_ref = I(x0)/I(a). The remaining 1D integral uses a 2^14-point
   midpoint rule (error < 1e-9 for the test geometry). */
static double edge_integrand(double x, double h, double a) {
    return 2.0 * h * asinh(a / sqrt(h * h + x * x));
}
static double edge_integral(double x0, double h, double a) {
    int m = 1 << 14;
    double lo = -a, hi = x0, sum = 0.0, dx;
    if (hi <= lo) return 0.0;
    dx = (hi - lo) / m;
    for (int i = 0; i < m; i++) sum += edge_integrand(lo + (i + 0.5) * dx, h, a);
    return sum * dx;
}

/* upper bound for |V_N - V_ref| for a straight vertical edge:
     |V_N - V| <= 2(S + E) / W
   S = cos_max * strip area (one cell wide) bounds the weight of cells
   straddling the edge; E is the midpoint-rule smoothing error of cos
   over the light; W is the total weight. */
static double edge_bound(double h, double a, int n, double wref) {
    double cell = 2.0 * a / n;
    double area = 4.0 * a * a;
    double cos_max = 1.0;              /* shading point below light center */
    double s = cos_max * cell * 2.0 * a;
    double fcurv = 4.0 / (h * h);      /* |f_xx|+|f_yy| <= 4/h^2, f = cos */
    double e = area * cell * cell * fcurv / 24.0;
    return 2.0 * (s + e) / wref;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    AreaLight L;
    L.c = (V3){ 0, 0, 10 }; L.n = (V3){ 0, 0, -1 };
    L.u = (V3){ 1, 0, 0 };  L.v = (V3){ 0, 1, 0 };
    L.hx = 0.5; L.hy = 0.5;
    V3 p = { 0, 0, 0 };

    /* grid regression: samples cover the interior, not just edges */
    int inside = 1;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) {
            V3 q = cell_center(&L, 4, i, j);
            if (fabs(q.x) >= 0.5 || fabs(q.y) >= 0.5) inside = 0;
        }
    check("centers_inside", inside);
    V3 qc = cell_center(&L, 3, 1, 1);
    check("center_covered", fabs(qc.x) < 1e-12 && fabs(qc.y) < 1e-12);

    /* analytic references: exact values */
    check("no_occluder", visibility(&L, p, 8, NULL, NULL) == 1.0);
    check("full_occluder", visibility(&L, p, 8, full_occ, NULL) == 0.0);
    EdgeCtx e0 = { 0.0 };
    check("half_exact", fabs(visibility(&L, p, 8, edge_occ, &e0) - 0.5) < 1e-12);

    /* off-center straight edge vs semi-analytic reference, with bound */
    double h = 10.0, a = 0.5, x0 = 0.17;
    double wref = edge_integral(a, h, a);
    double vref = edge_integral(x0, h, a) / wref;
    EdgeCtx ex = { x0 };
    double v8 = visibility(&L, p, 8, edge_occ, &ex);
    double v32 = visibility(&L, p, 32, edge_occ, &ex);
    double b8 = edge_bound(h, a, 8, wref);
    double b32 = edge_bound(h, a, 32, wref);
    printf("edge_ref=%.8f\n", vref);
    printf("edge_err8=%.8f\n", fabs(v8 - vref));
    printf("edge_bound8=%.8f\n", b8);
    printf("edge_err32=%.8f\n", fabs(v32 - vref));
    printf("edge_bound32=%.8f\n", b32);
    check("edge_within_bound8", fabs(v8 - vref) <= b8);
    check("edge_within_bound32", fabs(v32 - vref) <= b32);
    check("bound_shrinks", b32 < b8 / 3.0); /* bound is O(1/n) */
    check("grid_converges", fabs(v32 - vref) <= fabs(v8 - vref) + 1e-9);

    /* penumbra width: wall at z=5 halfway to the light, receiver slides
       along x. Analytic width w = s*d2/d1 = 1*5/5 = 1. Grid n quantizes
       the transition into ~n steps; edge position error <= w/n. */
    WallCtx wall = { 5.0 };
    int n = 16;
    double dx = 0.005, w = 1.0;
    double x0m = 0, x25 = 0, x75 = 0, x1m = 0;
    int g0 = 0, g25 = 0, g75 = 0, g1 = 0;
    for (double x = -0.8; x <= 0.8; x += dx) {
        V3 q = { x, 0, 0 };
        double v = visibility(&L, q, n, wall_occ, &wall);
        if (v > 1e-9 && !g0) { x0m = x; g0 = 1; }
        if (v >= 0.25 && !g25) { x25 = x; g25 = 1; }
        if (v >= 0.75 && !g75) { x75 = x; g75 = 1; }
        if (v >= 1.0 - 1e-9 && !g1) { x1m = x; g1 = 1; }
    }
    double tol = 2.0 * (w / n + dx);
    printf("penumbra_ref=%.8f\n", w);
    printf("penumbra_meas=%.8f\n", x1m - x0m);
    printf("penumbra_bound=%.8f\n", tol);
    check("penumbra_swept", g0 && g25 && g75 && g1);
    check("penumbra_full", fabs((x1m - x0m) - w) <= tol);
    check("penumbra_mid", fabs((x75 - x25) - 0.5 * w) <= tol);

    /* degenerate cases: explicit fallback, never a division by zero */
    AreaLight pl = L; pl.hx = 0.0; pl.hy = 0.0;
    double vlit = visibility(&pl, p, 8, NULL, NULL);
    double vshd = visibility(&pl, p, 8, full_occ, NULL);
    check("zero_area_lit", vlit == 1.0 && !isnan(vlit));
    check("zero_area_shadow", vshd == 0.0 && !isnan(vshd));
    EdgeCtx epos = { 0.17 }, eneg = { -0.1 };
    check("zero_area_hard_lit", visibility(&pl, p, 8, edge_occ, &epos) == 1.0);
    check("zero_area_hard_shd", visibility(&pl, p, 8, edge_occ, &eneg) == 0.0);
    check("zero_grid_fallback", visibility(&L, p, 0, NULL, NULL) == 1.0);
    AreaLight back = L; back.n = (V3){ 0, 0, 1 };
    double vb = visibility(&back, p, 8, NULL, NULL);
    check("backfacing_zero", vb == 0.0 && !isnan(vb));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
