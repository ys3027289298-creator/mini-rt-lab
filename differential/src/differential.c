
/* differential: ray differential propagation for texture-filter footprints.

   Why this exists: after one specular reflection a naive tracer either keeps
   the camera pixel footprint unchanged or rescales it by the raw direction
   derivative times travel distance. Both ignore the surface curvature term
   and the tangent-plane projection, so the texture LOD is wrong -- over-blurred
   on convex mirrors, aliased at grazing angles.

   This module implements the Igehy 1999 ray differential machinery:

     transfer to a surface  (Igehy eq. 12, surface-normal form)
       dt = -((dP + t dD).N) / (D.N)
       dP' = dP + t dD + dt D                 (tangent to the surface)

     reflection             (Igehy eq. 17)
       r  = D - 2 (D.N) N
       dr = dD - 2 [ (dD.N + D.dN) N + (D.N) dN ]

     footprint              (Igehy eq. 18-21, scalar filter width)
       w = 0.5 (|dP/dx| + |dP/dy|)

   Every returned differential is checked: non-finite values are rejected,
   negative squared lengths are rejected (defensive against corrupted
   externally-built RayDiff structs), extreme curvature is reported and the
   footprint truncated to a caller-supplied cap.

   The self-test compares propagated footprints against per-pixel finite
   differences and verifies an explicit upper bound on the discrepancy. */
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct { float x, y, z; } V3;

/* Status codes. 0 is success; every failure leaves outputs untouched. */
#define RD_OK         0
#define RD_ERR_NAN    1   /* non-finite (NaN/Inf) value in ray or derivative */
#define RD_ERR_NEG    2   /* negative squared differential length            */
#define RD_ERR_PARAM  3   /* t < 0, radius <= 0, cap <= 0, zero direction    */
#define RD_ERR_DEGEN  4   /* UV tangent basis not invertible                 */
#define RD_WARN_CLAMP 5   /* footprint exceeded the cap and was rescaled     */

typedef struct {
    V3 o;          /* ray origin            */
    V3 d;          /* ray direction (unit)  */
    V3 dx, dy;     /* d(origin)/d(pixel)    */
    V3 ddx, ddy;   /* d(direction)/d(pixel) */
} RayDiff;

/* ------------------------------------------------------------------ vec3 */

static V3 vadd(V3 a, V3 b) { V3 r = { a.x+b.x, a.y+b.y, a.z+b.z }; return r; }
static V3 vsub(V3 a, V3 b) { V3 r = { a.x-b.x, a.y-b.y, a.z-b.z }; return r; }
static V3 vmul(V3 a, float s) { V3 r = { a.x*s, a.y*s, a.z*s }; return r; }
static float vdot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static float vlen2(V3 a) { return vdot(a, a); }
static float vlen(V3 a) { return sqrtf(vlen2(a)); }
static V3 vnorm(V3 a) {
    float l = vlen(a);
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x/l; r.y = a.y/l; r.z = a.z/l; }
    return r;
}
static int vfinite(V3 a) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
}

/* ------------------------------------------------------------------ guard */

/* Reject negative squared lengths: a valid Euclidean derivative can never
   produce one, so seeing it means the RayDiff was corrupted upstream. */
static int neg_len(float s2) { return s2 < 0.0f; }

static int check_ray(const RayDiff* rd) {
    if (!vfinite(rd->o) || !vfinite(rd->d) ||
        !vfinite(rd->dx) || !vfinite(rd->dy) ||
        !vfinite(rd->ddx) || !vfinite(rd->ddy))
        return RD_ERR_NAN;
    if (neg_len(vlen2(rd->dx)) || neg_len(vlen2(rd->dy)) ||
        neg_len(vlen2(rd->ddx)) || neg_len(vlen2(rd->ddy)))
        return RD_ERR_NEG;
    if (vlen(rd->d) <= 0.0f) return RD_ERR_PARAM;
    return RD_OK;
}

/* -------------------------------------------------------------- transfer */

/* Propagate ray differentials onto a hit surface.
   N: unit surface normal (either sign), t: ray distance,
   max_width: hard cap on the output footprint (<=0 -> RD_ERR_PARAM).
   Outputs px/py are tangent-plane position differentials at the hit. */
static int rd_transfer(const RayDiff* rd, V3 N, float t,
                       float max_width, V3* px, V3* py, int* clamped) {
    int rc = check_ray(rd);
    if (rc != RD_OK) return rc;
    if (!vfinite(N)) return RD_ERR_NAN;
    if (t < 0.0f || max_width <= 0.0f) return RD_ERR_PARAM;
    float dn = vdot(rd->d, N);
    if (fabsf(dn) < 1e-20f) return RD_ERR_PARAM;       /* edge-on hit */

    /* Each pixel axis: dt keeps the differentiated hit point on the
       surface. Projecting the free-flight differential onto N (Igehy
       eq. 12, surface form) gives the general rule
         dt = -((dP + t dD).N) / (D.N);
       the equivalent D-form dt = -(dP.D + t dD.D)/(D.N) only holds when
       the incoming dP is already tangent, which fails for pinhole rays
       with zero origin differential. */
    V3 pin[2]  = { rd->dx, rd->dy };
    V3 ddin[2] = { rd->ddx, rd->ddy };
    V3 out[2];
    for (int i = 0; i < 2; i++) {
        V3 flight = vadd(pin[i], vmul(ddin[i], t));
        float dt = -vdot(flight, N) / dn;
        out[i] = vadd(flight, vmul(rd->d, dt));
    }

    float w = 0.5f * (vlen(out[0]) + vlen(out[1]));
    if (neg_len(vlen2(out[0])) || neg_len(vlen2(out[1]))) return RD_ERR_NEG;
    *clamped = 0;
    if (w > max_width) {
        /* Extreme curvature / long path: cap the footprint rather than
           poisoning the texture cache. Direction itself is never touched. */
        float s = max_width / w;
        out[0] = vmul(out[0], s);
        out[1] = vmul(out[1], s);
        *clamped = 1;
        w = max_width;
    }
    *px = out[0];
    *py = out[1];
    return *clamped ? RD_WARN_CLAMP : RD_OK;
}

/* -------------------------------------------------------------- reflect */

/* Specular reflection with differential.
   p: hit point, n: unit shading normal (either sign),
   outward: unit outward radial normal of the curved surface (centre -> p;
   pass n with radius=INFINITY for a flat surface), radius: curvature
   radius (sphere radius; INFINITY = plane), d: incident unit direction,
   px/py: tangent position differentials from rd_transfer.
   The dN sign is resolved against n so the result does not depend on
   which way the shading normal is flipped. */
static int rd_reflect(V3 p, V3 n, V3 outward, float radius, V3 d,
                      V3 px, V3 py, RayDiff* out) {
    if (!vfinite(p) || !vfinite(n) || !vfinite(outward) || !vfinite(d) ||
        !vfinite(px) || !vfinite(py))
        return RD_ERR_NAN;
    if (!(radius > 0.0f)) return RD_ERR_PARAM;

    V3 r = vsub(d, vmul(n, 2.0f * vdot(d, n)));

    /* Normal derivative. Sphere: |dN| = |dP|/R; sign follows the shading
       normal's orientation. Plane: dN = 0. */
    V3 dnx, dny;
    if (isinf(radius)) {
        V3 zero = { 0, 0, 0 };
        dnx = dny = zero;
    } else {
        /* Outward radial normal: dN = dP/R. The shading normal may be
           flipped, in which case its derivative flips too: dN = -dP/R.
           n.outward < 0 means the shading normal is the flipped one. */
        float sign = (vdot(n, outward) >= 0.0f) ? 1.0f : -1.0f;
        dnx = vmul(px, sign / radius);
        dny = vmul(py, sign / radius);
    }

    V3 dnin[2] = { dnx, dny };
    V3 dr[2];
    for (int i = 0; i < 2; i++) {
        /* dr = dD - 2[(dD.n + d.dN) n + (d.n) dN].
           The incident direction derivative is zero at a transferred hit
           (the tangent-plane projection absorbs it), so dD = 0 here. */
        float ddn = vdot(d, dnin[i]);
        V3 term = vadd(vmul(n, ddn), vmul(dnin[i], vdot(d, n)));
        dr[i] = vmul(term, -2.0f);
    }

    memset(out, 0, sizeof(*out));
    out->o = p;
    out->d = r;
    out->dx = px;
    out->dy = py;
    out->ddx = dr[0];
    out->ddy = dr[1];
    return check_ray(out);
}

/* -------------------------------------------------------------- measure */

/* Scalar footprint width and the two tangent axis lengths. */
static int rd_measure(V3 px, V3 py, float* width, float* ax, float* ay) {
    if (!vfinite(px) || !vfinite(py)) return RD_ERR_NAN;
    float lx = vlen2(px), ly = vlen2(py);
    if (neg_len(lx) || neg_len(ly)) return RD_ERR_NEG;
    *ax = sqrtf(lx);
    *ay = sqrtf(ly);
    *width = 0.5f * (*ax + *ay);
    return RD_OK;
}

/* Texture-space filter widths (texels) for a surface parameterised by two
   tangent vectors e1, e2 spanning one texel each. Solves the 2x2 Gram
   system instead of a divide-by-component that breaks at grazing angles. */
static int rd_uv_widths(V3 px, V3 py, V3 e1, V3 e2,
                        float* wu, float* wv) {
    if (!vfinite(px) || !vfinite(py) || !vfinite(e1) || !vfinite(e2))
        return RD_ERR_NAN;
    float a = vdot(e1, e1), b = vdot(e1, e2), c = vdot(e2, e2);
    float det = a * c - b * b;
    if (!(det > 1e-20f)) return RD_ERR_DEGEN;
    float ux, uy, vx, vy;
    float r1x = vdot(px, e1), r2x = vdot(px, e2);
    float r1y = vdot(py, e1), r2y = vdot(py, e2);
    ux = ( c * r1x - b * r2x) / det;
    vx = (-b * r1x + a * r2x) / det;
    uy = ( c * r1y - b * r2y) / det;
    vy = (-b * r1y + a * r2y) / det;
    *wu = 0.5f * (fabsf(ux) + fabsf(uy));
    *wv = 0.5f * (fabsf(vx) + fabsf(vy));
    return RD_OK;
}

/* ============================================================ self-test */

typedef struct { double x, y, z; } D3;
static D3 dadd(D3 a, D3 b) { return (D3){a.x+b.x,a.y+b.y,a.z+b.z}; }
static D3 dsub(D3 a, D3 b) { return (D3){a.x-b.x,a.y-b.y,a.z-b.z}; }
static D3 dmul(D3 a, double s) { return (D3){a.x*s,a.y*s,a.z*s}; }
static double ddot(D3 a, D3 b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
static double dlen(D3 a) { return sqrt(ddot(a,a)); }
static D3 dnorm(D3 a) { double l = dlen(a); return (D3){a.x/l,a.y/l,a.z/l}; }

/* Double-precision per-pixel ground truth: orthographic ray (u,0,v)->+y,
   sphere at (0,4,0), reflect, hit receiver plane y=2. */
static D3 ground_truth_Q(double u, double v, double R) {
    D3 o = { u, 0.0, v }, d = { 0, 1, 0 };
    D3 C = { 0, 4, 0 };
    D3 oc = dsub(o, C);
    double B = ddot(oc, d);
    double cc = ddot(oc, oc) - R * R;
    double disc = B * B - cc;
    double t = -B - sqrt(disc);
    D3 P = dadd(o, dmul(d, t));
    D3 N = dnorm(dsub(P, C));
    D3 r = dsub(d, dmul(N, 2.0 * ddot(d, N)));
    double s = (2.0 - P.y) / r.y;
    return dadd(P, dmul(r, s));
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {

    /* ---- camera: orthographic, one world unit per pixel ---------------- */
    RayDiff cam;
    memset(&cam, 0, sizeof(cam));
    cam.o = (V3){ 0, 0, 0 };
    cam.d = (V3){ 0, 1, 0 };
    cam.dx = (V3){ 1, 0, 0 };
    cam.dy = (V3){ 0, 0, 1 };

    /* ============ 1. flat mirror: curvature term must be zero ========== */
    /* plane y=2, normal facing the incoming ray */
    V3 pP = { 0, 2, 0 }, pN = { 0, 1, 0 };
    V3 ppx, ppy; int clamp_flag = 0;
    int rc = rd_transfer(&cam, pN, 2.0f, 1e9f, &ppx, &ppy, &clamp_flag);
    check("plane_transfer_ok", rc == RD_OK);
    RayDiff refl_plane;
    rc = rd_reflect(pP, pN, pN, INFINITY, cam.d, ppx, ppy, &refl_plane);
    check("plane_reflect_ok", rc == RD_OK);
    check("plane_dir_flipped", fabsf(refl_plane.d.y + 1.0f) < 1e-6f);
    check("plane_normal_diff_zero",
          vlen(refl_plane.ddx) < 1e-6f && vlen(refl_plane.ddy) < 1e-6f);
    /* receiver plane one unit behind the mirror */
    V3 rN = { 0, -1, 0 };
    V3 qpx, qpy;
    rc = rd_transfer(&refl_plane, rN, 1.0f, 1e9f,
                     &qpx, &qpy, &clamp_flag);
    float w_plane, ax, ay;
    rd_measure(qpx, qpy, &w_plane, &ax, &ay);
    printf("plane_width=%.6f (exact 1.0)\n", w_plane);
    check("plane_width_exact", fabsf(w_plane - 1.0f) < 1e-5f);

    /* ============ 2. unit sphere centre: analytic width is exactly 3 ==== */
    const float R = 1.0f;
    V3 sP = { 0, 3, 0 }, sN = { 0, -1, 0 };
    V3 spx, spy;
    rc = rd_transfer(&cam, sN, 3.0f, 1e9f, &spx, &spy, &clamp_flag);
    check("sphere_transfer_ok", rc == RD_OK);
    float w_surf;
    rd_measure(spx, spy, &w_surf, &ax, &ay);
    check("sphere_surface_footprint_1", fabsf(w_surf - 1.0f) < 1e-6f);

    RayDiff refl_sph;
    V3 sOut = { 0, -1, 0 };  /* sphere centre (0,4,0) -> hit P (0,3,0) */
    rc = rd_reflect(sP, sN, sOut, R, cam.d, spx, spy, &refl_sph);
    check("sphere_reflect_ok", rc == RD_OK);
    /* dD at centre: dr = -2 (D.N) dN = -2(-1)(1,0,0) = (2,0,0) */
    check("sphere_curvature_term",
          fabsf(refl_sph.ddx.x - 2.0f) < 1e-6f &&
          fabsf(refl_sph.ddy.z - 2.0f) < 1e-6f);

    V3 sqx, sqy;
    rc = rd_transfer(&refl_sph, rN, 1.0f, 1e9f,
                     &sqx, &sqy, &clamp_flag);
    float w_sph;
    rd_measure(sqx, sqy, &w_sph, &ax, &ay);
    printf("sphere_center_width=%.6f (exact 3.0, growth x%.3f)\n",
           w_sph, w_sph / w_surf);
    check("sphere_center_width_exact3", fabsf(w_sph - 3.0f) < 2e-5f);
    check("footprint_growth_reported", w_sph > w_surf * 2.99f);

    /* ============ 3. per-pixel FD ground truth + explicit error bound === */
    /* Central differences in double precision. Truncation error of a
       central difference is O(delta^2); Richardson: if C is the delta-step
       estimate and F the delta/2 estimate, |F - truth| <= |C-F|/3,
       hence |C - truth| <= 4/3 |C-F|. Float roundoff of the propagated
       result is bounded separately by B_fp. */
    const double delta = 2e-3;
    double wfd_c = 0.0, wfd_f = 0.0;
    for (int pass = 0; pass < 2; pass++) {
        double h = pass == 0 ? delta : delta / 2.0;
        D3 qxp = ground_truth_Q( h, 0, R);
        D3 qxm = ground_truth_Q(-h, 0, R);
        D3 qzp = ground_truth_Q(0,  h, R);
        D3 qzm = ground_truth_Q(0, -h, R);
        D3 gx = dmul(dsub(qxp, qxm), 1.0 / (2 * h));
        D3 gz = dmul(dsub(qzp, qzm), 1.0 / (2 * h));
        double w = 0.5 * (dlen(gx) + dlen(gz));
        if (pass == 0) wfd_c = w; else wfd_f = w;
    }
    double b_trunc = 4.0 / 3.0 * fabs(wfd_c - wfd_f);
    double b_fp = 64.0 * 1.2e-7 * (double)w_sph + 1e-6;
    double err = fabs((double)w_sph - wfd_f);
    printf("fd_coarse=%.8f fd_fine=%.8f analytic=%.8f\n",
           wfd_c, wfd_f, w_sph);
    printf("error=%.3e bound=%.3e (trunc=%.3e fp=%.3e)\n",
           err, b_trunc + b_fp, b_trunc, b_fp);
    check("error_within_bound", err <= b_trunc + b_fp);
    check("bound_is_tight", b_trunc > 0.0 && b_trunc < 1e-4);

    /* same check off-axis (u0 = 0.5), where the footprint is anisotropic */
    RayDiff cam2 = cam;
    cam2.o = (V3){ 0.5f, 0, 0 };
    double u0 = 0.5;
    double yp = 4.0 - sqrt(1.0 - u0 * u0);
    V3 oP = { 0.5f, (float)yp, 0 };
    /* Hit-side shading normal (faces the incoming ray); the outward
       radial normal oOut below is passed separately to fix the dN sign. */
    V3 oN = vnorm((V3){ 0.5f, (float)(4.0 - yp), 0 });
    V3 opx, opy;
    rd_transfer(&cam2, oN, (float)yp, 1e9f, &opx, &opy, &clamp_flag);
    RayDiff orefl;
    V3 oOut = vnorm((V3){ 0.5f, (float)(yp - 4.0), 0 }); /* centre -> P */
    rd_reflect(oP, oN, oOut, R, cam2.d, opx, opy, &orefl);
    V3 oqx, oqy;
    float ot = fabsf((2.0f - oP.y) / orefl.d.y);
    rd_transfer(&orefl, rN, ot, 1e9f, &oqx, &oqy, &clamp_flag);
    float w_off;
    rd_measure(oqx, oqy, &w_off, &ax, &ay);
    double wc = 0.0, wf = 0.0;
    for (int pass = 0; pass < 2; pass++) {
        double h = pass == 0 ? delta : delta / 2.0;
        D3 qxp = ground_truth_Q(u0 + h, 0, R);
        D3 qxm = ground_truth_Q(u0 - h, 0, R);
        D3 qzp = ground_truth_Q(u0,  h, R);
        D3 qzm = ground_truth_Q(u0, -h, R);
        D3 gx = dmul(dsub(qxp, qxm), 1.0 / (2 * h));
        D3 gz = dmul(dsub(qzp, qzm), 1.0 / (2 * h));
        double w = 0.5 * (dlen(gx) + dlen(gz));
        if (pass == 0) wc = w; else wf = w;
    }
    double bt2 = 4.0 / 3.0 * fabs(wc - wf);
    double bf2 = 64.0 * 1.2e-7 * (double)w_off + 1e-6;
    printf("offaxis analytic=%.8f fd_fine=%.8f bound=%.3e\n",
           w_off, wf, bt2 + bf2);
    check("offaxis_error_within_bound",
          fabs((double)w_off - wf) <= bt2 + bf2);

    /* ============ 4. naive distance scaling is wrong (the blur bug) ==== */
    /* Pinhole camera at (0,-3,0) looking at ground plane y=0 through
       screen y=-2. Exact footprint at the ground is 3 world units/pixel;
       the naive |dD|*t projection under-reports it off-axis. */
    RayDiff pin;
    memset(&pin, 0, sizeof(pin));
    pin.o = (V3){ 0, -3, 0 };
    V3 screen = { 0.5f, -2.0f, 0 };
    pin.d = vnorm(vsub(screen, pin.o));
    /* d(direction)/d(screen x) for d = (u,1,0)/sqrt(1+u^2), u=0.5:
         dd/du = (1,-u,0)/(1+u^2)^{3/2}.
       The x factor is 1/(1+u^2)^{3/2}, not 1/(1+u^2): getting it wrong
       leaves a derivative not orthogonal to d, which makes the tangent
       projection miss and the naive |dD|*t scaling diverge from truth.
       Along screen z, d = (u,1,w)/sqrt(1+u^2+w^2) gives at w=0:
         dd/dw = (0,0,1)/sqrt(1+u^2) -- no y component. */
    float u = 0.5f, u2 = 1.0f + u * u;
    float kx = 1.0f / (u2 * sqrtf(u2));
    float ky = -u * kx;
    pin.ddx = (V3){ kx, ky, 0 };
    pin.ddy = (V3){ 0, 0, 1.0f / sqrtf(u2) };
    float tg = 3.0f / pin.d.y;
    V3 gN = { 0, 1, 0 };
    V3 gpx, gpy;
    rd_transfer(&pin, gN, tg, 1e9f, &gpx, &gpy, &clamp_flag);
    float w_ground;
    rd_measure(gpx, gpy, &w_ground, &ax, &ay);
    float naive = 0.5f * (vlen(pin.ddx) + vlen(pin.ddy)) * tg;
    printf("ground exact=%.6f naive_scaling=%.6f ratio=%.3f\n",
           w_ground, naive, naive / w_ground);
    check("ground_exact3", fabsf(w_ground - 3.0f) < 2e-5f);
    check("naive_underestimates", naive < w_ground * 0.97f);

    /* ============ 5. extreme curvature: report and truncate ============ */
    V3 cx, cy; int was_clamped = 0;
    rc = rd_transfer(&refl_sph, rN, 1.0f, 0.5f,
                     &cx, &cy, &was_clamped);
    check("clamp_reported", rc == RD_WARN_CLAMP && was_clamped == 1);
    float w_cap;
    rd_measure(cx, cy, &w_cap, &ax, &ay);
    check("clamp_width_capped", fabsf(w_cap - 0.5f) < 1e-6f);
    /* ray direction must survive clamping untouched */
    check("clamp_preserves_dir",
          fabsf(refl_sph.d.x) < 1e-6f && fabsf(refl_sph.d.y + 1.0f) < 1e-6f);
    rc = rd_transfer(&refl_sph, rN, 1.0f, 0.0f,
                     &cx, &cy, &was_clamped);
    check("zero_cap_rejected", rc == RD_ERR_PARAM);

    /* ============ 6. negative / non-finite differentials rejected ===== */
    check("neg_guard_branch", neg_len(-1.0f) && !neg_len(0.0f));
    RayDiff bad = cam;
    bad.dx = (V3){ NAN, 0, 0 };
    rc = rd_transfer(&bad, sN, 1.0f, 1e9f, &cx, &cy, &was_clamped);
    check("nan_rejected", rc == RD_ERR_NAN);
    bad = cam; bad.ddy = (V3){ INFINITY, 0, 0 };
    rc = rd_transfer(&bad, sN, 1.0f, 1e9f, &cx, &cy, &was_clamped);
    check("inf_rejected", rc == RD_ERR_NAN);
    rc = rd_transfer(&cam, sN, -1.0f, 1e9f, &cx, &cy, &was_clamped);
    check("negative_t_rejected", rc == RD_ERR_PARAM);
    rc = rd_reflect(sP, sN, sOut, 0.0f, cam.d, spx, spy, &bad);
    check("zero_radius_rejected", rc == RD_ERR_PARAM);
    V3 zero3 = { 0, 0, 0 };
    float fu, fv;
    rc = rd_uv_widths(sqx, sqy, zero3, zero3, &fu, &fv);
    check("degenerate_basis_rejected", rc == RD_ERR_DEGEN);

    /* ============ 7. texture-space filter widths ======================= */
    /* receiver textured at 2 world units per texel in x and z */
    V3 e1 = { 2, 0, 0 }, e2 = { 0, 0, 2 };
    rc = rd_uv_widths(sqx, sqy, e1, e2, &fu, &fv);
    check("uv_ok", rc == RD_OK);
    printf("uv_filter_widths=%.6f %.6f (expect 0.75 0.75)\n", fu, fv);
    check("uv_width_correct",
          fabsf(fu - 0.75f) < 2e-5f && fabsf(fv - 0.75f) < 2e-5f);
    rd_uv_widths(qpx, qpy, e1, e2, &fu, &fv);
    check("uv_plane_correct", fabsf(fu - 0.25f) < 2e-5f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
