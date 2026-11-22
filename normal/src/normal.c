
/* normal: area-and-angle weighted vertex-normal reconstruction.
   Per-face geometry normals make a mesh look faceted. Each face normal is
   accumulated at its three corners with weight area * corner_angle, binned
   by welded position, then normalized. The result is compared against the
   supplied shading normals (with separate seam error) and against analytic
   normals (with a provable upper bound). Degenerate faces are listed and
   never participate in the average. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_VERT 128
#define MAX_FACE 192
#define EPS_DEGEN 1e-10f   /* cross-product length below which a face is degenerate */
#define EPS_WELD 1e-7f     /* position distance for welding two vertices */
/* float32 slop for angle compares: acos near 1 has a sqrt singularity, so a
   re-normalized unit vector compared with itself already differs by ~3e-4 rad */
#define EPS_ANG 1e-3f
#define PI 3.14159265358979323846f

typedef struct { float x, y, z; } V3;

static V3 vsub(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static V3 vadd(V3 a, V3 b) { V3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static V3 vscale(V3 a, float s) { V3 r = { a.x * s, a.y * s, a.z * s }; return r; }
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vcross(V3 a, V3 b) {
    V3 r = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return r;
}
static float vlen(V3 a) { return sqrtf(vdot(a, a)); }
static V3 vnorm(V3 a) {
    float l = vlen(a);
    V3 r = { 0, 0, 0 };
    if (l > 0.0f) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}
/* angle between unit vectors, radians */
static float vangle(V3 a, V3 b) {
    float d = vdot(a, b);
    if (d > 1.0f) d = 1.0f;
    if (d < -1.0f) d = -1.0f;
    return acosf(d);
}
/* identity check on vectors: acos has a sqrt singularity at 1, so even a
   vector compared with itself yields ~3e-4 rad in float32. Directional
   equality is checked on the components instead. */
static int vclose(V3 a, V3 b, float tol) {
    return fabsf(a.x - b.x) < tol && fabsf(a.y - b.y) < tol &&
           fabsf(a.z - b.z) < tol;
}

typedef struct { int v[3]; } Tri;
typedef struct {
    V3 pos[MAX_VERT];   /* positions */
    V3 sn[MAX_VERT];    /* supplied shading normals, one per vertex index */
    Tri face[MAX_FACE];
    int nv, nf;
} Mesh;

static int add_vert(Mesh* m, V3 p, V3 shade_n) {
    if (m->nv >= MAX_VERT) return -1;
    m->pos[m->nv] = p;
    m->sn[m->nv] = shade_n;
    return m->nv++;
}
static int add_tri(Mesh* m, int a, int b, int c) {
    if (m->nf >= MAX_FACE) return -1;
    m->face[m->nf].v[0] = a;
    m->face[m->nf].v[1] = b;
    m->face[m->nf].v[2] = c;
    return m->nf++;
}

/* one face's geometry. returns 0 for a degenerate face. */
static int face_geom(const Mesh* m, int f, V3* fn, float* area) {
    const Tri* t = &m->face[f];
    V3 e1 = vsub(m->pos[t->v[1]], m->pos[t->v[0]]);
    V3 e2 = vsub(m->pos[t->v[2]], m->pos[t->v[0]]);
    V3 c = vcross(e1, e2);
    float l = vlen(c);
    if (l < EPS_DEGEN) return 0;
    *fn = vscale(c, 1.0f / l);
    *area = 0.5f * l;
    return 1;
}

/* interior angle of face f at its corner k, radians */
static float corner_angle(const Mesh* m, int f, int k) {
    const Tri* t = &m->face[f];
    V3 e1 = vnorm(vsub(m->pos[t->v[k]], m->pos[t->v[(k + 1) % 3]]));
    V3 e2 = vnorm(vsub(m->pos[t->v[k]], m->pos[t->v[(k + 2) % 3]]));
    return vangle(e1, e2);
}

/* weld vertex indices by position; weld[i] is the position id of vertex i */
static void weld_build(const Mesh* m, int* weld, int* nweld) {
    *nweld = 0;
    for (int i = 0; i < m->nv; i++) {
        int id = -1;
        for (int k = 0; k < i; k++) {
            if (vlen(vsub(m->pos[i], m->pos[k])) < EPS_WELD) { id = weld[k]; break; }
        }
        if (id < 0) id = (*nweld)++;
        weld[i] = id;
    }
}

static int weld_members(const int* weld, int nv, int w, int* out) {
    int k = 0;
    for (int i = 0; i < nv; i++) if (weld[i] == w) out[k++] = i;
    return k;
}

/* reconstruct one unit normal per welded position.
   weight of face f at corner k = area(f) * angle(f, k).
   degenerate face ids go into deg[], excluded from every accumulation.
   welds touched only by degenerate faces keep a zero (undefined) normal. */
static void reconstruct(const Mesh* m, int* weld, V3* recon, int* deg, int* ndeg) {
    int nw = 0;
    weld_build(m, weld, &nw);
    *ndeg = 0;
    V3 acc[MAX_VERT];
    memset(acc, 0, sizeof(acc));
    for (int f = 0; f < m->nf; f++) {
        V3 fn;
        float area;
        if (!face_geom(m, f, &fn, &area)) { deg[(*ndeg)++] = f; continue; }
        for (int k = 0; k < 3; k++) {
            float theta = corner_angle(m, f, k);
            int vi = m->face[f].v[k];
            acc[weld[vi]] = vadd(acc[weld[vi]], vscale(fn, area * theta));
        }
    }
    for (int w = 0; w < nw; w++) recon[w] = vnorm(acc[w]);
}

typedef struct {
    int ndeg;
    int deg[MAX_FACE];
    int nweld, nseam, undefined;
    int seam_identity;      /* structural: seam copies share one weld id => one normal */
    float shade_err_max;     /* max angle vs shading normal, all vertices */
    float seam_shade_err;    /* max angle vs shading normal at seam vertices */
    float seam_shade_jump;   /* max spread between shading normals in one weld */
    float seam_recon_jump;   /* numeric spread between reconstructed normals in one weld */
    float analytic_err_max;  /* max beta: angle(recon, analytic) */
    float analytic_bound_max;/* max alpha: max angle(face normal, analytic) per weld */
} Report;

static void analyze(const Mesh* m, const V3* analytic, Report* r) {
    int weld[MAX_VERT];
    V3 recon[MAX_VERT];
    memset(r, 0, sizeof(*r));
    reconstruct(m, weld, recon, r->deg, &r->ndeg);

    int members[MAX_VERT];
    r->nweld = 0;
    for (int i = 0; i < m->nv; i++) if (weld[i] + 1 > r->nweld) r->nweld = weld[i] + 1;

    r->shade_err_max = 0.0f;
    r->seam_shade_err = 0.0f;
    r->seam_shade_jump = 0.0f;
    r->seam_recon_jump = 0.0f;
    r->seam_identity = 1;
    r->analytic_err_max = 0.0f;
    r->analytic_bound_max = 0.0f;

    for (int w = 0; w < r->nweld; w++) {
        int k = weld_members(weld, m->nv, w, members);
        int seam = k > 1;
        if (seam) r->nseam++;
        int defined = vlen(recon[w]) > 0.5f;
        if (!defined) { r->undefined += k; continue; }

        for (int a = 0; a < k; a++) {
            float e = vangle(recon[w], m->sn[members[a]]);
            if (e > r->shade_err_max) r->shade_err_max = e;
            if (seam && e > r->seam_shade_err) r->seam_shade_err = e;
            /* structural zero-jump guarantee: both copies index the same weld,
               so they literally return the same reconstructed normal */
            if (seam && weld[members[a]] != w) r->seam_identity = 0;
            for (int b = a + 1; b < k; b++) {
                float js = vangle(m->sn[members[a]], m->sn[members[b]]);
                if (js > r->seam_shade_jump) r->seam_shade_jump = js;
                /* numeric cross-check of the same fact; float noise only */
                float jr = vangle(recon[weld[members[a]]], recon[weld[members[b]]]);
                if (jr > r->seam_recon_jump) r->seam_recon_jump = jr;
            }
        }

        if (analytic) {
            float beta = vangle(recon[w], analytic[members[0]]);
            if (beta > r->analytic_err_max) r->analytic_err_max = beta;
            float alpha = 0.0f;
            for (int f = 0; f < m->nf; f++) {
                V3 fn;
                float area;
                if (!face_geom(m, f, &fn, &area)) continue;
                for (int c = 0; c < 3; c++) {
                    if (weld[m->face[f].v[c]] != w) continue;
                    float d = vangle(fn, analytic[m->face[f].v[c]]);
                    if (d > alpha) alpha = d;
                }
            }
            if (alpha > r->analytic_bound_max) r->analytic_bound_max = alpha;
        }
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* single face: normal direction, area, degenerate flag */
static void test_face_geometry(void) {
    Mesh m;
    memset(&m, 0, sizeof(m));
    V3 z = { 0, 0, 1 };
    int v0 = add_vert(&m, (V3){ 0, 0, 0 }, z);
    int v1 = add_vert(&m, (V3){ 1, 0, 0 }, z);
    int v2 = add_vert(&m, (V3){ 0, 1, 0 }, z);
    add_tri(&m, v0, v1, v2);
    V3 fn;
    float area = 0.0f;
    check("face_valid", face_geom(&m, 0, &fn, &area) == 1);
    check("face_area_half", fabsf(area - 0.5f) < 1e-6f);
    check("face_normal_pos_z", fabsf(fn.x) < 1e-7f && fabsf(fn.y) < 1e-7f &&
                               fabsf(fn.z - 1.0f) < 1e-7f);

    Mesh d;
    memset(&d, 0, sizeof(d));
    add_vert(&d, (V3){ 0, 0, 0 }, z);
    add_vert(&d, (V3){ 1, 0, 0 }, z);
    add_vert(&d, (V3){ 2, 0, 0 }, z); /* collinear */
    add_tri(&d, 0, 1, 2);
    check("face_degenerate", face_geom(&d, 0, &fn, &area) == 0);
}

/* hand-computed accumulation: two faces sharing vertex 0,
   T1 area 0.5 angle pi/2 normal +z, T2 area 0.5 angle pi/4 normal -y,
   sum is proportional to (0,-1,2), normalized (0,-1,2)/sqrt(5) */
static void test_weights(void) {
    Mesh m;
    memset(&m, 0, sizeof(m));
    V3 z = { 0, 0, 1 };
    add_vert(&m, (V3){ 0, 0, 0 }, z);   /* 0: shared vertex */
    add_vert(&m, (V3){ 1, 0, 0 }, z);   /* 1 */
    add_vert(&m, (V3){ 0, 1, 0 }, z);   /* 2 */
    add_vert(&m, (V3){ 1, 0, 0 }, z);   /* 3 */
    add_vert(&m, (V3){ 1, 0, 1 }, z);   /* 4 */
    add_tri(&m, 0, 1, 2);
    add_tri(&m, 0, 3, 4);

    int weld[MAX_VERT];
    V3 recon[MAX_VERT];
    int deg[MAX_FACE], ndeg = 0;
    reconstruct(&m, weld, recon, deg, &ndeg);
    V3 expect = vnorm((V3){ 0, -1, 2 });
    V3 got = recon[weld[0]];
    check("weight_area_angle", vclose(got, expect, 1e-6f));
    check("weight_no_degenerate", ndeg == 0);
}

/* good quad plus three degenerate faces, plus an orphan region whose vertices
   only a degenerate face references */
static void test_degenerate(void) {
    Mesh m;
    memset(&m, 0, sizeof(m));
    V3 z = { 0, 0, 1 };
    int p0 = add_vert(&m, (V3){ 0, 0, 0 }, z);
    int p1 = add_vert(&m, (V3){ 1, 0, 0 }, z);
    int p2 = add_vert(&m, (V3){ 1, 1, 0 }, z);
    int p3 = add_vert(&m, (V3){ 0, 1, 0 }, z);
    add_tri(&m, p0, p1, p2);
    add_tri(&m, p0, p2, p3);
    add_tri(&m, p0, p1, p1);                 /* repeated vertex */
    add_tri(&m, p0, p1, add_vert(&m, (V3){ 2, 0, 0 }, z)); /* collinear, adds vertex 6 */
    int o0 = add_vert(&m, (V3){ 5, 5, 5 }, z);
    int o1 = add_vert(&m, (V3){ 6, 5, 5 }, z);
    int o2 = add_vert(&m, (V3){ 7, 5, 5 }, z);
    add_tri(&m, o0, o1, o2);                 /* degenerate-only orphan region */

    Report r;
    analyze(&m, NULL, &r);
    printf("degenerate_faces=[");
    for (int i = 0; i < r.ndeg; i++) printf("%d%s", r.deg[i], i + 1 < r.ndeg ? "," : "");
    printf("]\n");
    check("degen_count", r.ndeg == 3);
    check("degen_ids", r.deg[0] == 2 && r.deg[1] == 3 && r.deg[2] == 4);
    /* four welds are touched only by degenerate faces: collinear vertex (2,0,0)
       plus the three orphan vertices at y=5 */
    check("degen_orphan_undefined", r.undefined == 4);

    /* the same mesh with degenerate faces removed must reconstruct identically */
    Mesh c;
    memset(&c, 0, sizeof(c));
    for (int i = 0; i < 4; i++) add_vert(&c, m.pos[i], m.sn[i]);
    add_tri(&c, 0, 1, 2);
    add_tri(&c, 0, 2, 3);
    int wm[MAX_VERT], wc[MAX_VERT];
    V3 nm[MAX_VERT], nc[MAX_VERT];
    int dm[MAX_FACE], dc[MAX_FACE], nm_deg = 0, nc_deg = 0;
    reconstruct(&m, wm, nm, dm, &nm_deg);
    reconstruct(&c, wc, nc, dc, &nc_deg);
    int same = 1;
    for (int i = 0; i < c.nv && same; i++) {
        if (vlen(vsub(nm[wm[i]], nc[wc[i]])) > 1e-7f) same = 0;
    }
    check("degen_excluded_from_average", same);
    check("clean_has_no_degenerate", nc_deg == 0);
}

/* tessellated cylinder side. analytic normal at every vertex is radial.
   irregular spacing makes the reconstructed normal miss the exact radial
   direction, which exercises the error bound instead of trivial equality. */
static void build_cylinder(Mesh* m, int n, const float* delta_deg) {
    memset(m, 0, sizeof(*m));
    float phi[MAX_VERT];
    phi[0] = 0.0f;
    for (int k = 0; k < n; k++) phi[k + 1] = phi[k] + delta_deg[k] * PI / 180.0f;
    phi[n] = 0.0f; /* close the circle: last column is a welded duplicate of the first */
    for (int k = 0; k <= n; k++) {
        V3 radial = { cosf(phi[k]), sinf(phi[k]), 0.0f };
        add_vert(m, (V3){ radial.x, radial.y, 0.0f }, radial);
        add_vert(m, (V3){ radial.x, radial.y, 1.0f }, radial);
    }
    for (int k = 0; k < n; k++) {
        int b0 = 2 * k, t0 = 2 * k + 1, b1 = 2 * (k + 1), t1 = 2 * (k + 1) + 1;
        add_tri(m, b0, b1, t0);
        add_tri(m, b1, t1, t0);
    }
}

static void test_cylinder(int irregular) {
    int n = 8;
    float regular[8], uneven[8] = { 20, 30, 40, 60, 50, 35, 45, 75 };
    for (int k = 0; k < n; k++) regular[k] = 360.0f / n;
    Mesh m;
    build_cylinder(&m, n, irregular ? uneven : regular);

    Report r;
    analyze(&m, m.sn, &r); /* shading normals equal the analytic radial normals */
    const char* tag = irregular ? "cyl_irregular" : "cyl_regular";
    printf("%s: analytic_err_deg=%.5f bound_deg=%.5f shade_err_deg=%.5f "
           "seams=%d recon_jump_deg=%.6f shade_jump_deg=%.6f\n",
           tag, r.analytic_err_max * 180.0f / PI, r.analytic_bound_max * 180.0f / PI,
           r.shade_err_max * 180.0f / PI, r.nseam,
           r.seam_recon_jump * 180.0f / PI, r.seam_shade_jump * 180.0f / PI);

    check("cyl_weld_seams", r.nseam == 2);
    check("cyl_recon_unit", r.undefined == 0);
    check("cyl_seam_identity", r.seam_identity == 1);
    check("cyl_seam_recon_jump_zero", r.seam_recon_jump < EPS_ANG);
    if (!irregular) {
        /* symmetric umbrella: area-angle average is radial up to float noise */
        check("cyl_regular_analytic_exact", r.analytic_err_max < EPS_ANG);
        check("cyl_regular_shade_exact", r.shade_err_max < EPS_ANG);
        check("cyl_wrap_smooth_shade_jump", r.seam_shade_jump < EPS_ANG);
    } else {
        /* uneven spacing must actually cost something, or the bound check is vacuous */
        check("cyl_irregular_error_present", r.analytic_err_max > 1e-4f);
    }
    /* upper bound: beta <= alpha for every weld (convex cone of unit vectors) */
    check("cyl_analytic_bound_holds",
          r.analytic_err_max <= r.analytic_bound_max + EPS_ANG);
}

/* cube exported the way faceted formats do it: every face owns four vertices
   and a flat shading normal, so each corner position is a three-way seam. */
static void test_cube(void) {
    Mesh m;
    memset(&m, 0, sizeof(m));
    V3 axes[6] = {
        { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 },
        { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 }
    };
    for (int s = 0; s < 6; s++) {
        V3 n = axes[s];
        V3 u = (fabsf(n.x) > 0.5f) ? (V3){ 0, 1, 0 } : (V3){ 1, 0, 0 };
        V3 v = vcross(n, u); /* u x v = n */
        V3 c = n;
        V3 q[4] = {
            vsub(vsub(c, u), v), vsub(vadd(c, u), v),
            vadd(vadd(c, u), v), vsub(vadd(c, v), u)
        };
        int base = m.nv;
        for (int i = 0; i < 4; i++) add_vert(&m, q[i], n);
        add_tri(&m, base, base + 1, base + 2);
        add_tri(&m, base, base + 2, base + 3);
    }

    Report r;
    analyze(&m, NULL, &r);
    printf("cube: seams=%d seam_shade_err_deg=%.5f shade_jump_deg=%.5f "
           "recon_jump_deg=%.6f\n",
           r.nseam, r.seam_shade_err * 180.0f / PI,
           r.seam_shade_jump * 180.0f / PI, r.seam_recon_jump * 180.0f / PI);

    check("cube_corner_seams", r.nseam == 8);
    check("cube_seam_identity", r.seam_identity == 1);
    check("cube_recon_jump_zero", r.seam_recon_jump < EPS_ANG);
    /* three orthogonal face normals, equal weights -> body diagonal */
    float expect_err = acosf(1.0f / sqrtf(3.0f));
    check("cube_seam_shade_error", fabsf(r.seam_shade_err - expect_err) < 1e-5f);
    check("cube_shade_jump_right_angle", fabsf(r.seam_shade_jump - PI / 2.0f) < 1e-5f);

    int weld[MAX_VERT];
    V3 recon[MAX_VERT];
    int deg[MAX_FACE], ndeg = 0;
    reconstruct(&m, weld, recon, deg, &ndeg);
    int diagonal = 1;
    for (int i = 0; i < m.nv; i++) {
        if (!vclose(recon[weld[i]], vnorm(m.pos[i]), 1e-6f)) diagonal = 0;
    }
    check("cube_recon_body_diagonal", diagonal);
    check("cube_no_degenerate", ndeg == 0);
}

/* direct stress test of the bound theorem, independent of any mesh:
   take unit vectors inside a cone of half-angle alpha around an axis,
   combine with arbitrary non-negative weights, and require the normalized
   result to stay within alpha of the axis. Runs many deterministic trials
   with skewed weights so a vacuous (equal-weight-only) check cannot pass. */
static unsigned lcg_state = 0x12345678u;
static float urand(void) {
    lcg_state = lcg_state * 1103515245u + 12345u;
    return (float)(lcg_state >> 8) / (float)(1u << 24);
}

static void test_bound_random(void) {
    int trials = 2000, violations = 0;
    float worst_margin = 0.0f;  /* largest (beta - alpha) seen; must stay <= slop */
    float best_margin = 1e30f;  /* smallest (alpha - beta); bound must be tight somewhere */
    for (int t = 0; t < trials; t++) {
        V3 axis = vnorm((V3){ urand() - 0.5f, urand() - 0.5f, urand() - 0.5f });
        V3 helper = (fabsf(axis.z) < 0.9f) ? (V3){ 0, 0, 1 } : (V3){ 1, 0, 0 };
        V3 t1 = vnorm(vcross(axis, helper));
        V3 t2 = vcross(axis, t1);
        float alpha = 0.05f + urand() * 1.2f; /* cone half-angle, radians */
        int cnt = 2 + (int)(urand() * 7.0f); /* 2..8 vectors */
        V3 acc = { 0, 0, 0 };
        for (int i = 0; i < cnt; i++) {
            float theta = urand() * alpha;
            float phi = urand() * 2.0f * PI;
            V3 dir = vadd(vscale(axis, cosf(theta)),
                      vadd(vscale(t1, sinf(theta) * cosf(phi)),
                           vscale(t2, sinf(theta) * sinf(phi))));
            float w = 0.01f + urand() * urand() * 100.0f; /* skewed weights */
            acc = vadd(acc, vscale(dir, w));
        }
        V3 mean = vnorm(acc);
        float beta = vangle(mean, axis);
        if (beta > alpha + EPS_ANG) violations++;
        if (beta - alpha > worst_margin) worst_margin = beta - alpha;
        if (alpha - beta < best_margin) best_margin = alpha - beta;
    }
    printf("bound_random: trials=%d violations=%d worst_margin_rad=%.2e "
           "best_margin_rad=%.2e\n", trials, violations, worst_margin, best_margin);
    check("bound_random_holds", violations == 0);
    /* some trial must approach the bound within a few degrees; otherwise the
       cone construction is wrong and the bound check is vacuous */
    check("bound_random_tight", best_margin < 0.05f);
}

int main(void) {
    test_face_geometry();
    test_weights();
    test_degenerate();
    test_cylinder(0);
    test_cylinder(1);
    test_cube();
    test_bound_random();

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
