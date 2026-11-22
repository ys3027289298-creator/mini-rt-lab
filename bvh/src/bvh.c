
#define _POSIX_C_SOURCE 199309L
/* bvh: median-split AABB tree over triangle bounds, nearest-hit traversal,
 * incremental refit for animated scenes, rebuild criterion, and audit.
 *
 * Animated-scene path:
 *   - build() fixes the leaf partition.
 *   - a moved prim only affects its leaf and that leaf's ancestors; refit()
 *     recomputes exactly that union, bottom-up.
 *   - with the prim set fixed, refit is exact, so a rebuild is never needed
 *     for correctness; rebuild is only needed when the partition itself is
 *     obsolete (prim set changed) or when traversal quality has drifted
 *     (SAH proxy). See README for the derived criterion.
 *   - audit() checks the containment invariant that makes traversal culling
 *     sound; raycast vs. brute force catches false culls produced by a stale
 *     box, e.g. a broken dirty set or a skipped ancestor. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_PRIM 4096
#define MAX_NODE (2 * MAX_PRIM)
#define MAX_STACK 128

#define CONTAINS_EPS 1e-5f
#define HIT_EPS 1e-4f
#define REBUILD_RATIO 1.5

typedef struct { float min[3], max[3]; } AABB;
typedef struct { AABB box; float cx, cy, cz; int id; float v[3][3]; int is_tri; } Prim;
typedef struct { AABB box; int parent, left, right; int first, count; } Node;

typedef struct {
    Prim p[MAX_PRIM];
    Node n[MAX_NODE];
    int nprim, nnode, root;
    int built_nprim;        /* nprim at the last build: topology-change test */
    int slot[MAX_PRIM];     /* prim id -> index in p[] (p is reordered by build) */
    int leaf_of[MAX_PRIM];  /* prim id -> leaf node */
    unsigned mark[MAX_NODE];
    unsigned stamp;
    double area_sum;        /* sum of internal node AABBs (SAH proxy numerator) */
    float q_build;          /* area_sum / root area at the last build */
    int nleaves, ninternal;
    int rebuilds, refits;
} Tree;

static void aabb_empty(AABB* b) {
    for (int i = 0; i < 3; i++) { b->min[i] = 1e30f; b->max[i] = -1e30f; }
}
static void aabb_grow(AABB* b, const AABB* o) {
    for (int i = 0; i < 3; i++) {
        if (o->min[i] < b->min[i]) b->min[i] = o->min[i];
        if (o->max[i] > b->max[i]) b->max[i] = o->max[i];
    }
}
static float aabb_area(const AABB* b) {
    float dx = b->max[0] - b->min[0];
    float dy = b->max[1] - b->min[1];
    float dz = b->max[2] - b->min[2];
    if (dx < 0 || dy < 0 || dz < 0) return 0.0f;
    return 2.0f * (dx * dy + dy * dz + dz * dx);
}
static int aabb_hit(const AABB* b, const float o[3], const float d[3], float tmax, float* tmin) {
    float t0 = 0.0f, t1 = tmax;
    for (int i = 0; i < 3; i++) {
        if (fabsf(d[i]) < 1e-8f) {
            if (o[i] < b->min[i] || o[i] > b->max[i]) return 0;
        } else {
            float inv = 1.0f / d[i];
            float ta = (b->min[i] - o[i]) * inv;
            float tb = (b->max[i] - o[i]) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    if (tmin) *tmin = t0;
    return 1;
}
static int aabb_contains(const AABB* outer, const AABB* inner, float eps) {
    for (int i = 0; i < 3; i++) {
        if (inner->min[i] < outer->min[i] - eps) return 0;
        if (inner->max[i] > outer->max[i] + eps) return 0;
    }
    return 1;
}

static int tree_add(Tree* t, float minx, float miny, float minz, float maxx, float maxy, float maxz) {
    if (t->nprim >= MAX_PRIM) return -1;
    Prim* p = &t->p[t->nprim];
    p->box.min[0] = minx; p->box.min[1] = miny; p->box.min[2] = minz;
    p->box.max[0] = maxx; p->box.max[1] = maxy; p->box.max[2] = maxz;
    p->cx = (minx + maxx) * 0.5f; p->cy = (miny + maxy) * 0.5f; p->cz = (minz + maxz) * 0.5f;
    p->id = t->nprim;
    p->is_tri = 0;
    memset(p->v, 0, sizeof(p->v));
    t->slot[t->nprim] = t->nprim;
    t->leaf_of[t->nprim] = -1;
    return t->nprim++;
}

static int add_tri(Tree* t, const float v[3][3]) {
    float mn[3] = { 1e30f, 1e30f, 1e30f }, mx[3] = { -1e30f, -1e30f, -1e30f };
    for (int k = 0; k < 3; k++) {
        for (int i = 0; i < 3; i++) {
            if (v[k][i] < mn[i]) mn[i] = v[k][i];
            if (v[k][i] > mx[i]) mx[i] = v[k][i];
        }
    }
    int id = tree_add(t, mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]);
    if (id < 0) return -1;
    Prim* p = &t->p[t->slot[id]];
    memcpy(p->v, v, sizeof(p->v));
    p->is_tri = 1;
    return id;
}

/* Move prim `id` to new vertices. Keeps box/centroid in sync. O(cost of one prim). */
static void apply_tri(Tree* t, int id, const float v[3][3]) {
    if (id < 0 || id >= t->nprim) return;
    int pos = t->slot[id];
    if (pos < 0 || pos >= t->nprim) return;
    Prim* p = &t->p[pos];
    memcpy(p->v, v, sizeof(p->v));
    aabb_empty(&p->box);
    for (int k = 0; k < 3; k++) {
        AABB vb;
        for (int i = 0; i < 3; i++) { vb.min[i] = v[k][i]; vb.max[i] = v[k][i]; }
        aabb_grow(&p->box, &vb);
    }
    p->cx = (p->box.min[0] + p->box.max[0]) * 0.5f;
    p->cy = (p->box.min[1] + p->box.max[1]) * 0.5f;
    p->cz = (p->box.min[2] + p->box.max[2]) * 0.5f;
}

static int build(Tree* t, int first, int count, int parent) {
    if (t->nnode >= MAX_NODE || count <= 0) return -1;
    int ni = t->nnode++;
    Node* n = &t->n[ni];
    n->parent = parent;
    n->left = n->right = -1;
    n->first = count;   /* provisional: leaf */
    n->count = 0;
    aabb_empty(&n->box);
    for (int i = first; i < first + count; i++) aabb_grow(&n->box, &t->p[i].box);

    if (count <= 2) {
        n->first = first; n->count = count;
        n->left = n->right = -1;
        for (int i = first; i < first + count; i++) {
            int id = t->p[i].id;
            t->slot[id] = i;
            t->leaf_of[id] = ni;
        }
        t->nleaves++;
        return ni;
    }
    /* split on widest axis at box-center */
    float w[3] = { n->box.max[0] - n->box.min[0], n->box.max[1] - n->box.min[1], n->box.max[2] - n->box.min[2] };
    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);
    float c[3] = { n->box.min[0] + w[0] * 0.5f, n->box.min[1] + w[1] * 0.5f, n->box.min[2] + w[2] * 0.5f };
    int mid = first;
    for (int i = first; i < first + count; i++) {
        float v = (axis == 0) ? t->p[i].cx : (axis == 1 ? t->p[i].cy : t->p[i].cz);
        if (v < c[axis]) {
            Prim tmp = t->p[i]; t->p[i] = t->p[mid]; t->p[mid] = tmp;
            mid++;
        }
    }
    if (mid == first || mid == first + count) mid = first + count / 2;
    /* A spatial split can be arbitrarily skewed (e.g. one outlier); bound it so
     * recursion depth and the traversal stack stay logarithmic. */
    if ((mid - first) * 7 < count || (first + count - mid) * 7 < count)
        mid = first + count / 2;
    n->first = -1; n->count = 0;
    int l = build(t, first, mid - first, ni);
    int r = build(t, mid, first + count - mid, ni);
    if (l < 0 || r < 0) return -1;
    n->left = l; n->right = r;
    aabb_empty(&n->box);
    aabb_grow(&n->box, &t->n[l].box);
    aabb_grow(&n->box, &t->n[r].box);
    t->ninternal++;
    t->area_sum += aabb_area(&n->box);
    return ni;
}

static double q_now(const Tree* t) {
    if (t->root < 0) return 0.0;
    float ra = aabb_area(&t->n[t->root].box);
    return (ra > 0.0f) ? t->area_sum / (double)ra : 0.0;
}

static void rebuild(Tree* t) {
    t->rebuilds++;
    t->nnode = 0;
    t->nleaves = t->ninternal = 0;
    t->area_sum = 0.0;
    t->root = build(t, 0, t->nprim, -1);
    t->built_nprim = t->nprim;
    t->q_build = (float)q_now(t);
}

/* ---- criteria ------------------------------------------------------------ */

/* Rebuild criterion, split into "must" and "quality":
 *   REBUILD_TOPOLOGY: the prim set changed. A leaf partition is a partition of
 *     the old prim set; refit only recomputes boxes over existing ranges, so a
 *     newly added prim is in no leaf and a removed one leaves a stale range.
 *     This is the only correctness-mandatory trigger.
 *   REBUILD_SAH: prim set fixed, but motion moved geometry between nodes, so
 *     the old median split no longer matches the layout. The SAH proxy is
 *     sum(internal AABB areas)/root AABB area (scale invariant), compared with
 *     its value at build time. Refit cannot repair the partition; a rebuild
 *     restores it. This is a profitability trigger, not a correctness one.
 * With a fixed prim set refit is *exact*, so hits are correct even when the
 * SAH drift is huge; ignoring the trigger only costs traversal time. */
#define REBUILD_TOPOLOGY 1
#define REBUILD_SAH 2

static int should_rebuild(const Tree* t) {
    int r = 0;
    if (t->nprim != t->built_nprim) r |= REBUILD_TOPOLOGY;
    if (t->q_build > 0.0f && q_now(t) > REBUILD_RATIO * (double)t->q_build) r |= REBUILD_SAH;
    return r;
}

/* ---- incremental refit --------------------------------------------------- */

/* Recompute the box of one dirty leaf (from its prims) and its ancestors.
 * The update set is exactly: dirty leaves + all their ancestors. Any node
 * outside that union has a subtree of untouched prims, so its box is still
 * valid. */
static int refit_fix(Tree* t, int ni, unsigned s, int* updated) {
    Node* n = &t->n[ni];
    if (n->count > 0) {
        if (t->mark[ni] != s) return 0;
        AABB b;
        aabb_empty(&b);
        for (int i = n->first; i < n->first + n->count; i++) aabb_grow(&b, &t->p[i].box);
        n->box = b;
        (*updated)++;
        return 1;
    }
    int ch = 0;
    if (n->left >= 0 && t->mark[n->left] == s) ch |= refit_fix(t, n->left, s, updated);
    if (n->right >= 0 && t->mark[n->right] == s) ch |= refit_fix(t, n->right, s, updated);
    if (ch) {
        t->area_sum -= aabb_area(&n->box);
        aabb_empty(&n->box);
        aabb_grow(&n->box, &t->n[n->left].box);
        aabb_grow(&n->box, &t->n[n->right].box);
        t->area_sum += aabb_area(&n->box);
        (*updated)++;
        return 1;
    }
    return 0;
}

/* Refit for prim ids that moved. Returns the number of nodes whose box was
 * recomputed. Duplicate ids are fine. */
static int refit(Tree* t, const int* ids, int nids) {
    t->refits++;
    if (++t->stamp == 0) { memset(t->mark, 0, sizeof(t->mark)); t->stamp = 1; }
    unsigned s = t->stamp;
    for (int k = 0; k < nids; k++) {
        int id = ids[k];
        if (id < 0 || id >= t->nprim) continue;
        int lf = t->leaf_of[id];
        if (lf < 0 || lf >= t->nnode) continue;
        if (t->mark[lf] != s) {
            t->mark[lf] = s;
            int p = t->n[lf].parent;
            while (p >= 0 && t->mark[p] != s) { t->mark[p] = s; p = t->n[p].parent; }
        }
    }
    int updated = 0;
    if (t->root >= 0) refit_fix(t, t->root, s, &updated);
    return updated;
}

/* ---- containment audit --------------------------------------------------- */

/* The invariant traversal relies on: every node box contains the boxes of all
 * prims in its subtree. Returns the number of violating nodes. A stale leaf or
 * a skipped ancestor shows up here immediately. */
static int audit_contains(const Tree* t, float eps) {
    int bad = 0;
    for (int i = 0; i < t->nnode; i++) {
        const Node* n = &t->n[i];
        if (n->count > 0) {
            for (int j = n->first; j < n->first + n->count; j++) {
                if (!aabb_contains(&n->box, &t->p[j].box, eps)) { bad++; break; }
            }
        } else if (n->left >= 0 && n->right >= 0) {
            if (!aabb_contains(&n->box, &t->n[n->left].box, eps) ||
                !aabb_contains(&n->box, &t->n[n->right].box, eps)) bad++;
        }
    }
    return bad;
}

/* ---- intersection -------------------------------------------------------- */

static int tri_hit(const float v[3][3], const float o[3], const float d[3], float tmax, float* tout) {
    float e1[3], e2[3], pv[3], tv[3], qv[3];
    for (int i = 0; i < 3; i++) {
        e1[i] = v[1][i] - v[0][i];
        e2[i] = v[2][i] - v[0][i];
        tv[i] = o[i] - v[0][i];
    }
    pv[0] = d[1] * e2[2] - d[2] * e2[1];
    pv[1] = d[2] * e2[0] - d[0] * e2[2];
    pv[2] = d[0] * e2[1] - d[1] * e2[0];
    float det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
    if (fabsf(det) < 1e-9f) return 0;
    float inv = 1.0f / det;
    float u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
    if (u < 0.0f || u > 1.0f) return 0;
    qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
    qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
    qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
    float w = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) * inv;
    if (w < 0.0f || u + w > 1.0f) return 0;
    float tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
    if (tt < HIT_EPS || tt > tmax) return 0;
    *tout = tt;
    return 1;
}

/* Box tests performed by the last raycast(). Used to show that a refit tree
 * traverses as cheaply as a freshly rebuilt one (same hierarchy, same boxes),
 * i.e. that the speedup comes from avoiding rebuild, not from worse queries. */
static long g_box_tests = 0;

static int raycast(const Tree* t, const float o[3], const float d[3], float tmax,
                   float* tout, int* pid) {
    int stack[MAX_STACK];
    int sp = 0;
    float bestt = tmax;
    int bestid = -1;
    if (t->root < 0) return 0;
    stack[sp++] = t->root;
    while (sp > 0) {
        int ni = stack[--sp];
        const Node* n = &t->n[ni];
        g_box_tests++;
        if (!aabb_hit(&n->box, o, d, bestt, NULL)) continue;
        if (n->count > 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                const Prim* p = &t->p[i];
                if (!p->is_tri) continue;
                float tt;
                if (tri_hit(p->v, o, d, bestt, &tt) && tt < bestt) {
                    bestt = tt; bestid = p->id;
                }
            }
        } else {
            if (n->left >= 0) stack[sp++] = n->left;
            if (n->right >= 0) stack[sp++] = n->right;
        }
    }
    if (bestid >= 0) { if (tout) *tout = bestt; if (pid) *pid = bestid; return 1; }
    return 0;
}

static int brute(const Tree* t, const float o[3], const float d[3], float tmax,
                 float* tout, int* pid) {
    float bestt = tmax;
    int bestid = -1;
    for (int i = 0; i < t->nprim; i++) {
        const Prim* p = &t->p[i];
        if (!p->is_tri) continue;
        float tt;
        if (tri_hit(p->v, o, d, bestt, &tt) && tt < bestt) { bestt = tt; bestid = p->id; }
    }
    if (bestid >= 0) { if (tout) *tout = bestt; if (pid) *pid = bestid; return 1; }
    return 0;
}

/* 1 if tree traversal and brute force disagree on hit/miss, prim id, or depth. */
static int ray_mismatch(const Tree* t, const float o[3], const float d[3]) {
    float t1 = 1e30f, t2 = 1e30f;
    int id1 = -1, id2 = -1;
    int h1 = raycast(t, o, d, 1e30f, &t1, &id1);
    int h2 = brute(t, o, d, 1e30f, &t2, &id2);
    if (h1 != h2) return 1;
    if (!h1) return 0;
    if (id1 != id2) return 1;
    float scale = (t2 > 1.0f) ? t2 : 1.0f;
    if (fabsf(t1 - t2) > 1e-3f * scale) return 1;
    return 0;
}

/* ---- misc ---------------------------------------------------------------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static unsigned rng_state = 12345u;
static unsigned rnd(void) {
    rng_state = rng_state * 1664525u + 1013904223u;
    return rng_state;
}
static float frand(float a, float b) {
    return a + (b - a) * ((float)(rnd() >> 8) / 16777216.0f);
}

/* Build a triangle around center c: random vertices in a unit-ish cube. */
static void make_tri(float v[3][3], const float c[3], float radius) {
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < 3; i++)
            v[k][i] = c[i] + frand(-radius, radius);
}

/* ---- tests --------------------------------------------------------------- */

static Tree g_a, g_b, g_s;  /* large; keep off the stack */

static void test_legacy(void) {
    Tree t;
    memset(&t, 0, sizeof(t));
    tree_add(&t, 1, 0, 0, 2, 1, 1);
    tree_add(&t, 5, 0, 0, 6, 1, 1);
    tree_add(&t, 3, 0, 0, 4, 1, 1);
    check("prim_count", t.nprim == 3);

    int root = build(&t, 0, t.nprim, -1);
    check("build_ok", root == 0 && t.nnode >= 2);

    AABB whole;
    aabb_empty(&whole);
    for (int i = 0; i < t.nprim; i++) aabb_grow(&whole, &t.p[i].box);
    check("root_covers", fabsf(t.n[0].box.min[0] - whole.min[0]) < 1e-5f &&
                         fabsf(t.n[0].box.max[0] - whole.max[0]) < 1e-5f);

    float o[3] = { 0, 0.5f, 0.5f }, d[3] = { 1, 0, 0 };
    float tm = 0;
    check("hit_box", aabb_hit(&t.n[0].box, o, d, 100.0f, &tm) == 1 && tm > 0.0f);

    float d2[3] = { -1, 0, 0 };
    check("miss_box", aabb_hit(&t.n[0].box, o, d2, 100.0f, &tm) == 0);

    float o3[3] = { 0, 5, 5 };
    check("miss_parallel", aabb_hit(&t.n[0].box, o3, d, 100.0f, &tm) == 0);

    check("area_positive", aabb_area(&t.n[0].box) > 0.0f);
    check("audit_after_build", audit_contains(&t, CONTAINS_EPS) == 0);

    Tree big;
    memset(&big, 0, sizeof(big));
    int rc = 0;
    for (int i = 0; i < MAX_PRIM; i++) rc = tree_add(&big, (float)i, 0, 0, (float)i + 0.5f, 1, 1);
    check("capacity", rc == MAX_PRIM - 1 && big.nprim == MAX_PRIM);
    check("build_full", build(&big, 0, big.nprim, -1) == 0);
}

/* Scope of a single-prim refit: must be exactly that prim's leaf plus its
 * ancestors, and every node outside the scope must be untouched. */
static void test_refit_scope(void) {
    const int NP = 256;
    memset(&g_s, 0, sizeof(g_s));
    float base[256][3][3];
    float c[3];
    for (int i = 0; i < NP; i++) {
        c[0] = frand(0, 64); c[1] = frand(0, 64); c[2] = frand(0, 64);
        make_tri(base[i], c, 1.0f);
        add_tri(&g_s, base[i]);
    }
    rebuild(&g_s);
    check("scope_tree_ok", g_s.root == 0 && audit_contains(&g_s, CONTAINS_EPS) == 0);

    int id = 42;
    int lf = g_s.leaf_of[id];
    unsigned char in_scope[MAX_NODE];
    memset(in_scope, 0, sizeof(in_scope));
    int expected = 0;
    for (int ni = lf; ni >= 0; ni = g_s.n[ni].parent) { in_scope[ni] = 1; expected++; }

    AABB before[MAX_NODE];
    for (int i = 0; i < g_s.nnode; i++) before[i] = g_s.n[i].box;

    float moved[3][3];
    memcpy(moved, base[id], sizeof(moved));
    for (int k = 0; k < 3; k++) moved[k][1] += 7.0f;
    apply_tri(&g_s, id, moved);
    int ids[1] = { id };
    int updated = refit(&g_s, ids, 1);

    check("refit_scope_exact", updated == expected);
    int outside_touched = 0;
    for (int i = 0; i < g_s.nnode; i++) {
        if (!in_scope[i] && memcmp(&before[i], &g_s.n[i].box, sizeof(AABB)) != 0) outside_touched = 1;
    }
    check("refit_scope_isolated", outside_touched == 0);
    check("refit_scope_audit", audit_contains(&g_s, CONTAINS_EPS) == 0);
    printf("scope_nodes=%d tree_nodes=%d\n", updated, g_s.nnode);
}

/* Animation: full refit vs. full rebuild each frame.
 * Rigid translation preserves every AABB, so the SAH proxy is unchanged and
 * the criterion must not ask for a rebuild; this is the refit-favourable case.
 * Every frame is checked against brute force and against a tree rebuilt from
 * scratch, so "refit == rebuild" is verified on hits, not assumed. */
static void test_animation_rigid(void) {
    const int NP = 1024;
    const int FRAMES = 96;
    const int NRAY = 256;

    Tree* A = &g_a;
    Tree* B = &g_b;
    memset(A, 0, sizeof(*A));
    memset(B, 0, sizeof(*B));

    static float base[1024][3][3];
    static float ray_dir[256][3], ray_aim[256][3];

    for (int i = 0; i < NP; i++) {
        int cluster = i % 8;
        float c[3] = { (float)(cluster % 4) * 16.0f + 8.0f,
                       (float)(cluster / 4) * 16.0f + 8.0f,
                       (float)(i / 8 % 8) * 4.0f + 4.0f };
        make_tri(base[i], c, 0.6f);
        add_tri(A, base[i]);
    }
    rebuild(A);
    *B = *A;
    rebuild(B);
    check("rigid_audit_build", audit_contains(A, CONTAINS_EPS) == 0);

    /* Rays aimed at a primitive's current centroid: guaranteed hits (targets
     * are precomputed, the aim point follows the per-frame offset), so the
     * comparison is not vacuously true. */
    for (int r = 0; r < NRAY; r++) {
        int target = (int)(rnd() % (unsigned)NP);
        for (int i = 0; i < 3; i++)
            ray_aim[r][i] = (base[target][0][i] + base[target][1][i] + base[target][2][i]) / 3.0f;
        float dir[3] = { frand(-1, 1), frand(-1, 1), frand(-1, 1) };
        float len = sqrtf(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        if (len < 1e-3f) { dir[0] = 0; dir[1] = 1; dir[2] = 0; len = 1; }
        for (int i = 0; i < 3; i++) ray_dir[r][i] = dir[i] / len;
    }

    float moved[3][3];
    double t_refit = 0.0, t_rebuild = 0.0;
    int mismatches = 0, audits_bad = 0, triggers = 0, hit_rays = 0;

    for (int f = 1; f <= FRAMES; f++) {
        float off[3] = { (float)f * 0.05f, (float)((f % 7) - 3) * 0.05f, (float)f * 0.02f };
        for (int i = 0; i < NP; i++) {
            for (int k = 0; k < 3; k++)
                for (int j = 0; j < 3; j++) moved[k][j] = base[i][k][j] + off[j];
            apply_tri(A, i, moved);
            apply_tri(B, i, moved);
        }
        int ids[1024];
        for (int i = 0; i < NP; i++) ids[i] = i;

        double a0 = now_sec();
        refit(A, ids, NP);
        double a1 = now_sec();
        rebuild(B);
        double b1 = now_sec();
        t_refit += a1 - a0;
        t_rebuild += b1 - a1;

        if (should_rebuild(A)) triggers++;
        audits_bad += audit_contains(A, CONTAINS_EPS);

        for (int r = 0; r < NRAY; r++) {
            float ro[3];
            for (int i = 0; i < 3; i++)
                ro[i] = ray_aim[r][i] + off[i] - ray_dir[r][i] * 70.0f;
            float t1 = 0, t2 = 0, t3 = 0;
            int id1 = -1, id2 = -1, id3 = -1;
            int h1 = raycast(A, ro, ray_dir[r], 1e30f, &t1, &id1);
            int h2 = brute(A, ro, ray_dir[r], 1e30f, &t2, &id2);
            int h3 = raycast(B, ro, ray_dir[r], 1e30f, &t3, &id3);
            if (h2) hit_rays++;
            int bad = 0;
            if (h1 != h2 || h3 != h2) bad = 1;
            else if (h1 && (id1 != id2 || id2 != id3)) bad = 1;
            else if (h1) {
                float s1 = t2 > 1.0f ? t2 : 1.0f;
                if (fabsf(t1 - t2) > 1e-3f * s1 || fabsf(t3 - t2) > 1e-3f * s1) bad = 1;
            }
            mismatches += bad;
        }
    }

    check("rigid_mismatch", mismatches == 0);
    check("rigid_audit", audits_bad == 0);
    check("rigid_criterion_quiet", triggers == 0);
    check("rigid_rays_hit", hit_rays > FRAMES * NRAY / 4);
    double refit_us = t_refit / FRAMES * 1e6;
    double rebuild_us = t_rebuild / FRAMES * 1e6;
    printf("rigid_frames=%d prims=%d rays_per_frame=%d hit_rays=%d\n", FRAMES, NP, NRAY, hit_rays);
    printf("refit_us_per_frame=%.2f rebuild_us_per_frame=%.2f speedup=%.2fx\n",
           refit_us, rebuild_us, rebuild_us / refit_us);
    check("refit_cheaper", refit_us < rebuild_us);
    check("refit_speedup_gt2", rebuild_us > 2.0 * refit_us);
}

/* Partial motion, the realistic animation case: a small fraction of prims
 * moves per frame. Refit touches only those leaves and their ancestors; the
 * per-frame rebuild pays for all prims. Same hit check as above. */
static void test_animation_partial(void) {
    const int NP = 4096;
    const int FRAMES = 64;
    const int NRAY = 128;
    const int MOVING = NP * 3 / 100;

    Tree* A = &g_a;
    Tree* B = &g_b;
    memset(A, 0, sizeof(*A));
    memset(B, 0, sizeof(*B));

    static float base[4096][3][3];
    static int tg[128];
    static float rdir[128][3];

    for (int i = 0; i < NP; i++) {
        float c[3] = { frand(0, 256), frand(0, 256), frand(0, 256) };
        make_tri(base[i], c, 0.5f);
        add_tri(A, base[i]);
    }
    rebuild(A);
    *B = *A;
    rebuild(B);
    check("partial_audit_build", audit_contains(A, CONTAINS_EPS) == 0);

    for (int r = 0; r < NRAY; r++) {
        tg[r] = (int)(rnd() % (unsigned)NP);
        float dir[3] = { frand(-1, 1), frand(-1, 1), frand(-1, 1) };
        float len = sqrtf(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        if (len < 1e-3f) { dir[0] = 0; dir[1] = 1; dir[2] = 0; len = 1; }
        for (int i = 0; i < 3; i++) rdir[r][i] = dir[i] / len;
    }

    static float cur[4096][3][3];
    memcpy(cur, base, sizeof(base));
    int ids[4096];
    for (int i = 0; i < NP; i++) ids[i] = i;

    float v[3][3];
    double t_refit = 0.0, t_rebuild = 0.0;
    long boxes_refit = 0, boxes_rebuild = 0;
    int mismatches = 0, audits_bad = 0, triggers = 0, hit_rays = 0, touched = 0;

    for (int f = 1; f <= FRAMES; f++) {
        float off[3] = { (float)(f % 5) * 0.1f, (float)(f % 3) * 0.1f, (float)(f % 7) * 0.1f };
        for (int k = 0; k < MOVING; k++) {
            int i = (int)((unsigned)k * 37u + (unsigned)f * 11u) % NP;
            for (int a = 0; a < 3; a++)
                for (int j = 0; j < 3; j++) v[a][j] = base[i][a][j] + off[j];
            memcpy(cur[i], v, sizeof(v));
            apply_tri(A, i, v);
            apply_tri(B, i, v);
            ids[k] = i;
        }

        double a0 = now_sec();
        touched += refit(A, ids, MOVING);
        double a1 = now_sec();
        rebuild(B);
        double b1 = now_sec();
        t_refit += a1 - a0;
        t_rebuild += b1 - a1;

        if (should_rebuild(A) & REBUILD_TOPOLOGY) triggers++;
        audits_bad += audit_contains(A, CONTAINS_EPS);

        for (int r = 0; r < NRAY; r++) {
            float c[3], ro[3];
            for (int i = 0; i < 3; i++)
                c[i] = (cur[tg[r]][0][i] + cur[tg[r]][1][i] + cur[tg[r]][2][i]) / 3.0f;
            for (int i = 0; i < 3; i++) ro[i] = c[i] - rdir[r][i] * 70.0f;
            float t1 = 0, t2 = 0, t3 = 0;
            int id1 = -1, id2 = -1, id3 = -1;
            long b0 = g_box_tests;
            int h1 = raycast(A, ro, rdir[r], 1e30f, &t1, &id1);
            long b1 = g_box_tests;
            int h2 = brute(A, ro, rdir[r], 1e30f, &t2, &id2);
            int h3 = raycast(B, ro, rdir[r], 1e30f, &t3, &id3);
            long b2 = g_box_tests;
            boxes_refit += b1 - b0;
            boxes_rebuild += b2 - b1;
            if (h2) hit_rays++;
            int bad = (h1 != h2 || h3 != h2);
            if (!bad && h2) {
                bad = (id1 != id2 || id2 != id3);
                float s1 = t2 > 1.0f ? t2 : 1.0f;
                if (!bad && (fabsf(t1 - t2) > 1e-3f * s1 || fabsf(t3 - t2) > 1e-3f * s1)) bad = 1;
            }
            mismatches += bad;
        }
    }

    check("partial_mismatch", mismatches == 0);
    check("partial_audit", audits_bad == 0);
    check("partial_no_topology_trigger", triggers == 0);
    check("partial_rays_hit", hit_rays > FRAMES * NRAY / 4);
    double refit_us = t_refit / FRAMES * 1e6;
    double rebuild_us = t_rebuild / FRAMES * 1e6;
    printf("partial_frames=%d prims=%d moving_per_frame=%d nodes=%d nodes_refit_per_frame=%.0f\n",
           FRAMES, NP, MOVING, A->nnode, (double)touched / FRAMES);
    printf("partial_refit_us=%.2f partial_rebuild_us=%.2f speedup=%.2fx\n",
           refit_us, rebuild_us, rebuild_us / refit_us);
    check("partial_refit_cheaper", refit_us < rebuild_us);
    check("partial_speedup_gt5", rebuild_us > 5.0 * refit_us);
    /* A refit tree must not be a worse query structure than a rebuilt one:
     * identical hierarchy and boxes mean identical traversal work. */
    check("partial_traversal_not_worse", boxes_refit <= boxes_rebuild);
    printf("partial_box_tests_refit=%ld rebuilt=%ld\n", boxes_refit, boxes_rebuild);
}

/* Adversarial layout: one far outlier plus a tight cluster. An unclamped
 * spatial split degrades to a linear chain here; depth must stay logarithmic
 * or the traversal stack (MAX_STACK) overflows. */
static void test_skewed_layout(void) {
    static Tree s;
    memset(&s, 0, sizeof(s));
    float c[3] = { 0, 0, 0 };
    float v[3][3];
    for (int i = 0; i < 1000; i++) {
        c[0] = frand(-1, 1); c[1] = frand(-1, 1); c[2] = frand(-1, 1);
        make_tri(v, c, 0.01f);
        add_tri(&s, v);
    }
    c[0] = 1e4f; c[1] = 1e4f; c[2] = 1e4f;
    make_tri(v, c, 0.01f);
    add_tri(&s, v);
    rebuild(&s);

    int max_depth = 0;
    for (int i = 0; i < s.nnode; i++) {
        int d = 0;
        for (int ni = i; ni >= 0; ni = s.n[ni].parent) d++;
        if (d > max_depth) max_depth = d;
    }
    check("skewed_audit", audit_contains(&s, CONTAINS_EPS) == 0);
    check("skewed_depth_bounded", max_depth <= 48);

    float o[3] = { 1e4f - 50, 1e4f - 50, 1e4f - 50 }, d[3] = { 1, 1, 1 };
    check("skewed_ray_ok", ray_mismatch(&s, o, d) == 0);
    printf("skewed_max_depth=%d nodes=%d\n", max_depth, s.nnode);
}

/* Incoherent motion: each prim is re-scattered every frame. The partition goes
 * stale fast, so the SAH proxy must cross the threshold and the criterion must
 * fire. First run refit-only to show the drift, then honour the criterion. */
static void test_animation_scatter(void) {
    const int NP = 1024;
    const int FRAMES = 32;
    const int NRAY = 128;

    Tree* T = &g_a;
    memset(T, 0, sizeof(*T));
    static float cur[1024][3][3];
    for (int i = 0; i < NP; i++) {
        float c[3] = { frand(0, 64), frand(0, 64), frand(0, 64) };
        make_tri(cur[i], c, 0.8f);
        add_tri(T, cur[i]);
    }
    rebuild(T);
    float q_build = T->q_build;

    Tree* R = &g_b;
    *R = *T;

    float v[3][3];
    int ids[1024];
    for (int i = 0; i < NP; i++) ids[i] = i;

    double peak_drift = 0.0;
    int audits_bad = 0;
    for (int f = 0; f < FRAMES; f++) {
        for (int i = 0; i < NP; i++) {
            float c[3] = { frand(0, 64), frand(0, 64), frand(0, 64) };
            make_tri(v, c, 0.8f);
            apply_tri(R, i, v);
        }
        refit(R, ids, NP);
        double d = q_now(R) / (double)q_build;
        if (d > peak_drift) peak_drift = d;
        audits_bad += audit_contains(R, CONTAINS_EPS);
    }
    check("scatter_drift_detected", peak_drift > (double)REBUILD_RATIO);
    check("scatter_drift_audit", audits_bad == 0);

    Tree* C = T;
    /* C is still the pristine build (R was a copy). */
    int rebuilds = 0, mismatches = 0;
    audits_bad = 0;
    float rays_o[128][3], rays_d[128][3];
    for (int r = 0; r < NRAY; r++) {
        int target = (int)(rnd() % (unsigned)NP);
        float tc[3] = { (cur[target][0][0] + cur[target][1][0] + cur[target][2][0]) / 3.0f,
                        (cur[target][0][1] + cur[target][1][1] + cur[target][2][1]) / 3.0f,
                        (cur[target][0][2] + cur[target][1][2] + cur[target][2][2]) / 3.0f };
        float dir[3] = { frand(-1, 1), frand(-1, 1), frand(-1, 1) };
        float len = sqrtf(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        if (len < 1e-3f) { dir[0] = 0; dir[1] = 1; dir[2] = 0; len = 1; }
        for (int i = 0; i < 3; i++) rays_d[r][i] = dir[i] / len;
        for (int i = 0; i < 3; i++) rays_o[r][i] = tc[i] - rays_d[r][i] * 70.0f;
    }

    for (int f = 0; f < FRAMES; f++) {
        for (int i = 0; i < NP; i++) {
            float c[3] = { frand(0, 64), frand(0, 64), frand(0, 64) };
            make_tri(v, c, 0.8f);
            apply_tri(C, i, v);
        }
        refit(C, ids, NP);
        if (should_rebuild(C) & REBUILD_SAH) { rebuild(C); rebuilds++; }
        audits_bad += audit_contains(C, CONTAINS_EPS);
        for (int r = 0; r < NRAY; r++) mismatches += ray_mismatch(C, rays_o[r], rays_d[r]);
    }
    check("scatter_rebuilds_triggered", rebuilds > 0);
    check("scatter_mismatch", mismatches == 0);
    check("scatter_audit", audits_bad == 0);
    printf("scatter_drift_peak=%.2f threshold=%.2f frames=%d rebuilds=%d\n",
           peak_drift, REBUILD_RATIO, FRAMES, rebuilds);
}

/* Failure detection: a stale box must not silently drop hits.
 * Two ways the update criterion can fail:
 *   1. dirty set misses the moved prim (stale prim + leaf),
 *   2. scope computation updates the leaf but not its ancestors.
 * Both are caught by audit_contains(); both produce a ray that brute force
 * hits and the stale traversal culls. After a correct refit both are clean. */
static void make_small_scene(Tree* t, float base[3][3][3]) {
    memset(t, 0, sizeof(*t));
    float c0[3] = { 2, 0.3f, 0.2f };
    float c1[3] = { 4, -0.3f, 0.2f };
    float c2[3] = { 6, 0.3f, 0.2f };
    make_tri(base[0], c0, 0.2f);
    make_tri(base[1], c1, 0.2f);
    make_tri(base[2], c2, 0.2f);
    for (int i = 0; i < 3; i++) add_tri(t, base[i]);
    rebuild(t);
}

static void move_prim_high(Tree* t, int id) {
    float v[3][3] = { { -5, 500, -5 }, { 5, 500, -5 }, { 0, 500, 5 } };
    apply_tri(t, id, v);
}

static void test_stale_detection(void) {
    float base[3][3][3];
    static Tree s;
    make_small_scene(&s, base);
    check("stale_scene_ok", audit_contains(&s, CONTAINS_EPS) == 0);

    float o[3] = { 0, 0, 0 }, d[3] = { 0, 1, 0 };
    float tb = 0;
    int idb = -1;
    check("stale_brute_hits", brute(&s, o, d, 1e30f, &tb, &idb) == 0);  /* nothing high yet */

    /* (1) Dirty-set failure: move the prim, recompute nothing. */
    static Tree stale;
    stale = s;
    move_prim_high(&stale, 0);
    int audit_bad = audit_contains(&stale, CONTAINS_EPS) > 0;
    float tt = 0;
    int idt = -1;
    int tree_hit = raycast(&stale, o, d, 1e30f, &tt, &idt);
    int brute_hit = brute(&stale, o, d, 1e30f, &tb, &idb);
    check("stale_audit_detects", audit_bad);
    check("stale_false_cull", brute_hit == 1 && tree_hit == 0);
    check("stale_brute_id", brute_hit == 1 && idb == 0);

    int ids[1] = { 0 };
    refit(&stale, ids, 1);
    check("stale_repaired", audit_contains(&stale, CONTAINS_EPS) == 0);
    check("stale_repaired_hit", ray_mismatch(&stale, o, d) == 0);

    /* (2) Scope failure: update the leaf box only, leave ancestors stale. */
    static Tree partial;
    partial = s;
    move_prim_high(&partial, 0);
    int lf = partial.leaf_of[0];
    aabb_empty(&partial.n[lf].box);
    for (int i = partial.n[lf].first; i < partial.n[lf].first + partial.n[lf].count; i++)
        aabb_grow(&partial.n[lf].box, &partial.p[i].box);
    int partial_audit = audit_contains(&partial, CONTAINS_EPS) > 0;
    int partial_tree = raycast(&partial, o, d, 1e30f, &tt, &idt);
    int partial_brute = brute(&partial, o, d, 1e30f, &tb, &idb);
    check("scope_audit_detects", partial_audit);
    check("scope_false_cull", partial_brute == 1 && partial_tree == 0);

    refit(&partial, ids, 1);
    check("scope_repaired", audit_contains(&partial, CONTAINS_EPS) == 0);
    check("scope_repaired_hit", ray_mismatch(&partial, o, d) == 0);
}

/* Topology change: appending a prim invalidates the partition. Refit cannot
 * place it, so the criterion must say TOPOLOGY and a rebuild must fix hits. */
static void test_topology_change(void) {
    float base[3][3][3];
    static Tree s;
    make_small_scene(&s, base);

    float o[3] = { 0, 0, 0 }, d[3] = { 0, 1, 0 };
    float v[3][3] = { { -5, 500, -5 }, { 5, 500, -5 }, { 0, 500, 5 } };
    add_tri(&s, v);
    check("topo_criterion", (should_rebuild(&s) & REBUILD_TOPOLOGY) != 0);

    float tb = 0;
    int idb = -1;
    int brute_hit = brute(&s, o, d, 1e30f, &tb, &idb);
    float tt = 0;
    int idt = -1;
    int tree_hit = raycast(&s, o, d, 1e30f, &tt, &idt);
    check("topo_wrong_before_rebuild", brute_hit == 1 && tree_hit == 0);

    rebuild(&s);
    check("topo_clears_criterion", (should_rebuild(&s) & REBUILD_TOPOLOGY) == 0);
    check("topo_ok_after_rebuild", ray_mismatch(&s, o, d) == 0);
    check("topo_audit_after", audit_contains(&s, CONTAINS_EPS) == 0);
}

int main(void) {
    test_legacy();
    test_refit_scope();
    test_animation_rigid();
    test_animation_partial();
    test_animation_scatter();
    test_skewed_layout();
    test_stale_detection();
    test_topology_change();

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
