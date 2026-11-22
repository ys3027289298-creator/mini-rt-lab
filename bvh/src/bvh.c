
/* bvh: median-split AABB tree over triangle bounds, nearest-hit traversal. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_PRIM 256
#define MAX_NODE 512

typedef struct { float min[3], max[3]; } AABB;
typedef struct { AABB box; float cx, cy, cz; int id; } Prim;
typedef struct { AABB box; int left, right; int first, count; } Node;

typedef struct {
    Prim p[MAX_PRIM];
    Node n[MAX_NODE];
    int nprim, nnode;
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

static int tree_add(Tree* t, float minx, float miny, float minz, float maxx, float maxy, float maxz) {
    if (t->nprim >= MAX_PRIM) return -1;
    Prim* p = &t->p[t->nprim];
    p->box.min[0] = minx; p->box.min[1] = miny; p->box.min[2] = minz;
    p->box.max[0] = maxx; p->box.max[1] = maxy; p->box.max[2] = maxz;
    p->cx = (minx + maxx) * 0.5f; p->cy = (miny + maxy) * 0.5f; p->cz = (minz + maxz) * 0.5f;
    p->id = t->nprim;
    return t->nprim++;
}

static int build(Tree* t, int first, int count) {
    if (t->nnode >= MAX_NODE) return -1;
    int ni = t->nnode++;
    Node* n = &t->n[ni];
    aabb_empty(&n->box);
    for (int i = first; i < first + count; i++) aabb_grow(&n->box, &t->p[i].box);

    if (count <= 2) {
        n->left = n->right = -1;
        n->first = first; n->count = count;
        return ni;
    }
    /* split on widest axis at centroid median */
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
    n->first = -1; n->count = 0;
    n->left = build(t, first, mid - first);
    n->right = build(t, mid, first + count - mid);
    return ni;
}

/* ---- traversal ----
   Three queries over the same tree:
   - occlude_full : enters every node and tests every primitive, with no
                    early out — the "each shadow ray walks the whole tree"
                    baseline and the ground-truth occlusion predicate.
   - occlude_early: returns on the first hit, skips whole subtrees whose
                    entry t exceeds the bound, tests the near child first so
                    the bound shrinks as fast as possible.
   - nearest_hit  : ordered descent with the current best t as the moving
                    upper bound; checks occlusion and t against brute force.

   `visits` counts nodes entered, identically for full and early, so the
   reduction ratio measures only what early-out saves. `probes`, passed to
   occlude_early only, counts the extra child-box tests it issues before
   pushing a child; the full traversal never issues those because it does not
   prune. Both are reported so the ratio is not flattering. */

typedef struct { int node; float entry; } StackEntry;

static int occlude_full(const Tree* t, int root, const float o[3], const float d[3],
                        float tmax, unsigned long* visits) {
    int stack[MAX_NODE], sp = 0;
    int any = 0;
    stack[sp++] = root;
    while (sp > 0) {
        const Node* n = &t->n[stack[--sp]];
        (*visits)++;
        if (n->left < 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                float tm;
                if (aabb_hit(&t->p[i].box, o, d, tmax, &tm)) any = 1;
            }
        } else {
            stack[sp++] = n->left;
            stack[sp++] = n->right;
        }
    }
    return any;
}

static int occlude_early(const Tree* t, int root, const float o[3], const float d[3],
                         float tmax, unsigned long* visits, unsigned long* probes) {
    StackEntry st[MAX_NODE];
    int sp = 0;
    float t0;
    (*visits)++;
    if (!aabb_hit(&t->n[root].box, o, d, tmax, &t0)) return 0;
    st[sp].node = root; st[sp].entry = t0; sp++;
    while (sp > 0) {
        StackEntry e = st[--sp];
        (*visits)++;
        if (e.entry > tmax) continue;
        const Node* n = &t->n[e.node];
        if (n->left < 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                float tm;
                if (aabb_hit(&t->p[i].box, o, d, tmax, &tm)) return 1;
            }
        } else {
            float tl = 0.0f, tr = 0.0f;
            (*probes) += 2;
            int hl = aabb_hit(&t->n[n->left].box, o, d, tmax, &tl);
            int hr = aabb_hit(&t->n[n->right].box, o, d, tmax, &tr);
            /* LIFO: push far child first so the near child pops first */
            if (hl && hr) {
                if (tl <= tr) {
                    st[sp].node = n->right; st[sp].entry = tr; sp++;
                    st[sp].node = n->left;  st[sp].entry = tl; sp++;
                } else {
                    st[sp].node = n->left;  st[sp].entry = tl; sp++;
                    st[sp].node = n->right; st[sp].entry = tr; sp++;
                }
            } else if (hl) {
                st[sp].node = n->left;  st[sp].entry = tl; sp++;
            } else if (hr) {
                st[sp].node = n->right; st[sp].entry = tr; sp++;
            }
        }
    }
    return 0;
}

static int nearest_hit(const Tree* t, int root, const float o[3], const float d[3],
                       float tmax, float* out_t, int* out_id, unsigned long* visits) {
    StackEntry st[MAX_NODE];
    int sp = 0, found = 0, best_id = -1;
    float best = tmax, t0;
    (*visits)++;
    if (!aabb_hit(&t->n[root].box, o, d, tmax, &t0)) return 0;
    st[sp].node = root; st[sp].entry = t0; sp++;
    while (sp > 0) {
        StackEntry e = st[--sp];
        (*visits)++;
        if (e.entry > best) continue;
        const Node* n = &t->n[e.node];
        if (n->left < 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                float tm;
                if (!aabb_hit(&t->p[i].box, o, d, best, &tm)) continue;
                if (!found || tm < best || (tm == best && t->p[i].id < best_id)) {
                    best = tm; best_id = t->p[i].id; found = 1;
                }
            }
        } else {
            float tl = 0.0f, tr = 0.0f;
            int hl = aabb_hit(&t->n[n->left].box, o, d, best, &tl);
            int hr = aabb_hit(&t->n[n->right].box, o, d, best, &tr);
            if (hl && hr) {
                if (tl <= tr) {
                    st[sp].node = n->right; st[sp].entry = tr; sp++;
                    st[sp].node = n->left;  st[sp].entry = tl; sp++;
                } else {
                    st[sp].node = n->left;  st[sp].entry = tl; sp++;
                    st[sp].node = n->right; st[sp].entry = tr; sp++;
                }
            } else if (hl) {
                st[sp].node = n->left;  st[sp].entry = tl; sp++;
            } else if (hr) {
                st[sp].node = n->right; st[sp].entry = tr; sp++;
            }
        }
    }
    if (found) { *out_t = best; *out_id = best_id; }
    return found;
}

/* O(n) reference used to validate nearest_hit, same predicate and tie-break. */
static int brute_nearest(const Tree* t, const float o[3], const float d[3],
                         float tmax, float* out_t, int* out_id) {
    int found = 0, best_id = -1;
    float best = tmax;
    for (int i = 0; i < t->nprim; i++) {
        float tm;
        if (!aabb_hit(&t->p[i].box, o, d, tmax, &tm)) continue;
        if (!found || tm < best || (tm == best && t->p[i].id < best_id)) {
            best = tm; best_id = t->p[i].id; found = 1;
        }
    }
    if (found) { *out_t = best; *out_id = best_id; }
    return found;
}

/* ---- small vector helpers for the test camera ---- */
static float v_dot(const float a[3], const float b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static void v_cross(const float a[3], const float b[3], float out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}
static void v_norm(float v[3]) {
    float l = sqrtf(v_dot(v, v));
    if (l > 0.0f) { v[0] /= l; v[1] /= l; v[2] /= l; }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Tree t;
    memset(&t, 0, sizeof(t));
    tree_add(&t, 1, 0, 0, 2, 1, 1);
    tree_add(&t, 5, 0, 0, 6, 1, 1);
    tree_add(&t, 3, 0, 0, 4, 1, 1);
    check("prim_count", t.nprim == 3);

    int root = build(&t, 0, t.nprim);
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

    Tree big;
    memset(&big, 0, sizeof(big));
    int rc = 0;
    for (int i = 0; i < MAX_PRIM; i++) rc = tree_add(&big, (float)i, 0, 0, (float)i + 0.5f, 1, 1);
    check("capacity", rc == MAX_PRIM - 1 && big.nprim == MAX_PRIM);
    check("build_full", build(&big, 0, big.nprim) == 0);

    /* ---- early-out traversal vs full traversal, pixel by pixel ---- */
    Tree grid;
    memset(&grid, 0, sizeof(grid));
    /* 8x8x2 staggered unit-ish boxes in [-8,8]^3, one box deliberately
       straddling the camera ray path so early exit triggers often. */
    for (int iz = 0; iz < 2; iz++)
        for (int iy = 0; iy < 8; iy++)
            for (int ix = 0; ix < 8; ix++) {
                float x = -7.0f + 2.0f * ix + 0.25f * ((iy + iz) & 1);
                float y = -7.0f + 2.0f * iy;
                float z = -4.0f + 6.0f * iz;
                tree_add(&grid, x, y, z, x + 0.9f, y + 0.9f, z + 0.9f);
            }
    check("grid_prims", grid.nprim == 128);
    int grid_root = build(&grid, 0, grid.nprim);
    check("grid_build", grid_root == 0);

    const int W = 128, H = 128;
    float fwd[3] = { 0, 0, 1 }, upv[3] = { 0, 1, 0 };
    float right[3];
    v_cross(fwd, upv, right);
    v_norm(right);
    const float fov = 0.9f, tmax = 100.0f;

    unsigned long vis_full = 0, vis_early = 0, probes_early = 0;
    int pixels = 0, occ_full = 0, occ_early = 0, early_fp = 0, early_fn = 0;
    int nearest_mismatch = 0, nearest_found = 0, nearest_vs_occl = 0;
    int early_visit_overflow = 0;

    for (int py = 0; py < H; py++) {
        for (int px = 0; px < W; px++) {
            float u = (2.0f * (px + 0.5f) / W - 1.0f);
            float v = (1.0f - 2.0f * (py + 0.5f) / H);
            float d[3] = { fwd[0] + right[0] * u * fov,
                           fwd[1] + right[1] * u * fov + upv[1] * v * fov,
                           fwd[2] + right[2] * u * fov + upv[2] * v * fov };
            v_norm(d);
            const float o[3] = { 0.13f, 0.07f, -9.0f };  /* slightly off-axis */

            unsigned long vf = 0, ve = 0, pe = 0;
            int a = occlude_full(&grid, grid_root, o, d, tmax, &vf);
            int b = occlude_early(&grid, grid_root, o, d, tmax, &ve, &pe);
            vis_full += vf;
            vis_early += ve;
            probes_early += pe;
            pixels++;
            if (a) occ_full++;
            if (b) occ_early++;
            if (!a && b) early_fp++;       /* early-out claims hit, full says miss */
            if (a && !b) early_fn++;       /* forbidden: false negative */
            if (ve > vf) early_visit_overflow++;  /* entering a node must be shared work */

            /* nearest hit must agree with the occlusion predicate, and with brute force */
            float nt = 0, bt = 0;
            int nid = -1, bid = -1;
            unsigned long vn = 0;
            int nh = nearest_hit(&grid, grid_root, o, d, tmax, &nt, &nid, &vn);
            int bh = brute_nearest(&grid, o, d, tmax, &bt, &bid);
            if (nh != bh || (nh && (nt != bt || nid != bid))) nearest_mismatch++;
            if (nh) nearest_found++;
            if (nh != a) nearest_vs_occl++;
        }
    }

    check("pixel_count", pixels == W * H);
    check("early_false_negatives_zero", early_fn == 0);
    check("early_false_positives_zero", early_fp == 0);
    check("early_matches_full", occ_full == occ_early);
    check("nearest_matches_brute", nearest_mismatch == 0);
    check("nearest_matches_occlusion", nearest_vs_occl == 0);
    check("scene_has_hits_and_misses", occ_full > 0 && occ_full < pixels);

    /* early-out must visit strictly fewer nodes overall on this scene */
    check("visits_reduced", vis_early < vis_full);
    check("early_visits_bounded_by_full", early_visit_overflow == 0);

    printf("occluded_px=%d/%d\n", occ_full, pixels);
    printf("visits_full=%lu\n", vis_full);
    printf("visits_early=%lu\n", vis_early);
    printf("visits_early_plus_probes=%lu\n", vis_early + probes_early);
    printf("visit_reduction=%.4f\n", vis_full ? 1.0 - (double)vis_early / (double)vis_full : 0.0);
    printf("work_reduction_incl_probes=%.4f\n",
           vis_full ? 1.0 - (double)(vis_early + probes_early) / (double)vis_full : 0.0);

    /* targeted cases: on-axis ray into the near box group must early-exit
       after a tiny fraction of the tree. */
    {
        const float o2[3] = { 1.5f, -0.5f, -9.0f }, d2[3] = { 0, 0, 1 };
        unsigned long vf = 0, ve = 0, pe = 0;
        int a = occlude_full(&grid, grid_root, o2, d2, tmax, &vf);
        int b = occlude_early(&grid, grid_root, o2, d2, tmax, &ve, &pe);
        check("axis_hit_agrees", a == 1 && b == 1);
        check("axis_hit_prunes", ve * 4 < vf);
        printf("axis_visits_full=%lu axis_visits_early=%lu\n", vf, ve);
    }
    {
        /* ray pointing away from the scene: full traversal still walks the
           whole tree (root box test fails first, so both are cheap here);
           verify both report miss. */
        const float o3[3] = { 0, 0, -9.0f }, d3[3] = { 0, 0, -1 };
        unsigned long vf = 0, ve = 0, pe = 0;
        check("away_miss_agrees",
              occlude_full(&grid, grid_root, o3, d3, tmax, &vf) == 0 &&
              occlude_early(&grid, grid_root, o3, d3, tmax, &ve, &pe) == 0);
    }
    /* t upper bound: the only boxes on the axis sit at z in [-4.0,-3.1] and
       [2.0,2.9], i.e. entry t = 5 and 11 from z=-9. A bound of 4 is before
       the first entry, so both queries must report a miss. */
    {
        const float o4[3] = { 1.5f, -0.5f, -9.0f }, d4[3] = { 0, 0, 1 };
        unsigned long vf = 0, ve = 0, pe = 0;
        int a = occlude_full(&grid, grid_root, o4, d4, 4.0f, &vf);
        int b = occlude_early(&grid, grid_root, o4, d4, 4.0f, &ve, &pe);
        check("tbound_excludes_both", a == 0 && b == 0);
        check("tbound_prunes_early", ve < vf);
        printf("tbound_visits_full=%lu tbound_visits_early=%lu\n", vf, ve);
    }
    /* nearest hit under the same bound must also miss; with a bound past the
       first box it must return that box's entry distance and agree with the
       early-out predicate. */
    {
        const float o5[3] = { 1.5f, -0.5f, -9.0f }, d5[3] = { 0, 0, 1 };
        float nt = 0;
        int nid = -1;
        unsigned long vn = 0;
        check("nearest_tbound_miss", nearest_hit(&grid, grid_root, o5, d5, 4.0f, &nt, &nid, &vn) == 0);
        check("nearest_tbound_hit",
              nearest_hit(&grid, grid_root, o5, d5, 6.0f, &nt, &nid, &vn) == 1 &&
              fabsf(nt - 5.0f) < 1e-4f);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
