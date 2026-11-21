
/* bvh: median-split AABB tree over triangle bounds, nearest-hit traversal.
   Nearest-hit accuracy is validated per-pixel against brute force:
   both sides sort/select candidates by (t, id) so coplanar and
   shared-edge ties resolve identically; degenerate triangles and
   t <= 0 hits are listed separately, never selected. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_PRIM 32
#define MAX_NODE 64

typedef struct { float min[3], max[3]; } AABB;
typedef struct {
    AABB box;
    float cx, cy, cz;
    int id;
    float v[3][3];   /* triangle vertices (all zero for box-only prims) */
    int degenerate;  /* zero-area triangle */
} Prim;
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
    memset(p->v, 0, sizeof(p->v));
    p->degenerate = 0;
    p->id = t->nprim;
    return t->nprim++;
}

static int tree_add_tri(Tree* t, const float v0[3], const float v1[3], const float v2[3]) {
    if (t->nprim >= MAX_PRIM) return -1;
    Prim* p = &t->p[t->nprim];
    memcpy(p->v[0], v0, sizeof(p->v[0]));
    memcpy(p->v[1], v1, sizeof(p->v[1]));
    memcpy(p->v[2], v2, sizeof(p->v[2]));
    aabb_empty(&p->box);
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < 3; i++) {
            if (p->v[k][i] < p->box.min[i]) p->box.min[i] = p->v[k][i];
            if (p->v[k][i] > p->box.max[i]) p->box.max[i] = p->v[k][i];
        }
    p->cx = (p->box.min[0] + p->box.max[0]) * 0.5f;
    p->cy = (p->box.min[1] + p->box.max[1]) * 0.5f;
    p->cz = (p->box.min[2] + p->box.max[2]) * 0.5f;
    float e1[3], e2[3];
    for (int i = 0; i < 3; i++) { e1[i] = v1[i] - v0[i]; e2[i] = v2[i] - v0[i]; }
    float nx = e1[1] * e2[2] - e1[2] * e2[1];
    float ny = e1[2] * e2[0] - e1[0] * e2[2];
    float nz = e1[0] * e2[1] - e1[1] * e2[0];
    p->degenerate = (nx * nx + ny * ny + nz * nz) < 1e-12f;
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

/* Möller–Trumbore. Returns 1 when the ray crosses the triangle at any t
   (sign included); *t_out receives the raw t so the caller can route
   t <= 0 hits to a separate list instead of the nearest-hit selection. */
static int tri_intersect(const Prim* p, const float o[3], const float d[3], float* t_out) {
    const float* v0 = p->v[0];
    const float* v1 = p->v[1];
    const float* v2 = p->v[2];
    float e1[3], e2[3], tv[3];
    for (int i = 0; i < 3; i++) { e1[i] = v1[i] - v0[i]; e2[i] = v2[i] - v0[i]; tv[i] = o[i] - v0[i]; }
    float px = d[1] * e2[2] - d[2] * e2[1];
    float py = d[2] * e2[0] - d[0] * e2[2];
    float pz = d[0] * e2[1] - d[1] * e2[0];
    float det = e1[0] * px + e1[1] * py + e1[2] * pz;
    if (fabsf(det) < 1e-8f) return 0;   /* ray parallel to plane, or degenerate triangle */
    float inv = 1.0f / det;
    float u = (tv[0] * px + tv[1] * py + tv[2] * pz) * inv;
    if (u < 0.0f || u > 1.0f) return 0;
    float qx = tv[1] * e1[2] - tv[2] * e1[1];
    float qy = tv[2] * e1[0] - tv[0] * e1[2];
    float qz = tv[0] * e1[1] - tv[1] * e1[0];
    float v = (d[0] * qx + d[1] * qy + d[2] * qz) * inv;
    if (v < 0.0f || u + v > 1.0f) return 0;
    *t_out = (e2[0] * qx + e2[1] * qy + e2[2] * qz) * inv;
    return 1;
}

typedef struct { int id; float t; } Hit;

/* deterministic candidate order: smaller t first, ties (coplanar overlap,
   shared edges) broken by smaller prim id — identical on both sides,
   which is what makes the tree-vs-brute-force comparison exactly zero */
static int hit_cmp(const void* a, const void* b) {
    const Hit* x = (const Hit*)a;
    const Hit* y = (const Hit*)b;
    if (x->t < y->t) return -1;
    if (x->t > y->t) return 1;
    return x->id - y->id;
}
static int hit_better(const Hit* h, const Hit* best) {
    return best->id < 0 || h->t < best->t || (h->t == best->t && h->id < best->id);
}

/* brute-force reference: test every triangle, sort all t > 0 candidates
   by (t, id), take the first. t <= 0 hits are counted into n_neg/neg_mask
   (per-prim bit) and never become candidates. */
static Hit brute_nearest(const Tree* tr, const float o[3], const float d[3],
                         Hit* cand, int* ncand, int* n_neg, unsigned* neg_mask) {
    Hit local[MAX_PRIM];
    if (!cand) cand = local;
    int nc = 0;
    for (int i = 0; i < tr->nprim; i++) {
        float t;
        if (!tri_intersect(&tr->p[i], o, d, &t)) continue;
        if (t <= 0.0f) {
            if (n_neg) (*n_neg)++;
            if (neg_mask) *neg_mask |= 1u << tr->p[i].id;
            continue;
        }
        cand[nc].id = tr->p[i].id;
        cand[nc].t = t;
        nc++;
    }
    qsort(cand, (size_t)nc, sizeof(Hit), hit_cmp);
    if (ncand) *ncand = nc;
    Hit best = { -1, 0.0f };
    if (nc > 0) best = cand[0];
    return best;
}

/* BVH nearest hit. Children are pushed far-first so the nearer subtree is
   explored first; a child is skipped only when its box entry is beyond the
   current best t, with a small slack so float rounding in the slab test can
   never prune the true winner (the overlapping-boxes failure mode).
   Leaf selection uses the same (t, id) order as brute_nearest. */
static Hit tree_nearest(const Tree* tr, const float o[3], const float d[3]) {
    Hit best = { -1, 0.0f };
    int stack[MAX_NODE], sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        const Node* n = &tr->n[stack[--sp]];
        if (n->left < 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                float t;
                if (!tri_intersect(&tr->p[i], o, d, &t)) continue;
                if (t <= 0.0f) continue;   /* listed separately via brute force */
                Hit h = { tr->p[i].id, t };
                if (hit_better(&h, &best)) best = h;
            }
        } else {
            float tmax = best.id < 0 ? 1e30f : best.t * 1.0001f + 1e-5f;
            float tl = 0.0f, trr = 0.0f;
            int hl = aabb_hit(&tr->n[n->left].box, o, d, tmax, &tl);
            int hr = aabb_hit(&tr->n[n->right].box, o, d, tmax, &trr);
            if (hl && hr) {
                if (tl <= trr) { stack[sp++] = n->right; stack[sp++] = n->left; }
                else           { stack[sp++] = n->left;  stack[sp++] = n->right; }
            } else if (hl) {
                stack[sp++] = n->left;
            } else if (hr) {
                stack[sp++] = n->right;
            }
        }
    }
    return best;
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

    /* ---- nearest-hit validation: BVH vs brute force ----
       Scene ids: 0 back wall z=5; 1 front tri z=3 (overlaps 0);
       2,3 quad halves z=4 (shared edge); 4,5 coplanar overlapping z=4.5;
       6 degenerate (collinear); 7 behind the camera (t < 0). */
    Tree s;
    memset(&s, 0, sizeof(s));
    tree_add_tri(&s, (float[]){ -3, -3, 5 }, (float[]){ 3, -3, 5 }, (float[]){ 0, 3, 5 });
    tree_add_tri(&s, (float[]){ -1, -1, 3 }, (float[]){ 1, -1, 3 }, (float[]){ 0, 1, 3 });
    tree_add_tri(&s, (float[]){ 1, 0, 4 }, (float[]){ 2, 0, 4 }, (float[]){ 2, 1, 4 });
    tree_add_tri(&s, (float[]){ 1, 0, 4 }, (float[]){ 2, 1, 4 }, (float[]){ 1, 1, 4 });
    tree_add_tri(&s, (float[]){ -2, 0, 4.5f }, (float[]){ -1, 0, 4.5f }, (float[]){ -1, 1, 4.5f });
    tree_add_tri(&s, (float[]){ -1.75f, -0.25f, 4.5f }, (float[]){ -0.75f, -0.25f, 4.5f }, (float[]){ -1.25f, 0.75f, 4.5f });
    tree_add_tri(&s, (float[]){ 0, 0, 6 }, (float[]){ 1, 1, 6 }, (float[]){ 2, 2, 6 });
    tree_add_tri(&s, (float[]){ -0.5f, -0.5f, -2 }, (float[]){ 0.5f, -0.5f, -2 }, (float[]){ 0, 0.5f, -2 });
    check("scene_prims", s.nprim == 8);
    check("scene_build", build(&s, 0, s.nprim) == 0);

    float zo[3] = { 0, 0, 0 };

    /* degenerate triangles: listed separately, never hit */
    int deg_ok = 1;
    for (int i = 0; i < s.nprim; i++)
        if (s.p[i].degenerate != (s.p[i].id == 6)) deg_ok = 0;
    check("degenerate_flag", deg_ok);
    {
        const Prim* deg = 0;
        for (int i = 0; i < s.nprim; i++) if (s.p[i].id == 6) deg = &s.p[i];
        float dd[3] = { 1, 1, 6 }, ttmp = 0;
        check("degenerate_no_hit", deg != 0 && !tri_intersect(deg, zo, dd, &ttmp));
        printf("degenerate_prims=");
        int first = 1;
        for (int i = 0; i < s.nprim; i++)
            if (s.p[i].degenerate) { printf("%s%d", first ? "" : ",", s.p[i].id); first = 0; }
        printf("\n");
    }

    /* overlapping objects: the front triangle must win in both */
    {
        float zd[3] = { 0, 0, 1 };
        Hit ob = brute_nearest(&s, zo, zd, 0, 0, 0, 0);
        Hit ot = tree_nearest(&s, zo, zd);
        check("overlap_front_wins", ob.id == 1 && ob.t == 3.0f);
        check("overlap_tree_eq_brute", ot.id == ob.id && ot.t == ob.t);
    }

    /* shared edge: ray through the quad's diagonal midpoint hits both
       halves with the same t; the tie resolves to the smaller id */
    {
        float de[3] = { 1.5f, 0.5f, 4.0f };
        Hit cand[MAX_PRIM];
        int nc = 0;
        Hit eb = brute_nearest(&s, zo, de, cand, &nc, 0, 0);
        Hit et = tree_nearest(&s, zo, de);
        check("edge_two_candidates", nc == 2 && cand[0].t == cand[1].t);
        check("edge_tie_by_id", eb.id == 2 && cand[1].id == 3);
        check("edge_tree_eq_brute", et.id == eb.id && et.t == eb.t);
    }

    /* coplanar overlap: both triangles at z=4.5, same plane, same t */
    {
        float dc[3] = { -1.25f, 0.25f, 4.5f };
        Hit cand[MAX_PRIM];
        int nc = 0;
        Hit cb = brute_nearest(&s, zo, dc, cand, &nc, 0, 0);
        Hit ct = tree_nearest(&s, zo, dc);
        check("coplanar_two_candidates", nc == 2 &&
              ((cand[0].id == 4 && cand[1].id == 5) || (cand[0].id == 5 && cand[1].id == 4)));
        check("coplanar_tree_eq_brute", ct.id == cb.id && ct.t == cb.t);
    }

    /* per-pixel sweep: tree and brute force must agree exactly */
    int mism = 0, neg_total = 0;
    unsigned neg_mask = 0;
    for (int py = 0; py < 64; py++) {
        for (int px = 0; px < 64; px++) {
            float dpx[3] = { -2.0f + 4.0f * (px + 0.5f) / 64.0f,
                             -2.0f + 4.0f * (py + 0.5f) / 64.0f, 1.0f };
            int nn = 0;
            unsigned nm = 0;
            Hit a = brute_nearest(&s, zo, dpx, 0, 0, &nn, &nm);
            Hit b = tree_nearest(&s, zo, dpx);
            neg_total += nn;
            neg_mask |= nm;
            if (a.id != b.id || a.t != b.t) mism++;
        }
    }
    printf("pixels=%d mismatches=%d\n", 64 * 64, mism);
    check("pixels_equal", mism == 0);

    /* t <= 0 hits: listed separately, excluded from nearest */
    {
        printf("negative_t_prims=");
        int first = 1;
        for (int i = 0; i < MAX_PRIM; i++)
            if (neg_mask & (1u << i)) { printf("%s%d", first ? "" : ",", i); first = 0; }
        printf("\n");
    }
    check("negative_t_only_7", neg_mask == (1u << 7));
    check("negative_t_count", neg_total > 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
