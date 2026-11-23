
/* bvh: median-split AABB tree over triangle bounds, nearest-hit traversal. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_PRIM 32
#define MAX_NODE 64

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
    int axis = (w[0] < w[1]) ? (w[0] < w[2] ? 0 : 2) : (w[1] < w[2] ? 1 : 2);
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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
