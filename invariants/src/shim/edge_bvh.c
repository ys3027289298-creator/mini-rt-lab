/*
 * edge_bvh.c -- observer shim for bvh/src/bvh.c.
 *
 * The module is compiled *unmodified*; only its entry point is renamed so its
 * self-test main() does not collide with ours.  Everything here is a pure
 * forwarder except the two constructors (tree_add/build), which are the
 * module's own.  All module symbols are static inside bvh.c, so the adapters
 * live in this one TU.
 */
#define main bvh_module_main
#define failures bvh_module_failures
#include "../../bvh/src/bvh.c"
#undef main
#undef failures

#include "inv_impl.h"

/* The module works in float; the checker's contract is double.  These adapters
 * are the only place the two meet, and they pass through unchanged -- no
 * nudging of inputs or outputs. */
int edge_box_hit(const double lo[3], const double hi[3],
                 const double o[3], const double d[3],
                 double tmax, double* t) {
    AABB b;
    float of[3], df[3];
    float tf = 0.0f;
    int hit;
    for (int i = 0; i < 3; i++) {
        b.min[i] = (float)lo[i];
        b.max[i] = (float)hi[i];
        of[i]    = (float)o[i];
        df[i]    = (float)d[i];
    }
    hit = aabb_hit(&b, of, df, (float)tmax, &tf);
    if (t) *t = (double)tf;
    return hit;
}

int edge_build(InvTree* out) {
    static Tree t;
    InvBox scene[INV_MAX_PRIM];
    int n = inv_scene_boxes(scene);
    memset(&t, 0, sizeof t);
    for (int i = 0; i < n; i++) {
        tree_add(&t,
                 (float)scene[i].lo[0], (float)scene[i].lo[1], (float)scene[i].lo[2],
                 (float)scene[i].hi[0], (float)scene[i].hi[1], (float)scene[i].hi[2]);
    }
    int root = build(&t, 0, t.nprim);
    if (root != 0 || t.nprim != n) return -1;
    out->nnodes = t.nnode;
    out->nprims = t.nprim;
    for (int i = 0; i < t.nnode && i < INV_MAX_NODE; i++) {
        out->node[i].box.lo[0] = t.n[i].box.min[0];
        out->node[i].box.lo[1] = t.n[i].box.min[1];
        out->node[i].box.lo[2] = t.n[i].box.min[2];
        out->node[i].box.hi[0] = t.n[i].box.max[0];
        out->node[i].box.hi[1] = t.n[i].box.max[1];
        out->node[i].box.hi[2] = t.n[i].box.max[2];
        out->node[i].left  = t.n[i].left;
        out->node[i].right = t.n[i].right;
        out->node[i].first = t.n[i].first;
        out->node[i].count = t.n[i].count;
    }
    for (int i = 0; i < t.nprim && i < INV_MAX_PRIM; i++) {
        out->prim[i].lo[0] = t.p[i].box.min[0];
        out->prim[i].lo[1] = t.p[i].box.min[1];
        out->prim[i].lo[2] = t.p[i].box.min[2];
        out->prim[i].hi[0] = t.p[i].box.max[0];
        out->prim[i].hi[1] = t.p[i].box.max[1];
        out->prim[i].hi[2] = t.p[i].box.max[2];
    }
    return 0;
}

int edge_bvh_failures(void) { return bvh_module_failures; }
