/*
 * bvh_shim.c -- compile the real bvh/src/bvh.c into the renderer.
 * BVH_SRC is supplied by the harness as a quoted absolute path (pristine or
 * mutant sandbox).  The module's own self-test main is renamed, not removed,
 * so a mutant cannot silently change what we link against.
 */
#include "iface.h"

#ifndef BVH_SRC
#error "BVH_SRC must be defined to the path of bvh.c"
#endif

#define main bvh_module_selftest_main
#include BVH_SRC
#undef main

static Tree g_tree;
static int g_root = -1;
static int g_ok = 0;

int rt_bvh_build(const float* boxes, int n) {
    memset(&g_tree, 0, sizeof(g_tree));
    g_root = -1;
    g_ok = 0;
    if (n < 0 || n > MAX_PRIM) return 0;
    for (int i = 0; i < n; i++) {
        const float* b = boxes + 6 * i;
        if (tree_add(&g_tree, b[0], b[1], b[2], b[3], b[4], b[5]) < 0) return 0;
    }
    if (g_tree.nprim == 0) { g_ok = 1; return 1; }
    g_root = build(&g_tree, 0, g_tree.nprim);
    if (g_root < 0) return 0;
    g_ok = 1;
    return 1;
}

int rt_bvh_prim_count(void) { return g_tree.nprim; }

/*
 * Nearest-hit traversal.  When two primitives are hit at *exactly* the same t
 * (coincident faces of overlapping boxes), the lower primitive id wins.  The
 * build permutes the primitive array, so without this rule the reported box
 * would depend on tree shape, and a shape-only change would look like a visual
 * difference.  An explicit stack is used instead of recursion so
 * the *visit order* is fixed (left subtree first, siblings pushed right-first)
 * and does not depend on the compiler's call layout.  The order is a pure
 * function of the tree, so two builds of this file agree bit for bit.
 */
static void walk(const float o[3], const float d[3], float* best, int* bp) {
    int stack[MAX_NODE];
    int sp = 0;
    stack[sp++] = g_root;
    while (sp > 0) {
        int ni = stack[--sp];
        if (ni < 0 || ni >= g_tree.nnode) continue;
        const Node* n = &g_tree.n[ni];
        float t;
        if (!aabb_hit(&n->box, o, d, *best, &t)) continue;
        if (n->left < 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                int id = g_tree.p[i].id;
                if (aabb_hit(&g_tree.p[i].box, o, d, *best, &t) &&
                    (t < *best || (t == *best && id < *bp))) {
                    *best = t;
                    *bp = id;
                }
            }
            continue;
        }
        /* push right before left so the left child is visited first */
        if (n->right >= 0 && sp < MAX_NODE) stack[sp++] = n->right;
        if (n->left >= 0 && sp < MAX_NODE) stack[sp++] = n->left;
    }
}

int rt_bvh_hit(const float o[3], const float d[3], float tmax, int* box, float* t) {
    if (!g_ok || g_tree.nprim == 0) return 0;
    float best = tmax;
    int bp = 1 << 30;
    walk(o, d, &best, &bp);
    if (bp == (1 << 30)) return 0;
    *box = bp;
    *t = best;
    return 1;
}

int rt_bvh_hit_brute(const float o[3], const float d[3], float tmax, int* box, float* t) {
    if (!g_ok || g_tree.nprim == 0) return 0;
    float best = tmax;
    int bp = 1 << 30;
    for (int i = 0; i < g_tree.nprim; i++) {
        int id = g_tree.p[i].id;
        float tt;
        if (aabb_hit(&g_tree.p[i].box, o, d, best, &tt) &&
            (tt < best || (tt == best && id < bp))) {
            best = tt; bp = id;
        }
    }
    if (bp == (1 << 30)) return 0;
    *box = bp; *t = best;
    return 1;
}
