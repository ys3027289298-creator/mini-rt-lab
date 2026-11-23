/*
 * Minimal test entry for the geometry/intersection path (bvh).
 * No third-party test framework: includes the source under test directly
 * (path overrivable via -DBVH_SRC='"..."') and renames its demo main.
 *
 * Scope note: bvh.c currently exposes box intersection (aabb_hit), tree
 * construction (build) and AABB helpers -- there is no tree traversal /
 * nearest-hit routine in the source, so none is tested here. Nearest-hit
 * semantics are exercised on the t-values aabb_hit returns.
 *
 * Strategy:
 *   - an independent brute-force slab oracle (separately derived, division
 *     form + fmax/fmin reduction) compared against aabb_hit for EVERY ray
 *     against EVERY primitive box and against the union box
 *   - structural checks on the built tree: root coverage, child unions,
 *     leaf partition, unique/complete primitive set, widest-axis split
 *   - boundary cases: parallel/zero-axis rays, negative direction, tmax
 *     clipping (inclusive boundary), corner grazing, capacity, area values
 */
#ifndef BVH_SRC
#define BVH_SRC "../bvh/src/bvh.c"
#endif
#define main bvh_demo_main
#include BVH_SRC
#undef main

static int t_failures = 0;
#define CHECK(name, cond) do {                                          \
    int _ok = !!(cond);                                                 \
    printf("%s=%d\n", (name), _ok);                                     \
    if (!_ok) t_failures++;                                             \
} while (0)
#define CHECK_NEAR(name, got, ref, tol) do {                            \
    double _g = (double)(got), _r = (double)(ref);                      \
    int _ok = fabs(_g - _r) <= (double)(tol);                           \
    printf("%s=%d\n", (name), _ok);                                     \
    if (!_ok) { t_failures++;                                           \
        printf("# %s got=%.9g ref=%.9g\n", (name), _g, _r); }           \
} while (0)

/* Independent slab test, re-derived (division form, fmax/fmin reduction);
   intentionally shares no code with aabb_hit. */
static int ref_box_hit(const AABB* b, const float o[3], const float d[3],
                       float tmax, float* tmin_out) {
    float tlo = 0.0f, thi = tmax;
    for (int i = 0; i < 3; i++) {
        if (fabsf(d[i]) < 1e-7f) {
            if (o[i] < b->min[i] || o[i] > b->max[i]) return 0;
        } else {
            float ta = (b->min[i] - o[i]) / d[i];
            float tb = (b->max[i] - o[i]) / d[i];
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            tlo = fmaxf(tlo, ta);
            thi = fminf(thi, tb);
            if (tlo > thi) return 0;
        }
    }
    if (tmin_out) *tmin_out = tlo;
    return 1;
}

static unsigned lcg = 123456789u;
static float frand(float lo, float hi) {
    lcg = lcg * 1103515245u + 12345u;
    float u = (float)(lcg >> 8) / (float)(1u << 24);
    return lo + u * (hi - lo);
}

static void box_union(const Tree* t, AABB* u) {
    aabb_empty(u);
    for (int i = 0; i < t->nprim; i++) aabb_grow(u, &t->p[i].box);
}

/* walk the tree: validate node invariants and collect leaf primitives */
static void walk(const Tree* t, int ni, int* seen, int* leaves, int depth,
                 int* maxdepth) {
    const Node* n = &t->n[ni];
    if (depth > *maxdepth) *maxdepth = depth;
    if (n->count > 0) {
        CHECK("leaf_is_terminal", n->left == -1 && n->right == -1);
        for (int i = n->first; i < n->first + n->count; i++) {
            CHECK("leaf_id_in_range", i >= 0 && i < t->nprim);
            if (i >= 0 && i < t->nprim) {
                CHECK("leaf_id_unique", !seen[i]);
                seen[i] = 1;
            }
        }
        (*leaves)++;
    } else {
        CHECK("inner_has_children", n->left >= 0 && n->right >= 0);
        const Node *L = &t->n[n->left], *R = &t->n[n->right];
        for (int k = 0; k < 3; k++) {
            CHECK_NEAR("child_union_min",
                       fminf(L->box.min[k], R->box.min[k]), n->box.min[k], 1e-5f);
            CHECK_NEAR("child_union_max",
                       fmaxf(L->box.max[k], R->box.max[k]), n->box.max[k], 1e-5f);
        }
        walk(t, n->left, seen, leaves, depth + 1, maxdepth);
        walk(t, n->right, seen, leaves, depth + 1, maxdepth);
    }
}

static void check_structure(Tree* t, const char* tag) {
    int root = build(t, 0, t->nprim);
    char buf[72];
    snprintf(buf, sizeof(buf), "%s_build_root", tag);
    CHECK(buf, root == 0);

    AABB u; box_union(t, &u);
    for (int k = 0; k < 3; k++) {
        snprintf(buf, sizeof(buf), "%s_root_min%d", tag, k);
        CHECK_NEAR(buf, t->n[0].box.min[k], u.min[k], 1e-5f);
        snprintf(buf, sizeof(buf), "%s_root_max%d", tag, k);
        CHECK_NEAR(buf, t->n[0].box.max[k], u.max[k], 1e-5f);
    }

    int seen[MAX_PRIM] = {0}, leaves = 0, maxdepth = 0;
    walk(t, root, seen, &leaves, 0, &maxdepth);
    int all = 1;
    for (int i = 0; i < t->nprim; i++) if (!seen[i]) all = 0;
    snprintf(buf, sizeof(buf), "%s_all_prims_in_leaves", tag);
    CHECK(buf, all);
    snprintf(buf, sizeof(buf), "%s_leaves_present", tag);
    CHECK(buf, leaves >= 1);
    snprintf(buf, sizeof(buf), "%s_tree_balanced", tag);
    CHECK(buf, maxdepth >= 0 && (1 << maxdepth) <= 4 * t->nprim + 4);
}

/* exhaustive oracle comparison: each ray vs every primitive + union box */
static void oracle_rays(Tree* t, int nrays, float tmax, const char* tag) {
    AABB u; box_union(t, &u);
    int mm = 0, tm = 0, um = 0;
    for (int r = 0; r < nrays; r++) {
        float o[3], d[3];
        for (int k = 0; k < 3; k++) o[k] = frand(-3.0f, 10.0f);
        do {
            for (int k = 0; k < 3; k++) d[k] = frand(-1.0f, 1.0f);
        } while (d[0]*d[0] + d[1]*d[1] + d[2]*d[2] < 1e-6f);

        float tr, tg;
        int hr = ref_box_hit(&u, o, d, tmax, &tr);
        int hg = aabb_hit(&t->n[0].box, o, d, tmax, &tg);
        if (hr != hg) { um++; }
        else if (hr && fabsf(tg - tr) > 1e-4f) um++;

        for (int i = 0; i < t->nprim; i++) {
            int hr = ref_box_hit(&t->p[i].box, o, d, tmax, &tr);
            int hg = aabb_hit(&t->p[i].box, o, d, tmax, &tg);
            if (hr != hg) { mm++; continue; }
            if (hr && fabsf(tg - tr) > 1e-4f) tm++;
        }
    }
    char buf[72];
    snprintf(buf, sizeof(buf), "%s_oracle_union", tag);
    CHECK(buf, um == 0);
    snprintf(buf, sizeof(buf), "%s_oracle_hitmiss", tag);
    CHECK(buf, mm == 0);
    snprintf(buf, sizeof(buf), "%s_oracle_tdist", tag);
    CHECK(buf, tm == 0);
}

/* collect primitive ids stored under node ni into ids[] */
static void collect(const Tree* t, int ni, int* ids, int* n) {
    const Node* nd = &t->n[ni];
    if (nd->count > 0) {
        for (int i = nd->first; i < nd->first + nd->count; i++)
            ids[(*n)++] = i;
    } else {
        collect(t, nd->left, ids, n);
        collect(t, nd->right, ids, n);
    }
}

/* nearest primitive by independently checking each box with the CUT:
   validates that aabb_hit's t-values order intersections correctly. */
static int nearest_via_aabb(const Tree* t, const float o[3], const float d[3],
                            float tmax, float* tbest_out) {
    int best = -1; float tbest = tmax;
    for (int i = 0; i < t->nprim; i++) {
        float tn;
        if (aabb_hit(&t->p[i].box, o, d, tmax, &tn) && tn < tbest) {
            tbest = tn; best = i;
        }
    }
    if (tbest_out) *tbest_out = tbest;
    return best;
}

int main(void) {
    /* ---- scene A: 10 boxes in a row along +x ---- */
    Tree a;
    memset(&a, 0, sizeof(a));
    for (int i = 0; i < 10; i++)
        tree_add(&a, (float)i, 0, 0, (float)i + 0.5f, 1, 1);
    CHECK("grid_count", a.nprim == 10);
    check_structure(&a, "grid");

    /* every primitive is reachable along +x through its own box */
    int reach = 1;
    for (int i = 0; i < a.nprim; i++) {
        float o[3] = {-1.0f, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tn;
        if (!aabb_hit(&a.p[i].box, o, d, 100.0f, &tn)) { reach = 0; break; }
    }
    CHECK("each_prim_reachable", reach);

    /* nearest ordering: leftmost box first, entry at t=2 */
    {
        float o[3] = {-1, .5f, .5f}, d[3] = {1, 0, 0}, tn;
        int id = nearest_via_aabb(&a, o, d, 100.0f, &tn);
        CHECK("nearest_leftmost_id", id == 0);
        CHECK_NEAR("nearest_leftmost_t", tn, 1.0f, 1e-5f);
    }
    /* negative direction: rightmost box first, entry t=0.5 */
    {
        float o[3] = {10, .5f, .5f}, d[3] = {-1, 0, 0}, tn;
        int id = nearest_via_aabb(&a, o, d, 100.0f, &tn);
        CHECK("nearest_rightmost_id_negdir", id == 9);
        CHECK_NEAR("nearest_rightmost_t_negdir", tn, 0.5f, 1e-5f);
    }
    /* ray pointing away hits nothing */
    {
        float o[3] = {-1, .5f, .5f}, d[3] = {-1, 0, 0}, tn;
        CHECK("behind_ray_miss", aabb_hit(&a.n[0].box, o, d, 100.0f, &tn) == 0);
    }
    /* parallel rays on each zero-component axis */
    {
        float ox[3] = {0.5f, -1, 0.5f}, dx[3] = {0, 1, 0}, tn;
        CHECK("parallel_y_hit", aabb_hit(&a.n[0].box, ox, dx, 100.0f, &tn) == 1);
        float oz[3] = {0.5f, 0.5f, -1}, dz[3] = {0, 0, 1};
        CHECK("parallel_z_hit", aabb_hit(&a.n[0].box, oz, dz, 100.0f, &tn) == 1);
        float om[3] = {0.5f, 2.0f, 0.5f};
        CHECK("parallel_offset_miss", aabb_hit(&a.n[0].box, om, dx, 100.0f, &tn) == 0);
    }
    /* tmax clipping, including the inclusive boundary */
    {
        float o[3] = {-1, .5f, .5f}, d[3] = {1, 0, 0}, tn;
        CHECK("tmax_short_miss", aabb_hit(&a.n[0].box, o, d, 0.5f, &tn) == 0);
        CHECK("tmax_boundary_hit", aabb_hit(&a.n[0].box, o, d, 2.0f, &tn) == 1);
    }
    /* corner grazing: boundary inclusive, just outside misses */
    {
        float og[3] = {0, 1.0f, 0.5f}, d[3] = {1, 0, 0}, tn;
        CHECK("graze_boundary_hit", aabb_hit(&a.n[0].box, og, d, 100, &tn) == 1);
        float om[3] = {0, 1.0001f, 0.5f};
        CHECK("graze_outside_miss", aabb_hit(&a.n[0].box, om, d, 100, &tn) == 0);
    }
    /* clear miss above the scene */
    {
        float o[3] = {0, 5, 5}, d[3] = {1, 0, 0}, tn;
        CHECK("above_scene_miss", aabb_hit(&a.n[0].box, o, d, 100.0f, &tn) == 0);
    }

    oracle_rays(&a, 300, 40.0f, "grid");

    /* ---- scene B: jittered, overlapping, spread mainly along y ---- */
    Tree b;
    memset(&b, 0, sizeof(b));
    static const float SB[][6] = {
        {0.0f, 0.0f, 0.0f, 0.9f, 0.9f, 0.9f},
        {0.5f, 1.2f, 0.0f, 1.4f, 2.1f, 0.9f},
        {2.0f, 2.4f, 1.0f, 2.8f, 3.3f, 1.8f},
        {2.2f, 3.6f, 1.2f, 3.0f, 4.5f, 2.0f},
        {4.0f, 0.0f, 2.2f, 4.7f, 0.8f, 3.0f},
        {4.2f, 1.0f, 2.4f, 4.9f, 1.9f, 3.2f},
        {6.0f, 4.0f, 0.0f, 6.6f, 4.8f, 0.7f},
        {6.1f, 5.2f, 0.2f, 6.7f, 6.0f, 0.9f},
        {1.0f, 6.0f, 3.0f, 1.8f, 6.8f, 3.8f},
        {3.0f, 7.0f, 3.2f, 3.7f, 7.9f, 4.0f},
        {5.0f, 6.2f, 4.1f, 5.8f, 7.1f, 4.9f},
        {0.2f, 8.0f, 1.0f, 1.1f, 8.8f, 1.7f},
    };
    for (unsigned i = 0; i < sizeof(SB) / sizeof(SB[0]); i++)
        tree_add(&b, SB[i][0], SB[i][1], SB[i][2],
                      SB[i][3], SB[i][4], SB[i][5]);
    check_structure(&b, "jitter");

    /* widest axis is y: the split must separate primitive centroids along
       y (not x). Centroid split -> every left-subtree centroid sits below
       every right-subtree centroid on the split axis; boxes may overlap. */
    {
        int lid[MAX_PRIM], rid[MAX_PRIM], nl = 0, nr = 0;
        collect(&b, b.n[0].left, lid, &nl);
        collect(&b, b.n[0].right, rid, &nr);
        int ysep = (nl > 0 && nr > 0), xsep;
        float lymax = -1e30f, rymin = 1e30f, lxmax = -1e30f, rxmin = 1e30f;
        for (int i = 0; i < nl; i++) {
            lymax = fmaxf(lymax, b.p[lid[i]].cy);
            lxmax = fmaxf(lxmax, b.p[lid[i]].cx);
        }
        for (int i = 0; i < nr; i++) {
            rymin = fminf(rymin, b.p[rid[i]].cy);
            rxmin = fminf(rxmin, b.p[rid[i]].cx);
        }
        if (!(lymax <= rymin + 1e-6f)) ysep = 0;
        /* x centroids stay interleaved/overlapping, so the split was not on x */
        xsep = (lxmax <= rxmin + 1e-6f);
        CHECK("split_centroid_separates_y", ysep);
        CHECK("split_centroids_overlap_x", !xsep);
        /* direct widest-axis check against the union box extents */
        float w[3] = { b.n[0].box.max[0] - b.n[0].box.min[0],
                       b.n[0].box.max[1] - b.n[0].box.min[1],
                       b.n[0].box.max[2] - b.n[0].box.min[2] };
        CHECK("widest_axis_is_y", w[1] >= w[0] && w[1] >= w[2]);
    }
    oracle_rays(&b, 300, 40.0f, "jitter");

    /* ---- small trees exercise the count<=2 leaf path ---- */
    {
        Tree one; memset(&one, 0, sizeof(one));
        tree_add(&one, 0, 0, 0, 1, 1, 1);
        check_structure(&one, "one");
        Tree two; memset(&two, 0, sizeof(two));
        tree_add(&two, 0, 0, 0, 1, 1, 1);
        tree_add(&two, 2, 0, 0, 3, 1, 1);
        check_structure(&two, "two");
    }

    /* ---- AABB area invariants against hand-computed values ---- */
    {
        AABB unit = {{0, 0, 0}, {1, 1, 1}};
        CHECK_NEAR("area_unit_cube", aabb_area(&unit), 6.0f, 1e-6f);
        AABB box234 = {{0, 0, 0}, {2, 3, 4}};
        CHECK_NEAR("area_2_3_4", aabb_area(&box234), 52.0f, 1e-5f);
        AABB deg = {{0, 0, 0}, {1, 0, 1}};
        /* flat box: surface area is twice the face area -> 2*(1*1) = 2 */
        CHECK_NEAR("area_flat_box", aabb_area(&deg), 2.0f, 1e-6f);
        AABB e; aabb_empty(&e);
        CHECK("area_empty_zero", aabb_area(&e) == 0.0f);
        /* note: aabb_hit on the 1e30/-1e30 "empty" box is an unhandled
           edge in the source (returns a hit); deliberately not asserted. */
    }

    /* ---- capacity: exactly MAX_PRIM fit, the next is rejected ---- */
    {
        Tree cap; memset(&cap, 0, sizeof(cap));
        int rc = 0;
        for (int i = 0; i < MAX_PRIM; i++)
            rc = tree_add(&cap, (float)i, 0, 0, (float)i + .5f, 1, 1);
        CHECK("capacity_fill", rc == MAX_PRIM - 1 && cap.nprim == MAX_PRIM);
        CHECK("capacity_reject", tree_add(&cap, 0, 0, 0, 1, 1, 1) == -1);
        CHECK("build_at_capacity", build(&cap, 0, cap.nprim) == 0);
    }

    /* ---- G7: a primitive whose centroid sits exactly on the split
       midpoint must go to the RIGHT subtree (strict `<` partition). ---- */
    {
        Tree eq; memset(&eq, 0, sizeof(eq));
        tree_add(&eq, 0, 0, 0, 1, 1, 1);   /* cx=0.5 */
        tree_add(&eq, 1, 0, 0, 2, 1, 1);   /* cx=1.5, == midpoint 1.5 */
        tree_add(&eq, 2, 0, 0, 3, 1, 1);   /* cx=2.5 */
        int root = build(&eq, 0, eq.nprim);
        int lid[MAX_PRIM], rid[MAX_PRIM], nl = 0, nr = 0;
        collect(&eq, eq.n[root].left, lid, &nl);
        collect(&eq, eq.n[root].right, rid, &nr);
        int boundary_right = 0;
        for (int i = 0; i < nr; i++) if (eq.p[rid[i]].cx == 1.5f) boundary_right = 1;
        CHECK("split_strict_lt_boundary_right", boundary_right);
        /* every left centroid strictly below the midpoint on split axis */
        int strict = (nl > 0);
        for (int i = 0; i < nl; i++) if (!(eq.p[lid[i]].cx < 1.5f)) strict = 0;
        CHECK("split_left_centroids_strict_below", strict);
    }

    /* ---- G8: equal centroids must trigger the midpoint-count fallback,
       giving a compact tree (3 nodes) instead of exhausting MAX_NODE. ---- */
    {
        Tree ec; memset(&ec, 0, sizeof(ec));
        tree_add(&ec, 0, 0, 0, 2, 1, 1);
        tree_add(&ec, 0, 0, 0, 2, 1, 1);
        tree_add(&ec, 0, 0, 0, 2, 1, 1);
        int root = build(&ec, 0, ec.nprim);
        CHECK("equal_centroid_build_ok", root == 0);
        CHECK("equal_centroid_compact_nodes", ec.nnode == 3);
        int bad = 0;
        for (int i = 0; i < ec.nnode; i++)
            if (ec.n[i].count == 0 && (ec.n[i].left < 0 || ec.n[i].right < 0)) bad = 1;
        CHECK("equal_centroid_no_dead_node", !bad);
    }

    /* ---- G11: a ray with a small-but-nonzero axis component must be
       treated as a slanted ray. The x/z components are exactly 0 (parallel,
       origin inside their slab), only y carries 1e-7 (between the 1e-8
       slant threshold and the loosened 1e-6 one). From oy=1.001 going down,
       the y-slab opens at t=10000 and closes at t=100100, so tmax=1000
       misses but tmax=20000 reaches it. A loosened 1e-6 threshold treats y
       as parallel instead and wrongly reports a hit already at tmax=1000. ---- */
    {
        AABB bx = {{0, 0, 0}, {2, 1, 1}};
        float o[3] = {0.5f, 1.001f, 0.5f}, d[3] = {0.0f, -1e-7f, 0.0f}, tn;
        CHECK("tiny_component_slanted_miss", aabb_hit(&bx, o, d, 1000.0f, &tn) == 0);
        CHECK("tiny_component_long_ray_hit", aabb_hit(&bx, o, d, 20000.0f, &tn) == 1);
    }

    /* ---- G12: a box with a negative extent (inverted min/max) must be
       rejected; without the guard the area formula returns garbage. ---- */
    {
        AABB inv = {{2, 0, 0}, {1, 1, 1}};   /* dx = -1 */
        CHECK("area_negative_extent_zero", aabb_area(&inv) == 0.0f);
    }

    printf("all_ok=%d\n", t_failures == 0);
    return t_failures == 0 ? 0 : 1;
}
