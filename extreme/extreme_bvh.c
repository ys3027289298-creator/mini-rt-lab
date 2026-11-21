/* extreme_bvh.c — bvh 模块极端场景矩阵。
   直接 #include 模块翻译单元（main 重命名），测试的就是仓库原始代码，
   仓库源文件不做任何修改。 */
#define main bvh_repo_selftest
#include "../bvh/src/bvh.c"
#undef main
#include "harness.h"

/* ---------- baseline：人工构造基准用例（已知期望值，用于比对） ---------- */

static void c_repo_selftest(void) {
    EXPECT_TRUE("repo_selftest_rc", bvh_repo_selftest() == 0);
}

static void c_unit_box_hit(void) {
    AABB b = {{0, 0, 0}, {1, 1, 1}};
    float o[3] = {-1, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 100.0f, &tm);
    EXPECT_TRUE("hit", hit == 1);
    EXPECT_TRUE("tmin_exact", fabsf(tm - 1.0f) < 1e-6f);
    EXPECT_FINITE("tmin", tm);
}

/* ---------- long_thin：超长细三角形 ----------
   触发条件：包围盒边长比 >= 1e12（长边 1e6，薄边 1e-6）；
   或薄边尺寸低于所在坐标量级的 float 相对精度（~1e-7）。 */

static void c_long_thin_build_hit(void) {
    Tree t; memset(&t, 0, sizeof t);
    tree_add(&t, 0, 0, 0, 1e6f, 1e-6f, 1e-6f);  /* 边长比 1e12 */
    tree_add(&t, 0, 0, 0, 1e-6f, 1e6f, 1e-6f);
    tree_add(&t, 0, 0, 0, 1e-6f, 1e-6f, 1e6f);
    int root = build(&t, 0, t.nprim);
    EXPECT_TRUE("build_ok", root == 0);
    EXPECT_FINITE("root_area", aabb_area(&t.n[0].box));

    float o[3] = {-1, 0.5e-6f, 0.5e-6f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&t.p[0].box, o, d, 1e9f, &tm);
    EXPECT_TRUE("hit_long_axis", hit == 1);
    EXPECT_FINITE("tmin", tm);

    float o2[3] = {-1, 2e-6f, 0.5e-6f}; /* 越过薄边必须 miss */
    EXPECT_TRUE("miss_thin", aabb_hit(&t.p[0].box, o2, d, 1e9f, &tm) == 0);
}

static void c_long_thin_collapse(void) {
    /* 1e-6 的薄边在 1e6 坐标量级下被 float 舍入吞掉（eps(1e6)≈0.0625） */
    Tree t; memset(&t, 0, sizeof t);
    tree_add(&t, 1e6f, 1e6f, 1e6f, 1e6f + 1.0f, 1e6f + 1e-6f, 1e6f + 1e-6f);
    EXPECT_TRUE("thin_dims_collapsed",
                t.p[0].box.min[1] == t.p[0].box.max[1]); /* 记录坍塌事实 */
    int root = build(&t, 0, t.nprim);
    EXPECT_TRUE("build_ok", root == 0);
    EXPECT_FINITE("centroid", t.p[0].cx);
    EXPECT_FINITE("area", aabb_area(&t.p[0].box));
    float o[3] = {0, 1e6f, 1e6f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&t.p[0].box, o, d, 1e9f, &tm);
    EXPECT_TRUE("hit", hit == 1);
    EXPECT_FINITE("tmin", tm);
}

/* ---------- near_far：极近 / 极远平面 ----------
   触发条件：tmax <= 1e-30（极近）或 tmax >= 1e30（极远）；
   射线原点距图元 >= 1e25。 */

static void c_tiny_tmax(void) {
    AABB b = {{0, 0, 0}, {1, 1, 1}};
    float o[3] = {-1, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 1e-30f, &tm);
    EXPECT_TRUE("no_hit_beyond_near", hit == 0);
    EXPECT_FINITE("tmin", tm);
}

static void c_zero_tmax(void) {
    AABB b = {{0, 0, 0}, {1, 1, 1}};
    float o[3] = {0.5f, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 0.0f, &tm); /* 近平面贴到原点 */
    EXPECT_TRUE("hit_inside", hit == 1);
    EXPECT_TRUE("tmin_zero", tm == 0.0f);
}

static void c_huge_tmax_far_origin(void) {
    AABB b = {{0, 0, 0}, {1, 1, 1}};
    float o[3] = {-1e25f, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 1e30f, &tm);
    EXPECT_TRUE("hit", hit == 1);
    EXPECT_TRUE("tmin_value", fabsf(tm - 1e25f) / 1e25f < 1e-6f);
    EXPECT_FINITE("tmin", tm);
}

static void c_origin_on_plane(void) {
    AABB b = {{0, 0, 0}, {1, 1, 1}};
    float o[3] = {0, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 100.0f, &tm);
    EXPECT_TRUE("hit", hit == 1);
    EXPECT_TRUE("tmin_zero", tm == 0.0f);
}

/* ---------- zero_area：零面积光源 / 退化图元 ----------
   触发条件：包围盒至少一维 min==max（面积或体积为零）。 */

static void c_point_box(void) {
    AABB b = {{1, 1, 1}, {1, 1, 1}};
    EXPECT_TRUE("area_zero", aabb_area(&b) == 0.0f);
    float o[3] = {0, 1, 1}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 100.0f, &tm);
    EXPECT_TRUE("hit_point", hit == 1);
    EXPECT_TRUE("tmin_one", fabsf(tm - 1.0f) < 1e-6f);
    EXPECT_FINITE("tmin", tm);
}

static void c_flat_box(void) {
    AABB b = {{0, 0, 0}, {1, 1, 0}}; /* 零厚度平板 */
    EXPECT_TRUE("area_finite", aabb_area(&b) == 2.0f);
    float o[3] = {0.5f, 0.5f, -1}, d[3] = {0, 0, 1}, tm = -1;
    int hit = aabb_hit(&b, o, d, 100.0f, &tm);
    EXPECT_TRUE("hit_flat", hit == 1);
    EXPECT_FINITE("tmin", tm);
}

static void c_identical_prims_build(void) {
    /* 32 个质心完全相同的点图元：分裂必然落空，走 count/2 兜底，
       验证不会无限递归 / 节点数组不越界 */
    Tree t; memset(&t, 0, sizeof t);
    for (int i = 0; i < MAX_PRIM; i++) tree_add(&t, 1, 1, 1, 1, 1, 1);
    int root = build(&t, 0, t.nprim);
    EXPECT_TRUE("build_terminates", root == 0);
    EXPECT_TRUE("node_budget", t.nnode <= MAX_NODE);
}

/* ---------- huge_coords：超大坐标 ----------
   触发条件：|坐标| >= 1e18（float 精度坍塌）；
   坐标 >= 1.7e38（质心 min+max 溢出）；边长 >= 1e19（面积溢出）。 */

static void c_coords_1e18(void) {
    Tree t; memset(&t, 0, sizeof t);
    tree_add(&t, 1e18f, 0, 0, 1e18f + 1, 1, 1); /* +1 被舍入吞掉 */
    int root = build(&t, 0, t.nprim);
    EXPECT_TRUE("build_ok", root == 0);
    EXPECT_FINITE("centroid", t.p[0].cx);
    float o[3] = {0, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&t.p[0].box, o, d, 1e30f, &tm);
    EXPECT_TRUE("hit", hit == 1);
    EXPECT_FINITE("tmin", tm);
}

static void c_coords_near_fltmax(void) {
    Tree t; memset(&t, 0, sizeof t);
    tree_add(&t, 3e38f, 3e38f, 3e38f, 3e38f, 3e38f, 3e38f);
    build(&t, 0, t.nprim);
    /* (3e38 + 3e38) * 0.5 -> 6e38 > FLT_MAX：质心计算溢出 */
    EXPECT_FINITE("centroid_x", t.p[0].cx);
    EXPECT_FINITE("centroid_y", t.p[0].cy);
    EXPECT_FINITE("centroid_z", t.p[0].cz);
}

static void c_area_overflow(void) {
    AABB b = {{0, 0, 0}, {1e20f, 1e20f, 1e20f}};
    /* 1e20 * 1e20 = 1e40 > FLT_MAX：表面积溢出 */
    EXPECT_FINITE("area", aabb_area(&b));
}

static void c_opposite_extremes(void) {
    /* 原点在 -3e38、盒在 +3e38：跨度 6e38 溢出为 inf，
       slab 测试退化为 t0=inf > tmax -> miss（结果凑巧正确，记录行为） */
    AABB b = {{3e38f, 0, 0}, {3.3e38f, 1, 1}};
    float o[3] = {-3e38f, 0.5f, 0.5f}, d[3] = {1, 0, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 1e30f, &tm);
    EXPECT_TRUE("miss_beyond_tmax", hit == 0);
    EXPECT_FINITE("tmin", tm);
}

static void c_denormal_dir(void) {
    AABB b = {{0, 0, 0}, {1, 1, 1}};
    float o[3] = {0.5f, 0.5f, 0.5f}, d[3] = {1e-40f, 1, 0}, tm = -1;
    int hit = aabb_hit(&b, o, d, 100.0f, &tm); /* 次法向方向分量按平行处理 */
    EXPECT_TRUE("hit_inside", hit == 1);
    EXPECT_FINITE("tmin", tm);
}

int main(void) {
    static const Case cases[] = {
        {"baseline",    "repo_selftest",        c_repo_selftest},
        {"baseline",    "unit_box_hit",         c_unit_box_hit},
        {"long_thin",   "build_and_hit",        c_long_thin_build_hit},
        {"long_thin",   "rounding_collapse",    c_long_thin_collapse},
        {"near_far",    "tiny_tmax",            c_tiny_tmax},
        {"near_far",    "zero_tmax",            c_zero_tmax},
        {"near_far",    "huge_tmax_far_origin", c_huge_tmax_far_origin},
        {"near_far",    "origin_on_plane",      c_origin_on_plane},
        {"zero_area",   "point_box",            c_point_box},
        {"zero_area",   "flat_box",             c_flat_box},
        {"zero_area",   "identical_prims_build",c_identical_prims_build},
        {"huge_coords", "coords_1e18",          c_coords_1e18},
        {"huge_coords", "near_fltmax",          c_coords_near_fltmax},
        {"huge_coords", "area_overflow",        c_area_overflow},
        {"huge_coords", "opposite_extremes",    c_opposite_extremes},
        {"huge_coords", "denormal_dir",         c_denormal_dir},
    };
    return run_matrix("bvh", cases, (int)(sizeof cases / sizeof cases[0]));
}
