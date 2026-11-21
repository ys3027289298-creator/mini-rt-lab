/* extreme_material.c — material 模块极端场景矩阵。
   直接 #include 模块翻译单元（main 重命名），仓库源文件不做任何修改。 */
#define main material_repo_selftest
#include "../material/src/material.c"
#undef main
#include "harness.h"

/* ---------- baseline：人工构造基准用例 ---------- */

static void c_repo_selftest(void) {
    EXPECT_TRUE("repo_selftest_rc", material_repo_selftest() == 0);
}

static void c_ggx_baseline(void) {
    /* 手算：rough=0.5 -> a2=0.0625, d=0.0625, ndf=1/(pi*0.0625)=5.0930 */
    float ndf = ggx_ndf(1.0f, 0.5f);
    EXPECT_TRUE("ndf_value", fabsf(ndf - 5.0930f) < 1e-2f);
    EXPECT_FINITE("ndf", ndf);
}

/* ---------- zero_area：零面积光源（delta 分布 / 掠射 pdf=0） ----------
   触发条件：roughness = 0（完美镜面 = 点/零面积光源的 BRDF 极限）；
   或 rough^4 < FLT_EPSILON/2（rough < ~0.016）：float 下 1-(1-a2) 舍入为 0，
   d 坍塌为 0 -> 除零（a2>0 得 Inf，a2 下溢为 0 得 NaN）；
   cos_theta = 0（掠射方向 pdf 为零）。 */

static void c_ggx_rough_zero(void) {
    float ndf = ggx_ndf(1.0f, 0.0f); /* d = 0 -> 0/0 */
    EXPECT_FINITE("ndf_rough0", ndf);
}

static void c_ggx_rough_0p01(void) {
    /* rough=0.01：常见材质粗糙度，a2=1e-8 < FLT_EPSILON -> d 坍塌为 0 */
    float ndf = ggx_ndf(1.0f, 0.01f);
    EXPECT_FINITE("ndf_rough0p01", ndf);
}

static void c_ggx_rough_underflow_nan(void) {
    /* rough=1e-20：a2=rough^4=1e-80 下溢为 0 -> 0/0 */
    float ndf = ggx_ndf(1.0f, 1e-20f);
    EXPECT_FINITE("ndf_rough1e-20", ndf);
}

static void c_ggx_rough_small_ok(void) {
    /* rough=0.1：a2=1e-4 > FLT_EPSILON，d 不坍塌，应有限 */
    float ndf = ggx_ndf(1.0f, 0.1f);
    EXPECT_FINITE("ndf_rough0p1", ndf);
}

static void c_cosine_pdf_grazing(void) {
    float p = cosine_pdf(0.0f); /* 掠射：pdf=0，下游 1/pdf 即 inf */
    EXPECT_FINITE("pdf", p);
    EXPECT_TRUE("pdf_zero", p == 0.0f);
}

static void c_cosine_sample_edges(void) {
    static const float us[][2] = {
        {0.0f, 0.0f}, {0.0f, 1.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.5f, 0.5f},
    };
    for (int i = 0; i < 5; i++) {
        V3 v = cosine_sample(us[i][0], us[i][1]);
        EXPECT_FINITE("x", v.x);
        EXPECT_FINITE("y", v.y);
        EXPECT_FINITE("z", v.z);
        float len = sqrtf(dot3(v, v));
        EXPECT_TRUE("unit_or_zero", len < 1e-4f || fabsf(len - 1.0f) < 1e-4f);
    }
}

static void c_cosine_sample_out_of_domain(void) {
    /* u1<0：sqrtf 得 NaN，norm3 的 l>0 守卫应吞掉并返回零向量 */
    V3 v = cosine_sample(-0.5f, 0.25f);
    EXPECT_FINITE("x", v.x);
    EXPECT_FINITE("y", v.y);
    EXPECT_FINITE("z", v.z);
}

static void c_norm_zero_vector(void) {
    V3 z = {0, 0, 0};
    V3 n = norm3(z); /* 零方向（零面积光源退化方向） */
    EXPECT_FINITE("x", n.x);
    EXPECT_FINITE("y", n.y);
    EXPECT_FINITE("z", n.z);
}

/* ---------- huge_coords：超大量级参数 ---------- */

static void c_lambert_huge_albedo(void) {
    float b = brdf_lambert(1e38f);
    EXPECT_FINITE("brdf", b);
}

static void c_ggx_grazing_nh(void) {
    float ndf = ggx_ndf(0.0f, 0.5f); /* nh=0：d=1，正常 */
    EXPECT_FINITE("ndf", ndf);
    EXPECT_TRUE("ndf_value", fabsf(ndf - 0.0625f / 3.14159265f) < 1e-4f);
}

int main(void) {
    static const Case cases[] = {
        {"baseline",    "repo_selftest",            c_repo_selftest},
        {"baseline",    "ggx_value",                c_ggx_baseline},
        {"zero_area",   "ggx_rough_zero",           c_ggx_rough_zero},
        {"zero_area",   "ggx_rough_0p01",           c_ggx_rough_0p01},
        {"zero_area",   "ggx_rough_underflow_nan",  c_ggx_rough_underflow_nan},
        {"zero_area",   "ggx_rough_small_ok",       c_ggx_rough_small_ok},
        {"zero_area",   "cosine_pdf_grazing",       c_cosine_pdf_grazing},
        {"zero_area",   "cosine_sample_edges",      c_cosine_sample_edges},
        {"zero_area",   "cosine_sample_out_domain", c_cosine_sample_out_of_domain},
        {"zero_area",   "norm_zero_vector",         c_norm_zero_vector},
        {"huge_coords", "lambert_huge_albedo",      c_lambert_huge_albedo},
        {"huge_coords", "ggx_grazing_nh",           c_ggx_grazing_nh},
    };
    return run_matrix("material", cases, (int)(sizeof cases / sizeof cases[0]));
}
