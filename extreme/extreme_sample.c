/* extreme_sample.c — sample 模块极端场景矩阵。
   直接 #include 模块翻译单元（main 重命名），仓库源文件不做任何修改。 */
#define main sample_repo_selftest
#include "../sample/src/sample.c"
#undef main
#include "harness.h"

static double f_huge(double u) { (void)u; return 1e308; }

/* ---------- baseline：人工构造基准用例 ---------- */

static void c_repo_selftest(void) {
    EXPECT_TRUE("repo_selftest_rc", sample_repo_selftest() == 0);
}

static void c_quad_mean_baseline(void) {
    Stat s = stratify(64, f_quad);
    EXPECT_TRUE("mean_third", fabs(s.mean - 1.0 / 3.0) < 1e-3);
    EXPECT_FINITE("mean", s.mean);
    EXPECT_FINITE("var", s.var);
}

static void c_strata_1(void) {
    Stat s = stratify(1, f_quad); /* 单样本：u=0.5 -> 0.25 */
    EXPECT_TRUE("mean_quarter", fabs(s.mean - 0.25) < 1e-15);
    EXPECT_TRUE("var_zero", s.var == 0.0);
    EXPECT_FINITE("mean", s.mean);
}

/* ---------- zero_area：零面积光源（采样数为零） ----------
   触发条件：strata = 0（光源上一个样本都没有）。 */

static void c_zero_strata(void) {
    Stat s = stratify(0, f_quad);
    /* sum/n = 0/0：均值 NaN */
    EXPECT_FINITE("mean", s.mean);
    EXPECT_FINITE("var", s.var);
}

/* ---------- huge_coords：超大量级被积函数 ----------
   触发条件：积分估计值超出 double 范围（f >= 1e308，累加溢出）。 */

static void c_exploding_integrand(void) {
    Stat s = stratify(64, f_huge);
    /* sum 溢出 -> mean = inf；var 的 NaN 被 var>0 守卫吞掉 */
    EXPECT_FINITE("mean", s.mean);
    EXPECT_FINITE("var", s.var);
}

/* ---------- 低差异序列极值索引 ---------- */

static void c_halton_max_index(void) {
    double a = halton2(4294967295u, 2);
    EXPECT_FINITE("halton2_max", a);
    EXPECT_TRUE("in_unit", a >= 0.0 && a < 1.0);
    double b = halton(4294967295u, 7); /* radical_inverse 溢出回绕路径 */
    EXPECT_FINITE("halton_max", b);
    double c = halton2(4294967295u, 13);
    EXPECT_FINITE("halton2_max_b13", c);
    EXPECT_TRUE("in_unit_b13", c >= 0.0 && c < 1.0);
}

int main(void) {
    static const Case cases[] = {
        {"baseline",    "repo_selftest",        c_repo_selftest},
        {"baseline",    "quad_mean",            c_quad_mean_baseline},
        {"baseline",    "strata_1",             c_strata_1},
        {"zero_area",   "zero_strata",          c_zero_strata},
        {"huge_coords", "exploding_integrand",  c_exploding_integrand},
        {"huge_coords", "halton_max_index",     c_halton_max_index},
    };
    return run_matrix("sample", cases, (int)(sizeof cases / sizeof cases[0]));
}
