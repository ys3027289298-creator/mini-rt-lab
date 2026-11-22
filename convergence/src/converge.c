
/* converge: 两种低差异序列 (Halton vs Sobol) 的收敛指数对比。
 *
 * 方法: 对每个序列、每个样本数 N, 用 Cranley-Patterson 随机移位重复
 * REPS 次, 以 RMSE 作为误差; 在 log-log 坐标下最小二乘拟合
 *     err = C * N^(-alpha)
 * 输出 alpha 及其 95% 置信区间、同样本下两序列的误差比例、拟合残差
 * 给出的误差上界 (误差 <= 拟合值 * exp(r_max)); 样本点不足或数据
 * 不可用时报告"不可拟合"及原因。
 *
 * 用法: ./converge [k]      k = 样本数档位个数 (默认 9, N = 2^4 .. 2^(3+k))
 *       ./converge --selftest
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define REPS 32
#define MAX_POINTS 16
#define MIN_FIT_POINTS 4

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- rng: splitmix64 ---------- */
static unsigned long long rng_state = 0x123456789abcdefULL;
static unsigned long long next_u64(void) {
    unsigned long long z = (rng_state += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static double next_double(void) { return (next_u64() >> 11) * (1.0 / 9007199254740992.0); }

/* ---------- 低差异序列: 返回 [0,1)^2 中第 i 个点 ---------- */
typedef void (*Seq)(unsigned i, double* x, double* y);

static double vdc(unsigned i, unsigned base) {
    double f = 1.0 / base, r = 0.0;
    while (i > 0) { r += f * (i % base); i /= base; f /= base; }
    return r;
}

static void halton_pt(unsigned i, double* x, double* y) { *x = vdc(i, 2); *y = vdc(i, 3); }

/* Sobol 2D: 两维都用 Gray 码 g(i)=i^(i>>1) 索引 (顺序必须一致, 否则
 * 维间配对错位、点退化到对角线上)。
 * dim0: v_k = 2^-k; dim1 用本原多项式 x+1, 方向数 m_k = (2*m_{k-1}) XOR
 * m_{k-1}, m_1 = 1 (即 1,3,5,15,17,51,85,255,...)。 */
static unsigned sobol_v[2][32];
static void sobol_init(void) {
    unsigned m = 1;
    for (int k = 1; k <= 32; k++) {
        sobol_v[0][k - 1] = 1u << (32 - k);
        sobol_v[1][k - 1] = m << (32 - k);
        m = (2u * m) ^ m;
    }
}
static void sobol_pt(unsigned i, double* x, double* y) {
    unsigned g = i ^ (i >> 1), v0 = 0, v1 = 0;
    for (int k = 0; k < 32; k++)
        if (g & (1u << k)) { v0 ^= sobol_v[0][k]; v1 ^= sobol_v[1][k]; }
    *x = (double)v0 / 4294967296.0;
    *y = (double)v1 / 4294967296.0;
}

/* ---------- 场景: 被积函数与精确值 ---------- */
typedef struct { const char* name; double (*f)(double, double); double exact; } Scene;

static double f_smooth(double x, double y) { return sin(M_PI * x) * sin(M_PI * y); }
static double f_disk(double x, double y) {
    double dx = x - 0.5, dy = y - 0.5;
    return (dx * dx + dy * dy < 0.16) ? 1.0 : 0.0;
}
static const Scene scenes[] = {
    { "smooth(sin*sin)", f_smooth, 4.0 / (M_PI * M_PI) },
    { "disk(r=0.4)",     f_disk,   M_PI * 0.16 },
};

/* ---------- 误差: REPS 次随机移位的 RMSE ---------- */
static double rmse_error(Seq seq, const Scene* sc, int n) {
    double se = 0.0;
    for (int r = 0; r < REPS; r++) {
        double sx = next_double(), sy = next_double(), sum = 0.0;
        for (int i = 0; i < n; i++) {
            double x, y;
            seq((unsigned)i, &x, &y);
            x += sx; if (x >= 1.0) x -= 1.0;
            y += sy; if (y >= 1.0) y -= 1.0;
            sum += sc->f(x, y);
        }
        double e = sum / n - sc->exact;
        se += e * e;
    }
    return sqrt(se / REPS);
}

/* ---------- 幂律拟合 ---------- */
typedef struct {
    int ok;
    char reason[64];
    int n;
    double alpha, alpha_lo, alpha_hi; /* 指数及 95% 置信区间 */
    double coef;                      /* C: err = C * N^-alpha */
    double resid_bound;               /* exp(max|log 残差|): 误差 <= 拟合值 * 此值 */
    double r2;
} Fit;

/* t_{0.975, df}, df = 1..30 查表, 之外渐近插值到 1.96 */
static double t975(int df) {
    static const double t[] = {
        0, 12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
        2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
        2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042
    };
    if (df < 1) return 0.0;
    if (df <= 30) return t[df];
    return 1.96 + (2.042 - 1.96) * 30.0 / df;
}

static Fit fit_power(const double* ns, const double* err, int n) {
    Fit ft; memset(&ft, 0, sizeof(ft)); ft.n = n;
    if (n < MIN_FIT_POINTS) {
        snprintf(ft.reason, sizeof(ft.reason), "样本点不足 (%d < %d)", n, MIN_FIT_POINTS);
        return ft;
    }
    double lx[MAX_POINTS], ly[MAX_POINTS], mx = 0.0, my = 0.0;
    for (int i = 0; i < n; i++) {
        if (err[i] <= 0.0) {
            snprintf(ft.reason, sizeof(ft.reason), "误差非正 (err[%d]=%.3g), 无法取对数", i, err[i]);
            return ft;
        }
        lx[i] = log(ns[i]); ly[i] = log(err[i]);
        mx += lx[i]; my += ly[i];
    }
    mx /= n; my /= n;
    double sxx = 0.0, sxy = 0.0, sst = 0.0;
    for (int i = 0; i < n; i++) {
        sxx += (lx[i] - mx) * (lx[i] - mx);
        sxy += (lx[i] - mx) * (ly[i] - my);
        sst += (ly[i] - my) * (ly[i] - my);
    }
    if (sxx <= 0.0) {
        snprintf(ft.reason, sizeof(ft.reason), "样本数无变化 (log N 方差为 0)");
        return ft;
    }
    double b = sxy / sxx, a = my - b * mx;
    double sse = 0.0, rmax = 0.0;
    for (int i = 0; i < n; i++) {
        double r = ly[i] - (a + b * lx[i]);
        sse += r * r;
        if (fabs(r) > rmax) rmax = fabs(r);
    }
    double se_b = (n > 2) ? sqrt(sse / (n - 2) / sxx) : 0.0;
    double tc = t975(n - 2);
    ft.ok = 1;
    ft.alpha = -b;
    ft.alpha_lo = -b - tc * se_b;
    ft.alpha_hi = -b + tc * se_b;
    ft.coef = exp(a);
    ft.resid_bound = exp(rmax);
    ft.r2 = (sst > 0.0) ? 1.0 - sse / sst : 1.0;
    return ft;
}

static double fit_eval(const Fit* ft, double n) { return ft->coef * pow(n, -ft->alpha); }

/* ---------- 报告 ---------- */
static void run_report(int k) {
    double ns[MAX_POINTS];
    for (int i = 0; i < k; i++) ns[i] = (double)(1 << (4 + i));
    const Seq seqs[2] = { halton_pt, sobol_pt };
    const char* names[2] = { "halton", "sobol" };

    for (unsigned s = 0; s < sizeof(scenes) / sizeof(scenes[0]); s++) {
        const Scene* sc = &scenes[s];
        printf("场景 %s, 精确值 = %.6f, 每点 %d 次随机移位取 RMSE\n", sc->name, sc->exact, REPS);
        Fit fits[2];
        for (int q = 0; q < 2; q++) {
            double err[MAX_POINTS];
            printf("  序列 %s:\n", names[q]);
            for (int i = 0; i < k; i++) {
                err[i] = rmse_error(seqs[q], sc, (int)ns[i]);
                printf("    N=%-6.0f rmse=%.6g\n", ns[i], err[i]);
            }
            fits[q] = fit_power(ns, err, k);
            if (!fits[q].ok) {
                printf("    不可拟合: %s\n", fits[q].reason);
                continue;
            }
            printf("    指数 alpha = %.3f, 95%% CI [%.3f, %.3f]\n",
                   fits[q].alpha, fits[q].alpha_lo, fits[q].alpha_hi);
            printf("    残差上界: 误差 <= 拟合值 x %.3f, R^2 = %.4f\n",
                   fits[q].resid_bound, fits[q].r2);
        }
        if (fits[0].ok && fits[1].ok) {
            double n_max = ns[k - 1];
            double e0 = fit_eval(&fits[0], n_max), e1 = fit_eval(&fits[1], n_max);
            int better = e1 < e0 ? 1 : 0;
            double ratio = better ? e0 / e1 : e1 / e0;
            printf("  对比 (N=%.0f): %s 误差更低, 为另一序列的 1/%.2f (误差比 %.2f:1);",
                   n_max, names[better], ratio, ratio);
            printf(" 指数 %s alpha=%.3f %s alpha=%.3f\n",
                   names[0], fits[0].alpha, names[1], fits[1].alpha);
            if (fits[0].alpha_lo < fits[1].alpha_hi && fits[1].alpha_lo < fits[0].alpha_hi)
                printf("  注意: 两序列指数的 95%% 置信区间重叠, 指数差异不显著\n");
        } else {
            printf("  对比不可用: 至少一个序列不可拟合\n");
        }
        printf("\n");
    }
}

/* ---------- 自检 ---------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int selftest(void) {
    sobol_init();
    double x, y;
    sobol_pt(1, &x, &y);
    check("sobol_pt1", fabs(x - 0.5) < 1e-12 && fabs(y - 0.5) < 1e-12);
    sobol_pt(2, &x, &y);
    check("sobol_pt2", fabs(x - 0.75) < 1e-12 && fabs(y - 0.25) < 1e-12);
    /* (0,2)-序列性质: 前 4 个点必须落在 4 个不同象限 */
    int quad = 0;
    for (unsigned i = 0; i < 4; i++) {
        sobol_pt(i, &x, &y);
        quad |= 1 << ((x >= 0.5) + 2 * (y >= 0.5));
    }
    check("sobol_quadrants", quad == 0xF);
    halton_pt(1, &x, &y);
    check("halton_pt1", fabs(x - 0.5) < 1e-12 && fabs(y - 1.0 / 3.0) < 1e-12);

    /* 合成数据 err = 2*N^-0.75, 应精确恢复 alpha = 0.75 */
    double ns[6], err[6];
    for (int i = 0; i < 6; i++) { ns[i] = (double)(1 << (4 + i)); err[i] = 2.0 * pow(ns[i], -0.75); }
    Fit ft = fit_power(ns, err, 6);
    check("fit_ok", ft.ok);
    check("fit_alpha", fabs(ft.alpha - 0.75) < 1e-9);
    check("fit_ci_contains", ft.alpha_lo <= 0.75 && ft.alpha_hi >= 0.75);
    check("fit_resid_small", ft.resid_bound < 1.0 + 1e-6);

    Fit few = fit_power(ns, err, 3);
    check("fit_few_not_ok", !few.ok);
    double bad[4] = { 1.0, 0.5, 0.0, 0.3 };
    Fit zero = fit_power(ns, bad, 4);
    check("fit_zero_not_ok", !zero.ok);
    double same_n[4] = { 64, 64, 64, 64 };
    Fit flat = fit_power(same_n, err, 4);
    check("fit_flat_not_ok", !flat.ok);

    /* 实际序列: sobol 在光滑场景下应随 N 收敛 */
    double e_lo = rmse_error(sobol_pt, &scenes[0], 64);
    double e_hi = rmse_error(sobol_pt, &scenes[0], 4096);
    check("sobol_converges", e_hi < e_lo);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    sobol_init();
    if (argc > 1 && strcmp(argv[1], "--selftest") == 0) return selftest();
    int k = (argc > 1) ? atoi(argv[1]) : 9;
    if (k < 1) k = 1;
    if (k > MAX_POINTS) k = MAX_POINTS;
    run_report(k);
    return 0;
}
