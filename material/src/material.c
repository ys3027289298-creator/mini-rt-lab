
/* material: Lambert + GGX BRDF evaluation and cosine-weighted sampling,
   with a PDF check against the analytic cosine distribution.
   Also: dielectric interface Fresnel (glass/water), rough-interface
   perturbation with energy compensation, and conservation self-checks. */
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct { float x, y, z; } V3;

static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}

/* Lambert: albedo / pi */
static float brdf_lambert(float albedo) { return albedo / 3.14159265f; }

/* GGX normal distribution */
static float ggx_ndf(float nh, float rough) {
    float a = rough * rough;
    float a2 = a * a;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / (3.14159265f * d * d);
}

/* cosine-weighted hemisphere: pdf = cos/pi */
static float cosine_pdf(float cos_theta) { return cos_theta / 3.14159265f; }

static V3 cosine_sample(float u1, float u2) {
    float r = sqrtf(u1);
    float phi = 2.0f * 3.14159265f * u2;
    V3 v = { r * cosf(phi), r * sinf(phi), sqrtf(1.0f - u1 > 0 ? 1.0f - u1 : 0) };
    return norm3(v);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---- 电介质交界:解析菲涅耳 + 粗糙扰动补偿 ----
   玻璃-水交界过亮的常见原因:反射/透射未按菲涅耳加权(R、T 各算各的
   或直接相加超过 1),或透射侧漏乘 (nt/ni)^2 的辐射度缩放。这里的 R/T
   是吞吐量权重,守恒条件为 R+T=1(界面无吸收)。 */
#define ETA_GLASS 1.50
#define ETA_WATER 1.33

/* 解析菲涅耳反射率(无吸收、非偏振)。ci = |cos theta_i|,全反射返回 1。 */
static double fresnel_R_exact(double ci, double ni, double nt) {
    if (ci <= 0.0) return 1.0;
    double s2 = ni * ni * (1.0 - ci * ci) / (nt * nt); /* sin^2 theta_t */
    if (s2 >= 1.0) return 1.0;                        /* 全反射 */
    double ct = sqrt(1.0 - s2);
    double rs = (ni * ci - nt * ct) / (ni * ci + nt * ct);
    double rp = (nt * ci - ni * ct) / (nt * ci + ni * ct);
    return 0.5 * (rs * rs + rp * rp);
}

/* 透射率的独立解析式(振幅系数 + 折射波束截面比),不取 1-R,
   用于验证能量守恒解析条件: Rs + Ts = 1, Rp + Tp = 1(逐偏振成立)。 */
static double fresnel_T_exact(double ci, double ni, double nt) {
    if (ci <= 0.0) return 0.0;
    double s2 = ni * ni * (1.0 - ci * ci) / (nt * nt);
    if (s2 >= 1.0) return 0.0;                        /* 全反射 */
    double ct = sqrt(1.0 - s2);
    double ts = 2.0 * ni * ci / (ni * ci + nt * ct);
    double tp = 2.0 * ni * ci / (nt * ci + ni * ct);
    double beam = nt * ct / (ni * ci);                /* 波束截面比 */
    return 0.5 * beam * (ts * ts + tp * tp);
}

#define ROUGH_STRATA 24  /* 24x24 确定性分层求积,无随机数 */

/* 粗糙交界:GGX 微法线扰动下对菲涅耳取平均(逆 CDF 采样,pdf ∝ D(m)(m·n))。
   单次散射模型中 wi·m <= 0 的背向微面不贡献,可见权重 W = P(wi·m > 0) < 1,
   不补偿时 R、T 同被低估一个因子 W(掠射角最明显,W 最低约 0.5)。
   补偿项:乘以 1/W 重归一化,使 R+T = 1 在全角度成立;
   visible_w 返回 W,补偿量即 (1-W)。
   rough == 0:微法线分布退化为 delta,直接用法线 —— 退化为镜面。 */
static double fresnel_avg_rough(double ci, double ni, double nt, double rough,
                                int want_transmission, double* visible_w) {
    double si = sqrt(1.0 - ci * ci > 0.0 ? 1.0 - ci * ci : 0.0);
    double a2 = rough * rough;
    double acc = 0.0;
    int vis = 0, tot = 0;
    for (int i = 0; i < ROUGH_STRATA; i++) {
        for (int j = 0; j < ROUGH_STRATA; j++) {
            double u1 = (i + 0.5) / ROUGH_STRATA;
            double u2 = (j + 0.5) / ROUGH_STRATA;
            double cm, sm, phi;
            if (rough <= 0.0) { cm = 1.0; sm = 0.0; phi = 0.0; } /* 退化:镜面 */
            else {
                cm = sqrt((1.0 - u1) / (1.0 + (a2 - 1.0) * u1));
                sm = sqrt(1.0 - cm * cm > 0.0 ? 1.0 - cm * cm : 0.0);
                phi = 2.0 * 3.14159265358979323846 * u2;
            }
            double cim = si * sm * cos(phi) + ci * cm; /* wi·m,wi 在 xz 平面 */
            tot++;
            if (cim <= 0.0) continue; /* 背向微面:单次散射不贡献 */
            vis++;
            acc += want_transmission ? fresnel_T_exact(cim, ni, nt)
                                     : fresnel_R_exact(cim, ni, nt);
        }
    }
    if (visible_w) *visible_w = (double)vis / tot;
    if (vis == 0) return want_transmission ? 0.0 : 1.0;
    return acc / vis; /* 含 1/W 补偿 */
}

int main(void) {
    check("lambert_unit", fabsf(brdf_lambert(3.14159265f) - 1.0f) < 1e-4f);
    check("lambert_zero", brdf_lambert(0.0f) == 0.0f);

    check("ggx_peak", ggx_ndf(1.0f, 0.5f) > ggx_ndf(0.2f, 0.5f));
    check("ggx_positive", ggx_ndf(0.5f, 0.3f) > 0.0f);

    check("pdf_half", fabsf(cosine_pdf(0.5f) - 0.5f / 3.14159265f) < 1e-6f);

    V3 s = cosine_sample(0.5f, 0.25f);
    check("sample_unit", fabsf(sqrtf(dot3(s, s)) - 1.0f) < 1e-4f);
    check("sample_z_nonneg", s.z >= 0.0f);

    V3 up = { 0, 0, 1 };
    check("sample_hemisphere", dot3(s, up) >= -1e-6f);

    /* mean of cos over samples ~ 2/3 for cosine sampling */
    double acc = 0.0; int n = 512;
    for (int i = 0; i < n; i++) {
        float u1 = (float)((i % 16) + 0.5f) / 16.0f;
        float u2 = (float)((i / 16) + 0.5f) / 32.0f;
        V3 v = cosine_sample(u1, u2);
        acc += v.z;
    }
    check("mean_cos_close", fabs(acc / n - 2.0 / 3.0) < 0.02);

    /* ---- 电介质交界:菲涅耳能量守恒与粗糙补偿 ---- */

    /* 解析守恒条件:R+T=1,T 由独立振幅公式算出(不取 1-R);
       玻璃->水 与 水->玻璃(含全反射区)两个方向 */
    double max_cons = 0.0;
    for (int k = 1; k <= 180; k++) {
        double ci = k / 180.0;
        double e1 = fresnel_R_exact(ci, ETA_GLASS, ETA_WATER) + fresnel_T_exact(ci, ETA_GLASS, ETA_WATER);
        double e2 = fresnel_R_exact(ci, ETA_WATER, ETA_GLASS) + fresnel_T_exact(ci, ETA_WATER, ETA_GLASS);
        if (fabs(e1 - 1.0) > max_cons) max_cons = fabs(e1 - 1.0);
        if (fabs(e2 - 1.0) > max_cons) max_cons = fabs(e2 - 1.0);
    }
    printf("report: exact_RplusT_max_dev=%.3e\n", max_cons);
    check("fresnel_energy_exact", max_cons < 1e-9);

    /* 正入射:R = ((ni-nt)/(ni+nt))^2,玻璃->水约 0.36% */
    double r0 = (ETA_GLASS - ETA_WATER) / (ETA_GLASS + ETA_WATER);
    r0 *= r0;
    check("fresnel_normal_glass_water",
          fabs(fresnel_R_exact(1.0, ETA_GLASS, ETA_WATER) - r0) < 1e-12);

    /* 玻璃->水(密到疏),临界角约 62.5 度,之外全反射:R=1, T=0 */
    check("fresnel_tir", fresnel_R_exact(0.3, ETA_GLASS, ETA_WATER) == 1.0 &&
                         fresnel_T_exact(0.3, ETA_GLASS, ETA_WATER) == 0.0);

    /* 粗糙度为零退化为镜面:粗糙平均必须逐角度等于解析菲涅耳 */
    double max_deg = 0.0;
    for (int k = 1; k <= 180; k++) {
        double ci = k / 180.0;
        double d1 = fabs(fresnel_avg_rough(ci, ETA_GLASS, ETA_WATER, 0.0, 0, NULL)
                         - fresnel_R_exact(ci, ETA_GLASS, ETA_WATER));
        double d2 = fabs(fresnel_avg_rough(ci, ETA_WATER, ETA_GLASS, 0.0, 0, NULL)
                         - fresnel_R_exact(ci, ETA_WATER, ETA_GLASS));
        if (d1 > max_deg) max_deg = d1;
        if (d2 > max_deg) max_deg = d2;
    }
    printf("report: rough_zero_max_dev=%.3e\n", max_deg);
    check("rough_zero_specular", max_deg < 1e-12);

    /* 粗糙扰动扫描:全角度总能量 R+T 不得超过 1,偏差与补偿量给出报告。
       误差上界:远离临界角为 O(a^2)(对称分布一阶项抵消),临界角附近因
       F ~= 1 - c*sqrt(theta_c-theta) 的平方根奇异性为 O(sqrt(a)),
       统一取 C*sqrt(a) 上界,C=2。 */
    double max_total_dev = 0.0, max_over = 0.0, max_comp = 0.0, max_err_ratio = 0.0;
    double roughs[] = { 0.05, 0.1, 0.2, 0.4 };
    for (int r = 0; r < 4; r++) {
        double a = roughs[r];
        double max_err = 0.0;
        for (int k = 1; k <= 180; k++) {
            double ci = k / 180.0;
            for (int d = 0; d < 2; d++) {
                double ni = d ? ETA_WATER : ETA_GLASS;
                double nt = d ? ETA_GLASS : ETA_WATER;
                double wr = 0.0;
                double R = fresnel_avg_rough(ci, ni, nt, a, 0, &wr);
                double T = fresnel_avg_rough(ci, ni, nt, a, 1, NULL);
                double total = R + T;
                if (fabs(total - 1.0) > max_total_dev) max_total_dev = fabs(total - 1.0);
                if (total - 1.0 > max_over) max_over = total - 1.0;
                if (1.0 - wr > max_comp) max_comp = 1.0 - wr;
                double err = fabs(R - fresnel_R_exact(ci, ni, nt));
                if (err > max_err) max_err = err;
            }
        }
        if (max_err / sqrt(a) > max_err_ratio) max_err_ratio = max_err / sqrt(a);
        printf("report: rough=%.2f max|R_rough-R_exact|=%.4f (bound %.4f)\n",
               a, max_err, 2.0 * sqrt(a));
    }
    printf("report: rough_total_max_dev=%.3e max_over_one=%.3e max_compensation=%.4f\n",
           max_total_dev, max_over, max_comp);
    check("rough_total_le_one", max_over <= 1e-6);
    check("rough_total_dev_small", max_total_dev < 1e-6);
    check("rough_err_bounded", max_err_ratio <= 2.0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
