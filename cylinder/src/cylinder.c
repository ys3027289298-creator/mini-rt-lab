
/* cylinder: finite open cylinder (side surface) ray intersection.
   二次方程求解 + 轴向 slab 裁剪；轴近平行/近垂直两条稳定分支；
   交叉积判别式防对消；long double 解析参照做误差上界核对。 */
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <float.h>

#ifndef SWEEP_N
#define SWEEP_N 200000
#endif

/* ---- 向量 ---- */
typedef struct { float x, y, z; } V3;

static V3 v3(float x, float y, float z) { V3 r = { x, y, z }; return r; }
static V3 vsub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static V3 vadd(V3 a, V3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static V3 vmul(V3 a, float s) { return v3(a.x * s, a.y * s, a.z * s); }
static float vdot(V3 a, V3 b) { return fmaf(a.x, b.x, fmaf(a.y, b.y, a.z * b.z)); }
static float vlen2(V3 a) { return vdot(a, a); }
static float vlen(V3 a) { return sqrtf(vlen2(a)); }
static V3 vcross(V3 a, V3 b) {
    return v3(a.y * b.z - a.z * b.y,
              a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}

/* ---- 几何对象与返回码 ---- */
typedef struct {
    V3 c;         /* 轴上基点：z=c.z=z0 处的圆心 */
    V3 axis;      /* 单位轴向 */
    float r;      /* 半径，必须 > 0 */
    float z0, z1; /* 轴向区间 [z0,z1]，z0 < z1 */
} Cylinder;

typedef struct {
    float t; /* 射线参数 */
    V3 p;    /* 命中点 */
    V3 n;    /* 朝外的法线（已取向使 dot(n,d) <= 0） */
} Hit;

enum {
    CYL_HIT = 0,
    CYL_MISS = 1,          /* 无命中 */
    CYL_ZERO_RADIUS = -1,  /* 半径为零或负：退化体，显式报告 */
    CYL_DEGENERATE = -2,   /* 轴为零向量或区间非法 */
    CYL_BAD_RAY = -3       /* 射线方向为零 */
};

/* 单舍入误差单位 u，与分支阈值（推导见 README） */
#define U_FLT      FLT_EPSILON        /* 2^-24 ≈ 1.1921e-7 */
#define EPS_PAR    1.0e-14            /* A/|d|^2 低于此走近平行分支（double 极限附近） */
#define DISC_GUARD (8.0 * 1.19209290e-16) /* double 判别式保护带（约 8*eps_dbl） */
#define PAR_GUARD  (8.0 * 1.19209290e-16) /* 平行分支 R^2 判定保护带 */
#define SIDE_EPS   (4.0 * 1.19209290e-16) /* 平行掠入点表面判定带 */

/* double 精度向量工具：输入输出仍是 float，仅在求交关键路径上提升精度，
   避免小半径/远起点/掠入射时投影与判别式的灾难性对消。 */
typedef struct { double x, y, z; } DV;
static DV dv(V3 a) { DV r = { a.x, a.y, a.z }; return r; }
static DV dsub(DV a, DV b) { DV r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static DV dadd(DV a, DV b) { DV r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static DV dmul(DV a, double s) { DV r = { a.x * s, a.y * s, a.z * s }; return r; }
static double ddot(DV a, DV b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static DV dcross(DV a, DV b) {
    DV r = { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
    return r;
}

/* 有限开口圆柱（仅侧面）与射线 o+t d 的最近命中。
   t 限制在 [tmin,tmax]，命中写入 *hit。返回码见 enum。 */
static int cyl_hit(const Cylinder* cy, V3 of, V3 df, float tmin, float tmax, Hit* hit) {
    if (!(cy->r > 0.0f)) return CYL_ZERO_RADIUS;
    if (vlen2(cy->axis) < 0.5f) return CYL_DEGENERATE;
    if (!(cy->z1 > cy->z0)) return CYL_DEGENERATE;
    if (vlen2(df) == 0.0f) return CYL_BAD_RAY;

    DV o = dv(of), d = dv(df), c = dv(cy->c), ax = dv(cy->axis);
    double r = cy->r, rr = r * r;

    /* 垂直于轴平面上的投影（double 累加）：w_perp = w - (w·axis) axis */
    double da = ddot(d, ax);
    DV dperp = dsub(d, dmul(ax, da));
    double A = ddot(dperp, dperp);          /* = |d|^2 - (d·axis)^2 */
    double d2 = ddot(d, d);

    DV oc = dsub(o, c);
    double oa = ddot(oc, ax);
    DV operp = dsub(oc, dmul(ax, oa));

    double B2 = ddot(dperp, operp);         /* 方程的 1/2 一次系数 */
    double C = ddot(operp, operp) - rr;

    /* 轴向 slab：z(t)=oa+t*da 落在 [z0,z1]，求 t 区间 [Tlo,Thi]（da 近退化稳定分支）。
       边界根（光线恰在端面处命中侧面）会因 z(t) 的舍入被推出区间而丢命中——
       这正是细长柱渲染缺口的来源。按 z(t) 的计算误差给 t 边界松弛，
       松弛只向内扩大 slab，不会凭空制造几何上不存在的命中。 */
    double Tlo = tmin, Thi = tmax;
    if (fabs(da) >= 1.0e-12 * sqrt(fmax(d2, 1.0e-30))) {
        double inv = 1.0 / da;
        double ta = cy->z0 * inv;
        double tb = cy->z1 * inv;
        if (ta > tb) { double tmp = ta; ta = tb; tb = tmp; }
        double span = fabs(tb - ta);
        double tslack = 4.0 * (double)U_FLT * (fabs(ta) + fabs(tb) + span);
        ta -= tslack;
        tb += tslack;
        if (ta > Tlo) Tlo = ta;
        if (tb < Thi) Thi = tb;
    } else {
        /* da≈0：轴向近退化，不做除法，仅判断起点轴向位置（含同样的相对松弛） */
        double zslack = 4.0 * (double)U_FLT
                        * fmax(fabs(oa), fmax(fabs((double)cy->z0), fabs((double)cy->z1)));
        if (oa < cy->z0 - zslack || oa > cy->z1 + zslack) return CYL_MISS;
    }
    if (Tlo > Thi) return CYL_MISS;

    if (A < (double)EPS_PAR * d2) {
        /* ---- 稳定分支①：射线与轴近平行（A→0，二次方程退化为线性）----
           不除 A、不除 sqrt(A)；直接用 operp 长度与 slab 内最近点判定，
           因此近退化时不会因对消或除零丢命中。 */
        double R2lo = C + 2.0 * B2 * Tlo + A * Tlo * Tlo; /* |operp(Tlo)|^2 - r^2 */
        if (R2lo > (double)PAR_GUARD * rr) return CYL_MISS; /* 整条 slab 都在柱外 */
        if (R2lo >= -(double)SIDE_EPS * rr) {
            /* 掠入：在 Tlo 处贴着侧面，记一次命中（精确平行光线唯一的侧面命中机会） */
            DV p = dadd(o, dmul(d, Tlo));
            DV radial = dsub(p, dadd(c, dmul(ax, ddot(dsub(p, c), ax))));
            double rl = sqrt(ddot(radial, radial));
            DV n = (rl > 0.0) ? dmul(radial, 1.0 / rl) : (DV){ 0, 0, 0 };
            if (ddot(n, d) > 0.0) n = dmul(n, -1.0);
            hit->t = (float)Tlo;
            hit->p = v3((float)p.x, (float)p.y, (float)p.z);
            hit->n = v3((float)n.x, (float)n.y, (float)n.z);
            return CYL_HIT;
        }
        /* 柱内平行：侧面在 slab 外才相交，开口圆柱侧面在 slab 内不命中 */
        return CYL_MISS;
    }

    /* ---- 一般分支：交叉积判别式（double 精度下的稳定符号测试）----
       B2^2 - A C 与  Ar^2 - |d_perp x o_perp|^2  恒等，两项均非负。
       double 累加把对消余量从 ~1e-7 提升到 ~1e-16；再叠加保护带吸收
       输入 float 量化与剩余舍入，带内负值视为切线钳到零，不丢命中。 */
    DV cr = dcross(dperp, operp);
    double Ar2 = A * rr;
    double cr2 = ddot(cr, cr);
    double disc2 = Ar2 - cr2;
    double disc_band = (double)DISC_GUARD * (Ar2 + cr2);
    if (disc2 < -disc_band) return CYL_MISS; /* 误差带之外仍为负：真无交 */
    if (disc2 < 0.0) disc2 = 0.0;            /* 带内当切线，不丢命中 */

    /* 稳定二次求根（q 形式，避免 -b±sqrt 的对消丢根），全程 double */
    double s = sqrt(disc2);
    double q = (B2 >= 0.0) ? -(B2 + s) : -(B2 - s);
    double ta, tb;
    if (q == 0.0) { ta = 0.0; tb = 0.0; }
    else {
        double a = q / A;
        double b = C / q;
        ta = (a <= b) ? a : b;
        tb = (a <= b) ? b : a;
    }

    /* slab 裁剪取最近侧面命中 */
    double tc = ta;
    if (tc < Tlo) tc = tb;
    if (tc < Tlo || tc > Thi) return CYL_MISS;

    DV p = dadd(o, dmul(d, tc));
    DV radial = dsub(p, dadd(c, dmul(ax, ddot(dsub(p, c), ax))));
    double rl = sqrt(ddot(radial, radial));
    DV n = (rl > 0.0) ? dmul(radial, 1.0 / rl) : (DV){ 0, 0, 0 };
    if (ddot(n, d) > 0.0) n = dmul(n, -1.0);
    hit->t = (float)tc;
    hit->p = v3((float)p.x, (float)p.y, (float)p.z);
    hit->n = v3((float)n.x, (float)n.y, (float)n.z);
    return CYL_HIT;
}

/* =====================================================================
   long double 解析参照：同样的几何，用 80 位长双精度独立重算，
   作为"解析解"核对 float 实现的误差上界。
   ===================================================================== */
typedef struct { long double x, y, z; } LV;

static LV lv(V3 a) { LV r = { a.x, a.y, a.z }; return r; }
static LV lsub(LV a, LV b) { LV r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static LV ladd(LV a, LV b) { LV r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static LV lmul(LV a, long double s) { LV r = { a.x * s, a.y * s, a.z * s }; return r; }
static long double ldot(LV a, LV b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static LV lcross(LV a, LV b) {
    LV r = { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
    return r;
}

/* 参照命中判定。返回 1 表示在 [tmin,tmax] 内有侧面命中，*tref 给出最近参数。
   *tangent 置 1 表示判别式落在 double 保护带内：这类情形产品端报 hit 或 miss
   都在数值误差允许范围内，不算漏/误命中。 */
static int ref_hit(const Cylinder* cy, V3 of, V3 df,
                   long double tmin, long double tmax, long double* tref, int* tangent) {
    LV o = lv(of), d = lv(df), c = lv(cy->c), ax = lv(cy->axis);
    long double r = cy->r, r2 = r * r;
    *tangent = 0;

    long double da = ldot(d, ax);
    LV dperp = lsub(d, lmul(ax, da));
    long double A = ldot(dperp, dperp);
    LV oc = lsub(o, c);
    long double oa = ldot(oc, ax);
    LV operp = lsub(oc, lmul(ax, oa));
    long double B2 = ldot(dperp, operp);
    long double C = ldot(operp, operp) - r2;

    long double Tlo = tmin, Thi = tmax;
    if (da != 0.0L) {
        long double ta = cy->z0 / da, tb = cy->z1 / da;
        if (ta > tb) { long double tmp = ta; ta = tb; tb = tmp; }
        /* 与产品端同一公式、按 float 舍入单位给的 slab 松弛，使边界根两侧判定一致 */
        long double span = fabsl(tb - ta);
        long double tslack = 4.0L * (long double)U_FLT * (fabsl(ta) + fabsl(tb) + span);
        ta -= tslack;
        tb += tslack;
        if (ta > Tlo) Tlo = ta;
        if (tb < Thi) Thi = tb;
    } else {
        long double zslack = 4.0L * (long double)U_FLT
                             * fmaxl(fabsl(oa), fmaxl(fabsl((long double)cy->z0),
                                                       fabsl((long double)cy->z1)));
        if (oa < cy->z0 - zslack || oa > cy->z1 + zslack) return 0;
    }
    if (Tlo > Thi) return 0;

    long double t1, t2;
    if (A == 0.0L) {
        /* 精确平行：仅当恰好贴侧面时在 slab 边界产生命中（测度零，参照给 miss） */
        return 0;
    }
    long double disc = B2 * B2 - A * C;
    /* 与产品端同一保护带公式：Ar^2 + |d_perp x o_perp|^2，
       带内的判别式符号不可信，产品端报 hit 或 miss 都不算硬不一致。 */
    LV lcr = lcross(dperp, operp);
    long double Ar2 = A * r2;
    long double cr2 = ldot(lcr, lcr);
    long double band = (long double)DISC_GUARD * (Ar2 + cr2);
    if (disc < -band) return 0;
    if (disc < 0.0L) { disc = 0.0L; *tangent = 1; }
    long double s = sqrtl(disc);
    t1 = (-B2 - s) / A;
    t2 = (-B2 + s) / A;
    long double tc = t1;
    if (tc < Tlo) tc = t2;
    if (tc < Tlo || tc > Thi) return 0;
    *tref = tc;
    return 1;
}

/* 后验残差上界（可在产品端直接计算的向后误差）。
   求交在 double 关键路径上完成、根回存 float，故残差界由两部分组成：
     gamma16_d * 多项式规模量：double 投影/点积/q 求根的舍入链；
     u_flt * 多项式规模量：根 t 回存 float 的量化（主导项）；
   以及 DISC_GUARD*r^2*|A| 的判别式保护带。 */
#define U_DBL 2.2204460492503130e-16
static double residual_bound(double A, double B2, double C, double rr, double t) {
    const double g16 = 16.0 * U_DBL / (1.0 - 16.0 * U_DBL);
    double at = fabs(t);
    double ops = fabs(A) * at * at + 2.0 * fabs(B2) * at + fabs(C);
    return (g16 + (double)U_FLT) * ops + (double)DISC_GUARD * rr * fabs(A);
}

/* =====================================================================
   自检
   ===================================================================== */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* 可复现的 LCG，种子可经 -DRNG_SEED=... 覆盖以做多轮差分 */
#ifndef RNG_SEED
#define RNG_SEED 0x123456789abcdef0ULL
#endif
static uint64_t rng_state = RNG_SEED;
static float rng_f(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (float)((rng_state >> 11) & 0xFFFFFF) / (float)0x1000000;
}
static float rng_range(float lo, float hi) { return lo + rng_f() * (hi - lo); }

int main(void) {
    /* 标准柱：轴 = z，基点 (0,0,0)，r=1，z∈[0,4] */
    Cylinder cy = { v3(0, 0, 0), v3(0, 0, 1), 1.0f, 0.0f, 4.0f };
    Hit h;

    /* 1. 垂直入射打柱身：解析 t=4，命中 (1,0,1)，法线 +x */
    int rc = cyl_hit(&cy, v3(5, 0, 1), v3(-1, 0, 0), 0.0f, 100.0f, &h);
    check("radial_hit", rc == CYL_HIT && fabsf(h.t - 4.0f) < 1e-6f);
    check("radial_point", fabsf(h.p.x - 1.0f) < 1e-6f && fabsf(h.p.y) < 1e-6f);
    check("radial_normal", fabsf(h.n.x - 1.0f) < 1e-6f && fabsf(vdot(h.n, v3(-1, 0, 0)) + 1.0f) < 1e-6f);

    /* 2. 切线：判别式恰为零，不许丢命中 */
    rc = cyl_hit(&cy, v3(0, 1, 1), v3(-1, 0, 0), 0.0f, 100.0f, &h);
    check("tangent_kept", rc == CYL_HIT && fabsf(h.t) < 1e-5f &&
                           fabsf(h.n.y - 1.0f) < 1e-5f);

    /* 3. 内部出发：取穿出点 t2（裁剪后不能误判 miss），解析 t=1。
       出口几何法线 +x，按着色约定翻转成朝向光线，故报告 -x。 */
    rc = cyl_hit(&cy, v3(0, 0, 1), v3(1, 0, 0), 1e-4f, 100.0f, &h);
    check("inside_exit", rc == CYL_HIT && fabsf(h.t - 1.0f) < 1e-6f &&
                         fabsf(h.n.x + 1.0f) < 1e-6f);

    /* 4. 命中点轴向超出 slab：裁剪后 miss */
    rc = cyl_hit(&cy, v3(5, 0, 9), v3(-1, 0, 0), 0.0f, 100.0f, &h);
    check("clip_miss", rc == CYL_MISS);

    /* 5. tmax 裁剪：真实命中 t=4，tmax=3 时被剔除 */
    rc = cyl_hit(&cy, v3(5, 0, 1), v3(-1, 0, 0), 0.0f, 3.0f, &h);
    check("tmax_cull", rc == CYL_MISS);

    /* 6. 反向（最近根为负） */
    rc = cyl_hit(&cy, v3(5, 0, 1), v3(1, 0, 0), 0.0f, 100.0f, &h);
    check("behind_miss", rc == CYL_MISS);

    /* 7. 精确平行于轴：柱内，开口圆柱侧面不命中 */
    rc = cyl_hit(&cy, v3(0.2f, 0, 0.5f), v3(0, 0, 1), 0.0f, 100.0f, &h);
    check("parallel_inside_miss", rc == CYL_MISS);

    /* 8. 精确平行于轴：柱外 */
    rc = cyl_hit(&cy, v3(2, 0, 0.5f), v3(0, 0, 1), 0.0f, 100.0f, &h);
    check("parallel_outside_miss", rc == CYL_MISS);

    /* 9. 近平行（轴近退化）掠入：距轴 r 的平行光线，保护带内不许丢命中 */
    rc = cyl_hit(&cy, v3(1.0f, 0, 0.5f), v3(0, 0, 1), 0.0f, 100.0f, &h);
    check("grazing_parallel_kept", rc == CYL_HIT && fabsf(h.t) < 1e-5f);

    /* 10. 近平行但略倾斜（parfrac≈1e-6，走二次方程分支），slab 中部进入柱身。
       解析：1.001-1e-3 t=1 → t=1，z(t)=2+t=3 仍在 [0,4] 内。缺口回归用例。 */
    rc = cyl_hit(&cy, v3(1.001f, 0, 2.0f), v3(-1e-3f, 0, 1), 0.0f, 100.0f, &h);
    check("near_parallel_entry", rc == CYL_HIT);
    if (rc == CYL_HIT) {
        check("near_parallel_t", fabsf(h.t - 1.0f) < 1e-3f);
    } else {
        check("near_parallel_t", 0);
    }

    /* 11. 近垂直于轴（da 近退化的 slab 分支）：光线沿轴分量极微但仍在 slab 内 */
    rc = cyl_hit(&cy, v3(5, 0, 2), v3(-1, 0, 1e-14f), 0.0f, 100.0f, &h);
    check("near_perp_hit", rc == CYL_HIT && fabsf(h.t - 4.0f) < 1e-4f);

    /* 12. 斜轴圆柱：轴 (0,1,0)，基点 (0,0,0)，z 坐标变成"垂直轴"分量，不参与 */
    Cylinder cyx = { v3(0, 0, 0), v3(0, 1, 0), 1.0f, 0.0f, 4.0f };
    rc = cyl_hit(&cyx, v3(5, 0, 0), v3(-1, 0, 0), 0.0f, 100.0f, &h);
    check("tilted_axis_hit", rc == CYL_HIT && fabsf(h.t - 4.0f) < 1e-6f);

    /* 13. 退化输入显式报告 */
    Cylinder bad_r = cy; bad_r.r = 0.0f;
    check("zero_radius_reported", cyl_hit(&bad_r, v3(5, 0, 1), v3(-1, 0, 0), 0, 100, &h) == CYL_ZERO_RADIUS);
    Cylinder neg_r = cy; neg_r.r = -1.0f;
    check("neg_radius_reported", cyl_hit(&neg_r, v3(5, 0, 1), v3(-1, 0, 0), 0, 100, &h) == CYL_ZERO_RADIUS);
    Cylinder bad_ax = cy; bad_ax.axis = v3(0, 0, 0);
    check("zero_axis_reported", cyl_hit(&bad_ax, v3(5, 0, 1), v3(-1, 0, 0), 0, 100, &h) == CYL_DEGENERATE);
    Cylinder bad_iv = cy; bad_iv.z0 = 4.0f; bad_iv.z1 = 4.0f;
    check("bad_interval_reported", cyl_hit(&bad_iv, v3(5, 0, 1), v3(-1, 0, 0), 0, 100, &h) == CYL_DEGENERATE);
    check("zero_dir_reported", cyl_hit(&cy, v3(5, 0, 1), v3(0, 0, 0), 0, 100, &h) == CYL_BAD_RAY);

    /* ---- 随机差分测试：float 实现 vs long double 解析参照 ---- */
    enum { N_SWEEP = SWEEP_N };
    int mismatch = 0;                 /* 保护带之外的命中/未命中硬不一致（真丢命中） */
    int guard_accepts = 0;            /* 参照在切线带内，产品按保护带报 hit（合法） */
    int n_hits = 0, n_near = 0;
    float max_fwd_rel = 0.0f;         /* 最大正向相对误差 |t-tref|/max(1,|tref|) */
    float max_fwd_ratio = 0.0f;       /* 正向误差 / 推导上界，要求 <=1 */
    float max_res_ratio = 0.0f;       /* 残差 / 后验残差上界，要求 <=1 */
    float max_nang = 0.0f;            /* 法线角度误差（弧度） */

    for (int i = 0; i < N_SWEEP; i++) {
        /* 随机柱：轴方向在球面上，半径与长度跨尺度 */
        float az = rng_range(0.0f, 6.2831853f);
        float el = asinf(rng_range(-1.0f, 1.0f));
        V3 axis = v3(cosf(el) * cosf(az), cosf(el) * sinf(az), sinf(el));
        float r = rng_range(1e-3f, 10.0f);
        float half = rng_range(0.1f, 50.0f);
        Cylinder rnd = { v3(rng_range(-5, 5), rng_range(-5, 5), rng_range(-5, 5)),
                         axis, r, -half, half };

        V3 o = v3(rng_range(-30, 30), rng_range(-30, 30), rng_range(-30, 30));
        V3 d;
        int bucket = i % 5;
        if (bucket == 0) {
            /* 近平行：d = axis*大 + 极小垂直分量 */
            V3 t1 = v3(1, 0, 0);
            if (fabsf(axis.x) > 0.9f) t1 = v3(0, 1, 0);
            V3 perp = vsub(t1, vmul(axis, vdot(t1, axis)));
            float pl = vlen(perp);
            perp = vmul(perp, 1.0f / pl);
            d = vadd(vmul(axis, rng_range(0.5f, 2.0f)),
                     vmul(perp, rng_range(-1e-5f, 1e-5f)));
        } else if (bucket == 1) {
            /* 近垂直于轴：da 极小 */
            V3 t1 = v3(1, 0, 0);
            if (fabsf(axis.x) > 0.9f) t1 = v3(0, 1, 0);
            V3 perp = vsub(t1, vmul(axis, vdot(t1, axis)));
            perp = vmul(perp, 1.0f / vlen(perp));
            d = vadd(perp, vmul(axis, rng_range(-1e-6f, 1e-6f)));
        } else if (bucket == 2) {
            /* 切线附近：先取柱身一点，方向使其擦过 */
            float phi = rng_range(0, 6.2831853f);
            V3 radial = vsub(v3(cosf(phi), sinf(phi), 0),
                             vmul(axis, vdot(v3(cosf(phi), sinf(phi), 0), axis)));
            radial = vmul(radial, 1.0f / vlen(radial));
            V3 p = vadd(rnd.c, vadd(vmul(radial, r), vmul(axis, rng_range(-half, half))));
            V3 tan = vcross(axis, radial);
            o = vadd(p, vmul(radial, rng_range(-0.05f, 0.05f) * r));
            d = vadd(tan, vmul(radial, rng_range(-0.02f, 0.02f)));
        } else if (bucket == 3) {
            /* 一般随机方向 */
            float da2 = rng_range(0, 6.2831853f), de2 = asinf(rng_range(-1, 1));
            d = v3(cosf(de2) * cosf(da2), cosf(de2) * sinf(da2), sinf(de2));
        } else {
            /* 远距离小目标（细长柱）：半径 1e-3 量级、起点很远 */
            d = v3(rng_range(-1, 1), rng_range(-1, 1), rng_range(-1, 1));
        }

        long double tref = 0.0L;
        int tangent = 0;
        int want = ref_hit(&rnd, o, d, 1e-4L, 1000.0L, &tref, &tangent);
        int got_rc = cyl_hit(&rnd, o, d, 1e-4f, 1000.0f, &h);
        int got = (got_rc == CYL_HIT);

        /* 近退化比例 parfrac=A/|d|^2（double 下计算），用于统计与近退化分桶 */
        double da_d = ddot(dv(d), dv(rnd.axis));
        double d2_d = ddot(dv(d), dv(d));
        double parfrac = fmax(1.0 - da_d * da_d / d2_d, 1.19209290e-16);
        n_near += (parfrac < 1.0e-6);

        if (want != got) {
            if (tangent && got) { guard_accepts++; continue; } /* 保护带救回切线，允许 */
            mismatch++;
            if (mismatch <= 10)
                printf("  mismatch b=%d want=%d got=%d pf=%g tangent=%d tref=%Lf\n",
                       bucket, want, got, (float)parfrac, tangent, tref);
            continue;
        }
        if (!got) continue;
        n_hits++;

        float fwd = fabsf(h.t - (float)tref) / fmaxf(1.0f, fabsf((float)tref));
        if (fwd > max_fwd_rel) max_fwd_rel = fwd;

        /* 与实现同路径的 double 系数，用于误差核对 */
        DV dd = dv(d), cdd = dv(rnd.c), axd = dv(rnd.axis);
        DV dperpd = dsub(dd, dmul(axd, ddot(dd, axd)));
        double A = ddot(dperpd, dperpd);
        DV ocd = dsub(dv(o), cdd);
        DV operpd = dsub(ocd, dmul(axd, ddot(ocd, axd)));
        double B2 = ddot(dperpd, operpd);
        double rr = (double)r * r;
        double C = ddot(operpd, operpd) - rr;
        double th = h.t;

        /* 残差（double 求值）与后验上界 */
        double resid = fabs(A * th * th + 2.0 * B2 * th + C);
        double rb = residual_bound(A, B2, C, rr, th);
        double ratio = (rb > 0.0) ? resid / rb : 0.0;
        if (ratio > max_res_ratio) max_res_ratio = (float)ratio;

        /* 正向误差推导上界。double 求交后系数对消已消除，主导项是：
           (a) 分离根一阶扰动 |Δt| <= rb / |P'(t)|；
           (b) 近重根平方根分裂 sqrt(rb/A)（切线，O(sqrt(eps_dbl))）；
           (c) 保护带钳位位移 sqrt(DISC_GUARD)*r/sqrt(A)；
           (d) 近零根符号翻转后两根间距 2|B2|/A；
           以及结果回存 float 的量化地板 ~u*max(1,|t|)。 */
        double deriv = fabs(2.0 * A * th + 2.0 * B2);
        double back = rb;
        double disc_dd = A * rr - ddot(dcross(dperpd, operpd), dcross(dperpd, operpd));
        double band_d = (double)DISC_GUARD * (A * rr + ddot(dcross(dperpd, operpd),
                                                              dcross(dperpd, operpd)));
        double clamp_shift = (disc_dd < 0.0 && disc_dd > -band_d)
                             ? sqrt((double)DISC_GUARD) * r / sqrt(A) : 0.0;
        double root_floor = 2.0 * fabs(B2) / A;
        double tol;
        if (deriv > 0.0) {
            tol = fmax(back / deriv, sqrt(back / A));
            tol = fmax(tol, clamp_shift);
            tol = fmax(tol, root_floor);
            tol = fmax(tol, (double)U_FLT * fmax(1.0, fabs(th)));
        } else {
            tol = fmax(sqrt(back / A), fmax(clamp_shift, root_floor));
        }
        float fratio = (float)(fabs(th - (double)tref) / tol);
        if (fratio > max_fwd_ratio) max_fwd_ratio = fratio;

        /* 法线与解析径向方向比较 */
        LV lp = ladd(lv(o), lmul(lv(d), tref));
        LV lrad = lsub(lp, ladd(lv(rnd.c), lmul(lv(rnd.axis),
                     ldot(lsub(lp, lv(rnd.c)), lv(rnd.axis)))));
        long double ll = sqrtl(ldot(lrad, lrad));
        LV nref = lmul(lrad, 1.0L / ll);
        float nd = fmaxf(-1.0f, fminf(1.0f,
                    h.n.x * (float)nref.x + h.n.y * (float)nref.y + h.n.z * (float)nref.z));
        float ang = acosf(fabsf(nd));
        if (ang > max_nang) max_nang = ang;
    }

    printf("sweep_hits=%d sweep_near_parallel=%d guard_accepts=%d\n", n_hits, n_near, guard_accepts);
    printf("sweep_mismatch=%d\n", mismatch);
    printf("max_forward_rel_err=%g\n", max_fwd_rel);
    printf("max_forward_over_bound=%g (要求 <=1)\n", max_fwd_ratio);
    printf("max_residual_over_bound=%g (要求 <=1)\n", max_res_ratio);
    printf("max_normal_angle_err_rad=%g\n", max_nang);

    check("no_missed_hits", mismatch == 0);
    check("forward_inside_bound", max_fwd_ratio <= 1.0f);
    check("residual_inside_bound", max_res_ratio <= 1.0f);
    /* 法线：常规命中远优于 1e-3；重根保护带内的最坏情形约 sqrt(8u) 量级，放宽到 5e-2 */
    check("normal_accuracy", max_nang < 5.0e-2f);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
