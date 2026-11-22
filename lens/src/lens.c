
/* lens: Brown-Conrady 径向+切向畸变正模型、牛顿反解、收敛阶与单调区间、
   解析反解(Cardano)比对与误差上界、系数可辨识性、未收敛像素计数。 */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846

typedef struct { double k1, k2, k3, p1, p2; } Dist;

/* ---- 正模型 ----
   r2 = x^2+y^2, rad = 1 + k1 r2 + k2 r2^2 + k3 r2^3
   xd = x*rad + 2 p1 x y + p2 (r2 + 2 x^2)
   yd = y*rad + p1 (r2 + 2 y^2) + 2 p2 x y                       */
static void distort(const Dist* c, double x, double y, double* xd, double* yd) {
    double r2 = x * x + y * y;
    double rad = 1.0 + c->k1 * r2 + c->k2 * r2 * r2 + c->k3 * r2 * r2 * r2;
    *xd = x * rad + 2.0 * c->p1 * x * y + c->p2 * (r2 + 2.0 * x * x);
    *yd = y * rad + c->p1 * (r2 + 2.0 * y * y) + 2.0 * c->p2 * x * y;
}

static void distort_jac(const Dist* c, double x, double y, double J[2][2]) {
    double r2 = x * x + y * y;
    double rad = 1.0 + c->k1 * r2 + c->k2 * r2 * r2 + c->k3 * r2 * r2 * r2;
    double dr = c->k1 + 2.0 * c->k2 * r2 + 3.0 * c->k3 * r2 * r2; /* d(rad)/d(r2) */
    J[0][0] = rad + 2.0 * dr * x * x + 2.0 * c->p1 * y + 6.0 * c->p2 * x;
    J[0][1] = 2.0 * dr * x * y + 2.0 * c->p1 * x + 2.0 * c->p2 * y;
    J[1][0] = J[0][1];
    J[1][1] = rad + 2.0 * dr * y * y + 6.0 * c->p1 * y + 2.0 * c->p2 * x;
}

/* ---- 实系数三次方程实根: a x^3 + b x^2 + c x + d = 0 ---- */
static int solve_cubic_real(double a, double b, double c, double d, double roots[3]) {
    if (fabs(a) < 1e-30) {
        if (fabs(b) < 1e-30) {
            if (fabs(c) < 1e-30) return 0;
            roots[0] = -d / c;
            return 1;
        }
        double disc = c * c - 4.0 * b * d;
        if (disc < 0.0) return 0;
        double sq = sqrt(disc);
        roots[0] = (-c - sq) / (2.0 * b);
        roots[1] = (-c + sq) / (2.0 * b);
        return 2;
    }
    double A = b / a, B = c / a, C = d / a;
    double p = B - A * A / 3.0;
    double q = 2.0 * A * A * A / 27.0 - A * B / 3.0 + C;
    double off = -A / 3.0;
    double disc = q * q / 4.0 + p * p * p / 27.0;
    if (disc > 0.0) {
        double sq = sqrt(disc);
        roots[0] = cbrt(-q / 2.0 + sq) + cbrt(-q / 2.0 - sq) + off;
        return 1;
    }
    if (fabs(p) < 1e-30) { roots[0] = off; return 1; }
    double m = 2.0 * sqrt(-p / 3.0);
    double arg = 3.0 * q / (2.0 * p) * sqrt(-3.0 / p);
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;
    double theta = acos(arg) / 3.0;
    for (int k = 0; k < 3; k++)
        roots[k] = m * cos(theta - 2.0 * PI * k / 3.0) + off;
    return 3;
}

/* ---- 单调区间 ----
   径向映射 g(r) = r (1 + k1 r^2 + k2 r^4 + k3 r^6) 可逆当且仅当
   g'(r) = 1 + 3 k1 s + 5 k2 s^2 + 7 k3 s^3 > 0 (s = r^2)。
   返回使 g' > 0 成立的最大 r；无正根则返回 HUGE_VAL。 */
static double monotonic_radius(const Dist* c) {
    double roots[3];
    int n = solve_cubic_real(7.0 * c->k3, 5.0 * c->k2, 3.0 * c->k1, 1.0, roots);
    double best = HUGE_VAL;
    for (int i = 0; i < n; i++)
        if (roots[i] > 1e-12 && roots[i] < best) best = roots[i];
    return best == HUGE_VAL ? HUGE_VAL : sqrt(best);
}

/* ---- 牛顿反解 ----
   初值取恒等 (x,y)=(xd,yd)。truth 非空时记录误差序列并估计收敛阶
   q: |e_{n+1}| ~ C |e_n|^q。返回 1=收敛。 */
static int undistort(const Dist* c, double xd, double yd, double* x, double* y,
                     int maxit, double tol, const double* truth,
                     double* order, int* iters) {
    double ux = xd, uy = yd;
    double hist[64];
    int nh = 0, conv = 0, it;
    for (it = 0; it < maxit; it++) {
        double fx, fy;
        distort(c, ux, uy, &fx, &fy);
        double rx = fx - xd, ry = fy - yd;
        if (truth && nh < 64) hist[nh++] = hypot(ux - truth[0], uy - truth[1]);
        if (hypot(rx, ry) < tol) { conv = 1; break; }
        double J[2][2];
        distort_jac(c, ux, uy, J);
        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (fabs(det) < 1e-14) break;
        double dx = (J[1][1] * rx - J[0][1] * ry) / det;
        double dy = (J[0][0] * ry - J[1][0] * rx) / det;
        ux -= dx;
        uy -= dy;
        if (!isfinite(ux) || !isfinite(uy) || fabs(ux) > 1e3 || fabs(uy) > 1e3) break;
    }
    if (truth && conv && nh < 64) hist[nh++] = hypot(ux - truth[0], uy - truth[1]);
    if (order) {
        /* 对 log|e_{n+1}| ~ q log|e_n| 做最小二乘, 只用远离舍入误差的点对 */
        double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
        int m = 0;
        for (int i = 0; i + 1 < nh; i++) {
            if (hist[i] > 1e-9 && hist[i + 1] > 0.0 && hist[i + 1] < hist[i]) {
                double lx = log(hist[i]), ly = log(hist[i + 1]);
                sx += lx; sy += ly; sxx += lx * lx; sxy += lx * ly;
                m++;
            }
        }
        *order = 0.0;
        if (m >= 2) {
            double den = m * sxx - sx * sx;
            if (fabs(den) > 1e-30) *order = (m * sxy - sx * sy) / den;
        }
    }
    *x = ux;
    *y = uy;
    if (iters) *iters = it;
    return conv;
}

/* ---- 解析反解(仅径向 k1): rd = r + k1 r^3, Cardano, 取单调支上的根 ---- */
static double undistort_radial_k1(double k1, double rd, double rmax) {
    if (fabs(k1) < 1e-30) return rd;
    double roots[3];
    int n = solve_cubic_real(k1, 0.0, 1.0, -rd, roots);
    double best = 0.0;
    int found = 0;
    for (int i = 0; i < n; i++) {
        double r = roots[i];
        if (((rd >= 0.0) == (r >= 0.0)) && fabs(r) <= rmax + 1e-9) {
            if (!found || fabs(r) < fabs(best)) { best = r; found = 1; }
        }
    }
    return found ? best : NAN;
}

/* ---- 可辨识性 ----
   畸变关于系数 (k1,k2,k3,p1,p2) 是线性的，观测雅可比与系数取值无关，
   只依赖采样构型。对 A = J^T J 做 Jacobi 特征分解，
   近零特征值对应的特征向量支撑集即不可辨识系数组。 */
static const char* PNAMES[5] = { "k1", "k2", "k3", "p1", "p2" };

static void jacobi5(const double Ain[5][5], double eval[5], double V[5][5]) {
    double A[5][5];
    memcpy(A, Ain, sizeof A);
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) V[i][j] = (i == j) ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 200; sweep++) {
        double off = 0.0;
        for (int p = 0; p < 5; p++)
            for (int q = p + 1; q < 5; q++) off += A[p][q] * A[p][q];
        if (off < 1e-24) break;
        for (int p = 0; p < 5; p++) {
            for (int q = p + 1; q < 5; q++) {
                if (fabs(A[p][q]) < 1e-300) continue;
                double tau = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                double t = (tau >= 0.0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
                double co = 1.0 / sqrt(1.0 + t * t), si = t * co;
                for (int k = 0; k < 5; k++) {
                    double akp = A[k][p], akq = A[k][q];
                    A[k][p] = co * akp - si * akq;
                    A[k][q] = si * akp + co * akq;
                }
                for (int k = 0; k < 5; k++) {
                    double apk = A[p][k], aqk = A[q][k];
                    A[p][k] = co * apk - si * aqk;
                    A[q][k] = si * apk + co * aqk;
                }
                for (int k = 0; k < 5; k++) {
                    double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = co * vkp - si * vkq;
                    V[k][q] = si * vkp + co * vkq;
                }
            }
        }
    }
    for (int i = 0; i < 5; i++) eval[i] = A[i][i];
}

static int ident_report(const char* name, const double* sx, const double* sy, int n) {
    double A[5][5] = { { 0 } };
    for (int i = 0; i < n; i++) {
        double x = sx[i], y = sy[i];
        double r2 = x * x + y * y, r4 = r2 * r2, r6 = r4 * r2;
        double jx[5] = { x * r2, x * r4, x * r6, 2.0 * x * y, r2 + 2.0 * x * x };
        double jy[5] = { y * r2, y * r4, y * r6, r2 + 2.0 * y * y, 2.0 * x * y };
        for (int a = 0; a < 5; a++)
            for (int b = 0; b < 5; b++)
                A[a][b] += jx[a] * jx[b] + jy[a] * jy[b];
    }
    double eval[5], V[5][5];
    jacobi5(A, eval, V);
    double lmax = 0.0;
    for (int i = 0; i < 5; i++) if (fabs(eval[i]) > lmax) lmax = fabs(eval[i]);
    printf("eig_%s=", name);
    for (int i = 0; i < 5; i++) printf("%s%.3e", i ? "," : "", eval[i]);
    printf("\n");
    int nulls = 0;
    int usupport[5] = { 0 };
    for (int i = 0; i < 5; i++) {
        if (fabs(eval[i]) <= 1e-10 * lmax) {
            nulls++;
            for (int j = 0; j < 5; j++)
                if (fabs(V[j][i]) > 0.2) usupport[j] = 1;
        }
    }
    if (nulls > 0) {
        printf("unidentifiable_%s={", name);
        int first = 1;
        for (int j = 0; j < 5; j++) {
            if (usupport[j]) {
                printf("%s%s", first ? "" : ",", PNAMES[j]);
                first = 0;
            }
        }
        printf("}\n");
    }
    printf("nullity_%s=%d\n", name, nulls);
    return nulls;
}

/* ---- 像素扫描: 统计未收敛(或落到单调区间外的折叠支)像素数 ---- */
static long sweep(const Dist* c, int W, int H, double f, double rgate) {
    long bad = 0;
    double cx = W * 0.5, cy = H * 0.5;
    for (int py = 0; py < H; py++) {
        for (int px = 0; px < W; px++) {
            double xd = (px + 0.5 - cx) / f, yd = (py + 0.5 - cy) / f;
            double x, y;
            int ok = undistort(c, xd, yd, &x, &y, 50, 1e-12, NULL, NULL, NULL);
            if (!ok || hypot(x, y) > rgate + 1e-6) bad++;
        }
    }
    return bad;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* 广角(非单调)与常规(单调)两组系数 */
    Dist wide = { -0.5, 0.1, 0.0, 1e-3, -1.5e-3 };
    Dist mild = { -0.1, 0.02, 0.0, 5e-4, -5e-4 };

    /* 1. 往返一致性 */
    double xd, yd, ux, uy;
    distort(&wide, 0.3, -0.2, &xd, &yd);
    int ok = undistort(&wide, xd, yd, &ux, &uy, 50, 1e-13, NULL, NULL, NULL);
    check("roundtrip", ok && fabs(ux - 0.3) < 1e-9 && fabs(uy + 0.2) < 1e-9);

    /* 2. 单调区间 */
    double rs_wide = monotonic_radius(&wide);
    double rs_mild = monotonic_radius(&mild);
    printf("monotonic_interval_wide=[0,%.6f]\n", rs_wide);
    printf("monotonic_interval_mild=[0,%s]\n", isinf(rs_mild) ? "inf" : "?");
    double s9 = 0.9 * rs_wide;
    double gp = 1.0 + 3.0 * wide.k1 * s9 * s9 + 5.0 * wide.k2 * s9 * s9 * s9 * s9;
    check("monotone_inside", gp > 0.0);
    check("monotone_boundary", fabs(1.0 + 3.0 * wide.k1 * rs_wide * rs_wide +
          5.0 * wide.k2 * rs_wide * rs_wide * rs_wide * rs_wide) < 1e-6);

    /* 3. 收敛阶(单调区间内, 真解已知) */
    double qsum = 0.0;
    int qn = 0;
    for (int i = 0; i < 3; i++) {
        double r = 0.2 + 0.15 * i, ang = 0.7;
        double tx = r * cos(ang), ty = r * sin(ang);
        double tdx, tdy, rx, ry, order;
        distort(&wide, tx, ty, &tdx, &tdy);
        double truth[2] = { tx, ty };
        if (undistort(&wide, tdx, tdy, &rx, &ry, 50, 1e-13, truth, &order, NULL) && order > 0.0) {
            printf("order_r%.2f=%.3f\n", r, order);
            qsum += order;
            qn++;
        }
    }
    double qavg = qn ? qsum / qn : 0.0;
    printf("convergence_order=%.3f\n", qavg);
    check("order_quadratic", qavg >= 1.95);

    /* 4. 与解析反解比对及误差上界(仅径向 k1 = -0.3) */
    Dist rad = { -0.3, 0.0, 0.0, 0.0, 0.0 };
    double rmax = monotonic_radius(&rad);
    double rtest = 0.9 * rmax;
    double gmin = 1.0 + 3.0 * rad.k1 * rtest * rtest; /* g' 在 [0,rtest] 上的下界 */
    double max_err = 0.0, max_bound = 0.0;
    int bound_ok = 1;
    for (int i = 1; i <= 20; i++) {
        double r = rtest * i / 20.0;
        double rd = r * (1.0 + rad.k1 * r * r);
        double rc = undistort_radial_k1(rad.k1, rd, rmax);
        double rn, ry0;
        double truth[2] = { r, 0.0 };
        if (!undistort(&rad, rd, 0.0, &rn, &ry0, 50, 1e-14, truth, NULL, NULL)) { bound_ok = 0; continue; }
        double fx, fy;
        distort(&rad, rn, ry0, &fx, &fy);
        double bound = hypot(fx - rd, fy) / gmin; /* 后验上界 |r-r*| <= |F|/min g' */
        double err = fabs(rn - rc);
        if (err > max_err) max_err = err;
        if (bound > max_bound) max_bound = bound;
        if (err > bound + 1e-12) bound_ok = 0;
    }
    printf("analytic_max_err=%.3e\n", max_err);
    printf("analytic_err_bound=%.3e\n", max_bound);
    check("analytic_within_bound", bound_ok);

    /* 5. 可辨识性: 不同采样构型下的不可辨识系数组 */
    double gx[256], gy[256];
    int n = 0;
    for (int i = 0; i < 16; i++)
        for (int j = 0; j < 16; j++) {
            gx[n] = -0.8 + 1.6 * i / 15.0;
            gy[n] = -0.8 + 1.6 * j / 15.0;
            n++;
        }
    check("ident_grid", ident_report("grid", gx, gy, n) == 0);

    double rx2[64], ry2[64];
    for (int i = 0; i < 64; i++) {
        rx2[i] = 0.5 * cos(2.0 * PI * i / 64.0);
        ry2[i] = 0.5 * sin(2.0 * PI * i / 64.0);
    }
    check("ident_ring", ident_report("ring", rx2, ry2, 64) == 2); /* {k1,k2,k3} 二维零空间 */

    double px1[1] = { 0.3 }, py1[1] = { 0.4 };
    check("ident_point", ident_report("point", px1, py1, 1) == 3);

    double sx2[25], sy2[25];
    n = 0;
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) {
            sx2[n] = -0.02 + 0.04 * i / 4.0;
            sy2[n] = -0.02 + 0.04 * j / 4.0;
            n++;
        }
    check("ident_smallr", ident_report("smallr", sx2, sy2, n) >= 1); /* k3 数值不可辨识 */

    /* 6. 未收敛像素计数 (1280x960, f=700) */
    long bad_wide = sweep(&wide, 1280, 960, 700.0, rs_wide);
    long bad_mild = sweep(&mild, 1280, 960, 700.0, rs_mild);
    printf("nonconverged_pixels_wide=%ld\n", bad_wide);
    printf("nonconverged_pixels_mild=%ld\n", bad_mild);
    check("mild_all_converge", bad_mild == 0);
    check("wide_corner_fail", bad_wide > 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
