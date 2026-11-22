
/* sample: stratified and low-discrepancy sampling with variance estimate. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N_MAX 256

typedef struct { double mean; double var; int n; } Stat;

static unsigned radical_inverse(unsigned i, unsigned base) {
    unsigned rev = 0, base_n = base;
    while (i > 0) { rev += (i % base) * (base_n / base); i /= base; base_n *= base; }
    return rev;
}

static double halton(unsigned i, unsigned base) { return (double)radical_inverse(i, base) / 4294967296.0 * base; }
static double halton2(unsigned i, unsigned base) {
    double f = 1.0, r = 0.0; unsigned k = i;
    while (k > 0) { f /= base; r += f * (k % base); k /= base; }
    return r;
}

/* estimate integral of f over [0,1) with stratified samples */
static Stat stratify(int strata, double (*f)(double)) {
    Stat s; memset(&s, 0, sizeof(s));
    double sum = 0.0, sum2 = 0.0;
    int n = 0;
    for (int i = 0; i < strata; i++) {
        double u = ((double)i + 0.5) / (double)strata;
        double v = f(u);
        sum += v; sum2 += v * v; n++;
    }
    s.n = n;
    s.mean = sum / n;
    double var = (n > 1) ? (sum2 / n - s.mean * s.mean) * n / (n - 1) : 0.0;
    s.var = var > 0 ? var : 0;
    return s;
}

static double f_quad(double u) { return u * u; }
static double f_const(double u) { (void)u; return 1.0; }

/* Light-emission power compensation.
   N lights, power P_i, each contributing c_i per unit power. Emission picks
   light i with probability q_i (target distribution, sum q_i = 1) and forms
   the estimator X = c_i * w_i with per-photon weight w_i = (P_i / P_total)/q_i,
   so E_q[X] = sum_i c_i P_i / P_total (unbiased for the power-normalized
   integral). Uniform emission uses q_u,i = 1/n with its own weights; the two
   variances are compared.

   Rejected: target not integrating to one, any negative weight. */

#define LIGHT_MAX 64

/* weights for sampling distribution q over powers p; 0 ok, -1 reject */
static int emit_weights(int n, const double p[LIGHT_MAX], const double q[LIGHT_MAX],
                        double total_power, double w[LIGHT_MAX]) {
    if (n <= 0 || n > LIGHT_MAX || total_power <= 0.0) return -1;
    double sum_q = 0.0;
    for (int i = 0; i < n; i++) {
        if (q[i] < 0.0) return -1;
        sum_q += q[i];
    }
    if (fabs(sum_q - 1.0) > 1e-6) return -1;   /* target must integrate to one */
    for (int i = 0; i < n; i++) {
        if (q[i] <= 0.0) {
            if (p[i] > 0.0) return -1;         /* power at a zero-probability light */
            w[i] = 0.0;
            continue;
        }
        double wi = (p[i] / total_power) / q[i];
        if (wi < 0.0) return -1;               /* reject negative weight */
        w[i] = wi;
    }
    return 0;
}

/* variance of estimator X = c_i w_i when i is drawn from q */
static double emit_variance(int n, const double c[LIGHT_MAX], const double w[LIGHT_MAX],
                            const double q[LIGHT_MAX]) {
    double mean = 0.0;
    for (int i = 0; i < n; i++) mean += q[i] * (c[i] * w[i]);
    double var = 0.0;
    for (int i = 0; i < n; i++) {
        double t = c[i] * w[i] - mean;
        var += q[i] * t * t;
    }
    return var;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    check("halton_in_unit", halton2(1, 2) >= 0.0 && halton2(1, 2) < 1.0);
    check("halton_first", fabs(halton2(0, 2)) < 1e-12);

    Stat a = stratify(64, f_quad);
    check("quad_mean_close", fabs(a.mean - 1.0 / 3.0) < 1e-3);

    Stat b = stratify(16, f_quad);
    Stat c = stratify(64, f_quad);
    check("variance_drops", c.var <= b.var + 1e-6);

    Stat d = stratify(32, f_const);
    check("const_exact", fabs(d.mean - 1.0) < 1e-12 && d.var < 1e-12);

    check("halton_distinct", halton2(1, 2) != halton2(2, 2));
    check("halton_monotone_base", halton2(1, 2) > halton2(0, 2));

    /* --- power-compensated emission --- */
    /* 4 lights, very different powers; the dim ones are the noise source at
       fixed photon count under uniform emission. */
    double p[LIGHT_MAX] = { 1.0, 10.0, 100.0, 1000.0 };
    int n = 4;
    double total = 0.0;
    for (int i = 0; i < n; i++) total += p[i];

    /* target distribution = power-normalized, integrates to one by construction */
    double q[LIGHT_MAX];
    for (int i = 0; i < n; i++) q[i] = p[i] / total;
    double qsum = 0.0;
    for (int i = 0; i < n; i++) qsum += q[i];
    check("target_integrates_one", fabs(qsum - 1.0) < 1e-12);

    /* integrand per light: geometric visibility to a surface, not constant, so
       neither scheme is degenerate and the variance drop is meaningful. */
    double ci[LIGHT_MAX] = { 0.9, 0.5, 0.2, 0.1 };

    double w[LIGHT_MAX], wu[LIGHT_MAX];
    check("weights_ok", emit_weights(n, p, q, total, w) == 0);

    double qu[LIGHT_MAX];
    for (int i = 0; i < n; i++) qu[i] = 1.0 / n;
    check("weights_uniform_ok", emit_weights(n, p, qu, total, wu) == 0);

    /* unbiasedness: E_q[c*w] must equal E_uniform[c*w] = sum c_i P_i / total */
    double ex_q = 0.0, ex_u = 0.0;
    for (int i = 0; i < n; i++) { ex_q += q[i] * (ci[i] * w[i]); ex_u += qu[i] * (ci[i] * wu[i]); }
    check("unbiased_match", fabs(ex_q - ex_u) < 1e-12);

    /* report the variance drop versus uniform emission */
    double vq = emit_variance(n, ci, w, q);
    double vu = emit_variance(n, ci, wu, qu);
    printf("var_target=%.12f\n", vq);
    printf("var_uniform=%.12f\n", vu);
    printf("var_ratio=%.12f\n", vu > 0.0 ? vq / vu : 0.0);
    check("variance_drops_vs_uniform", vq < vu);

    /* with constant integrand and power-proportional emission the estimator is
       exactly constant, i.e. the optimum of this scheme. */
    double cone[LIGHT_MAX] = { 1.0, 1.0, 1.0, 1.0 };
    double wc[LIGHT_MAX];
    check("weights_const_ok", emit_weights(n, p, q, total, wc) == 0);
    double vc = emit_variance(n, cone, wc, q);
    check("optimum_variance_zero", vc < 1e-12);

    /* reject: target does not integrate to one */
    double qbad[LIGHT_MAX];
    for (int i = 0; i < n; i++) qbad[i] = q[i];
    qbad[0] += 0.5;
    check("reject_target_not_one", emit_weights(n, p, qbad, total, w) != 0);

    /* reject: negative weight (a light emits with negative power) */
    double pneg[LIGHT_MAX];
    for (int i = 0; i < n; i++) pneg[i] = p[i];
    pneg[1] = -5.0;
    check("reject_negative_weight", emit_weights(n, pneg, q, total, w) != 0);

    /* reject: positive power placed on a zero-probability light */
    double qzero[LIGHT_MAX];
    double qz_rest = 0.0;
    for (int i = 0; i < n; i++) { qzero[i] = (i == 3) ? 0.0 : q[i]; qz_rest += qzero[i]; }
    for (int i = 0; i < n; i++) qzero[i] /= qz_rest;    /* others renormalized to 1 */
    qzero[3] = 0.0;
    check("reject_power_at_zero_q", emit_weights(n, p, qzero, total, w) != 0);

    /* --- power-compensated emission over an 8x8 light grid --- */
    /* Power falls off from one corner; target q is power-normalized. Integrand
       is a broad falloff in the other direction, so the optimum q is not simply
       power-proportional and the variance drop is a nontrivial number. */
    {
        int g = 8, gn = g * g;
        double gp[LIGHT_MAX], gq[LIGHT_MAX], gqu[LIGHT_MAX], gtm[LIGHT_MAX];
        double gw[LIGHT_MAX], gwu[LIGHT_MAX];
        double gtotal = 0.0;
        for (int j = 0; j < g; j++)
            for (int i = 0; i < g; i++) {
                int k = j * g + i;
                gp[k] = 1.0 / (1.0 + (double)(i * i + j * j));
                gtm[k] = 0.05 + 0.95 * exp(-(double)((g - 1 - i) * (g - 1 - i) + j * j) / 32.0);
                gtotal += gp[k];
            }
        double gsum = 0.0;
        for (int k = 0; k < gn; k++) gsum += gp[k];
        double gqsum = 0.0;
        for (int k = 0; k < gn; k++) { gq[k] = gp[k] / gsum; gqu[k] = 1.0 / gn; gqsum += gq[k]; }
        check("grid_target_integrates_one", fabs(gqsum - 1.0) < 1e-12);
        check("grid_weights_ok", emit_weights(gn, gp, gq, gtotal, gw) == 0);
        check("grid_uniform_ok", emit_weights(gn, gp, gqu, gtotal, gwu) == 0);

        double gex = 0.0, gexu = 0.0;
        for (int k = 0; k < gn; k++) {
            gex += gq[k] * (gtm[k] * gw[k]);
            gexu += gqu[k] * (gtm[k] * gwu[k]);
        }
        check("grid_unbiased_match", fabs(gex - gexu) < 1e-12);

        double gvq = emit_variance(gn, gtm, gw, gq);
        double gvu = emit_variance(gn, gtm, gwu, gqu);
        printf("grid_var_target=%.12f\n", gvq);
        printf("grid_var_uniform=%.12f\n", gvu);
        printf("grid_var_ratio=%.12f\n", gvu > 0.0 ? gvq / gvu : 0.0);
        check("grid_variance_drops", gvq < gvu);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
