
/* mis: multiple importance sampling with the balance heuristic.
   Two unbiased strategies over the hemisphere:
     A: BRDF-lobe importance sampling  pA(w) = (e1+1)/(2pi) (R.w)^e1
     B: light-lobe importance sampling pB(w) = (e2+1)/(2pi) (L.w)^e2
   Integrand (glossy reflection of a lobe-shaped light off a Lambert patch):
     f(w) = rho/pi * (R.w)^e1 * (L.w)^e2 * (n.w)
   Estimators:
     single:  X_i = f(w)/p_i(w)
     MIS:     X   = w_A f/pA + w_B f/pB,  w_i = p_i/(pA+pB) (balance heuristic)
   Unbiasedness is checked against (1) an analytic integral on a uniform-cap
   scene and (2) deterministic midpoint quadrature on the two-lobe scene.
   A weight-sum guard rejects any combine whose weights do not sum to 1. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

typedef struct { double x, y, z; } V3;

static V3 v3(double x, double y, double z) { V3 r = { x, y, z }; return r; }
static double dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 cross3(V3 a, V3 b) {
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static V3 norm3(V3 a) {
    double l = sqrt(dot3(a, a));
    return l > 0 ? v3(a.x / l, a.y / l, a.z / l) : v3(0, 0, 0);
}
/* rotate a direction given in a frame (t, b, n) into world coordinates */
static V3 from_frame(V3 t, V3 b, V3 n, V3 d) {
    return v3(t.x * d.x + b.x * d.y + n.x * d.z,
              t.y * d.x + b.y * d.y + n.y * d.z,
              t.z * d.x + b.z * d.y + n.z * d.z);
}
/* frame whose z-axis is c; x-axis lies in the plane of c and the normal z */
static void frame_around(V3 c, V3 *t, V3 *b, V3 *n) {
    V3 z0 = v3(0, 0, 1);
    V3 x = norm3(cross3(z0, c));
    if (dot3(x, x) == 0.0) x = v3(1, 0, 0);
    *t = x; *b = cross3(c, x); *n = c;
}

/* xorshift64* RNG, independent streams per trial via splitmix64 seeding */
static unsigned long long rng_state;
static unsigned long long rng_u64(void) {
    unsigned long long x = rng_state;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    rng_state = x;
    return x * 0x2545F4914F6CDD1DULL;
}
static double rng_uniform(void) { return (rng_u64() >> 11) * 0x1.0p-53; }
static unsigned long long splitmix64(unsigned long long *s) {
    unsigned long long z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* ---------------- scene ---------------- */

typedef struct {
    V3 n;      /* surface normal */
    V3 R;      /* BRDF lobe center (reflection direction) */
    V3 L;      /* light lobe center */
    double e1; /* BRDF lobe exponent */
    double e2; /* light lobe exponent */
    double rho;
} Scene;

static double integrand(const Scene *s, V3 w) {
    double dR = dot3(s->R, w), dL = dot3(s->L, w), dn = dot3(s->n, w);
    if (dR <= 0.0 || dL <= 0.0 || dn <= 0.0) return 0.0;
    return s->rho / PI * pow(dR, s->e1) * pow(dL, s->e2) * dn;
}
static double pdf_brdf(const Scene *s, V3 w) {
    double dR = dot3(s->R, w);
    if (dR <= 0.0) return 0.0;
    return (s->e1 + 1.0) / (2.0 * PI) * pow(dR, s->e1);
}
static double pdf_light(const Scene *s, V3 w) {
    double dL = dot3(s->L, w);
    if (dL <= 0.0) return 0.0;
    return (s->e2 + 1.0) / (2.0 * PI) * pow(dL, s->e2);
}

/* power-cosine lobe sampling around frame z-axis:
   cos theta = u1^(1/(e+1)), phi = 2 pi u2 */
static V3 sample_lobe(V3 center, double e, double u1, double u2) {
    double ct = pow(u1, 1.0 / (e + 1.0));
    double st = sqrt(fmax(0.0, 1.0 - ct * ct));
    double phi = 2.0 * PI * u2;
    V3 local = v3(st * cos(phi), st * sin(phi), ct);
    V3 t, b, n;
    frame_around(center, &t, &b, &n);
    return from_frame(t, b, n, local);
}

/* ---------------- balance-heuristic combine with weight-sum guard ---------------- */

typedef struct {
    double value;
    int ok;            /* weights summed to 1 and at least one pdf was nonzero */
    double weight_sum; /* actual weight sum, for diagnostics */
} Combine;

static Combine mis_combine(double f, double pA, double pB) {
    Combine c;
    c.value = 0.0;
    c.weight_sum = 0.0;
    double sum = pA + pB;
    if (sum <= 0.0) {
        /* unsampleable point: unbiased combine requires f == 0 here */
        c.ok = (f == 0.0);
        return c;
    }
    double wA = pA / sum, wB = pB / sum;
    c.weight_sum = wA + wB;
    c.ok = fabs(c.weight_sum - 1.0) < 1e-12;
    /* wA*f/pA + wB*f/pB simplifies to the closed form f/(pA+pB) */
    c.value = c.ok ? f / sum : 0.0;
    return c;
}

/* ---------------- deterministic quadrature reference ---------------- */

/* midpoint rule over the hemisphere:
   I = sum f(w) dw, dw = dcos(theta) dphi, w=(sin phi st, cos phi st, ct).
   JA/JB are the weighted partial integrals int f p_i/(pA+pB) dw, JA+JB = I. */
typedef struct {
    double mean;
    double m2A;
    double m2B;
    double m2mis; /* int f^2/(pA+pB) dw */
    double JA;
    double JB;
} Ref;

static Ref quadrature(const Scene *s, int n_ct, int n_phi) {
    Ref r; memset(&r, 0, sizeof(r));
    double dct = 1.0 / n_ct, dphi = 2.0 * PI / n_phi;
    double dw = dct * dphi;
    for (int i = 0; i < n_ct; i++) {
        double ct = (i + 0.5) * dct;
        double st = sqrt(fmax(0.0, 1.0 - ct * ct));
        for (int j = 0; j < n_phi; j++) {
            double phi = (j + 0.5) * dphi;
            V3 w = v3(st * cos(phi), st * sin(phi), ct);
            double f = integrand(s, w);
            double pA = pdf_brdf(s, w), pB = pdf_light(s, w);
            double sum = pA + pB;
            r.mean += f * dw;
            if (pA > 0.0) r.m2A += f * f / pA * dw;
            if (pB > 0.0) r.m2B += f * f / pB * dw;
            if (sum > 0.0) {
                r.m2mis += f * f / sum * dw;
                r.JA += f * pA / sum * dw;
                r.JB += f * pB / sum * dw;
            }
        }
    }
    return r;
}

/* ---------------- Monte Carlo experiment ---------------- */

#define TRIALS 4096

typedef struct {
    double mean;      /* estimator mean across trials */
    double variance;  /* variance of one trial estimate */
    double bias;      /* mean - reference */
    int guard_fails;  /* weight-sum violations */
} Result;

/* Equal total budget: every estimator spends 2m samples per trial.
   mode 0: strategy A only (2m draws); mode 1: strategy B only (2m draws);
   mode 2: balance-heuristic MIS (m paired draws). */
static Result run_experiment(const Scene *s, int m, int mode, double reference) {
    Result out; memset(&out, 0, sizeof(out));
    double sum = 0.0, sum2 = 0.0;
    int budget = 2 * m;
    unsigned long long seed = 0x12345678ULL ^ ((unsigned long long)m << 32) ^
                              ((unsigned long long)mode << 16);
    for (int t = 0; t < TRIALS; t++) {
        rng_state = splitmix64(&seed);
        double acc = 0.0;
        if (mode == 0) {
            for (int k = 0; k < budget; k++) {
                V3 w = sample_lobe(s->R, s->e1, rng_uniform(), rng_uniform());
                double f = integrand(s, w), p = pdf_brdf(s, w);
                acc += p > 0.0 ? f / p : 0.0;
            }
        } else if (mode == 1) {
            for (int k = 0; k < budget; k++) {
                V3 w = sample_lobe(s->L, s->e2, rng_uniform(), rng_uniform());
                double f = integrand(s, w), p = pdf_light(s, w);
                acc += p > 0.0 ? f / p : 0.0;
            }
        } else {
            for (int k = 0; k < m; k++) {
                V3 wa = sample_lobe(s->R, s->e1, rng_uniform(), rng_uniform());
                V3 wb = sample_lobe(s->L, s->e2, rng_uniform(), rng_uniform());
                Combine ca = mis_combine(integrand(s, wa), pdf_brdf(s, wa),
                                         pdf_light(s, wa));
                Combine cb = mis_combine(integrand(s, wb), pdf_brdf(s, wb),
                                         pdf_light(s, wb));
                if (!ca.ok || !cb.ok) out.guard_fails++;
                acc += ca.value + cb.value;
            }
        }
        /* single strategies average over their 2m draws; the Veach multi-sample
           MIS estimator sums weighted contributions over its m pairs (its
           expectation is already I), so it is normalized by the pair count. */
        double est = (mode == 2) ? acc / m : acc / budget;
        sum += est;
        sum2 += est * est;
    }
    out.mean = sum / TRIALS;
    double var = sum2 / TRIALS - out.mean * out.mean;
    out.variance = var * TRIALS / (TRIALS - 1);
    out.bias = out.mean - reference;
    return out;
}

/* ---------------- analytic reference: uniform light cap, Lambert BRDF ---------------- */

/* scene 1: f(w) = rho/pi * cos  inside a cap of half-angle tc centered at L.
   analytic integral: I = rho (1-cos^2 tc) cos gamma, gamma = angle(n, L).
   strategies: A cosine hemisphere pA=cos/pi, B uniform cap pB=1/Omega. */
typedef struct { double gamma, tc, rho; } CapScene;

static double cap_integrand(const CapScene *c, V3 w) {
    V3 n = v3(0, 0, 1);
    V3 L = v3(sin(c->gamma), 0, cos(c->gamma));
    if (dot3(L, w) < cos(c->tc) || dot3(n, w) <= 0.0) return 0.0;
    return c->rho / PI * dot3(n, w);
}
static double cap_pa(V3 w) {
    V3 n = v3(0, 0, 1);
    double dn = dot3(n, w);
    return dn > 0.0 ? dn / PI : 0.0;
}
static double cap_pb(const CapScene *c, V3 w) {
    V3 L = v3(sin(c->gamma), 0, cos(c->gamma));
    double omega = 2.0 * PI * (1.0 - cos(c->tc));
    return dot3(L, w) >= cos(c->tc) ? 1.0 / omega : 0.0;
}
static V3 cap_sample_light(const CapScene *c, double u1, double u2) {
    double ct = 1.0 - u1 * (1.0 - cos(c->tc));
    double st = sqrt(fmax(0.0, 1.0 - ct * ct));
    double phi = 2.0 * PI * u2;
    V3 local = v3(st * cos(phi), st * sin(phi), ct);
    V3 L = v3(sin(c->gamma), 0, cos(c->gamma));
    V3 t, b, n;
    frame_around(L, &t, &b, &n);
    return from_frame(t, b, n, local);
}
static V3 cap_sample_cosine(double u1, double u2) {
    double r = sqrt(u1), phi = 2.0 * PI * u2;
    return v3(r * cos(phi), r * sin(phi), sqrt(fmax(0.0, 1.0 - u1)));
}

static double cap_analytic(const CapScene *c) {
    return c->rho * (1.0 - cos(c->tc) * cos(c->tc)) * cos(c->gamma);
}

/* ---------------- test harness ---------------- */

static int failures = 0;
static void check(const char *name, int ok) {
    printf("  %s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* ---- scene 1: analytic unbiasedness, uniform cap + Lambert ---- */
    printf("scene1: uniform light cap, analytic integral\n");
    CapScene cap = { 60.0 * PI / 180.0, 25.0 * PI / 180.0, 0.8 };
    double ref1 = cap_analytic(&cap);
    int M1 = 256;
    unsigned long long seed = 0xABCDEFULL;
    double sA = 0.0, sB = 0.0, sM = 0.0;
    int cap_guard_fails = 0;
    for (int k = 0; k < M1; k++) {
        rng_state = splitmix64(&seed);
        V3 wa = cap_sample_cosine(rng_uniform(), rng_uniform());
        V3 wb = cap_sample_light(&cap, rng_uniform(), rng_uniform());
        double fa = cap_integrand(&cap, wa), pa = cap_pa(wa), pb = cap_pb(&cap, wa);
        double fb = cap_integrand(&cap, wb), qa = cap_pa(wb), qb = cap_pb(&cap, wb);
        sA += pa > 0.0 ? fa / pa : 0.0;
        sB += qb > 0.0 ? fb / qb : 0.0;
        Combine ca = mis_combine(fa, pa, pb), cb = mis_combine(fb, qa, qb);
        if (!ca.ok || !cb.ok) cap_guard_fails++;
        sM += ca.value + cb.value;
    }
    double estA = sA / M1, estB = sB / M1, estM = sM / M1;
    /* 5-sigma tolerances from conservative per-sample variance bounds */
    double tolA = 5.0 * sqrt(fmax(ref1, 1e-9) / M1);
    double tolB = 5.0 * sqrt(0.05 / M1);
    double tolM = 5.0 * sqrt(0.05 / (2 * M1));
    printf("  analytic I=%.6f  estA=%.6f  estB=%.6f  estMIS=%.6f\n",
           ref1, estA, estB, estM);
    check("cap_analytic_positive", ref1 > 0.0);
    check("cap_strategyA_unbiased", fabs(estA - ref1) < tolA);
    check("cap_strategyB_unbiased", fabs(estB - ref1) < tolB);
    check("cap_mis_unbiased", fabs(estM - ref1) < tolM);
    check("cap_weight_guard_clean", cap_guard_fails == 0);

    /* guard must fire on an unsampleable nonzero contribution */
    Combine bad = mis_combine(1.0, 0.0, 0.0);
    check("guard_reports_zero_pdf", bad.ok == 0);
    Combine good = mis_combine(0.3, 0.2, 0.8);
    check("guard_passes_unit_weights", good.ok == 1 && fabs(good.weight_sum - 1.0) < 1e-15);
    check("guard_mis_value", fabs(good.value - 0.3) < 1e-15);

    /* ---- scene 2: two-lobe product, quadrature truth, variance comparison ---- */
    printf("scene2: separated BRDF/light lobes, quadrature reference\n");
    Scene sc;
    sc.n = v3(0, 0, 1);
    sc.R = norm3(v3(sin(35.0 * PI / 180.0), 0, cos(35.0 * PI / 180.0)));
    sc.L = norm3(v3(sin(65.0 * PI / 180.0), 0, cos(65.0 * PI / 180.0)));
    sc.e1 = 8.0;
    sc.e2 = 16.0;
    sc.rho = 0.8;

    Ref ref = quadrature(&sc, 1024, 2048);
    double I = ref.mean;
    double vA_th = ref.m2A - I * I;
    double vB_th = ref.m2B - I * I;
    /* one paired BH draw X = va+vb has E[X]=JA+JB=I and
       E[X^2]=MM+2 JA JB, hence Var(X)=MM-JA^2-JB^2 (>=0 by Cauchy-Schwarz).
       A trial with m pairs spends 2m samples; its per-sample constant is
       2 Var(X), matching the singles' per-sample constants below. */
    double v_pair = ref.m2mis - ref.JA * ref.JA - ref.JB * ref.JB;
    double vM_th = 2.0 * v_pair;
    printf("  quadrature I=%.7f  (JA+JB=%.7f, partition check)\n",
           I, ref.JA + ref.JB);
    printf("  per-sample variance constants (equal budget 2m):\n");
    printf("    A=%.7f  B=%.7f  MIS=%.7f  MIS/best=%.4f\n",
           vA_th, vB_th, vM_th, vM_th / fmin(vA_th, vB_th));
    check("weight_partition_identity", fabs(ref.JA + ref.JB - I) < 1e-10);
    check("quadrature_converged", I > 0.0 && vA_th > 0.0 && vB_th > 0.0 && v_pair > 0.0);
    check("theory_mis_beats_both", vM_th < vA_th && vM_th < vB_th);

    /* empirical: equal total sample budget 2m for every estimator */
    printf("  empirical %d trials, equal budget 2m:\n", TRIALS);
    printf("  %6s %10s %10s %10s | %9s %9s %9s | %s\n",
           "m", "varA", "varB", "varMIS", "biasA", "biasB", "biasM", "guard");
    int all_unbiased = 1, mis_wins = 1;
    int ms[] = { 16, 64, 256 };
    for (unsigned mi = 0; mi < sizeof(ms) / sizeof(ms[0]); mi++) {
        int m = ms[mi];
        Result rA = run_experiment(&sc, m, 0, I);
        Result rB = run_experiment(&sc, m, 1, I);
        Result rM = run_experiment(&sc, m, 2, I);
        /* 5-sigma tolerance on the measured bias */
        double tolA = 5.0 * sqrt(rA.variance / TRIALS);
        double tolB = 5.0 * sqrt(rB.variance / TRIALS);
        double tolM = 5.0 * sqrt(rM.variance / TRIALS);
        int okA = fabs(rA.bias) < tolA;
        int okB = fabs(rB.bias) < tolB;
        int okM = fabs(rM.bias) < tolM;
        int win = rM.variance < rA.variance && rM.variance < rB.variance;
        all_unbiased &= okA && okB && okM;
        mis_wins &= win;
        printf("  %6d %10.7f %10.7f %10.7f | %9.2e %9.2e %9.2e | fails=%d %s\n",
               m, rA.variance, rB.variance, rM.variance,
               rA.bias, rB.bias, rM.bias, rM.guard_fails,
               win ? "MIS_LOWEST" : "NO_WIN");
        check("guard_clean_at_m", rM.guard_fails == 0);
    }
    check("all_estimators_unbiased_5sigma", all_unbiased);
    check("mis_variance_lowest_every_m", mis_wins);

    /* ---- theoretical bounds (all constants are per-sample, budget-normalized) ---- */
    printf("theory:\n");
    /* naive equal-weight MIS: weight 1/2 each, per-sample second moment (MA+MB)/2.
       Cauchy-Schwarz on q=(pA+pB)/2 gives 2*MM <= MA+MB, hence BH <= naive. */
    double v_naive = 0.5 * (ref.m2A + ref.m2B) - I * I;
    /* mixture sampling from q=(pA+pB)/2: per-sample constant 2MM-I^2.
       BH vs mixture reduces to JA^2+JB^2 >= 2 JA JB (C-S), equality iff JA=JB. */
    double v_mix = 2.0 * ref.m2mis - I * I;
    /* Veach's general balance-heuristic bound (n=2): V[BH] <= V_i + I^2. */
    double veach = fmin(vA_th, vB_th) + I * I;
    printf("  lower bound 0 <= V[MIS]=%.7f (zero variance iff p*=f/I sampled)\n", vM_th);
    printf("  BH vs naive equal-weight MIS: %.7f <= %.7f  [C-S: 2MM <= MA+MB]\n",
           vM_th, v_naive);
    printf("  BH vs mixture q=(pA+pB)/2:    %.7f <= %.7f  [C-S: JA^2+JB^2 >= 2JAJB]\n",
           vM_th, v_mix);
    printf("  Veach BH bound: %.7f <= min(V_A,V_B)+I^2=%.7f\n", vM_th, veach);
    check("nonnegative_lower_bound", vM_th >= 0.0);
    check("bh_beats_naive_equal_weight", vM_th <= v_naive + 1e-12);
    check("bh_not_worse_than_mixture", vM_th <= v_mix + 1e-12);
    check("veach_bound_holds", vM_th <= veach + 1e-12);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
