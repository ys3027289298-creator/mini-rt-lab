/* compensate: inverse-probability compensation for low-probability bright
   events, with a closed-form post-compensation variance and explicit bias
   accounting when the compensation weight is truncated.

   Model
   -----
   The domain is split into K cells.  Cell i has area fraction p_i (the
   reference measure: what a natural / area-proportional sampler draws) and
   carries radiance f_i.  A ray sample from the natural measure estimates

       I = sum_i p_i f_i                                        (mean radiance)

   with per-sample value f_i and variance

       Var_nat = sum_i p_i f_i^2 - I^2.

   A single tiny but very bright cell (p small, f huge) makes p_i f_i^2 large
   and wrecks Var_nat, even though its mean contribution is modest.

   Compensation
   ------------
   Draw cell i from a proposal q instead, and compensate each sample by the
   inverse-probability weight w_i = p_i / q_i:

       X_i = w_i f_i = (p_i / q_i) f_i.

   For any q with q_i > 0 wherever p_i f_i > 0,

       E_q[X] = sum_i q_i (p_i/q_i) f_i = sum_i p_i f_i = I,

   so the compensated estimator is unbiased under *every* proposal.  Its
   variance follows in closed form from E[X^2]:

       Var_q[X] = sum_i q_i (p_i f_i / q_i)^2 - I^2
                = sum_i p_i^2 f_i^2 / q_i - I^2.

   The variance is zero when q_i is proportional to p_i f_i (the energy
   proposal): then w_i f_i = I for every i and each sample returns I exactly.
   The factor 1/q_i is what keeps the estimate unbiased when the proposal is
   tilted; dropping it leaves the biased "raw" estimator whose mean is
   sum_i q_i f_i, not I.

   Truncation
   ----------
   w_i = p_i/q_i is unbounded as q_i -> 0, and a bad energy estimate (a pilot
   that underestimates the bright cell) pushes q_i there.  Clamping the weight
   at W replaces X_i by min(w_i, W) f_i, which is still unbiased for a
   different integral:

       E[X_clamped] = I - sum_{w_i > W} p_i f_i (1 - W/w_i) = I - lost(W).

   The bias is exactly the energy of the clipped events; this module computes
   that closed form, reports which cells were clipped and how much each
   contributes, so the bias is localised rather than silent.

   Everything below is evaluated as an exact sum over the discrete
   distribution (no sampling), except the two marked Monte Carlo checks. */

#include <stdio.h>
#include <math.h>
#include <string.h>

#define K 64
#define MC_N 200000L

typedef struct {
    double p[K];    /* reference (natural) probability of cell i */
    double f[K];    /* radiance carried by cell i */
} Model;

typedef struct {
    double mean;    /* E[estimator]                                   */
    double var;     /* Var[estimator]                                 */
    double mse;     /* var + (mean - I)^2                             */
    double lost;    /* sum p_i f_i (1 - W/w_i) over clipped cells     */
    double wmax;    /* largest weight p_i/q_i seen (before clamping)  */
    double wmax_used; /* largest weight actually applied             */
    int nclip;      /* number of cells with w_i > W                   */
    int worst;      /* clipped cell contributing the most lost energy */
    double worst_lost;
    int unbiased;   /* mean == I within tolerance                     */
} Moments;

/* ------------------------------------------------------------------ */

static double integral(const Model* m) {
    double I = 0.0;
    for (int i = 0; i < K; i++) I += m->p[i] * m->f[i];
    return I;
}

/* q_i proportional to p_i f_i.  Normalised so sum q = 1. */
static void energy_proposal(const Model* m, double q[K]) {
    double I = integral(m);
    for (int i = 0; i < K; i++) q[i] = (I > 0.0) ? m->p[i] * m->f[i] / I : 1.0 / K;
}

/* q = (1 - lam) * p + lam * energy.  lam = 0 is the natural sampler,
   lam = 1 the zero-variance energy sampler. */
static void mixed_proposal(const Model* m, double lam, double q[K]) {
    double qe[K], sp = 0.0;
    energy_proposal(m, qe);
    for (int i = 0; i < K; i++) sp += m->p[i];
    for (int i = 0; i < K; i++)
        q[i] = (1.0 - lam) * (m->p[i] / sp) + lam * qe[i];
}

/* Proposal from a measured energy estimate g: q_i proportional to p_i g_i.
   A pilot that underestimates the bright cell gives it too little mass and
   inflates the weight p_i/q_i. */
static void pilot_proposal(const Model* m, const double g[K], double q[K]) {
    double s = 0.0;
    for (int i = 0; i < K; i++) s += m->p[i] * g[i];
    for (int i = 0; i < K; i++) q[i] = (s > 0.0) ? m->p[i] * g[i] / s : 1.0 / K;
}

/* Compensated estimator: draw from q, average min(p_i/q_i, W) f_i.
   W <= 0 disables the clamp. */
static Moments compensated(const Model* m, const double q[K], double W) {
    Moments r;
    memset(&r, 0, sizeof(r));
    r.worst = -1;
    double I = integral(m), e1 = 0.0, e2 = 0.0;

    for (int i = 0; i < K; i++) {
        double pi = m->p[i], qi = q[i];
        if (pi <= 0.0 || qi <= 0.0 || m->f[i] == 0.0) continue;
        double w = pi / qi;
        if (w > r.wmax) r.wmax = w;
        double used = w;
        if (W > 0.0 && w > W) {
            double l = pi * m->f[i] * (1.0 - W / w);
            r.lost += l;
            r.nclip++;
            if (l > r.worst_lost) { r.worst_lost = l; r.worst = i; }
            used = W;
        }
        if (used > r.wmax_used) r.wmax_used = used;
        double x = used * m->f[i];
        e1 += qi * x;
        e2 += qi * x * x;
    }
    r.mean = e1;
    r.var = e2 - e1 * e1;
    if (r.var < 0.0) r.var = 0.0;
    r.mse = r.var + (r.mean - I) * (r.mean - I);
    r.unbiased = fabs(r.mean - I) <= 1e-12 * (fabs(I) + 1.0);
    return r;
}

/* Uncompensated ("raw") estimator: draw from q, average f_i with no weight.
   E = sum q_i f_i, which equals I only when q == p. */
static Moments uncompensated(const Model* m, const double q[K]) {
    Moments r;
    memset(&r, 0, sizeof(r));
    r.worst = -1;
    double I = integral(m), e1 = 0.0, e2 = 0.0;
    for (int i = 0; i < K; i++) {
        e1 += q[i] * m->f[i];
        e2 += q[i] * m->f[i] * m->f[i];
    }
    r.mean = e1;
    r.var = e2 - e1 * e1;
    if (r.var < 0.0) r.var = 0.0;
    r.mse = r.var + (r.mean - I) * (r.mean - I);
    r.unbiased = fabs(r.mean - I) <= 1e-12 * (fabs(I) + 1.0);
    return r;
}

/* Closed-form post-compensation variance, sum p_i^2 f_i^2 / q_i - I^2. */
static double var_closed_form(const Model* m, const double q[K]) {
    double I = integral(m), s = 0.0;
    for (int i = 0; i < K; i++)
        if (q[i] > 0.0) s += m->p[i] * m->p[i] * m->f[i] * m->f[i] / q[i];
    return s - I * I;
}

/* Share of the compensated estimator's second moment, sum p_i^2 f_i^2 / q_i,
   contributed by cell k.  This is the quantity whose blow-up at small q_k is
   the firefly problem. */
static double second_moment_share(const Model* m, const double q[K], int k) {
    double s = 0.0;
    for (int i = 0; i < K; i++)
        if (q[i] > 0.0) s += m->p[i] * m->p[i] * m->f[i] * m->f[i] / q[i];
    if (s <= 0.0) return 0.0;
    return m->p[k] * m->p[k] * m->f[k] * m->f[k] / q[k] / s;
}

/* Report the sources of truncation bias: the clipped cells with the largest
   lost energy, plus totals, so the bias can be attributed rather than just
   counted. */
static void clip_report(const Model* m, const double q[K], double W) {
    double I = integral(m), lost = 0.0, top[3] = { 0, 0, 0 };
    int idx[3] = { -1, -1, -1 }, n = 0;
    for (int i = 0; i < K; i++) {
        if (m->p[i] <= 0.0 || q[i] <= 0.0) continue;
        double w = m->p[i] / q[i];
        if (w > W) {
            double l = m->p[i] * m->f[i] * (1.0 - W / w);
            lost += l;
            n++;
            for (int s = 0; s < 3; s++) {
                if (l > top[s]) {
                    for (int t = 2; t > s; t--) { top[t] = top[t - 1]; idx[t] = idx[t - 1]; }
                    top[s] = l; idx[s] = i;
                    break;
                }
            }
        }
    }
    for (int s = 0; s < 3 && idx[s] >= 0; s++) {
        int i = idx[s];
        printf("  lost from i=%d p=%.6g f=%.6g w=%.6g share_of_lost=%.6g\n",
               i, m->p[i], m->f[i], m->p[i] / q[i], lost > 0 ? top[s] / lost : 0.0);
    }
    printf("  clipped=%d lost_total=%.12g lost_frac_of_I=%.6g\n",
           n, lost, I > 0.0 ? lost / I : 0.0);
}

/* ------------------------------------------------------------------ */
/* deterministic PRNG (splitmix64) and a small sampler, for the MC checks */

static unsigned long long rng_state;

static double next_u(void) {
    rng_state += 0x9E3779B97F4A7C15ULL;
    unsigned long long z = rng_state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return (double)((z ^ (z >> 31)) >> 11) * (1.0 / 9007199254740992.0);
}

static int sample_cell(const double q[K]) {
    double u = next_u(), acc = 0.0;
    for (int i = 0; i < K; i++) {
        acc += q[i];
        if (u < acc) return i;
    }
    return K - 1;
}

/* Monte Carlo the compensated estimator: sample from q, average X. */
static void monte_carlo(const Model* m, const double q[K], double W, long n,
                        double* out_mean, double* out_var, long* out_hits) {
    double s = 0.0, s2 = 0.0;
    long hits = 0;
    for (long k = 0; k < n; k++) {
        int i = sample_cell(q);
        double w = m->p[i] / q[i];
        if (W > 0.0 && w > W) w = W;
        double x = w * m->f[i];
        s += x;
        s2 += x * x;
        if (i == K - 1) hits++;
    }
    *out_mean = s / (double)n;
    *out_var = s2 / (double)n - (*out_mean) * (*out_mean);
    if (*out_var < 0.0) *out_var = 0.0;
    *out_hits = hits;
}

/* ------------------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static void check_close(const char* name, double got, double want, double tol) {
    double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
    int ok = fabs(got - want) <= tol * scale;
    printf("%s=%d got=%.12g want=%.12g\n", name, ok ? 1 : 0, got, want);
    if (!ok) failures++;
}

/* 63 dim cells (area 1/64 each share of the remainder) plus one bright cell
   with 1/1024 of the area and 16384x the radiance.  Natural sampling hits the
   bright cell once per 1024 samples, yet it carries 94% of I. */
static Model make_scene(void) {
    Model m;
    memset(&m, 0, sizeof(m));
    double pb = 1.0 / 1024.0;
    for (int i = 0; i < K - 1; i++) { m.p[i] = (1.0 - pb) / (K - 1); m.f[i] = 1.0; }
    m.p[K - 1] = pb;
    m.f[K - 1] = 16384.0;
    return m;
}

int main(void) {
    Model m = make_scene();
    double I = integral(&m);
    double q[K], g[K];
    int bright = K - 1;

    /* ---- 1. baseline: natural sampling, weight is 1 ------------------ */
    Moments nat = compensated(&m, m.p, 0.0);
    check("nat_unbiased", nat.unbiased);
    check_close("nat_mean", nat.mean, 16.9990234375, 1e-12);
    check_close("nat_var", nat.var, 261856.03222560883, 1e-9);
    check_close("nat_wmax", nat.wmax, 1.0, 1e-12);
    check_close("nat_p_bright", m.p[bright], 1.0 / 1024.0, 1e-15);
    check_close("nat_bright_energy_share",
                m.p[bright] * m.f[bright] / I, 0.941230539438157, 1e-9);
    /* the bright cell alone supplies almost all of E[X^2] */
    check_close("nat_bright_second_moment_share",
                second_moment_share(&m, m.p, bright), 0.999996189043, 1e-9);
    printf("nat_sd=%.12g sd_over_mean=%.6g\n",
           sqrt(nat.var), sqrt(nat.var) / I);

    /* ---- 2. energy proposal: unbiased, zero variance ----------------- */
    energy_proposal(&m, q);
    Moments en = compensated(&m, q, 0.0);
    check("energy_unbiased", en.unbiased);
    check("energy_zero_var", en.var < 1e-9);
    check_close("energy_q_bright", q[bright], 16.0 / I, 1e-12);
    check_close("energy_wmax", en.wmax, I, 1e-12);
    /* same proposal, weights dropped: hugely biased */
    Moments raw = uncompensated(&m, q);
    check("energy_raw_biased", !raw.unbiased);
    check_close("energy_raw_mean", raw.mean, 15421.179927615327, 1e-9);
    check_close("energy_raw_var", raw.var, 14846858.754091268, 1e-9);
    check("energy_comp_mse_better", en.mse < raw.mse);
    printf("energy_raw_rel_bias=%.6g\n", (raw.mean - I) / I);

    /* ---- 3. mixing natural -> energy: variance falls ---------------- */
    printf("lam,mean,var,var_ratio,drop_pct,wmax\n");
    double prev = nat.var, base = nat.var;
    for (int k = 0; k <= 4; k++) {
        double lam = 0.25 * k;
        mixed_proposal(&m, lam, q);
        Moments r = compensated(&m, q, 0.0);
        printf("mix lam=%.2f mean=%.12g var=%.12g ratio=%.9g drop=%.6f wmax=%.9g\n",
               lam, r.mean, r.var, r.var / base, 100.0 * (1.0 - r.var / base), r.wmax);
        check("mix_unbiased", r.unbiased);
        if (k > 0) check("mix_var_drops", r.var <= prev + 1e-9);
        prev = r.var;
    }
    mixed_proposal(&m, 0.25, q);
    check_close("mix_q1_var", compensated(&m, q, 0.0).var, 796.901294, 1e-7);
    mixed_proposal(&m, 0.50, q);
    check_close("mix_q2_var", compensated(&m, q, 0.0).var, 256.3251871, 1e-7);
    mixed_proposal(&m, 0.75, q);
    check_close("mix_q3_var", compensated(&m, q, 0.0).var, 76.95030967, 1e-7);

    /* ---- 4. closed-form variance agrees with the exact moments ----- */
    for (int k = 0; k <= 3; k++) {
        double lam = 0.25 * k;
        mixed_proposal(&m, lam, q);
        Moments r = compensated(&m, q, 0.0);
        check_close("var_formula_matches", var_closed_form(&m, q), r.var, 1e-9);
    }

    /* standard error at N = 1e6 samples */
    {
        mixed_proposal(&m, 0.75, q);
        double v75 = compensated(&m, q, 0.0).var;
        printf("se_at_1e6 nat=%.6g mix0.75=%.6g energy=%.6g\n",
               sqrt(nat.var / 1e6), sqrt(v75 / 1e6), 0.0);
    }

    /* ---- 5. Monte Carlo sanity -------------------------------------- */
    rng_state = 0;
    energy_proposal(&m, q);
    {   /* every sample equals I, so the mean is exact for any seed */
        double mm, vv; long hits;
        monte_carlo(&m, q, 0.0, MC_N, &mm, &vv, &hits);
        check_close("mc_energy_mean", mm, I, 1e-9);
        check_close("mc_energy_var", vv, 0.0, 1e-9);
        printf("mc_energy_hits=%ld/%ld\n", hits, MC_N);
    }
    {
        double mm, vv; long hits;
        mixed_proposal(&m, 0.75, q);
        Moments r = compensated(&m, q, 0.0);
        monte_carlo(&m, q, 0.0, MC_N, &mm, &vv, &hits);
        /* within 5 standard errors of I */
        check_close("mc_mix_mean_5se", mm, I, 5.0 * sqrt(r.var / (double)MC_N) / I);
        printf("mc_mix mean=%.12g se=%.6g var=%.12g (theory %.12g) hits=%ld\n",
               mm, sqrt(r.var / (double)MC_N), vv, r.var, hits);
    }

    /* ---- 6. bad pilot: the bright cell is under-sampled ------------ */
    /* Compare three estimators of the same I under the same mis-piloted
       proposal q: (a) compensated, (b) compensated with the weight dropped
       (i.e. a plain average of f under q), (c) direct sampling from p. */
    printf("pilot k,qbright,wbright,comp_var,raw_var,nat_var,"
           "comp_bias,raw_bias,comp_vs_nat,comp_vs_raw\n");
    for (int k = 1; k <= 30; k *= 2) {
        for (int i = 0; i < K; i++) g[i] = 1.0;
        g[bright] = 16384.0 / (double)(1 << k);   /* understate by 2^k */
        pilot_proposal(&m, g, q);
        Moments r = compensated(&m, q, 0.0);
        Moments u = uncompensated(&m, q);
        printf("pilot k=%d qbright=%.9g wbright=%.9g comp_var=%.12g raw_var=%.12g "
               "nat_var=%.12g comp_bias=%.3g raw_bias=%.6g comp_vs_nat=%.9g "
               "comp_vs_raw=%.9g\n",
               k, q[bright], r.wmax, r.var, u.var, nat.var,
               r.mean - I, u.mean - I, r.var / nat.var, r.var / u.var);
        /* the compensated estimator stays unbiased under every q ... */
        check("pilot_unbiased", r.unbiased);
        /* ... while the one without the weight does not */
        check("pilot_raw_biased", !u.unbiased);
        /* when the pilot is better than natural, compens-var beats nat-var */
        if (k <= 8) check("pilot_comp_var_below_nat", r.var < nat.var);
        /* Once the pilot under-samples the bright cell hard enough, dropping
           the weight gives a *lower* variance than compensating: that is the
           firefly blow-up.  The unweighted estimate is then simply wrong
           (bias = 12..14500), so we require the compensated estimator to
           win on MSE until the pilot is bad, and print the crossover. */
        if (k <= 8) check("pilot_comp_mse_below_raw", r.mse < u.mse);
        else printf("  note: at k=%d unweighted mse=%.6g < compensated mse=%.6g "
                    "but unweighted bias=%.3g (biased estimator, not comparable)\n",
                    k, u.mse, r.mse, u.mean - I);
    }

    /* m = 30: the weight is ~65472 and the variance is 1.7e10. */
    for (int i = 0; i < K; i++) g[i] = 1.0;
    g[bright] = 16384.0 / 1073741824.0;
    pilot_proposal(&m, g, q);
    Moments bad = compensated(&m, q, 0.0);
    check_close("bad_wmax", bad.wmax, 65472.0009765625, 1e-9);
    check_close("bad_var", bad.var, 1.716309194e10, 1e-8);
    check("bad_unbiased", bad.unbiased);
    /* 200k samples expect 200000*q_bright < 1 hit: unbiased, never seen */
    printf("bad_expected_hits=%ld qbright=%.9g\n",
           (long)(MC_N * q[bright]), q[bright]);
    check("bad_expected_hits_under_one", MC_N * q[bright] < 1.0);

    /* ---- 7. truncate the unbounded weight and account the bias ------ */
    /* Bias identity: I - E = lost, exactly.  On the energy proposal the
       bright cell has the *smallest* weight (I/f = 0.001), so clamping at 16
       clips the 63 dim cells instead; still a useful sanity check of the
       identity and of reporting. */
    energy_proposal(&m, q);
    Moments e16 = compensated(&m, q, 16.0);
    check("clip16_biased", !e16.unbiased);
    check_close("clip16_mean", e16.mean, 16.9403113689895, 1e-9);
    check_close("clip16_lost", e16.lost, 0.0587120685105, 1e-9);
    check_close("clip16_bias_equals_lost", I - e16.mean, e16.lost, 1e-12);
    check("clip16_clips_dim", e16.nclip == K - 1 && e16.worst != bright);
    Moments e17 = compensated(&m, q, 17.0);
    check("clip17_exact", e17.unbiased && e17.nclip == 0);

    /* On the mis-piloted proposal the bright weight reaches 65472.  Clamping
       is the only way to keep the variance finite; the bias is exactly the
       energy of the clipped cell and is reported per event. */
    for (int i = 0; i < K; i++) g[i] = 1.0;
    g[bright] = 16384.0 / 1073741824.0;
    pilot_proposal(&m, g, q);
    printf("clamp on bad pilot (unclamped wmax=%.9g):\n", bad.wmax);
    for (int k = 0; k <= 10; k += 2) {
        double W = pow(2.0, k);
        Moments r = compensated(&m, q, W);
        printf("  W=%-8g mean=%.12g bias=%.12g lost=%.12g identity=%.3g "
               "var=%.12g mse=%.12g nclip=%d worst=%d\n",
               W, r.mean, r.mean - I, r.lost, (r.mean - I) + r.lost,
               r.var, r.mse, r.nclip, r.worst);
        check("clip_bias_identity",
              fabs((r.mean - I) + r.lost) <= 1e-9 * (fabs(I) + 1.0));
        check("clip_source_is_bright", r.nclip == 1 && r.worst == bright);
    }
    printf("clip detail at W=1:\n");
    clip_report(&m, q, 1.0);
    Moments c1 = compensated(&m, q, 1.0);
    check_close("clip1_lost_frac", c1.lost / I, 0.9412161633609729, 1e-12);
    check_close("clip1_var", c1.var, 4.003421683, 1e-7);

    /* At W=1 the clamped estimator has 94% of I stripped out, so it is
       useless as an estimate even though its variance is tiny.  MSE alone
       would still prefer it to the natural sampler: MSE must not be used to
       bless a biased estimator. */
    check("clip1_mse_below_natural", c1.mse < nat.mse);
    printf("clip1_bias_over_I=%.6g mse=%.12g natural_mse=%.12g\n",
           (c1.mean - I) / I, c1.mse, nat.mse);

    printf("I=%.12g natural_sd=%.12g\n", I, sqrt(nat.var));
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
