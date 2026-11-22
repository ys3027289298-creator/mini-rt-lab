/* firefly: path-weight truncation for isolated bright pixels ("fireflies"),
   with a certified bias bound and an untruncated reference.

   Scene model (synthetic, deterministic):
       X_pix = L(x,y) + s(x) * W
   where L is a smooth low-albedo base, s a per-column scale, and W is the
   path weight.  With probability P_TRUE a sample carries a heavy Pareto
   tail; otherwise W = 0.  The tail is what produces isolated bright pixels.

   Truncation clamps the *weight*, not the pixel:  W_K = min(W, K).
   The per-pixel loss is therefore exactly s * E[(W-K)^+], and since the
   lost energy is non-negative the image bias is one-sided (darker).

   Certified bound.  We do not trust P_TRUE / A_TRUE; we use a declared
   envelope that must majorize the true survival function:
       P(W > w) <= P_ENV * (W_MIN / w)^A_ENV      for w >= W_MIN,
       P(W > 0) <= P_ENV.
   Since (W-K)^+ = int_K^inf 1{W>t} dt, tail dominance gives
       E[(W-K)^+] <= Delta_env(K)
   for every K >= 0, with Delta_env in closed form (delta_tail below).
   With s <= S_MAX this yields the certified per-pixel bound
       |bias| <= S_MAX * Delta_env(K).
   The renderer is validated against it empirically: the measured bias with
   a 99.7% (3-sigma) allowance must stay under the bound, for every K.

   The reference image is rendered from the same random streams (common
   random numbers), so the truncation bias is isolated from sampling noise. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define IMG_W 96
#define IMG_H 96
#define IMG_N (IMG_W * IMG_H)

#define SPP 16
#define NSEEDS 24

/* scene */
#define L0 0.18
#define L1 0.40
#define L_RIPPLE 0.03
#define S_MID 0.90
#define S_AMP 0.10
#define S_MAX (S_MID + S_AMP)

/* true weight generator */
#define P_TRUE 0.010
#define A_TRUE 2.5
#define W_MIN 4.0

/* declared conservative envelope (must majorize the true tail) */
#define P_ENV 0.012
#define A_ENV 2.4

/* isolated-bright detector: value > local_median + DET_T * local_sigma,
   where sigma is the neighbour MAD scaled to a standard deviation, and the
   value must also clear DET_FLOOR.  Thresholds are always evaluated against
   the *untruncated* reference so counts across K are comparable. */
#define DET_T 4.0
#define DET_FLOOR 1.0
#define SIG_FLOOR 0.02

/* acceptance: certified relative bound budget, and energy-loss report level */
#define TAU 0.050
#define EPS 0.035

static const double KS[] = { 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 32 };
#define NK ((int)(sizeof(KS) / sizeof(KS[0])))

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* --- RNG: splitmix64, independent stream per (seed, sample) ------------- */

static unsigned long long mix64(unsigned long long z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
typedef struct { unsigned long long s; } Rng;
static void rng_seed(Rng* r, unsigned long long seed) { r->s = mix64(seed); }
static double rng_next(Rng* r) {
    r->s = mix64(r->s);
    return (double)(r->s >> 11) * (1.0 / 9007199254740992.0); /* [0,1) */
}

/* --- scene fields ------------------------------------------------------- */

#define TWO_PI 6.28318530717958647692

static double field_base(int x, int y) {
    return L0 + (L1 - L0) * (double)x / (double)(IMG_W - 1)
              + L_RIPPLE * sin(TWO_PI * (double)y / (double)IMG_H);
}
static double field_scale(int x) {
    return S_MID + S_AMP * cos(TWO_PI * (double)x / (double)IMG_W);
}

/* exact mean of the untruncated image over pixels (no sampling) */
static double reference_mean_exact(void) {
    double sum = 0.0;
    double ew = P_TRUE * W_MIN * A_TRUE / (A_TRUE - 1.0);
    for (int y = 0; y < IMG_H; y++)
        for (int x = 0; x < IMG_W; x++)
            sum += field_base(x, y) + field_scale(x) * ew;
    return sum / (double)IMG_N;
}

/* --- weight tail and its closed-form loss ------------------------------- */

/* E[(W-K)^+] for P(W=0)=1-p and Pareto(xmin, a) otherwise */
static double delta_tail(double K, double p, double a, double xmin) {
    if (K <= 0.0) return p * xmin * a / (a - 1.0);
    if (K < xmin) return p * (xmin * a / (a - 1.0) - K);
    return p * pow(xmin, a) * pow(K, 1.0 - a) / (a - 1.0);
}
static double delta_true(double K) { return delta_tail(K, P_TRUE, A_TRUE, W_MIN); }
static double delta_env(double K)  { return delta_tail(K, P_ENV,  A_ENV,  W_MIN); }

/* certified per-pixel absolute bias bound */
static double bias_bound_abs(double K) { return S_MAX * delta_env(K); }

/* --- rendering (common random numbers) ---------------------------------- */

static double img_ref[IMG_N];
static double img_k[NK][IMG_N];

static void render_all(unsigned long long seed) {
    memset(img_ref, 0, sizeof(img_ref));
    memset(img_k, 0, sizeof(img_k));

    for (int s = 0; s < SPP; s++) {
        Rng r;
        rng_seed(&r, seed * 1000003ULL + (unsigned long long)s);
        for (int y = 0; y < IMG_H; y++) {
            for (int x = 0; x < IMG_W; x++) {
                double w = 0.0;
                if (rng_next(&r) < P_TRUE) {
                    double u = rng_next(&r);
                    if (u < 1e-300) u = 1e-300;
                    w = W_MIN * pow(u, -1.0 / A_TRUE);
                }
                double v = field_base(x, y);
                double sc = field_scale(x);
                int i = y * IMG_W + x;
                img_ref[i] += v + sc * w;
                for (int k = 0; k < NK; k++) {
                    double wc = (w > KS[k]) ? KS[k] : w;
                    img_k[k][i] += v + sc * wc;
                }
            }
        }
    }
    for (int i = 0; i < IMG_N; i++) {
        img_ref[i] /= (double)SPP;
        for (int k = 0; k < NK; k++) img_k[k][i] /= (double)SPP;
    }
}

/* --- isolated-bright detector ------------------------------------------- */

static int cmp_double(const void* a, const void* b) {
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

/* median and MAD-sigma of the 8 neighbours (3x3 minus centre), edge-clamped */
static void local_stats(const double* img, int x, int y, double* med, double* sig) {
    double nb[8], dev[8];
    int n = 0;
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (dx == 0 && dy == 0) continue;
            int xx = x + dx, yy = y + dy;
            if (xx < 0) xx = 0;
            if (xx >= IMG_W) xx = IMG_W - 1;
            if (yy < 0) yy = 0;
            if (yy >= IMG_H) yy = IMG_H - 1;
            nb[n++] = img[yy * IMG_W + xx];
        }
    }
    qsort(nb, 8, sizeof(double), cmp_double);
    *med = 0.5 * (nb[3] + nb[4]);
    for (int i = 0; i < 8; i++) dev[i] = fabs(nb[i] - *med);
    qsort(dev, 8, sizeof(double), cmp_double);
    *sig = 0.5 * (dev[3] + dev[4]) / 0.6745;
    if (*sig < SIG_FLOOR) *sig = SIG_FLOOR;
}

static int count_bright(const double* img) {
    int c = 0;
    for (int y = 0; y < IMG_H; y++) {
        for (int x = 0; x < IMG_W; x++) {
            double med, sig;
            local_stats(img_ref, x, y, &med, &sig);
            double v = img[y * IMG_W + x];
            if (v > med + DET_T * sig && v > DET_FLOOR) c++;
        }
    }
    return c;
}

static double image_mean(const double* img) {
    double s = 0.0;
    for (int i = 0; i < IMG_N; i++) s += img[i];
    return s / (double)IMG_N;
}

int main(void) {
    /* 1. the declared envelope really dominates the true tail */
    int env_ok = 1;
    for (int i = 0; i <= 20000; i++) {
        double w = W_MIN * pow(10.0, 4.0 * (double)i / 20000.0);
        if (P_TRUE * pow(W_MIN / w, A_TRUE) > P_ENV * pow(W_MIN / w, A_ENV) + 1e-18)
            env_ok = 0;
    }
    check("envelope_dominates_true_tail", env_ok);

    /* 2. the envelope bound really covers the generator's exact bias */
    int covers_ok = 1;
    for (int k = 0; k < NK; k++) {
        double exact = (S_MID) * delta_true(KS[k]); /* mean scale is S_MID */
        if (exact > bias_bound_abs(KS[k]) + 1e-15) covers_ok = 0;
    }
    check("bound_covers_exact_bias", covers_ok);

    const unsigned long long seed0 = 9001ULL;
    double mu_ref = reference_mean_exact();

    int cnt_ref[NSEEDS];
    int cnt_k[NK][NSEEDS];
    double bias[NK][NSEEDS];   /* mean_K - mean_ref, common random numbers */
    double loss[NK][NSEEDS];   /* measured relative energy loss vs reference */

    for (int s = 0; s < NSEEDS; s++) {
        render_all(seed0 + 7919ULL * (unsigned long long)s);
        double mr = image_mean(img_ref);
        cnt_ref[s] = count_bright(img_ref);
        for (int k = 0; k < NK; k++) {
            double mk = image_mean(img_k[k]);
            cnt_k[k][s] = count_bright(img_k[k]);
            bias[k][s] = mk - mr;
            loss[k][s] = (mr - mk) / mr;
        }
    }

    /* clamping is monotone in K, so with a detector fixed on the reference
       the bright count must be non-decreasing in K and never exceed ref */
    int mono_ok = 1;
    for (int s = 0; s < NSEEDS; s++) {
        if (cnt_k[0][s] > cnt_ref[s]) mono_ok = 0;
        if (cnt_k[NK - 1][s] > cnt_ref[s]) mono_ok = 0;
        for (int k = 1; k < NK; k++)
            if (cnt_k[k][s] < cnt_k[k - 1][s]) mono_ok = 0;
    }
    check("bright_counts_monotone_in_K", mono_ok);

    double ref_mean_count = 0.0;
    for (int s = 0; s < NSEEDS; s++) ref_mean_count += (double)cnt_ref[s];
    ref_mean_count /= (double)NSEEDS;

    printf("mu_ref_exact=%.6f\n", mu_ref);
    printf("ref_bright=%.2f\n", ref_mean_count);
    printf("K bright drop_pct bound_abs bound_pct meas_abs meas_pct se3_pct "
           "loss_pct flags\n");

    int bias_ok = 1;      /* measured bias stays under the certified bound */
    int loss_flag_any = 0;
    double best_K = -1.0;
    double best_bound_pct = 0.0, best_meas_pct = 0.0, best_loss_pct = 0.0;
    double best_count = 0.0;

    for (int k = 0; k < NK; k++) {
        double K = KS[k];
        double sum_c = 0.0, sum_b = 0.0;
        for (int s = 0; s < NSEEDS; s++) {
            sum_c += (double)cnt_k[k][s];
            sum_b += bias[k][s];
        }
        double mean_c = sum_c / (double)NSEEDS;
        double mean_b = sum_b / (double)NSEEDS;
        double sq = 0.0;
        for (int s = 0; s < NSEEDS; s++) {
            double d = bias[k][s] - mean_b;
            sq += d * d;
        }
        double se = sqrt(sq / (double)(NSEEDS - 1) / (double)NSEEDS);
        double sum_l = 0.0;
        for (int s = 0; s < NSEEDS; s++) sum_l += loss[k][s];
        double mean_l = sum_l / (double)NSEEDS;

        double bound_abs = bias_bound_abs(K);
        double bound_pct = 100.0 * bound_abs / mu_ref;
        double meas_pct = 100.0 * mean_b / mu_ref;
        double se3_pct = 100.0 * 3.0 * se / mu_ref;
        double loss_pct = 100.0 * mean_l;
        double drop_pct = 100.0 * (1.0 - mean_c / ref_mean_count);

        int under_bound = fabs(mean_b) <= bound_abs + 3.0 * se;
        if (!under_bound) bias_ok = 0;

        int loss_flag = mean_l > EPS;
        if (loss_flag) loss_flag_any = 1;

        if (best_K < 0.0 && (bound_abs / mu_ref) <= TAU && !loss_flag) {
            best_K = K;
            best_bound_pct = bound_pct;
            best_meas_pct = meas_pct;
            best_loss_pct = loss_pct;
            best_count = mean_c;
        }

        printf("%.0f %.2f %.1f %.5f %.2f %.5f %.2f %.2f %.2f %s\n",
               K, mean_c, drop_pct, bound_abs, bound_pct, mean_b, meas_pct,
               se3_pct, loss_pct,
               !under_bound ? "BOUND_VIOLATION" : (loss_flag ? "LOSS" : "ok"));
    }

    check("measured_bias_within_bound", bias_ok);
    check("small_threshold_loss_reported", loss_flag_any);
    check("budget_threshold_selected", best_K >= 0.0);
    if (best_K >= 0.0) {
        printf("selected_K=%.0f bound_pct=%.2f meas_pct=%.2f loss_pct=%.2f "
               "bright=%.2f drop_pct=%.1f\n",
               best_K, best_bound_pct, best_meas_pct, best_loss_pct,
               best_count, 100.0 * (1.0 - best_count / ref_mean_count));
        check("selected_threshold_removes_majority",
              best_count <= 0.5 * ref_mean_count);
        check("selected_loss_within_eps", best_loss_pct <= 100.0 * EPS);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
