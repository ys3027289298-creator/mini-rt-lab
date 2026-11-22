/* dither: 8-bit quantization error analysis — ordered dither vs error diffusion.
 *
 * Theory (q = 1 LSB):
 *   round:    |e| <= q/2, error is a deterministic sawtooth of the signal
 *             -> low-frequency flat runs on smooth gradients (banding).
 *   ordered:  y = Q(x+t), t in [-q/2, q/2)  ->  |e| <= q.
 *             error two-valued {-f, 1-f} given frac f: zero mean,
 *             var = f(1-f), average q^2/6; white floor + Bayer tones.
 *   diffusion: discarded fraction e_i = u_i - Q(u_i) is pushed forward,
 *             u_i = x_i + sum w_k e_{i-k}, |e_i| <= q/2,
 *             y_i - x_i = sum w_k e_{i-k} - e_i  ->  |e| <= (q/2)(1 + sum|w_k|).
 *             Conservation: sum(y-x) = (sum w_k - 1) sum e_i - buffered,
 *             so sum w_k = 1 keeps the running error bounded. Noise shaping:
 *             Y(z) = X(z) + E(z)(W(z) - 1); sum w_k = 1 nulls the DC noise
 *             gain (blue noise). sum w_k != 1 leaks DC and the accumulated
 *             error drifts -> reported by kernel_check.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define W 256
#define H 256
#define NK 4

static const double PI = 3.14159265358979323846;

static double ref[H][W]; /* high-precision reference image */
static int out[H][W];    /* quantized 8-bit output */

static const unsigned char BAYER[8][8] = {
    { 0,32, 8,40, 2,34,10,42},
    {48,16,56,24,50,18,58,26},
    {12,44, 4,36,14,46, 6,38},
    {60,28,52,20,62,30,54,22},
    { 3,35,11,43, 1,33, 9,41},
    {51,19,59,27,49,17,57,25},
    {15,47, 7,39,13,45, 5,37},
    {63,31,55,23,61,29,53,21}
};

typedef struct { int dx, dy; double w; } KernelTap;

/* Floyd–Steinberg */
static const KernelTap FS[NK] = {
    { 1, 0, 7.0/16 }, { -1, 1, 3.0/16 }, { 0, 1, 5.0/16 }, { 1, 1, 1.0/16 }
};

static double kernel_sum(const KernelTap* k) {
    double s = 0;
    for (int i = 0; i < NK; i++) s += k[i].w;
    return s;
}

static double kernel_sumabs(const KernelTap* k) {
    double s = 0;
    for (int i = 0; i < NK; i++) s += fabs(k[i].w);
    return s;
}

/* report kernel weight sum; diffusion drifts at DC unless it is exactly 1 */
static int kernel_check(const char* name, const KernelTap* k) {
    double s = kernel_sum(k);
    int ok = fabs(s - 1.0) < 1e-9;
    printf("kernel_%s_sum=%.9f\n", name, s);
    printf("kernel_%s_sum_ok=%d\n", name, ok);
    if (!ok)
        printf("warning: kernel %s weights sum to %.9f != 1 (DC leak, error drift)\n", name, s);
    return ok;
}

static void quant_round(void) {
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            out[y][x] = (int)lround(ref[y][x]);
}

static void quant_ordered(void) {
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double t = (BAYER[y & 7][x & 7] + 0.5) / 64.0 - 0.5;
            out[y][x] = (int)lround(ref[y][x] + t);
        }
}

static void quant_ed(const KernelTap* k) {
    static double err[H][W]; /* error pushed forward to each pixel */
    memset(err, 0, sizeof err);
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double u = ref[y][x] + err[y][x];
            int v = (int)lround(u);
            out[y][x] = v;
            double e = u - (double)v; /* discarded fraction, pushed forward */
            for (int i = 0; i < NK; i++) {
                int nx = x + k[i].dx, ny = y + k[i].dy;
                if (nx >= 0 && nx < W && ny < H) err[ny][nx] += k[i].w * e;
            }
        }
}

typedef struct {
    double maxe, mae, rmse, mean, var, lowfrac;
    long band;
} Stats;

static Stats analyze(void) {
    Stats s; memset(&s, 0, sizeof s);
    double sum = 0, sum2 = 0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double e = out[y][x] - ref[y][x];
            double ae = fabs(e);
            if (ae > s.maxe) s.maxe = ae;
            s.mae += ae; sum += e; sum2 += e * e;
        }
    long n = (long)W * H;
    s.mean = sum / n;
    s.mae /= n;
    s.rmse = sqrt(sum2 / n);
    s.var = sum2 / n - s.mean * s.mean;

    /* band pixels: members of horizontal flat runs too long to be dither
     * micro-structure (Bayer period is 8; threshold 16 gives margin) */
    for (int y = 0; y < H; y++) {
        int x = 0;
        while (x < W) {
            int x2 = x;
            while (x2 + 1 < W && out[y][x2 + 1] == out[y][x]) x2++;
            if (x2 - x + 1 >= 16) s.band += x2 - x + 1;
            x = x2 + 1;
        }
    }

    /* error spectrum: per-row DFT power, averaged; lowfreq = bins 1..W/16 */
    static double power[W / 2 + 1];
    memset(power, 0, sizeof power);
    for (int y = 0; y < H; y++) {
        static double e[W];
        double m = 0;
        for (int x = 0; x < W; x++) { e[x] = out[y][x] - ref[y][x]; m += e[x]; }
        m /= W;
        for (int x = 0; x < W; x++) e[x] -= m;
        for (int k = 1; k <= W / 2; k++) {
            double re = 0, im = 0;
            for (int x = 0; x < W; x++) {
                double a = -2.0 * PI * k * x / W;
                re += e[x] * cos(a); im += e[x] * sin(a);
            }
            power[k] += re * re + im * im;
        }
    }
    double lo = 0, tot = 0;
    for (int k = 1; k <= W / 2; k++) { tot += power[k]; if (k <= W / 16) lo += power[k]; }
    s.lowfrac = tot > 0 ? lo / tot : 0;
    return s;
}

static void report(const char* name, Stats s) {
    printf("max_%s=%.6f\n", name, s.maxe);
    printf("rmse_%s=%.6f\n", name, s.rmse);
    printf("mean_%s=%.6e\n", name, s.mean);
    printf("var_%s=%.6f\n", name, s.var);
    printf("band_%s=%ld\n", name, s.band);
    printf("lowfreq_%s=%.4f\n", name, s.lowfrac);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* smooth gradient: 8 levels over the width -> 32px flat bands when rounded */
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
            ref[y][x] = 100.0 + 6.0 * x / (W - 1) + 2.0 * y / (H - 1);

    /* kernel weight-sum report (sum != 1 => DC leak + drift) */
    check("kernel_fs_sum_ok", kernel_check("fs", FS));
    static const KernelTap BAD[NK] = {
        { 1, 0, 0.40 }, { -1, 1, 0.30 }, { 0, 1, 0.20 }, { 1, 1, 0.05 }
    };
    check("kernel_bad_detected", !kernel_check("bad", BAD));

    const double bound_round = 0.5;                 /* q/2 */
    const double bound_ord = 1.0;                   /* q/2 + q/2 */
    const double bound_ed = 0.5 * (1.0 + kernel_sumabs(FS)); /* (q/2)(1+sum|w|) */
    printf("bound_round=%.3f\nbound_ordered=%.3f\nbound_ed=%.3f\n",
           bound_round, bound_ord, bound_ed);

    quant_round();   Stats r = analyze(); report("round", r);
    quant_ordered(); Stats o = analyze(); report("ordered", o);
    quant_ed(FS);    Stats d = analyze(); report("ed", d);

    /* measured error must not exceed the theoretical bounds */
    check("bound_round_ok", r.maxe <= bound_round + 1e-9);
    check("bound_ordered_ok", o.maxe <= bound_ord + 1e-9);
    check("bound_ed_ok", d.maxe <= bound_ed + 1e-9);

    /* banding pixel counts and reduction */
    printf("band_reduction_ordered_pct=%.2f\n", 100.0 * (1.0 - (double)o.band / r.band));
    printf("band_reduction_ed_pct=%.2f\n", 100.0 * (1.0 - (double)d.band / r.band));
    check("band_ordered_less", o.band < r.band);
    check("band_ed_less", d.band < r.band);

    /* spectrum: diffusion/ordered push error energy out of the low band */
    check("lowfreq_ordered_less", o.lowfrac < r.lowfrac);
    check("lowfreq_ed_less", d.lowfrac < r.lowfrac);

    /* error variance: ramp sawtooth is uniform -> q^2/12; ordered-dither
     * error vs source is two-valued {-f, 1-f} -> average f(1-f) = q^2/6 */
    check("var_round_near_1_12", fabs(r.var - 1.0 / 12.0) < 0.01);
    check("var_ordered_near_1_6", fabs(o.var - 1.0 / 6.0) < 0.01);

    /* sum w = 1 keeps the running error sum bounded -> near-zero DC error */
    check("dc_ed_small", fabs(d.mean) < 0.01);

    /* a kernel with sum != 1 drifts: demonstrate with the bad kernel */
    quant_ed(BAD); Stats b = analyze();
    printf("drift_bad_mean=%.6e\n", b.mean);
    check("drift_when_sum_ne_1", fabs(b.mean) > 10.0 * fabs(d.mean));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
