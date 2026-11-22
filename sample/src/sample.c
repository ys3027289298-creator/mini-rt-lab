
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

/* ---------- correlated multi-jittered sampling (stacked tiles) ---------- */

#define CMJ_MAX_GRID 64

typedef struct { unsigned s; } Rng;
static unsigned rng_u32(Rng* r) {
    unsigned x = r->s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    r->s = x;
    return x;
}
static double rng_u01(Rng* r) { return (double)(rng_u32(r) >> 8) * (1.0 / 16777216.0); }

static int cmj_perm[CMJ_MAX_GRID * CMJ_MAX_GRID];
static double cmj_buf[2 * CMJ_MAX_GRID * CMJ_MAX_GRID];

/* One layer = grid*grid stratified 2D points, one per tile of a grid x grid
   tiling. Consecutive points are grouped into blocks of `block`; every point
   in a block shares the block's jitter vector (intra-block correlation), each
   block draws its own jitter (inter-block independence), and each `layer`
   reseeds the RNG (inter-layer independence). A point never leaves its tile
   for any block size, so the stratification bounds hold for all block >= 1.
   Returns 0 on success, -1 when block < 1 or block > grid*grid. */
static int cmj_layer(double* out, int grid, int block, unsigned layer) {
    if (grid < 1 || grid > CMJ_MAX_GRID) return -1;
    int n = grid * grid;
    if (block < 1 || block > n) return -1;

    Rng r = { 0x9e3779b9u + layer * 0x85ebca6bu };
    if (r.s == 0) r.s = 0x27d4eb2fu;

    for (int i = 0; i < n; i++) cmj_perm[i] = i;
    for (int i = n - 1; i > 0; i--) {
        int j = (int)(rng_u32(&r) % (unsigned)(i + 1));
        int t = cmj_perm[i]; cmj_perm[i] = cmj_perm[j]; cmj_perm[j] = t;
    }

    double jx = 0.0, jy = 0.0;
    for (int i = 0; i < n; i++) {
        if (i % block == 0) { jx = rng_u01(&r); jy = rng_u01(&r); }
        int sx = cmj_perm[i] % grid, sy = cmj_perm[i] / grid;
        out[2 * i]     = ((double)sx + jx) / (double)grid;
        out[2 * i + 1] = ((double)sy + jy) / (double)grid;
    }
    return 0;
}

/* max |count - n/(k*k)| / n over a k x k histogram of 2D points */
static double hist2d_err(const double* xy, int n, int k) {
    static int cnt[32 * 32];
    memset(cnt, 0, sizeof(int) * (size_t)(k * k));
    for (int i = 0; i < n; i++) {
        int bx = (int)(xy[2 * i] * k);     if (bx >= k) bx = k - 1;
        int by = (int)(xy[2 * i + 1] * k); if (by >= k) by = k - 1;
        cnt[by * k + bx]++;
    }
    double expect = (double)n / (double)(k * k);
    double err = 0.0;
    for (int i = 0; i < k * k; i++) {
        double d = fabs(cnt[i] - expect) / (double)n;
        if (d > err) err = d;
    }
    return err;
}

/* Deterministic upper bound on hist2d_err for any one-point-per-tile pattern
   on a grid x grid tiling: a histogram bin touches at most (ceil(g/k)+1)^2
   tiles and fully contains at least max(0, floor(g/k)-1)^2. If k divides g
   the bins align with tiles and the error is exactly 0. */
static double hist2d_bound(int grid, int k) {
    if (grid % k == 0) return 0.0;
    int n = grid * grid;
    double expect = (double)n / (double)(k * k);
    int rows_in = grid / k - 1; if (rows_in < 0) rows_in = 0;
    int rows_touch = (grid + k - 1) / k + 1;
    double hi = (double)(rows_touch * rows_touch) - expect;
    double lo = expect - (double)(rows_in * rows_in);
    double b = hi > lo ? hi : lo;
    return b / (double)n;
}

/* Pearson r of x-coordinates of consecutive point pairs (i, i+1),
   accumulated over `layers` independent layers. */
static double cmj_consecutive_corr(int grid, int block, int layers) {
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    long m = 0;
    int n = grid * grid;
    for (int L = 0; L < layers; L++) {
        cmj_layer(cmj_buf, grid, block, (unsigned)L);
        for (int i = 0; i + 1 < n; i++) {
            double a = cmj_buf[2 * i], b = cmj_buf[2 * (i + 1)];
            sx += a; sy += b; sxx += a * a; syy += b * b; sxy += a * b; m++;
        }
    }
    double cov = sxy / m - (sx / m) * (sy / m);
    double vx = sxx / m - (sx / m) * (sx / m);
    double vy = syy / m - (sy / m) * (sy / m);
    return cov / sqrt(vx * vy);
}

/* Pearson r between the sequences of layer 0 and layer 1 (x and y pooled). */
static double cmj_cross_layer_corr(int grid, int block) {
    static double other[2 * CMJ_MAX_GRID * CMJ_MAX_GRID];
    cmj_layer(cmj_buf, grid, block, 0);
    cmj_layer(other, grid, block, 1);
    int n = grid * grid;
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
    for (int i = 0; i < 2 * n; i++) {
        double a = cmj_buf[i], b = other[i];
        sx += a; sy += b; sxx += a * a; syy += b * b; sxy += a * b;
    }
    double m = 2.0 * n;
    double cov = sxy / m - (sx / m) * (sy / m);
    double vx = sxx / m - (sx / m) * (sx / m);
    double vy = syy / m - (sy / m) * (sy / m);
    return cov / sqrt(vx * vy);
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

    /* --- correlated multi-jittered sampling --- */

    /* block size zero (and negative / oversized) must be rejected */
    check("cmj_rejects_zero_block", cmj_layer(cmj_buf, 8, 0, 0) == -1);
    check("cmj_rejects_negative_block", cmj_layer(cmj_buf, 8, -4, 0) == -1);
    check("cmj_rejects_oversize_block", cmj_layer(cmj_buf, 8, 65, 0) == -1);
    check("cmj_accepts_valid_block", cmj_layer(cmj_buf, 8, 8, 0) == 0);

    /* 2D histogram uniformity vs pure random, with deterministic bound */
    int n = CMJ_MAX_GRID * CMJ_MAX_GRID;
    cmj_layer(cmj_buf, CMJ_MAX_GRID, n, 0);          /* worst case: block = n */
    double e8 = hist2d_err(cmj_buf, n, 8);
    double b8 = hist2d_bound(CMJ_MAX_GRID, 8);       /* 8 | 64 -> bound 0 */
    check("cmj_hist8_within_bound", e8 <= b8 + 1e-12);

    cmj_layer(cmj_buf, CMJ_MAX_GRID, 7, 1);          /* arbitrary block size */
    double e12 = hist2d_err(cmj_buf, n, 12);
    double b12 = hist2d_bound(CMJ_MAX_GRID, 12);
    check("cmj_hist12_within_bound", e12 <= b12 + 1e-12);

    Rng rr = { 0xdeadbeefu };
    for (int i = 0; i < 2 * n; i++) cmj_buf[i] = rng_u01(&rr);
    double r8 = hist2d_err(cmj_buf, n, 8);
    double r12 = hist2d_err(cmj_buf, n, 12);
    check("cmj_beats_random_hist", e8 <= r8 && e12 <= r12);
    printf("hist8_err=%.6f bound=%.6f random=%.6f\n", e8, b8, r8);
    printf("hist12_err=%.6f bound=%.6f random=%.6f\n", e12, b12, r12);

    /* correlation vs block size. Theory: with the per-layer stratum
       permutation, consecutive positions are sampled without replacement
       (antithetic term -1/(n-1)) and a shared block jitter adds +1/g^2 for
       pairs inside one block, so r(B) = (1-1/B)/g^2 - 1/(n-1). */
    printf("# consecutive-pair x correlation, grid=8, 20000 layers\n");
    double c_full = cmj_consecutive_corr(8, 64, 20000);
    double c_16   = cmj_consecutive_corr(8, 16, 20000);
    double c_4    = cmj_consecutive_corr(8, 4, 20000);
    double c_1    = cmj_consecutive_corr(8, 1, 20000);
    printf("corr_block_64=%.5f (expect %.5f)\n", c_full, 63.0 / 64 / 64 - 1.0 / 63);
    printf("corr_block_16=%.5f (expect %.5f)\n", c_16,  15.0 / 16 / 64 - 1.0 / 63);
    printf("corr_block_4=%.5f (expect %.5f)\n",  c_4,   3.0 / 4 / 64 - 1.0 / 63);
    printf("corr_block_1=%.5f (expect %.5f)\n",  c_1,   -1.0 / 63);
    check("cmj_intra_layer_correlated", c_full > c_1 + 0.010);
    check("cmj_corr_drops_with_block", c_full > c_4 && c_4 > c_1);
    check("cmj_block1_antithetic", fabs(c_1 - (-1.0 / 63)) < 0.004);

    /* layers are independent of each other */
    double cl = cmj_cross_layer_corr(CMJ_MAX_GRID, n);
    printf("corr_cross_layer=%.5f\n", cl);
    check("cmj_inter_layer_independent", fabs(cl) < 0.05);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
