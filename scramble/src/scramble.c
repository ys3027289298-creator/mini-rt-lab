/* scramble: randomized per-pixel sample sequences for reuse-free rendering.

   Reusing one sample sequence for every pixel and every frame leaves a
   fixed error pattern: it is correlated across neighbouring pixels and it
   never averages out over frames, so the picture shows regular structure
   ("pattern").  The fix is to randomize the sequence per pixel and per
   frame while keeping within-pixel stratification.

   Two randomization kernels are provided and can be combined:

     SC_ROT   Cranley-Patterson rotation: add a uniform random offset
              (mod 1) to every sample of the sequence.  Cheap, decorrelates
              sequences between pixels and frames, and leaves the sequence
              uniformly distributed inside the pixel.
     SC_PERM  Owen (nested uniform) digit scrambling: permute the binary
              digits of the sample index.  A bijection per dimension, so
              the (0,4,2)-net is preserved exactly and every frame is still
              perfectly stratified.

   Both kernels draw their randomness from a hash of
   (seed, frame?, pixel?, dimension, kernel), so any (seed, frame, pixel)
   reproduces the same stream on any run and on any machine: only integer
   hashing and IEEE-754 conversions are used.  SC_FRAME / SC_PIXEL control
   which of frame and pixel enter the hash, which allows the old
   "one sequence for everything" behaviour to be reproduced for
   comparison.

   The companion report (`./scramble bench`) renders a small analytic
   scene with each mode and reports:

     - sequence correlation between neighbouring pixels and between frames;
     - RMSE per frame and of the temporal mean;
     - persistence = rmse(mean) / rmse(frame);
     - pattern pixels: area whose error in the temporal mean (and in 4x4
       block means) exceeds 1% of the scene dynamic range, i.e. coherent
       error that survives averaging instead of per-pixel noise.
*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define SC_BITS 4u
#define SC_N (1u << SC_BITS)
#define SC_NMASK (SC_N - 1u)

enum {
    SC_FRAME = 1u,
    SC_PIXEL = 2u,
    SC_ROT   = 4u,
    SC_PERM  = 8u
};

typedef struct { double u, v; } sc_point;

/* ---- hashing ----------------------------------------------------------- */

static uint64_t mix64(uint64_t z) {
    z += 0x9E3779B97F4A7C15ull;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}

static uint64_t mix_key(uint64_t k, uint64_t v, uint64_t salt) {
    return mix64(k ^ mix64(v * 0x9E3779B97F4A7C15ull + salt));
}

/* uniform double in [0,1): 53 random bits, exactly representable */
static double hash_unit(uint64_t h) {
    return (double)(h >> 11) * 0x1.0p-53;
}

static uint32_t reverse_bits(uint32_t v) {
    uint32_t r = 0;
    for (unsigned i = 0; i < SC_BITS; i++) { r = (r << 1) | (v & 1u); v >>= 1; }
    return r;
}

/* ---- randomization kernels --------------------------------------------- */

/* nested uniform (Owen) scrambling of one SC_BITS-digit number; a
   bijection for every key, hence the net survives */
static uint32_t owen_permute(uint32_t x, uint64_t key) {
    uint32_t out = 0, prefix = 0;
    for (int lev = (int)SC_BITS - 1; lev >= 0; lev--) {
        uint32_t b = (x >> lev) & 1u;
        uint64_t h = mix64(key ^ mix64((uint64_t)prefix * 2u +
                                       (uint64_t)lev + 0xD1B54A32D192ED03ull));
        out |= (b ^ (uint32_t)(h >> 63)) << lev;
        prefix = (prefix << 1) | b;
    }
    return out;
}

/* ---- sequence generator ------------------------------------------------ */

static void sc_sequence(uint64_t seed, uint32_t frame, uint32_t pixel,
                        unsigned flags, sc_point out[SC_N]) {
    uint64_t k = mix64(seed ^ 0xA5A5A5A5A5A5A5A5ull);
    if (flags & SC_FRAME) k = mix_key(k, (uint64_t)frame + 1u, 0x1111111111111111ull);
    if (flags & SC_PIXEL) k = mix_key(k, (uint64_t)pixel + 1u, 0x2222222222222222ull);

    uint64_t kx = mix_key(k, 1u, 0x3333333333333333ull);
    uint64_t ky = mix_key(k, 2u, 0x4444444444444444ull);
    double rot_u = hash_unit(mix_key(kx, 3u, 0x5555555555555555ull));
    double rot_v = hash_unit(mix_key(ky, 4u, 0x6666666666666666ull));

    for (uint32_t i = 0; i < SC_N; i++) {
        /* base set: a (0,4,2)-net -- x in natural digit order, y bit-reversed
           Gray code, so every 4x4 cell and every 1x16 / 16x1 strip holds
           exactly one point */
        uint32_t x = i;
        uint32_t y = reverse_bits(i ^ (i >> 1));
        if (flags & SC_PERM) {
            x = owen_permute(x, kx);
            y = owen_permute(y, ky);
        }
        double u = (x + 0.5) / (double)SC_N;
        double v = (y + 0.5) / (double)SC_N;
        if (flags & SC_ROT) {
            /* Cranley-Patterson: toroidal shift of the whole set */
            u += rot_u; u -= floor(u);
            v += rot_v; v -= floor(v);
        }
        out[i].u = u;
        out[i].v = v;
    }
}

/* ---- sequence statistics ----------------------------------------------- */

/* index-matched Pearson correlation between two sequences (both dims) */
static double corr_points(const sc_point *a, const sc_point *b) {
    double ma = 0.0, mb = 0.0, va = 0.0, vb = 0.0, cab = 0.0;
    for (uint32_t i = 0; i < SC_N; i++) { ma += a[i].u + a[i].v; mb += b[i].u + b[i].v; }
    ma /= 2.0 * SC_N; mb /= 2.0 * SC_N;
    for (uint32_t i = 0; i < SC_N; i++) {
        double da0 = a[i].u - ma, da1 = a[i].v - ma;
        double db0 = b[i].u - mb, db1 = b[i].v - mb;
        va += da0 * da0 + da1 * da1;
        vb += db0 * db0 + db1 * db1;
        cab += da0 * db0 + da1 * db1;
    }
    if (va < 1e-30 || vb < 1e-30) return 0.0;
    return cab / sqrt(va * vb);
}

/* horizontally adjacent pixels, same frame, same sample index */
static double neighbor_corr(unsigned flags, uint64_t seed, uint32_t frame,
                            int w, int h, int step) {
    double sum = 0.0; int n = 0;
    for (int py = 0; py < h; py += step) {
        for (int px = 0; px + 1 < w; px += step) {
            uint32_t p = (uint32_t)(py * w + px);
            sc_point a[SC_N], b[SC_N];
            sc_sequence(seed, frame, p, flags, a);
            sc_sequence(seed, frame, p + 1u, flags, b);
            sum += corr_points(a, b); n++;
        }
    }
    return n ? sum / n : 0.0;
}

/* same pixel in two frames, same sample index */
static double frame_corr(unsigned flags, uint64_t seed, uint32_t fa, uint32_t fb,
                         int w, int h, int step) {
    double sum = 0.0; int n = 0;
    for (int py = 0; py < h; py += step) {
        for (int px = 0; px < w; px += step) {
            uint32_t p = (uint32_t)(py * w + px);
            sc_point a[SC_N], b[SC_N];
            sc_sequence(seed, fa, p, flags, a);
            sc_sequence(seed, fb, p, flags, b);
            sum += corr_points(a, b); n++;
        }
    }
    return n ? sum / n : 0.0;
}

/* ---- demo scene and artifact metric ------------------------------------ */

#define IMG_W 160
#define IMG_H 120
#define BLOCK 4
#define DIGEST_PX 1024

#define CELL 3.3

static double scene(double x, double y) {
    long cx = (long)floor(x / CELL), cy = (long)floor(y / CELL);
    double checker = ((cx + cy) & 1L) ? 0.0 : 1.0;   /* even cells are lit */
    double shade = 0.65 + 0.35 * (y / (double)IMG_H);   /* slow shading */
    return (0.15 + 0.70 * checker) * shade;
}

/* length of [x0, x0+1) covered by x cells of even parity */
static double even_cell_coverage(double x0) {
    double s = 0.0;
    double b = floor(x0 / CELL) * CELL;
    while (b < x0 + 1.0) {
        double lo = b > x0 ? b : x0;
        double hi = b + CELL < x0 + 1.0 ? b + CELL : x0 + 1.0;
        if (((long)floor((lo + 1e-9) / CELL) & 1L) == 0) s += hi - lo;
        b += CELL;
    }
    return s;
}

/* Analytic pixel average: the checker parity is separable per cell and the
   shading is linear in y, so integrating y cell by cell and using the x
   coverage per y parity is exact.  A supersampled reference would leave a
   numerical floor that every error metric below could not go under. */
static double reference_pixel(int px, int py) {
    double acc = 0.0;
    double x_even = even_cell_coverage((double)px);
    double y0 = (double)py, y1 = y0 + 1.0;
    double b = floor(y0 / CELL) * CELL;
    while (b < y1) {
        double lo = b > y0 ? b : y0;
        double hi = b + CELL < y1 ? b + CELL : y1;
        double ym = 0.5 * (lo + hi);               /* shade is linear in y */
        double shade = 0.65 + 0.35 * ym / (double)IMG_H;
        /* even y cell: on where x is even; odd y cell: on where x is odd */
        double frac = (((long)floor((lo + 1e-9) / CELL) & 1L) == 0)
                          ? x_even
                          : 1.0 - x_even;
        acc += (hi - lo) * (0.15 + 0.70 * frac) * shade;
        b += CELL;
    }
    return acc;
}

static void reference_image(double *ref, double *lo_out, double *hi_out) {
    double lo = 1e30, hi = -1e30;
    for (int py = 0; py < IMG_H; py++) {
        for (int px = 0; px < IMG_W; px++) {
            double v = reference_pixel(px, py);
            ref[py * IMG_W + px] = v;
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
    }
    *lo_out = lo; *hi_out = hi;
}

typedef struct {
    double rmse_frame;     /* mean per-frame RMSE */
    double rmse_mean;      /* RMSE of the temporal mean */
    double persist;        /* rmse_mean / rmse_frame; 1 = never averages */
    double err_fcorr;      /* error correlation first vs last frame */
    int pattern_px;        /* pixels above threshold in temporal mean error */
    int pattern_blk;       /* blocks above threshold, temporal mean */
    int total_px;
} sc_stats;

/* ---- pattern detection -------------------------------------------------
   The pattern is the part of the error that survives averaging:
   block-averaging removes per-pixel noise but not correlated structure, and
   averaging over frames removes independent frame noise but not the
   component reused frame after frame.  An error that exceeds TAU_FRAC of the
   scene dynamic range is counted as a pattern pixel/block.  */
#define TAU_FRAC 0.01

static int pattern_blocks(const double *err, double tau) {
    int hb = IMG_H / BLOCK, wb = IMG_W / BLOCK, hit = 0;
    for (int by = 0; by < hb; by++) {
        for (int bx = 0; bx < wb; bx++) {
            double s = 0.0;
            for (int y = 0; y < BLOCK; y++)
                for (int x = 0; x < BLOCK; x++)
                    s += err[(by * BLOCK + y) * IMG_W + bx * BLOCK + x];
            if (fabs(s / (BLOCK * BLOCK)) > tau) hit++;
        }
    }
    return hit * BLOCK * BLOCK;
}

/* same criterion per pixel, no block averaging: reported for context only,
   since single-pixel noise can cross the threshold by chance */
static int pattern_pixels(const double *err, double tau) {
    int n = IMG_W * IMG_H, hit = 0;
    for (int i = 0; i < n; i++)
        if (fabs(err[i]) > tau) hit++;
    return hit;
}

static void render_stats(unsigned flags, uint64_t seed, const double *ref,
                         double tau, int frames, sc_stats *st) {
    int n = IMG_W * IMG_H;
    double *acc = calloc((size_t)n, sizeof(double));
    double *first = malloc((size_t)n * sizeof(double));
    double *last = malloc((size_t)n * sizeof(double));
    double *err = malloc((size_t)n * sizeof(double));
    if (!acc || !first || !last || !err) { fprintf(stderr, "oom\n"); exit(2); }

    double rms_frames = 0.0;
    for (int f = 0; f < frames; f++) {
        double rms = 0.0;
        for (int py = 0; py < IMG_H; py++) {
            for (int px = 0; px < IMG_W; px++) {
                uint32_t pix = (uint32_t)(py * IMG_W + px);
                sc_point pts[SC_N];
                sc_sequence(seed, (uint32_t)f, pix, flags, pts);
                double s = 0.0;
                for (uint32_t i = 0; i < SC_N; i++)
                    s += scene(px + pts[i].u, py + pts[i].v);
                double v = s / SC_N;
                acc[pix] += v;
                double e = v - ref[pix];
                rms += e * e;
                if (f == 0) first[pix] = e;
                if (f == frames - 1) last[pix] = e;
            }
        }
        rms_frames += sqrt(rms / n);
    }
    st->rmse_frame = rms_frames / frames;

    double rms_mean = 0.0;
    for (int i = 0; i < n; i++) {
        err[i] = acc[i] / frames - ref[i];
        rms_mean += err[i] * err[i];
    }
    st->rmse_mean = sqrt(rms_mean / n);
    st->persist = st->rmse_frame > 0 ? st->rmse_mean / st->rmse_frame : 0.0;

    st->pattern_px = pattern_pixels(err, tau);
    st->pattern_blk = pattern_blocks(err, tau);
    st->total_px = n;

    double ma = 0.0, mb = 0.0, va = 0.0, vb = 0.0, cab = 0.0;
    for (int i = 0; i < n; i++) { ma += first[i]; mb += last[i]; }
    ma /= n; mb /= n;
    for (int i = 0; i < n; i++) {
        double da = first[i] - ma, db = last[i] - mb;
        va += da * da; vb += db * db; cab += da * db;
    }
    st->err_fcorr = (va > 1e-30 && vb > 1e-30) ? cab / sqrt(va * vb) : 0.0;

    free(acc); free(first); free(last); free(err);
}

/* ---- image dump -------------------------------------------------------- */

/* write a grayscale PGM of the temporal-mean error, rescaled so that zero
   error is mid-gray; optional: pattern enhancement is not applied */
static int write_error_pgm(const char *path, unsigned flags, uint64_t seed,
                           const double *ref, int frames) {
    int n = IMG_W * IMG_H;
    double *acc = calloc((size_t)n, sizeof(double));
    if (!acc) return 0;
    for (int f = 0; f < frames; f++) {
        for (int py = 0; py < IMG_H; py++) {
            for (int px = 0; px < IMG_W; px++) {
                uint32_t pix = (uint32_t)(py * IMG_W + px);
                sc_point pts[SC_N];
                sc_sequence(seed, (uint32_t)f, pix, flags, pts);
                double s = 0.0;
                for (uint32_t i = 0; i < SC_N; i++)
                    s += scene(px + pts[i].u, py + pts[i].v);
                acc[pix] += s / SC_N;
            }
        }
    }
    double peak = 0.0;
    for (int i = 0; i < n; i++) {
        acc[i] = acc[i] / frames - ref[i];
        double a = fabs(acc[i]);
        if (a > peak) peak = a;
    }
    if (peak < 1e-12) peak = 1e-12;
    FILE *f = fopen(path, "wb");
    if (!f) { free(acc); return 0; }
    fprintf(f, "P5\n%d %d\n255\n", IMG_W, IMG_H);
    for (int i = 0; i < n; i++) {
        int v = (int)(128.0 + 127.0 * acc[i] / peak);
        if (v < 0) v = 0;
        if (v > 255) v = 255;
        fputc(v, f);
    }
    fclose(f);
    free(acc);
    return 1;
}

/* ---- deterministic digest ---------------------------------------------- */

static uint64_t frame_digest(unsigned flags, uint64_t seed, uint32_t frame,
                             uint32_t pixels) {
    uint64_t h = 1469598103934665603ull;
    for (uint32_t p = 0; p < pixels; p++) {
        sc_point pts[SC_N];
        sc_sequence(seed, frame, p, flags, pts);
        const unsigned char *b = (const unsigned char *)pts;
        for (size_t i = 0; i < sizeof(pts); i++) { h ^= b[i]; h *= 1099511628211ull; }
    }
    return h;
}

static uint64_t checksum(unsigned flags, uint64_t seed, uint32_t frames, uint32_t pixels) {
    uint64_t h = 1469598103934665603ull;
    for (uint32_t f = 0; f < frames; f++) {
        uint64_t d = frame_digest(flags, seed, f, pixels);
        h ^= d;
        h *= 1099511628211ull;
    }
    return h;
}

/* ---- self checks ------------------------------------------------------- */

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* a (0,4,2)-net: each 4x4 cell, each column strip and each row strip holds
   exactly one point -- this is what "stratified inside the pixel" means.
   Note this holds for SC_PERM; SC_ROT is a continuous toroidal shift and
   deliberately trades strict stratification for decorrelation. */
static int is_net(const sc_point pts[SC_N]) {
    unsigned cell[4][4], row[SC_N], col[SC_N];
    memset(cell, 0, sizeof(cell));
    memset(row, 0, sizeof(row));
    memset(col, 0, sizeof(col));
    for (uint32_t i = 0; i < SC_N; i++) {
        unsigned x = (unsigned)(pts[i].u * SC_N), y = (unsigned)(pts[i].v * SC_N);
        if (x >= SC_N || y >= SC_N) return 0;
        cell[y / 4][x / 4]++;
        col[x]++; row[y]++;
    }
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
            if (cell[a][b] != 1) return 0;
    for (uint32_t i = 0; i < SC_N; i++)
        if (row[i] != 1 || col[i] != 1) return 0;
    return 1;
}

static void run_checks(void) {
    const unsigned all = SC_FRAME | SC_PIXEL | SC_ROT | SC_PERM;
    sc_point a[SC_N], b[SC_N];

    /* golden digest: if the stream changes, every previously rendered frame
       changes, so this pins the generator.  Recompute only on purpose. */
    check("golden_seed42", checksum(all, 42, 8, DIGEST_PX) == 0xd65062c7604ae784ull);

    /* permutation alone preserves the net for every seed/frame/pixel */
    int net_ok = 1;
    for (uint64_t seed = 0; seed < 8; seed++)
        for (uint32_t pix = 0; pix < 64; pix += 9)
            for (uint32_t f = 0; f < 3; f++) {
                sc_sequence(seed, f, pix, SC_PERM, a);
                if (!is_net(a)) net_ok = 0;
                sc_sequence(seed, f, pix, SC_PERM | SC_FRAME | SC_PIXEL, a);
                if (!is_net(a)) net_ok = 0;
            }
    check("perm_preserves_net", net_ok);

    /* Owen scrambling is a bijection: no sample lands on another's cell */
    int bij = 1;
    for (uint64_t seed = 0; seed < 8; seed++) {
        sc_sequence(seed, 0, 0, SC_PERM, a);
        unsigned seen[SC_N];
        memset(seen, 0, sizeof(seen));
        for (uint32_t i = 0; i < SC_N; i++) {
            unsigned x = (unsigned)(a[i].u * SC_N);
            if (x >= SC_N || seen[x]++) bij = 0;
        }
    }
    check("perm_bijection", bij);

    /* rotation keeps every sample inside the pixel and changes frames */
    sc_sequence(5, 0, 0, SC_ROT, a);
    int in_unit = 1;
    for (uint32_t i = 0; i < SC_N; i++)
        if (!(a[i].u >= 0.0 && a[i].u < 1.0 && a[i].v >= 0.0 && a[i].v < 1.0))
            in_unit = 0;
    check("rotation_in_unit", in_unit);

    /* reuse mode: identical sequence for every pixel and frame */
    sc_sequence(123, 0, 0, 0, a);
    sc_sequence(123, 7, 555, 0, b);
    check("reuse_identical", memcmp(a, b, sizeof(a)) == 0);

    /* randomization mode: changes across pixels and across frames */
    sc_sequence(123, 0, 0, all, a);
    sc_sequence(123, 0, 1, all, b);
    int pixel_changes = memcmp(a, b, sizeof(a)) != 0;
    sc_sequence(123, 1, 0, all, b);
    check("pixel_changes", pixel_changes);
    check("frame_changes", memcmp(a, b, sizeof(a)) != 0);

    /* fixed seed reproduces the same bytes, and order does not matter */
    uint64_t c0 = checksum(all, 42, 8, 1024);
    uint64_t c1 = checksum(all, 42, 8, 1024);
    check("deterministic", c0 == c1);
    check("seed_matters", c0 != checksum(all, 43, 8, 1024));

    sc_point f2[SC_N], f2b[SC_N];
    sc_sequence(7, 2, 10, all, f2);
    (void)checksum(all, 7, 8, 64);          /* evaluate other frames first */
    sc_sequence(7, 2, 10, all, f2b);
    check("frame_reproducible", memcmp(f2, f2b, sizeof(f2)) == 0);

    /* sequence correlation */
    double r_neigh = neighbor_corr(all, 5, 0, 96, 72, 7);
    double r_neigh_reuse = neighbor_corr(0, 5, 0, 96, 72, 7);
    check("neighbor_decorrelated", fabs(r_neigh) < 0.10);
    check("neighbor_corr_reuse", r_neigh_reuse > 0.999);

    double r_frame = frame_corr(all, 5, 0, 1, 96, 72, 7);
    double r_frame_pix = frame_corr(SC_PIXEL | SC_ROT | SC_PERM, 5, 0, 1, 96, 72, 7);
    check("frame_decorrelated", fabs(r_frame) < 0.10);
    check("pixel_mode_static_frames", r_frame_pix > 0.999);

    /* artifact metric: randomization cuts the persistent pattern */
    double *ref = malloc(sizeof(double) * IMG_W * IMG_H);
    if (!ref) { fprintf(stderr, "oom\n"); exit(2); }
    double lo, hi;
    reference_image(ref, &lo, &hi);
    double tau = TAU_FRAC * (hi - lo);
    sc_stats sr, sf;
    render_stats(0, 20261122, ref, tau, 16, &sr);
    render_stats(all, 20261122, ref, tau, 16, &sf);
    check("pattern_drops", sf.pattern_blk * 4 < sr.pattern_blk);
    check("persistence_drops", sf.persist < 0.7 * sr.persist);
    free(ref);
}

/* ---- report ------------------------------------------------------------ */

static const struct { const char *name; unsigned flags; } MODES[] = {
    { "reuse",        0 },
    { "frame_rot",    SC_FRAME | SC_ROT },
    { "frame_perm",   SC_FRAME | SC_PERM },
    { "frame_both",   SC_FRAME | SC_ROT | SC_PERM },
    { "pixel_rot",    SC_PIXEL | SC_ROT },
    { "pixel_perm",   SC_PIXEL | SC_PERM },
    { "pixel_both",   SC_PIXEL | SC_ROT | SC_PERM },
    { "full",         SC_FRAME | SC_PIXEL | SC_ROT | SC_PERM },
};
#define NMODES (int)(sizeof(MODES) / sizeof(MODES[0]))

static void bench(uint64_t seed, int frames) {
    double *ref = malloc(sizeof(double) * IMG_W * IMG_H);
    if (!ref) { fprintf(stderr, "oom\n"); exit(2); }
    double lo, hi;
    reference_image(ref, &lo, &hi);
    double tau = TAU_FRAC * (hi - lo);

    printf("seed=%llu frames=%d samples=%u image=%dx%d block=%d tau=%.0f%% of range %.3f\n\n",
           (unsigned long long)seed, frames, SC_N, IMG_W, IMG_H, BLOCK,
           TAU_FRAC * 100.0, hi - lo);

    printf("| mode | seq corr neighbor | seq corr frame 0-1 | err corr first-last | rmse/frame | rmse/mean | persistence | pattern px | pattern blk | pattern %% |\n");
    printf("|---|---|---|---|---|---|---|---|---|---|\n");
    for (int i = 0; i < NMODES; i++) {
        sc_stats st;
        render_stats(MODES[i].flags, seed, ref, tau, frames, &st);
        double rn = neighbor_corr(MODES[i].flags, seed, 0, IMG_W, IMG_H, 11);
        double rf = frame_corr(MODES[i].flags, seed, 0, 1, IMG_W, IMG_H, 11);
        printf("| %s | %+.3f | %+.3f | %+.3f | %.5f | %.5f | %.3f | %d | %d | %.2f%% |\n",
               MODES[i].name, rn, rf, st.err_fcorr, st.rmse_frame, st.rmse_mean,
               st.persist, st.pattern_px, st.pattern_blk,
               100.0 * st.pattern_blk / st.total_px);
    }

    printf("\n## persistent pattern vs frame count (temporal mean error)\n\n");
    printf("| mode | frames | rmse/mean | persistence | pattern px | pattern blk | pattern %% |\n");
    printf("|---|---|---|---|---|---|---|\n");
    int fs[] = { 1, 4, 16, 64 };
    unsigned picks[] = { 0u, SC_FRAME | SC_PIXEL | SC_ROT | SC_PERM };
    const char *names[] = { "reuse", "full" };
    for (int m = 0; m < 2; m++) {
        for (int k = 0; k < 4; k++) {
            sc_stats st;
            render_stats(picks[m], seed, ref, tau, fs[k], &st);
            printf("| %s | %d | %.5f | %.3f | %d | %d | %.2f%% |\n", names[m], fs[k],
                   st.rmse_mean, st.persist, st.pattern_px, st.pattern_blk,
                   100.0 * st.pattern_blk / st.total_px);
        }
    }

    printf("\n## frame digests (seed=%llu, pixels=%d; rerun must match)\n\n",
           (unsigned long long)seed, DIGEST_PX);
    printf("| frame | digest |\n|---|---|\n");
    for (uint32_t f = 0; f < 8; f++)
        printf("| %u | %016llx |\n", f,
               (unsigned long long)frame_digest(SC_FRAME | SC_PIXEL | SC_ROT | SC_PERM,
                                                seed, f, DIGEST_PX));
    printf("\nfull checksum=%016llx\n",
           (unsigned long long)checksum(SC_FRAME | SC_PIXEL | SC_ROT | SC_PERM,
                                        seed, 8, DIGEST_PX));
    free(ref);
}

static void usage(const char *prog) {
    fprintf(stderr,
            "usage: %s                      run self checks\n"
            "       %s bench [seed] [frames] artifact / correlation table\n"
            "       %s checksum [seed]       one reproducible digest\n"
            "       %s dump [seed] [frames]  error PGMs for reuse and full\n",
            prog, prog, prog, prog);
}

int main(int argc, char **argv) {
    if (argc < 2 || strcmp(argv[1], "check") == 0) {
        run_checks();
        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }
    uint64_t seed = 20261122ull;
    if (argc >= 3) seed = strtoull(argv[2], NULL, 10);
    if (strcmp(argv[1], "bench") == 0) {
        int frames = argc >= 4 ? atoi(argv[3]) : 16;
        if (frames < 1) frames = 1;
        bench(seed, frames);
        return 0;
    }
    if (strcmp(argv[1], "checksum") == 0) {
        printf("checksum=%016llx\n",
               (unsigned long long)checksum(SC_FRAME | SC_PIXEL | SC_ROT | SC_PERM,
                                            seed, 8, DIGEST_PX));
        return 0;
    }
    if (strcmp(argv[1], "dump") == 0) {
        int frames = argc >= 4 ? atoi(argv[3]) : 16;
        if (frames < 1) frames = 1;
        double *ref = malloc(sizeof(double) * IMG_W * IMG_H);
        if (!ref) { fprintf(stderr, "oom\n"); exit(2); }
        double lo, hi;
        reference_image(ref, &lo, &hi);
        int a = write_error_pgm("err_reuse.pgm", 0, seed, ref, frames);
        int b = write_error_pgm("err_full.pgm",
                                SC_FRAME | SC_PIXEL | SC_ROT | SC_PERM,
                                seed, ref, frames);
        free(ref);
        if (!a || !b) { fprintf(stderr, "cannot write pgm\n"); return 2; }
        printf("wrote err_reuse.pgm err_full.pgm (frames=%d seed=%llu)\n",
               frames, (unsigned long long)seed);
        return 0;
    }
    usage(argv[0]);
    return 2;
}
