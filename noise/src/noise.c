/* noise: per-pixel noise estimation from a local pixel neighbourhood.
 *
 * The renderer stores, per pixel and channel, the mean colour and how many
 * samples produced it. What makes the image shimmer is the variance of that
 * mean, which no single pixel can reveal. It has to be inferred from how much
 * nearby pixels disagree.
 *
 * Third differences. Along a row take four consecutive pixels and form
 *     T = m[x] - 3 m[x+1] + 3 m[x+2] - m[x+3].
 * T is exactly zero for any linear or quadratic signal, so a smooth image
 * leaves only noise. With independent per-pixel noise of variance s2,
 *     Var(T) = (1 + 9 + 9 + 1) s2 = 20 s2.
 * Over n mutually disjoint quadruples, S = sum(T^2) satisfies
 *     S / (20 s2) ~ chi2(n),
 * so est = S/(20 n) is unbiased and S is sufficient for s2. Vertical
 * quadruples give a second, independent set of the same kind.
 *
 * The bound. S/(20 s2) is chi2(n), so S/20 has variance 2 n s2^2 and each
 * orientation's estimate estH = S/(20 n) has standard error sqrt(2 estH^2/n).
 * The two orientations are averaged, so
 *     se = 0.5 * sqrt(2 estH^2/nH + 2 estV^2/nV),
 * and the reported bound is z * se. The orientations share no quadruples but
 * do share a few pixels, which only makes the true se smaller; z = 3 therefore
 * gives a conservative ~99% guarantee per pixel under Gaussian noise with a
 * locally cubic signal. The bound widens as the window shrinks, which is the
 * honest behaviour of any variance estimate.
 *
 * What it targets. The estimate is a weighted window average of per-pixel
 * variance, with the quadruple weights 1, 9, 9, 1. The self-test compares it
 * against that same weighted truth, not the centre pixel. A plain window
 * average over the square would be a different target and mismatch wherever
 * the noise level varies within the window.
 *
 * Not estimable. A window clipped by the image border or made small by a
 * short radius may hold fewer than MIN_QUAD quadruples in either orientation.
 * There the noise cannot be separated from curvature with any confidence, so
 * the pixel is flagged, given no bound, and assigned the mean weight by the
 * allocator rather than a value the estimator cannot support.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NOISE_K 3
#define MIN_QUAD 10        /* quadruples per orientation needed to report */
#define NOISE_RADIUS 6     /* default halfwidth: 13x13 window */
#define NOISE_PI 3.14159265358979323846

typedef struct {
    double v[NOISE_K];
} Pixel;

typedef struct {
    int ok;                /* 0 => neighbourhood too small, no estimate */
    int nh, nv;            /* disjoint quadruples per orientation */
    double est[NOISE_K];   /* point estimate of per-pixel variance */
    double val[NOISE_K];   /* max(est, 0): the noise-map value */
    double lo[NOISE_K];    /* lower end of the confidence bound */
    double hi[NOISE_K];    /* upper end */
    double bound[NOISE_K]; /* z standard errors of `est` */
} NoisePixel;

/* ---- estimator ---------------------------------------------------------- */

/* Sum of squared third differences over disjoint quadruples in one
   orientation, per channel. Returns the number of quadruples. */
static int quad_scan(const Pixel* img, int w, int x0, int x1, int y0, int y1,
                     int horizontal, double S[NOISE_K]) {
    for (int c = 0; c < NOISE_K; c++) S[c] = 0.0;
    int n = 0;
    if (horizontal) {
        for (int y = y0; y <= y1; y++)
            for (int x = x0; x + 3 <= x1; x += 4) {
                size_t p0 = (size_t)y * w + x;
                for (int c = 0; c < NOISE_K; c++) {
                    double t = img[p0].v[c] - 3.0 * img[p0 + 1].v[c] +
                               3.0 * img[p0 + 2].v[c] - img[p0 + 3].v[c];
                    S[c] += t * t;
                }
                n++;
            }
        return n;
    }
    for (int x = x0; x <= x1; x++)
        for (int y = y0; y + 3 <= y1; y += 4) {
            size_t p0 = (size_t)y * w + x;
            for (int c = 0; c < NOISE_K; c++) {
                double t = img[p0].v[c] - 3.0 * img[p0 + w].v[c] +
                           3.0 * img[p0 + 2 * w].v[c] - img[p0 + 3 * w].v[c];
                S[c] += t * t;
            }
            n++;
        }
    return n;
}

/* Estimate the noise map. radius is the window halfwidth; z multiplies the
   standard error to form the bound. Returns the number of estimated pixels. */
static int noise_estimate(const Pixel* img, int w, int h, int k, int radius,
                          double z, NoisePixel* out) {
    int estimable = 0;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            NoisePixel* np = &out[(size_t)y * w + x];
            memset(np, 0, sizeof(*np));
            int x0 = x - radius < 0 ? 0 : x - radius;
            int x1 = x + radius >= w ? w - 1 : x + radius;
            int y0 = y - radius < 0 ? 0 : y - radius;
            int y1 = y + radius >= h ? h - 1 : y + radius;

            double SH[NOISE_K], SV[NOISE_K];
            int nh = quad_scan(img, w, x0, x1, y0, y1, 1, SH);
            int nv = quad_scan(img, w, x0, x1, y0, y1, 0, SV);
            np->nh = nh;
            np->nv = nv;
            np->ok = (nh >= MIN_QUAD && nv >= MIN_QUAD);
            if (!np->ok) continue;   /* lo/hi/bound stay zero: no claim */

            for (int c = 0; c < k; c++) {
                double estH = SH[c] / (20.0 * (double)nh);
                double estV = SV[c] / (20.0 * (double)nv);
                double est = 0.5 * (estH + estV);
                /* S/(20 s2) is chi2(n), so Var(S/20) = 2 n s2^2 and
                   Var(estH) = 2 estH^2 / n. The orientations are averaged. */
                double vH = 2.0 * estH * estH / (double)nh;
                double vV = 2.0 * estV * estV / (double)nv;
                double bound = z * 0.5 * sqrt(vH + vV);
                np->est[c] = est;
                np->val[c] = est > 0.0 ? est : 0.0;
                np->lo[c] = est > bound ? est - bound : 0.0;
                np->hi[c] = est + bound;
                np->bound[c] = bound;
            }
            estimable++;
        }
    }
    return estimable;
}

/* ---- adaptive sampling --------------------------------------------------- */

static double noise_weight(const NoisePixel* np, int k) {
    double m = 0.0;
    for (int c = 0; c < k; c++) if (np->val[c] > m) m = np->val[c];
    return m;
}

/* Sample counts proportional to the estimated variance. The residual variance
   of the mean after N samples is var/N, so N proportional to var leaves a
   roughly constant residual level across the image. Pixels that could not be
   estimated take the mean weight so a missing estimate does not starve them. */
static void noise_allocate(const NoisePixel* map, int n, int k, int base,
                           int min_n, int max_n, int* alloc) {
    double sum = 0.0;
    int okn = 0;
    for (int i = 0; i < n; i++)
        if (map[i].ok) { sum += noise_weight(&map[i], k); okn++; }
    double mean_w = okn ? sum / (double)okn : 1.0;
    if (mean_w <= 0.0) mean_w = 1.0;
    for (int i = 0; i < n; i++) {
        double w = map[i].ok ? noise_weight(&map[i], k) : mean_w;
        int a = (int)floor((double)base * w / mean_w + 0.5);
        if (a < min_n) a = min_n;
        if (a > max_n) a = max_n;
        alloc[i] = a;
    }
}

/* ---- noise map output ---------------------------------------------------- */

/* Text rendering of one channel, block-averaged by `step`, on a fixed scale
   `hi` so two maps can be compared. Pass scale = 0 to fit the map's own max.
   '?' marks pixels the estimator refused to characterise. */
static void print_noise_map(const NoisePixel* map, int w, int h, int ch, int step,
                            double scale) {
    static const char* ramp = " .:-=+*#%@";
    double hi = scale;
    if (hi <= 0.0)
        for (int i = 0; i < w * h; i++)
            if (map[i].ok && map[i].val[ch] > hi) hi = map[i].val[ch];
    if (hi <= 0.0) hi = 1.0;
    for (int y = 0; y < h; y += step) {
        for (int x = 0; x < w; x += step) {
            double s = 0.0;
            int n = 0, any = 0;
            for (int yy = y; yy < y + step && yy < h; yy++)
                for (int xx = x; xx < x + step && xx < w; xx++) {
                    const NoisePixel* np = &map[(size_t)yy * w + xx];
                    if (np->ok) { s += np->val[ch]; n++; any = 1; }
                }
            if (!any) { putchar('?'); continue; }
            double v = s / n;
            int idx = (int)(v / hi * 9.0 + 0.5);
            if (idx < 0) idx = 0;
            if (idx > 9) idx = 9;
            putchar(ramp[idx]);
        }
        putchar('\n');
    }
}

/* ---- deterministic synthetic scenes (tests only) ------------------------ */

static unsigned wang_hash(unsigned x) {
    x = (x ^ 61u) ^ (x >> 16); x *= 9u;
    x ^= x >> 4; x *= 0x27d4eb2du;
    return x ^ (x >> 15);
}

static double gauss_at(unsigned key) {
    unsigned s = wang_hash(key ^ 0x9e3779b9u);
    double u1 = ((double)(s & 0xFFFFFFu) + 1.0) / 16777217.0;
    s = wang_hash(s);
    double u2 = ((double)(s & 0xFFFFFFu) + 1.0) / 16777217.0;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * NOISE_PI * u2);
}

static double smoothstep(double t) { return t * t * (3.0 - 2.0 * t); }

static double lattice(int ix, int iy, unsigned seed) {
    unsigned h = wang_hash(((unsigned)ix * 73856093u) ^
                           ((unsigned)iy * 19349663u) ^ seed);
    return (double)(h & 0xFFFFu) / 65535.0;
}

static double value_noise(double x, double y, unsigned seed) {
    int ix = (int)floor(x), iy = (int)floor(y);
    double fx = smoothstep(x - (double)ix), fy = smoothstep(y - (double)iy);
    double a = lattice(ix, iy, seed) +
               (lattice(ix + 1, iy, seed) - lattice(ix, iy, seed)) * fx;
    double b = lattice(ix, iy + 1, seed) +
               (lattice(ix + 1, iy + 1, seed) - lattice(ix, iy + 1, seed)) * fx;
    return a + (b - a) * fy;
}

typedef struct {
    int w, h, k;
    double* signal;
    double* var;      /* true per-pixel, per-channel variance of the mean */
} Scene;

static int scene_alloc(Scene* s, int w, int h, int k) {
    if (k < 1 || k > NOISE_K) return 0;
    s->w = w; s->h = h; s->k = k;
    s->signal = (double*)calloc((size_t)w * h * k, sizeof(double));
    s->var = (double*)calloc((size_t)w * h * k, sizeof(double));
    return s->signal && s->var;
}

static void scene_free(Scene* s) { free(s->signal); free(s->var); }

static void scene_render(const Scene* sc, unsigned seed, Pixel* out) {
    for (int y = 0; y < sc->h; y++)
        for (int x = 0; x < sc->w; x++) {
            size_t p = (size_t)y * sc->w + x;
            for (int c = 0; c < sc->k; c++) {
                size_t i = p * sc->k + c;
                unsigned key = wang_hash((unsigned)p * 2654435761u ^
                                         (unsigned)(c + 1) * 2246822519u ^ seed);
                double v = sc->var[i] > 0.0 ? sqrt(sc->var[i]) * gauss_at(key) : 0.0;
                out[p].v[c] = sc->signal[i] + v;
            }
        }
}

/* ---- self-test ---------------------------------------------------------- */

#define SW 25
#define SH 25

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void fill_uniform(Scene* sc, double level, double var) {
    for (int i = 0; i < sc->w * sc->h * sc->k; i++) {
        sc->signal[i] = level * (1.0 + 0.1 * (double)(i % sc->k));
        sc->var[i] = var;
    }
}

/* A slowly varying truth, so a window average is comparable to the pixel. */
static void fill_graded(Scene* sc, unsigned seed) {
    for (int y = 0; y < sc->h; y++)
        for (int x = 0; x < sc->w; x++) {
            double s = value_noise(x / 11.0, y / 11.0, seed) - 0.5;
            double v = 0.06 + 0.60 * value_noise(x / 10.0, y / 10.0, seed + 7u);
            for (int c = 0; c < sc->k; c++) {
                size_t i = ((size_t)y * sc->w + x) * sc->k + c;
                sc->signal[i] = 1.0 + 0.4 * s;
                sc->var[i] = v * (1.0 + 0.2 * c);
            }
        }
}

/* The variance the estimator targets is not a plain window average: E[T^2]
   weights the four pixels of a quadruple by 1, 9, 9, 1. Mirror the geometry of
   quad_scan to accumulate, for one window, each pixel's total weight, and take
   the weighted mean of the per-pixel variance. This is the exact truth a
   noiseless estimate would return. */
static void quad_truth(int x0, int x1, int y0, int y1, int bw,
                       double* weight, double* wsum) {
    static const double c[4] = { 1, 9, 9, 1 };
    for (int i = 0; i < bw * (y1 - y0 + 1); i++) weight[i] = 0.0;
    *wsum = 0.0;
    for (int y = y0; y <= y1; y++)
        for (int x = x0; x + 3 <= x1; x += 4) {
            size_t p = (size_t)(y - y0) * bw + (x - x0);
            for (int i = 0; i < 4; i++) { weight[p + i] += c[i]; *wsum += c[i]; }
        }
    for (int x = x0; x <= x1; x++)
        for (int y = y0; y + 3 <= y1; y += 4)
            for (int i = 0; i < 4; i++) {
                size_t p = (size_t)(y + i - y0) * bw + (x - x0);
                weight[p] += c[i]; *wsum += c[i];
            }
}

/* Fill out[] (w*h) with the stencil-weighted window variance of `var`. */
static void stencil_truth(const double* var, int w, int h, int radius, double* out) {
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++) {
            int x0 = x - radius < 0 ? 0 : x - radius;
            int x1 = x + radius >= w ? w - 1 : x + radius;
            int y0 = y - radius < 0 ? 0 : y - radius;
            int y1 = y + radius >= h ? h - 1 : y + radius;
            int bw = x1 - x0 + 1, bh = y1 - y0 + 1;
            static double weight[SW * SH];
            double wsum = 0.0;
            quad_truth(x0, x1, y0, y1, bw, weight, &wsum);
            double acc = 0.0;
            for (int yy = 0; yy < bh; yy++)
                for (int xx = 0; xx < bw; xx++)
                    acc += weight[(size_t)yy * bw + xx] * var[(size_t)(y0 + yy) * w + x0 + xx];
            out[(size_t)y * w + x] = wsum > 0.0 ? acc / wsum : 0.0;
        }
}

/* Average the estimator over R independent renderings of the same scene. Each
   estimate is unbiased, so the average suppresses the estimator's own sampling
   noise by sqrt(R) and exposes how well the map follows the true variance. */
static void estimate_avg(const Scene* sc, int radius, double z, int k, int R,
                         NoisePixel* out, Pixel* img, NoisePixel* scratch) {
    int n = sc->w * sc->h;
    for (int i = 0; i < n; i++) {
        memset(&out[i], 0, sizeof(out[i]));
        out[i].ok = 1;
    }
    for (int r = 0; r < R; r++) {
        scene_render(sc, 0x5000u + (unsigned)r * 104729u, img);
        noise_estimate(img, sc->w, sc->h, k, radius, z, scratch);
        for (int i = 0; i < n; i++) {
            if (!scratch[i].ok) { out[i].ok = 0; continue; }
            for (int c = 0; c < k; c++) out[i].est[c] += scratch[i].est[c];
        }
    }
    for (int i = 0; i < n; i++) {
        if (!out[i].ok) continue;
        for (int c = 0; c < k; c++) {
            out[i].est[c] /= (double)R;
            out[i].val[c] = out[i].est[c] > 0.0 ? out[i].est[c] : 0.0;
        }
    }
}

static double pearson(const double* a, const double* b, const int* use, int n) {
    double ma = 0.0, mb = 0.0; int m = 0;
    for (int i = 0; i < n; i++) if (use[i]) { ma += a[i]; mb += b[i]; m++; }
    if (m < 3) return 0.0;
    ma /= m; mb /= m;
    double sab = 0.0, saa = 0.0, sbb = 0.0;
    for (int i = 0; i < n; i++) if (use[i]) {
        double da = a[i] - ma, db = b[i] - mb;
        sab += da * db; saa += da * da; sbb += db * db;
    }
    return (saa > 0.0 && sbb > 0.0) ? sab / sqrt(saa * sbb) : 0.0;
}

/* Coefficient of variation over usable pixels. */
static double cv(const double* v, const int* use, int n) {
    double m = 0.0, m2 = 0.0; int c = 0;
    for (int i = 0; i < n; i++) if (use[i]) { m += v[i]; m2 += v[i] * v[i]; c++; }
    if (c < 2) return 0.0;
    m /= c; m2 /= c;
    double var = m2 - m * m;
    return (var > 0.0 && m > 0.0) ? sqrt(var) / m : 0.0;
}

int main(void) {
    static Pixel img[SW * SH];
    static NoisePixel map[SW * SH];
    static int alloc[SW * SH], use[SW * SH];
    static double truth[SW * SH], estv[SW * SH], resid[SW * SH];
    Scene sc;

    /* 1. constant signal, constant noise: unbiased, and the 95% bound covers
          the truth at close to the advertised rate over realisations. */
    if (!scene_alloc(&sc, SW, SH, 1)) return 2;
    double tvar = 0.25;
    fill_uniform(&sc, 1.0, tvar);
    int covered = 0, total = 0, ordered = 1;
    double sum_est = 0.0, max_bound = 0.0;
    for (unsigned r = 0; r < 6; r++) {
        scene_render(&sc, 0x1000u + r * 7919u, img);
        noise_estimate(img, SW, SH, 1, NOISE_RADIUS, 3.0, map);
        for (int i = 0; i < SW * SH; i++) {
            if (!map[i].ok) continue;
            total++;
            if (fabs(map[i].est[0] - tvar) <= map[i].bound[0]) covered++;
            if (!(map[i].lo[0] >= 0.0 && map[i].lo[0] <= map[i].est[0] &&
                  map[i].est[0] <= map[i].hi[0])) ordered = 0;
            sum_est += map[i].est[0];
            if (map[i].bound[0] > max_bound) max_bound = map[i].bound[0];
        }
    }
    check("estimable_some", total > 0);
    check("interval_ordered", ordered);
    check("coverage_at_least_90", covered >= (int)(0.90 * total));
    check("unbiased", fabs(sum_est / total - tvar) < 0.02);
    check("bound_finite", max_bound < 4.0 * tvar);
    printf("const_coverage=%.3f est_mean=%.4f truth=%.4f max_bound=%.4f\n",
           (double)covered / total, sum_est / total, tvar, max_bound);

    /* 2. the bound is not free: more quadruples in a bigger window must make
          it smaller. */
    fill_uniform(&sc, 1.0, 0.25);
    scene_render(&sc, 0x2222u, img);
    double b4 = 0.0, b8 = 0.0; int n4 = 0, n8 = 0;
    noise_estimate(img, SW, SH, 1, 4, 3.0, map);
    for (int i = 0; i < SW * SH; i++) if (map[i].ok) { b4 += map[i].bound[0]; n4++; }
    noise_estimate(img, SW, SH, 1, 8, 3.0, map);
    for (int i = 0; i < SW * SH; i++) if (map[i].ok) { b8 += map[i].bound[0]; n8++; }
    check("bound_shrinks_with_window", n4 > 0 && n8 > 0 && b8 / n8 < 0.8 * b4 / n4);
    printf("mean_bound_r4=%.4f r8=%.4f\n", b4 / n4, b8 / n8);

    /* 3. smooth signal: the third difference annihilates its curvature, so
          the estimate stays near truth and coverage does not collapse. */
    fill_uniform(&sc, 0.0, 0.20);
    for (int y = 0; y < SH; y++)
        for (int x = 0; x < SW; x++)
            sc.signal[y * SW + x] = 1.0 + 0.8 * (value_noise(x / 18.0, y / 18.0, 0x333u) - 0.5);
    scene_render(&sc, 0x4444u, img);
    noise_estimate(img, SW, SH, 1, NOISE_RADIUS, 3.0, map);
    double ssum = 0.0; int scov = 0, sn = 0;
    for (int i = 0; i < SW * SH; i++) {
        if (!map[i].ok) continue;
        ssum += map[i].est[0]; sn++;
        if (fabs(map[i].est[0] - 0.20) <= map[i].bound[0]) scov++;
    }
    check("smooth_mean_close", fabs(ssum / sn - 0.20) < 0.02);
    check("smooth_coverage", scov >= (int)(0.90 * sn));
    printf("smooth_mean=%.4f coverage=%.3f n=%d\n", ssum / sn, (double)scov / sn, sn);

    /* 4. too small a neighbourhood: nothing estimable, no invented bound. */
    int e = noise_estimate(img, SW, SH, 1, 1, 3.0, map);
    int none = 1;
    for (int i = 0; i < SW * SH; i++)
        if (map[i].ok || map[i].bound[0] != 0.0) none = 0;
    check("tiny_window_not_estimable", e == 0 && none);
    check("radius2_not_estimable", noise_estimate(img, SW, SH, 1, 2, 3.0, map) == 0);
    printf("tiny_window_estimable=%d\n", e);

    /* 5. at radius 4 a clipped corner window holds too few quadruples and is
          flagged, while the same radius in the interior is estimable. */
    noise_estimate(img, SW, SH, 1, 4, 3.0, map);
    check("corner_window_not_estimable", map[0].ok == 0);
    check("interior_window_estimable", map[(SH / 2) * SW + SW / 2].ok == 1);
    printf("corner_nh=%d nv=%d interior_nh=%d nv=%d\n",
           map[0].nh, map[0].nv, map[(SH / 2) * SW + SW / 2].nh,
           map[(SH / 2) * SW + SW / 2].nv);

    /* 6. tiny image: every window is clipped, everything is flagged. */
    {
        static Pixel small[7 * 7];
        static NoisePixel sm[7 * 7];
        Scene ss;
        if (!scene_alloc(&ss, 7, 7, 1)) return 2;
        fill_uniform(&ss, 1.0, 0.25);
        scene_render(&ss, 0x55u, small);
        (void)noise_estimate(small, 7, 7, 1, 3, 3.0, sm);
        int all_flagged = 1;
        for (int i = 0; i < 7 * 7; i++) if (sm[i].ok) all_flagged = 0;
        check("tiny_image_not_estimable", all_flagged);
        scene_free(&ss);
    }

    /* 7. graded noise, multi-channel: the map tracks the truth independently
          per channel, and allocation is monotone in the map. */
    scene_free(&sc);
    if (!scene_alloc(&sc, SW, SH, 3)) return 2;
    fill_graded(&sc, 0xABCDu);
    scene_render(&sc, 0x6666u, img);
    static NoisePixel avg[SW * SH], scratch[SW * SH];
    int ne = noise_estimate(img, SW, SH, 3, NOISE_RADIUS, 3.0, map);
    check("graded_estimable", ne > 0);
    estimate_avg(&sc, NOISE_RADIUS, 3.0, 3, 12, avg, img, scratch);
    static double winvar[SW * SH];
    for (int c = 0; c < 3; c++) {
        /* the estimate averages over the window, so compare it with the
           window-average of the true variance, not the centre pixel */
        for (int i = 0; i < SW * SH; i++) winvar[i] = sc.var[i * 3 + c];
        stencil_truth(winvar, SW, SH, NOISE_RADIUS, truth);
        for (int i = 0; i < SW * SH; i++) {
            estv[i] = avg[i].val[c];
            use[i] = avg[i].ok;
        }
        double r = pearson(estv, truth, use, SW * SH);
        printf("graded_corr_ch%d=%.4f\n", c, r);
        check("graded_corr", r > 0.9);
    }
    noise_allocate(map, SW * SH, 3, 32, 4, 256, alloc);
    int in_range = 1, monotone = 1;
    for (int i = 0; i < SW * SH; i++)
        if (alloc[i] < 4 || alloc[i] > 256) in_range = 0;
    for (int i = 0; i < SW * SH; i++) {
        if (!map[i].ok) continue;
        for (int j = 0; j < SW * SH; j++) {
            if (!map[j].ok) continue;
            if (noise_weight(&map[i], 3) > noise_weight(&map[j], 3) &&
                alloc[i] < alloc[j]) monotone = 0;
        }
    }
    check("alloc_in_range", in_range);
    check("alloc_monotone_in_map", monotone);

    /* 8. allocation consistency. Counts are allocated from the map, the
          scene is re-rendered with variance scaled by base/alloc, and the
          residual is re-estimated. If the allocation follows the map, the
          residual is flat: its spread should collapse well below the spread
          the same scene shows under a uniform allocation. Estimates are
          averaged over several realisations so this measures the map, not the
          estimator's own sampling noise. */
    scene_free(&sc);
    if (!scene_alloc(&sc, SW, SH, 1)) return 2;
    fill_graded(&sc, 0xFEEDu);
    estimate_avg(&sc, NOISE_RADIUS, 3.0, 1, 12, avg, img, scratch);
    noise_allocate(avg, SW * SH, 1, 32, 4, 256, alloc);
    stencil_truth(sc.var, SW, SH, NOISE_RADIUS, truth);
    for (int i = 0; i < SW * SH; i++) use[i] = avg[i].ok;

    /* uniform allocation: spread of the unmodified truth */
    double cv_uniform = cv(truth, use, SW * SH);

    /* map-driven allocation: spread of the scaled truth */
    static double scaled[SW * SH];
    for (int i = 0; i < SW * SH; i++) scaled[i] = sc.var[i] * 32.0 / (double)alloc[i];
    stencil_truth(scaled, SW, SH, NOISE_RADIUS, truth);
    double cv_alloc = cv(truth, use, SW * SH);
    printf("cv_uniform=%.3f cv_alloc_in_truth=%.3f\n", cv_uniform, cv_alloc);
    check("alloc_evens_out_truth", cv_alloc < 0.5 * cv_uniform);

    /* end to end: re-render at the allocated counts and re-estimate */
    for (int i = 0; i < SW * SH; i++) sc.var[i] = scaled[i];
    estimate_avg(&sc, NOISE_RADIUS, 3.0, 1, 12, avg, img, scratch);
    for (int i = 0; i < SW * SH; i++) {
        resid[i] = avg[i].val[0];
        use[i] = avg[i].ok;
    }
    double cv_resid = cv(resid, use, SW * SH);
    check("alloc_flattens_residual", cv_resid < 0.5 * cv_uniform);
    printf("cv_residual=%.3f\n", cv_resid);

    /* 9. noise map and residual map, both from averaged estimates of the same
          graded scene so the before/after contrast reflects the allocation. */
    {
        static NoisePixel before[SW * SH];
        Scene g;
        if (!scene_alloc(&g, SW, SH, 1)) return 2;
        fill_graded(&g, 0xFEEDu);
        estimate_avg(&g, NOISE_RADIUS, 3.0, 1, 12, before, img, scratch);
        noise_allocate(before, SW * SH, 1, 32, 4, 256, alloc);
        for (int i = 0; i < SW * SH; i++)
            g.var[i] *= 32.0 / (double)alloc[i];
        estimate_avg(&g, NOISE_RADIUS, 3.0, 1, 12, map, img, scratch);
        double scale = 0.0;
        for (int i = 0; i < SW * SH; i++)
            if (before[i].ok && before[i].val[0] > scale) scale = before[i].val[0];
        printf("noise_map_before (block 2x2, same scale, '?'=not estimable):\n");
        print_noise_map(before, SW, SH, 0, 2, scale);
        printf("noise_map_after  (block 2x2, same scale):\n");
        print_noise_map(map, SW, SH, 0, 2, scale);
        scene_free(&g);
    }

    scene_free(&sc);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
