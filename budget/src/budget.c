/* budget: block-wise error statistics for a stochastic renderer.

   Four things, each with a self-check in main():

   1. Error heatmap. Per-block error is written to a PPM, one solid cell
      per block. A cell's intensity is a linear ramp of its error
      normalised by the maximum, so the brightest cell and the numerically
      worst block are the same location by construction. The test does not
      trust that: it reads the PPM back and scans for its brightest pixel,
      then compares that block index against the numeric maximum. A layout
      or packing bug fails the test instead of hiding a mismatch.

   2. Worst-block ranking by quantile. Blocks are ordered by error and the
      error at a quantile boundary is reported (q=0 worst, q=1 best), plus
      the list of blocks in the worst decile.

   3. Budget allocation. Each block's estimate has variance var_k / N_k,
      and total error is sum_k var_k / N_k. Minimising that under
      sum_k N_k = B gives N_k proportional to sqrt(var_k), the usual
      variance-optimal split. It is compared against uniform allocation on
      the same budget and the ratio is printed.

   4. Unmet blocks. After spending the whole budget, blocks still above a
      per-block error target are listed, worst first. This is the case the
      uniform split hides: the budget ran out before every block met the
      target.

   Per-block variance is synthetic (scene_var), shaped like a render: a
   flat floor, a smooth hot region, and a speckle region. Blocks differ
   enough that the worst block is unambiguous. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define W_MAX 256
#define H_MAX 256
#define MAX_BLOCK (W_MAX * H_MAX)

typedef struct { int x, y, w, h; } Block;

typedef struct {
    int w, h;                    /* image size in pixels */
    int bw, bh;                  /* blocks per axis */
    int nb;                      /* number of blocks */
    Block b[MAX_BLOCK];
    double var[MAX_BLOCK];       /* per-block single-sample variance */
    double err[MAX_BLOCK];       /* per-block error at the current N */
} Field;

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static double scene_var(int px, int py, int w, int h) {
    double base = 0.04;
    double cx = w * 0.72, cy = h * 0.38;
    double dx = px - cx, dy = py - cy;
    double sigma = w * 0.10;
    double hot = 0.9 * exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma));
    double speck = ((px * 7 + py * 13) % 11 == 0) ? 0.25 * (double)((px * 3 + py) % 5) : 0.0;
    return base + hot + speck;
}

static int field_init(Field* f, int w, int h, int bs) {
    if (w <= 0 || h <= 0 || w > W_MAX || h > H_MAX) return -1;
    memset(f, 0, sizeof(*f));
    f->w = w; f->h = h;
    f->bw = (w + bs - 1) / bs;
    f->bh = (h + bs - 1) / bs;
    f->nb = f->bw * f->bh;
    if (f->nb > MAX_BLOCK) return -1;
    for (int by = 0; by < f->bh; by++) {
        for (int bx = 0; bx < f->bw; bx++) {
            int i = by * f->bw + bx;
            Block* b = &f->b[i];
            b->x = bx * bs; b->y = by * bs;
            b->w = (b->x + bs > w) ? w - b->x : bs;
            b->h = (b->y + bs > h) ? h - b->y : bs;
            double s = 0.0;
            for (int py = b->y; py < b->y + b->h; py++)
                for (int px = b->x; px < b->x + b->w; px++)
                    s += scene_var(px, py, w, h);
            f->var[i] = s / (double)(b->w * b->h);
        }
    }
    return 0;
}

/* Error of every block for a given per-block sample count: var / N. */
static void field_eval(Field* f, const int* n) {
    for (int i = 0; i < f->nb; i++) {
        int k = n[i] > 0 ? n[i] : 1;
        f->err[i] = f->var[i] / (double)k;
    }
}

static double field_total(const Field* f) {
    double s = 0.0;
    for (int i = 0; i < f->nb; i++) s += f->err[i];
    return s;
}

/* Largest-error block; strict > keeps ties at the first row-major block,
   the same order the heatmap scan uses. */
static int field_worst(const Field* f) {
    int wi = -1; double wv = -1.0;
    for (int i = 0; i < f->nb; i++)
        if (f->err[i] > wv) { wv = f->err[i]; wi = i; }
    return wi;
}

static void alloc_uniform(int* n, int nb, long budget) {
    long per = budget / nb, rem = budget % nb;
    for (int i = 0; i < nb; i++) n[i] = (int)per + (i < rem ? 1 : 0);
}

/* N_k proportional to sqrt(var_k). One sample of floor for every block so
   none is left unestimated, then largest-remainder rounding of the rest. */
static void alloc_optimal(int* n, const double* var, int nb, long budget) {
    if (nb <= 0) return;
    if (budget < nb) { for (int i = 0; i < nb; i++) n[i] = 1; return; }
    double wsum = 0.0;
    for (int i = 0; i < nb; i++) wsum += sqrt(var[i] > 0 ? var[i] : 0);
    if (wsum <= 0) { alloc_uniform(n, nb, budget); return; }
    long room = budget - nb;
    long used = 0;
    for (int i = 0; i < nb; i++) {
        double w = sqrt(var[i] > 0 ? var[i] : 0);
        n[i] = 1 + (int)floor((double)room * w / wsum);
        used += n[i];
    }
    long left = budget - used;          /* 0 .. nb-1 */
    while (left > 0) {
        int bi = -1; double best = -1.0;
        for (int i = 0; i < nb; i++) {
            double ideal = (double)room * sqrt(var[i] > 0 ? var[i] : 0) / wsum;
            double frac = ideal - floor(ideal);
            if (n[i] == 1 + (int)floor(ideal) && frac > best) { best = frac; bi = i; }
        }
        if (bi < 0) break;
        n[bi]++;
        left--;
    }
}

/* Indices ordered by error descending, worst first. */
static int rank_by_error(const Field* f, const double* v, int* out) {
    int idx[MAX_BLOCK];
    for (int i = 0; i < f->nb; i++) idx[i] = i;
    for (int i = 1; i < f->nb; i++) {
        int key = idx[i], j = i - 1;
        while (j >= 0 && v[idx[j]] < v[key]) { idx[j + 1] = idx[j]; j--; }
        idx[j + 1] = key;
    }
    if (out) for (int i = 0; i < f->nb; i++) out[i] = idx[i];
    return f->nb;
}

/* Nearest-rank value at quantile q of the descending error list
   (q = 0 is the worst block, q = 1 the best). */
static double quantile_error(const Field* f, const double* v, double q) {
    int rank[MAX_BLOCK];
    int n = rank_by_error(f, v, rank);
    if (n <= 0) return 0.0;
    if (q < 0) q = 0;
    if (q > 1) q = 1;
    int k = (int)floor(q * (double)(n - 1));
    if (k < 0) k = 0;
    if (k >= n) k = n - 1;
    return v[rank[k]];
}

/* Blocks still above tol once the budget is spent; worst first. */
static int unmet_blocks(const Field* f, double tol, int* out, int cap) {
    int rank[MAX_BLOCK];
    int n = rank_by_error(f, f->err, rank);
    int m = 0;
    for (int i = 0; i < n && m < cap; i++)
        if (f->err[rank[i]] > tol) out[m++] = rank[i];
    return m;
}

/* Per-block error to PPM. One cell per block; intensity is error / max
   mapped linearly to 0..255. */
static int heatmap_write(const Field* f, const char* path, int cell) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return -1;
    fprintf(fp, "P6\n%d %d\n255\n", f->bw * cell, f->bh * cell);
    double mx = 0.0;
    for (int i = 0; i < f->nb; i++) if (f->err[i] > mx) mx = f->err[i];
    for (int by = 0; by < f->bh; by++) {
        for (int yy = 0; yy < cell; yy++) {
            for (int bx = 0; bx < f->bw; bx++) {
                double t = (mx > 0) ? f->err[by * f->bw + bx] / mx : 0.0;
                if (t < 0) t = 0;
                if (t > 1) t = 1;
                unsigned char px[3] = { (unsigned char)lround(t * 255.0),
                                        (unsigned char)lround(t * 255.0),
                                        (unsigned char)lround(t * 255.0) };
                for (int xx = 0; xx < cell; xx++) fwrite(px, 1, 3, fp);
            }
        }
    }
    if (fclose(fp) != 0) return -1;
    return 0;
}

/* Read back a PPM heatmap and return the block index of its brightest
   pixel, or -1 if the file cannot be parsed. */
static int heatmap_peak_block(const char* path, int bw, int cell) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return -1;
    int iw, ih, maxv;
    if (fscanf(fp, "P6 %d %d %d", &iw, &ih, &maxv) != 3) { fclose(fp); return -1; }
    if (fgetc(fp) == EOF) { fclose(fp); return -1; }
    int best = -1, bestv = -1;
    for (int by = 0; by < ih / cell; by++) {
        for (int yy = 0; yy < cell; yy++) {
            for (int bx = 0; bx < bw; bx++) {
                for (int xx = 0; xx < cell; xx++) {
                    unsigned char px[3];
                    if (fread(px, 1, 3, fp) != 3) { fclose(fp); return -1; }
                    if ((int)px[0] > bestv) { bestv = px[0]; best = by * bw + bx; }
                }
            }
        }
    }
    fclose(fp);
    return best;
}

int main(void) {
    const int bs = 16;
    Field f, g;
    int n[MAX_BLOCK];

    check("field_init", field_init(&f, 160, 96, bs) == 0);
    check("block_count", f.nb == 10 * 6);
    check("block_edges",
          f.b[f.bw - 1].x + f.b[f.bw - 1].w == f.w &&
          f.b[(f.bh - 1) * f.bw + (f.bw - 1)].y +
              f.b[(f.bh - 1) * f.bw + (f.bw - 1)].h == f.h);

    /* 64 samples per block, uniformly */
    long budget = 64L * f.nb;
    alloc_uniform(n, f.nb, budget);
    long ub = 0;
    for (int i = 0; i < f.nb; i++) ub += n[i];
    check("uniform_uses_budget", ub == budget);
    field_eval(&f, n);

    int worst = field_worst(&f);
    check("worst_valid", worst >= 0 && worst < f.nb);

    /* heatmap: peak cell must equal the numeric worst block */
    const char* heat = "error_heatmap_uniform.ppm";
    check("heatmap_write", heatmap_write(&f, heat, 8) == 0);
    int peak = heatmap_peak_block(heat, f.bw, 8);
    check("peak_matches_worst", peak == worst);

    /* quantile ranking */
    int rank[MAX_BLOCK];
    check("rank_count", rank_by_error(&f, f.err, rank) == f.nb);
    check("rank_sorted", f.err[rank[0]] >= f.err[rank[f.nb - 1]]);
    check("rank_top_is_worst", rank[0] == worst);
    double q10 = quantile_error(&f, f.err, 0.10);
    double q50 = quantile_error(&f, f.err, 0.50);
    double q90 = quantile_error(&f, f.err, 0.90);
    check("quantile_monotone", q10 >= q50 && q50 >= q90);
    check("quantile_worst_ge_top", f.err[worst] >= q10 - 1e-12);

    /* worst decile, by the same ranking */
    int decile = (f.nb + 9) / 10;
    printf("worst_decile_count=%d\n", decile);
    for (int i = 0; i < decile; i++)
        printf("worst[%d]=%d,%d err=%.9f\n", i,
               rank[i] % f.bw, rank[i] / f.bw, f.err[rank[i]]);

    double total_uniform = field_total(&f);
    int worst_u = worst;
    (void)worst_u;

    /* variance-optimal allocation on the same budget */
    int no[MAX_BLOCK];
    alloc_optimal(no, f.var, f.nb, budget);
    long sum = 0;
    int minn = 1 << 30, maxn = 0;
    for (int i = 0; i < f.nb; i++) {
        sum += no[i];
        if (no[i] < minn) minn = no[i];
        if (no[i] > maxn) maxn = no[i];
    }
    check("alloc_uses_budget", sum == budget);
    check("alloc_nonempty", minn >= 1);
    check("alloc_varies", maxn > minn);

    memcpy(&g, &f, sizeof(g));
    field_eval(&g, no);
    double total_opt = field_total(&g);
    int worst_o = field_worst(&g);
    check("alloc_reduces_total", total_opt < total_uniform - 1e-12);

    printf("total_uniform=%.9f\n", total_uniform);
    printf("total_opt=%.9f\n", total_opt);
    printf("opt_over_uniform=%.6f\n", total_opt / total_uniform);
    printf("worst_block_uniform=%d,%d\n", worst_u % f.bw, worst_u / f.bw);
    printf("worst_block_opt=%d,%d\n", worst_o % g.bw, worst_o / g.bw);
    printf("q10=%.9f\nq50=%.9f\nq90=%.9f\n", q10, q50, q90);

    /* heatmap of the allocated frame; peak must track its own worst */
    const char* heat2 = "error_heatmap_opt.ppm";
    check("heatmap_write_opt", heatmap_write(&g, heat2, 8) == 0);
    int peak2 = heatmap_peak_block(heat2, g.bw, 8);
    check("peak_matches_worst_opt", peak2 == worst_o);

    /* budget exhausted with unmet blocks: target is the error level of
       the 10th-percentile block, so the worst ~10% cannot meet it */
    double tol = quantile_error(&g, g.err, 0.10);
    int got[MAX_BLOCK];
    int m = unmet_blocks(&g, tol, got, MAX_BLOCK);
    printf("unmet_tol=%.9f\n", tol);
    printf("unmet_count=%d\n", m);
    for (int i = 0; i < m; i++)
        printf("unmet[%d]=%d,%d err=%.9f n=%d\n", i,
               got[i] % g.bw, got[i] / g.bw, g.err[got[i]], no[got[i]]);
    check("unmet_nonempty", m > 0);
    check("unmet_worst_first", m > 0 && got[0] == worst_o);
    int sorted_ok = 1, all_over = 1;
    for (int i = 0; i < m; i++) {
        if (g.err[got[i]] <= tol) all_over = 0;
        if (i && g.err[got[i - 1]] < g.err[got[i]]) sorted_ok = 0;
    }
    check("unmet_all_over_tol", all_over == 1);
    check("unmet_sorted", sorted_ok == 1);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
