/* moire: the three-tier subdivision experiment.

   For each of three subdivision densities it renders the same facade twice --
   once with the fixed stratified grid, once with the low-discrepancy offset --
   measures the dominant stripe spacing of each against the analytic reference,
   and reports lambda*D. A lambda*D that stays constant across the three tiers
   means the offset did not decouple the sample lattice from the panel lattice.

   It also writes, per tier and per sampler, the actual sample positions used
   over a small tile of the wall (sample_positions.csv), plus the measured
   stripe profile (stripe_profile.csv), so the sample distribution can be
   plotted next to the stripes it produced. */
#include "glass.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define W 384
#define H 192
#define SPP 4
#define NTIER 3

static const int tiers[NTIER] = { 708, 368, 249 };

/* Dump the absolute sample positions over a tile of the wall. Absolute (not
   pixel-local) coordinates are what make the phase relation to the panel
   lattice visible: with the fixed grid every pixel repeats the same lattice,
   so relative to the panel edges the samples drift by frac(p) per pixel. */
static void dump_samples(const GlassWall* w, GlassSampleMode mode, const char* path) {
    FILE* f = fopen(path, "w");
    if (!f) return;
    double s[GLASS_MAX_SPP * GLASS_MAX_SPP][2];
    fprintf(f, "px,py,x,y\n");
    int nx = 40, ny = 4;              /* a 40x4-pixel strip */
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nx; x++) {
            int n = glass_pixel_samples(w, x, y, mode, s);
            for (int i = 0; i < n; i++)
                fprintf(f, "%d,%d,%.6f,%.6f\n", x, y,
                        x + s[i][0], y + s[i][1]);
        }
    }
    fclose(f);
}

/* The stripe profile for one render is the column-averaged relative residual
   after the same detrending the spacing measurement uses. Recompute it here
   so the CSV lines up with the reported spacing. */
static void dump_profile(const double* img, const double* ref, const GlassWall* w,
                         const char* path) {
    int Wd = w->width;
    double* prof = (double*)malloc(sizeof(double) * (size_t)Wd);
    double mean = 0.0;
    for (int x = 0; x < Wd; x++) {
        double acc = 0.0;
        for (int y = 0; y < w->height; y++) {
            long i = (long)y * Wd + x;
            acc += (img[i] - ref[i]) / ref[i];
        }
        prof[x] = acc / w->height;
        mean += prof[x];
    }
    mean /= Wd;
    int half = 16;
    FILE* f = fopen(path, "w");
    if (!f) { free(prof); return; }
    fprintf(f, "x,residual\n");
    for (int x = 0; x < Wd; x++) {
        double s = 0.0;
        int c = 0;
        for (int j = -half; j <= half; j++) {
            int k = x + j;
            if (k >= 0 && k < Wd) { s += prof[k] - mean; c++; }
        }
        fprintf(f, "%d,%.8f\n", x, (prof[x] - mean) - s / c);
    }
    fclose(f);
    free(prof);
}

int main(void) {
    double* grid = (double*)malloc(sizeof(double) * W * H);
    double* off = (double*)malloc(sizeof(double) * W * H);
    double* ref = (double*)malloc(sizeof(double) * W * H);
    if (!grid || !off || !ref) { fprintf(stderr, "alloc failed\n"); return 1; }

    printf("wall=%dx%d spp=%d (grid %dx%d) panes per row D = %d %d %d\n",
           W, H, SPP * SPP, SPP, SPP, tiers[0], tiers[1], tiers[2]);
    printf("%-5s %-7s %-7s %-8s %-8s %-9s %-9s %-9s %-9s %-9s\n",
           "D", "pitch", "frac", "pred_lam", "grid_lam", "grid_lamD",
           "band_grid", "band_off", "px_grid", "px_off");

    FILE* summary = fopen("out/summary.csv", "w");
    if (summary)
        fprintf(summary,
                "D,pitch,frac_x_lattice,pred_lam,grid_lam,grid_lamD,off_lam,"
                "off_lamD,band_grid,band_off,pixel_grid,pixel_off,band_ratio\n");

    double grid_lamD[NTIER], off_lamD[NTIER], band_ratio[NTIER];
    for (int t = 0; t < NTIER; t++) {
        GlassWall w = glass_wall(W, H, tiers[t], SPP);
        glass_reference(&w, ref);

        char path[160];
        glass_render(&w, GLASS_SAMPLE_GRID, grid);
        double pg = 0, ag = 0;
        double lg = glass_fringe_spacing(grid, ref, &w, &pg, &ag);
        snprintf(path, sizeof(path), "out/samples_D%d_grid.csv", tiers[t]);
        dump_samples(&w, GLASS_SAMPLE_GRID, path);
        snprintf(path, sizeof(path), "out/stripe_D%d_grid.csv", tiers[t]);
        dump_profile(grid, ref, &w, path);

        glass_render(&w, GLASS_SAMPLE_OFFSET, off);
        double po = 0, ao = 0;
        double lo = glass_fringe_spacing(off, ref, &w, &po, &ao);
        snprintf(path, sizeof(path), "out/samples_D%d_offset.csv", tiers[t]);
        dump_samples(&w, GLASS_SAMPLE_OFFSET, path);
        snprintf(path, sizeof(path), "out/stripe_D%d_offset.csv", tiers[t]);
        dump_profile(off, ref, &w, path);

        /* displacement of the sampling lattice against the panel lattice, as
           a fraction of a panel: frac(S*p), folded into [0, 0.5] */
        double p = glass_panel_pitch(&w);
        double frac = (double)SPP * p - floor((double)SPP * p);
        if (frac > 0.5) frac = 1.0 - frac;
        grid_lamD[t] = lg * tiers[t];
        off_lamD[t] = lo * tiers[t];
        (void)ag; (void)ao;
        double bg = glass_band_rms(grid, ref, &w);
        double bo = glass_band_rms(off, ref, &w);
        double xg = glass_pixel_rms(grid, ref, &w);
        double xo = glass_pixel_rms(off, ref, &w);
        band_ratio[t] = bo > 0 ? bg / bo : 0.0;
        printf("%-5d %-7.4f %-7.4f %-8.2f %-8.2f %-9.0f %-9.5f %-9.5f %-9.5f %-9.5f\n",
               tiers[t], p, frac, glass_predicted_stripe(&w), lg,
               grid_lamD[t], bg, bo, xg, xo);
        if (summary)
            fprintf(summary, "%d,%.6f,%.6f,%.4f,%.4f,%.1f,%g,%g,%.5f,%.5f,%.5f,%.5f,%.2f\n",
                    tiers[t], p, frac, glass_predicted_stripe(&w), lg,
                    grid_lamD[t], lo, off_lamD[t], bg, bo, xg, xo, band_ratio[t]);
    }
    if (summary) fclose(summary);

    /* constancy of lambda*D across tiers: ~1 means the offset did nothing */
    double gmin = grid_lamD[0], gmax = grid_lamD[0];
    double omin = off_lamD[0], omax = off_lamD[0];
    for (int t = 1; t < NTIER; t++) {
        if (grid_lamD[t] < gmin) gmin = grid_lamD[t];
        if (grid_lamD[t] > gmax) gmax = grid_lamD[t];
        if (off_lamD[t] < omin) omin = off_lamD[t];
        if (off_lamD[t] > omax) omax = off_lamD[t];
    }
    double gspread = gmin > 0 ? gmax / gmin : 0;
    double ospread = omin > 0 ? omax / omin : 0;
    printf("\nlambda*D constancy (max/min, near 1 = no effect):\n");
    printf("  grid  : %.3f\n", gspread);
    printf("  offset: %.3f\n", ospread);

    free(grid); free(off); free(ref);
    return 0;
}
