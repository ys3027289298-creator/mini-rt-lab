/* diffimg: per-pixel difference between the two LDR outputs.
   Operates on tone-mapped sRGB bytes (what a viewer actually compares),
   emits global stats, a luminance-difference histogram (CSV), an ASCII
   histogram and a blue->red spatial heat map. Optional PFM pair adds
   linear-space RMSE. */
#include "math3.h"
#include "image.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NBINS 32

static float lum(unsigned char r, unsigned char g, unsigned char b) {
    return 0.2126f * r + 0.7152f * g + 0.0722f * b;
}

static int cmp_float(const void* a, const void* b) {
    float x = *(const float*)a, y = *(const float*)b;
    return (x > y) - (x < y);
}

static unsigned char heat_byte(float t) { return (unsigned char)clampf(t * 255.0f + 0.5f, 0, 255); }

static double percentile(const float* sorted, size_t n, double q) {
    return sorted[(size_t)(q * (n - 1))];
}

static double frac_le(const long* hist, size_t n, float t) {
    int upto = (int)(t / 256.0f * NBINS);
    long acc = 0;
    for (int j = 0; j <= upto && j < NBINS; j++) acc += hist[j];
    return (double)acc / n;
}

int main(int argc, char** argv) {
    const char* a_path = "out/rt.ppm";
    const char* b_path = "out/rast.ppm";
    const char* csv_path = "out/diff_hist.csv";
    const char* heat_path = "out/diff_heat.ppm";
    const char* hdr_a = NULL, *hdr_b = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-a") && i + 1 < argc) a_path = argv[++i];
        else if (!strcmp(argv[i], "-b") && i + 1 < argc) b_path = argv[++i];
        else if (!strcmp(argv[i], "-csv") && i + 1 < argc) csv_path = argv[++i];
        else if (!strcmp(argv[i], "-heat") && i + 1 < argc) heat_path = argv[++i];
        else if (!strcmp(argv[i], "-hf") && i + 2 < argc) { hdr_a = argv[++i]; hdr_b = argv[++i]; }
    }

    unsigned char *pa, *pb;
    int wa, ha, wb, hb;
    if (!image_read_ppm(&pa, &wa, &ha, a_path) || !image_read_ppm(&pb, &wb, &hb, b_path)) {
        fprintf(stderr, "cannot read inputs\n"); return 1;
    }
    if (wa != wb || ha != hb) { fprintf(stderr, "size mismatch\n"); return 1; }
    size_t n = (size_t)wa * ha;

    long hist[NBINS] = {0};
    float* dl = (float*)malloc(n * sizeof(float));
    double sum_abs = 0.0, sum2 = 0.0, maxd = 0.0;
    double channel_mae[3] = {0};
    unsigned char* heat = (unsigned char*)malloc(3 * n);

    for (size_t i = 0; i < n; i++) {
        float la = lum(pa[i * 3], pa[i * 3 + 1], pa[i * 3 + 2]);
        float lb = lum(pb[i * 3], pb[i * 3 + 1], pb[i * 3 + 2]);
        float d = la - lb;                 /* signed, in sRGB levels [0,255] */
        float ad = fabsf(d);
        dl[i] = ad;
        sum_abs += ad; sum2 += d * d;
        if (ad > maxd) maxd = ad;
        for (int c = 0; c < 3; c++)
            channel_mae[c] += abs((int)pa[i * 3 + c] - (int)pb[i * 3 + c]);
        int bin = (int)(ad / 256.0f * NBINS);
        if (bin >= NBINS) bin = NBINS - 1;
        hist[bin]++;
        float t = ad / 64.0f;              /* heat saturates at 64/255 */
        heat[i * 3 + 0] = heat_byte(t);
        heat[i * 3 + 1] = heat_byte(1.0f - fabsf(2.0f * t - 1.0f));
        heat[i * 3 + 2] = heat_byte(1.0f - t);
    }

    qsort(dl, n, sizeof(float), cmp_float);
    float mean = (float)(sum_abs / n);
    float rmse = (float)sqrt(sum2 / n);

    FILE* csv = fopen(csv_path, "w");
    fprintf(csv, "bin,lo,hi,count,fraction\n");
    for (int k = 0; k < NBINS; k++)
        fprintf(csv, "%d,%.1f,%.1f,%ld,%.6f\n", k,
                k * 256.0 / NBINS, (k + 1) * 256.0 / NBINS, hist[k], (double)hist[k] / n);
    fclose(csv);
    image_write_ppm(heat, wa, ha, heat_path);

    printf("== per-pixel LDR diff (sRGB levels 0..255, %zux%zu=%zu px) ==\n",
           (size_t)wa, (size_t)ha, n);
    printf("MAE(lum)=%.3f  RMSE(lum)=%.3f  max|dL|=%.2f\n", mean, rmse, maxd);
    printf("MAE channels: R=%.3f G=%.3f B=%.3f\n",
           channel_mae[0] / n, channel_mae[1] / n, channel_mae[2] / n);
    printf("percentiles |dL|: p50=%.2f p90=%.2f p95=%.2f p99=%.2f\n",
           percentile(dl, n, 0.50), percentile(dl, n, 0.90),
           percentile(dl, n, 0.95), percentile(dl, n, 0.99));
    printf("fraction |dL|<=4: %.4f   <=8: %.4f   <=16: %.4f   <=32: %.4f\n",
           frac_le(hist, n, 4), frac_le(hist, n, 8),
           frac_le(hist, n, 16), frac_le(hist, n, 32));

    printf("\nASCII histogram of |dL| (bar length scaled to peak bin):\n");
    long peak = 1;
    for (int k = 0; k < NBINS; k++) if (hist[k] > peak) peak = hist[k];
    for (int k = 0; k < NBINS; k++) {
        int bars = (int)(40.0 * hist[k] / peak);
        printf("[%5.1f,%5.1f) %6ld %6.2f%% |", k * 256.0 / NBINS, (k + 1) * 256.0 / NBINS,
               hist[k], 100.0 * hist[k] / n);
        for (int j = 0; j < bars; j++) putchar('#');
        putchar('\n');
    }

    if (hdr_a && hdr_b) {
        HdrImage ha2, hb2;
        if (image_read_pfm(&ha2, hdr_a) && image_read_pfm(&hb2, hdr_b) &&
            ha2.w == hb2.w && ha2.h == hb2.h) {
            double s2 = 0;
            for (size_t i = 0; i < (size_t)ha2.w * ha2.h; i++) {
                V3 e = vsub(ha2.px[i], hb2.px[i]);
                s2 += vdot(e, e) / 3.0;
            }
            printf("\nlinear HDR RMSE(radiance)=%.5f\n", sqrt(s2 / ((size_t)ha2.w * ha2.h)));
            image_free(&ha2); image_free(&hb2);
        }
    }

    free(pa); free(pb); free(dl); free(heat);
    return 0;
}
