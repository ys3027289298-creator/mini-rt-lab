/*
 * Acceptance verifier for the sampling module.
 *
 * Pulls the scene-sized number of draws from the unified sampling stream and
 * bins them into a histogram. Emits:
 *   <outdir>/histogram.csv  per-bin counts with the scene header
 * and runs a chi-square uniformity acceptance check on the spot.
 *
 * All randomness comes from shared/rng.c — no other random source exists.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rng.h"
#include "scene.h"
#include "src/sample.c"

#define NBINS 16

#undef main
int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <outdir>\n", argv[0]);
        return 2;
    }
    const char* outdir = argv[1];
    Rng r;
    rng_seed(&r, SCENE_SEED, STREAM_SAMPLE);

    int draws = SCENE_WIDTH * SCENE_HEIGHT * SCENE_SPP;
    long counts[NBINS];
    memset(counts, 0, sizeof(counts));
    double sum = 0.0;
    for (int i = 0; i < draws; i++) {
        double u = rng_next_double(&r);
        int b = (int)(u * NBINS);
        if (b >= NBINS) b = NBINS - 1;
        counts[b]++;
        sum += u;
    }
    double expected = (double)draws / NBINS;
    double chi2 = 0.0;
    double max_rel = 0.0;
    for (int b = 0; b < NBINS; b++) {
        double d = (double)counts[b] - expected;
        chi2 += d * d / expected;
        double rel = d / expected;
        if (rel < 0) rel = -rel;
        if (rel > max_rel) max_rel = rel;
    }
    double mean = sum / draws;
    int uniform_ok = (chi2 < 40.0) && (max_rel < 0.05);

    char path[512];
    snprintf(path, sizeof(path), "%s/histogram.csv", outdir);
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); return 2; }
    fprintf(f, "# scene_id=%s\n", SCENE_ID);
    fprintf(f, "# seed=%llu\n", (unsigned long long)SCENE_SEED);
    fprintf(f, "# stream=%llu\n", (unsigned long long)STREAM_SAMPLE);
    fprintf(f, "# draws=%d,nbins=%d\n", draws, NBINS);
    fprintf(f, "bin,lo,hi,count\n");
    for (int b = 0; b < NBINS; b++) {
        fprintf(f, "%d,%.6f,%.6f,%ld\n",
                b, (double)b / NBINS, (double)(b + 1) / NBINS, counts[b]);
    }
    fclose(f);

    printf("scene_id=%s\n", SCENE_ID);
    printf("seed=%llu\n", (unsigned long long)SCENE_SEED);
    printf("stream=%llu\n", (unsigned long long)STREAM_SAMPLE);
    printf("draws=%d\n", draws);
    printf("mean=%.9f\n", mean);
    printf("chi2=%.6f\n", chi2);
    printf("max_rel=%.6f\n", max_rel);
    printf("hist_uniform=%d\n", uniform_ok);
    printf("artifact=%s\n", path);
    return uniform_ok ? 0 : 1;
}
