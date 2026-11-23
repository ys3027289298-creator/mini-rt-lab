/*
 * Acceptance verifier for the material module.
 *
 * Draws the scene-sized number of cosine-weighted hemisphere directions from
 * the unified material stream, estimates the sample variance of cos(theta)=z
 * (analytic value 1/18 for cosine-weighted sampling) and emits the data as
 * JSON for the shared scene manifest.
 *
 * All randomness comes from shared/rng.c — no other random source exists.
 */
#include <stdio.h>
#include <math.h>
#include "rng.h"
#include "scene.h"
#include "src/material.c"

#undef main
int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <outdir>\n", argv[0]);
        return 2;
    }
    const char* outdir = argv[1];
    Rng r;
    rng_seed(&r, SCENE_SEED, STREAM_MATERIAL);

    int n = SCENE_WIDTH * SCENE_HEIGHT * SCENE_SPP;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; i++) {
        float u1 = rng_next_float(&r);
        float u2 = rng_next_float(&r);
        V3 v = cosine_sample(u1, u2);
        double z = (double)v.z;
        sum += z;
        sum2 += z * z;
    }
    double mean = sum / n;
    double var = (sum2 / n) - mean * mean;
    double analytic_var = 1.0 / 18.0; /* Var(cos) under pdf = cos/pi */
    double analytic_mean = 2.0 / 3.0;
    int variance_ok = fabs(var - analytic_var) < 0.02 &&
                      fabs(mean - analytic_mean) < 0.02;

    char path[512];
    snprintf(path, sizeof(path), "%s/variance.json", outdir);
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); return 2; }
    fprintf(f, "{\n");
    fprintf(f, "  \"scene_id\": \"%s\",\n", SCENE_ID);
    fprintf(f, "  \"seed\": %llu,\n", (unsigned long long)SCENE_SEED);
    fprintf(f, "  \"stream\": %llu,\n", (unsigned long long)STREAM_MATERIAL);
    fprintf(f, "  \"n\": %d,\n", n);
    fprintf(f, "  \"statistic\": \"cos_theta_sample_variance\",\n");
    fprintf(f, "  \"mean\": %.12f,\n", mean);
    fprintf(f, "  \"variance\": %.12f,\n", var);
    fprintf(f, "  \"analytic_mean\": %.12f,\n", analytic_mean);
    fprintf(f, "  \"analytic_variance\": %.12f,\n", analytic_var);
    fprintf(f, "  \"pass\": %s\n", variance_ok ? "true" : "false");
    fprintf(f, "}\n");
    fclose(f);

    printf("scene_id=%s\n", SCENE_ID);
    printf("seed=%llu\n", (unsigned long long)SCENE_SEED);
    printf("stream=%llu\n", (unsigned long long)STREAM_MATERIAL);
    printf("n=%d\n", n);
    printf("mean=%.9f (analytic %.9f)\n", mean, analytic_mean);
    printf("variance=%.9f (analytic %.9f)\n", var, analytic_var);
    printf("variance_ok=%d\n", variance_ok);
    printf("artifact=%s\n", path);
    return variance_ok ? 0 : 1;
}
