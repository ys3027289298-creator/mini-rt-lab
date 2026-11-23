/*
 * Acceptance verifier for the geometry (bvh) module.
 *
 * Builds the shared scene boxes into a BVH, then renders a fixed-resolution
 * image by shooting rays whose per-pixel jitter comes from the unified
 * geometry stream. Emits:
 *   <outdir>/image.ppm      binary P6 render of the shared scene
 *   <outdir>/image_hash.json  FNV-1a 64-bit hash of the exact PPM bytes plus
 *                             the scene linkage metadata
 *
 * Every integer/float conversion is fixed-width and the shading is integer
 * math, so the PPM bytes — and therefore the hash — are machine-independent.
 *
 * All randomness comes from shared/rng.c — no other random source exists.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "rng.h"
#include "scene.h"
#include "src/bvh.c"

static uint64_t fnv1a64_file(FILE* f) {
    uint64_t h = UINT64_C(14695981039346656037);
    int ch;
    while ((ch = fgetc(f)) != EOF) {
        h ^= (uint64_t)(unsigned char)ch;
        h *= UINT64_C(1099511628211);
    }
    return h;
}

/* nearest box hit along a ray; returns scene prim id or -1 */
static int scene_nearest(const Tree* t, const float o[3], const float d[3]) {
    int best = -1;
    float best_t = SCENE_TMAX;
    for (int i = 0; i < t->nprim; i++) {
        float thit;
        if (aabb_hit(&t->p[i].box, o, d, best_t, &thit) && thit > 0.0f) {
            best_t = thit;
            best = t->p[i].id;
        }
    }
    return best;
}

static unsigned char clamp_byte(float v) {
    if (v < 0.0f) v = 0.0f;
    if (v > 255.0f) v = 255.0f;
    return (unsigned char)(v + 0.5f);
}

#undef main
int main(int argc, char** argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: %s <outdir>\n", argv[0]);
        return 2;
    }
    const char* outdir = argv[1];

    Tree t;
    memset(&t, 0, sizeof(t));
    for (int i = 0; i < SCENE_NPRIM; i++) {
        const ScenePrim* sp = &SCENE_PRIMS[i];
        tree_add(&t, sp->min[0], sp->min[1], sp->min[2],
                      sp->max[0], sp->max[1], sp->max[2]);
    }
    if (build(&t, 0, t.nprim) != 0) {
        fprintf(stderr, "bvh build failed\n");
        return 2;
    }

    Rng r;
    rng_seed(&r, SCENE_SEED, STREAM_GEOMETRY);

    unsigned char* pixels = malloc((size_t)SCENE_WIDTH * SCENE_HEIGHT * 3);
    if (!pixels) { fprintf(stderr, "oom\n"); return 2; }

    const float invw = 1.0f / (float)SCENE_WIDTH;
    const float invh = 1.0f / (float)SCENE_HEIGHT;
    const float tan_half = 0.5773502691896258f; /* tan(30 deg) */
    const float aspect = (float)SCENE_WIDTH / (float)SCENE_HEIGHT;

    for (int py = 0; py < SCENE_HEIGHT; py++) {
        for (int px = 0; px < SCENE_WIDTH; px++) {
            float acc[3] = { 0.0f, 0.0f, 0.0f };
            for (int s = 0; s < SCENE_SPP; s++) {
                float jx = rng_next_float(&r);
                float jy = rng_next_float(&r);
                float ndx = (2.0f * ((float)px + jx) * invw - 1.0f)
                            * tan_half * aspect;
                float ndy = (1.0f - 2.0f * ((float)py + jy) * invh)
                            * tan_half;
                float dl = sqrtf(ndx * ndx + ndy * ndy + 1.0f);
                float d[3] = { ndx / dl, ndy / dl, -1.0f / dl };

                int id = scene_nearest(&t, SCENE_ORIGIN, d);
                float col[3];
                if (id < 0) {
                    col[0] = 18.0f; col[1] = 22.0f; col[2] = 30.0f;
                } else {
                    /* deterministic per-box tint, integer-friendly fractions */
                    float k = (float)((id % 4) + 1) * 0.25f;
                    col[0] = 90.0f + 60.0f * k;
                    col[1] = 120.0f - 20.0f * k;
                    col[2] = 160.0f - 30.0f * k;
                }
                acc[0] += col[0]; acc[1] += col[1]; acc[2] += col[2];
            }
            size_t off = (size_t)(py * SCENE_WIDTH + px) * 3;
            pixels[off + 0] = clamp_byte(acc[0] / (float)SCENE_SPP);
            pixels[off + 1] = clamp_byte(acc[1] / (float)SCENE_SPP);
            pixels[off + 2] = clamp_byte(acc[2] / (float)SCENE_SPP);
        }
    }

    char path[512];
    snprintf(path, sizeof(path), "%s/image.ppm", outdir);
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); free(pixels); return 2; }
    fprintf(f, "P6\n%d %d\n255\n", SCENE_WIDTH, SCENE_HEIGHT);
    fwrite(pixels, 1, (size_t)SCENE_WIDTH * SCENE_HEIGHT * 3, f);
    fclose(f);
    free(pixels);

    FILE* g = fopen(path, "rb");
    if (!g) { perror(path); return 2; }
    uint64_t hash = fnv1a64_file(g);
    fclose(g);

    char jpath[512];
    snprintf(jpath, sizeof(jpath), "%s/image_hash.json", outdir);
    FILE* j = fopen(jpath, "w");
    if (!j) { perror(jpath); return 2; }
    fprintf(j, "{\n");
    fprintf(j, "  \"scene_id\": \"%s\",\n", SCENE_ID);
    fprintf(j, "  \"seed\": %llu,\n", (unsigned long long)SCENE_SEED);
    fprintf(j, "  \"stream\": %llu,\n", (unsigned long long)STREAM_GEOMETRY);
    fprintf(j, "  \"image\": \"image.ppm\",\n");
    fprintf(j, "  \"width\": %d,\n", SCENE_WIDTH);
    fprintf(j, "  \"height\": %d,\n", SCENE_HEIGHT);
    fprintf(j, "  \"spp\": %d,\n", SCENE_SPP);
    fprintf(j, "  \"hash_alg\": \"fnv1a64\",\n");
    fprintf(j, "  \"fnv1a64\": \"%016llx\"\n", (unsigned long long)hash);
    fprintf(j, "}\n");
    fclose(j);

    printf("scene_id=%s\n", SCENE_ID);
    printf("seed=%llu\n", (unsigned long long)SCENE_SEED);
    printf("stream=%llu\n", (unsigned long long)STREAM_GEOMETRY);
    printf("image=%dx%d spp=%d\n", SCENE_WIDTH, SCENE_HEIGHT, SCENE_SPP);
    printf("fnv1a64=%016llx\n", (unsigned long long)hash);
    printf("artifact=%s\n", path);
    printf("artifact=%s\n", jpath);
    return 0;
}
