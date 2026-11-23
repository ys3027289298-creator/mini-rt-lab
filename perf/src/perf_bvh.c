/*
 * bvh rendering path for the performance harness.
 *
 * The module source is compiled in UNMODIFIED (only its self-test entry point
 * is renamed); this file adds the workload around it.  Rays are jittered by a
 * deterministic PCG stream and intersected with the module's own AABB test, so
 * the measured path is the repository's geometry code, not a reimplementation.
 *
 * The work counter is incremented once per traced ray.  The harness requires
 * it to equal width*height*spp, which makes the work budget enforceable
 * independently of how fast the machine happens to be.
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define main bvh_module_selftest_main
#include "../../bvh/src/bvh.c"
#undef main

#include "../include/perf_render.h"
#include "../include/perf_rng.h"

/* 12 boxes in one row, camera at (0,2,6) looking down -z. */
#define PERF_NPRIM 12
static const float kBoxMin[PERF_NPRIM][3] = {
    { -3.30f, 0.0f,  0.00f }, { -2.70f, 0.0f, -0.10f },
    { -2.10f, 0.0f, -0.20f }, { -1.50f, 0.0f, -0.30f },
    { -0.90f, 0.0f,  0.00f }, { -0.30f, 0.0f, -0.10f },
    {  0.30f, 0.0f, -0.20f }, {  0.90f, 0.0f, -0.30f },
    {  1.50f, 0.0f,  0.00f }, {  2.10f, 0.0f, -0.10f },
    {  2.70f, 0.0f, -0.20f }, {  3.30f, 0.0f, -0.30f }
};
static const float kBoxMax[PERF_NPRIM][3] = {
    { -2.90f, 1.0f,  0.50f }, { -2.30f, 1.0f,  0.40f },
    { -1.70f, 1.0f,  0.30f }, { -1.10f, 1.0f,  0.20f },
    { -0.50f, 1.0f,  0.50f }, {  0.10f, 1.0f,  0.40f },
    {  0.70f, 1.0f,  0.30f }, {  1.30f, 1.0f,  0.20f },
    {  1.90f, 1.0f,  0.50f }, {  2.50f, 1.0f,  0.40f },
    {  3.10f, 1.0f,  0.30f }, {  3.70f, 1.0f,  0.20f }
};

static uint64_t g_bvh_work;

static unsigned char clamp_byte(float v) {
    if (!(v > 0.0f)) return 0;       /* also catches NaN */
    if (v > 255.0f) return 255;
    return (unsigned char)(v + 0.5f);
}

void perf_bvh_render(const PerfRenderParams* p, unsigned char* rgb) {
    g_bvh_work = 0;
    Tree t;
    memset(&t, 0, sizeof(t));
    for (int i = 0; i < PERF_NPRIM; i++)
        tree_add(&t, kBoxMin[i][0], kBoxMin[i][1], kBoxMin[i][2],
                     kBoxMax[i][0], kBoxMax[i][1], kBoxMax[i][2]);
    if (build(&t, 0, t.nprim) != 0) {
        /* The tree is at capacity only if the scene grew past MAX_PRIM;
           leave the buffer black and let the harness report a failure. */
        memset(rgb, 0, (size_t)p->width * (size_t)p->height * 3u);
        return;
    }

    PerfRng rng;
    perf_rng_seed(&rng, p->seed, p->stream);

    const float invw = 1.0f / (float)p->width;
    const float invh = 1.0f / (float)p->height;
    const float tan_half = 0.5773502691896258f;   /* tan(30 deg) */
    const float aspect = (float)p->width / (float)p->height;

    for (int py = 0; py < p->height; py++) {
        for (int px = 0; px < p->width; px++) {
            float acc[3] = { 0.0f, 0.0f, 0.0f };
            for (int s = 0; s < p->spp; s++) {
                float jx = perf_rng_float(&rng);
                float jy = perf_rng_float(&rng);
                float ndx = (2.0f * ((float)px + jx) * invw - 1.0f)
                            * tan_half * aspect;
                float ndy = (1.0f - 2.0f * ((float)py + jy) * invh) * tan_half;
                float dl = sqrtf(ndx * ndx + ndy * ndy + 1.0f);
                float o[3] = { 0.0f, 2.0f, 6.0f };
                float d[3] = { ndx / dl, ndy / dl, -1.0f / dl };

                int best_id = -1;
                float best_t = 20.0f;
                for (int i = 0; i < t.nprim; i++) {
                    float hit;
                    if (aabb_hit(&t.p[i].box, o, d, best_t, &hit) &&
                        hit > 0.0f && hit < best_t) {
                        best_t = hit;
                        best_id = t.p[i].id;
                    }
                }
                g_bvh_work++;

                float col[3];
                if (best_id < 0) {
                    col[0] = 20.0f; col[1] = 24.0f; col[2] = 32.0f;
                } else {
                    float k = (float)((best_id % 4) + 1) * 0.25f;
                    col[0] = 90.0f + 60.0f * k;
                    col[1] = 120.0f - 20.0f * k;
                    col[2] = 160.0f - 30.0f * k;
                }
                acc[0] += col[0]; acc[1] += col[1]; acc[2] += col[2];
            }
            size_t off = ((size_t)py * (size_t)p->width + (size_t)px) * 3u;
            float inv = 1.0f / (float)p->spp;
            rgb[off + 0] = clamp_byte(acc[0] * inv);
            rgb[off + 1] = clamp_byte(acc[1] * inv);
            rgb[off + 2] = clamp_byte(acc[2] * inv);
        }
    }
}

/* The harness reads the counter through this accessor; it exists so the
   unit's work can be asserted against the budget, not just timed. */
uint64_t perf_bvh_work(void) { return g_bvh_work; }
