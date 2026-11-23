#ifndef PERF_RENDER_H
#define PERF_RENDER_H

/*
 * The three rendering paths under a performance budget.
 *
 * Each path owns one image buffer of width*height RGB bytes.  Its scheduled
 * work is exactly width*height*spp units (one ray, one stratified draw, or one
 * material sample per unit); the harness checks that the measured workload is
 * the scheduled workload, so a path that quietly does extra work fails even
 * if the machine is fast.
 *
 * All three paths draw from the same scene seed and use disjoint streams, in
 * the same spirit as the repository's shared-scene contract.
 */
#include <stdint.h>

#define PERF_SCENE_ID "mini-rt-lab:perf-scene-v1"
#define PERF_SCENE_SEED UINT64_C(0x9E3779B97F4A7C15)

#define PERF_STREAM_BVH      UINT64_C(1)
#define PERF_STREAM_SAMPLE   UINT64_C(2)
#define PERF_STREAM_MATERIAL UINT64_C(3)

/* Timed repeats per path; the CLI refuses anything outside this range. */
#define PERF_MAX_REPEATS 16

typedef struct {
    int width;
    int height;
    int spp;        /* rays / draws / samples per pixel, per path */
    uint64_t seed;
    uint64_t stream;
} PerfRenderParams;

typedef void (*PerfRenderFn)(const PerfRenderParams* p, unsigned char* rgb);

/* bvh/: pinhole rays through a median-split AABB tree built from 12 boxes. */
void perf_bvh_render(const PerfRenderParams* p, unsigned char* rgb);

/* sample/: per-pixel stratified estimate f(u)=u^k plus a Halton spatial term. */
void perf_sample_render(const PerfRenderParams* p, unsigned char* rgb);

/* material/: GGX NDF evaluated from cosine-weighted hemisphere samples. */
void perf_material_render(const PerfRenderParams* p, unsigned char* rgb);

#endif /* PERF_RENDER_H */
