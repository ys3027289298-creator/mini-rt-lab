#ifndef MINI_RT_SCENE_H
#define MINI_RT_SCENE_H

/*
 * The single scene shared by all three acceptance artifacts:
 *   sample/   -> histogram of the sampling stream
 *   material/ -> variance of the material stream
 *   bvh/      -> image hash of the geometry stream
 *
 * Every artifact embeds SCENE_ID and SCENE_SEED so the harness can prove the
 * three datasets were produced from the same scene. Stream ids keep the three
 * draw sequences disjoint while deriving from the one seed.
 */

#define SCENE_ID "mini-rt-lab:scene-v1"
#define SCENE_SEED UINT64_C(0x9E3779B97F4A7C15)

/* Image / integration resolution shared by the verifiers. */
#define SCENE_WIDTH  64
#define SCENE_HEIGHT 48
#define SCENE_SPP    256

/* Geometry: 12 axis-aligned boxes, one row along x, deterministic layout. */
#define SCENE_NPRIM 12

/* Per-consumer stream selectors for rng_seed(). */
#define STREAM_SAMPLE   UINT64_C(1)
#define STREAM_MATERIAL UINT64_C(2)
#define STREAM_GEOMETRY UINT64_C(3)

typedef struct {
    float min[3];
    float max[3];
} ScenePrim;

static const ScenePrim SCENE_PRIMS[SCENE_NPRIM] = {
    { { -3.30f, 0.0f,  0.00f }, { -2.90f, 1.0f, 0.50f } },
    { { -2.70f, 0.0f, -0.10f }, { -2.30f, 1.0f, 0.40f } },
    { { -2.10f, 0.0f, -0.20f }, { -1.70f, 1.0f, 0.30f } },
    { { -1.50f, 0.0f, -0.30f }, { -1.10f, 1.0f, 0.20f } },
    { { -0.90f, 0.0f,  0.00f }, { -0.50f, 1.0f, 0.50f } },
    { { -0.30f, 0.0f, -0.10f }, {  0.10f, 1.0f, 0.40f } },
    { {  0.30f, 0.0f, -0.20f }, {  0.70f, 1.0f, 0.30f } },
    { {  0.90f, 0.0f, -0.30f }, {  1.30f, 1.0f, 0.20f } },
    { {  1.50f, 0.0f,  0.00f }, {  1.90f, 1.0f, 0.50f } },
    { {  2.10f, 0.0f, -0.10f }, {  2.50f, 1.0f, 0.40f } },
    { {  2.70f, 0.0f, -0.20f }, {  3.10f, 1.0f, 0.30f } },
    { {  3.30f, 0.0f, -0.30f }, {  3.70f, 1.0f, 0.20f } }
};

/* Pinhole camera: origin, view direction -z, half-FOV 30 degrees. */
static const float SCENE_ORIGIN[3] = { 0.0f, 2.0f, 6.0f };
#define SCENE_TMAX 20.0f

#endif /* MINI_RT_SCENE_H */
