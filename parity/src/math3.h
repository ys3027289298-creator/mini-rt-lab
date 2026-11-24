/* math3: float Vec3 / Ray helpers shared by every back end. */
#ifndef PARITY_MATH3_H
#define PARITY_MATH3_H

#include <math.h>

#define PI 3.14159265358979323846f
#define INV_PI 0.31830988618379067154f

typedef struct { float x, y, z; } V3;
typedef struct { V3 o, d; } Ray;

static inline V3 v3(float x, float y, float z) { V3 r = { x, y, z }; return r; }
static inline V3 vadd(V3 a, V3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline V3 vsub(V3 a, V3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline V3 vmul(V3 a, V3 b) { return v3(a.x * b.x, a.y * b.y, a.z * b.z); }
static inline V3 vscl(V3 a, float s) { return v3(a.x * s, a.y * s, a.z * s); }
static inline float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline V3 vcross(V3 a, V3 b) {
    return v3(a.y * b.z - a.z * b.y,
              a.z * b.x - a.x * b.z,
              a.x * b.y - a.y * b.x);
}
static inline float vlen2(V3 a) { return vdot(a, a); }
static inline float vlen(V3 a) { return sqrtf(vlen2(a)); }
static inline V3 vnorm(V3 a) {
    float l = vlen(a);
    return l > 1e-12f ? vscl(a, 1.0f / l) : v3(0, 0, 0);
}
static inline V3 vlerp(V3 a, V3 b, float t) { return vadd(a, vscl(vsub(b, a), t)); }
static inline V3 vmax(V3 a, V3 b) {
    return v3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z));
}
static inline V3 vmin(V3 a, V3 b) {
    return v3(fminf(a.x, b.x), fminf(a.y, b.y), fminf(a.z, b.z));
}
static inline V3 reflect_v(V3 d, V3 n) { return vsub(d, vscl(n, 2.0f * vdot(d, n))); }

static inline float clampf(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}
static inline V3 vclamp(V3 a, float lo, float hi) {
    return v3(clampf(a.x, lo, hi), clampf(a.y, lo, hi), clampf(a.z, lo, hi));
}
static inline float lerpf(float a, float b, float t) { return a + (b - a) * t; }

static inline float rng_next(unsigned* state) {
    /* xorshift32 */
    unsigned x = *state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    *state = x;
    return (x >> 8) * (1.0f / 16777216.0f);
}

#endif
