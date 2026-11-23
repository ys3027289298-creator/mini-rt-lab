/*
 * rt.c -- minimal path tracer used as the observation oracle for mutation
 * testing.  It contains NO copy of the modules under test: all intersection
 * and BRDF work goes through the rt_bvh_* / rt_* entry points implemented by
 * bvh_shim.c and mat_shim.c, which #include the real bvh/src/bvh.c and
 * material/src/material.c.  A mutant is observed by compiling this renderer
 * against a patched copy of those sources.
 *
 * Usage:  rt <scene> <width> <height> <spp> <seed>
 * Output: stdout, "RT1\0" | u32 w | u32 h | u32 spp | u32 fnv1a(pixels) |
 *         w*h*3 raw RGB8, row major, top row first.
 * Exit 0 if the frame was produced, 2 on bad usage/allocation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include "iface.h"

#define MAX_BOXES 24
#define MAX_SPHERES 4
#define MAX_LIGHTS 1
#define BOUNCES 3

/* ------------------------------------------------------------------ vectors */
typedef struct { float x, y, z; } Vec;
static Vec v3(float x, float y, float z) { Vec v = { x, y, z }; return v; }
static Vec vsub(Vec a, Vec b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static Vec vmul(Vec a, float s) { return v3(a.x * s, a.y * s, a.z * s); }
static Vec vadd(Vec a, Vec b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static float vdot(Vec a, Vec b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static Vec vnorm(Vec a) { float o[3]; rt_norm3((const float[]){ a.x, a.y, a.z }, o); return v3(o[0], o[1], o[2]); }

/* ---------------------------------------------------------------- scene data */
typedef struct {
    float lo[3], hi[3];     /* BVH primitive bounds */
    float albedo, rough;
} Box;

typedef struct { Vec c; float r, albedo, rough; } Sphere;

typedef struct {
    const char* name;
    Box box[MAX_BOXES];
    int nbox;
    Sphere sph[MAX_SPHERES];
    int nsph;
    Vec cam;
    Vec light;              /* direction *towards* the light */
    float fov;
} Scene;

static Scene S;
static int g_brute = 0;   /* RT_BRUTE=1: bypass traversal, test every box */

static void add_box(Scene* s, float x0, float y0, float z0, float x1, float y1, float z1,
                    float albedo, float rough) {
    Box* b = &s->box[s->nbox++];
    b->lo[0] = x0; b->lo[1] = y0; b->lo[2] = z0;
    b->hi[0] = x1; b->hi[1] = y1; b->hi[2] = z1;
    b->albedo = albedo; b->rough = rough;
}
static void add_sphere(Scene* s, float x, float y, float z, float r,
                       float albedo, float rough) {
    Sphere* p = &s->sph[s->nsph++];
    p->c = v3(x, y, z); p->r = r; p->albedo = albedo; p->rough = rough;
}

/*
 * Every scene is built so that no two surfaces are exactly coplanar or
 * exactly tangent: a tie in hit distance would let traversal order (and thus
 * tree shape) pick a different primitive, which would make tree-shape-only
 * mutants look visible.  The four scenes together exercise:
 *
 * "grid":    wall + floor + 4x3 boxes at staggered positions/depths + spheres.
 *            Splitting, occlusion, primitive ids.
 * "rough":   roughness ramp 0.03/0.40/0.85 on spheres against a wall.
 *            GGX shape, Lambert/pdf estimator, sampler.
 * "occlude": ten narrow boxes in a row along the view axis, overlapping in
 *            screen space.  Nearest-hit ordering and split partitioning.
 * "grazing": a slab whose front face is 0.005 in front of the camera, so
 *            primary rays start essentially on it (slab entry parameter), plus
 *            thin posts and a floor near-parallel to the rays.
 * "warm":    dark warm-albedo ramp, so small albedo errors are not clipped
 *            away by the 8-bit transfer curve.
 */
static int scene_build(Scene* s, const char* name) {
    memset(s, 0, sizeof(*s));
    s->name = name;
    s->cam = v3(0.0f, 0.35f, -6.0f);
    s->light = vnorm(v3(0.35f, 0.85f, -0.4f));
    s->fov = 0.72f;

    if (!strcmp(name, "grid")) {
        add_box(s, -5.0f, -1.2f, 6.0f, 5.0f, 5.0f, 6.4f, 0.42f, 0.6f);   /* back wall */
        add_box(s, -5.0f, -1.2f, -1.0f, 5.0f, -0.90f, 6.4f, 0.30f, 0.7f); /* floor */
        int k = 0;
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 4; col++) {
                float x = -2.1f + 1.4f * col;
                float h = 0.8f + 0.25f * row;
                float z = 1.2f + 0.45f * row + 0.2f * col;
                float hue = (float)k / 12.0f;
                add_box(s, x - 0.35f, -0.87f, z, x + 0.35f, -0.87f + h, z + 0.5f,
                        0.15f + 0.6f * hue, 0.08f + 0.12f * row + 0.05f * col);
                k++;
            }
        }
        add_sphere(s, 2.15f, -0.33f, 0.55f, 0.55f, 0.75f, 0.25f);
        add_sphere(s, -2.4f, -0.52f, -0.8f, 0.35f, 0.55f, 0.05f);
        return 1;
    }
    if (!strcmp(name, "warm")) s->light = vnorm(v3(-0.2f, 0.45f, -0.8f));
    if (!strcmp(name, "rough")) {
        add_box(s, -7.0f, -2.0f, 6.0f, 7.0f, 5.0f, 6.5f, 0.35f, 0.9f);
        add_box(s, -7.0f, -2.0f, -2.0f, 7.0f, -1.60f, 6.2f, 0.20f, 0.95f);
        add_sphere(s, -2.05f, -0.52f, 1.9f, 1.05f, 0.72f, 0.03f);
        add_sphere(s, 0.0f, -0.52f, 1.9f, 1.05f, 0.72f, 0.40f);
        add_sphere(s, 2.05f, -0.52f, 1.9f, 1.05f, 0.72f, 0.85f);
        add_box(s, -1.2f, -1.55f, 1.0f, 1.2f, -0.50f, 1.35f, 0.30f, 0.6f);
        return 1;
    }
    if (!strcmp(name, "occlude")) {
        add_box(s, -6.0f, -2.0f, 8.0f, 6.0f, 5.0f, 8.4f, 0.4f, 0.7f);
        add_box(s, -3.0f, -1.1f, -1.5f, 3.0f, -0.98f, 8.2f, 0.3f, 0.8f);
        for (int i = 0; i < 10; i++) {
            float z = 0.3f + 0.5f * i;
            float x = (i & 1) ? 0.6f : -0.6f;
            add_box(s, x - 0.8f, -0.95f, z, x + 0.8f, 0.5f, z + 0.12f,
                    0.12f + 0.06f * i, 0.1f + 0.05f * i);
        }
        return 1;
    }
    if (!strcmp(name, "grazing")) {
        /* front face 0.005 in front of the camera: primary rays enter the slab
           at t ~ 0.005, which is inside any small positive starting t0 */
        add_box(s, -2.4f, -1.2f, -5.995f, -0.35f, 1.7f, -5.35f, 0.55f, 0.35f);
        add_box(s, -8.0f, -1.9f, -6.0f, 8.0f, -1.80f, 8.0f, 0.28f, 0.85f);
        add_box(s, -8.0f, -1.9f, 7.6f, 8.0f, 5.0f, 8.0f, 0.45f, 0.6f);
        for (int i = 0; i < 6; i++) {
            float x = -2.0f + 0.8f * i;
            add_box(s, x - 0.07f, -1.75f, -1.0f + 0.2f * i, x + 0.07f, 1.2f, -0.82f + 0.2f * i,
                    0.55f, 0.25f);
        }
        add_sphere(s, 1.2f, -1.20f, 2.0f, 0.55f, 0.7f, 0.2f);
        return 1;
    }
    if (!strcmp(name, "warm")) {
        add_box(s, -7.0f, -2.0f, 5.0f, 7.0f, 5.0f, 5.5f, 0.16f, 0.9f);
        add_box(s, -7.0f, -2.0f, -2.0f, 7.0f, -1.70f, 5.2f, 0.10f, 0.95f);
        add_sphere(s, -1.9f, -0.66f, 1.6f, 1.0f, 0.18f, 0.05f);
        add_sphere(s, 0.0f, -0.66f, 1.6f, 1.0f, 0.14f, 0.45f);
        add_sphere(s, 1.9f, -0.66f, 1.6f, 1.0f, 0.20f, 0.9f);
        add_box(s, -3.0f, -1.65f, 3.2f, 3.0f, -0.40f, 3.5f, 0.12f, 0.7f);
        return 1;
    }
    return 0;
}

/* ---------------------------------------------------------------- surfaces */
static void box_normal(const Box* b, Vec p, Vec* n) {
    const float e = 1e-3f;
    if (fabsf(p.x - b->lo[0]) < e) *n = v3(-1, 0, 0);
    else if (fabsf(p.x - b->hi[0]) < e) *n = v3(1, 0, 0);
    else if (fabsf(p.y - b->lo[1]) < e) *n = v3(0, -1, 0);
    else if (fabsf(p.y - b->hi[1]) < e) *n = v3(0, 1, 0);
    else if (fabsf(p.z - b->lo[2]) < e) *n = v3(0, 0, -1);
    else *n = v3(0, 0, 1);
}

typedef struct { int kind; int idx; float t; } Hit;  /* kind 0 box, 1 sphere */

static int sphere_intersect(const Sphere* s, Vec o, Vec d, float tmax, float* tout) {
    Vec oc = vsub(o, s->c);
    float b = vdot(oc, d);
    float c = vdot(oc, oc) - s->r * s->r;
    float disc = b * b - c;
    if (disc < 0.0f) return 0;
    float sq = sqrtf(disc);
    float t = -b - sq;
    if (t < 1e-4f) t = -b + sq;
    if (t < 1e-4f || t > tmax) return 0;
    *tout = t;
    return 1;
}

static int intersect(Vec o, Vec d, float tmax, Hit* h) {
    h->t = tmax; h->kind = -1; h->idx = -1;
    int bi = -1;
    float bt = tmax;
    const float oo[3] = { o.x, o.y, o.z }, dd[3] = { d.x, d.y, d.z };
    int (*hit_fn)(const float*, const float*, float, int*, float*) =
        g_brute ? rt_bvh_hit_brute : rt_bvh_hit;
    if (hit_fn(oo, dd, tmax, &bi, &bt) && bt < h->t) {
        h->t = bt; h->kind = 0; h->idx = bi;
    }
    for (int i = 0; i < S.nsph; i++) {
        float t;
        if (sphere_intersect(&S.sph[i], o, d, h->t, &t) && t < h->t) {
            h->t = t; h->kind = 1; h->idx = i;
        }
    }
    return h->kind >= 0;
}

static Vec sky(Vec d) {
    float y = 0.5f * (d.y + 1.0f);
    if (y < 0.0f) y = 0.0f;
    if (y > 1.0f) y = 1.0f;
    return v3(0.10f + 0.22f * y, 0.16f + 0.26f * y, 0.38f + 0.34f * y);
}

/* ------------------------------------------------------- deterministic RNG */
static uint32_t rng_s, rng_i;
static void rng_seed(uint32_t seed) {
    rng_s = seed + 0x9e3779b9u;
    rng_i = 0xda3e39cbu;
    rng_s = rng_s * 747796405u + 2891336453u;
    uint32_t w = ((rng_s >> ((rng_s >> 28u) + 4u)) ^ rng_s) * 277803737u;
    rng_s = (w >> 22u) ^ w;
}
static uint32_t rng_u32(void) {
    uint32_t old = rng_s;
    rng_s = old * 747796405u + 2891336453u;
    uint32_t w = ((old >> ((old >> 28u) + 4u)) ^ old) * 277803737u;
    return (w >> 22u) ^ w;
}
static float rng_f(void) { return (float)(rng_u32() >> 8) * (1.0f / 16777216.0f); }

/* --------------------------------------------------------------- path tracer */
typedef struct { Vec p, n, albedo; float rough; } Surf;

static int shade_point(Hit h, Vec o, Vec d, Vec* p, Surf* s) {
    p->x = o.x + d.x * h.t; p->y = o.y + d.y * h.t; p->z = o.z + d.z * h.t;
    if (h.kind == 0) {
        const Box* b = &S.box[h.idx];
        box_normal(b, *p, &s->n);
        s->albedo = v3(b->albedo, b->albedo, b->albedo);
        s->rough = b->rough;
    } else {
        const Sphere* sp = &S.sph[h.idx];
        s->n = vnorm(v3((p->x - sp->c.x) / sp->r,
                        (p->y - sp->c.y) / sp->r,
                        (p->z - sp->c.z) / sp->r));
        s->albedo = v3(sp->albedo, sp->albedo, sp->albedo);
        s->rough = sp->rough;
    }
    s->p = *p;
    if (vdot(s->n, d) > 0.0f) s->n = vmul(s->n, -1.0f);
    return 1;
}

static Vec trace(Vec o, Vec d) {
    Vec acc = v3(0, 0, 0), thr = v3(1, 1, 1);
    for (int depth = 0; depth < BOUNCES; depth++) {
        Hit h;
        if (!intersect(o, d, 1e30f, &h)) {
            Vec s = sky(d);
            return vadd(acc, v3(thr.x * s.x, thr.y * s.y, thr.z * s.z));
        }
        Vec p; Surf sf;
        shade_point(h, o, d, &p, &sf);

        /* shadow ray towards the single directional light */
        Vec l = S.light;
        Hit sh;
        float lit;
        if (intersect(vadd(p, vmul(sf.n, 1e-3f)), l, 1e30f, &sh)) lit = 0.0f;
        else lit = vdot(sf.n, l) > 0.0f ? vdot(sf.n, l) : 0.0f;

        /* direct: Lambert diffuse + a small GGX specular lobe */
        float lam = rt_brdf_lambert(sf.albedo.x);
        Vec vdir = vmul(d, -1.0f);
        Vec hv = vnorm(vadd(l, vdir));
        float nh = vdot(sf.n, hv);
        float spec = rt_ggx_ndf(nh, sf.rough) * (sf.rough * sf.rough) * lit * 0.25f;
        float direct = lam * lit;
        if (!(direct == direct) || !(spec == spec) || direct < 0.0f || spec < 0.0f) {
            return v3(1, 1, 1);              /* poison flag: keep it visible */
        }
        acc = vadd(acc, v3(thr.x * (direct + spec),
                           thr.y * (direct + spec),
                           thr.z * (direct + spec)));

        /* cosine-weighted bounce; estimator weight via brdf and pdf */
        float u1 = rng_f(), u2 = rng_f();
        float ls[3];
        rt_cosine_sample(u1, u2, ls);
        Vec lv = v3(ls[0], ls[1], ls[2]);
        float cos_t = lv.z;
        float pdf = rt_cosine_pdf(cos_t);
        if (!(pdf > 0.0f) || !(cos_t == cos_t)) return v3(1, 1, 1);
        float w = lam * cos_t / pdf;
        if (!(w == w) || w < 0.0f) return v3(1, 1, 1);

        /* tangent frame around sf.n (Duff et al.) */
        Vec n = sf.n;
        float sg = n.z >= 0.0f ? 1.0f : -1.0f;
        float aa = -1.0f / (sg + n.z);
        float cc = n.x * n.y * aa;
        Vec tan = v3(1.0f + sg * n.x * n.x * aa, sg * cc, -sg * n.x);
        Vec bit = v3(cc, sg + n.y * n.y * aa, -n.y);
        Vec wd = vnorm(vadd(vadd(vmul(tan, lv.x), vmul(bit, lv.y)), vmul(n, lv.z)));

        thr = vmul(thr, w);
        if (thr.x < 1e-4f && thr.y < 1e-4f && thr.z < 1e-4f) return acc;
        o = vadd(p, vmul(sf.n, 1e-3f));
        d = wd;
    }
    return acc;
}

/* ------------------------------------------------------------------ output */
static uint32_t fnv1a(const unsigned char* p, size_t n, uint32_t h) {
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 16777619u; }
    return h;
}
static unsigned char to8(float c) {
    if (!(c == c)) return 255;
    if (c < 0.0f) c = 0.0f;
    if (c > 1.0f) c = 1.0f;
    float g = powf(c, 1.0f / 2.2f);
    int i = (int)(255.0f * g + 0.5f);
    if (i < 0) i = 0;
    if (i > 255) i = 255;
    return (unsigned char)i;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: %s <grid|rough|occlude> <w> <h> <spp> <seed>\n", argv[0]);
        return 2;
    }
    if (!scene_build(&S, argv[1])) {
        fprintf(stderr, "unknown scene '%s'\n", argv[1]);
        return 2;
    }
    g_brute = getenv("RT_BRUTE") && getenv("RT_BRUTE")[0] == '1';
    int W = atoi(argv[2]), H = atoi(argv[3]), spp = atoi(argv[4]);
    uint32_t seed = (uint32_t)strtoul(argv[5], NULL, 0);
    if (W <= 0 || H <= 0 || spp <= 0 || W > 4096 || H > 4096 || spp > 4096) return 2;

    float* boxes = NULL;
    if (S.nbox > 0) {
        boxes = (float*)malloc(sizeof(float) * 6 * (size_t)S.nbox);
        if (!boxes) return 2;
        for (int i = 0; i < S.nbox; i++) {
            memcpy(boxes + 6 * i, S.box[i].lo, 3 * sizeof(float));
            memcpy(boxes + 6 * i + 3, S.box[i].hi, 3 * sizeof(float));
        }
        if (!rt_bvh_build(boxes, S.nbox)) {
            fprintf(stderr, "bvh build failed (%d prims)\n", S.nbox);
            return 2;
        }
        if (rt_bvh_prim_count() != S.nbox) {
            fprintf(stderr, "bvh prim count %d != %d (mutant dropped primitives)\n",
                    rt_bvh_prim_count(), S.nbox);
            return 3;
        }
    }

    size_t npx = (size_t)W * H;
    unsigned char* img = (unsigned char*)malloc(npx * 3);
    if (!img) return 2;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            double r = 0, g = 0, b = 0;
            for (int s = 0; s < spp; s++) {
                uint32_t pix = (uint32_t)(y * W + x);
                rng_seed(seed ^ (pix * 2654435761u) ^ ((uint32_t)s * 0x9e3779b9u));
                float jx = rng_f(), jy = rng_f();
                float ndcx = ((float)x + jx) / (float)W * 2.0f - 1.0f;
                float ndcy = 1.0f - ((float)y + jy) / (float)H * 2.0f;
                float aspect = (float)W / (float)H;
                Vec d = vnorm(v3(ndcx * S.fov * aspect, ndcy * S.fov, 1.0f));
                Vec c = trace(S.cam, d);
                if (!(c.x == c.x) || !(c.y == c.y) || !(c.z == c.z)) c = v3(1, 1, 1);
                r += c.x; g += c.y; b += c.z;
            }
            float inv = 1.0f / (float)spp;
            float cr = (float)(r * inv), cg = (float)(g * inv), cb = (float)(b * inv);
            size_t o = (size_t)(y * W + x) * 3;
            img[o + 0] = to8(cr);
            img[o + 1] = to8(cg);
            img[o + 2] = to8(cb);
        }
    }

    uint32_t hv = fnv1a(img, npx * 3, 2166136261u);
    uint32_t hdr[4];
    hdr[0] = (uint32_t)W; hdr[1] = (uint32_t)H;
    hdr[2] = (uint32_t)spp; hdr[3] = hv;
    fwrite("RT1\0", 1, 4, stdout);
    fwrite(hdr, sizeof(hdr), 1, stdout);
    fwrite(img, 1, npx * 3, stdout);
    fflush(stdout);
    free(img);
    free(boxes);
    return 0;
}
