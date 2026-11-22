/*
 * render: multi-threaded software path tracer with a strict determinism
 * contract.
 *
 * The same parameters must produce a bit-identical frame on every run and
 * with every thread count. The contract is:
 *
 *   1. No mutable global state on the deterministic path. The scene is built
 *      once in main() and only read after worker threads start.
 *   2. Each pixel owns a private PCG32 RNG, seeded from (x, y, seed) via a
 *      hash. There is no shared RNG, so the draw sequence for a pixel does
 *      not depend on which thread runs it or when.
 *   3. Each framebuffer pixel is written exactly once, by the single thread
 *      that owns its row band; no read-modify-write sharing, no cross-thread
 *      reduction.
 *   4. Per-pixel sample accumulation follows a fixed order (s = 0..spp-1),
 *      so floating-point summation is independent of scheduling.
 *
 * --rng=global is deliberate fault injection: a process-wide RNG under a
 * mutex, seeded from the clock. Draw order then depends on thread scheduling
 * and the seed changes every run, reproducing the "same frame renders
 * differently twice" bug class. It exists only to exercise the consistency
 * harness and must never be used for real output.
 *
 * Output is a binary RGB Portable FloatMap (PFM, little-endian float32),
 * i.e. the raw framebuffer with no tone mapping or quantization.
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

#define PI 3.14159265358979323846f
#define MAX_DEPTH 1          /* one indirect bounce beyond the primary ray */
#define TMIN 1e-4f

typedef struct { float x, y, z; } V3;

static V3 vadd(V3 a, V3 b) { return (V3){ a.x + b.x, a.y + b.y, a.z + b.z }; }
static V3 vsub(V3 a, V3 b) { return (V3){ a.x - b.x, a.y - b.y, a.z - b.z }; }
static V3 vmul(V3 a, V3 b) { return (V3){ a.x * b.x, a.y * b.y, a.z * b.z }; }
static V3 vscale(V3 a, float k) { return (V3){ a.x * k, a.y * k, a.z * k }; }
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 vcross(V3 a, V3 b) {
    return (V3){ a.y * b.z - a.z * b.y,
                 a.z * b.x - a.x * b.z,
                 a.x * b.y - a.y * b.x };
}
static V3 vnorm(V3 a) {
    float l = sqrtf(vdot(a, a));
    return l > 0.0f ? vscale(a, 1.0f / l) : a;
}

/* ---- PCG32: per-pixel private RNG, no shared instance ------------------ */

typedef struct { uint64_t state; uint64_t inc; } Rng;

static uint32_t pcg_next(Rng* r) {
    uint64_t old = r->state;
    r->state = old * 6364136223846793005ULL + r->inc;
    uint32_t xors = (uint32_t)(((old >> 18u) ^ old) >> 27u);
    uint32_t rot = (uint32_t)(old >> 59u);
    return (xors >> rot) | (xors << ((-rot) & 31u));
}

static void rng_seed(Rng* r, uint64_t s0, uint64_t s1) {
    r->state = 0U;
    r->inc = (s1 << 1u) | 1u;
    pcg_next(r);
    r->state += s0;
    pcg_next(r);
}

static float rng_float(Rng* r) {
    /* 24-bit mantissa float in [0, 1) */
    return (float)(pcg_next(r) >> 8) * (1.0f / 16777216.0f);
}

static uint64_t mix64(uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* ---- Fault injection: one mutable global RNG (see file header) --------- */

static int g_rng_global = 0;
static Rng g_rng;
static pthread_mutex_t g_rng_lock = PTHREAD_MUTEX_INITIALIZER;

/* Single draw entry point used everywhere in the renderer. */
static float rf(Rng* local) {
    if (g_rng_global) {
        pthread_mutex_lock(&g_rng_lock);
        float u = rng_float(&g_rng);
        pthread_mutex_unlock(&g_rng_lock);
        return u;
    }
    return rng_float(local);
}

/* ---- Scene: built once before threads launch, read-only afterwards ----- */

typedef struct { V3 c; float r; V3 albedo; float emit; } Sphere;

#define MAX_SPHERES 8

typedef struct {
    Sphere sp[MAX_SPHERES];
    int nsp;
} Scene;

/* The scene instance is a non-static global pointer set once in main();
 * worker threads only ever read through it. */
static const Scene* g_scene;

static void scene_add(Scene* s, V3 c, float r, V3 albedo, float emit) {
    if (s->nsp >= MAX_SPHERES) { fprintf(stderr, "scene full\n"); exit(2); }
    s->sp[s->nsp++] = (Sphere){ c, r, albedo, emit };
}

static int sphere_hit(const Sphere* s, V3 o, V3 d, float tmin, float* t) {
    V3 oc = vsub(o, s->c);
    float b = vdot(oc, d);
    float cc = vdot(oc, oc) - s->r * s->r;
    float disc = b * b - cc;
    if (disc < 0.0f) return 0;
    float sq = sqrtf(disc);
    float t0 = -b - sq;
    if (t0 > tmin) { *t = t0; return 1; }
    t0 = -b + sq;
    if (t0 > tmin) { *t = t0; return 1; }
    return 0;
}

/* Returns sphere id in *id, or 0 on a miss. */
static int scene_hit(const Scene* s, V3 o, V3 d, float tmin,
                     float* t, int* id) {
    float best = 1e30f;
    int hit = 0;
    for (int i = 0; i < s->nsp; i++) {
        float ti;
        if (sphere_hit(&s->sp[i], o, d, tmin, &ti) && ti < best) {
            best = ti; *id = i; hit = 1;
        }
    }
    if (hit) *t = best;
    return hit;
}

static V3 sky(V3 d) {
    float a = 0.5f * (d.y + 1.0f);
    V3 white = { 1.0f, 1.0f, 1.0f };
    V3 blue = { 0.5f, 0.7f, 1.0f };
    return vadd(vscale(white, 1.0f - a), vscale(blue, a));
}

static V3 cosine_hemisphere(V3 n, Rng* r) {
    float u1 = rf(r);
    float u2 = rf(r);
    float rr = sqrtf(u1);
    float phi = 2.0f * PI * u2;
    V3 helper = fabsf(n.x) > 0.9f ? (V3){ 0, 1, 0 } : (V3){ 1, 0, 0 };
    V3 t = vnorm(vcross(helper, n));
    V3 b = vcross(n, t);
    V3 dir = vadd(vadd(vscale(t, rr * cosf(phi)),
                       vscale(b, rr * sinf(phi))),
                  vscale(n, sqrtf(1.0f - u1)));
    return vnorm(dir);
}

static V3 trace(V3 o, V3 d, int depth, Rng* r) {
    float t;
    int id;
    if (!scene_hit(g_scene, o, d, TMIN, &t, &id)) return sky(d);

    const Sphere* sp = &g_scene->sp[id];
    V3 p = vadd(o, vscale(d, t));
    V3 n = vnorm(vsub(p, sp->c));
    if (sp->emit > 0.0f) return vscale(sp->albedo, sp->emit);
    if (depth <= 0) return (V3){ 0, 0, 0 };

    V3 nd = cosine_hemisphere(n, r);
    V3 ind = trace(p, nd, depth - 1, r);
    return vmul(sp->albedo, ind);
}

/* ---- Camera ------------------------------------------------------------ */

typedef struct {
    V3 origin, lower_left, horiz, vert;
} Camera;

static Camera camera_make(int w, int h) {
    Camera cam;
    cam.origin = (V3){ 0.0f, 1.0f, 3.0f };
    V3 look = vnorm((V3){ 0.0f, -0.15f, -1.0f });
    V3 up = (V3){ 0, 1, 0 };
    V3 right = vnorm(vcross(look, up));
    V3 true_up = vcross(right, look);
    float half_h = tanf(0.5f * 60.0f * PI / 180.0f);
    float half_w = half_h * ((float)w / (float)h);
    V3 center = vadd(cam.origin, look);
    cam.horiz = vscale(right, 2.0f * half_w);
    cam.vert = vscale(true_up, 2.0f * half_h);
    cam.lower_left = vsub(vsub(center, vscale(cam.horiz, 0.5f)),
                          vscale(cam.vert, 0.5f));
    return cam;
}

static V3 camera_ray(const Camera* cam, int x, int y, int w, int h, Rng* r) {
    /* sub-pixel jitter: fixed two draws per sample, always in this order */
    float jx = rf(r);
    float jy = rf(r);
    float u = ((float)x + jx) / (float)w;
    float v = 1.0f - ((float)y + jy) / (float)h; /* y=0 is the top row */
    V3 target = vadd(cam->lower_left,
              vadd(vscale(cam->horiz, u), vscale(cam->vert, v)));
    return vnorm(vsub(target, cam->origin));
}

/* ---- Rendering --------------------------------------------------------- */

typedef struct {
    const Camera* cam;
    float* fb;
    int w, h, spp, seed;
    int y0, y1;
} Job;

static void* worker(void* arg) {
    Job* job = arg;
    const int w = job->w, h = job->h, spp = job->spp;
    for (int y = job->y0; y < job->y1; y++) {
        for (int x = 0; x < w; x++) {
            Rng local;
            uint64_t idx = (uint64_t)y * (uint64_t)w + (uint64_t)x;
            rng_seed(&local,
                     mix64(idx + (uint64_t)job->seed),
                     mix64(~idx ^ (uint64_t)job->seed));

            V3 acc = { 0, 0, 0 };
            for (int s = 0; s < spp; s++) {
                V3 d = camera_ray(job->cam, x, y, w, h, &local);
                V3 c = trace(job->cam->origin, d, MAX_DEPTH, &local);
                acc = vadd(acc, c);
            }
            float* px = &job->fb[((size_t)y * w + x) * 3];
            px[0] = acc.x / (float)spp;
            px[1] = acc.y / (float)spp;
            px[2] = acc.z / (float)spp;
        }
    }
    return NULL;
}

/* PFM: "PF" RGB float32, rows stored bottom-to-top, negative scale means
 * little-endian (the native endianness of the target platform). */
static int write_pfm(const char* path, int w, int h, const float* fb) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return 1; }
    if (fprintf(f, "PF\n%d %d\n-1.0\n", w, h) < 0) { fclose(f); return 1; }
    for (int y = h - 1; y >= 0; y--) {
        if (fwrite(&fb[(size_t)y * w * 3], sizeof(float),
                   (size_t)w * 3, f) != (size_t)w * 3) {
            perror(path);
            fclose(f);
            return 1;
        }
    }
    if (fclose(f) != 0) { perror(path); return 1; }
    return 0;
}

static void usage(const char* prog) {
    fprintf(stderr,
        "usage: %s -w WIDTH -h HEIGHT -s SPP -t THREADS [-o OUT.pfm]\n"
        "          [--rng det|global] [--seed N]\n", prog);
}

int main(int argc, char** argv) {
    int w = 0, h = 0, spp = 0, threads = 0, seed = 1;
    const char* out = NULL;
    const char* rng_mode = "det";

    static struct option long_opts[] = {
        {"rng", required_argument, 0, 1000},
        {"seed", required_argument, 0, 1001},
        {0, 0, 0, 0}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "w:h:s:t:o:", long_opts, NULL))
           != -1) {
        switch (opt) {
            case 'w': w = atoi(optarg); break;
            case 'h': h = atoi(optarg); break;
            case 's': spp = atoi(optarg); break;
            case 't': threads = atoi(optarg); break;
            case 'o': out = optarg; break;
            case 1000: rng_mode = optarg; break;
            case 1001: seed = atoi(optarg); break;
            default: usage(argv[0]); return 2;
        }
    }
    if (w <= 0 || h <= 0 || spp <= 0 || threads <= 0 || !out) {
        usage(argv[0]);
        return 2;
    }
    if (strcmp(rng_mode, "det") == 0) {
        g_rng_global = 0;
    } else if (strcmp(rng_mode, "global") == 0) {
        g_rng_global = 1;
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        uint64_t s = mix64((uint64_t)ts.tv_nsec
                         ^ ((uint64_t)getpid() << 17)
                         ^ ((uint64_t)ts.tv_sec << 33));
        rng_seed(&g_rng, s, 0xabcdefULL);
    } else {
        fprintf(stderr, "unknown --rng mode: %s\n", rng_mode);
        return 2;
    }

    Scene scene;
    memset(&scene, 0, sizeof(scene));
    scene_add(&scene, (V3){ 0, -100.5f, -1 }, 100.0f, (V3){ 0.8f, 0.8f, 0.7f }, 0);
    scene_add(&scene, (V3){ 0, 0, -1 }, 0.5f, (V3){ 0.8f, 0.2f, 0.2f }, 0);
    scene_add(&scene, (V3){ -1.1f, 0, -1.4f }, 0.4f, (V3){ 0.2f, 0.8f, 0.3f }, 0);
    scene_add(&scene, (V3){ 1.1f, -0.1f, -1.6f }, 0.35f, (V3){ 0.3f, 0.4f, 0.9f }, 0);
    scene_add(&scene, (V3){ 0.6f, 1.2f, -2.2f }, 0.25f, (V3){ 1.0f, 0.9f, 0.7f }, 6.0f);
    g_scene = &scene;

    Camera cam = camera_make(w, h);

    float* fb = malloc((size_t)w * h * 3 * sizeof(float));
    if (!fb) { perror("malloc"); return 1; }

    /* Split rows into contiguous bands so every row is owned by exactly
     * one thread for any thread count, including counts that do not divide
     * h (and counts larger than h). */
    int nthreads = threads < h ? threads : h;
    Job* jobs = malloc((size_t)nthreads * sizeof(Job));
    pthread_t* tids = malloc((size_t)nthreads * sizeof(pthread_t));
    if (!jobs || !tids) { perror("malloc"); return 1; }
    int chunk = h / nthreads;
    int rem = h % nthreads;
    int y = 0;
    for (int i = 0; i < nthreads; i++) {
        int count = chunk + (i < rem ? 1 : 0);
        jobs[i] = (Job){ &cam, fb, w, h, spp, seed, y, y + count };
        y += count;
        if (pthread_create(&tids[i], NULL, worker, &jobs[i]) != 0) {
            perror("pthread_create");
            return 1;
        }
    }
    for (int i = 0; i < nthreads; i++) pthread_join(tids[i], NULL);

    int rc = write_pfm(out, w, h, fb);
    free(tids);
    free(jobs);
    free(fb);
    if (rc == 0) {
        fprintf(stderr,
                "rendered %dx%d spp=%d threads=%d rng=%s -> %s\n",
                w, h, spp, threads, rng_mode, out);
    }
    return rc;
}
