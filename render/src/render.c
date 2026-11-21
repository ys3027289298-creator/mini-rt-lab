
/* render: deterministic tiled multithreaded renderer.

   Determinism contract:
     - The image is partitioned into a tiles_x by tiles_y grid of whole
       pixels. A tile owns its pixels exclusively, so writes never race and
       the result cannot depend on thread scheduling.
     - Each pixel accumulates its samples in one fixed order: sample index
       0..spp-1 ascending. The per-(pixel, sample) random stream is seeded
       from (px, py, sample), so sample k has the same value regardless of
       spp, tiling or thread count.
     - Accumulation uses a Neumaier compensated sum; its result is a pure
       function of the fed sequence. With a fixed sequence the pixel value is
       bit-identical for every thread/tile configuration.
     - Samples are normalized to +0.0f, so -0.0f from one lane cannot differ
       from +0.0f in another.
     - The build disables fast-math and FP contraction (see Makefile) so the
       arithmetic is the same for every configuration and every run.
   Tiles only decide which thread evaluates a pixel, never how. */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#define PI_F 3.14159265358979323846f
#define MAX_SPHERES 4
#define MAX_THREADS 64

typedef struct { float x, y, z; } V3;
typedef struct { float r, g, b; } Color;
typedef struct { float ox, oy, oz, dx, dy, dz, tmin, tmax; } Ray;
typedef struct { float x, y, z, r; } Sphere;
typedef struct { double sum, comp; int n; } Acc;

typedef struct {
    int width, height, spp;
    int tiles_x, tiles_y;
    int nspheres;
    Sphere sph[MAX_SPHERES];
    Color* buf;
} Scene;

typedef struct { const Scene* sc; int id, ntiles, stride; } Worker;
typedef struct { int code; int tiles_x, tiles_y, threads; } RenderResult;

static void acc_add(Acc* a, double v) {
    double t = a->sum + v;
    if (fabs(a->sum) >= fabs(v)) a->comp += (a->sum - t) + v;
    else                         a->comp += (v - t) + a->sum;
    a->sum = t;
    a->n++;
}
static double acc_value(const Acc* a) { return a->sum + a->comp; }

static inline float poszero(float x) { return (x == 0.0f) ? 0.0f : x; }
static inline float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    if (l > 0.0f) { a.x /= l; a.y /= l; a.z /= l; }
    return a;
}

static uint32_t pcg_hash(uint32_t v) {
    v ^= v >> 16; v *= 0x7feb352du;
    v ^= v >> 15; v *= 0x846ca68bu;
    v ^= v >> 16;
    return v;
}
/* 24 random bits -> [0,1]; the gap is >= 1/2^24 so it cannot round to 1 */
static inline float u32_to_unit(uint32_t x) {
    return (float)(x >> 8) * (1.0f / 16777216.0f);
}
static uint32_t sample_seed(int px, int py, int s) {
    uint32_t h = (uint32_t)px * 0x9e3779b1u;
    h ^= (uint32_t)py * 0x85ebca77u;
    h ^= (uint32_t)s  * 0xc2b2ae3du;
    return pcg_hash(h);
}

/* cosine-weighted direction in the hemisphere around n */
static V3 cosine_hemisphere(V3 n, float u1, float u2) {
    float phi = 2.0f * PI_F * u2;
    float r = sqrtf(u1);
    V3 local = { r * cosf(phi), r * sinf(phi), sqrtf(1.0f - u1) };
    V3 up = (fabsf(n.z) < 0.999f) ? (V3){ 0.0f, 0.0f, 1.0f } : (V3){ 1.0f, 0.0f, 0.0f };
    V3 b1 = norm3((V3){ n.y * up.z - n.z * up.y,
                        n.z * up.x - n.x * up.z,
                        n.x * up.y - n.y * up.x });
    V3 b2 = (V3){ n.y * b1.z - n.z * b1.y,
                  n.z * b1.x - n.x * b1.z,
                  n.x * b1.y - n.y * b1.x };
    return norm3((V3){ b1.x * local.x + b2.x * local.y + n.x * local.z,
                       b1.y * local.x + b2.y * local.y + n.y * local.z,
                       b1.z * local.x + b2.z * local.y + n.z * local.z });
}

/* nearest hit against the scene spheres; returns sphere index, or -1 */
static int trace(const Scene* sc, const Ray* r, float* t_out) {
    float best = r->tmax;
    int hit = -1;
    for (int i = 0; i < sc->nspheres; i++) {
        V3 d  = { r->dx, r->dy, r->dz };
        V3 oc = { r->ox - sc->sph[i].x, r->oy - sc->sph[i].y, r->oz - sc->sph[i].z };
        float a = dot3(d, d);
        float b = 2.0f * dot3(oc, d);
        float c = dot3(oc, oc) - sc->sph[i].r * sc->sph[i].r;
        float disc = b * b - 4.0f * a * c;
        if (disc < 0.0f) continue;
        float sq = sqrtf(disc);
        float t = (-b - sq) / (2.0f * a);
        if (t < r->tmin) t = (-b + sq) / (2.0f * a);
        if (t >= r->tmin && t < best) { best = t; hit = i; }
    }
    if (hit >= 0) *t_out = best;
    return hit;
}

static void shade_pixel(const Scene* sc, int px, int py, Color* out) {
    Acc ar = { 0, 0, 0 }, ag = { 0, 0, 0 }, ab = { 0, 0, 0 };
    Ray cam = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1e-3f, 1e30f };

    for (int s = 0; s < sc->spp; s++) {
        uint32_t seed = sample_seed(px, py, s);
        float jx = u32_to_unit(pcg_hash(seed ^ 0x68bc21ebu));
        float jy = u32_to_unit(pcg_hash(seed ^ 0x02e5be93u));
        float u1 = u32_to_unit(pcg_hash(seed ^ 0x967a889bu));
        float u2 = u32_to_unit(pcg_hash(seed ^ 0x2b6b3f4du));

        Ray r = cam;
        r.ox = 2.0f * ((float)px + jx) / (float)sc->width - 1.0f;
        r.oy = 2.0f * ((float)py + jy) / (float)sc->height - 1.0f;

        Color c = { 0.0f, 0.0f, 0.0f };
        float t;
        int hit = trace(sc, &r, &t);
        if (hit >= 0) {
            V3 p = { r.ox + t * r.dx, r.oy + t * r.dy, r.oz + t * r.dz };
            V3 n = norm3((V3){ (p.x - sc->sph[hit].x) / sc->sph[hit].r,
                               (p.y - sc->sph[hit].y) / sc->sph[hit].r,
                               (p.z - sc->sph[hit].z) / sc->sph[hit].r });
            V3 bounce = cosine_hemisphere(n, u1, u2);
            float ndl = dot3(n, bounce);
            if (ndl < 0.0f) ndl = 0.0f;
            float tint = (float)(px + py) / (float)(sc->width + sc->height);
            c.r = 0.6f * ndl + 0.4f * (n.x * 0.5f + 0.5f) + 0.05f * tint;
            c.g = 0.6f * ndl + 0.4f * (n.y * 0.5f + 0.5f) + 0.05f * u2;
            c.b = 0.6f * ndl + 0.4f * (n.z * 0.5f + 0.5f) + 0.05f * (1.0f - tint);
        }
        acc_add(&ar, (double)c.r);
        acc_add(&ag, (double)c.g);
        acc_add(&ab, (double)c.b);
    }

    out->r = poszero((float)(acc_value(&ar) / (double)sc->spp));
    out->g = poszero((float)(acc_value(&ag) / (double)sc->spp));
    out->b = poszero((float)(acc_value(&ab) / (double)sc->spp));
}

/* A tile is the rectangle [x0,x1) x [y0,y1) of pixels, computed by
   proportional partition so the tiles cover the image exactly. */
static void render_tile(const Scene* sc, int ti) {
    int tx = ti % sc->tiles_x;
    int ty = ti / sc->tiles_x;
    int x0 = (int)((int64_t)sc->width  * tx       / sc->tiles_x);
    int x1 = (int)((int64_t)sc->width  * (tx + 1) / sc->tiles_x);
    int y0 = (int)((int64_t)sc->height * ty       / sc->tiles_y);
    int y1 = (int)((int64_t)sc->height * (ty + 1) / sc->tiles_y);
    for (int y = y0; y < y1; y++)
        for (int x = x0; x < x1; x++)
            shade_pixel(sc, x, y, &sc->buf[(size_t)y * sc->width + x]);
}

static void* worker_main(void* arg) {
    Worker* w = (Worker*)arg;
    /* static tile distribution: tile ti always goes to the same thread id,
       and since tiles write disjoint pixels, completion order is irrelevant */
    for (int ti = w->id; ti < w->ntiles; ti += w->stride)
        render_tile(w->sc, ti);
    return NULL;
}

/* error codes: 1 bad image/spp, 2 bad tiling (zero/negative tiles),
   3 bad thread count, 4 tile count overflow, 5 thread start failure */
static RenderResult render(const Scene* sc, int threads) {
    RenderResult res = { 0, sc->tiles_x, sc->tiles_y, threads };
    if (sc->width <= 0 || sc->height <= 0 || sc->spp <= 0 || sc->buf == NULL) {
        res.code = 1; return res;
    }
    if (sc->tiles_x <= 0 || sc->tiles_y <= 0) { res.code = 2; return res; }
    if (threads <= 0) { res.code = 3; return res; }
    int64_t ntiles64 = (int64_t)sc->tiles_x * sc->tiles_y;
    if (ntiles64 > (int64_t)INT32_MAX) { res.code = 4; return res; }
    int ntiles = (int)ntiles64;

    int nthreads = threads < ntiles ? threads : ntiles;
    if (nthreads > MAX_THREADS) nthreads = MAX_THREADS;
    res.threads = nthreads;
    pthread_t th[MAX_THREADS];
    Worker wa[MAX_THREADS];
    memset(th, 0, sizeof(th));

    int started = 0;
    for (int k = 1; k < nthreads; k++) {
        wa[k] = (Worker){ sc, k, ntiles, nthreads };
        if (pthread_create(&th[k], NULL, worker_main, &wa[k]) != 0) {
            res.code = 5;
            for (int j = 1; j <= started; j++) pthread_join(th[j], NULL);
            return res;
        }
        started = k;
    }
    wa[0] = (Worker){ sc, 0, ntiles, nthreads };
    worker_main(&wa[0]);
    for (int k = 1; k <= started; k++) pthread_join(th[k], NULL);
    return res;
}

/* ---- self test ---- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static Scene make_scene(Color* buf) {
    Scene sc;
    memset(&sc, 0, sizeof(sc));
    sc.width = 97; sc.height = 61; sc.spp = 13;
    sc.tiles_x = 1; sc.tiles_y = 1;
    sc.buf = buf;
    sc.nspheres = 3;
    sc.sph[0] = (Sphere){ 0.0f, 0.0f, 3.0f, 1.0f };
    sc.sph[1] = (Sphere){ 1.3f, 0.4f, 4.0f, 0.6f };
    sc.sph[2] = (Sphere){ -1.2f, -0.5f, 3.4f, 0.7f };
    return sc;
}

static int same_image(const Color* a, const Color* b, int n) {
    return memcmp(a, b, (size_t)n * sizeof(Color)) == 0;
}
static int all_finite(const Color* img, int n) {
    for (int i = 0; i < n; i++)
        if (!isfinite(img[i].r) || !isfinite(img[i].g) || !isfinite(img[i].b)) return 0;
    return 1;
}
static int image_varies(const Color* img, int n) {
    float lo = img[0].r, hi = img[0].r;
    for (int i = 1; i < n; i++) {
        if (img[i].r < lo) lo = img[i].r;
        if (img[i].r > hi) hi = img[i].r;
    }
    return hi - lo > 1e-3f;
}
static int image_has_negative_zero(const Color* img, int n) {
    for (int i = 0; i < n; i++) {
        if (signbit(img[i].r) || signbit(img[i].g) || signbit(img[i].b)) return 1;
    }
    return 0;
}
/* FNV-1a over the raw pixel bytes: any bit difference changes the hash, so
   the printed value can be diffed across runs, builds and machines. */
static uint64_t image_hash(const Color* img, int n) {
    const unsigned char* p = (const unsigned char*)img;
    uint64_t h = 1469598103934665603ull;
    for (size_t i = 0; i < (size_t)n * sizeof(Color); i++) {
        h ^= p[i];
        h *= 1099511628211ull;
    }
    return h;
}

int main(void) {
    int npix = 97 * 61;
    Color* ref  = malloc(sizeof(Color) * (size_t)npix);
    Color* test = malloc(sizeof(Color) * (size_t)npix);
    Color* work = malloc(sizeof(Color) * (size_t)npix);
    if (!ref || !test || !work) { printf("alloc_failed\n"); return 1; }

    /* reference: single thread, single tile */
    Scene sc = make_scene(ref);
    check("reference_ok", render(&sc, 1).code == 0);
    check("reference_finite", all_finite(ref, npix));
    check("reference_varies", image_varies(ref, npix));
    check("no_negative_zero", !image_has_negative_zero(ref, npix));

    /* zero / negative tiling and thread counts are rejected */
    Scene bad = make_scene(test);
    bad.tiles_x = 0; bad.tiles_y = 4;
    check("tiles_x_zero_rejected", render(&bad, 4).code == 2);
    bad.tiles_x = 4; bad.tiles_y = 0;
    check("tiles_y_zero_rejected", render(&bad, 4).code == 2);
    bad.tiles_x = -1; bad.tiles_y = 2;
    check("tiles_negative_rejected", render(&bad, 4).code == 2);
    bad = make_scene(test);
    check("threads_zero_rejected", render(&bad, 0).code == 3);
    bad.width = 0;
    check("width_zero_rejected", render(&bad, 2).code == 1);
    bad = make_scene(test);
    bad.spp = 0;
    check("spp_zero_rejected", render(&bad, 2).code == 1);

    /* the tiling must cover every pixel exactly (sentinel fill) */
    for (int i = 0; i < npix; i++)
        work[i] = (Color){ -12345.0f, -12345.0f, -12345.0f };
    Scene cover = make_scene(work);
    cover.tiles_x = 7; cover.tiles_y = 5;
    check("cover_render_ok", render(&cover, 3).code == 0);
    int uncovered = 0;
    for (int i = 0; i < npix; i++)
        if (work[i].r == -12345.0f) uncovered++;
    check("covers_all_pixels", uncovered == 0);

    /* tile count and thread count must not change a single bit */
    struct { int tx, ty, th; const char* name; } cases[] = {
        {  1,  1,  1, "tile_1x1_t1" },
        {  1,  1,  8, "tile_1x1_t8" },
        {  3,  2,  1, "tile_3x2_t1" },
        {  3,  2,  4, "tile_3x2_t4" },
        {  7,  5,  3, "tile_7x5_t3" },
        {  7,  5,  8, "tile_7x5_t8" },
        { 16,  9,  2, "tile_16x9_t2" },
        { 64, 64,  8, "tile_64x64_t8" },
        { 97, 61,  6, "tile_perpixel_t6" },
        { 64, 64, 64, "tile_64x64_t64" },
        {1000, 997, 7, "tile_oversub_t7" },
    };
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        Scene c = make_scene(test);
        c.tiles_x = cases[i].tx;
        c.tiles_y = cases[i].ty;
        memset(test, 0, sizeof(Color) * (size_t)npix);
        RenderResult rr = render(&c, cases[i].th);
        check(cases[i].name, rr.code == 0 && same_image(ref, test, npix));
        /* every configuration must also hash to the single-thread image */
        if (!(rr.code == 0 && image_hash(ref, npix) == image_hash(test, npix)))
            failures++;
    }

    /* repeated runs with the same config are identical */
    int repeat_ok = 1;
    for (int k = 0; k < 3; k++) {
        Scene c = make_scene(test);
        c.tiles_x = 16; c.tiles_y = 9;
        memset(test, 0, sizeof(Color) * (size_t)npix);
        render(&c, 5);
        if (!same_image(ref, test, npix)) repeat_ok = 0;
    }
    check("repeat_runs_identical", repeat_ok);

    /* more threads than tiles still matches */
    Scene over = make_scene(test);
    over.tiles_x = 2; over.tiles_y = 2;
    memset(test, 0, sizeof(Color) * (size_t)npix);
    check("oversubscribed_threads", render(&over, 32).code == 0 && same_image(ref, test, npix));

    /* the parallel result must also be finite and cover the image */
    check("parallel_finite", all_finite(test, npix));
    check("parallel_varies", image_varies(test, npix));

    printf("image_hash=%016llx\n", (unsigned long long)image_hash(ref, npix));
    printf("all_ok=%d\n", failures == 0);
    free(ref); free(test); free(work);
    return failures == 0 ? 0 : 1;
}
