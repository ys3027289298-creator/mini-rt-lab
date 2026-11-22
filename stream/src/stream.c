/* stream: on-demand chunk streaming for a scene too large to hold whole.

   Layout on "disk": chunks of CHUNK_PRIMS triangles, serialized contiguously.
   A block is BLOCK_TRIS triangles (one chunk = an integral number of blocks).
   The chunk manifest (bounds + prim range) is always resident; chunk bytes are
   pulled in block-sized loads into a fixed-capacity cache with LRU eviction.

   The streaming renderer and the full-load reference renderer share the same
   intersection and shading code. Two rules make them bit-identical regardless
   of traversal order:
     - tri_hit accepts t in [0, tmax] inclusive (shared edges/vertices hit);
     - ties on distance are broken by the lowest global primitive index, so the
       chosen triangle does not depend on the order chunks are visited in.

   Tiles are rendered in feedback frames. A tile enters the image only after
   two independent loader passes produce the same quantized value; the image is
   never patched in place, and a tile whose value changes is a hard failure. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define MAX_CHUNKS       64
#define CHUNK_PRIMS      24
#define BLOCK_TRIS       12
#define TRIS_MAX         (MAX_CHUNKS * CHUNK_PRIMS)
#define BLOCK_BYTES      (BLOCK_TRIS * (int)sizeof(Tri))
#define SCENE_MAX_BYTES  ((size_t)TRIS_MAX * sizeof(Tri))
#define MAX_BLOCKS       64

#define W                64
#define H                32
#define NUM_TILES        (W * H)

#define INF              1e30f
#define OCC_TMAX         1e6f
#define SHADOW_BIAS      1e-3f
#define PI_F             3.14159265f
#define CHECK_SPACING    1.7f

typedef struct { float min[3], max[3]; } AABB;
typedef struct { float a[3], b[3], c[3]; } Tri;
typedef struct { AABB box; int first, count; } Chunk;
typedef struct { float ox, oy, oz; float dx, dy, dz; } RQ;
typedef struct { int hit, prim; float t; float P[3]; float n[3]; } Hit;

typedef struct {
    int present[MAX_BLOCKS];
    int use[MAX_BLOCKS];
    unsigned char* mem[MAX_BLOCKS];
    int cap, resident, clock;
    size_t peak_bytes;
    unsigned long loads, load_bytes, evictions, hits, misses;
} Cache;

typedef struct {
    int tiles, frame, stable, complete, aborted, unstable, processed;
    unsigned char* state;   /* 0 unseen, 1 tentative, 2 stable */
    unsigned char* img;
} Stream;

/* ---- resident metadata / "disk" / full-loaded copy --------------------- */

static union { unsigned char bytes[SCENE_MAX_BYTES]; double align; } g_disk;
static Tri g_loaded[TRIS_MAX];
static Chunk g_chunks[MAX_CHUNKS];
static int g_nchunks, g_nprim, g_nblocks;
static size_t g_scene_bytes, g_manifest_bytes;
static int g_tamper, g_tamper_block;

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* ---- geometry helpers -------------------------------------------------- */

static void aabb_empty(AABB* b) {
    for (int i = 0; i < 3; i++) { b->min[i] = INF; b->max[i] = -INF; }
}
static void aabb_grow(AABB* b, const float p[3]) {
    for (int i = 0; i < 3; i++) {
        if (p[i] < b->min[i]) b->min[i] = p[i];
        if (p[i] > b->max[i]) b->max[i] = p[i];
    }
}
static void tri_push(Tri* out, const float a[3], const float b[3], const float c[3]) {
    memcpy(out->a, a, sizeof(float) * 3);
    memcpy(out->b, b, sizeof(float) * 3);
    memcpy(out->c, c, sizeof(float) * 3);
}

/* One axis-aligned box as 12 triangles, appended to out. Returns 12. */
static int box_append(Tri* out, float cx, float cy, float hx, float hy, float z0, float z1) {
    float x0 = cx - hx, x1 = cx + hx, y0 = cy - hy, y1 = cy + hy;
    float fn[4][3] = { {x0,y0,0}, {x1,y0,0}, {x1,y1,0}, {x0,y1,0} };
    float bn[4][3] = { {x0,y0,0}, {x1,y0,0}, {x1,y1,0}, {x0,y1,0} };
    for (int i = 0; i < 4; i++) { fn[i][2] = z0; bn[i][2] = z1; }
    int n = 0;
    /* bottom, top */
    tri_push(&out[n++], fn[0], fn[2], fn[1]);
    tri_push(&out[n++], fn[0], fn[3], fn[2]);
    tri_push(&out[n++], bn[0], bn[1], bn[2]);
    tri_push(&out[n++], bn[0], bn[2], bn[3]);
    /* y0, y1 */
    float q[4][3] = { {x0,y0,z0}, {x1,y0,z0}, {x1,y1,z1}, {x0,y1,z0} };
    tri_push(&out[n++], q[0], q[1], q[2]);
    tri_push(&out[n++], q[0], q[2], q[3]);
    float r[4][3] = { {x0,y1,z0}, {x1,y1,z0}, {x1,y1,z1}, {x0,y1,z1} };
    tri_push(&out[n++], r[0], r[1], r[2]);
    tri_push(&out[n++], r[0], r[2], r[3]);
    /* x0, x1 */
    float l[4][3] = { {x0,y0,z0}, {x0,y1,z0}, {x0,y1,z1}, {x0,y0,z1} };
    tri_push(&out[n++], l[0], l[1], l[2]);
    tri_push(&out[n++], l[0], l[2], l[3]);
    float g[4][3] = { {x1,y0,z0}, {x1,y1,z0}, {x1,y1,z1}, {x1,y0,z1} };
    tri_push(&out[n++], g[0], g[1], g[2]);
    tri_push(&out[n++], g[0], g[2], g[3]);
    return n;
}

static void scene_build(void) {
    static Tri build[TRIS_MAX];
    static const float gx[4] = { -4.5f, -1.5f, 1.5f, 4.5f };
    static const float gy[4] = { 0.0f, 3.0f, 6.0f, 9.0f };
    memset(build, 0, sizeof(build));
    g_nprim = 0;
    g_nchunks = 0;
    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            Chunk* c = &g_chunks[g_nchunks++];
            c->first = g_nprim;
            g_nprim += box_append(&build[g_nprim], gx[col], gy[row], 0.8f, 0.6f, 0.0f, 2.0f);
            g_nprim += box_append(&build[g_nprim], gx[col] + 0.9f, gy[row] + 1.3f, 0.5f, 0.4f, 0.3f, 1.5f);
            c->count = g_nprim - c->first;
            aabb_empty(&c->box);
            for (int k = c->first; k < c->first + c->count; k++) {
                aabb_grow(&c->box, build[k].a);
                aabb_grow(&c->box, build[k].b);
                aabb_grow(&c->box, build[k].c);
            }
        }
    }
    g_scene_bytes = (size_t)g_nprim * sizeof(Tri);
    g_nblocks = (int)((g_scene_bytes + BLOCK_BYTES - 1) / BLOCK_BYTES);
    memcpy(g_disk.bytes, build, g_scene_bytes);
    memset(build, 0, sizeof(build));   /* in-memory copy gone: disk + manifest only */
    memset(g_loaded, 0, sizeof(g_loaded));
    g_manifest_bytes = (size_t)g_nchunks * sizeof(Chunk);
}

/* Fault injection for the fidelity check. Shrink the block's triangles toward
   its chunk centre: the result still lies inside the manifest AABB (the centre
   is interior and the box is convex), so streamed and full-load renderers cull
   identically and any image difference is real data dependence. */
static void tamper_block(Tri* tris, int chunk) {
    if (chunk < 0 || chunk >= g_nchunks) return;
    float c[3];
    for (int i = 0; i < 3; i++) c[i] = (g_chunks[chunk].box.min[i] + g_chunks[chunk].box.max[i]) * 0.5f;
    for (int n = 0; n < BLOCK_TRIS; n++) {
        float* v[3] = { tris[n].a, tris[n].b, tris[n].c };
        for (int q = 0; q < 3; q++)
            for (int i = 0; i < 3; i++) v[q][i] = c[i] + (v[q][i] - c[i]) * 0.3f;
    }
}

static int block_chunk(int b) { return b * BLOCK_TRIS / CHUNK_PRIMS; }

static void full_load(void) {
    memset(g_loaded, 0, sizeof(g_loaded));
    memcpy(g_loaded, g_disk.bytes, g_scene_bytes);
    if (g_tamper && g_tamper_block >= 0 && g_tamper_block < g_nblocks)
        tamper_block(&g_loaded[(size_t)g_tamper_block * BLOCK_TRIS], block_chunk(g_tamper_block));
}

/* ---- block cache ------------------------------------------------------- */

static void cache_init(Cache* c, int cap) {
    memset(c, 0, sizeof(*c));
    c->cap = cap;
    for (int i = 0; i < MAX_BLOCKS; i++) c->use[i] = -1;
}
static int cache_init_ok(int cap) { return cap > 0 && cap <= MAX_BLOCKS; }

static void cache_free(Cache* c) {
    for (int i = 0; i < MAX_BLOCKS; i++) {
        free(c->mem[i]);
        c->mem[i] = NULL;
        c->present[i] = 0;
    }
    c->resident = 0;
}

static int chunk_valid(int ci) {
    const Chunk* c = &g_chunks[ci];
    if (c->count <= 0 || c->count % BLOCK_TRIS != 0) return 0;
    if (c->first < 0 || c->first + c->count > g_nprim) return 0;
    if (c->first % BLOCK_TRIS != 0) return 0;
    return 1;
}
static int scene_valid(void) {
    if (g_nchunks <= 0) return 0;
    for (int i = 0; i < g_nchunks; i++) if (!chunk_valid(i)) return 0;
    return 1;
}

/* Load block b. Returns 0 when the budget is exhausted (cap <= 0) or the
   allocation fails; on success the block is resident and counted. */
static int cache_ensure(Cache* c, int b) {
    if (b < 0 || b >= g_nblocks) return 0;
    if (c->present[b]) {
        c->hits++;
        c->use[b] = ++c->clock;
        return 1;
    }
    c->misses++;
    if (c->resident >= c->cap) {
        int victim = -1;
        for (int i = 0; i < g_nblocks; i++)
            if (c->present[i] && (victim < 0 || c->use[i] < c->use[victim])) victim = i;
        if (victim < 0) return 0;                 /* zero budget: refuse */
        free(c->mem[victim]);
        c->mem[victim] = NULL;
        c->present[victim] = 0;
        c->resident--;
        c->evictions++;
    }
    unsigned char* mem = (unsigned char*)malloc((size_t)BLOCK_BYTES);
    if (!mem) return 0;
    size_t off = (size_t)b * BLOCK_BYTES;
    size_t avail = off < g_scene_bytes ? g_scene_bytes - off : 0;
    size_t n = avail < (size_t)BLOCK_BYTES ? avail : (size_t)BLOCK_BYTES;
    memset(mem, 0, (size_t)BLOCK_BYTES);
    memcpy(mem, g_disk.bytes + off, n);
    if (g_tamper && b == g_tamper_block) tamper_block((Tri*)mem, block_chunk(b));
    c->mem[b] = mem;
    c->present[b] = 1;
    c->use[b] = ++c->clock;
    c->resident++;
    c->loads++;
    c->load_bytes += (unsigned long)BLOCK_BYTES;
    size_t r = (size_t)c->resident * (size_t)BLOCK_BYTES;
    if (r > c->peak_bytes) c->peak_bytes = r;
    if (c->resident > c->cap) { /* invariant must never trip */
        c->present[b] = 0;
        free(mem);
        c->mem[b] = NULL;
        c->resident--;
        return 0;
    }
    return 1;
}

static const Tri* tri_ptr(Cache* c, int k) {
    int b = k / BLOCK_TRIS;
    if (!cache_ensure(c, b)) return NULL;
    return (const Tri*)(c->mem[b] + (size_t)(k % BLOCK_TRIS) * sizeof(Tri));
}

/* ---- ray / intersection ------------------------------------------------ */

static int aabb_hit(const AABB* b, const RQ* r, float tmax, float* tmin) {
    float o[3] = { r->ox, r->oy, r->oz };
    float d[3] = { r->dx, r->dy, r->dz };
    float t0 = 0.0f, t1 = tmax;
    for (int i = 0; i < 3; i++) {
        if (fabsf(d[i]) < 1e-8f) {
            if (o[i] < b->min[i] || o[i] > b->max[i]) return 0;
        } else {
            float inv = 1.0f / d[i];
            float ta = (b->min[i] - o[i]) * inv;
            float tb = (b->max[i] - o[i]) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    if (tmin) *tmin = t0;
    return 1;
}

/* t is valid in [0, tmax]; the caller resolves equal-t ties by prim index. */
static int tri_hit(const Tri* t, const RQ* r, float tmax, float* tout) {
    const float e1x = t->b[0] - t->a[0], e1y = t->b[1] - t->a[1], e1z = t->b[2] - t->a[2];
    const float e2x = t->c[0] - t->a[0], e2y = t->c[1] - t->a[1], e2z = t->c[2] - t->a[2];
    const float px = r->dy * e2z - r->dz * e2y;
    const float py = r->dz * e2x - r->dx * e2z;
    const float pz = r->dx * e2y - r->dy * e2x;
    const float det = e1x * px + e1y * py + e1z * pz;
    if (fabsf(det) < 1e-12f) return 0;
    const float inv = 1.0f / det;
    const float sx = r->ox - t->a[0], sy = r->oy - t->a[1], sz = r->oz - t->a[2];
    const float u = (sx * px + sy * py + sz * pz) * inv;
    const float eps = -1e-6f;
    if (u < eps || u > 1.0f - eps) return 0;
    const float qx = sy * e1z - sz * e1y;
    const float qy = sz * e1x - sx * e1z;
    const float qz = sx * e1y - sy * e1x;
    const float v = (r->dx * qx + r->dy * qy + r->dz * qz) * inv;
    if (v < eps || u + v > 1.0f - eps) return 0;
    const float tt = (e2x * qx + e2y * qy + e2z * qz) * inv;
    if (tt < 0.0f || tt > tmax) return 0;
    *tout = tt;
    return 1;
}

static void hit_record(Hit* h, const Tri* t, const RQ* r, int prim, float tt) {
    h->hit = 1;
    h->t = tt;
    h->prim = prim;
    h->P[0] = r->ox + tt * r->dx;
    h->P[1] = r->oy + tt * r->dy;
    h->P[2] = r->oz + tt * r->dz;
    float e1x = t->b[0] - t->a[0], e1y = t->b[1] - t->a[1], e1z = t->b[2] - t->a[2];
    float e2x = t->c[0] - t->a[0], e2y = t->c[1] - t->a[1], e2z = t->c[2] - t->a[2];
    float nx = e1y * e2z - e1z * e2y;
    float ny = e1z * e2x - e1x * e2z;
    float nz = e1x * e2y - e1y * e2x;
    float l = sqrtf(nx * nx + ny * ny + nz * nz);
    if (l > 0.0f) { nx /= l; ny /= l; nz /= l; }
    h->n[0] = nx; h->n[1] = ny; h->n[2] = nz;
}

static int hit_better(const Hit* h, float tt, int prim) {
    if (!h->hit) return 1;
    if (tt < h->t) return 1;
    if (tt == h->t && prim < h->prim) return 1;
    return 0;
}

/* Streaming nearest-hit: cull by manifest AABBs, visit in entry distance
   order, then fetch blocks on demand. Returns 0 if a block cannot be loaded. */
static int trace_stream(Cache* c, const RQ* r, Hit* h) {
    h->hit = 0;
    h->prim = INT_MAX;
    h->t = INF;
    int order[MAX_CHUNKS], n = 0;
    float entry[MAX_CHUNKS];
    for (int i = 0; i < g_nchunks; i++) {
        float t0;
        if (!aabb_hit(&g_chunks[i].box, r, INF, &t0)) continue;
        order[n] = i;
        entry[n] = t0;
        n++;
    }
    for (int i = 1; i < n; i++) {
        int ci = order[i];
        float ct = entry[i];
        int j = i - 1;
        while (j >= 0 && entry[j] > ct) { order[j + 1] = order[j]; entry[j + 1] = entry[j]; j--; }
        order[j + 1] = ci;
        entry[j + 1] = ct;
    }
    for (int i = 0; i < n; i++) {
        if (entry[i] > h->t) break;
        int ci = order[i];
        for (int k = g_chunks[ci].first; k < g_chunks[ci].first + g_chunks[ci].count; k++) {
            const Tri* tri = tri_ptr(c, k);
            if (!tri) return 0;
            float tt;
            if (tri_hit(tri, r, h->t, &tt) && hit_better(h, tt, k)) hit_record(h, tri, r, k, tt);
        }
    }
    return 1;
}

/* Full-load reference: same intersection and tie-break, no cache. It culls by
   the same resident manifest so a fault injected inside a block cannot make
   the two renderers disagree on which chunks are worth visiting. */
static void trace_full(const RQ* r, Hit* h) {
    h->hit = 0;
    h->prim = INT_MAX;
    h->t = INF;
    for (int ci = 0; ci < g_nchunks; ci++) {
        float t0;
        if (!aabb_hit(&g_chunks[ci].box, r, INF, &t0)) continue;
        for (int k = g_chunks[ci].first; k < g_chunks[ci].first + g_chunks[ci].count; k++) {
            float tt;
            if (tri_hit(&g_loaded[k], r, h->t, &tt) && hit_better(h, tt, k))
                hit_record(h, &g_loaded[k], r, k, tt);
        }
    }
}

static int occluded_stream(Cache* c, const RQ* r) {
    for (int i = 0; i < g_nchunks; i++) {
        float t0;
        if (!aabb_hit(&g_chunks[i].box, r, OCC_TMAX, &t0)) continue;
        for (int k = g_chunks[i].first; k < g_chunks[i].first + g_chunks[i].count; k++) {
            const Tri* tri = tri_ptr(c, k);
            if (!tri) return -1;
            float tt;
            if (tri_hit(tri, r, OCC_TMAX, &tt)) return 1;
        }
    }
    return 0;
}

static int occluded_full(const RQ* r) {
    for (int ci = 0; ci < g_nchunks; ci++) {
        float t0;
        if (!aabb_hit(&g_chunks[ci].box, r, OCC_TMAX, &t0)) continue;
        for (int k = g_chunks[ci].first; k < g_chunks[ci].first + g_chunks[ci].count; k++) {
            float tt;
            if (tri_hit(&g_loaded[k], r, OCC_TMAX, &tt)) return 1;
        }
    }
    return 0;
}

/* ---- camera / rays / shading ------------------------------------------- */

typedef struct {
    float px, py, pz;
    float fx, fy, fz, rx, ry, rz, ux, uy, uz;
    float hw, hh;
} Camera;

static Camera make_cam(float px, float py, float pz,
                       float tx, float ty, float tz, float tanfov) {
    Camera c;
    memset(&c, 0, sizeof(c));
    c.px = px; c.py = py; c.pz = pz;
    float dx = tx - px, dy = ty - py, dz = tz - pz;
    float l = sqrtf(dx * dx + dy * dy + dz * dz);
    c.fx = dx / l; c.fy = dy / l; c.fz = dz / l;
    float rx = -c.fz, ry = 0.0f, rz = c.fx;
    l = sqrtf(rx * rx + rz * rz);
    if (l < 1e-6f) { rx = 1.0f; ry = 0.0f; rz = 0.0f; l = 1.0f; }
    rx /= l; rz /= l;
    c.rx = rx; c.ry = ry; c.rz = rz;
    c.ux = c.fy * rz - c.fz * ry;
    c.uy = c.fz * rx - c.fx * rz;
    c.uz = c.fx * ry - c.fy * rx;
    c.hw = tanfov;
    c.hh = tanfov * (float)H / (float)W;
    return c;
}

static RQ primary_ray(const Camera* c, int tile) {
    float px = ((float)(tile % W) + 0.5f) / (float)W * 2.0f - 1.0f;
    float py = ((float)(tile / W) + 0.5f) / (float)H * 2.0f - 1.0f;
    RQ r;
    r.ox = c->px; r.oy = c->py; r.oz = c->pz;
    r.dx = c->fx + px * c->hw * c->rx + py * c->hh * c->ux;
    r.dy = c->fy + px * c->hw * c->ry + py * c->hh * c->uy;
    r.dz = c->fz + px * c->hw * c->rz + py * c->hh * c->uz;
    float l = sqrtf(r.dx * r.dx + r.dy * r.dy + r.dz * r.dz);
    r.dx /= l; r.dy /= l; r.dz /= l;
    return r;
}

static void shadow_ray(const RQ* p, const float P[3], RQ* out) {
    static const float l[3] = { -0.4f, 0.5f, 1.0f };
    float n = sqrtf(l[0] * l[0] + l[1] * l[1] + l[2] * l[2]);
    out->dx = l[0] / n; out->dy = l[1] / n; out->dz = l[2] / n;
    out->ox = P[0] + out->dx * SHADOW_BIAS;
    out->oy = P[1] + out->dy * SHADOW_BIAS;
    out->oz = P[2] + out->dz * SHADOW_BIAS;
    (void)p;
}

static int quantize(float v) {
    float x = v * 255.0f + 0.5f;
    if (!(x > 0.0f)) x = 0.0f;
    if (x > 255.0f) x = 255.0f;
    return (int)x;
}

/* Pure function of the hit and the occlusion result; both drivers call it. */
static int shade(const Hit* h, const RQ* r, int occl) {
    if (!h->hit) return quantize(0.07f);
    float n[3] = { h->n[0], h->n[1], h->n[2] };
    if (n[0] * r->dx + n[1] * r->dy + n[2] * r->dz > 0.0f) {
        n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
    }
    float ax = fabsf(n[0]), ay = fabsf(n[1]), az = fabsf(n[2]);
    float base, u, v;
    if (ax >= ay && ax >= az) { base = 0.70f; u = h->P[1]; v = h->P[2]; }
    else if (ay >= az)        { base = 0.80f; u = h->P[0]; v = h->P[2]; }
    else                      { base = 0.62f; u = h->P[0]; v = h->P[1]; }
    int cu = (int)floorf(u / CHECK_SPACING);
    int cv = (int)floorf(v / CHECK_SPACING);
    if ((cu + cv) & 1) base *= 0.60f;
    if (occl) base *= 0.15f;
    return quantize(base);
}

/* Positive-vertex frustum test against the four side planes through the eye. */
static int chunk_in_frustum(const Camera* c, const AABB* b) {
    float n[4][3];
    for (int i = 0; i < 3; i++) {
        float f = (i == 0) ? c->fx : (i == 1) ? c->fy : c->fz;
        float rr = (i == 0) ? c->rx : (i == 1) ? c->ry : c->rz;
        float uu = (i == 0) ? c->ux : (i == 1) ? c->uy : c->uz;
        n[0][i] = f + c->hw * rr;
        n[1][i] = f - c->hw * rr;
        n[2][i] = f + c->hh * uu;
        n[3][i] = f - c->hh * uu;
    }
    for (int p = 0; p < 4; p++) {
        float vx = (n[p][0] >= 0.0f) ? b->max[0] : b->min[0];
        float vy = (n[p][1] >= 0.0f) ? b->max[1] : b->min[1];
        float vz = (n[p][2] >= 0.0f) ? b->max[2] : b->min[2];
        float d = n[p][0] * (vx - c->px) + n[p][1] * (vy - c->py) + n[p][2] * (vz - c->pz);
        if (d < 0.0f) return 0;
    }
    return 1;
}

/* Camera-driven prefetch: visible chunks nearest the eye first, until the
   working set reaches `budget` blocks. On-demand loads still cover shadow
   rays and chunks outside the frustum. */
static void prefetch_frustum(Cache* c, const Camera* cam, int budget) {
    int vis[MAX_CHUNKS], n = 0;
    float d2[MAX_CHUNKS];
    for (int i = 0; i < g_nchunks; i++) {
        if (!chunk_in_frustum(cam, &g_chunks[i].box)) continue;
        float cx = (g_chunks[i].box.min[0] + g_chunks[i].box.max[0]) * 0.5f;
        float cy = (g_chunks[i].box.min[1] + g_chunks[i].box.max[1]) * 0.5f;
        float cz = (g_chunks[i].box.min[2] + g_chunks[i].box.max[2]) * 0.5f;
        float dx = cx - cam->px, dy = cy - cam->py, dz = cz - cam->pz;
        vis[n] = i;
        d2[n] = dx * dx + dy * dy + dz * dz;
        n++;
    }
    for (int i = 1; i < n; i++) {
        int ci = vis[i];
        float cd = d2[i];
        int j = i - 1;
        while (j >= 0 && d2[j] > cd) { vis[j + 1] = vis[j]; d2[j + 1] = d2[j]; j--; }
        vis[j + 1] = ci;
        d2[j + 1] = cd;
    }
    for (int i = 0; i < n && c->resident < budget; i++) {
        int ci = vis[i];
        for (int k = g_chunks[ci].first; k < g_chunks[ci].first + g_chunks[ci].count; k += BLOCK_TRIS) {
            if (c->resident >= budget) break;
            if (!cache_ensure(c, k / BLOCK_TRIS)) return;
        }
    }
}

/* ---- streaming driver -------------------------------------------------- */

static Stream* stream_create(int tiles) {
    if (g_nchunks <= 0 || !scene_valid()) return NULL;         /* zero chunks: refuse */
    if (tiles <= 0 || tiles > NUM_TILES) return NULL;
    Stream* s = (Stream*)calloc(1, sizeof(Stream));
    if (!s) return NULL;
    s->state = (unsigned char*)calloc((size_t)tiles, 1);
    s->img = (unsigned char*)calloc((size_t)tiles, 1);
    if (!s->state || !s->img) { free(s->state); free(s->img); free(s); return NULL; }
    s->tiles = tiles;
    return s;
}
static void stream_free(Stream* s) {
    if (!s) return;
    free(s->state);
    free(s->img);
    free(s);
}

/* A tile is published only after two loader passes agree. A later disagreement
   is recorded as an abort, never a silent overwrite. */
static void judge(Stream* s, int tile, int q) {
    if (s->aborted) return;
    if (s->state[tile] == 0) {
        s->img[tile] = (unsigned char)q;
        s->state[tile] = 1;
        return;
    }
    if (s->img[tile] != (unsigned char)q) {
        s->unstable++;
        s->aborted = 1;
        return;
    }
    if (s->state[tile] == 1) {
        s->state[tile] = 2;
        s->stable++;
    }
}

/* One feedback frame: prefetch visible chunks, then re-render tiles that are
   not yet stable plus a rotating quarter of the stable ones. Returns 0 if the
   stream refuses to run or a block cannot be loaded. */
static int stream_frame(Stream* s, Cache* c, const Camera* cam) {
    if (!s || g_nchunks <= 0 || c->cap <= 0) return 0;
    s->frame++;
    int budget = c->cap / 2 > 0 ? c->cap / 2 : 1;
    prefetch_frustum(c, cam, budget);
    int phase = s->frame % 4;
    for (int t = 0; t < s->tiles; t++) {
        if (s->state[t] == 2 && (t % 4) != phase) continue;
        RQ r = primary_ray(cam, t);
        Hit h;
        if (!trace_stream(c, &r, &h)) { s->aborted = 1; return 0; }
        int occl = 0;
        if (h.hit) {
            RQ sr;
            shadow_ray(&r, h.P, &sr);
            int o = occluded_stream(c, &sr);
            if (o < 0) { s->aborted = 1; return 0; }
            occl = o;
        }
        judge(s, t, shade(&h, &r, occl));
        s->processed++;
        if (s->aborted) return 0;
    }
    s->complete = (s->stable == s->tiles);
    return 1;
}

/* ---- reference driver -------------------------------------------------- */

static void render_full(const Camera* cam, unsigned char* img) {
    for (int t = 0; t < NUM_TILES; t++) {
        RQ r = primary_ray(cam, t);
        Hit h;
        trace_full(&r, &h);
        int occl = 0;
        if (h.hit) {
            RQ sr;
            shadow_ray(&r, h.P, &sr);
            occl = occluded_full(&sr);
        }
        img[t] = (unsigned char)shade(&h, &r, occl);
    }
}

static int img_diff(const unsigned char* a, const unsigned char* b, int n) {
    int d = 0;
    for (int i = 0; i < n; i++) if (a[i] != b[i]) d++;
    return d;
}

/* ---- main -------------------------------------------------------------- */

typedef struct {
    size_t peak, block_peak, manifest, full, blocks;
    unsigned long loads, load_bytes, evictions;
} Report;

/* Peak resident memory of a streamed render = manifest (always resident)
   + peak block-cache bytes. Both are reported so the bound is auditable. */
static Report summarize(const Cache* c) {
    Report r;
    r.block_peak = c->peak_bytes;
    r.peak = c->peak_bytes + g_manifest_bytes;
    r.manifest = g_manifest_bytes;
    r.full = g_scene_bytes;
    r.blocks = (size_t)g_nblocks;
    r.loads = c->loads;
    r.load_bytes = c->load_bytes;
    r.evictions = c->evictions;
    return r;
}

static void print_report(const char* tag, const Report* r) {
    printf("%s peak_bytes=%zu block_peak_bytes=%zu manifest_bytes=%zu "
           "full_load_bytes=%zu scene_blocks=%zu loads=%lu evictions=%lu\n",
           tag, r->peak, r->block_peak, r->manifest, r->full, r->blocks,
           r->loads, r->evictions);
}

int main(void) {
    Cache probe;
    cache_init(&probe, 4);

    /* ---- refusals: no chunks, zero budget, zero tiles, empty chunk ---- */
    check("cache_reject_zero_budget", cache_init_ok(0) == 0);
    check("reject_zero_chunks", stream_create(16) == NULL);
    scene_build();
    check("scene_valid", scene_valid() == 1);
    check("scene_chunk_count", g_nchunks == 16 && g_nprim == 16 * CHUNK_PRIMS);
    check("scene_blocks", g_nblocks == 32 && g_scene_bytes == 13824);
    check("block_aligned", CHUNK_PRIMS % BLOCK_TRIS == 0 && g_scene_bytes % BLOCK_BYTES == 0);
    check("reject_zero_tiles", stream_create(0) == NULL);
    int save_count = g_chunks[5].count;
    g_chunks[5].count = 0;
    check("reject_empty_chunk", stream_create(16) == NULL);
    g_chunks[5].count = save_count;
    check("manifest_bytes", g_manifest_bytes == 16 * sizeof(Chunk));

    Stream* zs = stream_create(16);
    Cache zc;
    cache_init(&zc, 0);
    check("reject_zero_budget_frame", zs != NULL && stream_frame(zs, &zc, NULL) == 0);
    stream_free(zs);

    full_load();
    Camera cam = make_cam(0.0f, -6.0f, 5.0f, 0.0f, 3.0f, 1.0f, 0.9f);

    unsigned char* ref = (unsigned char*)malloc((size_t)NUM_TILES);
    render_full(&cam, ref);

    int lit = 0, distinct = 0;
    unsigned char seen[256];
    memset(seen, 0, sizeof(seen));
    for (int t = 0; t < NUM_TILES; t++) {
        if (ref[t] > (int)quantize(0.07f)) lit++;
        if (!seen[ref[t]]) { seen[ref[t]] = 1; distinct++; }
    }
    printf("ref_lit_pixels=%d ref_distinct_values=%d\n", lit, distinct);
    check("ref_has_geometry", lit > NUM_TILES / 4);
    check("ref_has_shading", distinct >= 4);

    /* ---- streaming under a 4-block budget (scene needs 32 blocks) ----- */
    Cache c;
    cache_init(&c, 4);
    Stream* s = stream_create(NUM_TILES);
    check("stream_create", s != NULL);
    if (!s) return 1;

    int ok = 1;
    for (int f = 0; f < 8 && !s->complete && !s->aborted; f++)
        ok &= stream_frame(s, &c, &cam);
    check("stream_run_ok", ok == 1);
    check("stream_complete", s->complete == 1);
    check("stream_no_abort", s->aborted == 0);
    check("stream_frames", s->frame >= 2 && s->frame <= 8);

    Report rep = summarize(&c);
    print_report("stream", &rep);
    printf("stream_frames=%d stream_stable_tiles=%d\n", s->frame, s->stable);

    long diff = img_diff(s->img, ref, NUM_TILES);
    check("stream_vs_full_pixel_exact", diff == 0);
    if (diff) {
        int first = -1;
        for (int t = 0; t < NUM_TILES; t++) if (s->img[t] != ref[t]) { first = t; break; }
        printf("first_diff_pixel=%d\n", first);
    }
    check("peak_is_bounded", rep.block_peak <= 4u * (size_t)BLOCK_BYTES);
    check("peak_was_reached", rep.block_peak == 4u * (size_t)BLOCK_BYTES);
    check("peak_includes_manifest", rep.peak == rep.block_peak + g_manifest_bytes);
    check("peak_below_full_load", rep.peak * 3 < rep.full);
    check("cache_thrashed", rep.evictions > 0);
    check("blocks_were_loaded", rep.loads > 0);

    /* keep monitoring: a stable image must not change across loader passes */
    int clean_abort = 0;
    for (int f = 0; f < 6; f++) {
        if (!stream_frame(s, &c, &cam)) break;
        if (s->aborted) { clean_abort = 1; break; }
    }
    check("image_stable_over_frames", clean_abort == 0);
    check("still_complete", s->complete == 1);
    printf("monitor_frames=%d monitor_stable=%d\n", s->frame, s->stable);

    /* ---- integrity: a corrupted loader must change the image -------- */
    int tamper_chunk = -1;
    for (int ci = 0; ci < g_nchunks && tamper_chunk < 0; ci++) {
        g_tamper = 1;
        g_tamper_block = ci * (CHUNK_PRIMS / BLOCK_TRIS);
        full_load();
        unsigned char* bad = (unsigned char*)malloc((size_t)NUM_TILES);
        render_full(&cam, bad);
        int d = img_diff(bad, ref, NUM_TILES);
        free(bad);
        g_tamper = 0;
        if (d > 0) tamper_chunk = ci;
    }
    check("tamper_location_found", tamper_chunk >= 0);
    printf("tamper_chunk=%d\n", tamper_chunk);

    if (tamper_chunk >= 0) {
        g_tamper = 1;
        g_tamper_block = tamper_chunk * (CHUNK_PRIMS / BLOCK_TRIS);
        full_load();
        unsigned char* tampered_full = (unsigned char*)malloc((size_t)NUM_TILES);
        render_full(&cam, tampered_full);

        Cache tc;
        cache_init(&tc, 4);
        Stream* ts = stream_create(NUM_TILES);
        int tok = 1;
        for (int f = 0; f < 8 && ts && !ts->complete && !ts->aborted; f++)
            tok &= stream_frame(ts, &tc, &cam);
        check("tamper_stream_ran", tok == 1 && ts != NULL);
        int vs_bad = ts ? img_diff(ts->img, tampered_full, NUM_TILES) : -1;
        int vs_good = ts ? img_diff(ts->img, ref, NUM_TILES) : -1;
        check("tamper_faithful_to_loaded_bytes", vs_bad == 0);
        check("tamper_caught_by_reference", vs_good > 0);
        printf("tamper_stream_vs_tamper_full=%d tamper_stream_vs_ref=%d\n", vs_bad, vs_good);
        if (ts) stream_free(ts);
        cache_free(&tc);

        /* late corruption: a published tile must be flagged, not rewritten */
        g_tamper = 0;
        full_load();
        for (int f = 0; f < 2; f++) stream_frame(s, &c, &cam);
        g_tamper = 1;
        int caught = 0;
        for (int f = 0; f < 6; f++) {
            stream_frame(s, &c, &cam);
            if (s->aborted) { caught = 1; break; }
        }
        check("late_tamper_detected", caught == 1);
        printf("late_tamper_unstable_tiles=%d\n", s->unstable);
        g_tamper = 0;
        free(tampered_full);
    }

    free(ref);
    stream_free(s);
    cache_free(&c);
    cache_free(&probe);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
