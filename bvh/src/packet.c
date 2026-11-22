
/* packet: direction-grouped ray packet traversal vs. scalar reference.
 *
 * Coherent rays are bucketed by direction into packets of PACKET_K rays.
 * A packet walks the BVH once with an active-ray mask:
 *   - the node's bounds are loaded a single time and shared by every active
 *     ray (one node fetch per packet instead of one per ray);
 *   - a slab-interval envelope over the packet gives an exact all-miss
 *     proof, so a subtree can be pruned for the whole packet without any
 *     per-ray predicate;
 *   - otherwise each active ray runs the identical slab arithmetic used by
 *     the scalar path; rays that miss are masked out of the subtree.
 *
 * Exactness: pruning only ever removes tests whose own exact predicate says
 * miss, and every ray keeps private tmax / hit state.  Nearest-hit id and
 * hit distance are therefore bit-identical to scalar traversal and do not
 * depend on packet composition or input ray order.
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MAX_PRIM 256
#define MAX_NODE 600
#define MAX_RAYS 1024
#define PACKET_K 8
#define NBUCKETS 4                     /* transverse slope buckets per axis */
#define NBINS (6 * NBUCKETS * NBUCKETS)
#define STACK_MAX 64

typedef struct { float min[3], max[3]; } AABB;
typedef struct { AABB box; float cx, cy, cz; int id; } Prim;
typedef struct { AABB box; int left, right; int first, count; } Node;

typedef struct {
    Prim p[MAX_PRIM];
    Node n[MAX_NODE];
    int nprim, nnode;
} Tree;

typedef struct {
    float o[3], d[3];
    int pixel;                     /* image-space id, stable across reorder */
    int hit;                       /* nearest prim id, -1 on miss */
    float t;                       /* nearest hit distance */
} Ray;

typedef struct {
    long node_fetches;             /* node boxes loaded/tested */
    long ray_box_tests;            /* exact per-ray box predicates evaluated */
    long prim_tests;               /* per-ray primitive tests */
} Stats;

/* ---------- AABB math, shared verbatim by both traversal paths ---------- */

static void aabb_empty(AABB* b) {
    for (int i = 0; i < 3; i++) { b->min[i] = 1e30f; b->max[i] = -1e30f; }
}
static void aabb_grow(AABB* b, const AABB* o) {
    for (int i = 0; i < 3; i++) {
        if (o->min[i] < b->min[i]) b->min[i] = o->min[i];
        if (o->max[i] > b->max[i]) b->max[i] = o->max[i];
    }
}

/* Raw per-slab ray/box intervals, unclipped by [0, tmax].
 * Returns 1 when a parallel slab proves the ray misses the box. */
static int slab_intervals(const AABB* b, const float o[3], const float d[3],
                          float ta[3], float tb[3], int* parallel) {
    *parallel = 0;
    for (int i = 0; i < 3; i++) {
        if (fabsf(d[i]) < 1e-8f) {
            if (o[i] < b->min[i] || o[i] > b->max[i]) return 1;
            ta[i] = -INFINITY; tb[i] = INFINITY; *parallel = 1;
        } else {
            float inv = 1.0f / d[i];
            float u = (b->min[i] - o[i]) * inv;
            float v = (b->max[i] - o[i]) * inv;
            ta[i] = u < v ? u : v;
            tb[i] = u < v ? v : u;
        }
    }
    return 0;
}

/* Combine slab intervals with the ray's [0, tmax] window. */
static int slab_combine(const float ta[3], const float tb[3],
                        float tmax, float* tmin) {
    float t0 = 0.0f, t1 = tmax;
    for (int i = 0; i < 3; i++) {
        if (ta[i] > t0) t0 = ta[i];
        if (tb[i] < t1) t1 = tb[i];
        if (t0 > t1) return 0;
    }
    if (tmin) *tmin = t0;
    return 1;
}

/* Scalar entry point: identical predicate semantics to bvh.c's aabb_hit. */
static int aabb_hit(const AABB* b, const float o[3], const float d[3],
                    float tmax, float* tmin) {
    float ta[3], tb[3]; int parallel;
    if (slab_intervals(b, o, d, ta, tb, &parallel)) return 0;
    return slab_combine(ta, tb, tmax, tmin);
}

/* ---------- median-split build (same algorithm as bvh.c) ---------- */

static int tree_add(Tree* t, float minx, float miny, float minz,
                    float maxx, float maxy, float maxz) {
    if (t->nprim >= MAX_PRIM) return -1;
    Prim* p = &t->p[t->nprim];
    p->box.min[0] = minx; p->box.min[1] = miny; p->box.min[2] = minz;
    p->box.max[0] = maxx; p->box.max[1] = maxy; p->box.max[2] = maxz;
    p->cx = (minx + maxx) * 0.5f; p->cy = (miny + maxy) * 0.5f;
    p->cz = (minz + maxz) * 0.5f;
    p->id = t->nprim;
    return t->nprim++;
}

static int build(Tree* t, int first, int count) {
    if (t->nnode >= MAX_NODE) return -1;
    int ni = t->nnode++;
    Node* n = &t->n[ni];
    aabb_empty(&n->box);
    for (int i = first; i < first + count; i++) aabb_grow(&n->box, &t->p[i].box);

    if (count <= 2) {
        n->left = n->right = -1;
        n->first = first; n->count = count;
        return ni;
    }
    float w[3] = { n->box.max[0] - n->box.min[0],
                   n->box.max[1] - n->box.min[1],
                   n->box.max[2] - n->box.min[2] };
    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2)
                             : (w[1] > w[2] ? 1 : 2);
    float c = n->box.min[axis] + w[axis] * 0.5f;
    int mid = first;
    for (int i = first; i < first + count; i++) {
        float v = (axis == 0) ? t->p[i].cx : (axis == 1 ? t->p[i].cy : t->p[i].cz);
        if (v < c) {
            Prim tmp = t->p[i]; t->p[i] = t->p[mid]; t->p[mid] = tmp;
            mid++;
        }
    }
    if (mid == first || mid == first + count) mid = first + count / 2;
    n->first = -1; n->count = 0;
    n->left = build(t, first, mid - first);
    n->right = build(t, mid, first + count - mid);
    return ni;
}

/* ---------- scalar reference: one ray at a time, nearest hit ---------- */

static void trace_scalar(const Tree* t, Ray* rays, int nrays, Stats* s) {
    for (int r = 0; r < nrays; r++) {
        Ray* ray = &rays[r];
        ray->hit = -1; ray->t = INFINITY;
        int stack[STACK_MAX];
        int sp = 0;
        stack[sp++] = 0;
        while (sp > 0) {
            const Node* nd = &t->n[stack[--sp]];
            s->node_fetches++;
            s->ray_box_tests++;
            float tmin;
            if (!aabb_hit(&nd->box, ray->o, ray->d, ray->t, &tmin)) continue;
            if (nd->count > 0) {
                for (int i = 0; i < nd->count; i++) {
                    const Prim* pr = &t->p[nd->first + i];
                    s->prim_tests++;
                    float tt;
                    if (aabb_hit(&pr->box, ray->o, ray->d, ray->t, &tt) &&
                        tt < ray->t) {
                        ray->t = tt;
                        ray->hit = pr->id;
                    }
                }
            } else {
                stack[sp++] = nd->right;
                stack[sp++] = nd->left;   /* left popped first, fixed order */
            }
        }
    }
}

/* ---------- direction binning: 6 faces x NBUCKETS^2 slope grid ---------- */

static int dir_bin(const float d[3]) {
    float ax = fabsf(d[0]), ay = fabsf(d[1]), az = fabsf(d[2]);
    int face;
    if (ax >= ay && ax >= az) face = d[0] >= 0.0f ? 0 : 1;
    else if (ay >= az)        face = d[1] >= 0.0f ? 2 : 3;
    else                      face = d[2] >= 0.0f ? 4 : 5;
    int face_axis = face >> 1;
    int axes[2]; int k = 0;
    for (int i = 0; i < 3; i++) if (i != face_axis) axes[k++] = i;
    float denom = fabsf(d[face_axis]);
    int bucket[2];
    for (int i = 0; i < 2; i++) {
        float q = d[axes[i]] / denom;              /* in [-1, 1] */
        int b = (int)((q + 1.0f) * 0.5f * (float)NBUCKETS);
        if (b < 0) b = 0;
        if (b >= NBUCKETS) b = NBUCKETS - 1;
        bucket[i] = b;
    }
    return face * (NBUCKETS * NBUCKETS) + bucket[0] * NBUCKETS + bucket[1];
}

/* Counting-sort ray indices into direction bins; output in order[].
 * Rays inside each bin are additionally sorted by pixel id, so packet
 * composition (and every traversal statistic) is canonical and independent
 * of the input ray order. */
static void group_by_direction(const Ray* rays, int nrays, int* order) {
    int count[NBINS];
    int offset[NBINS];
    int start[NBINS];
    memset(count, 0, sizeof(count));
    for (int r = 0; r < nrays; r++) count[dir_bin(rays[r].d)]++;
    int acc = 0;
    for (int b = 0; b < NBINS; b++) {
        start[b] = acc; offset[b] = acc; acc += count[b];
    }
    for (int r = 0; r < nrays; r++) {
        int b = dir_bin(rays[r].d);
        order[offset[b]++] = r;
    }
    /* insertion-sort each bin block by pixel id (blocks are small) */
    for (int b = 0; b < NBINS; b++) {
        for (int i = start[b] + 1; i < offset[b]; i++) {
            int key = order[i];
            int j = i - 1;
            while (j >= start[b] && rays[order[j]].pixel > rays[key].pixel) {
                order[j + 1] = order[j];
                j--;
            }
            order[j + 1] = key;
        }
    }
}

/* ---------- packet traversal: one shared walk per coherent batch ---------- */

typedef struct { int node; uint32_t mask; } Frame;

static void trace_packets(const Tree* t, Ray* rays, int nrays,
                          const int* order, Stats* s) {
    for (int r = 0; r < nrays; r++) { rays[r].hit = -1; rays[r].t = INFINITY; }

    int pos = 0;
    while (pos < nrays) {
        int bin = dir_bin(rays[order[pos]].d);
        int k = 0;
        int idx[PACKET_K];
        while (pos < nrays && k < PACKET_K &&
               dir_bin(rays[order[pos]].d) == bin) {
            idx[k++] = order[pos++];
        }

        Frame stack[STACK_MAX];
        int sp = 0;
        stack[sp++] = (Frame){ 0, (1u << k) - 1u };
        while (sp > 0) {
            Frame f = stack[--sp];
            const Node* nd = &t->n[f.node];
            s->node_fetches++;   /* box loaded once, shared by the packet */

            float ta[PACKET_K][3], tb[PACKET_K][3];
            int bad[PACKET_K] = { 0 }, parallel = 0;
            float lo[3] = { INFINITY, INFINITY, INFINITY };
            float hi[3] = { -INFINITY, -INFINITY, -INFINITY };
            for (int j = 0; j < k; j++) {
                if (!(f.mask & (1u << j))) continue;
                const Ray* ray = &rays[idx[j]];
                if (slab_intervals(&nd->box, ray->o, ray->d,
                                   ta[j], tb[j], &parallel)) {
                    bad[j] = 1;
                } else if (!parallel) {
                    for (int i = 0; i < 3; i++) {
                        if (ta[j][i] < lo[i]) lo[i] = ta[j][i];
                        if (tb[j][i] > hi[i]) hi[i] = tb[j][i];
                    }
                }
            }

            uint32_t hitmask = 0;
            if (!parallel) {
                /* Exact all-miss proof: envelope empty on some slab pair
                 * implies every per-ray interval is disjoint (see README). */
                float L = fmaxf(fmaxf(lo[0], lo[1]), lo[2]);
                float H = fminf(fminf(hi[0], hi[1]), hi[2]);
                if (L > H) continue;
            }
            for (int j = 0; j < k; j++) {
                if (!(f.mask & (1u << j)) || bad[j]) continue;
                s->ray_box_tests++;
                Ray* ray = &rays[idx[j]];
                float tmin;
                if (slab_combine(ta[j], tb[j], ray->t, &tmin))
                    hitmask |= (1u << j);
            }
            if (!hitmask) continue;

            if (nd->count > 0) {
                for (int j = 0; j < k; j++) {
                    if (!(hitmask & (1u << j))) continue;
                    Ray* ray = &rays[idx[j]];
                    for (int i = 0; i < nd->count; i++) {
                        const Prim* pr = &t->p[nd->first + i];
                        s->prim_tests++;
                        float tt;
                        if (aabb_hit(&pr->box, ray->o, ray->d, ray->t, &tt) &&
                            tt < ray->t) {
                            ray->t = tt;
                            ray->hit = pr->id;
                        }
                    }
                }
            } else {
                stack[sp++] = (Frame){ nd->right, hitmask };
                stack[sp++] = (Frame){ nd->left,  hitmask };
            }
        }
    }
}

/* ---------- test scene and ray sets ---------- */

static void build_grid_scene(Tree* t) {
    memset(t, 0, sizeof(*t));
    for (int x = 0; x < 8; x++)
        for (int y = 0; y < 8; y++)
            for (int z = 0; z < 2; z++) {
                float cx = x * 1.1f, cy = y * 1.1f, cz = z * 1.1f;
                tree_add(t, cx - 0.45f, cy - 0.45f, cz - 0.45f,
                            cx + 0.45f, cy + 0.45f, cz + 0.45f);
            }
}

/* Pinhole camera: neighboring pixels have nearly equal directions. */
static int make_coherent_rays(Ray* rays, int side) {
    int n = 0;
    float o[3] = { 3.85f, 3.85f, -9.0f };
    for (int y = 0; y < side; y++)
        for (int x = 0; x < side; x++) {
            Ray* r = &rays[n];
            float tx = ((float)x + 0.5f) / (float)side * 8.3f - 0.3f;
            float ty = ((float)y + 0.5f) / (float)side * 8.3f - 0.3f;
            float tz = 0.55f;
            r->o[0] = o[0]; r->o[1] = o[1]; r->o[2] = o[2];
            float v[3] = { tx - o[0], ty - o[1], tz - o[2] };
            float inv = 1.0f / sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            r->d[0] = v[0] * inv; r->d[1] = v[1] * inv; r->d[2] = v[2] * inv;
            r->pixel = n; r->hit = -1; r->t = INFINITY;
            n++;
        }
    return n;
}

static unsigned lcg(unsigned* state) {
    *state = *state * 1103515245u + 12345u;
    return *state;
}

/* Incoherent control set: uniform-ish random directions from one origin. */
static int make_incoherent_rays(Ray* rays, int n, unsigned seed) {
    unsigned state = seed;
    float o[3] = { 3.85f, 3.85f, -4.0f };
    for (int r = 0; r < n; r++) {
        float v[3];
        do {
            for (int i = 0; i < 3; i++)
                v[i] = (float)(lcg(&state) & 0xffff) / 32767.5f - 1.0f;
        } while (v[0] * v[0] + v[1] * v[1] + v[2] * v[2] > 1.0f ||
                 v[0] * v[0] + v[1] * v[1] + v[2] * v[2] < 1e-6f);
        float inv = 1.0f / sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        rays[r].o[0] = o[0]; rays[r].o[1] = o[1]; rays[r].o[2] = o[2];
        rays[r].d[0] = v[0] * inv; rays[r].d[1] = v[1] * inv;
        rays[r].d[2] = v[2] * inv;
        rays[r].pixel = r; rays[r].hit = -1; rays[r].t = INFINITY;
    }
    return n;
}

static int results_equal(const Ray* a, const Ray* b, int n);
static void shuffle_rays(Ray* rays, int n, unsigned seed);

static int stress_scene(unsigned seed) {
    Tree t;
    unsigned st = seed;
    memset(&t, 0, sizeof(t));
    int nprims = 1 + (int)(lcg(&st) % 60u);
    for (int i = 0; i < nprims; i++) {
        float cx = (float)(lcg(&st) % 1600u) / 100.0f - 8.0f;
        float cy = (float)(lcg(&st) % 1600u) / 100.0f - 8.0f;
        float cz = (float)(lcg(&st) % 1600u) / 100.0f - 8.0f;
        float hx = 0.05f + (float)(lcg(&st) % 90u) / 100.0f;
        float hy = 0.05f + (float)(lcg(&st) % 90u) / 100.0f;
        float hz = 0.05f + (float)(lcg(&st) % 90u) / 100.0f;
        if (tree_add(&t, cx - hx, cy - hy, cz - hz,
                         cx + hx, cy + hy, cz + hz) < 0) return 0;
    }
    if (build(&t, 0, t.nprim) != 0) return 0;

    Ray a[MAX_RAYS], b[MAX_RAYS];
    int order[MAX_RAYS];
    int nrays = 96;
    for (int r = 0; r < nrays; r++) {
        Ray* ray = &a[r];
        for (int i = 0; i < 3; i++)
            ray->o[i] = (float)(lcg(&st) % 2000u) / 100.0f - 10.0f;
        float v[3];
        for (int i = 0; i < 3; i++)
            v[i] = (float)(lcg(&st) % 2001u) / 1000.0f - 1.0f;
        if (lcg(&st) % 4u == 0u) v[lcg(&st) % 3u] = 0.0f;  /* parallel slab */
        if (lcg(&st) % 7u == 0u && r > 0) {                 /* duplicate ray */
            float od[3];
            memcpy(ray->o, a[r - 1].o, sizeof(od));
            memcpy(od, a[r - 1].d, sizeof(od));
            v[0] = od[0]; v[1] = od[1]; v[2] = od[2];
        }
        float len = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (len < 1e-6f) { v[0] = 1.0f; v[1] = 0.0f; v[2] = 0.0f; len = 1.0f; }
        ray->d[0] = v[0] / len; ray->d[1] = v[1] / len; ray->d[2] = v[2] / len;
        ray->pixel = r; ray->hit = -1; ray->t = INFINITY;
    }

    trace_scalar(&t, a, nrays, &(Stats){ 0 });

    memcpy(b, a, sizeof(Ray) * (size_t)nrays);
    group_by_direction(b, nrays, order);
    trace_packets(&t, b, nrays, order, &(Stats){ 0 });
    if (!results_equal(a, b, nrays)) return 0;

    memcpy(b, a, sizeof(Ray) * (size_t)nrays);
    shuffle_rays(b, nrays, st ^ 0xdeadbeefu);
    group_by_direction(b, nrays, order);
    trace_packets(&t, b, nrays, order, &(Stats){ 0 });
    if (!results_equal(a, b, nrays)) return 0;

    return 1;
}

/* ---------- checks ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static uint32_t fbits(float f) { uint32_t u; memcpy(&u, &f, 4); return u; }

/* Compare two result sets keyed by pixel id: id and raw t bits must match. */
static int results_equal(const Ray* a, const Ray* b, int n) {
    for (int px = 0; px < n; px++) {
        const Ray* x = NULL; const Ray* y = NULL;
        for (int i = 0; i < n; i++) {
            if (a[i].pixel == px) x = &a[i];
            if (b[i].pixel == px) y = &b[i];
        }
        if (!x || !y || x->hit != y->hit || fbits(x->t) != fbits(y->t))
            return 0;
    }
    return 1;
}

static void shuffle_rays(Ray* rays, int n, unsigned seed) {
    unsigned state = seed;
    for (int i = n - 1; i > 0; i--) {
        int j = (int)(lcg(&state) % (unsigned)(i + 1));
        Ray tmp = rays[i]; rays[i] = rays[j]; rays[j] = tmp;
    }
}

static void reverse_rays(Ray* rays, int n) {
    for (int i = 0, j = n - 1; i < j; i++, j--) {
        Ray tmp = rays[i]; rays[i] = rays[j]; rays[j] = tmp;
    }
}

static int packets_are_bin_pure(const Ray* rays, const int* order, int nrays) {
    int pos = 0;
    while (pos < nrays) {
        int bin = dir_bin(rays[order[pos]].d);
        int k = 0;
        while (pos < nrays && k < PACKET_K &&
               dir_bin(rays[order[pos]].d) == bin) { pos++; k++; }
        if (k == 0) return 0;
    }
    return 1;
}

static void print_stats(const char* label, const Stats* st, int nrays) {
    printf("  %-10s fetches=%6ld (%6.2f/ray) box_preds=%6ld (%5.2f/ray) "
           "prim_tests=%6ld (%5.2f/ray)\n",
           label, st->node_fetches, (double)st->node_fetches / nrays,
           st->ray_box_tests, (double)st->ray_box_tests / nrays,
           st->prim_tests, (double)st->prim_tests / nrays);
}

static long pct(long base, long v) {
    return base == 0 ? 0 : (base - v) * 100L / base;
}

int main(void) {
    Tree t;
    build_grid_scene(&t);
    check("prim_count", t.nprim == 128);
    check("build_ok", build(&t, 0, t.nprim) == 0 && t.nnode > 2);

    Ray ref[MAX_RAYS], run[MAX_RAYS];
    int order[MAX_RAYS];

    /* ---- coherent camera rays ---- */
    int ncoh = make_coherent_rays(ref, 16);
    memcpy(run, ref, sizeof(Ray) * (size_t)ncoh);

    Stats sc = { 0 }, pc = { 0 };
    trace_scalar(&t, ref, ncoh, &sc);
    group_by_direction(run, ncoh, order);
    check("coherent_bin_pure", packets_are_bin_pure(run, order, ncoh));
    trace_packets(&t, run, ncoh, order, &pc);
    check("coherent_exact_match", results_equal(ref, run, ncoh));

    printf("coherent rays=%d (packets of %d):\n", ncoh, PACKET_K);
    print_stats("scalar", &sc, ncoh);
    print_stats("packet", &pc, ncoh);
    printf("  reduction: fetches -%ld%%, box_preds -%ld%%, prim_tests -%ld%%\n",
           pct(sc.node_fetches, pc.node_fetches),
           pct(sc.ray_box_tests, pc.ray_box_tests),
           pct(sc.prim_tests, pc.prim_tests));
    check("coherent_fetches_drop", pc.node_fetches < sc.node_fetches);
    check("coherent_preds_no_increase", pc.ray_box_tests <= sc.ray_box_tests);
    check("coherent_prim_no_increase", pc.prim_tests <= sc.prim_tests);

    /* shuffled input order: regroup, retrace, per-pixel must be identical */
    memcpy(run, ref, sizeof(Ray) * (size_t)ncoh);
    shuffle_rays(run, ncoh, 987654321u);
    group_by_direction(run, ncoh, order);
    Stats ps = { 0 };
    trace_packets(&t, run, ncoh, order, &ps);
    check("shuffle_exact_match", results_equal(ref, run, ncoh));
    check("shuffle_fetches_same", ps.node_fetches == pc.node_fetches);

    /* reversed input order */
    memcpy(run, ref, sizeof(Ray) * (size_t)ncoh);
    reverse_rays(run, ncoh);
    group_by_direction(run, ncoh, order);
    Stats pr = { 0 };
    trace_packets(&t, run, ncoh, order, &pr);
    check("reverse_exact_match", results_equal(ref, run, ncoh));
    check("reverse_fetches_same", pr.node_fetches == pc.node_fetches);

    /* ---- incoherent control rays: correctness still exact, gains smaller ---- */
    int ninc = make_incoherent_rays(ref, 256, 424242u);
    memcpy(run, ref, sizeof(Ray) * (size_t)ninc);
    Stats si = { 0 }, pi = { 0 };
    trace_scalar(&t, ref, ninc, &si);
    group_by_direction(run, ninc, order);
    trace_packets(&t, run, ninc, order, &pi);
    check("incoherent_exact_match", results_equal(ref, run, ninc));

    printf("incoherent rays=%d (packets of %d):\n", ninc, PACKET_K);
    print_stats("scalar", &si, ninc);
    print_stats("packet", &pi, ninc);
    printf("  reduction: fetches -%ld%%, box_preds -%ld%%, prim_tests -%ld%%\n",
           pct(si.node_fetches, pi.node_fetches),
           pct(si.ray_box_tests, pi.ray_box_tests),
           pct(si.prim_tests, pi.prim_tests));
    check("incoherent_fetches_drop", pi.node_fetches < si.node_fetches);
    check("incoherent_preds_no_increase", pi.ray_box_tests <= si.ray_box_tests);

    memcpy(run, ref, sizeof(Ray) * (size_t)ninc);
    shuffle_rays(run, ninc, 13579u);
    group_by_direction(run, ninc, order);
    Stats pish = { 0 };
    trace_packets(&t, run, ninc, order, &pish);
    check("incoherent_shuffle_match", results_equal(ref, run, ninc));
    check("incoherent_shuffle_fetches_same", pish.node_fetches == pi.node_fetches);

    /* miss coverage: a ray aimed away from the scene must report no hit */
    Ray miss[1] = {{ { 50.0f, 50.0f, 50.0f }, { 0, 1, 0 }, 0, -1, INFINITY }};
    Stats pm = { 0 };
    int mo[1] = { 0 };
    trace_packets(&t, miss, 1, mo, &pm);
    check("miss_stays_miss", miss[0].hit == -1);

    /* randomized stress: many scenes x random rays, bit-exact every time */
    int stress_ok = 1;
    for (unsigned seed = 1; seed <= 20u; seed++)
        stress_ok &= stress_scene(seed * 7919u + 13u);
    check("stress_20_scenes_bitexact", stress_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
