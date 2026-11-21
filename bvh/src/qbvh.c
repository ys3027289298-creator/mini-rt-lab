
/* qbvh: quantized narrow-node BVH, compared against the wide binary baseline.
 *
 * Wide layout (baseline, same as bvh.c): one 40 B node per tree node, own
 * box + two child indices, traversal chases one node at a time.
 *
 * Narrow layout: fixed-size groups of up to QN_ARITY children. One parent
 * box per group is the dequantization frame; each child packs a 6-byte
 * quantized box (8-bit offsets into the group box, rounded outward so no
 * hit can be lost) plus a 4-byte ref (leaf flag | prim range, or child
 * group index). A traversal step fetches one whole group (two cache
 * lines) and tests all its children from that single fetch.
 *
 * Self-test reports compression ratio, traversal cost (visits, box tests,
 * ns/ray), modeled 32 KiB L1 misses for both layouts, verifies hit results
 * match the wide tree on random rays, and rejects out-of-bounds packed refs.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define MAX_PRIM   16384
#define MAX_NODE   32767
#define MAX_QNODE  2100
#define QN_ARITY   8
#define LEAF_MAX   4
#define QLEAF      0x80000000u
#define NPRIM      16384
#define RAYS       20000

typedef struct { float min[3], max[3]; } AABB;
typedef struct { AABB box; float cx, cy, cz; int id; } Prim;
typedef struct { AABB box; int left, right; int first, count; } Node;   /* 40 B, wide */

typedef struct {
    uint8_t  qmin[3], qmax[3];   /* child box as 8-bit offsets into group box */
    uint16_t flags;              /* reserved, kept 0 */
    uint32_t ref;                /* QLEAF|(count<<16)|first, else group index */
} QChild;                        /* 12 B */

typedef struct {
    float    pmin[3], pmax[3];   /* group box = dequantization frame */
    QChild   ch[QN_ARITY];
    uint32_t nchild;
    uint32_t pad;
} QNode;                         /* 128 B = two cache lines, one fetch */

typedef struct { Prim p[MAX_PRIM]; Node n[MAX_NODE]; int nprim, nnode; } Tree;
typedef struct { QNode n[MAX_QNODE]; int nnode; } QTree;

static void aabb_empty(AABB* b) {
    for (int i = 0; i < 3; i++) { b->min[i] = 1e30f; b->max[i] = -1e30f; }
}
static void aabb_grow(AABB* b, const AABB* o) {
    for (int i = 0; i < 3; i++) {
        if (o->min[i] < b->min[i]) b->min[i] = o->min[i];
        if (o->max[i] > b->max[i]) b->max[i] = o->max[i];
    }
}
static int aabb_hit(const AABB* b, const float o[3], const float d[3], float tmax, float* tmin) {
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

static int tree_add(Tree* t, float minx, float miny, float minz, float maxx, float maxy, float maxz) {
    if (t->nprim >= MAX_PRIM) return -1;
    Prim* p = &t->p[t->nprim];
    p->box.min[0] = minx; p->box.min[1] = miny; p->box.min[2] = minz;
    p->box.max[0] = maxx; p->box.max[1] = maxy; p->box.max[2] = maxz;
    p->cx = (minx + maxx) * 0.5f; p->cy = (miny + maxy) * 0.5f; p->cz = (minz + maxz) * 0.5f;
    p->id = t->nprim;
    return t->nprim++;
}

/* split range on widest axis at centroid median; returns split position */
static int median_split(Prim* p, int first, int count) {
    AABB box;
    aabb_empty(&box);
    for (int i = first; i < first + count; i++) aabb_grow(&box, &p[i].box);
    float w[3] = { box.max[0] - box.min[0], box.max[1] - box.min[1], box.max[2] - box.min[2] };
    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);
    float c = box.min[axis] + w[axis] * 0.5f;
    int mid = first;
    for (int i = first; i < first + count; i++) {
        float v = (axis == 0) ? p[i].cx : (axis == 1 ? p[i].cy : p[i].cz);
        if (v < c) {
            Prim tmp = p[i]; p[i] = p[mid]; p[mid] = tmp;
            mid++;
        }
    }
    if (mid == first || mid == first + count) mid = first + count / 2;
    return mid;
}

/* wide binary build (same rule as bvh.c) */
static int build_wide(Tree* t, int first, int count) {
    if (t->nnode >= MAX_NODE) return -1;
    int ni = t->nnode++;
    Node* n = &t->n[ni];
    aabb_empty(&n->box);
    for (int i = first; i < first + count; i++) aabb_grow(&n->box, &t->p[i].box);
    if (count <= LEAF_MAX) {
        n->left = n->right = -1;
        n->first = first; n->count = count;
        return ni;
    }
    int mid = median_split(t->p, first, count);
    n->first = -1; n->count = 0;
    n->left = build_wide(t, first, mid - first);
    n->right = build_wide(t, mid, first + count - mid);
    return ni;
}

/* ---- narrow (quantized) layout ---- */

static void quantize(const QNode* qn, const AABB* b, uint8_t* qmin, uint8_t* qmax) {
    for (int k = 0; k < 3; k++) {
        float ext = qn->pmax[k] - qn->pmin[k];
        if (ext <= 0.0f) { qmin[k] = 0; qmax[k] = 0; continue; }
        float inv = 255.0f / ext;
        int lo = (int)floorf((b->min[k] - qn->pmin[k]) * inv);   /* round outward */
        int hi = (int)ceilf((b->max[k] - qn->pmin[k]) * inv);
        if (lo < 0) lo = 0; else if (lo > 255) lo = 255;
        if (hi < 0) hi = 0; else if (hi > 255) hi = 255;
        qmin[k] = (uint8_t)lo;
        qmax[k] = (uint8_t)hi;
    }
}

static void dequant(const QNode* qn, const QChild* c, AABB* b) {
    for (int k = 0; k < 3; k++) {
        float ext = qn->pmax[k] - qn->pmin[k];
        float sc = ext * (1.0f / 255.0f);
        float slack = ext * 1e-5f + 1e-30f;   /* absorb float rounding, stay conservative */
        b->min[k] = qn->pmin[k] + (float)c->qmin[k] * sc - slack;
        b->max[k] = qn->pmin[k] + (float)c->qmax[k] * sc + slack;
    }
}

/* narrow build: sort by centroid on widest axis, split into up to QN_ARITY
 * equal-count chunks so groups and leaves stay full */
static int qn_axis;
static int prim_cmp(const void* a, const void* b) {
    const Prim* x = a;
    const Prim* y = b;
    float u = (qn_axis == 0) ? x->cx : (qn_axis == 1) ? x->cy : x->cz;
    float v = (qn_axis == 0) ? y->cx : (qn_axis == 1) ? y->cy : y->cz;
    return (u > v) - (u < v);
}

static int qbuild(QTree* q, Prim* p, int first, int count) {
    if (q->nnode >= MAX_QNODE) return -1;
    int qi = q->nnode++;
    QNode* qn = &q->n[qi];
    AABB box;
    aabb_empty(&box);
    for (int i = first; i < first + count; i++) aabb_grow(&box, &p[i].box);
    for (int k = 0; k < 3; k++) { qn->pmin[k] = box.min[k]; qn->pmax[k] = box.max[k]; }
    qn->nchild = 0;
    qn->pad = 0;

    if (count <= LEAF_MAX) {
        QChild* qc = &qn->ch[0];
        qc->flags = 0;
        quantize(qn, &box, qc->qmin, qc->qmax);
        qc->ref = QLEAF | ((uint32_t)count << 16) | (uint32_t)first;
        qn->nchild = 1;
        return qi;
    }

    float w[3] = { box.max[0] - box.min[0], box.max[1] - box.min[1], box.max[2] - box.min[2] };
    qn_axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);
    qsort(p + first, (size_t)count, sizeof(Prim), prim_cmp);

    int nchunks = (count + LEAF_MAX - 1) / LEAF_MAX;
    if (nchunks > QN_ARITY) nchunks = QN_ARITY;

    for (int i = 0; i < nchunks; i++) {
        int f = first + (int)((long long)count * i / nchunks);
        int c = first + (int)((long long)count * (i + 1) / nchunks) - f;
        if (c <= 0) continue;
        AABB cb;
        aabb_empty(&cb);
        for (int j = f; j < f + c; j++) aabb_grow(&cb, &p[j].box);
        QChild* qc = &qn->ch[qn->nchild++];
        qc->flags = 0;
        quantize(qn, &cb, qc->qmin, qc->qmax);
        if (c <= LEAF_MAX) {
            qc->ref = QLEAF | ((uint32_t)c << 16) | (uint32_t)f;
        } else {
            int sub = qbuild(q, p, f, c);
            if (sub < 0) return -1;
            qc->ref = (uint32_t)sub;
        }
    }
    return qi;
}

/* reject any packed ref that points outside its array; 0 = ok, else -(group+1) */
static int qvalidate(const QTree* q, int nprim) {
    for (int i = 0; i < q->nnode; i++) {
        const QNode* qn = &q->n[i];
        if (qn->nchild < 1 || qn->nchild > QN_ARITY) return -(i + 1);
        for (int j = 0; j < (int)qn->nchild; j++) {
            uint32_t r = qn->ch[j].ref;
            if (r & QLEAF) {
                int first = (int)(r & 0xffffu);
                int count = (int)((r >> 16) & 0x7fffu);
                if (count < 1 || first + count > nprim) return -(i + 1);
            } else {
                if (r >= (uint32_t)q->nnode) return -(i + 1);
            }
        }
    }
    return 0;
}

/* ---- traversal, with optional cache-line touch instrumentation ---- */

typedef struct {
    uint64_t nodes, box_tests, prim_tests, touches;
    uint64_t misses, node_misses, prim_misses;
} Stats;

#define CACHE_SETS  64           /* 64 sets x 8 ways x 64 B = 32 KiB */
#define CACHE_WAYS  8
#define LINE_SHIFT  6

typedef struct {
    uint64_t tag[CACHE_SETS][CACHE_WAYS];
    uint64_t use[CACHE_SETS][CACHE_WAYS];
} SimCache;

static uint64_t sim_tick;

static void sim_line(SimCache* c, uint64_t line, Stats* s, int is_prim) {
    uint64_t set = line % CACHE_SETS;
    uint64_t tag = line / CACHE_SETS;
    s->touches++;
    for (int w = 0; w < CACHE_WAYS; w++) {
        if (c->use[set][w] && c->tag[set][w] == tag) {
            c->use[set][w] = ++sim_tick;
            return;
        }
    }
    s->misses++;
    if (is_prim) s->prim_misses++; else s->node_misses++;
    int v = 0;
    for (int w = 1; w < CACHE_WAYS; w++)
        if (c->use[set][w] < c->use[set][v]) v = w;
    c->tag[set][v] = tag;
    c->use[set][v] = ++sim_tick;
}

static void sim_touch(SimCache* c, const void* p, size_t n, Stats* s, int is_prim) {
    if (!c) return;
    uintptr_t a = (uintptr_t)p;
    uintptr_t b = a + n - 1;
    for (uint64_t line = a >> LINE_SHIFT; line <= b >> LINE_SHIFT; line++)
        sim_line(c, line, s, is_prim);
}

static int traverse_wide(const Tree* t, const float o[3], const float d[3],
                         SimCache* c, Stats* s, int* out_id, float* out_t) {
    int stack[128];
    int sp = 0;
    float best = 1e30f;
    int bid = -1;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        const Node* n = &t->n[ni];
        sim_touch(c, n, sizeof(*n), s, 0);
        s->nodes++;
        float t0;
        s->box_tests++;
        if (!aabb_hit(&n->box, o, d, best, &t0)) continue;
        if (n->left < 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                sim_touch(c, &t->p[i], sizeof(Prim), s, 1);
                s->prim_tests++;
                float tp;
                if (aabb_hit(&t->p[i].box, o, d, best, &tp) &&
                    (tp < best || (tp == best && t->p[i].id < bid))) {
                    best = tp;
                    bid = t->p[i].id;
                }
            }
        } else {
            if (sp + 2 > (int)(sizeof(stack) / sizeof(stack[0]))) return -1;
            stack[sp++] = n->left;
            stack[sp++] = n->right;
        }
    }
    *out_id = bid;
    *out_t = best;
    return 0;
}

/* returns -1 when a packed ref is out of bounds (corrupt tree rejected) */
static int traverse_narrow(const QTree* q, const Prim* prims, int nprim,
                           const float o[3], const float d[3],
                           SimCache* c, Stats* s, int* out_id, float* out_t) {
    int stack[128];
    int sp = 0;
    float best = 1e30f;
    int bid = -1;
    uint64_t steps = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        if (++steps > (1u << 20)) return -1;   /* cycle guard */
        int qi = stack[--sp];
        if (qi < 0 || qi >= q->nnode) return -1;
        const QNode* qn = &q->n[qi];
        sim_touch(c, qn, sizeof(*qn), s, 0);   /* one whole group per fetch */
        s->nodes++;
        float mn[3], sc[3], sl[3];             /* dequant frame, hoisted per group */
        for (int k = 0; k < 3; k++) {
            float ext = qn->pmax[k] - qn->pmin[k];
            mn[k] = qn->pmin[k];
            sc[k] = ext * (1.0f / 255.0f);
            sl[k] = ext * 1e-5f + 1e-30f;
        }
        for (int i = 0; i < (int)qn->nchild; i++) {
            AABB b;
            for (int k = 0; k < 3; k++) {
                b.min[k] = mn[k] + (float)qn->ch[i].qmin[k] * sc[k] - sl[k];
                b.max[k] = mn[k] + (float)qn->ch[i].qmax[k] * sc[k] + sl[k];
            }
            float t0;
            s->box_tests++;
            if (!aabb_hit(&b, o, d, best, &t0)) continue;
            uint32_t r = qn->ch[i].ref;
            if (r & QLEAF) {
                int first = (int)(r & 0xffffu);
                int count = (int)((r >> 16) & 0x7fffu);
                if (count < 1 || first + count > nprim) return -1;
                for (int j = first; j < first + count; j++) {
                    sim_touch(c, &prims[j], sizeof(Prim), s, 1);
                    s->prim_tests++;
                    float tp;
                    if (aabb_hit(&prims[j].box, o, d, best, &tp) &&
                        (tp < best || (tp == best && prims[j].id < bid))) {
                        best = tp;
                        bid = prims[j].id;
                    }
                }
            } else {
                if (r >= (uint32_t)q->nnode) return -1;
                if (sp + 1 > (int)(sizeof(stack) / sizeof(stack[0]))) return -1;
                stack[sp++] = (int)r;
            }
        }
    }
    *out_id = bid;
    *out_t = best;
    return 0;
}

/* ---- self test ---- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static uint32_t rng_s = 0x12345678u;
static float frand(void) {
    uint32_t x = rng_s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_s = x;
    return (float)(x >> 8) * (1.0f / 16777216.0f);
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(void) {
    check("sizes", sizeof(Node) == 40 && sizeof(QChild) == 12 && sizeof(QNode) == 128);

    /* quantization must round outward: dequantized box contains original */
    {
        QNode u;
        memset(&u, 0, sizeof(u));
        for (int k = 0; k < 3; k++) { u.pmin[k] = 0.0f; u.pmax[k] = 10.0f; }
        int ok = 1;
        for (int it = 0; it < 1000 && ok; it++) {
            AABB b, dq;
            for (int k = 0; k < 3; k++) {
                float a = frand() * 10.0f, e = frand() * (10.0f - a);
                b.min[k] = a; b.max[k] = a + e;
            }
            QChild qc;
            memset(&qc, 0, sizeof(qc));
            quantize(&u, &b, qc.qmin, qc.qmax);
            dequant(&u, &qc, &dq);
            for (int k = 0; k < 3; k++)
                if (dq.min[k] > b.min[k] || dq.max[k] < b.max[k]) ok = 0;
        }
        check("quant_conservative", ok);
    }

    Tree* t = malloc(sizeof(Tree));
    QTree* q = malloc(sizeof(QTree));
    if (!t || !q) { printf("all_ok=0\n"); return 1; }
    memset(t, 0, sizeof(*t));
    memset(q, 0, sizeof(*q));

    for (int i = 0; i < NPRIM; i++) {
        float c[3], h[3];
        for (int k = 0; k < 3; k++) {
            c[k] = -10.0f + 20.0f * frand();
            h[k] = 0.05f + 0.95f * frand();
        }
        tree_add(t, c[0] - h[0], c[1] - h[1], c[2] - h[2],
                 c[0] + h[0], c[1] + h[1], c[2] + h[2]);
    }
    check("prim_count", t->nprim == NPRIM);
    check("build_wide_ok", build_wide(t, 0, t->nprim) == 0);

    /* narrow build reorders prims independently; keep a private copy */
    QTree* q2 = q;
    Prim* qp = malloc(sizeof(Prim) * MAX_PRIM);
    if (!qp) { printf("all_ok=0\n"); return 1; }
    memcpy(qp, t->p, sizeof(Prim) * (size_t)t->nprim);
    q2->nnode = 0;
    check("build_narrow_ok", qbuild(q2, qp, 0, t->nprim) == 0);
    check("validate_ok", qvalidate(q2, t->nprim) == 0);

    size_t wide_bytes = (size_t)t->nnode * sizeof(Node);
    size_t narrow_bytes = (size_t)q2->nnode * sizeof(QNode);
    uint64_t entries = 0;
    for (int i = 0; i < q2->nnode; i++) entries += q2->n[i].nchild;
    printf("wide_nodes=%d\n", t->nnode);
    printf("narrow_groups=%d\n", q2->nnode);
    printf("avg_arity=%.2f\n", (double)entries / (double)q2->nnode);
    printf("wide_bytes=%zu\n", wide_bytes);
    printf("narrow_bytes=%zu\n", narrow_bytes);
    printf("compression=%.2f\n", (double)wide_bytes / (double)narrow_bytes);
    check("compresses", narrow_bytes * 2 <= wide_bytes);

    /* rays: odd ones aimed at a random prim, even ones random */
    float(*ro)[3] = malloc(sizeof(float[3]) * RAYS);
    float(*rd)[3] = malloc(sizeof(float[3]) * RAYS);
    if (!ro || !rd) { printf("all_ok=0\n"); return 1; }
    for (int r = 0; r < RAYS; r++) {
        for (int k = 0; k < 3; k++) ro[r][k] = -15.0f + 30.0f * frand();
        if (r & 1) {
            int pi = (int)(frand() * (float)t->nprim);
            if (pi >= t->nprim) pi = t->nprim - 1;
            const Prim* p = &qp[pi];
            float c[3] = { p->cx, p->cy, p->cz };
            for (int k = 0; k < 3; k++) rd[r][k] = c[k] - ro[r][k];
        } else {
            for (int k = 0; k < 3; k++) rd[r][k] = -1.0f + 2.0f * frand();
        }
    }

    /* consistency + traversal cost + modeled cache misses */
    SimCache* cw = calloc(1, sizeof(SimCache));
    SimCache* cn = calloc(1, sizeof(SimCache));
    if (!cw || !cn) { printf("all_ok=0\n"); return 1; }
    Stats sw, sn;
    memset(&sw, 0, sizeof(sw));
    memset(&sn, 0, sizeof(sn));
    int mism = 0;
    for (int r = 0; r < RAYS; r++) {
        int idw, idn;
        float tw, tn;
        if (traverse_wide(t, ro[r], rd[r], cw, &sw, &idw, &tw) != 0) mism++;
        if (traverse_narrow(q2, qp, t->nprim, ro[r], rd[r], cn, &sn, &idn, &tn) != 0) mism++;
        if (idw != idn || tw != tn) mism++;
    }
    printf("rays=%d\n", RAYS);
    printf("hit_mismatches=%d\n", mism);
    check("hit_consistent", mism == 0);

    printf("wide_nodes_per_ray=%.1f\n", (double)sw.nodes / RAYS);
    printf("narrow_groups_per_ray=%.1f\n", (double)sn.nodes / RAYS);
    printf("wide_box_tests_per_ray=%.1f\n", (double)sw.box_tests / RAYS);
    printf("narrow_box_tests_per_ray=%.1f\n", (double)sn.box_tests / RAYS);
    printf("wide_prim_tests_per_ray=%.2f\n", (double)sw.prim_tests / RAYS);
    printf("narrow_prim_tests_per_ray=%.2f\n", (double)sn.prim_tests / RAYS);
    printf("wide_line_touches=%llu\n", (unsigned long long)sw.touches);
    printf("narrow_line_touches=%llu\n", (unsigned long long)sn.touches);
    printf("wide_node_misses=%llu\n", (unsigned long long)sw.node_misses);
    printf("narrow_node_misses=%llu\n", (unsigned long long)sn.node_misses);
    printf("node_miss_reduction=%.2f\n", (double)sw.node_misses / (double)sn.node_misses);
    printf("wide_prim_misses=%llu\n", (unsigned long long)sw.prim_misses);
    printf("narrow_prim_misses=%llu\n", (unsigned long long)sn.prim_misses);
    printf("wide_misses=%llu\n", (unsigned long long)sw.misses);
    printf("narrow_misses=%llu\n", (unsigned long long)sn.misses);
    printf("miss_reduction=%.2f\n", (double)sw.misses / (double)sn.misses);
    check("fewer_misses", sn.misses < sw.misses);
    check("fewer_node_misses", sn.node_misses < sw.node_misses);

    /* wall-clock traversal cost, instrumentation off */
    {
        Stats z;
        memset(&z, 0, sizeof(z));
        int sink = 0, id;
        float tt;
        double t0 = now_s();
        for (int r = 0; r < RAYS; r++) {
            traverse_wide(t, ro[r], rd[r], NULL, &z, &id, &tt);
            sink += id;
        }
        double t1 = now_s();
        for (int r = 0; r < RAYS; r++) {
            traverse_narrow(q2, qp, t->nprim, ro[r], rd[r], NULL, &z, &id, &tt);
            sink += id;
        }
        double t2 = now_s();
        printf("wide_ns_per_ray=%.0f\n", (t1 - t0) * 1e9 / RAYS);
        printf("narrow_ns_per_ray=%.0f\n", (t2 - t1) * 1e9 / RAYS);
        printf("sink=%d\n", sink);
    }

    /* out-of-bounds packed refs must be rejected, by validate and by traversal */
    {
        Stats z;
        memset(&z, 0, sizeof(z));
        int id;
        float tt;
        AABB cb;
        dequant(&q2->n[0], &q2->n[0].ch[0], &cb);
        float ao[3] = {
            (cb.min[0] + cb.max[0]) * 0.5f - 50.0f,
            (cb.min[1] + cb.max[1]) * 0.5f + 0.3f,
            (cb.min[2] + cb.max[2]) * 0.5f - 0.2f
        };
        float ad[3] = { 50.0f, -0.3f, 0.2f };   /* passes through child 0 box */

        uint32_t saved = q2->n[0].ch[0].ref;

        q2->n[0].ch[0].ref = (uint32_t)q2->nnode + 1000;   /* group index past end */
        check("reject_internal_validate", qvalidate(q2, t->nprim) != 0);
        check("reject_internal_traverse",
              traverse_narrow(q2, qp, t->nprim, ao, ad, NULL, &z, &id, &tt) == -1);

        q2->n[0].ch[0].ref = QLEAF | (2u << 16) | (uint32_t)(t->nprim - 1); /* prim range past end */
        check("reject_leaf_validate", qvalidate(q2, t->nprim) != 0);
        check("reject_leaf_traverse",
              traverse_narrow(q2, qp, t->nprim, ao, ad, NULL, &z, &id, &tt) == -1);

        q2->n[0].ch[0].ref = saved;
        check("validate_restored", qvalidate(q2, t->nprim) == 0);
    }

    free(ro); free(rd); free(qp); free(cw); free(cn); free(q); free(t);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
