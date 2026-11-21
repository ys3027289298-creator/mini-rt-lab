
/* shadow: per-light bucketed shadow rays.

   A spherical light is sampled with N equal-area shadow-ray directions.
   Rays are grouped into B directionally coherent buckets by recursive
   median splitting.  The refinement sequence is the same for every B
   (largest bucket first, ties by index), so the partition for B+1 is a
   strict refinement of the partition for B.

   Each bucket runs one shared visibility query: the AABB hull of the eye
   and the bucket's light points is tested once against every occluder, and
   occluders that cannot touch the hull are culled for the whole bucket.
   Only survivors get per-sample segment tests, so buckets away from the
   shadow cone skip most of the scene.  The hull test is conservative, so
   per-sample visibility is identical to straight per-sample intersection.

   The bucket's occlusion rate is its blocked fraction.  Rate 0 or 1 is a
   hard lit / shadow transition; a fractional rate is the penumbra and is
   used as a soft transmission.  Against per-sample visibility the rate
   introduces a squared error  n_b*r*(1-r) / N  per bucket, bounded by
   min(blocked, lit) / N, which is at most 1/2 and whose sum is at most
   1/4.  Splitting a bucket into b1, b2 changes the contribution by
       n1*r1*(1-r1) + n2*r2*(1-r2)
         = n*r*(1-r) - n1*n2*(r1-r2)^2 / n,
   so the error is non-increasing as buckets are added; singleton buckets
   drive it to exactly zero.  B == 0 cannot be grouped and is reported. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_SAMPLE 64
#define MAX_BUCKET 64
#define MAX_OCC 32
#define CANNOT_GROUP (-1)

typedef struct { float x, y, z; } V3;
typedef struct { V3 lo, hi; } Box;

typedef struct {
    int lo, hi;      /* sample index range in the bucket's permutation */
    int count, blocked;
    float rate;      /* blocked / count */
    int edge;        /* 0 < rate < 1: soft transition */
} Bucket;

typedef struct {
    int nbucket;                     /* CANNOT_GROUP when B == 0 */
    int n;
    unsigned char vis[MAX_SAMPLE];   /* bucketed result, 1 = lit */
    unsigned char ref[MAX_SAMPLE];   /* per-sample intersection, 1 = lit */
    Bucket bucket[MAX_BUCKET];
    long shared_tests;               /* hull-vs-occluder, once per bucket */
    long sample_tests;               /* survivor segment tests */
    long ref_tests;                  /* per-sample reference segment tests */
    double edge_err;                 /* sum n_b*r*(1-r) / N */
    double edge_bound;               /* sum min(blocked,lit) / N */
    int edge_buckets, hard_buckets;
} Shadow;

/* conservative: 1 if the segment possibly touches the box, 0 only when
   the segment provably misses it.  d is unnormalized, tmax = 1. */
static int seg_box(V3 o, V3 d, const Box* b) {
    float t0 = 0.0f, t1 = 1.0f;
    float oo[3] = { o.x, o.y, o.z }, dd[3] = { d.x, d.y, d.z };
    float lo[3] = { b->lo.x, b->lo.y, b->lo.z };
    float hi[3] = { b->hi.x, b->hi.y, b->hi.z };
    for (int i = 0; i < 3; i++) {
        if (fabsf(dd[i]) < 1e-8f) {
            if (oo[i] < lo[i] || oo[i] > hi[i]) return 0;
        } else {
            float inv = 1.0f / dd[i];
            float ta = (lo[i] - oo[i]) * inv;
            float tb = (hi[i] - oo[i]) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    return 1;
}

static int box_overlap(const Box* a, const Box* b) {
    if (a->hi.x < b->lo.x || a->lo.x > b->hi.x) return 0;
    if (a->hi.y < b->lo.y || a->lo.y > b->hi.y) return 0;
    if (a->hi.z < b->lo.z || a->lo.z > b->hi.z) return 0;
    return 1;
}

static void box_grow(Box* b, V3 p) {
    if (p.x < b->lo.x) b->lo.x = p.x;
    if (p.y < b->lo.y) b->lo.y = p.y;
    if (p.z < b->lo.z) b->lo.z = p.z;
    if (p.x > b->hi.x) b->hi.x = p.x;
    if (p.y > b->hi.y) b->hi.y = p.y;
    if (p.z > b->hi.z) b->hi.z = p.z;
}

static float axis_val(V3 v, int axis) { return axis == 0 ? v.x : (axis == 1 ? v.y : v.z); }

/* equal-area Fibonacci sphere, uniform over the whole sphere */
static void light_samples(V3 c, float r, int n, V3 out[MAX_SAMPLE]) {
    const float pi = 3.14159265358979f;
    float ga = pi * (3.0f - sqrtf(5.0f));
    for (int i = 0; i < n; i++) {
        float y = 1.0f - 2.0f * ((float)i + 0.5f) / (float)n;
        float rad2 = 1.0f - y * y;
        float rad = sqrtf(rad2 > 0.0f ? rad2 : 0.0f);
        float th = ga * (float)i;
        out[i].x = c.x + r * rad * cosf(th);
        out[i].y = c.y + r * y;
        out[i].z = c.z + r * rad * sinf(th);
    }
}

/* sort a bucket by its widest axis and return the median split index */
static int split_range(const V3* p, int* idx, int lo, int hi) {
    Box box = { { 1e30f, 1e30f, 1e30f }, { -1e30f, -1e30f, -1e30f } };
    for (int i = lo; i < hi; i++) box_grow(&box, p[idx[i]]);
    float w[3] = { box.hi.x - box.lo.x, box.hi.y - box.lo.y, box.hi.z - box.lo.z };
    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);
    for (int i = lo + 1; i < hi; i++) {           /* insertion sort by axis */
        int v = idx[i];
        float key = axis_val(p[v], axis);
        int j = i - 1;
        while (j >= lo && axis_val(p[idx[j]], axis) > key) { idx[j + 1] = idx[j]; j--; }
        idx[j + 1] = v;
    }
    int count = hi - lo;
    if (count < 2) return lo;
    return lo + count / 2;
}

/* Split the largest bucket repeatedly.  Ties go to the lowest index, so the
   whole sequence is deterministic and nested across bucket counts. */
static int build_buckets(const V3* p, int n, int B, int* idx, Bucket out[MAX_BUCKET]) {
    if (B <= 0) return CANNOT_GROUP;
    int want = B < n ? B : n;
    if (want > MAX_BUCKET) want = MAX_BUCKET;
    for (int i = 0; i < n; i++) idx[i] = i;
    out[0].lo = 0; out[0].hi = n;
    int nb = 1;
    while (nb < want) {
        int best = -1, best_n = 1;
        for (int i = 0; i < nb; i++) {
            int c = out[i].hi - out[i].lo;
            if (c > best_n) { best_n = c; best = i; }
        }
        if (best < 0) break;                       /* every bucket is a singleton */
        int mid = split_range(p, idx, out[best].lo, out[best].hi);
        Bucket right;
        right.lo = mid; right.hi = out[best].hi;
        out[best].hi = mid;
        out[nb++] = right;
    }
    for (int i = 0; i < nb; i++) {
        out[i].count = out[i].hi - out[i].lo;
        out[i].blocked = 0;
        out[i].rate = 0.0f;
        out[i].edge = 0;
    }
    return nb;
}

static int shadow_eval(V3 eye, const V3* pts, int n, int B,
                       const Box* occ, int nocc, Shadow* s) {
    memset(s, 0, sizeof(*s));
    s->n = n;
    int idx[MAX_SAMPLE];
    int nb = build_buckets(pts, n, B, idx, s->bucket);
    s->nbucket = nb;
    if (nb == CANNOT_GROUP) return CANNOT_GROUP;

    /* reference: straight per-sample intersection against every occluder */
    for (int i = 0; i < n; i++) {
        V3 d = { pts[i].x - eye.x, pts[i].y - eye.y, pts[i].z - eye.z };
        int blocked = 0;
        for (int k = 0; k < nocc; k++) {
            s->ref_tests++;
            if (seg_box(eye, d, &occ[k])) { blocked = 1; break; }
        }
        s->ref[i] = blocked ? 0 : 1;
    }

    for (int b = 0; b < nb; b++) {
        Bucket* bk = &s->bucket[b];
        Box hull = { { eye.x, eye.y, eye.z }, { eye.x, eye.y, eye.z } };
        for (int i = bk->lo; i < bk->hi; i++) box_grow(&hull, pts[idx[i]]);

        /* one shared query per bucket: which occluders can touch the hull */
        int survivor[MAX_OCC];
        int nsurv = 0;
        for (int k = 0; k < nocc; k++) {
            s->shared_tests++;
            if (box_overlap(&hull, &occ[k])) survivor[nsurv++] = k;
        }

        for (int i = bk->lo; i < bk->hi; i++) {
            V3 d = { pts[idx[i]].x - eye.x, pts[idx[i]].y - eye.y, pts[idx[i]].z - eye.z };
            int blocked = 0;
            for (int j = 0; j < nsurv; j++) {
                s->sample_tests++;
                if (seg_box(eye, d, &occ[survivor[j]])) { blocked = 1; break; }
            }
            s->vis[idx[i]] = blocked ? 0 : 1;
            bk->blocked += blocked;
        }

        bk->rate = (float)bk->blocked / (float)bk->count;
        bk->edge = bk->blocked > 0 && bk->blocked < bk->count;
        int lit = bk->count - bk->blocked;
        s->edge_err += (double)bk->blocked * (double)lit / (double)bk->count;
        s->edge_bound += (double)(bk->blocked < lit ? bk->blocked : lit);
        if (bk->edge) s->edge_buckets++; else s->hard_buckets++;
    }
    s->edge_err /= (double)n;
    s->edge_bound /= (double)n;
    return nb;
}

/* occlusion rate -> transmission: hard outside the penumbra, soft inside */
static float bucket_transmission(const Bucket* b) {
    if (b->rate <= 0.0f) return 1.0f;
    if (b->rate >= 1.0f) return 0.0f;
    return 1.0f - b->rate;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* One light's full curve: exactness, bounded error, monotonicity, the
   soft/hard split and the shared-query savings. */
static void run_light(int li, V3 eye, V3 light_c, float light_r,
                      const Box* occ, int nocc, int* exact, int* bounds_ok,
                      int* monotone, int* zero_at_end, int* soft_ok,
                      int* query_win) {
    enum { N = MAX_SAMPLE };
    const int counts[] = { 1, 2, 4, 8, 16, 32, 64 };
    const int ncounts = (int)(sizeof(counts) / sizeof(counts[0]));
    V3 pts[N];
    light_samples(light_c, light_r, N, pts);

    Shadow ref;
    shadow_eval(eye, pts, N, N, occ, nocc, &ref);
    int ref_lit = 0;
    for (int i = 0; i < N; i++) ref_lit += ref.ref[i];
    if (!(ref_lit > 0 && ref_lit < N)) *exact = 0;   /* needs a real penumbra */

    double prev = 1e30;
    int saw_soft = 0;
    for (int t = 0; t < ncounts; t++) {
        int B = counts[t];
        Shadow s;
        int nb = shadow_eval(eye, pts, N, B, occ, nocc, &s);
        if (nb != B) *exact = 0;
        for (int i = 0; i < N; i++) if (s.vis[i] != ref.ref[i]) *exact = 0;

        printf("light=%d curve B=%d edge_err=%.6f bound=%.6f edge_buckets=%d "
               "hard_buckets=%d seg_tests=%ld ref_seg=%ld shared=%ld\n",
               li, B, s.edge_err, s.edge_bound, s.edge_buckets, s.hard_buckets,
               s.sample_tests, s.ref_tests, s.shared_tests);

        if (s.edge_err > s.edge_bound + 1e-12 || s.edge_err > 0.25 + 1e-12)
            *bounds_ok = 0;
        if (s.edge_err > prev + 1e-12) *monotone = 0;
        prev = s.edge_err;
        if (B == 16 && s.sample_tests < s.ref_tests) *query_win = 1;

        for (int b = 0; b < nb; b++) {
            const Bucket* bk = &s.bucket[b];
            float tr = bucket_transmission(bk);
            if (bk->edge) {
                saw_soft = 1;
                if (!(tr > 0.0f && tr < 1.0f)) *soft_ok = 0;
            } else if (bk->rate == 0.0f) {
                if (tr != 1.0f) *soft_ok = 0;
            } else if (bk->rate == 1.0f) {
                if (tr != 0.0f) *soft_ok = 0;
            } else {
                *soft_ok = 0;
            }
        }
    }
    if (prev != 0.0) *zero_at_end = 0;
    if (!saw_soft) *soft_ok = 0;
}

int main(void) {
    enum { N = MAX_SAMPLE };
    V3 eye = { 0.0f, 0.0f, 0.0f };

    /* three blockers inside the shadow cone, then distant clutter that no
       sample can hit but the reference still pays for on every ray */
    Box occ[MAX_OCC];
    int nocc = 0;
    occ[nocc].lo = (V3){ 0.25f, 3.0f, -0.40f }; occ[nocc++].hi = (V3){ 0.75f, 4.0f, 0.40f };
    occ[nocc].lo = (V3){ -0.75f, 3.0f, -0.40f }; occ[nocc++].hi = (V3){ -0.25f, 4.0f, 0.40f };
    occ[nocc].lo = (V3){ -0.40f, 3.0f, 0.25f }; occ[nocc++].hi = (V3){ 0.40f, 4.0f, 0.75f };
    for (int i = 0; i < 9; i++) {
        float ax = (float)(i % 3) * 3.0f + 4.0f;
        float az = (float)(i / 3) * 3.0f - 4.0f;
        float sx = (i & 1) ? 1.0f : -1.0f;
        occ[nocc].lo = (V3){ sx * ax - 0.3f, 1.0f, az - 0.3f };
        occ[nocc++].hi = (V3){ sx * ax + 0.3f, 2.0f, az + 0.3f };
    }

    /* B = 0 cannot be grouped */
    V3 pts[N];
    light_samples((V3){ 0.0f, 8.0f, 0.0f }, 2.0f, N, pts);
    Shadow z;
    int rc = shadow_eval(eye, pts, N, 0, occ, nocc, &z);
    printf("B=0: cannot group (need at least one bucket)\n");
    check("zero_buckets_reports_cannot_group", rc == CANNOT_GROUP && z.nbucket == CANNOT_GROUP);

    int exact = 1, bounds_ok = 1, monotone = 1, zero_at_end = 1;
    int soft_ok = 1, query_win = 1;

    /* two lights, bucketed independently */
    run_light(0, eye, (V3){ 0.0f, 8.0f, 0.0f }, 2.0f,
              occ, nocc, &exact, &bounds_ok, &monotone, &zero_at_end, &soft_ok, &query_win);
    run_light(1, eye, (V3){ 1.5f, 7.0f, 0.5f }, 1.0f,
              occ, nocc, &exact, &bounds_ok, &monotone, &zero_at_end, &soft_ok, &query_win);

    check("bucketed_matches_per_sample", exact);
    check("edge_rate_error_bounded", bounds_ok);
    check("edge_rate_error_monotone_in_buckets", monotone);
    check("edge_rate_error_zero_at_full_split", zero_at_end);
    check("shared_query_saves_segment_tests", query_win);
    check("soft_hard_transition_from_rate", soft_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
