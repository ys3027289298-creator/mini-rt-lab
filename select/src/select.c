
/* select: pick an acceleration structure (brute / uniform grid / bvh) and its
   parameter from geometry distribution + query pattern, then validate the
   heuristic against an exhaustive sweep of all candidates. Reports the
   deviation bound vs the exhaustive optimum and lists the scenarios where
   the criterion fails. Empty scenes/queries are rejected.

   Cost unit is a deterministic proxy for time: node visit = C_NODE,
   primitive test = C_PRIM. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_PRIM 64
#define MAX_NODE (2 * MAX_PRIM)
#define MAX_RAY 256
#define MAX_GRID_RES 16
#define MAX_CELL (MAX_GRID_RES * MAX_GRID_RES * MAX_GRID_RES)
#define MAX_ENTRY 8192

#define C_NODE 1.0
#define C_PRIM 2.0
#define FAIL_RATIO 1.5   /* heuristic/optimal above this counts as criterion failure */

typedef struct { float min[3], max[3]; } AABB;
typedef struct { AABB box[MAX_PRIM]; int n; } Scene;
typedef struct { float o[3], d[3]; } Ray;
typedef struct { Ray r[MAX_RAY]; int n, coherent; } Query;
typedef struct { int kind, param; } Choice; /* kind: 0 brute, 1 grid(res), 2 bvh(leaf) */

/* ---- rng ---- */
static unsigned long long rs = 88172645463325252ULL;
static void rseed(unsigned long long s) { rs = s ? s : 1; }
static unsigned long long rnext(void) { rs ^= rs << 13; rs >>= 7; rs ^= rs << 17; return rs; }
static float rf(void) { return (float)(rnext() % 0x1000000) / (float)0x1000000; }

/* ---- aabb ---- */
static void aabb_empty(AABB* b) {
    for (int i = 0; i < 3; i++) { b->min[i] = 1e30f; b->max[i] = -1e30f; }
}
static void aabb_grow(AABB* b, const AABB* o) {
    for (int i = 0; i < 3; i++) {
        if (o->min[i] < b->min[i]) b->min[i] = o->min[i];
        if (o->max[i] > b->max[i]) b->max[i] = o->max[i];
    }
}
static float aabb_area(const AABB* b) {
    float dx = b->max[0] - b->min[0];
    float dy = b->max[1] - b->min[1];
    float dz = b->max[2] - b->min[2];
    if (dx < 0 || dy < 0 || dz < 0) return 0.0f;
    return 2.0f * (dx * dy + dy * dz + dz * dx);
}
static int aabb_hit(const AABB* b, const Ray* r, float* t0p, float* t1p) {
    float t0 = 0.0f, t1 = 1e30f;
    for (int i = 0; i < 3; i++) {
        if (fabsf(r->d[i]) < 1e-8f) {
            if (r->o[i] < b->min[i] || r->o[i] > b->max[i]) return 0;
        } else {
            float inv = 1.0f / r->d[i];
            float ta = (b->min[i] - r->o[i]) * inv, tb = (b->max[i] - r->o[i]) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    *t0p = t0; *t1p = t1;
    return 1;
}

/* ---- scene / query generation (deterministic) ---- */
static void gen_scene(Scene* s, int kind, int n, unsigned long long seed) {
    rseed(seed);
    s->n = n;
    for (int i = 0; i < n; i++) {
        float c[3], h;
        if (kind == 0) {                    /* uniform in unit cube */
            c[0] = rf(); c[1] = rf(); c[2] = rf();
            h = 0.02f + 0.06f * rf();
        } else if (kind == 1) {             /* clustered around 3 centers */
            static const float cc[3][3] = { {0.25f,0.25f,0.25f}, {0.75f,0.6f,0.3f}, {0.4f,0.7f,0.8f} };
            int k = (int)(rf() * 3.0f); if (k > 2) k = 2;
            for (int a = 0; a < 3; a++) c[a] = cc[k][a] + (rf() - 0.5f) * 0.2f;
            h = 0.02f + 0.03f * rf();
        } else if (kind == 2) {             /* thin line along the diagonal */
            float t = (float)i / (float)(n > 1 ? n - 1 : 1);
            for (int a = 0; a < 3; a++) c[a] = 0.1f + 0.8f * t + (rf() - 0.5f) * 0.04f;
            h = 0.02f;
        } else if (kind == 3) {             /* one huge prim + many tiny */
            if (i == 0) { c[0] = c[1] = c[2] = 0.5f; h = 0.45f; }
            else { c[0] = rf(); c[1] = rf(); c[2] = rf(); h = 0.01f; }
        } else if (kind == 4) {             /* bimodal: half big, half tiny */
            c[0] = rf(); c[1] = rf(); c[2] = rf();
            h = (i % 2 == 0) ? 0.15f : 0.01f;
        } else if (kind == 5) {             /* outlier: dense corner + sparse rest */
            if (i < n * 7 / 8)
                for (int a = 0; a < 3; a++) c[a] = 0.05f + rf() * 0.15f;
            else
                for (int a = 0; a < 3; a++) c[a] = rf();
            h = 0.02f;
        } else if (kind == 6) {             /* stacked: one lone + one deep stack */
            if (i == 0) { c[0] = c[1] = c[2] = 0.1f; }
            else {
                for (int a = 0; a < 3; a++) c[a] = 0.8f + (rf() - 0.5f) * 0.01f;
            }
            h = 0.02f;
        } else {                            /* giants: 2 huge hidden among tiny */
            if (i < 2) { c[0] = rf(); c[1] = rf(); c[2] = rf(); h = 0.4f; }
            else { c[0] = rf(); c[1] = rf(); c[2] = rf(); h = 0.005f; }
        }
        for (int a = 0; a < 3; a++) { s->box[i].min[a] = c[a] - h; s->box[i].max[a] = c[a] + h; }
    }
}

static void gen_query(Query* q, const Scene* s, int n, int coherent, unsigned long long seed) {
    rseed(seed);
    AABB w; aabb_empty(&w);
    for (int i = 0; i < s->n; i++) aabb_grow(&w, &s->box[i]);
    float ext[3], diag = 0.0f;
    for (int a = 0; a < 3; a++) { ext[a] = w.max[a] - w.min[a]; diag += ext[a] * ext[a]; }
    diag = sqrtf(diag) + 1e-6f;
    q->n = n; q->coherent = coherent;
    if (coherent == 1) {                    /* parallel +z rays on a regular grid */
        int side = (int)ceilf(sqrtf((float)n));
        for (int i = 0; i < n; i++) {
            float u = ((float)(i % side) + 0.5f) / (float)side;
            float v = ((float)(i / side) + 0.5f) / (float)side;
            q->r[i].o[0] = w.min[0] + u * ext[0];
            q->r[i].o[1] = w.min[1] + v * ext[1];
            q->r[i].o[2] = w.min[2] - diag;
            q->r[i].d[0] = 0.0f; q->r[i].d[1] = 0.0f; q->r[i].d[2] = 1.0f;
        }
    } else {                                /* random rays aimed at the scene */
        for (int i = 0; i < n; i++) {
            float t[3];
            if (coherent == 2)              /* spot: all rays at one hot prim */
                for (int a = 0; a < 3; a++)
                    t[a] = (s->box[s->n - 1].min[a] + s->box[s->n - 1].max[a]) * 0.5f
                         + (rf() - 0.5f) * 0.05f * ext[a];
            else
                for (int a = 0; a < 3; a++) t[a] = w.min[a] + rf() * ext[a];
            float x, y, z, l2;
            do { x = 2*rf()-1; y = 2*rf()-1; z = 2*rf()-1; l2 = x*x + y*y + z*z; }
            while (l2 > 1.0f || l2 < 1e-6f);
            float l = sqrtf(l2);
            for (int a = 0; a < 3; a++) {
                float oc = (w.min[a] + w.max[a]) * 0.5f;
                float off = (a == 0 ? x : a == 1 ? y : z) / l * diag;
                q->r[i].o[a] = oc + off;
                q->r[i].d[a] = t[a] - q->r[i].o[a];
            }
            float dl = sqrtf(q->r[i].d[0]*q->r[i].d[0] + q->r[i].d[1]*q->r[i].d[1] + q->r[i].d[2]*q->r[i].d[2]);
            for (int a = 0; a < 3; a++) q->r[i].d[a] /= dl;
        }
    }
}

/* ---- brute force ---- */
static double brute_query(const Scene* s, const Ray* r, int* id) {
    double cost = 0.0; float best = 1e30f; *id = -1;
    for (int i = 0; i < s->n; i++) {
        float t0, t1;
        cost += C_PRIM;
        if (aabb_hit(&s->box[i], r, &t0, &t1) && t0 < best) { best = t0; *id = i; }
    }
    return cost;
}

/* ---- uniform grid ---- */
typedef struct {
    int res, ncell, nentry, ok;
    int start[MAX_CELL + 1];
    int entry[MAX_ENTRY];
    AABB bounds;
    float cell[3];
} Grid;

static void grid_build(Grid* g, const Scene* s, int res) {
    memset(g, 0, sizeof(*g));
    g->res = res; g->ncell = res * res * res; g->ok = 1;
    aabb_empty(&g->bounds);
    for (int i = 0; i < s->n; i++) aabb_grow(&g->bounds, &s->box[i]);
    for (int a = 0; a < 3; a++) {
        g->cell[a] = (g->bounds.max[a] - g->bounds.min[a]) / (float)res;
        if (g->cell[a] < 1e-9f) g->cell[a] = 1e-9f;
    }
    for (int i = 0; i < s->n; i++) {        /* count entries per cell */
        int lo[3], hi[3];
        for (int a = 0; a < 3; a++) {
            lo[a] = (int)((s->box[i].min[a] - g->bounds.min[a]) / g->cell[a]);
            hi[a] = (int)((s->box[i].max[a] - g->bounds.min[a]) / g->cell[a]);
            if (lo[a] < 0) lo[a] = 0;
            if (hi[a] >= res) hi[a] = res - 1;
        }
        for (int z = lo[2]; z <= hi[2]; z++)
            for (int y = lo[1]; y <= hi[1]; y++)
                for (int x = lo[0]; x <= hi[0]; x++) {
                    g->start[x + res * (y + res * z) + 1]++;
                    g->nentry++;
                }
    }
    if (g->nentry > MAX_ENTRY) { g->ok = 0; return; }
    for (int c = 0; c < g->ncell; c++) g->start[c + 1] += g->start[c];
    int cur[MAX_CELL];
    memcpy(cur, g->start, sizeof(cur));
    for (int i = 0; i < s->n; i++) {        /* fill */
        int lo[3], hi[3];
        for (int a = 0; a < 3; a++) {
            lo[a] = (int)((s->box[i].min[a] - g->bounds.min[a]) / g->cell[a]);
            hi[a] = (int)((s->box[i].max[a] - g->bounds.min[a]) / g->cell[a]);
            if (lo[a] < 0) lo[a] = 0;
            if (hi[a] >= res) hi[a] = res - 1;
        }
        for (int z = lo[2]; z <= hi[2]; z++)
            for (int y = lo[1]; y <= hi[1]; y++)
                for (int x = lo[0]; x <= hi[0]; x++)
                    g->entry[cur[x + res * (y + res * z)]++] = i;
    }
}

static double grid_query(const Grid* g, const Scene* s, const Ray* r, int* id) {
    *id = -1;
    float t0, t1;
    if (!aabb_hit(&g->bounds, r, &t0, &t1)) return C_NODE;
    double cost = C_NODE;
    int ix[3], step[3]; float tn[3], td[3];
    for (int a = 0; a < 3; a++) {
        float p = r->o[a] + t0 * r->d[a] - g->bounds.min[a];
        int i = (int)(p / g->cell[a]);
        if (i < 0) i = 0;
        if (i >= g->res) i = g->res - 1;
        ix[a] = i;
        if (fabsf(r->d[a]) < 1e-8f) { step[a] = 0; tn[a] = 1e30f; td[a] = 1e30f; }
        else if (r->d[a] > 0.0f) {
            step[a] = 1; td[a] = g->cell[a] / r->d[a];
            tn[a] = t0 + ((float)(i + 1) * g->cell[a] - p) / r->d[a];
        } else {
            step[a] = -1; td[a] = -g->cell[a] / r->d[a];
            tn[a] = t0 + ((float)i * g->cell[a] - p) / r->d[a];
        }
    }
    float best = 1e30f;
    for (;;) {                              /* 3D DDA walk */
        int c = ix[0] + g->res * (ix[1] + g->res * ix[2]);
        cost += C_NODE;
        for (int e = g->start[c]; e < g->start[c + 1]; e++) {
            float a0, a1;
            cost += C_PRIM;
            if (aabb_hit(&s->box[g->entry[e]], r, &a0, &a1) && a0 < best) { best = a0; *id = g->entry[e]; }
        }
        int a = 0;
        if (tn[1] < tn[a]) a = 1;
        if (tn[2] < tn[a]) a = 2;
        if (tn[a] > t1 || tn[a] > best) break;
        ix[a] += step[a];
        if (ix[a] < 0 || ix[a] >= g->res) break;
        tn[a] += td[a];
    }
    return cost;
}

/* ---- bvh (median split, leaf size parameter) ---- */
typedef struct { AABB box; int left, right, first, count; } BNode;
typedef struct { BNode n[MAX_NODE]; int idx[MAX_PRIM]; int nnode, leaf; } Bvh;

static float prim_c(const Scene* s, int i, int a) {
    return (s->box[i].min[a] + s->box[i].max[a]) * 0.5f;
}

static int bvh_build_rec(Bvh* b, const Scene* s, int first, int count) {
    if (b->nnode >= MAX_NODE) return -1;
    int ni = b->nnode++;
    BNode* nd = &b->n[ni];
    aabb_empty(&nd->box);
    for (int i = first; i < first + count; i++) aabb_grow(&nd->box, &s->box[b->idx[i]]);
    if (count <= b->leaf) {
        nd->left = nd->right = -1; nd->first = first; nd->count = count;
        return ni;
    }
    float w[3];
    for (int a = 0; a < 3; a++) w[a] = nd->box.max[a] - nd->box.min[a];
    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);
    float mid_c = nd->box.min[axis] + w[axis] * 0.5f;
    int mid = first;
    for (int i = first; i < first + count; i++) {
        if (prim_c(s, b->idx[i], axis) < mid_c) {
            int tmp = b->idx[i]; b->idx[i] = b->idx[mid]; b->idx[mid] = tmp;
            mid++;
        }
    }
    if (mid == first || mid == first + count) mid = first + count / 2;
    nd->first = -1; nd->count = 0;
    nd->left = bvh_build_rec(b, s, first, mid - first);
    nd->right = bvh_build_rec(b, s, mid, first + count - mid);
    return ni;
}

static void bvh_build(Bvh* b, const Scene* s, int leaf) {
    memset(b, 0, sizeof(*b));
    b->leaf = leaf;
    for (int i = 0; i < s->n; i++) b->idx[i] = i;
    bvh_build_rec(b, s, 0, s->n);
}

static double bvh_query(const Bvh* b, const Scene* s, const Ray* r, int* id) {
    *id = -1;
    double cost = 0.0; float best = 1e30f;
    int stack[MAX_NODE], sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        const BNode* nd = &b->n[stack[--sp]];
        float t0, t1;
        cost += C_NODE;
        if (!aabb_hit(&nd->box, r, &t0, &t1) || t0 > best) continue;
        if (nd->left < 0) {
            for (int i = nd->first; i < nd->first + nd->count; i++) {
                float a0, a1;
                cost += C_PRIM;
                if (aabb_hit(&s->box[b->idx[i]], r, &a0, &a1) && a0 < best) { best = a0; *id = b->idx[i]; }
            }
        } else {
            stack[sp++] = nd->left;
            stack[sp++] = nd->right;
        }
    }
    return cost;
}

/* ---- candidate sweep ---- */
static const Choice CANDS[] = {
    {0, 0}, {1, 4}, {1, 8}, {1, 16}, {2, 1}, {2, 2}, {2, 4}
};
#define NCAND 7

static void choice_name(Choice ch, char* buf) {
    if (ch.kind == 0) sprintf(buf, "brute");
    else if (ch.kind == 1) sprintf(buf, "grid%d", ch.param);
    else sprintf(buf, "bvh%d", ch.param);
}

/* measured cost of one candidate over the whole query set; hits may be NULL */
static double run_candidate(const Scene* s, const Query* q, Choice ch, int* hits) {
    double cost = 0.0;
    if (ch.kind == 0) {
        for (int i = 0; i < q->n; i++) cost += brute_query(s, &q->r[i], hits ? &hits[i] : &(int){0});
    } else if (ch.kind == 1) {
        static Grid g;
        grid_build(&g, s, ch.param);
        if (!g.ok) return 1e30;
        cost = (g.ncell + g.nentry) * C_NODE;
        for (int i = 0; i < q->n; i++) cost += grid_query(&g, s, &q->r[i], hits ? &hits[i] : &(int){0});
    } else {
        Bvh b;
        bvh_build(&b, s, ch.param);
        cost = b.nnode * 4.0 * C_NODE;
        for (int i = 0; i < q->n; i++) cost += bvh_query(&b, s, &q->r[i], hits ? &hits[i] : &(int){0});
    }
    return cost;
}

/* ---- heuristic criterion: estimate cost from scene stats + query pattern ---- */
typedef struct {
    int n, q, coherent;
    float ext[3];    /* scene extent per axis */
    float msize[3];  /* mean prim extent per axis */
    float rho;       /* child-box overlap of a root median split: (SA(L)+SA(R))/SA(root) */
    float hit_frac;  /* fraction of query rays hitting the root box */
} Stats;

static void scene_stats(const Scene* s, const Query* q, Stats* st) {
    AABB w; aabb_empty(&w);
    for (int i = 0; i < s->n; i++) aabb_grow(&w, &s->box[i]);
    st->n = s->n; st->q = q->n; st->coherent = q->coherent;
    for (int a = 0; a < 3; a++) {
        st->ext[a] = w.max[a] - w.min[a];
        if (st->ext[a] < 1e-6f) st->ext[a] = 1e-6f;
        st->msize[a] = 0.0f;
        for (int i = 0; i < s->n; i++) st->msize[a] += s->box[i].max[a] - s->box[i].min[a];
        st->msize[a] /= (float)s->n;
    }
    /* root median split on the widest axis, child surface-area ratio */
    {
        int axis = 0;
        for (int a = 1; a < 3; a++) if (st->ext[a] > st->ext[axis]) axis = a;
        float mid_c = w.min[axis] + st->ext[axis] * 0.5f;
        AABB bl, br; aabb_empty(&bl); aabb_empty(&br);
        for (int i = 0; i < s->n; i++) {
            float c = (s->box[i].min[axis] + s->box[i].max[axis]) * 0.5f;
            aabb_grow(c < mid_c ? &bl : &br, &s->box[i]);
        }
        float root_sa = aabb_area(&w);
        st->rho = root_sa > 0.0f ? (aabb_area(&bl) + aabb_area(&br)) / root_sa : 1.0f;
    }
    /* query pattern: fraction of rays that reach the scene at all */
    {
        int hit = 0;
        for (int i = 0; i < q->n; i++) {
            float t0, t1;
            if (aabb_hit(&w, &q->r[i], &t0, &t1)) hit++;
        }
        st->hit_frac = (float)hit / (float)q->n;
    }
}

static double est_cost(const Stats* st, Choice ch) {
    if (ch.kind == 0) return (double)st->q * st->n * C_PRIM;
    if (ch.kind == 2) {
        double d = ceil(log2((double)st->n / ch.param));
        if (d < 1.0) d = 1.0;
        /* expected traversal: rays missing the root pay one node test; rays
           entering branch with child-box ratio rho per level (geometric) */
        double visits = 0.0, reach = 1.0;
        for (int i = 0; i < (int)d; i++) { reach *= st->rho; visits += reach; }
        double per_ray = C_NODE + st->hit_frac * (visits * C_NODE + reach * ch.param * C_PRIM);
        double nodes = 2.0 * st->n / ch.param;  /* full binary tree over n/leaf leaves */
        return nodes * 4.0 * C_NODE + st->q * per_ray;
    }
    int r = ch.param;
    double cells = (double)r * r * r;
    double dup = 1.0;                       /* mean cells one prim overlaps */
    for (int a = 0; a < 3; a++) {
        double cs = st->ext[a] / r;
        dup *= 1.0 + st->msize[a] / (cs > 1e-9 ? cs : 1e-9);
    }
    double path = st->coherent == 1 ? r : 1.5 * r;   /* cells a ray crosses */
    return (cells + st->n * dup) * C_NODE
         + st->q * (C_NODE + st->hit_frac * (path * C_NODE + (st->n * dup / cells) * path * 0.5 * C_PRIM));
}

/* returns -1 on empty scene/query, 0 and the chosen candidate otherwise */
static int select_structure(const Scene* s, const Query* q, Choice* out) {
    if (!s || s->n <= 0 || !q || q->n <= 0) return -1;
    Stats st;
    scene_stats(s, q, &st);
    double bestc = 1e300;
    int bi = 0;
    for (int i = 0; i < NCAND; i++) {
        double c = est_cost(&st, CANDS[i]);
        if (c < bestc) { bestc = c; bi = i; }
    }
    *out = CANDS[bi];
    return 0;
}

/* ---- self test ---- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Scene s;
    Query q;
    Choice ch;

    gen_scene(&s, 0, 32, 7);
    gen_query(&q, &s, 32, 0, 8);

    /* empty input must be rejected */
    Scene es; es.n = 0;
    Query eq; eq.n = 0;
    check("reject_empty_scene", select_structure(&es, &q, &ch) == -1);
    check("reject_empty_query", select_structure(&s, &eq, &ch) == -1);
    check("accept_valid", select_structure(&s, &q, &ch) == 0);

    /* all structures must agree on the hit result */
    {
        int ref[MAX_RAY], got[MAX_RAY];
        int consistent = 1;
        run_candidate(&s, &q, CANDS[0], ref);
        for (int c = 1; c < NCAND; c++) {
            double co = run_candidate(&s, &q, CANDS[c], got);
            if (co >= 1e29) continue;       /* infeasible candidate */
            for (int i = 0; i < q.n; i++)
                if (got[i] != ref[i]) consistent = 0;
        }
        check("hit_consistent", consistent);
    }

    /* criterion sanity: few queries -> brute, many queries -> an index */
    {
        Query one; gen_query(&one, &s, 1, 0, 11);
        check("few_queries_brute", select_structure(&s, &one, &ch) == 0 && ch.kind == 0);
        Scene big; gen_scene(&big, 0, 48, 12);
        Query many; gen_query(&many, &big, 128, 0, 13);
        check("many_queries_indexed", select_structure(&big, &many, &ch) == 0 && ch.kind != 0);
    }

    /* scenario suite: heuristic pick vs exhaustive optimum */
    {
        static const char* skind[8] = { "uniform", "cluster", "line", "mixed", "bimodal", "outlier", "stacked", "giants" };
        static const char* qkind[3] = { "rnd", "coh", "spot" };
        static const int ns[] = { 16, 48 };
        static const int qs[] = { 8, 128 };
        double bound = 1.0;
        int nfail = 0, nscn = 0;
        for (int k = 0; k < 8; k++)
            for (int ni = 0; ni < 2; ni++)
                for (int qi = 0; qi < 2; qi++)
                    for (int coh = 0; coh < 3; coh++) {
                        Scene sc; Query qu; Choice pick;
                        gen_scene(&sc, k, ns[ni], 100 + k * 7 + ni * 3 + coh);
                        gen_query(&qu, &sc, qs[qi], coh, 200 + k * 5 + qi);
                        if (select_structure(&sc, &qu, &pick) != 0) continue;
                        double opt = 1e300, mine = 1e300;
                        for (int c = 0; c < NCAND; c++) {
                            double co = run_candidate(&sc, &qu, CANDS[c], NULL);
                            if (co < opt) opt = co;
                            if (CANDS[c].kind == pick.kind && CANDS[c].param == pick.param) mine = co;
#ifdef DEBUG_COSTS
                            { char cn[16]; choice_name(CANDS[c], cn); printf("  %s=%.0f", cn, co); }
#endif
                        }
#ifdef DEBUG_COSTS
                        { Stats st; scene_stats(&sc, &qu, &st);
                          for (int c = 0; c < NCAND; c++) printf("  est%d=%.0f", c, est_cost(&st, CANDS[c]));
                          printf("  rho=%.2f hit=%.2f scn=%s_n%d_q%d_%s\n", st.rho, st.hit_frac,
                                 skind[k], ns[ni], qs[qi], qkind[coh]); }
#endif
                        double ratio = mine / opt;
                        if (ratio > bound) bound = ratio;
                        nscn++;
                        if (ratio > FAIL_RATIO) {
                            char pn[16];
                            choice_name(pick, pn);
                            printf("fail=%s_n%d_q%d_%s pick=%s ratio=%.2f\n",
                                   skind[k], ns[ni], qs[qi], qkind[coh], pn, ratio);
                            nfail++;
                        }
                    }
        printf("scenarios=%d bound=%.3f fails=%d\n", nscn, bound, nfail);
        check("bound_valid", bound >= 1.0 && bound < 5.0);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
