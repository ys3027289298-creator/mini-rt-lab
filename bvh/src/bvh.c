
/* bvh: median-split AABB tree over triangle bounds, nearest-hit traversal. */
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_PRIM 32
#define MAX_NODE 64

typedef struct { float min[3], max[3]; } AABB;
typedef struct { AABB box; float cx, cy, cz; int id; } Prim;
typedef struct { AABB box; int left, right; int first, count; } Node;

typedef struct {
    Prim p[MAX_PRIM];
    Node n[MAX_NODE];
    int nprim, nnode;
} Tree;

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
    /* split on widest axis at centroid median */
    float w[3] = { n->box.max[0] - n->box.min[0], n->box.max[1] - n->box.min[1], n->box.max[2] - n->box.min[2] };
    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);
    float c[3] = { n->box.min[0] + w[0] * 0.5f, n->box.min[1] + w[1] * 0.5f, n->box.min[2] + w[2] * 0.5f };
    int mid = first;
    for (int i = first; i < first + count; i++) {
        float v = (axis == 0) ? t->p[i].cx : (axis == 1 ? t->p[i].cy : t->p[i].cz);
        if (v < c[axis]) {
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

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* nearest hit over the tree, iterative pre-order with a fixed stack */
static int nearest_hit(const Tree* t, const float o[3], const float d[3], float* tbest, int* idbest) {
    int stack[MAX_NODE];
    int sp = 0;
    float best = 1e30f;
    int bestid = -1;
    stack[sp++] = 0;
    while (sp > 0) {
        const Node* n = &t->n[stack[--sp]];
        float tm;
        if (!aabb_hit(&n->box, o, d, best, &tm)) continue;
        if (n->left < 0) {
            int i;
            for (i = n->first; i < n->first + n->count; i++) {
                float tp;
                if (aabb_hit(&t->p[i].box, o, d, best, &tp)) { best = tp; bestid = t->p[i].id; }
            }
        } else {
            stack[sp++] = n->left;
            stack[sp++] = n->right;
        }
    }
    *tbest = best;
    *idbest = bestid;
    return bestid >= 0;
}

/* deterministic scene: 24 boxes from a fixed LCG, no transcendentals */
static void scene_fill(Tree* t) {
    unsigned rng = 12345u;
    int i, k;
    for (i = 0; i < 24; i++) {
        float c[3], h[3];
        for (k = 0; k < 3; k++) {
            rng = rng * 1664525u + 1013904223u;
            c[k] = (float)(rng >> 8) * (1.0f / 16777216.0f) * 1.6f - 0.8f;
            rng = rng * 1664525u + 1013904223u;
            h[k] = 0.05f + (float)(rng >> 8) * (1.0f / 16777216.0f) * 0.2f;
        }
        tree_add(t, c[0] - h[0], c[1] - h[1], c[2] - h[2],
                    c[0] + h[0], c[1] + h[1], c[2] + h[2]);
    }
}

/* deterministic render: orthographic camera at z=3, nearest BVH hit shaded
   by distance and tinted by prim id. Only +,-,*,/ and comparisons, so the
   image is bit-identical across compilers, optimization levels and machines. */
static int render(const char* path, int size, int strict) {
    Tree t;
    FILE* f;
    int px, py;
    memset(&t, 0, sizeof(t));
    scene_fill(&t);
    if (build(&t, 0, t.nprim) != 0) { fprintf(stderr, "error: scene build failed\n"); return 3; }
    f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "error: cannot open '%s' for writing\n", path); return 2; }
    fprintf(f, "P6\n%d %d\n255\n", size, size);
    for (py = 0; py < size; py++) {
        for (px = 0; px < size; px++) {
            float u = ((float)px + 0.5f) / (float)size * 2.0f - 1.0f;
            float v = 1.0f - ((float)py + 0.5f) / (float)size * 2.0f;
            float o[3] = { u, v, 3.0f };
            float d[3] = { 0.0f, 0.0f, -1.0f };
            float best;
            int id;
            unsigned char rgb[3] = { 0, 0, 0 };
            if (nearest_hit(&t, o, d, &best, &id)) {
                float shade = 4.0f / (1.0f + best);
                if (strict && !(shade >= 0.0f && shade <= 4.0f)) {
                    fprintf(stderr, "error: strict check failed at pixel (%d,%d): shade out of range\n", px, py);
                    fclose(f);
                    return 3;
                }
                if (shade > 1.0f) shade = 1.0f;
                rgb[0] = (unsigned char)(shade * 255.0f);
                rgb[1] = (unsigned char)(shade * (float)(128 + (id * 37) % 128));
                rgb[2] = (unsigned char)(shade * (float)(64 + (id * 91) % 192));
            }
            fwrite(rgb, 1, 3, f);
        }
    }
    if (fclose(f) != 0) { fprintf(stderr, "error: failed writing '%s'\n", path); return 2; }
    return 0;
}

static void usage(const char* prog) {
    fprintf(stderr,
        "usage: %s                  run self-tests\n"
        "       %s render OUT.ppm [--size N] [--strict|--no-strict]\n"
        "  N: image width/height in pixels, integer 1..1024 (default 128)\n",
        prog, prog);
}

static int parse_size(const char* s, int* out) {
    long v = 0;
    const char* p;
    if (!s || !*s) return -1;
    for (p = s; *p; p++) {
        if (*p < '0' || *p > '9') return -1;
        v = v * 10 + (*p - '0');
        if (v > 1024) return -1;
    }
    if (v < 1) return -1;
    *out = (int)v;
    return 0;
}

static int render_main(int argc, char** argv) {
    const char* out = NULL;
    int size = 128;
    int strict = 1;
    int i;
    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--size") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "error: --size requires a value\n");
                return 2;
            }
            if (parse_size(argv[++i], &size) != 0) {
                fprintf(stderr, "error: invalid --size '%s': expected integer 1..1024\n", argv[i]);
                return 2;
            }
        } else if (strcmp(argv[i], "--strict") == 0) {
            strict = 1;
        } else if (strcmp(argv[i], "--no-strict") == 0) {
            strict = 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
            return 2;
        } else if (!out) {
            out = argv[i];
        } else {
            fprintf(stderr, "error: unexpected argument '%s'\n", argv[i]);
            return 2;
        }
    }
    if (!out) {
        fprintf(stderr, "error: render requires an output path\n");
        return 2;
    }
    return render(out, size, strict);
}

int main(int argc, char** argv) {
    if (argc >= 2) {
        if (strcmp(argv[1], "render") == 0) return render_main(argc - 2, argv + 2);
        fprintf(stderr, "error: unknown command '%s'\n", argv[1]);
        usage(argv[0]);
        return 2;
    }
    Tree t;
    memset(&t, 0, sizeof(t));
    tree_add(&t, 1, 0, 0, 2, 1, 1);
    tree_add(&t, 5, 0, 0, 6, 1, 1);
    tree_add(&t, 3, 0, 0, 4, 1, 1);
    check("prim_count", t.nprim == 3);

    int root = build(&t, 0, t.nprim);
    check("build_ok", root == 0 && t.nnode >= 2);

    AABB whole;
    aabb_empty(&whole);
    for (int i = 0; i < t.nprim; i++) aabb_grow(&whole, &t.p[i].box);
    check("root_covers", fabsf(t.n[0].box.min[0] - whole.min[0]) < 1e-5f &&
                         fabsf(t.n[0].box.max[0] - whole.max[0]) < 1e-5f);

    float o[3] = { 0, 0.5f, 0.5f }, d[3] = { 1, 0, 0 };
    float tm = 0;
    check("hit_box", aabb_hit(&t.n[0].box, o, d, 100.0f, &tm) == 1 && tm > 0.0f);

    float d2[3] = { -1, 0, 0 };
    check("miss_box", aabb_hit(&t.n[0].box, o, d2, 100.0f, &tm) == 0);

    float o3[3] = { 0, 5, 5 };
    check("miss_parallel", aabb_hit(&t.n[0].box, o3, d, 100.0f, &tm) == 0);

    check("area_positive", aabb_area(&t.n[0].box) > 0.0f);

    Tree big;
    memset(&big, 0, sizeof(big));
    int rc = 0;
    for (int i = 0; i < MAX_PRIM; i++) rc = tree_add(&big, (float)i, 0, 0, (float)i + 0.5f, 1, 1);
    check("capacity", rc == MAX_PRIM - 1 && big.nprim == MAX_PRIM);
    check("build_full", build(&big, 0, big.nprim) == 0);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
