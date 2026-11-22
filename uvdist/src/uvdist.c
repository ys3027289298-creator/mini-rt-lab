/* uvdist: per-triangle UV distortion (analytic area ratio) and
   distortion-aware sample-rate allocation, with heatmap export.

   Metric
   ------
   For triangle i with world vertices p0,p1,p2 and UVs q0,q1,q2:
       A_geo_i = 0.5 * |(p1-p0) x (p2-p0)|
       A_uv_i  = 0.5 * |cross2(q1-q0, q2-q0)|
       density_i = A_uv_i / A_geo_i
   Fix a reference density `ref` (texel density at which the asset is meant
   to be sampled). The signed log2 distortion and the unsigned area-ratio
   distortion are
       s_i = log2(density_i / ref)
       D_i = 2^{|s_i|} = max(r_i, 1/r_i),  r_i = density_i / ref
   so D >= 1 always, and D == 1 exactly when the triangle matches the target.
   |s| = 1 means the UV area is off by 2x, |s| = 2 by 4x.
       s > 0 : UV area too large -> texels crammed in (over-dense)
       s < 0 : UV area too small -> texture stretched/magnified (sparse)

   Every per-triangle area is also recomputed by a second, independent
   formula (Heron in 3D, shoelace in UV); the export reports the maximum
   deviation, so the distortions are cross-checked against the analytic
   area ratio rather than merely asserted.

   Anisotropy
   ----------
   The 2x2 tangent map dP/du, dP/dv has singular values 0 < smin <= smax.
   Their ratio (condition number) detects shear-like stretch that preserves
   area - and therefore leaves D at 1. It is exported for inspection but does
   not drive allocation, which follows the area ratio.

   Sample allocation
   -----------------
   `budget` samples are split by weight
       w_i = A_geo_i * 2^{k |s_i|}                (k = 1 by default)
   A triangle with zero distortion (s_i snapped to exactly 0) has
   w_i = A_geo_i, so it can never receive more than its geometric-area share
   of the budget: zero-distortion regions are not oversampled. Integer counts
   use largest-remainder rounding, so the budget is spent exactly and the
   split is deterministic.

   Export (`--export DIR`) writes:
       stats.csv         per-triangle metrics, region, samples
       flagged.txt       indices with D > threshold
       summary.json      totals, cross-check residual, zero-distortion audit
       uv_heatmap.png    D in texture (UV) space
       world_heatmap.png D drawn on the geometry in world space
       checker.png       checker texture sampled through UV, in world space
       sample.png        samples per unit world area, in world space
       uv_cond.png       condition number in UV space (shear detector)
       report.html       self-contained report with inline images + tables

   Self-test (`--selftest`) validates against closed-form maps. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPS 1e-12
#define SNAP_S 1e-9            /* |s| below this counts as exactly zero */
#define MAX_SAMPLES_PER_TRI 1000000
#define DEFAULT_THRESHOLD 2.0  /* flag triangles whose area ratio is off > 2x */
#define RAMP_STEPS 11

typedef struct { double x, y, z; } V3;
typedef struct { double u, v; } V2;
typedef struct { V3 p[3]; V2 uv[3]; } Tri;

typedef struct {
    double a_geo, a_uv, density, r, s, D;
    double stretch_min, stretch_max, cond;
    int valid;
} TriStat;

/* ------------------------------------------------------------------ */
/* vector helpers                                                      */
/* ------------------------------------------------------------------ */
static V3 v3sub(V3 a, V3 b) { V3 r = { a.x-b.x, a.y-b.y, a.z-b.z }; return r; }
static V3 v3cross(V3 a, V3 b) {
    V3 r = { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
    return r;
}
static double v3dot(V3 a, V3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
static double v3len(V3 a) { return sqrt(v3dot(a, a)); }

/* Primary area formulas. */
static double area_geo(const Tri* t) {
    V3 e1 = v3sub(t->p[1], t->p[0]);
    V3 e2 = v3sub(t->p[2], t->p[0]);
    return 0.5 * v3len(v3cross(e1, e2));
}
static double area_uv(const Tri* t) {
    double x1 = t->uv[1].u - t->uv[0].u, y1 = t->uv[1].v - t->uv[0].v;
    double x2 = t->uv[2].u - t->uv[0].u, y2 = t->uv[2].v - t->uv[0].v;
    return 0.5 * fabs(x1 * y2 - x2 * y1);
}

/* Independent area formulas, used only for the consistency cross-check. */
static double area_geo_heron(const Tri* t) {
    double a = v3len(v3sub(t->p[1], t->p[0]));
    double b = v3len(v3sub(t->p[2], t->p[1]));
    double c = v3len(v3sub(t->p[0], t->p[2]));
    double sp = 0.5 * (a + b + c);
    double h = sp * (sp - a) * (sp - b) * (sp - c);
    return h > 0.0 ? sqrt(h) : 0.0;
}
static double area_uv_shoelace(const Tri* t) {
    double sum = t->uv[0].u * (t->uv[1].v - t->uv[2].v)
               + t->uv[1].u * (t->uv[2].v - t->uv[0].v)
               + t->uv[2].u * (t->uv[0].v - t->uv[1].v);
    return 0.5 * fabs(sum);
}

/* Singular values of the UV -> world tangent map. */
static void jacobian_singular(const Tri* t, double* smin, double* smax) {
    V3 e1 = v3sub(t->p[1], t->p[0]);
    V3 e2 = v3sub(t->p[2], t->p[0]);
    double u1 = t->uv[1].u - t->uv[0].u, v1 = t->uv[1].v - t->uv[0].v;
    double u2 = t->uv[2].u - t->uv[0].u, v2 = t->uv[2].v - t->uv[0].v;
    double det = u1 * v2 - u2 * v1;
    *smin = *smax = 1.0;
    if (fabs(det) < EPS) return;
    /* [e1 e2] = [P_u P_v] * [u1 u2; v1 v2]^T
       P_u = (e1*v2 - e2*v1)/det, P_v = (e2*u1 - e1*u2)/det */
    V3 pu = { (e1.x*v2 - e2.x*v1)/det, (e1.y*v2 - e2.y*v1)/det, (e1.z*v2 - e2.z*v1)/det };
    V3 pv = { (e2.x*u1 - e1.x*u2)/det, (e2.y*u1 - e1.y*u2)/det, (e2.z*u1 - e1.z*u2)/det };
    double a = v3dot(pu, pu), b = v3dot(pu, pv), c = v3dot(pv, pv);
    double tr = a + c;
    double disc = sqrt(fmax(0.0, tr*tr - 4.0*(a*c - b*b)));
    *smax = sqrt(fmax(0.0, 0.5 * (tr + disc)));
    *smin = sqrt(fmax(0.0, 0.5 * (tr - disc)));
}

/* ------------------------------------------------------------------ */
/* analysis                                                            */
/* ------------------------------------------------------------------ */
static void compute_areas(const Tri* tris, int n, TriStat* st) {
    for (int i = 0; i < n; i++) {
        st[i].a_geo = area_geo(&tris[i]);
        st[i].a_uv = area_uv(&tris[i]);
        st[i].valid = (st[i].a_geo > EPS && st[i].a_uv > EPS);
        st[i].density = st[i].valid ? st[i].a_uv / st[i].a_geo : 0.0;
        st[i].r = 0.0;
        st[i].s = 0.0;
        st[i].D = 1.0;
        jacobian_singular(&tris[i], &st[i].stretch_min, &st[i].stretch_max);
        st[i].cond = (st[i].stretch_min > EPS)
                   ? st[i].stretch_max / st[i].stretch_min : 1.0;
    }
}

static double area_weighted_density(const TriStat* st, int n) {
    double auv = 0.0, ageo = 0.0;
    for (int i = 0; i < n; i++) {
        if (!st[i].valid) continue;
        auv += st[i].a_uv;
        ageo += st[i].a_geo;
    }
    return (ageo > EPS) ? auv / ageo : 0.0;
}

/* Assign distortion against an explicit reference density. ref <= 0 falls
   back to the area-weighted mean (total UV area / total world area). */
static double assign_distortion(TriStat* st, int n, double ref) {
    if (ref <= 0.0) ref = area_weighted_density(st, n);
    for (int i = 0; i < n; i++) {
        if (!st[i].valid || ref <= EPS) continue;
        st[i].r = st[i].density / ref;
        double s = log2(st[i].r);
        if (fabs(s) < SNAP_S) s = 0.0;   /* zero distortion is exact */
        st[i].s = s;
        st[i].D = (s == 0.0) ? 1.0 : pow(2.0, fabs(s));
    }
    return ref;
}

/* ------------------------------------------------------------------ */
/* sample allocation                                                  */
/* ------------------------------------------------------------------ */
static int allocate_samples(const TriStat* st, int n, long budget,
                            double strength, int* out) {
    double* w = (double*)malloc((size_t)n * sizeof(double));
    double* frac = (double*)calloc((size_t)n, sizeof(double));
    if (!w || !frac) { free(w); free(frac); return -1; }
    double wsum = 0.0;
    for (int i = 0; i < n; i++) {
        double wi = 0.0;
        if (st[i].valid) {
            double mult = (st[i].s == 0.0) ? 1.0 : pow(2.0, strength * fabs(st[i].s));
            if (mult > (double)MAX_SAMPLES_PER_TRI) mult = (double)MAX_SAMPLES_PER_TRI;
            wi = st[i].a_geo * mult;
        }
        w[i] = wi;
        wsum += wi;
    }
    if (wsum <= EPS) { free(w); free(frac); return -1; }

    long assigned = 0;
    for (int i = 0; i < n; i++) {
        double exact = (double)budget * w[i] / wsum;
        out[i] = (int)floor(exact);
        frac[i] = exact - (double)out[i];
        assigned += out[i];
    }
    while (assigned < budget) {
        int best = -1;
        for (int i = 0; i < n; i++) {
            if (out[i] >= MAX_SAMPLES_PER_TRI) continue;
            if (best < 0 || frac[i] > frac[best] + 1e-15 ||
                (fabs(frac[i] - frac[best]) <= 1e-15 && i < best))
                best = i;
        }
        if (best < 0) break;
        out[best]++;
        frac[best] -= 1.0;
        assigned++;
    }
    free(w); free(frac);
    return 0;
}

/* ------------------------------------------------------------------ */
/* planted sheet mesh                                                  */
/* ------------------------------------------------------------------ */
/* A flat 12x8-cell sheet. UV layout has three column bands:
       columns 0-1  : UV widened  4x in u -> over-dense, s = +2, D = 4
       columns 2-7  : UV identity          -> s =  0, D = 1
       columns 8-11 : UV narrowed 0.25x    -> sparse,     s = -2, D = 4
   A global UV shear of 0.20*v rides on top. It is a det=1 linear map, so it
   leaves every triangle's UV area - hence every s - unchanged, while pushing
   the condition number above 1. That separates the two phenomena: area
   distortion drives sampling; anisotropy is only reported. */
#define SHEET_CX 12
#define SHEET_CY 8
#define SHEET_NVX (SHEET_CX + 1)
#define SHEET_NVY (SHEET_CY + 1)
#define SHEET_TRI (2 * SHEET_CX * SHEET_CY)
#define SHEET_REF_TRI (2 * 3)      /* cell (3,0): an identity column */

static void sheet_uv(int i, int j, V2* out) {
    double u = 0.0;
    for (int c = 0; c < i; c++) {
        double du = 1.0 / SHEET_CX;
        if (c < 2) du *= 4.0;
        else if (c >= 8) du *= 0.25;
        u += du;
    }
    double v = (double)j / SHEET_CY;
    out->u = u + 0.20 * v;
    out->v = v;
}

static void build_sheet(Tri* tris) {
    V3 pos[SHEET_NVX][SHEET_NVY];
    V2 uvs[SHEET_NVX][SHEET_NVY];
    for (int i = 0; i < SHEET_NVX; i++)
        for (int j = 0; j < SHEET_NVY; j++) {
            pos[i][j].x = 0.2 * i;
            pos[i][j].y = 0.2 * j;
            pos[i][j].z = 0.0;
            sheet_uv(i, j, &uvs[i][j]);
        }
    int ti = 0;
    for (int j = 0; j < SHEET_CY; j++)
        for (int i = 0; i < SHEET_CX; i++) {
            Tri* a = &tris[ti++];
            a->p[0] = pos[i][j];     a->p[1] = pos[i+1][j];   a->p[2] = pos[i][j+1];
            a->uv[0] = uvs[i][j];    a->uv[1] = uvs[i+1][j];  a->uv[2] = uvs[i][j+1];
            Tri* b = &tris[ti++];
            b->p[0] = pos[i+1][j];   b->p[1] = pos[i+1][j+1]; b->p[2] = pos[i][j+1];
            b->uv[0] = uvs[i+1][j];  b->uv[1] = uvs[i+1][j+1]; b->uv[2] = uvs[i][j+1];
        }
}

static int sheet_col(int tri) { return (tri / 2) % SHEET_CX; }
static int sheet_row(int tri) { return (tri / 2) / SHEET_CX; }
static const char* sheet_region(int tri) {
    int col = sheet_col(tri);
    if (col < 2) return "overdense";
    if (col < 8) return "identity";
    return "sparse";
}

/* ------------------------------------------------------------------ */
/* image output (PPM, P6)                                             */
/* ------------------------------------------------------------------ */
typedef struct { int w, h; unsigned char* px; } Image;

static Image* image_new(int w, int h, const unsigned char bg[3]) {
    Image* im = (Image*)malloc(sizeof(Image));
    if (!im) return NULL;
    im->w = w; im->h = h;
    im->px = (unsigned char*)malloc((size_t)w * h * 3);
    if (!im->px) { free(im); return NULL; }
    for (int i = 0; i < w * h; i++) {
        im->px[i*3+0] = bg[0];
        im->px[i*3+1] = bg[1];
        im->px[i*3+2] = bg[2];
    }
    return im;
}
static void image_free(Image* im) { if (!im) return; free(im->px); free(im); }
static void image_put(Image* im, int x, int y, const unsigned char rgb[3]) {
    if (x < 0 || y < 0 || x >= im->w || y >= im->h) return;
    int o = (y * im->w + x) * 3;
    im->px[o+0] = rgb[0];
    im->px[o+1] = rgb[1];
    im->px[o+2] = rgb[2];
}
static int image_save(const Image* im, const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return -1; }
    fprintf(f, "P6\n%d %d\n255\n", im->w, im->h);
    fwrite(im->px, 1, (size_t)im->w * im->h * 3, f);
    fclose(f);
    return 0;
}

/* Minimal PNG writer: zlib stream built from stored (uncompressed) deflate
   blocks, so the module needs no zlib dependency. Produces valid PNGs that
   any viewer or browser can open. */
static unsigned long crc_table[256];
static int crc_table_ready = 0;
static void crc_build(void) {
    for (unsigned long n = 0; n < 256; n++) {
        unsigned long c = n;
        for (int k = 0; k < 8; k++)
            c = (c & 1) ? 0xedb88320UL ^ (c >> 1) : c >> 1;
        crc_table[n] = c;
    }
    crc_table_ready = 1;
}
static unsigned long crc32_buf(const unsigned char* p, size_t n) {
    if (!crc_table_ready) crc_build();
    unsigned long c = 0xffffffffUL;
    for (size_t i = 0; i < n; i++) c = crc_table[(c ^ p[i]) & 0xff] ^ (c >> 8);
    return c ^ 0xffffffffUL;
}
static unsigned long adler32_buf(const unsigned char* p, size_t n) {
    unsigned long a = 1, b = 0;
    for (size_t i = 0; i < n; i++) { a = (a + p[i]) % 65521; b = (b + a) % 65521; }
    return (b << 16) | a;
}
static void be32(unsigned char* o, unsigned long v) {
    o[0] = (unsigned char)(v >> 24); o[1] = (unsigned char)(v >> 16);
    o[2] = (unsigned char)(v >> 8);  o[3] = (unsigned char)v;
}
static int png_chunk(FILE* f, const char type[4], const unsigned char* data,
                     size_t len) {
    unsigned char hdr[8];
    be32(hdr, (unsigned long)len);
    memcpy(hdr + 4, type, 4);
    fwrite(hdr, 1, 8, f);
    if (len) fwrite(data, 1, len, f);
    unsigned char* crcbuf = (unsigned char*)malloc(len + 4);
    if (!crcbuf) return -1;
    memcpy(crcbuf, type, 4);
    if (len) memcpy(crcbuf + 4, data, len);
    unsigned char crc[4];
    be32(crc, crc32_buf(crcbuf, len + 4));
    free(crcbuf);
    fwrite(crc, 1, 4, f);
    return 0;
}
/* --- fixed-Huffman deflate (RFC 1951) ---------------------------------- */
typedef struct {
    unsigned char* buf;
    size_t len, cap;
    unsigned bitbuf;
    int bitcnt;
    int oom;
} BitW;

static void bw_byte(BitW* w, unsigned char b) {
    if (w->len == w->cap) {
        size_t nc = w->cap ? w->cap * 2 : 4096;
        unsigned char* nb = (unsigned char*)realloc(w->buf, nc);
        if (!nb) { w->oom = 1; return; }
        w->buf = nb; w->cap = nc;
    }
    w->buf[w->len++] = b;
}
static void bw_put(BitW* w, unsigned v, int n) {
    w->bitbuf |= (v & ((1u << n) - 1u)) << w->bitcnt;
    w->bitcnt += n;
    while (w->bitcnt >= 8) {
        bw_byte(w, (unsigned char)(w->bitbuf & 0xff));
        w->bitbuf >>= 8;
        w->bitcnt -= 8;
    }
}
static unsigned revbits(unsigned v, int n) {
    unsigned r = 0;
    for (int i = 0; i < n; i++) { r = (r << 1) | (v & 1u); v >>= 1; }
    return r;
}
static void bw_huff(BitW* w, unsigned code, int n) { bw_put(w, revbits(code, n), n); }

static void emit_lit(BitW* w, int sym) {
    if (sym <= 143)       bw_huff(w, 0x30u + (unsigned)sym, 8);
    else if (sym <= 255)  bw_huff(w, 0x190u + (unsigned)(sym - 144), 9);
    else if (sym <= 279)  bw_huff(w, (unsigned)(sym - 256), 7);
    else                  bw_huff(w, 0xC0u + (unsigned)(sym - 280), 8);
}

static const int LEN_BASE[29] = {
    3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,
    67,83,99,115,131,163,195,227,258 };
static const int LEN_EXTRA[29] = {
    0,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,0 };
static const int DIST_BASE[30] = {
    1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,
    1025,1537,2049,3073,4097,6145,8193,12289,16385,24577 };
static const int DIST_EXTRA[30] = {
    0,0,0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13 };

static void emit_match(BitW* w, int len, int dist) {
    int lc = 0;
    while (lc < 28 && LEN_BASE[lc + 1] <= len) lc++;
    emit_lit(w, 257 + lc);
    if (LEN_EXTRA[lc]) bw_put(w, (unsigned)(len - LEN_BASE[lc]), LEN_EXTRA[lc]);
    int dc = 0;
    while (dc < 29 && DIST_BASE[dc + 1] <= dist) dc++;
    bw_huff(w, (unsigned)dc, 5);
    if (DIST_EXTRA[dc]) bw_put(w, (unsigned)(dist - DIST_BASE[dc]), DIST_EXTRA[dc]);
}

/* Compress into a zlib stream (header + fixed-Huffman deflate + adler32).
   A single-slot hash table over 3-byte sequences with greedy matching is
   enough: row-repeated image data compresses many-fold. */
#define DH_HASH 32768
static unsigned char* zlib_deflate_fixed(const unsigned char* src, size_t n,
                                         size_t* out_len) {
    BitW w;
    memset(&w, 0, sizeof(w));
    bw_put(&w, 0x78, 8);          /* CMF */
    bw_put(&w, 0x01, 8);          /* FLG: fixed Huffman, check bit */
    bw_put(&w, 1, 1);             /* BFINAL = 1 */
    bw_put(&w, 1, 2);             /* BTYPE = 01 (fixed Huffman) */

    int* head = (int*)malloc(sizeof(int) * DH_HASH);
    if (!head) { free(w.buf); return NULL; }
    for (int i = 0; i < DH_HASH; i++) head[i] = -1;

    size_t i = 0;
    while (i < n) {
        int best_len = 0, best_dist = 0;
        if (i + 2 < n) {
            unsigned h = ((unsigned)src[i] * 0x1e35a7bdu +
                          (unsigned)src[i+1] * 0x9e3779b1u +
                          (unsigned)src[i+2]) & (DH_HASH - 1);
            int cand = head[h];
            head[h] = (int)i;
            if (cand >= 0 && (size_t)cand < i && i - (size_t)cand <= 32768) {
                size_t maxlen = n - i;
                if (maxlen > 258) maxlen = 258;
                size_t l = 0;
                while (l < maxlen && src[cand + l] == src[i + l]) l++;
                if (l >= 3) { best_len = (int)l; best_dist = (int)(i - cand); }
            }
            /* seed hashes for the bytes the match covers */
            size_t upto = i + (best_len > 0 ? (size_t)best_len : 1);
            for (size_t k = i + 1; k < upto && k + 2 < n; k++) {
                unsigned hh = ((unsigned)src[k] * 0x1e35a7bdu +
                               (unsigned)src[k+1] * 0x9e3779b1u +
                               (unsigned)src[k+2]) & (DH_HASH - 1);
                head[hh] = (int)k;
            }
        }
        if (best_len >= 3) {
            emit_match(&w, best_len, best_dist);
            i += (size_t)best_len;
        } else {
            emit_lit(&w, src[i]);
            i++;
        }
    }
    emit_lit(&w, 256);            /* end of block */
    if (w.bitcnt > 0) bw_byte(&w, (unsigned char)(w.bitbuf & 0xff));
    free(head);
    if (w.oom) { free(w.buf); return NULL; }

    unsigned long ad = adler32_buf(src, n);
    bw_byte(&w, (unsigned char)(ad >> 24));
    bw_byte(&w, (unsigned char)(ad >> 16));
    bw_byte(&w, (unsigned char)(ad >> 8));
    bw_byte(&w, (unsigned char)ad);
    if (w.oom) { free(w.buf); return NULL; }
    *out_len = w.len;
    return w.buf;
}

static int image_save_png(const Image* im, const char* path) {
    size_t stride = (size_t)im->w * 3 + 1;
    size_t rawlen = stride * (size_t)im->h;
    unsigned char* raw = (unsigned char*)malloc(rawlen);
    if (!raw) return -1;
    for (int y = 0; y < im->h; y++) {
        raw[(size_t)y * stride] = 0;   /* filter: none */
        memcpy(raw + (size_t)y * stride + 1, im->px + (size_t)y * im->w * 3,
               (size_t)im->w * 3);
    }
    size_t zlen = 0;
    unsigned char* z = zlib_deflate_fixed(raw, rawlen, &zlen);
    free(raw);
    if (!z) return -1;

    FILE* f = fopen(path, "wb");
    if (!f) { free(z); perror(path); return -1; }
    static const unsigned char sig[8] = { 0x89, 'P', 'N', 'G', 0x0d, 0x0a, 0x1a, 0x0a };
    fwrite(sig, 1, 8, f);
    unsigned char ihdr[13];
    be32(ihdr, (unsigned long)im->w);
    be32(ihdr + 4, (unsigned long)im->h);
    ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
    png_chunk(f, "IHDR", ihdr, 13);
    png_chunk(f, "IDAT", z, zlen);
    png_chunk(f, "IEND", NULL, 0);
    fclose(f);
    free(z);
    return 0;
}
/* Encode an image to a PNG byte buffer in memory (for the HTML report). */
static unsigned char* image_png_bytes(const Image* im, size_t* out_len) {
    size_t stride = (size_t)im->w * 3 + 1;
    size_t rawlen = stride * (size_t)im->h;
    unsigned char* raw = (unsigned char*)malloc(rawlen);
    if (!raw) return NULL;
    for (int y = 0; y < im->h; y++) {
        raw[(size_t)y * stride] = 0;
        memcpy(raw + (size_t)y * stride + 1, im->px + (size_t)y * im->w * 3,
               (size_t)im->w * 3);
    }
    size_t zlen = 0;
    unsigned char* z = zlib_deflate_fixed(raw, rawlen, &zlen);
    free(raw);
    if (!z) return NULL;

    size_t cap = zlen + 128;
    unsigned char* out = (unsigned char*)malloc(cap);
    if (!out) { free(z); return NULL; }
    size_t o = 0;
    static const unsigned char sig[8] = { 0x89, 'P', 'N', 'G', 0x0d, 0x0a, 0x1a, 0x0a };
    memcpy(out + o, sig, 8); o += 8;
    unsigned char ihdr[13];
    be32(ihdr, (unsigned long)im->w);
    be32(ihdr + 4, (unsigned long)im->h);
    ihdr[8] = 8; ihdr[9] = 2; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
    /* chunks into memory */
    struct { const char* type; const unsigned char* data; size_t len; } chunks[3] = {
        { "IHDR", ihdr, 13 }, { "IDAT", z, zlen }, { "IEND", NULL, 0 }
    };
    for (int c = 0; c < 3; c++) {
        unsigned char hdr[8];
        be32(hdr, (unsigned long)chunks[c].len);
        memcpy(hdr + 4, chunks[c].type, 4);
        memcpy(out + o, hdr, 8); o += 8;
        if (chunks[c].len) { memcpy(out + o, chunks[c].data, chunks[c].len); o += chunks[c].len; }
        unsigned char* crcbuf = (unsigned char*)malloc(chunks[c].len + 4);
        memcpy(crcbuf, chunks[c].type, 4);
        if (chunks[c].len) memcpy(crcbuf + 4, chunks[c].data, chunks[c].len);
        unsigned char crc[4];
        be32(crc, crc32_buf(crcbuf, chunks[c].len + 4));
        free(crcbuf);
        memcpy(out + o, crc, 4); o += 4;
    }
    free(z);
    *out_len = o;
    return out;
}
/* Pick PNG for .png paths, PPM otherwise. */
static int image_write(const Image* im, const char* path) {
    size_t n = strlen(path);
    if (n >= 4 && !strcmp(path + n - 4, ".png")) return image_save_png(im, path);
    return image_save(im, path);
}

/* Sequential blue ramp (validated palette, steps 100..700). */
static const unsigned char BLUE[RAMP_STEPS][3] = {
    { 0xcd, 0xe2, 0xfb }, { 0xb7, 0xd3, 0xf6 }, { 0x9e, 0xc5, 0xf4 },
    { 0x86, 0xb6, 0xef }, { 0x6d, 0xa7, 0xec }, { 0x55, 0x98, 0xe7 },
    { 0x39, 0x87, 0xe5 }, { 0x2a, 0x78, 0xd6 }, { 0x25, 0x6a, 0xbf },
    { 0x1c, 0x5c, 0xab }, { 0x0d, 0x36, 0x6b },
};
static void ramp(double x, double vmin, double vmax, unsigned char out[3]) {
    double f = (vmax > vmin) ? (x - vmin) / (vmax - vmin) : 0.0;
    if (!(f > 0.0)) f = 0.0;
    if (f > 1.0) f = 1.0;
    /* Continuous interpolation between the documented steps. Equal values
       therefore always render to the same colour; snapping to the nearest
       discrete step made float-equal values straddle a step boundary and
       produced a two-tone ripple across uniform regions. */
    double pos = f * (RAMP_STEPS - 1);
    int i0 = (int)floor(pos);
    if (i0 >= RAMP_STEPS - 1) i0 = RAMP_STEPS - 2;
    if (i0 < 0) i0 = 0;
    double t = pos - (double)i0;
    for (int c = 0; c < 3; c++) {
        double v = (double)BLUE[i0][c] * (1.0 - t) + (double)BLUE[i0 + 1][c] * t;
        out[c] = (unsigned char)floor(v + 0.5);
    }
}

/* Choose pixel dimensions that match the rect's aspect so the raster is
   not anisotropically scaled. */
static void fit_dims(const double rect[4], int maxw, int maxh, int* w, int* h) {
    double rw = rect[2] - rect[0], rh = rect[3] - rect[1];
    if (rw < 1e-12) rw = 1.0;
    if (rh < 1e-12) rh = 1.0;
    double aspect = rw / rh;
    if (aspect >= (double)maxw / (double)maxh) {
        *w = maxw;
        *h = (int)(maxw / aspect + 0.5);
    } else {
        *h = maxh;
        *w = (int)(maxh * aspect + 0.5);
    }
    if (*w < 64) *w = 64;
    if (*h < 64) *h = 64;
}

typedef void (*PixelFn)(void* payload, double b0, double b1, double b2,
                        int x, int y);

/* Rasterise one triangle over an orthographic rect. `gap` > 0 insets every
   edge by a fraction of a pixel so adjacent fills stay separable without a
   drawn border. */
static void raster_tri(Image* im, const double xy[3][2],
                       const double rect[4], PixelFn cb, void* payload,
                       int gap) {
    double rw = rect[2] - rect[0], rh = rect[3] - rect[1];
    double xd[3], yd[3];
    for (int k = 0; k < 3; k++) {
        xd[k] = (xy[k][0] - rect[0]) / rw * (im->w - 1);
        yd[k] = (1.0 - (xy[k][1] - rect[1]) / rh) * (im->h - 1);
    }
    int bx0 = (int)floor(fmin(xd[0], fmin(xd[1], xd[2]))) - 1;
    int bx1 = (int)ceil (fmax(xd[0], fmax(xd[1], xd[2]))) + 1;
    int by0 = (int)floor(fmin(yd[0], fmin(yd[1], yd[2]))) - 1;
    int by1 = (int)ceil (fmax(yd[0], fmax(yd[1], yd[2]))) + 1;
    if (bx0 < 0) bx0 = 0;
    if (by0 < 0) by0 = 0;
    if (bx1 >= im->w) bx1 = im->w - 1;
    if (by1 >= im->h) by1 = im->h - 1;

    double den = (yd[1] - yd[2]) * (xd[0] - xd[2]) + (xd[2] - xd[1]) * (yd[0] - yd[2]);
    if (fabs(den) < 1e-12) return;
    double margin = (gap > 0) ? (0.004 + 0.0015 * gap) : 0.0;
    for (int py = by0; py <= by1; py++) {
        for (int px = bx0; px <= bx1; px++) {
            double x = px + 0.5, y = py + 0.5;
            double b0 = ((yd[1] - yd[2]) * (x - xd[2]) + (xd[2] - xd[1]) * (y - yd[2])) / den;
            double b1 = ((yd[2] - yd[0]) * (x - xd[2]) + (xd[0] - xd[2]) * (y - yd[2])) / den;
            double b2 = 1.0 - b0 - b1;
            if (b0 < margin || b1 < margin || b2 < margin) continue;
            cb(payload, b0, b1, b2, px, py);
        }
    }
}

typedef struct {
    Image* im;
    unsigned char rgb[3];
} ConstPaint;
static void const_paint(void* payload, double b0, double b1, double b2,
                        int x, int y) {
    (void)b0; (void)b1; (void)b2;
    ConstPaint* p = (ConstPaint*)payload;
    image_put(p->im, x, y, p->rgb);
}

typedef struct {
    Image* im;
    const double (*uv)[2];
    int checks;
    unsigned char a[3], b[3];
} CheckerPaint;
static void checker_paint(void* payload, double b0, double b1, double b2,
                          int x, int y) {
    CheckerPaint* p = (CheckerPaint*)payload;
    double u = b0 * p->uv[0][0] + b1 * p->uv[1][0] + b2 * p->uv[2][0];
    double v = b0 * p->uv[0][1] + b1 * p->uv[1][1] + b2 * p->uv[2][1];
    long ci = (long)floor(u * p->checks);
    long cj = (long)floor(v * p->checks);
    int odd = (int)(((ci + cj) % 2 + 2) % 2);
    image_put(p->im, x, y, odd ? p->b : p->a);
}

/* Paint every triangle with one constant colour, given a per-triangle value
   and a colour map. `space` selects UV (0) or world (1) coordinates. */
static void draw_const_map(Image* im, const Tri* tris, const double* vals,
                           int n, int space, double vmin, double vmax,
                           const double rect[4], int gap) {
    for (int i = 0; i < n; i++) {
        double xy[3][2];
        for (int k = 0; k < 3; k++) {
            if (space == 0) { xy[k][0] = tris[i].uv[k].u; xy[k][1] = tris[i].uv[k].v; }
            else            { xy[k][0] = tris[i].p[k].x;  xy[k][1] = tris[i].p[k].y;  }
        }
        ConstPaint p;
        p.im = im;
        ramp(vals[i], vmin, vmax, p.rgb);
        raster_tri(im, xy, rect, const_paint, &p, gap);
    }
}

/* ------------------------------------------------------------------ */
/* export                                                              */
/* ------------------------------------------------------------------ */
static void bounds_uv(const Tri* tris, int n, double rect[4]);
static void bounds_world(const Tri* tris, int n, double rect[4]);

static void write_report(const char* dir, const TriStat* st,
                         const int* alloc, int n, double thresh, long budget,
                         double dmax, double condmax, double smax,
                         Image* img_uv, Image* img_world, Image* img_cond,
                         Image* img_sample, Image* img_checker,
                         long total, double total_area, int flagged,
                         double flagged_area, long flagged_samples,
                         int zero_count, double zero_area, long zero_samples,
                         double max_zero_ratio, double area_diff);

static int export_all(const char* dir, const Tri* tris, const TriStat* st,
                      const int* alloc, long budget, int n, double thresh) {
    char path[1024];

    /* ---- stats.csv ---- */
    snprintf(path, sizeof(path), "%s/stats.csv", dir);
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); return -1; }
    fprintf(f, "tri,col,row,region,a_geo,a_uv,density,ratio_to_ref,"
               "log2_stretch,distortion,stretch_min,stretch_max,cond,samples\n");
    for (int i = 0; i < n; i++)
        fprintf(f, "%d,%d,%d,%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d\n",
                i, sheet_col(i), sheet_row(i), sheet_region(i),
                st[i].a_geo, st[i].a_uv, st[i].density, st[i].r,
                st[i].s, st[i].D, st[i].stretch_min, st[i].stretch_max,
                st[i].cond, alloc ? alloc[i] : 0);
    fclose(f);

    /* ---- flagged.txt ---- */
    snprintf(path, sizeof(path), "%s/flagged.txt", dir);
    f = fopen(path, "w");
    if (!f) { perror(path); return -1; }
    fprintf(f, "# triangles with distortion D > %.6g (area ratio off by more than that)\n", thresh);
    fprintf(f, "# tri  region  distortion  log2_stretch  samples\n");
    for (int i = 0; i < n; i++)
        if (st[i].valid && st[i].D > thresh)
            fprintf(f, "%d  %s  %.6g  %.6g  %d\n",
                    i, sheet_region(i), st[i].D, st[i].s, alloc ? alloc[i] : 0);
    fclose(f);

    /* ---- per-triangle value arrays ---- */
    double* dvals = (double*)malloc(sizeof(double) * n);
    double* cvals = (double*)malloc(sizeof(double) * n);
    double* svals = (double*)malloc(sizeof(double) * n);
    if (!dvals || !cvals || !svals) { free(dvals); free(cvals); free(svals); return -1; }
    double dmax = 1.0, condmax = 1.0, smax = 0.0;
    for (int i = 0; i < n; i++) {
        dvals[i] = st[i].valid ? st[i].D : 1.0;
        cvals[i] = st[i].valid ? st[i].cond : 1.0;
        svals[i] = (alloc && st[i].valid && st[i].a_geo > EPS)
                 ? (double)alloc[i] / st[i].a_geo : 0.0;
        if (dvals[i] > dmax) dmax = dvals[i];
        if (cvals[i] > condmax) condmax = cvals[i];
        if (svals[i] > smax) smax = svals[i];
    }
    if (smax <= 0.0) smax = 1.0;

    const unsigned char surface[3] = { 0xfc, 0xfc, 0xfb };
    double urect[4], wrect[4];
    bounds_uv(tris, n, urect);
    bounds_world(tris, n, wrect);

    /* Render every view once; save as PNG and keep the image for the HTML
       report. */
    Image* img_uv = NULL;
    Image* img_world = NULL;
    Image* img_cond = NULL;
    Image* img_sample = NULL;
    Image* img_checker = NULL;
    {
        int w, h;
        fit_dims(urect, 900, 900, &w, &h);
        img_uv = image_new(w, h, surface);
        if (img_uv) {
            draw_const_map(img_uv, tris, dvals, n, 0, 1.0, dmax, urect, 1);
            snprintf(path, sizeof(path), "%s/uv_heatmap.png", dir);
            image_write(img_uv, path);
        }
    }
    {
        int w, h;
        fit_dims(wrect, 1000, 1000, &w, &h);
        img_world = image_new(w, h, surface);
        if (img_world) {
            draw_const_map(img_world, tris, dvals, n, 1, 1.0, dmax, wrect, 1);
            snprintf(path, sizeof(path), "%s/world_heatmap.png", dir);
            image_write(img_world, path);
        }
    }
    {
        int w, h;
        fit_dims(urect, 900, 900, &w, &h);
        img_cond = image_new(w, h, surface);
        if (img_cond) {
            draw_const_map(img_cond, tris, cvals, n, 0, 1.0, condmax, urect, 1);
            snprintf(path, sizeof(path), "%s/uv_cond.png", dir);
            image_write(img_cond, path);
        }
    }
    {
        int w, h;
        fit_dims(wrect, 1000, 1000, &w, &h);
        img_sample = image_new(w, h, surface);
        if (img_sample) {
            draw_const_map(img_sample, tris, svals, n, 1, 0.0, smax, wrect, 1);
            snprintf(path, sizeof(path), "%s/sample.png", dir);
            image_write(img_sample, path);
        }
    }
    {
        int w, h;
        fit_dims(wrect, 1000, 1000, &w, &h);
        const unsigned char chk_a[3] = { 0xf0, 0xef, 0xec };
        const unsigned char chk_b[3] = { 0x52, 0x51, 0x4e };
        img_checker = image_new(w, h, chk_a);
        if (img_checker) {
            for (int i = 0; i < n; i++) {
                double xy[3][2], uvp[3][2];
                for (int k = 0; k < 3; k++) {
                    xy[k][0] = tris[i].p[k].x;
                    xy[k][1] = tris[i].p[k].y;
                    uvp[k][0] = tris[i].uv[k].u;
                    uvp[k][1] = tris[i].uv[k].v;
                }
                CheckerPaint p;
                p.im = img_checker;
                p.uv = (const double (*)[2])uvp;
                p.checks = 16;
                memcpy(p.a, chk_a, 3);
                memcpy(p.b, chk_b, 3);
                raster_tri(img_checker, xy, wrect, checker_paint, &p, 0);
            }
            snprintf(path, sizeof(path), "%s/checker.png", dir);
            image_write(img_checker, path);
        }
    }
    free(dvals); free(cvals); free(svals);

    /* ---- audits ---- */
    int valid = 0, flagged = 0, zero_count = 0;
    long total = 0, flagged_samples = 0, zero_samples = 0;
    double zero_area = 0.0, flagged_area = 0.0, total_area = 0.0;
    double max_geo_diff = 0.0, max_uv_diff = 0.0, max_zero_ratio = 0.0;
    for (int i = 0; i < n; i++) {
        if (!st[i].valid) continue;
        valid++;
        total += alloc ? alloc[i] : 0;
        total_area += st[i].a_geo;
        double dg = fabs(area_geo(&tris[i]) - area_geo_heron(&tris[i]));
        double du = fabs(area_uv(&tris[i]) - area_uv_shoelace(&tris[i]));
        if (dg > max_geo_diff) max_geo_diff = dg;
        if (du > max_uv_diff) max_uv_diff = du;
        if (st[i].D > thresh) {
            flagged++;
            flagged_area += st[i].a_geo;
            flagged_samples += alloc ? alloc[i] : 0;
        }
        if (st[i].s == 0.0) {
            zero_count++;
            zero_area += st[i].a_geo;
            zero_samples += alloc ? alloc[i] : 0;
        }
    }
    if (alloc && total > 0 && total_area > 0)
        for (int i = 0; i < n; i++) {
            if (!st[i].valid || st[i].s != 0.0) continue;
            double area_share = st[i].a_geo / total_area;
            double sample_share = (double)alloc[i] / (double)total;
            if (area_share <= EPS) continue;
            double ratio = sample_share / area_share;
            if (ratio > max_zero_ratio) max_zero_ratio = ratio;
        }

    snprintf(path, sizeof(path), "%s/summary.json", dir);
    f = fopen(path, "w");
    if (!f) { perror(path); return -1; }
    fprintf(f, "{\n");
    fprintf(f, "  \"mesh\": \"sheet-%dx%d-cells\",\n", SHEET_CX, SHEET_CY);
    fprintf(f, "  \"triangles\": %d,\n", n);
    fprintf(f, "  \"valid_triangles\": %d,\n", valid);
    fprintf(f, "  \"reference\": \"density of identity-cell triangle %d\",\n", SHEET_REF_TRI);
    fprintf(f, "  \"reference_density\": %.12g,\n", st[SHEET_REF_TRI].density);
    fprintf(f, "  \"threshold_distortion\": %.9g,\n", thresh);
    fprintf(f, "  \"flagged_count\": %d,\n", flagged);
    fprintf(f, "  \"flagged_area_fraction\": %.9g,\n",
            total_area > EPS ? flagged_area / total_area : 0.0);
    fprintf(f, "  \"flagged_sample_fraction\": %.9g,\n",
            total > 0 ? (double)flagged_samples / (double)total : 0.0);
    fprintf(f, "  \"sample_budget\": %ld,\n", budget);
    fprintf(f, "  \"samples_total\": %ld,\n", total);
    fprintf(f, "  \"zero_distortion_triangles\": %d,\n", zero_count);
    fprintf(f, "  \"zero_distortion_area_fraction\": %.9g,\n",
            total_area > EPS ? zero_area / total_area : 0.0);
    fprintf(f, "  \"zero_distortion_sample_fraction\": %.9g,\n",
            total > 0 ? (double)zero_samples / (double)total : 0.0);
    fprintf(f, "  \"zero_distortion_max_sample_over_area\": %.9g,\n", max_zero_ratio);
    fprintf(f, "  \"zero_distortion_oversampled\": %d,\n", max_zero_ratio > 1.0 + 1e-6);
    fprintf(f, "  \"analytic_area_max_abs_diff\": %.3g,\n",
            max_geo_diff > max_uv_diff ? max_geo_diff : max_uv_diff);
    fprintf(f, "  \"max_distortion\": %.9g,\n", dmax);
    fprintf(f, "  \"max_cond\": %.9g,\n", condmax);
    fprintf(f, "  \"flagged_indices\": [");
    {
        int first = 1;
        for (int i = 0; i < n; i++) {
            if (!st[i].valid || st[i].D <= thresh) continue;
            fprintf(f, "%s%d", first ? "" : ",", i);
            first = 0;
        }
    }
    fprintf(f, "]\n}\n");
    fclose(f);

    /* ---- self-contained HTML report ---- */
    write_report(dir, st, alloc, n, thresh, budget, dmax, condmax, smax,
                 img_uv, img_world, img_cond, img_sample, img_checker,
                 total, total_area, flagged, flagged_area, flagged_samples,
                 zero_count, zero_area, zero_samples, max_zero_ratio,
                 max_geo_diff > max_uv_diff ? max_geo_diff : max_uv_diff);

    image_free(img_uv); image_free(img_world); image_free(img_cond);
    image_free(img_sample); image_free(img_checker);
    return 0;
}

static void bounds_uv(const Tri* tris, int n, double rect[4]) {
    double umin = 1e30, umax = -1e30, vmin = 1e30, vmax = -1e30;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++) {
            double u = tris[i].uv[k].u, v = tris[i].uv[k].v;
            if (u < umin) umin = u;
            if (u > umax) umax = u;
            if (v < vmin) vmin = v;
            if (v > vmax) vmax = v;
        }
    double w = umax - umin, h = vmax - vmin;
    double pad = 0.04 * (w > h ? w : h);
    if (pad <= 0) pad = 0.05;
    rect[0] = umin - pad; rect[1] = vmin - pad;
    rect[2] = umax + pad; rect[3] = vmax + pad;
}
static void bounds_world(const Tri* tris, int n, double rect[4]) {
    double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
    for (int i = 0; i < n; i++)
        for (int k = 0; k < 3; k++) {
            if (tris[i].p[k].x < xmin) xmin = tris[i].p[k].x;
            if (tris[i].p[k].x > xmax) xmax = tris[i].p[k].x;
            if (tris[i].p[k].y < ymin) ymin = tris[i].p[k].y;
            if (tris[i].p[k].y > ymax) ymax = tris[i].p[k].y;
        }
    double w = xmax - xmin, h = ymax - ymin;
    double pad = 0.04 * (w > h ? w : h);
    if (pad <= 0) pad = 0.05;
    rect[0] = xmin - pad; rect[1] = ymin - pad;
    rect[2] = xmax + pad; rect[3] = ymax + pad;
}

/* ------------------------------------------------------------------ */
/* HTML report                                                         */
/* ------------------------------------------------------------------ */
/* Data URI for an image, encoded as a PNG in memory. */
static void html_img_uri(FILE* f, const Image* im) {
    if (!im) return;
    size_t n = 0;
    unsigned char* png = image_png_bytes(im, &n);
    if (!png) return;
    static const char* b64 =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    fputs("data:image/png;base64,", f);
    for (size_t i = 0; i < n; i += 3) {
        unsigned v = (unsigned)png[i] << 16;
        if (i + 1 < n) v |= (unsigned)png[i+1] << 8;
        if (i + 2 < n) v |= (unsigned)png[i+2];
        fputc(b64[(v >> 18) & 63], f);
        fputc(b64[(v >> 12) & 63], f);
        fputc(i + 1 < n ? b64[(v >> 6) & 63] : '=', f);
        fputc(i + 2 < n ? b64[v & 63] : '=', f);
    }
    free(png);
}

static void html_scale_legend(FILE* f, double vmin, double vmax, const char* label) {
    fprintf(f, "<div class=\"scale\"><span class=\"scale-label\">%s</span>"
               "<span class=\"scale-min\">%.3g</span>"
               "<span class=\"scale-bar\">", label, vmin);
    for (int i = 0; i < RAMP_STEPS; i++)
        fprintf(f, "<i style=\"background:#%02x%02x%02x\"></i>",
                BLUE[i][0], BLUE[i][1], BLUE[i][2]);
    fprintf(f, "</span><span class=\"scale-max\">%.3g</span></div>\n", vmax);
}

/* Self-contained, offline HTML: all images inline, plus a table view of
   every triangle and the flagged list. Follows the same palette as the
   PNGs. The table is the WCAG-clean equivalent of the heatmaps, so no
   value is reachable only through colour. */
static void write_report(const char* dir, const TriStat* st,
                         const int* alloc, int n, double thresh, long budget,
                         double dmax, double condmax, double smax,
                         Image* img_uv, Image* img_world, Image* img_cond,
                         Image* img_sample, Image* img_checker,
                         long total, double total_area, int flagged,
                         double flagged_area, long flagged_samples,
                         int zero_count, double zero_area, long zero_samples,
                         double max_zero_ratio, double area_diff) {
    char path[1024];
    snprintf(path, sizeof(path), "%s/report.html", dir);
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); return; }

    fprintf(f,
"<!doctype html>\n<html lang=\"zh\"><head><meta charset=\"utf-8\">\n"
"<title>UV 失真报告</title>\n"
"<style>\n"
":root{--surface:#fcfcfb;--plane:#f9f9f7;--ink:#0b0b0b;--ink2:#52514e;"
"--muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;--blue:#2a78d6;--blue-deep:#0d366b;"
"--critical:#d03b3b;--border:rgba(11,11,11,.10)}\n"
"*{box-sizing:border-box}\n"
"body{margin:0;background:var(--plane);color:var(--ink);"
"font:14px/1.5 system-ui,-apple-system,'Segoe UI',sans-serif;padding:32px}\n"
"header{max-width:1200px;margin:0 auto 24px}\n"
"h1{font-size:24px;font-weight:650;margin:0 0 4px}\n"
"h2{font-size:16px;font-weight:600;margin:28px 0 12px}\n"
".sub{color:var(--ink2);margin:0}\n"
"main{max-width:1200px;margin:0 auto}\n"
".tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));"
"gap:12px;margin:20px 0}\n"
".tile{background:var(--surface);border:1px solid var(--border);"
"border-radius:10px;padding:14px 16px}\n"
".tile .k{color:var(--muted);font-size:12px;letter-spacing:.02em}\n"
".tile .v{font-size:26px;font-weight:650;margin-top:2px}\n"
".tile .n{color:var(--ink2);font-size:12px;margin-top:2px}\n"
".tile.bad .v{color:var(--critical)}\n"
".grid{display:grid;grid-template-columns:1fr 1fr;gap:20px}\n"
"@media(max-width:900px){.grid{grid-template-columns:1fr}}\n"
"figure{background:var(--surface);border:1px solid var(--border);"
"border-radius:10px;margin:0;padding:14px}\n"
"figure img{width:100%%;height:auto;display:block;border-radius:4px}\n"
"figcaption{color:var(--ink2);font-size:12.5px;margin-top:10px}\n"
".scale{display:flex;align-items:center;gap:8px;margin-top:12px;"
"font-size:12px;color:var(--ink2)}\n"
".scale-bar{display:flex;flex:1;height:10px;border-radius:5px;overflow:hidden}\n"
".scale-bar i{flex:1}\n"
".scale-min,.scale-max{font-variant-numeric:tabular-nums;min-width:36px}\n"
".scale-max{text-align:right}\n"
"table{width:100%%;border-collapse:collapse;background:var(--surface);"
"border:1px solid var(--border);border-radius:10px;overflow:hidden;"
"font-variant-numeric:tabular-nums;font-size:12.5px}\n"
"th,td{padding:7px 10px;text-align:right;border-bottom:1px solid var(--grid)}\n"
"th{color:var(--muted);font-weight:600;font-size:11.5px;text-align:right;"
"background:var(--surface);position:sticky;top:0}\n"
"th:first-child,td:first-child,th:nth-child(2),td:nth-child(2){text-align:left}\n"
"tbody tr:last-child td{border-bottom:none}\n"
".scroll{max-height:520px;overflow:auto;border-radius:10px;"
"border:1px solid var(--border)}\n"
".scroll table{border:none;border-radius:0}\n"
".tag{display:inline-block;padding:1px 8px;border-radius:999px;"
"font-size:11px;font-weight:600;border:1px solid var(--border)}\n"
".tag.bad{color:var(--critical)}\n"
".tag.ok{color:var(--ink2)}\n"
"code{font:12px ui-monospace,SFMono-Regular,Menlo,monospace;color:var(--ink2)}\n"
"</style></head>\n<body>\n");

    fprintf(f, "<header><h1>UV 失真报告</h1>\n<p class=\"sub\">"
               "%d 个三角形 · 参考密度 %.4g · 阈值 D &gt; %.3g · 采样预算 %ld"
               "</p></header>\n<main>\n", n,
            st[SHEET_REF_TRI].density, thresh, budget);

    /* headline tiles */
    fprintf(f, "<section class=\"tiles\">\n");
    fprintf(f, "<div class=\"tile\"><div class=\"k\">超阈值三角形</div>"
               "<div class=\"v\">%d</div><div class=\"n\">占面积 %.0f%% · 占采样 %.0f%%</div></div>\n",
            flagged, total_area > EPS ? 100.0 * flagged_area / total_area : 0.0,
            total > 0 ? 100.0 * (double)flagged_samples / (double)total : 0.0);
    fprintf(f, "<div class=\"tile\"><div class=\"k\">零失真三角形</div>"
               "<div class=\"v\">%d</div><div class=\"n\">占面积 %.0f%%</div></div>\n",
            zero_count, total_area > EPS ? 100.0 * zero_area / total_area : 0.0);
    fprintf(f, "<div class=\"tile\"><div class=\"k\">零失真采样占比</div>"
               "<div class=\"v\">%.0f%%</div>"
               "<div class=\"n\">面积占比 %.0f%% · 采样/面积比 %.2f%s</div></div>\n",
            total > 0 ? 100.0 * (double)zero_samples / (double)total : 0.0,
            total_area > EPS ? 100.0 * zero_area / total_area : 0.0, max_zero_ratio,
            max_zero_ratio > 1.0 + 1e-6 ? " · 超出" : " · 未超面积份额");
    fprintf(f, "<div class=\"tile\"><div class=\"k\">对拍残差</div>"
               "<div class=\"v\">%.1e</div><div class=\"n\">面积比 / SVD 独立公式</div></div>\n",
            area_diff);
    fprintf(f, "</section>\n");

    /* images in a 2-column grid */
    fprintf(f, "<div class=\"grid\">\n");

    if (img_uv) {
        fprintf(f, "<figure><img alt=\"UV 空间失真热图\" src=\"");
        html_img_uri(f, img_uv);
        fprintf(f, "\"><figcaption><b>UV 空间失真 D</b> —— 左侧 UV 过密、"
                   "右侧稀疏；这就是贴图被拉伸的位置。</figcaption>");
        html_scale_legend(f, 1.0, dmax, "D");
        fprintf(f, "</figure>\n");
    }
    if (img_checker) {
        fprintf(f, "<figure><img alt=\"棋盘纹理在几何上的外观\" src=\"");
        html_img_uri(f, img_checker);
        fprintf(f, "\"><figcaption><b>棋盘纹理</b> —— 同一棋盘经 UV 映射到几何："
                   "格子变窄处过密，变宽处被拉伸。</figcaption></figure>\n");
    }
    if (img_world) {
        fprintf(f, "<figure><img alt=\"世界空间失真热图\" src=\"");
        html_img_uri(f, img_world);
        fprintf(f, "\"><figcaption><b>世界空间失真 D</b> —— 按三角形着色，"
                   "与面积占比对照。</figcaption>");
        html_scale_legend(f, 1.0, dmax, "D");
        fprintf(f, "</figure>\n");
    }
    if (img_sample) {
        fprintf(f, "<figure><img alt=\"采样密度热图\" src=\"");
        html_img_uri(f, img_sample);
        fprintf(f, "\"><figcaption><b>采样密度</b> —— 每单位世界面积分到的采样数；"
                   "失真大的地方被提高，零失真区域保持面积份额。</figcaption>");
        html_scale_legend(f, 0.0, smax, "samples/area");
        fprintf(f, "</figure>\n");
    }
    if (img_cond) {
        fprintf(f, "<figure><img alt=\"条件数热图\" src=\"");
        html_img_uri(f, img_cond);
        fprintf(f, "\"><figcaption><b>各向异性（条件数）</b> —— 检测保面积的剪切失真；"
                   "不参与采样分配，仅报告。</figcaption>");
        html_scale_legend(f, 1.0, condmax, "σmax/σmin");
        fprintf(f, "</figure>\n");
    }
    fprintf(f, "</div>\n");

    /* allocation table by region */
    fprintf(f, "<h2>按区域的采样分配</h2>\n<div class=\"scroll\"><table>\n<thead><tr>"
               "<th>区域</th><th>三角形</th><th>面积占比</th>"
               "<th>采样占比</th><th>面积比</th></tr></thead><tbody>\n");
    {
        const char* regions[3] = { "overdense", "identity", "sparse" };
        const char* labels[3] = { "过密 (s=+2)", "零失真 (s=0)", "稀疏 (s=-2)" };
        for (int r = 0; r < 3; r++) {
            int cnt = 0;
            double a = 0.0;
            long sm = 0;
            for (int i = 0; i < n; i++) {
                if (!st[i].valid || strcmp(sheet_region(i), regions[r])) continue;
                cnt++; a += st[i].a_geo; sm += alloc ? alloc[i] : 0;
            }
            if (!cnt) continue;
            fprintf(f, "<tr><td>%s</td><td>%d</td><td>%.1f%%</td><td>%.1f%%</td>"
                       "<td>%.2f</td></tr>\n", labels[r], cnt,
                    total_area > EPS ? 100.0 * a / total_area : 0.0,
                    total > 0 ? 100.0 * (double)sm / (double)total : 0.0,
                    a > EPS ? ((double)sm / (double)total) / (a / total_area) : 0.0);
        }
    }
    fprintf(f, "</tbody></table></div>\n");

    /* flagged table */
    fprintf(f, "<h2>超阈值三角形（D &gt; %.3g）</h2>\n", thresh);
    fprintf(f, "<div class=\"scroll\"><table>\n<thead><tr>"
               "<th>三角形</th><th>区域</th><th>面积比 r</th><th>log2 拉伸 s</th>"
               "<th>失真 D</th><th>条件数</th><th>采样</th></tr></thead><tbody>\n");
    for (int i = 0; i < n; i++) {
        if (!st[i].valid || st[i].D <= thresh) continue;
        fprintf(f, "<tr><td>%d</td><td>%s</td><td>%.3f</td><td>%+.3f</td>"
                   "<td>%.3f</td><td>%.3f</td><td>%d</td></tr>\n",
                i, sheet_region(i), st[i].r, st[i].s, st[i].D, st[i].cond,
                alloc ? alloc[i] : 0);
    }
    fprintf(f, "</tbody></table></div>\n");

    /* full per-triangle table (table view = WCAG-clean twin of the maps) */
    fprintf(f, "<h2>全部三角形（表格视图）</h2>\n");
    fprintf(f, "<div class=\"scroll\"><table>\n<thead><tr>"
               "<th>tri</th><th>区域</th><th>a_geo</th><th>a_uv</th>"
               "<th>ratio</th><th>log2 s</th><th>D</th><th>cond</th>"
               "<th>采样</th></tr></thead><tbody>\n");
    for (int i = 0; i < n; i++) {
        if (!st[i].valid) continue;
        const char* cls = st[i].D > thresh ? "bad" : "ok";
        fprintf(f, "<tr><td>%d</td><td><span class=\"tag %s\">%s</span></td>"
                   "<td>%.4g</td><td>%.4g</td><td>%.4g</td><td>%+.3f</td>"
                   "<td>%.4g</td><td>%.3g</td><td>%d</td></tr>\n",
                i, cls, sheet_region(i), st[i].a_geo, st[i].a_uv, st[i].r,
                st[i].s, st[i].D, st[i].cond, alloc ? alloc[i] : 0);
    }
    fprintf(f, "</tbody></table></div>\n");
    fprintf(f, "</main></body></html>\n");
    fclose(f);
}

/* ------------------------------------------------------------------ */
/* self-test                                                           */
/* ------------------------------------------------------------------ */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static int near(double a, double b, double tol) { return fabs(a - b) <= tol; }

/* Unit square split into two triangles with a controllable UV rectangle
   and shear. */
static void make_quad(Tri out[2], double u, double v, double skew) {
    V3 p0 = { 0, 0, 0 }, p1 = { 1, 0, 0 }, p2 = { 1, 1, 0 }, p3 = { 0, 1, 0 };
    V2 a = { 0, 0 }, b = { u, skew }, c = { u, v + skew }, d = { 0, v };
    out[0].p[0] = p0; out[0].p[1] = p1; out[0].p[2] = p2;
    out[0].uv[0] = a; out[0].uv[1] = b; out[0].uv[2] = c;
    out[1].p[0] = p0; out[1].p[1] = p2; out[1].p[2] = p3;
    out[1].uv[0] = a; out[1].uv[1] = c; out[1].uv[2] = d;
}

static void test_identity(void) {
    Tri t[2]; TriStat st[2];
    make_quad(t, 1.0, 1.0, 0.0);
    compute_areas(t, 2, st);
    assign_distortion(st, 2, 0.0);
    check("identity_a_geo", near(st[0].a_geo, 0.5, 1e-12) && near(st[1].a_geo, 0.5, 1e-12));
    check("identity_a_uv", near(st[0].a_uv, 0.5, 1e-12) && near(st[1].a_uv, 0.5, 1e-12));
    check("identity_density_one", near(st[0].density, 1.0, 1e-12));
    check("identity_zero_distortion", st[0].s == 0.0 && st[1].s == 0.0);
    check("identity_D_one", st[0].D == 1.0 && st[1].D == 1.0);
    check("identity_isotropic", near(st[0].cond, 1.0, 1e-12));
    check("identity_area_paths_agree",
          near(area_geo(&t[0]), area_geo_heron(&t[0]), 1e-12) &&
          near(area_uv(&t[0]), area_uv_shoelace(&t[0]), 1e-12));
}

static void test_area_ratio(void) {
    Tri t[2]; TriStat st[2];
    make_quad(t, 2.0, 1.0, 0.0);          /* UV 2x1 -> density 2 */
    compute_areas(t, 2, st);
    assign_distortion(st, 2, 1.0);
    check("stretch2_ratio", near(st[0].r, 2.0, 1e-12));
    check("stretch2_s_plus_one", near(st[0].s, 1.0, 1e-12));
    check("stretch2_D_two", near(st[0].D, 2.0, 1e-12));
    check("stretch2_singular", near(st[0].stretch_max, 1.0, 1e-12) &&
                               near(st[0].stretch_min, 0.5, 1e-12));

    make_quad(t, 0.5, 0.5, 0.0);          /* UV 0.5x0.5 -> density 0.25 */
    compute_areas(t, 2, st);
    assign_distortion(st, 2, 1.0);
    check("shrink4_ratio", near(st[0].r, 0.25, 1e-12));
    check("shrink4_s_minus_two", near(st[0].s, -2.0, 1e-12));
    check("shrink4_D_four", near(st[0].D, 4.0, 1e-12));

    /* The distortion equals the area ratio recomputed by the independent
       Heron/shoelace path, and the density equals 1/|det J|, which the
       tangent SVD supplies through smax*smin. */
    double r_analytic = area_uv_shoelace(&t[0]) / area_geo_heron(&t[0]) / 1.0;
    check("analytic_ratio_matches", near(st[0].r, r_analytic, 1e-12));
    check("ratio_matches_svd_det",
          near(st[0].density, 1.0 / (st[0].stretch_max * st[0].stretch_min), 1e-12));
}

static void test_anisotropy(void) {
    Tri t[2]; TriStat st[2];
    make_quad(t, 2.0, 0.5, 0.0);          /* area ratio 1, cond 4 */
    compute_areas(t, 2, st);
    assign_distortion(st, 2, 1.0);
    check("aniso_ratio_one", near(st[0].r, 1.0, 1e-12));
    check("aniso_s_zero", st[0].s == 0.0);
    check("aniso_cond_four", near(st[0].cond, 4.0, 1e-9));

    make_quad(t, 1.0, 1.0, 0.5);          /* det-1 shear, area unchanged */
    compute_areas(t, 2, st);
    assign_distortion(st, 2, 1.0);
    check("shear_area_one", near(st[0].a_uv, 0.5, 1e-12));
    check("shear_s_zero", st[0].s == 0.0);
    check("shear_cond_raised", st[0].cond > 1.05);
}

static void test_reference_semantics(void) {
    /* Same geometry, UV areas 1:4 -> s = 0 and +2 against ref = 1. */
    Tri t[2];
    V3 p0 = { 0, 0, 0 }, p1 = { 1, 0, 0 }, p2 = { 1, 1, 0 }, p3 = { 0, 1, 0 };
    V2 a = { 0, 0 }, b = { 1, 0 }, c = { 1, 1 };
    V2 e = { 2, 0 }, f = { 2, 2 };
    t[0].p[0] = p0; t[0].p[1] = p1; t[0].p[2] = p2;
    t[0].uv[0] = a; t[0].uv[1] = b; t[0].uv[2] = c;
    t[1].p[0] = p0; t[1].p[1] = p2; t[1].p[2] = p3;
    t[1].uv[0] = a; t[1].uv[1] = e; t[1].uv[2] = f;
    TriStat st[2];
    compute_areas(t, 2, st);
    assign_distortion(st, 2, 1.0);
    check("ref_density_one", near(st[0].density, 1.0, 1e-12));
    check("ref_density_four", near(st[1].density, 4.0, 1e-12));
    check("ref_s_zero", st[0].s == 0.0);
    check("ref_s_plus_two", near(st[1].s, 2.0, 1e-12));
    check("ref_D_four", near(st[1].D, 4.0, 1e-12));

    /* Default reference is the area-weighted mean. */
    TriStat st2[2];
    compute_areas(t, 2, st2);
    double ref = assign_distortion(st2, 2, 0.0);
    check("ref_default_total", near(ref, (st2[0].a_uv + st2[1].a_uv) /
                                        (st2[0].a_geo + st2[1].a_geo), 1e-12));
}

static void test_allocation(void) {
    /* All-identity mesh: exact area split, nothing oversampled. */
    {
        Tri t[2]; TriStat st[2];
        make_quad(t, 1.0, 1.0, 0.0);
        compute_areas(t, 2, st);
        assign_distortion(st, 2, 1.0);
        int a[2];
        check("alloc_uniform_ok", allocate_samples(st, 2, 1000, 1.0, a) == 0);
        check("alloc_uniform_split", a[0] == 500 && a[1] == 500);
        check("alloc_uniform_total", a[0] + a[1] == 1000);
    }
    /* One identity, one 4x: 20% / 80%, deterministic. */
    {
        Tri t[2];
        V3 p0 = { 0, 0, 0 }, p1 = { 1, 0, 0 }, p2 = { 1, 1, 0 }, p3 = { 0, 1, 0 };
        V2 a = { 0, 0 }, b = { 1, 0 }, c = { 1, 1 };
        V2 e = { 2, 0 }, f = { 2, 2 };
        t[0].p[0] = p0; t[0].p[1] = p1; t[0].p[2] = p2;
        t[0].uv[0] = a; t[0].uv[1] = b; t[0].uv[2] = c;
        t[1].p[0] = p0; t[1].p[1] = p2; t[1].p[2] = p3;
        t[1].uv[0] = a; t[1].uv[1] = e; t[1].uv[2] = f;
        TriStat st[2];
        compute_areas(t, 2, st);
        assign_distortion(st, 2, 1.0);
        int al[2], al2[2];
        check("alloc_skew_ok", allocate_samples(st, 2, 1000, 1.0, al) == 0);
        allocate_samples(st, 2, 1000, 1.0, al2);
        check("alloc_skew_expected", al[0] == 200 && al[1] == 800);
        check("alloc_skew_total", al[0] + al[1] == 1000);
        check("alloc_deterministic", al[0] == al2[0] && al[1] == al2[1]);
        check("alloc_zero_not_oversampled", al[0] <= 500);
    }
    /* Three triangles: two identity (equal area 0.5) and one 4x-density
       (area 0.5). Weights are 0.5, 0.5, 2.0 -> shares 1/6, 1/6, 2/3.
       Largest-remainder gives 167, 167, 666. */
    {
        Tri t[3];
        V3 p0 = { 0, 0, 0 }, p1 = { 1, 0, 0 }, p2 = { 1, 1, 0 };
        V3 p3 = { 0, 1, 0 }, p4 = { 2, 1, 0 };
        V2 a = { 0, 0 }, b = { 1, 0 }, c = { 1, 1 }, d = { 0, 1 };
        V2 e = { 3, 0 }, f = { 3, 2 };
        t[0].p[0] = p0; t[0].p[1] = p1; t[0].p[2] = p2;
        t[0].uv[0] = a; t[0].uv[1] = b; t[0].uv[2] = c;
        t[1].p[0] = p0; t[1].p[1] = p2; t[1].p[2] = p3;
        t[1].uv[0] = a; t[1].uv[1] = c; t[1].uv[2] = d;
        t[2].p[0] = p1; t[2].p[1] = p4; t[2].p[2] = p2;
        t[2].uv[0] = b; t[2].uv[1] = e; t[2].uv[2] = f;
        TriStat st[3];
        compute_areas(t, 3, st);
        assign_distortion(st, 3, 1.0);
        check("tri3_valid", st[0].valid && st[1].valid && st[2].valid);
        check("tri3_zero_pair", st[0].s == 0.0 && st[1].s == 0.0);
        check("tri3_distorted", near(st[2].s, 2.0, 1e-12));
        check("tri3_distorted_density", near(st[2].density, 4.0, 1e-12));
        int al[3];
        allocate_samples(st, 3, 1000, 1.0, al);
        long zero = al[0] + al[1];
        check("tri3_zero_share", zero == 334);
        check("tri3_distorted_share", al[2] == 666);
        check("tri3_budget", zero + al[2] == 1000);
        check("tri3_zero_below_area_share", zero < 2.0 * 1000.0 / 3.0);
    }
}

static void test_sheet(void) {
    Tri* tris = (Tri*)malloc(sizeof(Tri) * SHEET_TRI);
    TriStat* st = (TriStat*)malloc(sizeof(TriStat) * SHEET_TRI);
    int* alloc = (int*)malloc(sizeof(int) * SHEET_TRI);
    if (!tris || !st || !alloc) {
        check("sheet_alloc_ok", 0);
        free(tris); free(st); free(alloc);
        return;
    }
    build_sheet(tris);
    compute_areas(tris, SHEET_TRI, st);
    double ref = st[SHEET_REF_TRI].density;
    assign_distortion(st, SHEET_TRI, ref);

    int invalid = 0, bad_identity = 0, bad_band = 0;
    int identity_n = 0, overdense_n = 0, sparse_n = 0;
    for (int i = 0; i < SHEET_TRI; i++) {
        if (!st[i].valid) { invalid++; continue; }
        const char* reg = sheet_region(i);
        if (!strcmp(reg, "identity")) {
            identity_n++;
            if (st[i].s != 0.0) bad_identity++;
        } else if (!strcmp(reg, "overdense")) {
            overdense_n++;
            if (!near(st[i].s, 2.0, 1e-6)) bad_band++;
        } else {
            sparse_n++;
            if (!near(st[i].s, -2.0, 1e-6)) bad_band++;
        }
    }
    check("sheet_all_valid", invalid == 0);
    check("sheet_identity_zero", bad_identity == 0 && identity_n == 2 * 6 * SHEET_CY);
    check("sheet_bands_exact", bad_band == 0 && overdense_n == 2 * 2 * SHEET_CY &&
                               sparse_n == 2 * 4 * SHEET_CY);
    check("sheet_ref_is_zero", st[SHEET_REF_TRI].s == 0.0);

    check("sheet_alloc_ok", allocate_samples(st, SHEET_TRI, 96000, 1.0, alloc) == 0);
    long total = 0, zero_samples = 0;
    double zero_area = 0.0, total_area = 0.0;
    for (int i = 0; i < SHEET_TRI; i++) {
        total += alloc[i];
        total_area += st[i].a_geo;
        if (st[i].s == 0.0) { zero_samples += alloc[i]; zero_area += st[i].a_geo; }
    }
    check("sheet_budget", total == 96000);
    check("sheet_zero_area_half", near(zero_area / total_area, 0.5, 1e-12));
    check("sheet_zero_sample_share", near((double)zero_samples / (double)total, 0.2, 1e-9));
    check("sheet_zero_not_oversampled",
          (double)zero_samples / (double)total < zero_area / total_area);

    /* Every area matches the independent formula path, and every density
       matches 1/(smax*smin) from the tangent SVD. */
    double maxdiff = 0.0, maxsvd = 0.0;
    for (int i = 0; i < SHEET_TRI; i++) {
        double dg = fabs(area_geo(&tris[i]) - area_geo_heron(&tris[i]));
        double du = fabs(area_uv(&tris[i]) - area_uv_shoelace(&tris[i]));
        if (dg > maxdiff) maxdiff = dg;
        if (du > maxdiff) maxdiff = du;
        double from_svd = 1.0 / (st[i].stretch_max * st[i].stretch_min);
        double d = fabs(st[i].density - from_svd);
        if (d > maxsvd) maxsvd = d;
    }
    check("sheet_area_paths_agree", maxdiff < 1e-12);
    check("sheet_density_matches_svd", maxsvd < 1e-9);

    /* Anisotropy is present (global shear) even where s == 0. */
    double cond_id = 0.0;
    for (int i = 0; i < SHEET_TRI; i++)
        if (st[i].s == 0.0 && st[i].cond > cond_id) cond_id = st[i].cond;
    check("sheet_shear_anisotropy", cond_id > 1.0);

    free(tris); free(st); free(alloc);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char** argv) {
    int self = 0, do_export = 0;
    const char* dir = ".";
    double thresh = DEFAULT_THRESHOLD;
    long budget = 96000;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--selftest")) self = 1;
        else if (!strcmp(argv[i], "--export")) {
            do_export = 1;
            if (i + 1 < argc && argv[i + 1][0] != '-') dir = argv[++i];
        } else if (!strcmp(argv[i], "--threshold") && i + 1 < argc) {
            thresh = atof(argv[++i]);
        } else if (!strcmp(argv[i], "--budget") && i + 1 < argc) {
            budget = atol(argv[++i]);
        }
    }

    if (self) {
        test_identity();
        test_area_ratio();
        test_anisotropy();
        test_reference_semantics();
        test_allocation();
        test_sheet();
        printf("all_ok=%d\n", failures == 0);
        return failures == 0 ? 0 : 1;
    }

    if (do_export) {
        Tri* tris = (Tri*)malloc(sizeof(Tri) * SHEET_TRI);
        TriStat* st = (TriStat*)malloc(sizeof(TriStat) * SHEET_TRI);
        int* alloc = (int*)malloc(sizeof(int) * SHEET_TRI);
        if (!tris || !st || !alloc) { fprintf(stderr, "out of memory\n"); return 1; }
        build_sheet(tris);
        compute_areas(tris, SHEET_TRI, st);
        double ref = st[SHEET_REF_TRI].density;
        assign_distortion(st, SHEET_TRI, ref);
        if (allocate_samples(st, SHEET_TRI, budget, 1.0, alloc) != 0) {
            fprintf(stderr, "allocation failed\n");
            return 1;
        }
        int rc = export_all(dir, tris, st, alloc, budget, SHEET_TRI, thresh);
        free(tris); free(st); free(alloc);
        return rc == 0 ? 0 : 1;
    }

    fprintf(stderr,
            "usage: %s --selftest | --export <dir> [--threshold D] [--budget N]\n",
            argv[0]);
    return 2;
}
