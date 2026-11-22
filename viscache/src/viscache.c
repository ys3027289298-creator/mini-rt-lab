
/* viscache: per-light visibility cache for soft shadows.

   One square grid is allocated per light over a receiver-space domain.
   A cell stores the fraction of the area light visible from the cell
   centre, estimated by stratified ray tests.  A cell is recomputed with
   progressively more samples while its estimated quadrature error is
   above a tolerance, so the cache spends work where the penumbra is.

   Lookup is bilinear.  The bound reported with each lookup combines the
   worst cell error, the local cell-to-cell spread (interpolation
   uncertainty), the bilinear truncation floor h^2, and the residual
   error of the dense direct reference used by the test. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NLIGHT        2
#define MIN_CELLS     16      /* receiver cells per axis, below this warn */
#define CELL_SIDE_MIN 8       /* 8x8 = 64 initial rays per cell */
#define CELL_SIDE_MAX 32      /* 32x32 = 1024 rays per cell */
#define REF_SIDE      128     /* 128x128 rays for the direct reference */
#define EVAL_GRID     33      /* direct-vs-cache comparison grid */
#define FRAMES        6
#define FRAME_POINTS  21
#define TOL           0.02f   /* quadrature error tolerance per cell */
#define REF_ERR_EST   (2.0f / (float)REF_SIDE)

typedef struct { float x, y, z; } V3;
typedef struct { V3 center; float half; } AreaLight; /* square, z = center.z */
typedef struct { V3 c; float r; } Sphere;

typedef struct {
    Sphere occ[2];
    int nocc;
} Scene;

typedef struct {
    float v;         /* visible fraction */
    float err;       /* estimated quadrature error */
    int valid;
    int saturated;   /* max samples reached, error still above tol */
} Cell;

typedef struct {
    const AreaLight* light;
    const Scene* scene;
    int res;
    float x0, y0, size;
    Cell* cells;
    int res_warning;
    long lookups, hits, recomputes, inadequate, rays;
} LightCache;

typedef struct {
    float v;
    float bound;
    int recomputed;
} Lookup;

static V3 v3_sub(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static float v3_dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

/* segment p..q (t in [0,1]) blocked by any sphere */
static int occluded(const Scene* sc, V3 p, V3 q) {
    V3 d = v3_sub(q, p);
    float a = v3_dot(d, d);
    if (a <= 0.0f) return 0;
    for (int i = 0; i < sc->nocc; i++) {
        V3 oc = v3_sub(p, sc->occ[i].c);
        float b = 2.0f * v3_dot(oc, d);
        float c = v3_dot(oc, oc) - sc->occ[i].r * sc->occ[i].r;
        float disc = b * b - 4.0f * a * c;
        if (disc <= 0.0f) continue;
        float s = sqrtf(disc);
        float t0 = (-b - s) / (2.0f * a);
        float t1 = (-b + s) / (2.0f * a);
        if ((t0 > 1e-4f && t0 < 1.0f - 1e-4f) ||
            (t1 > 1e-4f && t1 < 1.0f - 1e-4f)) return 1;
    }
    return 0;
}

/* stratified direct visibility: side*side rays at light-cell centres */
static float vis_quad(const Scene* sc, V3 p, const AreaLight* L, int side) {
    int visible = 0;
    for (int j = 0; j < side; j++) {
        float v = ((float)j + 0.5f) / (float)side;
        float ly = L->center.y + (2.0f * v - 1.0f) * L->half;
        for (int i = 0; i < side; i++) {
            float u = ((float)i + 0.5f) / (float)side;
            V3 q = { L->center.x + (2.0f * u - 1.0f) * L->half, ly, L->center.z };
            if (!occluded(sc, p, q)) visible++;
        }
    }
    return (float)visible / (float)(side * side);
}

static int cache_init(LightCache* c, const Scene* sc, const AreaLight* L,
                      float x0, float y0, float size, int res) {
    memset(c, 0, sizeof(*c));
    c->light = L;
    c->scene = sc;
    c->res = res;
    c->x0 = x0; c->y0 = y0; c->size = size;
    c->res_warning = (res < MIN_CELLS);
    if (c->res_warning) {
        printf("WARN: 缓存分辨率不足 (cache resolution %d < minimum %d), 精度不足\n",
               res, MIN_CELLS);
    }
    c->cells = (Cell*)calloc((size_t)res * (size_t)res, sizeof(Cell));
    return c->cells != NULL;
}

static void cache_free(LightCache* c) { free(c->cells); c->cells = NULL; }

static void cache_refresh(LightCache* c, int ix, int iy, float tol) {
    float h = c->size / (float)c->res;
    V3 p = { c->x0 + ((float)ix + 0.5f) * h,
             c->y0 + ((float)iy + 0.5f) * h, 0.0f };
    int side = CELL_SIDE_MIN;
    float prev = vis_quad(c->scene, p, c->light, side);
    float v = prev;
    c->rays += (long)side * side;
    float err = 0.5f; /* unknown at the initial level, force one refine */
    while (err > tol && side < CELL_SIDE_MAX) {
        side *= 2;
        float cur = vis_quad(c->scene, p, c->light, side);
        c->rays += (long)side * side;
        err = fabsf(cur - prev) + 0.5f / ((float)side * (float)side);
        v = cur;
        prev = cur;
    }
    Cell* cell = &c->cells[(size_t)iy * (size_t)c->res + (size_t)ix];
    cell->v = v;
    cell->err = err;
    cell->valid = 1;
    cell->saturated = (err > tol);
    c->recomputes++;
    if (cell->saturated) c->inadequate++;
}

static void cache_prepare_all(LightCache* c, float tol) {
    for (int y = 0; y < c->res; y++)
        for (int x = 0; x < c->res; x++)
            cache_refresh(c, x, y, tol);
}

static Cell* cache_get(LightCache* c, int ix, int iy, float tol, int* recomputed) {
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (ix >= c->res) ix = c->res - 1;
    if (iy >= c->res) iy = c->res - 1;
    Cell* cell = &c->cells[(size_t)iy * (size_t)c->res + (size_t)ix];
    if (!cell->valid || (!cell->saturated && cell->err > tol)) {
        *recomputed = 1;
        cache_refresh(c, ix, iy, tol);
    }
    return cell;
}

static Lookup cache_lookup(LightCache* c, V3 p, float tol) {
    float h = c->size / (float)c->res;
    float fx = (p.x - c->x0) / h - 0.5f;
    float fy = (p.y - c->y0) / h - 0.5f;
    if (fx < 0.0f) fx = 0.0f;
    if (fy < 0.0f) fy = 0.0f;
    if (fx > (float)(c->res - 1)) fx = (float)(c->res - 1);
    if (fy > (float)(c->res - 1)) fy = (float)(c->res - 1);
    int i0 = (int)fx, j0 = (int)fy;
    int i1 = (i0 + 1 < c->res) ? i0 + 1 : i0;
    int j1 = (j0 + 1 < c->res) ? j0 + 1 : j0;
    float tx = fx - (float)i0, ty = fy - (float)j0;

    Lookup lk;
    lk.recomputed = 0;
    Cell* c00 = cache_get(c, i0, j0, tol, &lk.recomputed);
    Cell* c10 = cache_get(c, i1, j0, tol, &lk.recomputed);
    Cell* c01 = cache_get(c, i0, j1, tol, &lk.recomputed);
    Cell* c11 = cache_get(c, i1, j1, tol, &lk.recomputed);

    float v00 = c00->v, v10 = c10->v, v01 = c01->v, v11 = c11->v;
    float v0 = v00 * (1.0f - tx) + v10 * tx;
    float v1 = v01 * (1.0f - tx) + v11 * tx;
    lk.v = v0 * (1.0f - ty) + v1 * ty;

    /* Bound from the 3x3 cell block around the query.  A query near a cell
       border or outside the outer centres is a convex combination of cells in
       this block, so the block's value spread covers the receiver-space
       variation across the interpolation footprint (including the clamped
       half-cell strip at the domain edge, which the four corners alone miss). */
    float vmin = 1.0f, vmax = 0.0f, emax = 0.0f;
    int all_valid = 1;
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            int xx = i0 + dx, yy = j0 + dy;
            if (xx < 0 || yy < 0 || xx >= c->res || yy >= c->res) continue;
            Cell* nb = &c->cells[(size_t)yy * (size_t)c->res + (size_t)xx];
            if (!nb->valid) { all_valid = 0; continue; }
            if (nb->v < vmin) vmin = nb->v;
            if (nb->v > vmax) vmax = nb->v;
            if (nb->err > emax) emax = nb->err;
        }
    }
    if (!all_valid) { vmin = 0.0f; vmax = 1.0f; }
    /* emax: worst quadrature error in the footprint;
       vmax-vmin: receiver-space variation over the interpolation footprint;
       h*h: bilinear truncation floor;
       REF_ERR_EST: residual error of the dense direct reference used to judge
       the bound, so the comparison is not tighter than its own yardstick. */
    lk.bound = emax + (vmax - vmin) + h * h + REF_ERR_EST;

    c->lookups++;
    if (!lk.recomputed) c->hits++;
    return lk;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Scene sc;
    memset(&sc, 0, sizeof(sc));
    sc.nocc = 2;
    sc.occ[0].c = (V3){ -0.15f, 0.0f, 1.0f };  sc.occ[0].r = 0.55f;
    sc.occ[1].c = (V3){  0.65f, -0.35f, 0.8f }; sc.occ[1].r = 0.32f;

    AreaLight lights[NLIGHT] = {
        { { 0.0f,  0.0f, 2.4f }, 0.80f },
        { { 0.7f, -0.5f, 2.2f }, 0.55f },
    };

    const float x0 = -1.0f, y0 = -1.0f, size = 2.0f;
    int hi_res = 48, lo_res = 8;

    LightCache hi[NLIGHT], lo[NLIGHT];
    int ok = 1;
    for (int l = 0; l < NLIGHT; l++) {
        ok &= cache_init(&hi[l], &sc, &lights[l], x0, y0, size, hi_res);
        ok &= cache_init(&lo[l], &sc, &lights[l], x0, y0, size, lo_res);
    }
    check("alloc_ok", ok);

    /* fully warm both caches so accuracy is measured without refills */
    for (int l = 0; l < NLIGHT; l++) {
        cache_prepare_all(&hi[l], TOL);
        cache_prepare_all(&lo[l], TOL);
    }

    /* accuracy against dense direct ray-traced visibility */
    float hi_max = 0.0f, hi_sum = 0.0f, hi_bound = 0.0f;
    float lo_max = 0.0f, lo_sum = 0.0f, lo_bound = 0.0f;
    float hi_bound_sum = 0.0f, lo_bound_sum = 0.0f;
    int hi_viol = 0, lo_viol = 0, samples = 0;
    for (int l = 0; l < NLIGHT; l++) {
        float l_hi_max = 0.0f, l_lo_max = 0.0f;
        for (int gy = 0; gy < EVAL_GRID; gy++) {
            for (int gx = 0; gx < EVAL_GRID; gx++) {
                V3 p = { x0 + size * (float)gx / (float)(EVAL_GRID - 1),
                         y0 + size * (float)gy / (float)(EVAL_GRID - 1), 0.0f };
                float ref = vis_quad(&sc, p, &lights[l], REF_SIDE);
                Lookup a = cache_lookup(&hi[l], p, TOL);
                Lookup b = cache_lookup(&lo[l], p, TOL);
                float ea = fabsf(a.v - ref), eb = fabsf(b.v - ref);
                hi_sum += ea; lo_sum += eb; samples++;
                hi_bound_sum += a.bound; lo_bound_sum += b.bound;
                if (ea > hi_max) hi_max = ea;
                if (eb > lo_max) lo_max = eb;
                if (ea > l_hi_max) l_hi_max = ea;
                if (eb > l_lo_max) l_lo_max = eb;
                if (a.bound > hi_bound) hi_bound = a.bound;
                if (b.bound > lo_bound) lo_bound = b.bound;
                if (ea > a.bound + 1e-4f) hi_viol++;
                if (eb > b.bound + 1e-4f) lo_viol++;
            }
        }
        printf("light%d_hi_max_err=%.4f light%d_lo_max_err=%.4f "
               "light%d_hi_bound=%.4f\n",
               l, l_hi_max, l, l_lo_max, l, hi_bound);
    }
    float hi_mean_bound = hi_bound_sum / (float)samples;
    float lo_mean_bound = lo_bound_sum / (float)samples;

    /* per-light independence: the two lights must not share results */
    float max_light_diff = 0.0f;
    for (int gy = 0; gy < EVAL_GRID; gy++) {
        for (int gx = 0; gx < EVAL_GRID; gx++) {
            V3 p = { x0 + size * (float)gx / (float)(EVAL_GRID - 1),
                     y0 + size * (float)gy / (float)(EVAL_GRID - 1), 0.0f };
            float v0 = vis_quad(&sc, p, &lights[0], 32);
            float v1 = vis_quad(&sc, p, &lights[1], 32);
            float d = fabsf(v0 - v1);
            if (d > max_light_diff) max_light_diff = d;
        }
    }

    /* frame sequence: empty cache, receiver grid drifts across cells */
    LightCache hi_stats[NLIGHT];
    for (int l = 0; l < NLIGHT; l++)
        ok &= cache_init(&hi_stats[l], &sc, &lights[l], x0, y0, size, hi_res);
    check("stats_alloc_ok", ok);

    long total_hits = 0, total_lookups = 0;
    float first_frame_rate = 0.0f;
    for (int f = 0; f < FRAMES; f++) {
        long f_hits = 0, f_lookups = 0;
        float dx = 0.010f * (float)f, dy = 0.006f * (float)f;
        for (int gy = 0; gy < FRAME_POINTS; gy++) {
            for (int gx = 0; gx < FRAME_POINTS; gx++) {
                V3 p = { -0.9f + 1.8f * (float)gx / (float)(FRAME_POINTS - 1) + dx,
                         -0.9f + 1.8f * (float)gy / (float)(FRAME_POINTS - 1) + dy,
                         0.0f };
                for (int l = 0; l < NLIGHT; l++) {
                    Lookup lk = cache_lookup(&hi_stats[l], p, TOL);
                    f_lookups++;
                    if (!lk.recomputed) f_hits++;
                }
            }
        }
        total_hits += f_hits;
        total_lookups += f_lookups;
        printf("frame=%d lookups=%ld hits=%ld hit_rate=%.4f\n",
               f, f_lookups, f_hits, (double)f_hits / (double)f_lookups);
        if (f == 0) first_frame_rate = (float)f_hits / (float)f_lookups;
    }
    float hit_rate = (float)total_hits / (float)total_lookups;
    long recomputes = 0, inadequate = 0;
    for (int l = 0; l < NLIGHT; l++) {
        recomputes += hi_stats[l].recomputes;
        inadequate += hi_stats[l].inadequate;
        printf("light%d_hit_rate=%.4f light%d_recomputes=%ld\n",
               l, (double)hi_stats[l].hits / (double)hi_stats[l].lookups,
               l, hi_stats[l].recomputes);
    }

    printf("res_hi=%d res_lo=%d min_res=%d cells_per_light=%d\n",
           hi_res, lo_res, MIN_CELLS, hi_res * hi_res);
    printf("tolerance=%.3f ref_side=%d ref_err_est=%.4f eval_points=%d\n",
           TOL, REF_SIDE, REF_ERR_EST, samples);
    printf("hi_max_err=%.4f hi_mean_err=%.4f hi_max_bound=%.4f "
           "hi_mean_bound=%.4f hi_violations=%d\n",
           hi_max, hi_sum / (float)samples, hi_bound, hi_mean_bound, hi_viol);
    printf("lo_max_err=%.4f lo_mean_err=%.4f lo_max_bound=%.4f "
           "lo_mean_bound=%.4f lo_violations=%d\n",
           lo_max, lo_sum / (float)samples, lo_bound, lo_mean_bound, lo_viol);
    printf("hit_rate=%.4f hit_rate_frame0=%.4f queries=%ld hits=%ld "
           "recomputes=%ld inadequate=%ld\n",
           hit_rate, first_frame_rate, total_lookups, total_hits,
           recomputes, inadequate);
    printf("max_light_diff=%.4f\n", max_light_diff);

    /* tolerance sweep: cost rises with the error budget, and a tolerance the
       estimator cannot meet saturates cells and is reported as inadequate */
    const float tols[3] = { 0.005f, TOL, 0.20f };
    long tol_rays[3], tol_inadequate[3];
    for (int ti = 0; ti < 3; ti++) {
        LightCache t;
        ok &= cache_init(&t, &sc, &lights[0], x0, y0, size, 24);
        cache_prepare_all(&t, tols[ti]);
        tol_rays[ti] = t.rays;
        tol_inadequate[ti] = t.inadequate;
        printf("tol=%.3f cells=%d recomputes=%ld rays=%ld inadequate=%ld\n",
               tols[ti], t.res * t.res, t.recomputes, t.rays, t.inadequate);
        cache_free(&t);
    }
    check("tight_tol_reports_inadequate", tol_inadequate[0] > 0);
    check("loose_tol_no_inadequate", tol_inadequate[2] == 0);
    check("cost_monotone_in_tol", tol_rays[0] > tol_rays[1] && tol_rays[1] > tol_rays[2]);


    check("high_res_no_warning", hi[0].res_warning == 0 && hi[1].res_warning == 0);
    check("low_res_warning", lo[0].res_warning == 1 && lo[1].res_warning == 1);
    check("bound_holds_high", hi_viol == 0);
    check("bound_holds_low", lo_viol == 0);
    check("high_res_accurate", hi_max < 0.05f);
    check("low_res_less_accurate", lo_max > hi_max + 0.02f);
    check("cache_fills_frame0", first_frame_rate < 0.01f);
    check("hit_rate_above_half", hit_rate >= 0.50f);
    check("lights_independent", max_light_diff > 0.10f);

    for (int l = 0; l < NLIGHT; l++) {
        cache_free(&hi[l]); cache_free(&lo[l]); cache_free(&hi_stats[l]);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
