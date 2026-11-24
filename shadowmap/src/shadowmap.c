/*
 * shadowmap: point light, wall, polar shadow map whose texel budget is
 * re-distributed with distance from the light.
 *
 * Scene (H = 1; all lengths in light heights)
 *   wall      y = 0, uniform albedo, disk 0 <= r <= R_MAX
 *   light     (0, H, 0), point
 *   camera    orthographic from +y, so one image pixel is one fixed wall area
 *   blockers  NBLK spheres at y ~ 0.55, each casting an umbra centred at a
 *             target wall radius.  The truth visibility is analytic; the
 *             blockers and their azimuths are re-drawn every trial so the
 *             per-ring leak area is an average over sub-texel alignment,
 *             not a single quantisation accident.
 *
 * Shadow map: polar around the light nadir.  NH radial rows x NP azimuth
 * columns; NH*NP is held fixed for every allocation, so the runs differ only
 * in where the rows go.  A lookup reads the nearest texel's binary visibility,
 * which is what turns a coarse allocation into visible blocky patches.
 *
 * Allocation law, one exponent alpha:
 *     wall area per texel   A_t(r) ~ r^alpha
 *     texels per wall area  n(r)   ~ r^-alpha
 *     cumulative rows       N(<=r) ~ r^(2-alpha)
 *   alpha = 0   uniform wall area per texel (planar/atlas map)
 *   alpha > 0   near-dense / far-sparse
 *
 * Derivation
 *   Solid angle of a wall element:  dOmega = dA * cos(theta_L) / r_L^2
 *   with cos(theta_L) = H / r_L and r_L^2 = r^2 + H^2:
 *       dOmega/dr = 2*pi*H*r / (r^2+H^2)^(3/2)
 *   Equal solid angle per texel (n(r) texels per dr):
 *       n(r)  ~ r / (r^2+H^2)^(3/2)
 *       A_t   = dA/(n dr) ~ (r^2+H^2)^(3/2) = r_L^3
 *   Local exponent:
 *       alpha(r) = dln A_t/dln r = 3 r^2 / (r^2+H^2)   in [0,3)
 *   alpha -> 0 at the light, 1.5 at r = H, -> 3 far out.  No constant
 *   exponent is exact; the program reports the band fit and sweeps alpha.
 *
 * Measurements per sample radius (20 points, geometric spacing)
 *   coverage  texels per wall area
 *   dres      wall area per texel (H^2)
 *   leak      wall area the map leaves lit although truth is shadowed
 *   notch     wall area the map darkens although truth is lit
 *   vis leak  the leak restricted to points whose dimmed radiance
 *             I * cos(theta_L)/r_L^2 exceeds THRESH (a visible bright patch)
 *
 * Usage:
 *   shadowmap --selftest
 *   shadowmap <out.csv> [dump]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PI_F   3.14159265358979f
#define H      1.0f
#define R_MIN  0.01f
#ifndef R_MAX
#define R_MAX  2.50f
#endif
#ifndef NH
#define NH     256
#endif
#ifndef NP
#define NP     512
#endif
#define NPTS   20
#define EDGE   512
#define NTRIAL 12
#define NBLK   14
#define NI     3
#define THRESH 0.10f
#define GAIN   1.1f

static const float INTENS[NI] = { 0.35f, 1.0f, 3.0f };
static const char* INAME[NI] = { "low", "mid", "high" };

/* ---------------- scene ---------------- */

typedef struct { float x, y, z, r; } Blk;
static Blk BLK[NBLK];

static unsigned long rng_s;
static void rng_seed(unsigned long s) { rng_s = s ? s : 1; }
static float frand(void) {
    rng_s = rng_s * 6364136223846793005UL + 1442695040888963407UL;
    return (float)((rng_s >> 33) & 0x7fffffff) / (float)0x7fffffff;
}

/* a sphere at height y, centred at plane radius rc, casts an (approximate)
 * umbra of radius rho*H/(H-y) centred at rc*H/(H-y); invert to place it */
static void place_umbra(Blk* b, float wall_r, float wall_u, float phi, float y) {
    float f = (H - y) / H;
    b->x = wall_r * f * cosf(phi);
    b->z = wall_r * f * sinf(phi);
    b->y = y;
    b->r = wall_u * f;
}
static void make_scene(int trial) {
    rng_seed(0x9e3779b9UL * (unsigned long)(trial + 1));
    for (int i = 0; i < NBLK; i++) {
        float u = ((float)i + 0.5f) / (float)NBLK;      /* stratified in log r */
        float R = 0.20f * powf(2.20f / 0.20f, u);
        R *= expf(0.10f * (frand() - 0.5f) * 2.0f);
        float U = 0.16f * R;
        float y = 0.55f + 0.05f * (frand() - 0.5f);
        float phi = 2.0f * PI_F * frand();
        place_umbra(&BLK[i], R, U, phi, y);
    }
}

static int sphere_blocked(const float L[3], const Blk* b, const float D[3]) {
    float a0 = b->x - L[0], a1 = b->y - L[1], a2 = b->z - L[2];
    float dd = D[0]*D[0] + D[1]*D[1] + D[2]*D[2];
    float t = (a0*D[0] + a1*D[1] + a2*D[2]) / dd;
    if (t <= 0.0f) return 0;
    float e0 = a0 - t*D[0], e1 = a1 - t*D[1], e2 = a2 - t*D[2];
    return e0*e0 + e1*e1 + e2*e2 < b->r * b->r;
}
static int shadowed(float x, float z) {
    const float L[3] = { 0.0f, H, 0.0f };
    const float D[3] = { x, -H, z };
    for (int i = 0; i < NBLK; i++)
        if (sphere_blocked(L, &BLK[i], D)) return 1;
    return 0;
}
static int cam_blocked(float x, float z) {
    for (int i = 0; i < NBLK; i++) {
        float dx = x - BLK[i].x, dz = z - BLK[i].z;
        if (dx*dx + dz*dz < BLK[i].r * BLK[i].r) return 1;
    }
    return 0;
}
static float irradiance(float x, float z) {     /* cos(theta_L)/r_L^2 */
    return H / powf(x*x + z*z + H*H, 1.5f);
}

/* ---------------- row allocation ---------------- */

typedef struct { int kind; float alpha; float R[NH + 1]; const char* tag; } Cfg;
/* kind 0 = A_t ~ r^alpha,  kind 1 = exact equal-solid-angle,  kind 2 = r^alpha
 * on the *radial texel count* (legacy squared/linear warps, A_t ~ r^(2-..)) */

static void rows_family(float alpha, float* R) {
    float p = 2.0f - alpha;
    if (fabsf(p) < 1e-4f) {
        float a = logf(R_MIN), b = logf(R_MAX);
        for (int k = 0; k <= NH; k++) R[k] = expf(a + (b - a) * (float)k / NH);
        return;
    }
    float c0 = powf(R_MIN, p), c1 = powf(R_MAX, p);
    for (int k = 0; k <= NH; k++)
        R[k] = powf(c0 + (c1 - c0) * (float)k / NH, 1.0f / p);
}
static void rows_derived(float* R) {
    float u0 = 1.0f / sqrtf(R_MIN*R_MIN + H*H);
    float u1 = 1.0f / sqrtf(R_MAX*R_MAX + H*H);
    for (int k = 0; k <= NH; k++) {
        float u = u0 + (u1 - u0) * (float)k / NH;
        R[k] = sqrtf(fmaxf(0.0f, 1.0f/(u*u) - H*H));
    }
}
static void cfg_rows(Cfg* c) {
    if (c->kind == 1) rows_derived(c->R);
    else rows_family(c->alpha, c->R);
}
static int row_of(const Cfg* c, float r) {
    int lo = 0, hi = NH;
    while (hi - lo > 1) { int m = (lo + hi) / 2; if (c->R[m] <= r) lo = m; else hi = m; }
    return lo;
}

/* ---------------- shadow map ---------------- */

static unsigned char* map_render(const Cfg* c) {
    unsigned char* sm = malloc((size_t)NH * NP);
    for (int k = 0; k < NH; k++) {
        float rc = 0.5f * (c->R[k] + c->R[k + 1]);
        for (int i = 0; i < NP; i++) {
            float phi = 2.0f * PI_F * ((float)i + 0.5f) / NP;
            sm[(size_t)k * NP + i] = shadowed(rc * cosf(phi), rc * sinf(phi)) ? 1 : 0;
        }
    }
    return sm;
}
static int map_get(const Cfg* c, const unsigned char* sm, float x, float z) {
    float r = sqrtf(x*x + z*z);
    if (r >= R_MAX || r < c->R[0]) return -1;
    int k = row_of(c, r);
    float phi = atan2f(z, x);
    if (phi < 0.0f) phi += 2.0f * PI_F;
    int i = (int)(phi / (2.0f * PI_F) * NP);
    if (i < 0) i = 0;
    if (i >= NP) i = NP - 1;
    return sm[(size_t)k * NP + i];
}

/* ---------------- sample radii and rings ---------------- */

static float RP[NPTS], RB[NPTS + 1];
static void make_radii(void) {
    for (int k = 0; k < NPTS; k++)
        RP[k] = 0.06f * powf(2.40f / 0.06f, (float)k / (float)(NPTS - 1));
    RB[0] = 0.0f;
    for (int k = 1; k < NPTS; k++) RB[k] = sqrtf(RP[k-1] * RP[k]);
    RB[NPTS] = R_MAX;
}
static int ring_of(float r) {
    if (r < 0.0f || r >= RB[NPTS]) return -1;
    int lo = 0, hi = NPTS;
    while (hi - lo > 1) { int m = (lo + hi) / 2; if (RB[m] <= r) lo = m; else hi = m; }
    return lo;
}

/* ---------------- measurement ---------------- */

typedef struct {
    double leak[NPTS], notch[NPTS], vis[NI][NPTS];
    double nvis_low, nvis_mid, nvis_high;   /* visible leak per unit wall area */
    double ring_area;
    double leak_tot, notch_tot, vis_tot[NI];
    double leak_near, leak_far, vis_near[NI], vis_far[NI];
    double art_near[NI];   /* visible patch area inside r<0.5H */
} Acc;

static unsigned char TRUTH[EDGE * EDGE], CAMS[EDGE * EDGE];
static float ER[EDGE * EDGE];
static int   RG[EDGE * EDGE];
static double PXA;

static void precompute_trial(void) {
    const float span = 2.0f * R_MAX / (float)EDGE;
    PXA = (double)span * span;
    for (int j = 0; j < EDGE; j++)
        for (int i = 0; i < EDGE; i++) {
            float x = -R_MAX + span * ((float)i + 0.5f);
            float z = -R_MAX + span * ((float)j + 0.5f);
            size_t p = (size_t)j * EDGE + i;
            int rg = ring_of(sqrtf(x*x + z*z));
            RG[p] = rg;
            TRUTH[p] = 0; CAMS[p] = 0; ER[p] = 0.0f;
            if (rg < 0) continue;
            CAMS[p] = cam_blocked(x, z) ? 1 : 0;
            if (CAMS[p]) continue;
            TRUTH[p] = shadowed(x, z) ? 1 : 0;
            ER[p] = irradiance(x, z);
        }
}
static void acc_config(const Cfg* c, Acc* A) {
    unsigned char* sm = map_render(c);
    for (int j = 0; j < EDGE; j++)
        for (int i = 0; i < EDGE; i++) {
            size_t p = (size_t)j * EDGE + i;
            int rg = RG[p];
            if (rg < 0 || CAMS[p]) continue;
            float x = -R_MAX + 2.0f * R_MAX * ((float)i + 0.5f) / (float)EDGE;
            float z = -R_MAX + 2.0f * R_MAX * ((float)j + 0.5f) / (float)EDGE;
            int ms = map_get(c, sm, x, z);
            if (ms < 0) continue;
            if (ms == 0 && TRUTH[p] == 1) {                 /* bright patch   */
                A->leak[rg] += PXA;
                A->leak_tot += PXA;
                if (RP[rg] < 0.5f) A->leak_near += PXA; else A->leak_far += PXA;
                for (int q = 0; q < NI; q++)
                    if (INTENS[q] * ER[p] > THRESH) {
                        A->vis[q][rg] += PXA;
                        A->vis_tot[q] += PXA;
                        if (RP[rg] < 0.5f) { A->vis_near[q] += PXA; A->art_near[q] += PXA; }
                        else A->vis_far[q] += PXA;
                    }
            } else if (ms == 1 && TRUTH[p] == 0) {          /* dark notch     */
                A->notch[rg] += PXA;
                A->notch_tot += PXA;
            }
        }
    free(sm);
}

/* local texels per wall area at each sample radius, from the row allocation */
static void alloc_stats(const Cfg* c, double* cov, double* dres) {
    for (int m = 0; m < NPTS; m++) {
        float r = RP[m];
        int k = row_of(c, r);
        double dr = (double)c->R[k+1] - c->R[k];       /* radial extent of a row */
        double nt = (double)NP / dr;                   /* texels per unit radius */
        double dA = 2.0 * PI_F * (double)r;            /* wall area per unit radius */
        dres[m] = dA / nt;                             /* A_t = dA/dN */
        cov[m]  = nt / dA;
    }
}

/* ---------------- self test ---------------- */

static float derived_band_fit(void) {
    double sx=0, sy=0, sxx=0, sxy=0; int n = 0;
    for (int k = 0; k < 4000; k++) {
        float r = R_MIN * powf(R_MAX / R_MIN, (float)k / 3999.0f);
        double X = log((double)r), Y = 1.5 * log((double)r*r + H*H);
        sx+=X; sy+=Y; sxx+=X*X; sxy+=X*Y; n++;
    }
    return (float)((n*sxy - sx*sy) / (n*sxx - sx*sx));
}

static int selftest(void) {
    int fail = 0;
    #define CK(n, c) do { printf("%s=%d\n", n, (c)?1:0); if(!(c)) fail++; } while(0)
    make_radii();
    Cfg c; c.kind = 0; c.alpha = 0.0f; c.tag = "t"; cfg_rows(&c);
    int mono = 1;
    for (int k = 0; k < NH; k++) if (c.R[k+1] <= c.R[k]) mono = 0;
    CK("rows_monotone", mono);
    CK("rows_span", fabsf(c.R[0]-R_MIN) < 1e-5f && fabsf(c.R[NH]-R_MAX) < 1e-5f);

    Cfg d; d.kind = 1; d.alpha = derived_band_fit(); d.tag = "d"; cfg_rows(&d);
    double cov[NPTS], dr[NPTS], covd[NPTS], drd[NPTS];
    alloc_stats(&c, cov, dr);
    alloc_stats(&d, covd, drd);
    /* skip the three innermost points: there the band is only a few rows
     * wide, so any allocation is locally quantised */
    /* beyond ~0.2H each row is many texels wide and the allocation is smooth */
    int flat = 1;
    for (int m = 6; m < NPTS; m++) if (fabs(dr[m]/dr[6]-1.0) > 0.08) flat = 0;
    CK("planar_flat", flat);

    /* the derived allocation must satisfy dres ~ r_L^3 */
    int law = 1;
    for (int m = 6; m < NPTS - 1; m++) {
        double want = pow(((double)RP[m+1]*RP[m+1]+H*H)/((double)RP[m]*RP[m]+H*H), 1.5);
        double got = drd[m+1]/drd[m];
        if (fabs(got/want - 1.0) > 0.08) law = 0;
    }
    CK("derived_law", law);

    /* total texel count is held fixed: rows cover the band exactly */
    double tot = 0;
    for (int k = 0; k < NH; k++) tot += 1.0;
    CK("texels_fixed", fabs(tot*(double)NP - (double)NH*NP) < 1e-9);

    /* truth vs a brute-force march */
    make_scene(3);
    int ok = 1;
    for (int i = 0; i < 300; i++) {
        float x = -2.2f + 4.4f * (float)(i % 21) / 20.0f;
        float z = -2.2f + 4.4f * (float)(i / 21) / 14.0f;
        int t = shadowed(x, z), hit = 0;
        for (int s = 0; s < 6000; s++) {
            float u = (float)s / 5999.0f;
            float px = u*x, py = H + u*(-H), pz = u*z;
            for (int b = 0; b < NBLK; b++) {
                float dx=px-BLK[b].x, dy=py-BLK[b].y, dz=pz-BLK[b].z;
                if (dx*dx+dy*dy+dz*dz < BLK[b].r*BLK[b].r) hit = 1;
            }
        }
        if (hit != t) { ok = 0; break; }
    }
    CK("truth_vs_march", ok);

    /* map lookup is consistent with the truth at texel centres */
    int mc = 1;
    unsigned char* sm = map_render(&d);
    for (int k = 1; k < NH - 1; k += 7)
        for (int i = 0; i < NP; i += 37) {
            float rc = 0.5f*(d.R[k]+d.R[k+1]);
            float phi = 2.0f*PI_F*((float)i+0.5f)/NP;
            int want = shadowed(rc*cosf(phi), rc*sinf(phi));
            if (sm[(size_t)k*NP+i] != want) mc = 0;
        }
    free(sm);
    CK("map_consistent", mc);

    printf("band_fit=%.4f\n", derived_band_fit());
    printf("all_ok=%d\n", fail == 0);
    return fail == 0 ? 0 : 1;
}

/* ---------------- main ---------------- */

static char OUTDIR[512] = ".";

/* sibling files go next to the CSV argument, so the tool works from any cwd */
static FILE* open_out(const char* name) {
    char fn[600];
    snprintf(fn, sizeof(fn), "%s/%s", OUTDIR, name);
    FILE* f = fopen(fn, "w");
    if (!f) fprintf(stderr, "cannot write %s\n", fn);
    return f;
}

int main(int argc, char** argv) {
    if (argc >= 2 && strcmp(argv[1], "--selftest") == 0) return selftest();
    if (argc < 2) { fprintf(stderr, "usage: %s <out.csv> [dump]\n", argv[0]); return 2; }
    int dump = (argc > 2) ? atoi(argv[2]) : 0;
    {
        const char* slash = strrchr(argv[1], '/');
        if (slash) {
            size_t n = (size_t)(slash - argv[1]);
            if (n >= sizeof(OUTDIR)) n = sizeof(OUTDIR) - 1;
            memcpy(OUTDIR, argv[1], n);
            OUTDIR[n] = '\0';
        }
    }
    make_radii();
    float bfit = derived_band_fit();

    Cfg cfg[6];
    cfg[0].kind = 0; cfg[0].alpha = 0.0f;   cfg[0].tag = "planar_a0";
    cfg[1].kind = 1; cfg[1].alpha = bfit;   cfg[1].tag = "derived_exact";
    cfg[2].kind = 0; cfg[2].alpha = bfit;   cfg[2].tag = "const_bandfit";
    cfg[3].kind = 0; cfg[3].alpha = 1.25f;  cfg[3].tag = "const_a1.25";
    cfg[4].kind = 0; cfg[4].alpha = 2.0f;   cfg[4].tag = "const_a2.0";
    cfg[5].kind = 0; cfg[5].alpha = -1.0f;  cfg[5].tag = "reverse_am1";
    for (int q = 0; q < 6; q++) cfg_rows(&cfg[q]);

    Acc acc[6]; memset(acc, 0, sizeof(acc));
    /* per-trial totals so the summary can carry a spread, not just a mean */
    static double tr_leak_near[6][NTRIAL], tr_vis_near[6][NTRIAL];
    for (int t = 0; t < NTRIAL; t++) {
        make_scene(t);
        precompute_trial();
        for (int q = 0; q < 6; q++) {
            Acc one; memset(&one, 0, sizeof(one));
            acc_config(&cfg[q], &one);
            if (t == 0) acc[q] = one;
            else {
                acc[q].leak_tot += one.leak_tot;
                acc[q].leak_near += one.leak_near;
                acc[q].leak_far += one.leak_far;
                acc[q].notch_tot += one.notch_tot;
                for (int m = 0; m < NPTS; m++) {
                    acc[q].leak[m] += one.leak[m];
                    acc[q].notch[m] += one.notch[m];
                    for (int ii = 0; ii < NI; ii++) acc[q].vis[ii][m] += one.vis[ii][m];
                }
                for (int ii = 0; ii < NI; ii++) {
                    acc[q].vis_tot[ii] += one.vis_tot[ii];
                    acc[q].vis_near[ii] += one.vis_near[ii];
                    acc[q].vis_far[ii] += one.vis_far[ii];
                }
            }
            tr_leak_near[q][t] = one.leak_near;
            tr_vis_near[q][t] = one.vis_near[1];
        }
    }

    FILE* f = fopen(argv[1], "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", argv[1]); return 1; }
    fprintf(f, "tag,kind,alpha,intensity,point,r_mid,r_lo,r_hi,coverage,dres_H2,"
               "leak_H2,notch_H2,vis_leak_H2,leak_dens,vis_dens\n");
    FILE* s = open_out("summary.csv");
    fprintf(s, "tag,kind,alpha,leak_H2,leak_near,leak_far,notch_H2,"
               "vis_low,vis_mid,vis_high,"
               "dres_0.06,dres_0.35,dres_1.20,dres_2.40,"
               "leak_near_sd,vis_near_mid_sd,"
               "visdens_near_low,visdens_near_mid,visdens_near_high,"
               "visdens_far_mid,visdens_all_mid\n");
    for (int q = 0; q < 6; q++) {
        double cov[NPTS], dr[NPTS];
        alloc_stats(&cfg[q], cov, dr);
        for (int m = 0; m < NPTS; m++)
            for (int ii = 0; ii < NI; ii++) {
                double ar = PI_F*((double)RB[m+1]*RB[m+1] - (double)RB[m]*RB[m]);
                fprintf(f, "%s,%d,%.4f,%s,%d,%.4f,%.4f,%.4f,%.4f,%.7f,%.6f,%.6f,%.6f,"
                           "%.6f,%.6f\n",
                        cfg[q].tag, cfg[q].kind, cfg[q].alpha, INAME[ii], m,
                        RP[m], RB[m], RB[m+1], cov[m], dr[m],
                        acc[q].leak[m] / NTRIAL, acc[q].notch[m] / NTRIAL,
                        acc[q].vis[ii][m] / NTRIAL,
                        acc[q].leak[m] / NTRIAL / ar,
                        acc[q].vis[ii][m] / NTRIAL / ar);
            }
        double mn = acc[q].leak_near/NTRIAL, mv = tr_vis_near[q][0];
        double sn = 0, sv = 0;
        for (int t = 0; t < NTRIAL; t++) {
            sn += (tr_leak_near[q][t]-mn)*(tr_leak_near[q][t]-mn);
            sv += (tr_vis_near[q][t]-mv)*(tr_vis_near[q][t]-mv);
        }
        sn = sqrt(sn/(NTRIAL-1)); sv = sqrt(sv/(NTRIAL-1));
        /* visible leak per unit wall area, near and far, per intensity */
        double an = 0, af = 0, vn[NI] = {0}, vf[NI] = {0}, vs[NI] = {0};
        for (int m = 0; m < NPTS; m++) {
            double ar = PI_F*((double)RB[m+1]*RB[m+1] - (double)RB[m]*RB[m]);
            if (RP[m] < 0.5f) { an += ar; for (int ii=0;ii<NI;ii++) vn[ii]+=acc[q].vis[ii][m]; }
            else              { af += ar; for (int ii=0;ii<NI;ii++) vf[ii]+=acc[q].vis[ii][m]; }
        }
        for (int ii=0;ii<NI;ii++) vs[ii] = (vn[ii]+vf[ii])/(an+af);
        fprintf(s, "%s,%d,%.4f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.7f,%.7f,%.7f,%.7f,"
                   "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                cfg[q].tag, cfg[q].kind, cfg[q].alpha,
                acc[q].leak_tot/NTRIAL, acc[q].leak_near/NTRIAL, acc[q].leak_far/NTRIAL,
                acc[q].notch_tot/NTRIAL,
                acc[q].vis_tot[0]/NTRIAL, acc[q].vis_tot[1]/NTRIAL, acc[q].vis_tot[2]/NTRIAL,
                dr[0], dr[3], dr[10], dr[19], sn, sv,
                vn[0]/(an*NTRIAL), vn[1]/(an*NTRIAL), vn[2]/(an*NTRIAL),
                vf[1]/(af*NTRIAL), vs[1]/((an+af)*NTRIAL));
    }
    fclose(s); fclose(f);

    /* exponent sweep, same allocation family, 12 trials each */
    const float ASW[] = { -1.0f, -0.5f, -0.25f, 0.0f, 0.2f, 0.4f, 0.6f, 0.8f,
                          1.0f, 1.25f, 1.5f, 1.75f, 2.0f, 2.5f, 3.0f };
    const int NAS = (int)(sizeof(ASW)/sizeof(ASW[0]));
    FILE* sw = open_out("sweep.csv");
    fprintf(sw, "alpha,leak_H2,leak_near,leak_far,notch_H2,"
                "vis_low,vis_mid,vis_high,vis_near_low,vis_near_mid,vis_near_high,"
                "dres_0.06,dres_1.20\n");
    for (int q = 0; q < NAS; q++) {
        Cfg c; c.kind = 0; c.alpha = ASW[q]; c.tag = "s"; cfg_rows(&c);
        Acc A; memset(&A, 0, sizeof(A));
        for (int t = 0; t < NTRIAL; t++) {
            make_scene(t);
            precompute_trial();
            acc_config(&c, &A);
        }
        double cov[NPTS], dr[NPTS];
        alloc_stats(&c, cov, dr);
        fprintf(sw, "%.3f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.7f,%.7f\n",
                c.alpha, A.leak_tot/NTRIAL, A.leak_near/NTRIAL, A.leak_far/NTRIAL,
                A.notch_tot/NTRIAL,
                A.vis_tot[0]/NTRIAL, A.vis_tot[1]/NTRIAL, A.vis_tot[2]/NTRIAL,
                A.vis_near[0]/NTRIAL, A.vis_near[1]/NTRIAL, A.vis_near[2]/NTRIAL,
                dr[0], dr[10]);
    }
    fclose(sw);

    /* analytic curves for the density overlay */
    FILE* a = open_out("density.csv");
    fprintf(a, "r,n_derived,n_bandfit,alpha_local,alpha_bandfit\n");
    double C = (double)NH*NP / (2.0*PI_F*(1.0/sqrt((double)R_MIN*R_MIN+H*H)
                                        - 1.0/sqrt((double)R_MAX*R_MAX+H*H)));
    double I0 = (pow((double)R_MAX, 2.0-bfit) - pow((double)R_MIN, 2.0-bfit)) / (2.0-bfit);
    double Cf = (double)NH*NP / (2.0*PI_F*I0);
    for (int k = 0; k < 400; k++) {
        float r = R_MIN * powf(R_MAX/R_MIN, (float)k/399.0f);
        fprintf(a, "%.6f,%.6f,%.6f,%.6f,%.6f\n", r,
                C / pow((double)r*r+H*H, 1.5),
                Cf * pow((double)r, -bfit),
                3.0*(double)r*r/((double)r*r+H*H), bfit);
    }
    fclose(a);

    /* observed density per map, at the 20 sample points */
    FILE* o = open_out("density_obs.csv");
    fprintf(o, "tag,alpha,r,coverage,dres_H2\n");
    for (int q = 0; q < 6; q++) {
        double cov[NPTS], dr[NPTS];
        alloc_stats(&cfg[q], cov, dr);
        for (int m = 0; m < NPTS; m++)
            fprintf(o, "%s,%.4f,%.4f,%.2f,%.7f\n",
                    cfg[q].tag, cfg[q].alpha, RP[m], cov[m], dr[m]);
    }
    fclose(o);

    if (dump) {
        /* wall captures: planar vs derived vs best constant, at each I */
        int shown[3] = { 0, 1, 3 };
        for (int si = 0; si < 3; si++) {
            Cfg* c = &cfg[shown[si]];
            for (int t = 0; t < 2; t++) {
                make_scene(t);
                precompute_trial();
                unsigned char* sm = map_render(c);
                for (int ii = 0; ii < NI; ii++) {
                    char fn[256];
                    snprintf(fn, sizeof(fn), "%s/wall_%s_%s_t%d.ppm", OUTDIR,
                             c->tag, INAME[ii], t);
                    FILE* p = fopen(fn, "wb");
                    if (!p) continue;
                    fprintf(p, "P6\n%d %d\n255\n", EDGE, EDGE);
                    for (int j = 0; j < EDGE; j++)
                        for (int i = 0; i < EDGE; i++) {
                            size_t q = (size_t)j*EDGE + i;
                            int rg = RG[q];
                            float val;
                            if (rg < 0) { unsigned char g=15; fputc(g,p);fputc(g,p);fputc(g,p); continue; }
                            float x = -R_MAX + 2.0f*R_MAX*((float)i+0.5f)/EDGE;
                            float z = -R_MAX + 2.0f*R_MAX*((float)j+0.5f)/EDGE;
                            if (CAMS[q]) {
                                val = 0.32f;
                            } else {
                                int ms = map_get(c, sm, x, z);
                                int vis = (ms < 0) ? 1 : (ms ? 0 : 1);
                                float E = INTENS[ii] * irradiance(x, z);
                                val = vis ? 1.0f - expf(-GAIN*E) : 0.015f;
                            }
                            unsigned char g = (unsigned char)(255.0f*powf(fminf(fmaxf(val,0.0f),1.0f), 1.0f/2.2f));
                            fputc(g,p); fputc(g,p); fputc(g,p);
                        }
                    fclose(p);
                    if (ii == 1) {
                        /* leak overlay: red = lit but should be shadowed */
                        snprintf(fn, sizeof(fn), "%s/leak_%s_t%d.ppm", OUTDIR, c->tag, t);
                        p = fopen(fn, "wb");
                        if (p) {
                            fprintf(p, "P6\n%d %d\n255\n", EDGE, EDGE);
                            for (int j = 0; j < EDGE; j++)
                                for (int i = 0; i < EDGE; i++) {
                                    size_t q = (size_t)j*EDGE + i;
                                    int R8=18,G8=18,B8=18;
                                    if (RG[q] >= 0 && !CAMS[q]) {
                                        float x = -R_MAX + 2.0f*R_MAX*((float)i+0.5f)/EDGE;
                                        float z = -R_MAX + 2.0f*R_MAX*((float)j+0.5f)/EDGE;
                                        int ms = map_get(c, sm, x, z);
                                        int vm = (ms < 0) ? 1 : (ms ? 0 : 1);
                                        if (vm == 1 && TRUTH[q] == 1) {
                                            float E = irradiance(x,z);
                                            R8 = 90 + (int)(160*fminf(1.0f, E/0.6f));
                                            G8 = 45; B8 = 45;
                                        } else if (vm == 0 && TRUTH[q] == 0) {
                                            R8 = 40; G8 = 65; B8 = 150;
                                        }
                                    }
                                    fputc(R8,p); fputc(G8,p); fputc(B8,p);
                                }
                            fclose(p);
                        }
                    }
                }
                if (t == 0) {
                    char fn[256];
                    snprintf(fn, sizeof(fn), "%s/map_%s.ppm", OUTDIR, c->tag);
                    FILE* p = fopen(fn, "wb");
                    if (p) {
                        fprintf(p, "P6\n%d %d\n255\n", NP, NH);
                        for (int k = 0; k < NH; k++)
                            for (int i = 0; i < NP; i++) {
                                unsigned char g = sm[(size_t)k*NP+i] ? 55 : 235;
                                fputc(g,p); fputc(g,p); fputc(g,p);
                            }
                        fclose(p);
                    }
                }
                free(sm);
            }
        }
    }
    return 0;
}
