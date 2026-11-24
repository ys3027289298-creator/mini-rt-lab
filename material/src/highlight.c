/* highlight: grazing-angle GGX highlight experiment.

   Scene: an infinite GGX plane (z=0) under an analytic environment light
   (gaussian "sun" near the horizon + faint ambient), so the highlight sits
   at grazing incidence near the horizon line.

   Estimators compared:
     old      cosine-weighted hemisphere sampling (the repo's original)
     vndf     visible-normal sampling (Heitz 2018)
     vndf+tb  VNDF plus a Taylor-modelled grazing band: inside
              t = cos(theta_i) < band the MC contribution is replaced by
              the deterministic quadrature of the second-order Taylor
              expansion of q(t) = pdf*brdf. This is a biased estimator —
              the bias is exactly the truncation error of the expansion
              integrated over the band, and is reported separately; the
              variance drops because in-band sample noise is removed
              entirely. The band width is chosen per pixel as the largest
              band whose measured truncation error stays under 10%.

   Outputs (to argv[1], default ../artifacts):
     highlight_old.pfm / highlight_new.pfm   full frames (linear, grayscale)
     stats.json                              crop window, per-roughness
                                             variances, truncation errors */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ggx.h"

#define W 256
#define H 144
#define SPP 64
#define REPS 16
#define CROP_W 160
#define CROP_H 56
#define T0_MAX 0.1 /* widest grazing band considered, in t = cos(theta_i) */
#define F0 0.04

/* ---------------------------------------------------------------- rng */
static unsigned long long rng_state;
static void rng_seed(unsigned long long s) {
    rng_state = s ? s : 0x9e3779b97f4a7c15ull;
}
static double rnd(void) { /* splitmix64 */
    unsigned long long z = (rng_state += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    z = z ^ (z >> 31);
    return (double)(z >> 11) * (1.0 / 9007199254740992.0);
}

/* ---------------------------------------------------------------- env */
typedef struct { V3d ldir; double sigma, L0, amb; } Env;

static double env_Li(const Env *e, V3d w) {
    if (w.z <= 0.0) return 0.0;
    double c = ggx_dot(w, e->ldir);
    if (c > 1.0) c = 1.0;
    double ang = acos(c);
    return e->amb + e->L0 * exp(-ang * ang / (2.0 * e->sigma * e->sigma));
}

/* --------------------------------------------------- camera / geometry */
static V3d cam_pos, cam_fwd, cam_right, cam_up;
static const double TAN_HALF = 0.36;

static void camera_init(void) {
    cam_pos = (V3d){ 0.0, -3.0, 1.0 };
    cam_fwd = ggx_norm((V3d){ 0.0, 1.0, -0.10 });
    V3d zup = { 0.0, 0.0, 1.0 };
    V3d r = { cam_fwd.y * zup.z - cam_fwd.z * zup.y,
              cam_fwd.z * zup.x - cam_fwd.x * zup.z,
              cam_fwd.x * zup.y - cam_fwd.y * zup.x };
    cam_right = ggx_norm(r);
    cam_up = (V3d){ cam_right.y * cam_fwd.z - cam_right.z * cam_fwd.y,
                    cam_right.z * cam_fwd.x - cam_right.x * cam_fwd.z,
                    cam_right.x * cam_fwd.y - cam_right.y * cam_fwd.x };
}

static V3d pixel_ray(int px, int py) {
    double nx = ((px + 0.5) / W * 2.0 - 1.0) * TAN_HALF;
    double ny = (1.0 - (py + 0.5) / H * 2.0) * TAN_HALF * ((double)H / W);
    V3d d = { cam_fwd.x + nx * cam_right.x + ny * cam_up.x,
              cam_fwd.y + nx * cam_right.y + ny * cam_up.y,
              cam_fwd.z + nx * cam_right.z + ny * cam_up.z };
    return ggx_norm(d);
}

/* ------------------------------------------------- grazing-band model */

/* largest band [0, t] whose Taylor-model truncation error stays < 10%,
   probed on a coarse (phi, t) grid; 0 disables the band replacement */
static double select_band(V3d wo, double alpha) {
    static const double cand[] = { 0.1, 0.05, 0.02, 0.01, 0.005, 0.002 };
    for (unsigned ci = 0; ci < sizeof(cand) / sizeof(cand[0]); ci++) {
        double tb = cand[ci];
        int ok = 1;
        for (int j = 0; j < 8 && ok; j++) {
            double phi = (j + 0.5) * 2.0 * GGX_PI / 8.0;
            double c0, c1, c2;
            q_taylor(wo, phi, alpha, F0, &c0, &c1, &c2);
            for (int k = 1; k <= 3; k++) {
                double t = tb * k / 3.0;
                double qe = q_product(wo, t, phi, alpha, F0);
                double qt = q_taylor_eval(c0, c1, c2, t);
                double denom = qe > 1e-30 ? qe : 1e-30;
                if (fabs(qe - qt) / denom > 0.10) { ok = 0; break; }
            }
        }
        if (ok) return tb;
    }
    return 0.0;
}

/* modelled band contribution: ∫_band q_tilde Li t / p dω, dω = dt dphi */
static double quad_band(const Env *e, V3d wo, double alpha, double band) {
    const int NP = 48, NT = 12;
    double sum = 0.0;
    for (int j = 0; j < NP; j++) {
        double phi = (j + 0.5) * 2.0 * GGX_PI / NP;
        double c0, c1, c2;
        q_taylor(wo, phi, alpha, F0, &c0, &c1, &c2);
        for (int k = 0; k < NT; k++) {
            double t = (k + 0.5) * band / NT;
            double st = sqrt(1.0 - t * t);
            V3d wi = { st * cos(phi), st * sin(phi), t };
            double p = vndf_pdf(wo, wi, alpha);
            sum += q_taylor_eval(c0, c1, c2, t) * env_Li(e, wi) * t / p;
        }
    }
    return sum * (2.0 * GGX_PI / NP) * (band / NT);
}

/* ------------------------------------------------------------ estimators */
static double mirror_L(const Env *e, V3d wo) {
    V3d wr = { -wo.x, -wo.y, wo.z };
    return ggx_schlick(wo.z, F0) * env_Li(e, wr);
}

static double est_cosine(const Env *e, V3d wo, double alpha) {
    if (alpha < GGX_ALPHA_MIN) return mirror_L(e, wo);
    double acc = 0.0;
    for (int s = 0; s < SPP; s++) {
        double u1 = rnd(), u2 = rnd();
        double r = sqrt(u1), phi = 2.0 * GGX_PI * u2;
        V3d wi = { r * cos(phi), r * sin(phi), sqrt(1.0 - u1 > 0 ? 1.0 - u1 : 0) };
        acc += GGX_PI * brdf_ggx(wo, wi, alpha, F0) * env_Li(e, wi);
    }
    return acc / SPP;
}

static double est_vndf(const Env *e, V3d wo, double alpha, double band,
                       int replace, double *band_part) {
    if (alpha < GGX_ALPHA_MIN) return mirror_L(e, wo);
    double acc = 0.0, bacc = 0.0;
    double G1o = ggx_G1(wo.z, alpha);
    for (int s = 0; s < SPP; s++) {
        V3d wi = vndf_sample(wo, alpha, rnd(), rnd());
        int inband = wi.z >= 0.0 && wi.z < band;
        if (inband && replace)
            continue; /* in-band contribution is replaced by quadrature */
        V3d h = ggx_norm(ggx_add(wo, wi));
        double w = ggx_schlick(ggx_dot(wi, h), F0)
                 * ggx_G2(wi.z, wo.z, alpha) / G1o;
        double x = w * env_Li(e, wi);
        if (inband) bacc += x; else acc += x;
    }
    if (band_part) *band_part = bacc / SPP;
    double L = (acc + (replace ? 0.0 : bacc)) / SPP;
    if (replace && band > 0.0) L += quad_band(e, wo, alpha, band);
    return L;
}

/* ------------------------------------------------------------- output */
static void write_pfm(const char *path, const double *img) {
    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); exit(1); }
    fprintf(f, "PF\n%d %d\n-1.0\n", W, H);
    for (int y = H - 1; y >= 0; y--)
        for (int x = 0; x < W; x++) {
            float v = (float)img[y * W + x];
            fwrite(&v, 4, 1, f); fwrite(&v, 4, 1, f); fwrite(&v, 4, 1, f);
        }
    fclose(f);
}

/* reference radiance by dense quadrature (for relative bias reporting) */
static double ref_L(const Env *e, V3d wo, double alpha) {
    if (alpha < GGX_ALPHA_MIN) return mirror_L(e, wo);
    const int NP = 256, NT = 96;
    double sum = 0.0;
    for (int j = 0; j < NP; j++)
        for (int k = 0; k < NT; k++) {
            double phi = (j + 0.5) * 2.0 * GGX_PI / NP;
            double t = (k + 0.5) / NT;
            double st = sqrt(1.0 - t * t);
            V3d wi = { st * cos(phi), st * sin(phi), t };
            sum += brdf_ggx(wo, wi, alpha, F0) * env_Li(e, wi) * t;
        }
    return sum * (2.0 * GGX_PI / NP) / NT;
}

/* truncation error of the Taylor model inside the chosen band, plus the
   bias it would induce if q_tilde were used without the E[g] correction */
static void truncation(const Env *e, V3d wo, double alpha, double band,
                       double *abs_max, double *rel_max, double *bias) {
    *abs_max = *rel_max = *bias = 0.0;
    if (band <= 0.0) return;
    const int NP = 64, NT = 24;
    for (int j = 0; j < NP; j++) {
        double phi = (j + 0.5) * 2.0 * GGX_PI / NP;
        double c0, c1, c2;
        q_taylor(wo, phi, alpha, F0, &c0, &c1, &c2);
        for (int k = 0; k < NT; k++) {
            double t = (k + 0.5) * band / NT;
            double st = sqrt(1.0 - t * t);
            V3d wi = { st * cos(phi), st * sin(phi), t };
            double qe = q_product(wo, t, phi, alpha, F0);
            double qt = q_taylor_eval(c0, c1, c2, t);
            double d = fabs(qe - qt);
            if (d > *abs_max) *abs_max = d;
            double rel = d / (qe > 1e-30 ? qe : 1e-30);
            if (rel > *rel_max) *rel_max = rel;
            double p = vndf_pdf(wo, wi, alpha);
            *bias += (qe - qt) * env_Li(e, wi) * t / p;
        }
    }
    *bias = fabs(*bias) * (2.0 * GGX_PI / NP) * (band / NT);
}

int main(int argc, char **argv) {
    const char *outdir = argc > 1 ? argv[1] : "../artifacts";
    const double roughs[] = { 0.3, 0.1, 0.02, 0.0 };
    const int NG = 4;
    Env env;
    env.ldir = ggx_norm((V3d){ 0.0, 0.995, 0.10 });
    env.sigma = 0.05;
    env.L0 = 300.0;
    env.amb = 0.03;
    camera_init();

    /* locate the highlight: coarse reference pass at roughness 0.1 */
    int hx = W / 2, hy = H / 2;
    double best = -1.0;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            V3d d = pixel_ray(x, y);
            if (d.z >= 0.0) continue;
            V3d wo = { -d.x, -d.y, -d.z };
            rng_seed(1234);
            double L = est_vndf(&env, wo, 0.01, 0.0, 0, NULL);
            if (L > best) { best = L; hx = x; hy = y; }
        }
    int cx = hx - CROP_W / 2, cy = hy - CROP_H / 2;
    if (cx < 0) cx = 0;
    if (cx + CROP_W > W) cx = W - CROP_W;
    if (cy < 0) cy = 0;
    if (cy + CROP_H > H) cy = H - CROP_H;
    fprintf(stderr, "highlight at (%d,%d), crop (%d,%d)+(%d,%d)\n",
            hx, hy, cx, cy, CROP_W, CROP_H);

    /* full frames at roughness 0.1: old vs new */
    double *img_old = calloc(W * H, sizeof(double));
    double *img_new = calloc(W * H, sizeof(double));
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            V3d d = pixel_ray(x, y);
            if (d.z >= 0.0) {
                img_old[y * W + x] = img_new[y * W + x] = env_Li(&env, d);
                continue;
            }
            V3d wo = { -d.x, -d.y, -d.z };
            double band = select_band(wo, 0.01);
            rng_seed(1000003ull * (y * W + x) + 7);
            img_old[y * W + x] = est_cosine(&env, wo, 0.01);
            rng_seed(1000003ull * (y * W + x) + 7);
            img_new[y * W + x] = est_vndf(&env, wo, 0.01, band, 1, NULL);
        }

    char path[512];
    snprintf(path, sizeof path, "%s/highlight_old.pfm", outdir);
    write_pfm(path, img_old);
    snprintf(path, sizeof path, "%s/highlight_new.pfm", outdir);
    write_pfm(path, img_new);

    /* per-group variance over the crop + truncation error at crop centre */
    FILE *js = NULL;
    snprintf(path, sizeof path, "%s/stats.json", outdir);
    js = fopen(path, "w");
    if (!js) { perror(path); return 1; }
    fprintf(js, "{\n \"crop\": {\"x\": %d, \"y\": %d, \"w\": %d, \"h\": %d},\n",
            cx, cy, CROP_W, CROP_H);
    fprintf(js, " \"spp\": %d, \"reps\": %d, \"t0_max\": %g,\n", SPP, REPS, T0_MAX);
    fprintf(js, " \"groups\": [\n");

    V3d dc = pixel_ray(hx, hy);
    V3d woc = { -dc.x, -dc.y, -dc.z };

    for (int g = 0; g < NG; g++) {
        double rough = roughs[g];
        double alpha = rough * rough;
        double sv[3] = { 0, 0, 0 }; /* summed per-pixel variances */
        double sband = 0.0;         /* summed in-band variance (removed) */
        int npix = 0;
        for (int y = cy; y < cy + CROP_H; y++)
            for (int x = cx; x < cx + CROP_W; x++) {
                V3d d = pixel_ray(x, y);
                if (d.z >= 0.0) continue;
                V3d wo = { -d.x, -d.y, -d.z };
                double band = alpha < GGX_ALPHA_MIN ? 0.0 : select_band(wo, alpha);
                double est[3][REPS], bp[REPS];
                for (int r = 0; r < REPS; r++) {
                    rng_seed(977ull * (y * W + x) + 31ull * r + 5);
                    est[0][r] = est_cosine(&env, wo, alpha);
                    rng_seed(977ull * (y * W + x) + 31ull * r + 5);
                    est[1][r] = est_vndf(&env, wo, alpha, band, 0, &bp[r]);
                    rng_seed(977ull * (y * W + x) + 31ull * r + 5);
                    est[2][r] = est_vndf(&env, wo, alpha, band, 1, NULL);
                }
                for (int k = 0; k < 3; k++) {
                    double m = 0.0;
                    for (int r = 0; r < REPS; r++) m += est[k][r];
                    m /= REPS;
                    double v = 0.0;
                    for (int r = 0; r < REPS; r++) v += (est[k][r] - m) * (est[k][r] - m);
                    sv[k] += v / (REPS - 1);
                }
                /* variance of the grazing-band contribution itself: this
                   is what the Taylor replacement removes entirely */
                double m = 0.0;
                for (int r = 0; r < REPS; r++) m += bp[r];
                m /= REPS;
                double v = 0.0;
                for (int r = 0; r < REPS; r++) v += (bp[r] - m) * (bp[r] - m);
                sband += v / (REPS - 1);
                npix++;
            }
        double vc = sv[0] / npix, vv = sv[1] / npix, vcv = sv[2] / npix;
        double vband = sband / npix;

        double band_c = alpha < GGX_ALPHA_MIN ? 0.0 : select_band(woc, alpha);
        double tabs = 0, trel = 0, tbias = 0, lref = 0;
        if (alpha >= GGX_ALPHA_MIN) {
            truncation(&env, woc, alpha, band_c, &tabs, &trel, &tbias);
            lref = ref_L(&env, woc, alpha);
        } else {
            lref = mirror_L(&env, woc);
        }

        fprintf(js,
                "  {\"rough\": %g, \"band\": %g, \"var_cos\": %.6e, "
                "\"var_vndf\": %.6e, \"var_vndf_tb\": %.6e, "
                "\"var_band\": %.6e, "
                "\"taylor_abs\": %.6e, \"taylor_rel\": %.6e, "
                "\"taylor_bias\": %.6e, \"L_ref\": %.6e}%s\n",
                rough, band_c, vc, vv, vcv, vband, tabs, trel, tbias, lref,
                g == NG - 1 ? "" : ",");
        fprintf(stderr,
                "rough=%.2f band=%.4f var cos=%.3e vndf=%.3e vndf+tb=%.3e "
                "(x%.1f) band_var=%.3e trunc rel=%.2e bias=%.2e\n",
                rough, band_c, vc, vv, vcv, vc / (vcv > 0 ? vcv : 1e-30),
                vband, trel, tbias);
    }
    fprintf(js, " ]\n}\n");
    fclose(js);
    free(img_old);
    free(img_new);
    return 0;
}
