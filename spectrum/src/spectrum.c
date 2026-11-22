
/* spectrum: wavelength-importance-sampled spectral rendering into CIE XYZ and
   sRGB, with a verified error bound against a dense reference.

   Something that "has no rainbow" almost always samples the spectrum at a few
   fixed wavelengths. Those samples miss narrow bands (an emission line hits
   the floor of a 5 nm grid and the colour reads wrong by tens of deltaE), and
   they never resolve a grating's angular dispersion. This module instead
   estimates the whole spectral integral

     X = int P(lambda) xbar(lambda) dlambda      CIE 1931 2-deg observer
     Y = int P(lambda) ybar(lambda) dlambda      (360..830 nm, 1 nm table)
     Z = int P(lambda) zbar(lambda) dlambda

   with wavelength importance sampling:

     P(lambda) is laid on K uniform bins, the target population is
       p_k ~ 0.98 * w_k/sum(w) + 0.02/K,   w_k = |P| * (xbar+ybar+zbar) at
     the bin centre, and samples inside a bin are uniform. The 2% floor keeps
     a bin the target missed from ever carrying unbounded weight. This is the
     three-channel generalisation of the optimal P*cmf density; the same
     machinery gives blue through cyan, green, yellow, orange and red for a
     grating because each direction is integrated over its own band.

   Then the sampled XYZ is mapped to display sRGB and compared with a dense,
   deterministic reference over the full observer range (uniform trapezoid,
   refined and checked in the self-test). The colour difference is given
   upper bounds, all enforced by the self-test:

     statistical   3*se per channel (CLT, 99.7%); se from the sample variance
     truncation    int |P|*xbar,|P|*ybar,|P|*zbar outside [lmin,lmax]: the
                   colour a truncated render box silently drops
     dE            CIE76 bound from the corner extremes of the XYZ box
     display       the XYZ bound box mapped through the RGB matrix; gamma'<=1
                   and clamping means 255*gamma'(channel)*|dRGB| levels

   The wavelength range is validated up front. A non-finite, reversed,
   outside-table or too-narrow range is rejected with the legal interval
   (360 nm <= lmin < lmax <= 830 nm) in the message; the estimator then
   refuses to produce numbers.

   Usage:  spectrum            run the self-test
           spectrum out.ppm    also write a spectral-grating strip (P6)
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "cie1931.h"

/* ------------------------------------------------------------- constants */

#define LAM_MIN    360.0
#define LAM_MAX    830.0
#define YBAR_INT   106.856914917   /* int ybar dlambda over the 1 nm table  */
#define NORM       (1.0 / YBAR_INT) /* equal-energy radiance 1 -> Y = 1      */
#define STEP_MIN   0.005          /* reference / quadrature restart step   */
#define REF_M      8              /* Simpson subintervals per 1 nm cell     */
#define MAX_BINS   512
#define PI_        3.14159265358979323846

/* ------------------------------------------------------------- colour */

/* piecewise-linear lookup of the CIE 1931 2-degree observer */
static double cmf(int ch, double lambda) {
    double t = (lambda - CIE_FIRST) / CIE_STEP;
    int k = (int)floor(t);
    if (k < 0) k = 0;
    if (k > CIE_N - 2) k = CIE_N - 2;
    double f = t - (double)k;
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;
    const double* tab = ch == 0 ? CIE_X : ch == 1 ? CIE_Y : CIE_Z;
    return tab[k] + (tab[k + 1] - tab[k]) * f;
}

/* sRGB / D65: XYZ -> linear RGB */
static const double XYZ2RGB[3][3] = {
    {  3.2404542, -1.5371385, -0.4985314 },
    { -0.9692660,  1.8760108,  0.0415560 },
    {  0.0556434, -0.2040259,  1.0572252 }
};

/* D65 relative spectral power, 380..780 nm at 5 nm (CIE 15:2004) */
static const double D65[81] = {
    49.9755, 52.3118, 54.6482, 68.7015, 82.7549, 87.1204, 91.4860, 92.4589,
    93.4318, 90.0570, 86.6823, 95.7736, 104.865, 110.936, 117.008, 117.410,
    117.812, 116.336, 114.861, 115.392, 115.923, 112.367, 108.811, 109.082,
    109.354, 108.578, 107.802, 106.296, 104.790, 106.239, 107.689, 106.047,
    104.405, 104.225, 104.046, 102.023, 100.000, 98.1671, 96.3342, 96.0611,
    95.7880, 92.2368, 88.6856, 89.3459, 90.0062, 89.8026, 89.5991, 88.6489,
    87.6987, 85.4936, 83.2886, 83.4939, 83.6992, 81.8630, 80.0268, 80.1207,
    80.2146, 81.2462, 82.2778, 80.2810, 78.2842, 74.0027, 69.7213, 70.6652,
    71.6091, 72.9790, 74.3490, 67.9765, 61.6040, 65.7448, 69.8856, 72.4863,
    75.0870, 69.3398, 63.5927, 55.0054, 46.4182, 56.6118, 66.8054, 65.0941,
    63.3828
};
static double d65(double lambda) {
    double t = (lambda - 380.0) / 5.0;
    int k = (int)floor(t);
    if (k < 0 || k > 79) return 0.0;
    double f = t - (double)k;
    return D65[k] + (D65[k + 1] - D65[k]) * f;
}

static double gamma_encode(double v) {
    if (!(v > 0.0)) return 0.0;
    if (v <= 0.0031308) return 12.92 * v;
    return 1.055 * pow(v, 1.0 / 2.4) - 0.055;
}
/* derivative of the clamped display transfer */
static double gamma_slope(double v) {
    if (v <= 0.0) return 12.92;
    if (v >= 1.0) return 0.0;
    if (v <= 0.0031308) return 12.92;
    return 1.055 / 2.4 * pow(v, 1.0 / 2.4 - 1.0);
}

static void xyz_to_linear(const double xyz[3], double lin[3]) {
    for (int i = 0; i < 3; i++)
        lin[i] = XYZ2RGB[i][0] * xyz[0] + XYZ2RGB[i][1] * xyz[1] + XYZ2RGB[i][2] * xyz[2];
}
/* CIE L*a*b* (D65 white) and the CIE76 difference */
static void xyz_to_lab(const double xyz[3], double lab[3]) {
    const double wp[3] = { 0.95047, 1.0, 1.08883 };
    double f[3];
    for (int i = 0; i < 3; i++) {
        double t = xyz[i] / wp[i];
        f[i] = t > 0.008856 ? cbrt(t) : 7.787 * t + 16.0 / 116.0;
    }
    lab[0] = 116.0 * f[1] - 16.0;
    lab[1] = 500.0 * (f[0] - f[1]);
    lab[2] = 200.0 * (f[1] - f[2]);
}
static double delta_e(const double a[3], const double b[3]) {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) +
                (a[1] - b[1]) * (a[1] - b[1]) +
                (a[2] - b[2]) * (a[2] - b[2]));
}
/* LCh hue in degrees */
static double lab_hue(const double lab[3]) {
    double h = atan2(lab[2], lab[1]) * 180.0 / PI_;
    return h < 0.0 ? h + 360.0 : h;
}

/* --------------------------------------------------------------- scenes */

typedef struct Scene {
    const char* name;
    double (*P)(const void*, double lambda);
    const void* param;
    double lmin, lmax;
} Scene;

/* Lambertian reflector of a given illuminant: P = albedo * L / pi */
typedef struct { double albedo; } Reflector;
static double flat_P(const void* p, double lambda) {
    (void)lambda;
    return ((const Reflector*)p)->albedo / PI_;
}
static double d65_P(const void* p, double lambda) {
    return ((const Reflector*)p)->albedo * d65(lambda) / PI_;
}

/* Planckian radiator, normalised to its own peak (T is the colour knob) */
static double planck_raw(double lam_nm, double T) {
    const double h = 6.62607015e-34, c = 2.99792458e8, k = 1.380649e-23;
    double l = lam_nm * 1e-9;
    return 1.0 / (pow(l, 5.0) * expm1(h * c / (l * k * T)));
}
typedef struct { double T; } Black;
static double black_P(const void* p, double lambda) {
    double T = ((const Black*)p)->T;
    double lp = 2.897771955e-3 / T * 1e9;       /* Wien peak, nm */
    return planck_raw(lambda, T) / planck_raw(lp, T);
}

/* a narrow spectral band: an emission line or one direction of a grating */
typedef struct { double center, sigma; } Band;
static double band_P(const void* p, double lambda) {
    const Band* b = p;
    double g = (lambda - b->center) / b->sigma;
    return exp(-0.5 * g * g);
}

/* ------------------------------------------------------------- validation */

typedef enum {
    RANGE_OK = 0,
    RANGE_NOT_FINITE,
    RANGE_ORDER,
    RANGE_OUTSIDE,
    RANGE_TOO_NARROW
} RangeStatus;

static RangeStatus check_range(double lmin, double lmax, char* msg, size_t n) {
    if (!isfinite(lmin) || !isfinite(lmax)) {
        snprintf(msg, n, "[%g, %g] is not finite; legal range: %.0f <= lmin < lmax <= %.0f nm",
                 lmin, lmax, LAM_MIN, LAM_MAX);
        return RANGE_NOT_FINITE;
    }
    if (!(lmin < lmax)) {
        snprintf(msg, n, "empty range [%g, %g]; legal range: %.0f <= lmin < lmax <= %.0f nm",
                 lmin, lmax, LAM_MIN, LAM_MAX);
        return RANGE_ORDER;
    }
    if (lmin < LAM_MIN || lmax > LAM_MAX) {
        snprintf(msg, n, "[%g, %g] is outside the observer table; legal range: "
                         "%.0f <= lmin < lmax <= %.0f nm", lmin, lmax, LAM_MIN, LAM_MAX);
        return RANGE_OUTSIDE;
    }
    if (lmax - lmin < 2.0 * STEP_MIN) {
        snprintf(msg, n, "[%g, %g] is narrower than the %.3f nm step; legal range: "
                         "%.0f <= lmin < lmax <= %.0f nm", lmin, lmax, 2.0 * STEP_MIN,
                 LAM_MIN, LAM_MAX);
        return RANGE_TOO_NARROW;
    }
    return RANGE_OK;
}

/* -------------------------------------------------------------- reference */

/* Deterministic reference. Between observer-table nodes (1 nm) the CMF is
   linear; between illuminant nodes (5 nm) D65 is linear; their product is
   then quadratic and Simpson's rule is exact on each 1 nm cell. Smooth
   spectra (Planck, Gaussians) converge quickly under the same subdivision.
   Each 1 nm cell is split into `m` subintervals, m even. */
static void integrate_box(const Scene* s, double lo, double hi, int m, double out[3]) {
    for (int ch = 0; ch < 3; ch++) out[ch] = 0.0;
    if (!(hi > lo) || !(s->P != NULL)) return;
    /* subdivide on the 1 nm observer grid plus the two interval ends */
    double cuts[1024];
    int nc = 0;
    cuts[nc++] = lo;
    double g = floor(lo) + 1.0;
    while (g < hi && nc < 1022) {
        if (g > lo) cuts[nc++] = g;
        g += 1.0;
    }
    cuts[nc++] = hi;
    for (int seg = 0; seg + 1 < nc; seg++) {
        int per = m;
        for (int k = 0; k < per; k += 2) {
            double a = cuts[seg] + (cuts[seg + 1] - cuts[seg]) * (double)k / (double)per;
            double b = cuts[seg] + (cuts[seg + 1] - cuts[seg]) * (double)(k + 1) / (double)per;
            double c = cuts[seg] + (cuts[seg + 1] - cuts[seg]) * (double)(k + 2) / (double)per;
            for (int ch = 0; ch < 3; ch++) {
                double fa = s->P(s->param, a) * cmf(ch, a);
                double fb = s->P(s->param, b) * cmf(ch, b);
                double fc = s->P(s->param, c) * cmf(ch, c);
                out[ch] += (fa + 4.0 * fb + fc) * ((c - a) / 6.0);
            }
        }
    }
    for (int ch = 0; ch < 3; ch++) out[ch] *= NORM;
}

/* Rigorous upper bound on the energy the render box leaves outside
   [lmin,lmax]: int |P|*CMF over the two side intervals, computed with
   Simpson at a fine fixed step and inflated for its own quadrature error. */
static double simpson_abs(const Scene* s, double lo, double hi, int ch) {
    if (!(hi > lo)) return 0.0;
    double step = 0.002;
    long n = (long)ceil((hi - lo) / step);
    if (n < 2) n = 2;
    if (n & 1) n++;
    double h = (hi - lo) / (double)n;
    double acc = 0.0;
    double f0 = fabs(s->P(s->param, lo) * cmf(ch, lo));
    for (long k = 0; k < n; k += 2) {
        double x1 = lo + (double)(k + 1) * h;
        double x2 = lo + (double)(k + 2) * h;
        double f1 = fabs(s->P(s->param, x1) * cmf(ch, x1));
        double f2 = x2 >= hi ? fabs(s->P(s->param, hi) * cmf(ch, hi))
                             : fabs(s->P(s->param, x2) * cmf(ch, x2));
        acc += (f0 + 4.0 * f1 + f2) * h / 3.0;
        f0 = f2;
    }
    return acc * NORM * (1.0 + 1e-6);
}
static void integrate_outside(const Scene* s, double out[3]) {
    for (int ch = 0; ch < 3; ch++)
        out[ch] = simpson_abs(s, LAM_MIN, s->lmin, ch) + simpson_abs(s, s->lmax, LAM_MAX, ch);
}

/* --------------------------------------------------------------- sampler */

typedef struct { unsigned long long s; } Rng;
static unsigned long long rng_next(Rng* r) {
    r->s = r->s * 6364136223846793005ULL + 1442695040888963407ULL;
    return r->s;
}
static double rng_unit(Rng* r) {
    return (double)((rng_next(r) >> 11) & ((1ULL << 53) - 1)) / 9007199254740992.0;
}

typedef struct {
    double xyz[3];     /* estimate                                  */
    double se[3];      /* standard error                            */
    double trunc[3];   /* deterministic outside-box term            */
    double unif[3];    /* uniform weights on the same draws         */
    double bound[3];   /* 3*se + trunc                              */
    long n;
    int ok;            /* range/argument validation                 */
} Estimate;

/* Wavelength importance sampling of int P*CMF over [lmin,lmax].
   K uniform bins, per-bin sampling probability from the centre values of
   |P|*(xbar+ybar+zbar), 2% floor, samples uniform inside the chosen bin,
   stratified over the K*N independent uniform variates. */
static Estimate estimate(const Scene* s, int K, long N, unsigned seed) {
    Estimate e;
    memset(&e, 0, sizeof e);
    char msg[256];
    if (s == NULL || check_range(s->lmin, s->lmax, msg, sizeof msg) != RANGE_OK ||
        K < 1 || K > MAX_BINS || N < 2) {
        return e;
    }
    e.ok = 1;

    double lo = s->lmin, hi = s->lmax;
    double bw = (hi - lo) / (double)K;

    /* per-bin target population, then the sampling probability */
    double w[MAX_BINS], prob[MAX_BINS], cum[MAX_BINS];
    double wsum = 0.0;
    for (int k = 0; k < K; k++) {
        double mid = lo + ((double)k + 0.5) * bw;
        double cm = cmf(0, mid) + cmf(1, mid) + cmf(2, mid);
        w[k] = fabs(s->P(s->param, mid)) * cm;
        if (!(w[k] >= 0.0) || !isfinite(w[k])) w[k] = 0.0;
        wsum += w[k];
    }
    if (!(wsum > 0.0)) {                    /* P == 0 everywhere */
        e.ok = 1;
        e.n = N;
        return e;
    }
    double acc = 0.0;
    for (int k = 0; k < K; k++) {
        prob[k] = 0.98 * w[k] / wsum + 0.02 / (double)K;
        acc += prob[k];
        cum[k] = acc;
    }
    cum[K - 1] = 1.0;

    /* the same wavelengths reweighted by a uniform density: the variance
       reduction the importance map buys, measured paired */
    double span = hi - lo;
    double pdfu = 1.0 / span;

    /* truncation: the colour this box throws away (deterministic) */
    integrate_outside(s, e.trunc);

    /* draw: stratified u chooses the bin, uniform position inside it */
    Rng rng = { 0x243F6A8885A308D3ULL ^ ((unsigned long long)seed << 32) ^ (unsigned)K };
    double sum[3] = { 0, 0, 0 }, sum2[3] = { 0, 0, 0 };
    double psum[3] = { 0, 0, 0 }, psum2[3] = { 0, 0, 0 };
    for (long i = 0; i < N; i++) {
        double u = ((double)i + rng_unit(&rng)) / (double)N;
        int k = 0;
        while (k < K - 1 && u > cum[k]) k++;
        double lambda = lo + ((double)k + rng_unit(&rng)) * bw;
        double pdf = prob[k] / bw;
        double p = s->P(s->param, lambda);
        for (int ch = 0; ch < 3; ch++) {
            double f = p * cmf(ch, lambda) * NORM;
            double T = f / pdf;
            sum[ch] += T;
            sum2[ch] += T * T;
            double Tu = f / pdfu;
            psum[ch] += Tu;
            psum2[ch] += Tu * Tu;
        }
    }
    e.n = N;
    for (int ch = 0; ch < 3; ch++) {
        double mean = sum[ch] / (double)N;
        double var = (sum2[ch] / (double)N - mean * mean) * (double)N / (double)(N - 1);
        if (!(var > 0.0)) var = 0.0;
        e.xyz[ch] = mean;
        e.se[ch] = sqrt(var / (double)N);
        e.bound[ch] = 3.0 * e.se[ch] + e.trunc[ch];
        double mu = psum[ch] / (double)N;
        double vu = (psum2[ch] / (double)N - mu * mu) * (double)N / (double)(N - 1);
        e.unif[ch] = sqrt((vu > 0.0 ? vu : 0.0) / (double)N);
    }
    return e;
}

/* fixed-wavelength estimator, the thing that "has no rainbow": m samples at
   the same wavelengths for every band, equal weights */
static void estimate_fixed(const Scene* s, const double* lams, int m, double out[3]) {
    for (int ch = 0; ch < 3; ch++) out[ch] = 0.0;
    for (int i = 0; i < m; i++) {
        double p = s->P(s->param, lams[i]);
        for (int ch = 0; ch < 3; ch++)
            out[ch] += p * cmf(ch, lams[i]) * (s->lmax - s->lmin) / (double)m * NORM;
    }
}

/* display levels implied by an XYZ error box around a reference colour */
static double display_levels(const double ref[3], const double dv[3]) {
    double lin[3], drgb[3];
    xyz_to_linear(ref, lin);
    double worst = 0.0;
    for (int i = 0; i < 3; i++) {
        drgb[i] = fabs(XYZ2RGB[i][0]) * dv[0] + fabs(XYZ2RGB[i][1]) * dv[1] +
                  fabs(XYZ2RGB[i][2]) * dv[2];
        double lvl = 255.0 * gamma_slope(lin[i]) * drgb[i];
        if (lvl > worst) worst = lvl;
    }
    return worst;
}

/* Upper bound on CIE76 dE for any XYZ inside the per-channel box
   ref[ch] +- dv[ch].  Each Lab component is monotone in X, Y and Z
   separately, so its extreme lies on a box corner; the norm of the
   componentwise extremes dominates the norm at every point of the box. */
static double lab_bound(const double ref[3], const double dv[3], const double labr[3]) {
    double span[3] = { 0, 0, 0 };
    for (int c = 0; c < 8; c++) {
        double xyz[3], lab[3];
        for (int ch = 0; ch < 3; ch++)
            xyz[ch] = ref[ch] + ((c >> ch) & 1 ? dv[ch] : -dv[ch]);
        xyz_to_lab(xyz, lab);
        for (int i = 0; i < 3; i++) {
            double d = fabs(lab[i] - labr[i]);
            if (d > span[i]) span[i] = d;
        }
    }
    return sqrt(span[0] * span[0] + span[1] * span[1] + span[2] * span[2]);
}

/* exposure that maps the reference to a chosen peak linear channel */
static double exposure_to(const double ref[3], double peak) {
    double lin[3];
    xyz_to_linear(ref, lin);
    double m = fmax(lin[0], fmax(lin[1], lin[2]));
    return m > 1e-30 ? peak / m : 0.0;
}

/* ---------------------------------------------------------- test harness */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* bracket one scene at one sample count, return 1 if inside the bound */
static int bracket_scene(const Scene* s, int K, long N, unsigned seed) {
    Estimate e = estimate(s, K, N, seed);
    double ref[3];
    /* full observer range: any colour the render box drops is an error too */
    integrate_box(s, LAM_MIN, LAM_MAX, REF_M, ref);
    double labr[3], labe[3];
    xyz_to_lab(ref, labr);
    xyz_to_lab(e.xyz, labe);
    double de = delta_e(labe, labr);
    double zmax = 0.0, err_rel = 0.0;
    int inside = e.ok;
    for (int ch = 0; ch < 3; ch++) {
        double err = fabs(e.xyz[ch] - ref[ch]);
        double scale = fabs(ref[ch]) > 1e-12 ? fabs(ref[ch]) : 1e-12;
        err_rel = fmax(err_rel, err / scale);
        if (e.se[ch] > 0.0) zmax = fmax(zmax, err / e.se[ch]);
        if (err > e.bound[ch] + 1e-12 * scale) inside = 0;
    }
    double labb = lab_bound(ref, e.bound, labr);
    double A = exposure_to(ref, 0.7);           /* display-level test at mid grey */
    double dvx[3] = { A * e.bound[0], A * e.bound[1], A * e.bound[2] };
    double refx[3] = { A * ref[0], A * ref[1], A * ref[2] };
    double lvl = display_levels(refx, dvx);
    printf("scene=%-10s N=%-6ld K=%-3d rel_err=%.3e z=%.2f dE=%.4f dE_bound=%.4f levels<=%.2f trunc=%.2e ok=%d\n",
           s->name, N, K, err_rel, zmax, de, labb, lvl, e.trunc[1], inside);
    return inside;
}

int main(int argc, char** argv) {
    char msg[256];
    Scene scenes[4];
    Reflector refl = { 0.9 };
    Band line = { 555.0, 5.0 };
    Band narrow = { 620.0, 12.0 };

    scenes[0] = (Scene){ "flat",   flat_P,    &refl,  LAM_MIN, LAM_MAX };
    scenes[1] = (Scene){ "d65",    d65_P,     &refl,  LAM_MIN, LAM_MAX };
    scenes[2] = (Scene){ "bb5800", black_P,   &(Black){ 5800.0 }, LAM_MIN, LAM_MAX };
    scenes[3] = (Scene){ "cut600", band_P,    &narrow, LAM_MIN, 600.0 }; /* box clips red */
    Scene line_scene = { "line555", band_P, &line, LAM_MIN, LAM_MAX };

    /* ---- observer / colour tables ---- */
    check("cie_y_peak_555", fabs(cmf(1, 555.0) - 1.0) < 1e-12);
    check("cie_table_ends", cmf(1, CIE_FIRST) < 1e-4 && cmf(1, CIE_FIRST + CIE_STEP * (CIE_N - 1)) < 1e-4);
    check("cie_channels_differ", fabs(cmf(0, 600.0) - cmf(2, 600.0)) > 0.5);
    {
        double xyz[3];
        integrate_box(&scenes[0], LAM_MIN, LAM_MAX, REF_M, xyz); /* albedo .9 flat */
        double ratio[3] = { xyz[0] / xyz[1], xyz[1] / xyz[1], xyz[2] / xyz[1] };
        check("equal_energy_near_equal", fabs(ratio[0] - 1.0) < 1e-3 && fabs(ratio[2] - 1.0) < 1e-3);
    }
    {
        double xyz[3];
        integrate_box(&scenes[1], LAM_MIN, LAM_MAX, REF_M, xyz);
        check("d65_whitepoint_srgb", fabs(xyz[0] / xyz[1] - 0.95047) < 2e-3 &&
                                      fabs(xyz[2] / xyz[1] - 1.08883) < 2e-3);
        double lin[3];
        xyz_to_linear((double[3]){ xyz[0] / xyz[1], 1.0, xyz[2] / xyz[1] }, lin);
        check("d65_maps_to_neutral", fabs(lin[0] - lin[1]) < 2e-3 && fabs(lin[2] - lin[1]) < 2e-3);
    }

    /* ---- range validation reports the legal bounds ---- */
    double a = 900.0, b = 950.0;
    check("range_rejects_outside", check_range(a, b, msg, sizeof msg) == RANGE_OUTSIDE);
    check("range_msg_names_360_830", strstr(msg, "360") != NULL && strstr(msg, "830") != NULL);
    a = 500.0; b = 500.0;
    check("range_rejects_empty", check_range(a, b, msg, sizeof msg) == RANGE_ORDER);
    a = NAN; b = 600.0;
    check("range_rejects_nan", check_range(a, b, msg, sizeof msg) == RANGE_NOT_FINITE);
    a = 300.0; b = 700.0;
    check("range_rejects_low", check_range(a, b, msg, sizeof msg) == RANGE_OUTSIDE);
    a = 500.0; b = 500.005;
    check("range_rejects_narrow", check_range(a, b, msg, sizeof msg) == RANGE_TOO_NARROW);
    a = 400.0; b = 700.0;
    check("range_accepts_legal", check_range(a, b, msg, sizeof msg) == RANGE_OK);

    /* the estimator refuses an illegal range instead of returning numbers */
    {
        Scene bad = { "bad", flat_P, &refl, 900.0, 950.0 };
        Estimate e = estimate(&bad, 32, 1024, 1u);
        check("bad_range_no_estimate", !e.ok && e.n == 0);
    }

    /* ---- reference convergence: halving the step must not move it ---- */
    {
        double r1[3], r2[3];
        integrate_box(&scenes[1], LAM_MIN, LAM_MAX, REF_M, r1);
        integrate_box(&scenes[1], LAM_MIN, LAM_MAX, 2 * REF_M, r2);
        double d = 0.0;
        for (int ch = 0; ch < 3; ch++)
            d = fmax(d, fabs(r1[ch] - r2[ch]) / fabs(r1[ch]));
        printf("reference_refinement_rel_diff=%.3e\n", d);
        check("reference_is_dense", d < 1e-6);
    }

    /* ---- importance map actually reduces variance on every scene ---- */
    {
        printf("--- variance reduction (importance vs uniform, same draws) ---\n");
        int all_better = 1;
        for (int i = 0; i < 4; i++) {
            Estimate e = estimate(&scenes[i], 32, 8192, 4242u);
            double worst = 0.0;
            for (int ch = 0; ch < 3; ch++) {
                double r = e.unif[ch] > 0.0 ? e.se[ch] / e.unif[ch] : 1.0;
                if (r > 0.9) all_better = 0;
                if (r > worst) worst = r;
            }
            printf("scene=%-10s best_se_ratio=%.3f\n", scenes[i].name, worst);
        }
        check("importance_beats_uniform", all_better);
    }

    /* ---- bracket: observed error inside 3*se + truncation, several seeds -- */
    {
        int all_inside = 1;
        printf("--- bound vs dense reference ---\n");
        for (int i = 0; i < 4; i++) {
            for (int seed = 0; seed < 6; seed++) {
                int ok = bracket_scene(&scenes[i], 64, 4096, 100u + (unsigned)seed);
                if (!ok) all_inside = 0;
            }
        }
        all_inside &= bracket_scene(&line_scene, 64, 4096, 7u);
        check("bracket_holds", all_inside);
    }

    /* ---- truncation bound is real: the boxed scene loses red energy ---- */
    {
        double ref_box[3], ref_full[3], trunc[3];
        integrate_box(&scenes[3], scenes[3].lmin, scenes[3].lmax, REF_M, ref_box);
        integrate_box(&scenes[3], LAM_MIN, LAM_MAX, REF_M, ref_full);
        integrate_outside(&scenes[3], trunc);
        double miss = 0.0;
        for (int ch = 0; ch < 3; ch++) {
            double d = fabs(ref_full[ch] - ref_box[ch]);
            if (d > trunc[ch] * (1.0 + 1e-9) + 1e-12) miss = 1.0;
        }
        check("truncation_covers_lost_colour", miss == 0.0);
        check("truncation_is_significant", trunc[0] > 1e-3 * ref_box[0]);
    }

    /* ---- convergence: more wavelength samples -> smaller dE and bound ---- */
    {
        Scene* s = &scenes[1];
        double ref[3], labr[3];
        integrate_box(s, s->lmin, s->lmax, REF_M, ref);
        xyz_to_lab(ref, labr);
        long counts[] = { 64, 256, 1024, 4096, 16384, 65536 };
        double mean_de[6], mean_bound[6];
        printf("--- convergence (d65 scene, 12 seeds) ---\n");
        for (int k = 0; k < 6; k++) {
            int K = (int)(counts[k] / 8) < 8 ? 8 : (int)(counts[k] / 8);
            if (K > 256) K = 256;
            double de_sum = 0.0, bd_sum = 0.0, zmax = 0.0;
            for (int seed = 0; seed < 12; seed++) {
                Estimate e = estimate(s, K, counts[k], 500u + (unsigned)seed);
                double labe[3];
                xyz_to_lab(e.xyz, labe);
                de_sum += delta_e(labe, labr);
                bd_sum += lab_bound(ref, e.bound, labr);
                for (int ch = 0; ch < 3; ch++)
                    if (e.se[ch] > 0.0)
                        zmax = fmax(zmax, fabs(e.xyz[ch] - ref[ch]) / e.se[ch]);
            }
            mean_de[k] = de_sum / 12.0;
            mean_bound[k] = bd_sum / 12.0;
            printf("N=%-6ld mean_dE=%.5f mean_dE_bound=%.5f max_z=%.2f\n",
                   counts[k], mean_de[k], mean_bound[k], zmax);
            if (zmax > 4.0) { failures++; printf("z_exceeds_4sigma=%d\n", 1); }
        }
        int decreasing = 1;
        for (int k = 1; k < 6; k++)
            if (!(mean_de[k] <= mean_de[k - 1] * 1.05 + 1e-6)) decreasing = 0;
        check("dE_halves_with_sampling", decreasing);
        check("dE_converges_sub1", mean_de[5] < 0.2);
        double shrink = mean_bound[0] / (mean_bound[5] > 1e-12 ? mean_bound[5] : 1e-12);
        check("bound_shrinks_with_sqrt_n", shrink > 20.0);
        check("bound_covers_error", mean_bound[5] + 1e-6 >= mean_de[5]);
    }

    /* ---- the rainbow: hue must run monotonically blue -> green -> red ---- */
    {
        const int nb = 54;
        double hue[54], lab[54][3], xyz[54][3];
        int bad_imp = 0;
        printf("--- spectral bands 430..695 nm, 5 nm steps ---\n");
        for (int i = 0; i < nb; i++) {
            Band b = { 430.0 + 5.0 * i, 6.0 };
            Scene sc = { "band", band_P, &b, LAM_MIN, LAM_MAX };
            integrate_box(&sc, LAM_MIN, LAM_MAX, REF_M, xyz[i]);
            xyz_to_lab(xyz[i], lab[i]);
            hue[i] = lab_hue(lab[i]);
            /* importance estimator at a modest sample count vs the reference */
            Estimate e = estimate(&sc, 64, 4096, 3u);
            double labe[3];
            xyz_to_lab(e.xyz, labe);
            /* the same bracket the other scenes use, in dE units */
            if (delta_e(labe, lab[i]) > lab_bound(xyz[i], e.bound, lab[i])) bad_imp++;
        }
        check("importance_bands_accurate", bad_imp == 0);
        /* unwrap and require a monotone descent */
        double un[54];
        un[0] = hue[0];
        for (int i = 1; i < nb; i++) {
            double h = hue[i];
            while (h > un[i - 1] + 180.0) h -= 360.0;
            while (h < un[i - 1] - 180.0) h += 360.0;
            un[i] = h;
        }
        int mono = 1;
        for (int i = 1; i < nb; i++)
            if (!(un[i] <= un[i - 1] + 1e-6)) mono = 0;
        check("hue_monotone_rainbow", mono);
        check("hue_span_wide", un[0] - un[nb - 1] > 280.0);
        check("violet_band_hue", hue[3] > 300.0 && hue[3] < 330.0);      /* 445 nm */
        check("cyan_band_hue", hue[12] > 150.0 && hue[12] < 230.0);      /* 490 nm */
        check("green_band_hue", hue[25] > 100.0 && hue[25] < 145.0);     /* 555 nm */
        check("red_band_hue", hue[41] > 20.0 && hue[41] < 45.0);         /* 635 nm */

        /* the fixed-wavelength estimator cannot follow the band: same dE test */
        const double fixed[3] = { 450.0, 550.0, 600.0 };
        int naive_bad = 0;
        double naive_worst = 0.0;
        for (int i = 0; i < nb; i++) {
            Band b = { 430.0 + 5.0 * i, 6.0 };
            Scene sc = { "band", band_P, &b, LAM_MIN, LAM_MAX };
            double xyz_n[3], lab_n[3];
            estimate_fixed(&sc, fixed, 3, xyz_n);
            xyz_to_lab(xyz_n, lab_n);
            double de = delta_e(lab_n, lab[i]);
            if (de > 10.0) naive_bad++;
            if (de > naive_worst) naive_worst = de;
        }
        printf("fixed-wavelength control: bands with dE>10: %d/%d, worst dE=%.1f\n",
               naive_bad, nb, naive_worst);
        check("fixed_sampling_fails_bands", naive_bad > nb / 2);
        check("fixed_sampling_worst_large", naive_worst > 30.0);
    }

    /* ---- optional image: a dispersive grating, one column per direction ---- */
    if (argc > 1) {
        const int W = 256, H = 64;
        const int MS = 4;               /* rows per column: mean over MS draws */
        double lin_col[256][3], lum[256];
        FILE* f = fopen(argv[1], "wb");
        if (f) {
            fprintf(f, "P6\n%d %d\n255\n", W, H);
            for (int x = 0; x < W; x++) {
                /* grating: deflection(lambda) = angle0 + D*(lambda - 520); the
                   image column at angle0 collects a band centred on lambda. */
                double center = 400.0 + (700.0 - 400.0) * (double)x / (double)(W - 1);
                Band b = { center, 6.0 };
                Scene sc = { "grating", band_P, &b, LAM_MIN, LAM_MAX };
                double acc[3] = { 0, 0, 0 }, accl = 0.0;
                for (int m = 0; m < MS; m++) {
                    Estimate e = estimate(&sc, 128, 8192, 31u + (unsigned)(x * MS + m));
                    double lin[3];
                    xyz_to_linear(e.xyz, lin);
                    for (int ch = 0; ch < 3; ch++) acc[ch] += lin[ch];
                    accl += e.xyz[1];
                }
                for (int ch = 0; ch < 3; ch++) lin_col[x][ch] = acc[ch] / MS;
                lum[x] = accl / MS;
            }
            /* the band at the ends of the visible range is almost black; smooth
               the strip with a small box kernel so the MC noise there does not
               show up as dither, then find the exposure peak */
            double lpv[256][3] = { { 0 } };
            for (int x = 0; x < W; x++) {
                int lo = x > 4 ? x - 4 : 0, hi = x + 4 < W - 1 ? x + 4 : W - 1;
                for (int j = lo; j <= hi; j++)
                    for (int ch = 0; ch < 3; ch++) lpv[x][ch] += lin_col[j][ch];
                for (int ch = 0; ch < 3; ch++) lin_col[x][ch] = lpv[x][ch] / (double)(hi - lo + 1);
            }
            double lum_peak = 0.0;
            for (int x = 0; x < W; x++) if (lum[x] > lum_peak) lum_peak = lum[x];
            for (int y = 0; y < H; y++) {
                for (int x = 0; x < W; x++) {
                    /* one exposure for the whole strip: relative luminance
                       comes from XYZ Y, so the spectral locus keeps its own
                       brightness profile (dark at the ends) */
                    double k = 0.75 / lum_peak;
                    unsigned char px[3];
                    for (int ch = 0; ch < 3; ch++) {
                        double v = gamma_encode(lin_col[x][ch] * k);
                        if (v < 0.0) v = 0.0;
                        if (v > 1.0) v = 1.0;
                        px[ch] = (unsigned char)(255.0 * v + 0.5);
                    }
                    fwrite(px, 1, 3, f);
                }
            }
            fclose(f);
            fprintf(stderr, "wrote %s (%dx%d)\n", argv[1], W, H);
        }
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
