
/* dof: thin-lens circle of confusion with pixel-diagonal clamp, bokeh self-check.
 *
 * model (all lengths in mm):
 *   r_lens = f^2 |s-d| / (2 N d (s-f))     thin-lens CoC radius
 *   r_eff  = sqrt(r_lens^2 + r_min^2)      clamp, r_min = p_diag/2
 * p_diag is the pixel diagonal; r_min is the radius of the disc reconstruction
 * filter, so convolution second moments add and the two limits are
 *   d -> s   : r_eff -> p_diag/2            (never collapses to zero)
 *   d -> inf : r_eff -> sqrt(r_inf^2 + r_min^2), r_inf = f^2 / (2 N (s-f))
 *
 * measurement: blobs are rendered on a sub-pixel grid (pitch p/8) so the
 * discrete second moment matches the continuous one; the display image itself
 * is rendered at pitch p with the same filter radius.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

/* ---------- analytic model ---------- */

static double coc_lens(double f, double N, double s, double d) {
    return f * f * fabs(s - d) / (2.0 * N * d * (s - f));
}

static double coc_effective(double f, double N, double s, double d, double r_min) {
    double r = coc_lens(f, N, s, d);
    return sqrt(r * r + r_min * r_min);
}

/* aperture bounds: 1 <= N <= 64, diameter f/N <= 100mm, object beyond focal length */
static int aperture_check(double f, double N, double d, const char** why) {
    if (!(N >= 1.0 && N <= 64.0)) { *why = "f-number N outside [1,64]"; return 0; }
    if (f / N > 100.0) { *why = "aperture diameter f/N > 100mm"; return 0; }
    if (d <= f) { *why = "object distance <= focal length, no real image"; return 0; }
    return 1;
}

/* ---------- renderer: thin-lens camera, small emissive spheres ---------- */

typedef struct { double x, z, r; } Sphere; /* centers on y=0 plane */

static unsigned long long rng_s;
static void rng_seed(unsigned long long s) {
    rng_s = s ? s : 0x123456789ULL;
    rng_s ^= rng_s >> 12; rng_s ^= rng_s << 25; rng_s ^= rng_s >> 27;
}
static double rnd(void) {
    rng_s ^= rng_s >> 12; rng_s ^= rng_s << 25; rng_s ^= rng_s >> 27;
    return (double)((rng_s * 0x2545F4914F6CDD1DULL) >> 11) * (1.0 / 9007199254740992.0);
}

/* render a square sensor window: center (cx,cy) mm, half-size `half` mm, n x n
 * samples, disc reconstruction filter radius r_f. lens at origin, sensor at
 * z=-si, scene at z>0. */
static void render_window(double* img, int n, double cx, double cy, double half,
                          double f, double N, double s, double r_f,
                          const Sphere* sp, int nsp, int grid) {
    double pitch = 2.0 * half / n;
    double si = f * s / (s - f);
    double aper = f / N;
    for (int j = 0; j < n; j++) for (int i = 0; i < n; i++) {
        double bx = cx - half + (i + 0.5) * pitch;
        double by = cy - half + (j + 0.5) * pitch;
        rng_seed(0x9e3779b97f4a7c15ULL + (unsigned long long)(j * n + i) * 0xbf58476d1ce4e5b9ULL);
        double acc = 0.0;
        for (int gj = 0; gj < grid; gj++) for (int gi = 0; gi < grid; gi++) {
            double u1 = (gi + rnd()) / grid, u2 = (gj + rnd()) / grid;
            double pr = r_f * sqrt(u1), pa = 2.0 * PI * u2;
            double px = bx + pr * cos(pa), py = by + pr * sin(pa);
            double ar = 0.5 * aper * sqrt(rnd()), aa = 2.0 * PI * rnd();
            double ax = ar * cos(aa), ay = ar * sin(aa);
            /* ray from aperture point through the pixel's focal-plane point */
            double fx = -px * s / si, fy = -py * s / si;
            double dx = fx - ax, dy = fy - ay, dz = s;
            double il = 1.0 / sqrt(dx * dx + dy * dy + dz * dz);
            dx *= il; dy *= il; dz *= il;
            for (int k = 0; k < nsp; k++) {
                double ox = ax - sp[k].x, oy = ay, oz = -sp[k].z;
                double b = ox * dx + oy * dy + oz * dz;
                double c = ox * ox + oy * oy + oz * oz - sp[k].r * sp[k].r;
                double h = b * b - c;
                if (h > 0.0 && -b - sqrt(h) > 0.0) { acc += 1.0; break; }
            }
        }
        img[j * n + i] = acc / (grid * grid);
    }
}

/* second-moment blob radius in a window centered on the blob.
 * disc PSF: <r^2> = r^2/2, convolution adds second moments, so
 * r_meas = sqrt(2*m2 - r_geo^2) recovers sqrt(r_lens^2 + r_min^2). */
static double measure_window(const double* img, int n, double half,
                             double r_geo, double* out_mass) {
    double pitch = 2.0 * half / n;
    double sum = 0.0, mx = 0.0, my = 0.0;
    for (int j = 0; j < n; j++) for (int i = 0; i < n; i++) {
        double w = img[j * n + i];
        double x = -half + (i + 0.5) * pitch, y = -half + (j + 0.5) * pitch;
        sum += w; mx += w * x; my += w * y;
    }
    if (sum <= 0.0) { *out_mass = 0.0; return -1.0; }
    mx /= sum; my /= sum;
    double m2 = 0.0;
    for (int j = 0; j < n; j++) for (int i = 0; i < n; i++) {
        double w = img[j * n + i];
        double x = -half + (i + 0.5) * pitch - mx, y = -half + (j + 0.5) * pitch - my;
        m2 += w * (x * x + y * y);
    }
    m2 /= sum;
    *out_mass = sum;
    double r2 = 2.0 * m2 - r_geo * r_geo;
    return r2 > 0.0 ? sqrt(r2) : 0.0;
}

/* measure one blob at sensor position (sx,sy): render a tight window at
 * sub-pixel pitch p/8 and take its second-moment radius */
static double measure_blob(double sx, double sy, double r_geo, double r_expect,
                           double f, double N, double s, double r_min, double p,
                           const Sphere* sp, int nsp, double* mass) {
    double half = 3.0 * r_expect + r_geo + 2.0 * p;
    int n = (int)ceil(2.0 * half / (p / 8.0)) + 1;
    double* img = malloc(sizeof(double) * (size_t)n * n);
    if (!img) { *mass = 0.0; return -1.0; }
    render_window(img, n, sx, sy, half, f, N, s, r_min, sp, nsp, 8);
    double rm = measure_window(img, n, half, r_geo, mass);
    free(img);
    return rm;
}

static int write_ppm(const char* path, const double* img, int R) {
    FILE* fp = fopen(path, "wb");
    if (!fp) return -1;
    double mx = 1e-30;
    for (int i = 0; i < R * R; i++) if (img[i] > mx) mx = img[i];
    fprintf(fp, "P6\n%d %d\n255\n", R, R);
    for (int i = 0; i < R * R; i++) {
        double v = sqrt(img[i] / mx);
        unsigned char b = (unsigned char)(v * 255.0 + 0.5);
        fputc(b, fp); fputc(b, fp); fputc(b, fp);
    }
    fclose(fp);
    return 0;
}

/* ---------- bokeh: three focal lengths, five defocus distances each ---------- */

static int bokeh(const char* outdir) {
    const double s = 2000.0, N = 2.0, W = 36.0;
    const int R = 400;
    const double fs[3] = {25.0, 50.0, 85.0};
    const double ds[5] = {1000.0, 1400.0, 1700.0, 2600.0, 5000.0};
    const double fracs[5] = {-0.7, -0.35, 0.0, 0.35, 0.7};
    const double p = W / R;
    const double r_min = p * sqrt(2.0) / 2.0;
    double* img = malloc(sizeof(double) * R * R);
    if (!img) return 1;
    for (int fi = 0; fi < 3; fi++) {
        double f = fs[fi];
        const char* why;
        if (!aperture_check(f, N, s, &why)) {
            printf("f=%.0f N=%.1f skipped: %s\n", f, N, why);
            continue;
        }
        double si = f * s / (s - f);
        Sphere sp[5]; double sx[5], rgeo[5];
        for (int k = 0; k < 5; k++) {
            double d = ds[k];
            sp[k].x = -fracs[k] * (W / 2.0) * d / si;
            sp[k].z = d;
            sp[k].r = 0.002 * d;
            sx[k] = fracs[k] * (W / 2.0);
            rgeo[k] = sp[k].r * (f * d / (d - f)) / d;
        }
        render_window(img, R, 0.0, 0.0, W / 2.0, f, N, s, r_min, sp, 5, 12);
        char path[512];
        snprintf(path, sizeof path, "%s/bokeh_f%02d.ppm", outdir, (int)f);
        if (write_ppm(path, img, R) != 0) { printf("write failed: %s\n", path); free(img); return 1; }
        printf("f=%.0fmm N=%.1f s=%.0fmm r_min=p_diag/2=%.4fmm -> %s\n", f, N, s, r_min, path);
        printf("  %8s %9s %9s %9s %8s\n", "d(mm)", "r_lens", "r_eff", "r_meas", "err(%)");
        for (int k = 0; k < 5; k++) {
            double d = ds[k];
            double rl = coc_lens(f, N, s, d);
            double re = coc_effective(f, N, s, d, r_min);
            double mass = 0.0;
            double rm = measure_blob(sx[k], 0.0, rgeo[k], re, f, N, s, r_min, p, sp, 5, &mass);
            printf("  %8.0f %9.4f %9.4f %9.4f %8.2f%s\n",
                   d, rl, re, rm, rm > 0.0 ? 100.0 * (rm - re) / re : 0.0,
                   rl < r_min ? " clamped" : "");
        }
    }
    /* aperture out-of-bounds probes: reported separately, never rendered */
    printf("aperture out-of-bounds (reported, not rendered):\n");
    const double pf[4] = {85.0, 400.0, 50.0, 50.0};
    const double pn[4] = {0.8, 2.0, 128.0, 2.0};
    const double pd[4] = {1000.0, 3000.0, 1000.0, 30.0};
    int n_oob = 0;
    for (int i = 0; i < 4; i++) {
        const char* why;
        if (!aperture_check(pf[i], pn[i], pd[i], &why)) {
            printf("  f=%.0fmm N=%.1f d=%.0fmm: %s\n", pf[i], pn[i], pd[i], why);
            n_oob++;
        }
    }
    if (!n_oob) printf("  none\n");
    free(img);
    return 0;
}

/* ---------- self-test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int run_tests(void) {
    const double f = 50.0, N = 2.0, s = 2000.0;
    check("coc_value", fabs(coc_lens(f, N, s, 1000.0) - 0.320513) < 1e-5);

    /* near limit: at the focal plane the CoC converges to p_diag/2, not zero */
    double rmin = 0.1591;
    check("clamp_at_focus", fabs(coc_effective(f, N, s, s, rmin) - rmin) < 1e-12);

    /* far limit: r_lens -> f^2/(2N(s-f)), r_eff -> quadrature with r_min */
    double rinf = f * f / (2.0 * N * (s - f));
    check("far_limit_analytic", fabs(coc_lens(f, N, s, 1e15) / rinf - 1.0) < 1e-9);
    check("far_limit_clamped",
          fabs(coc_effective(f, N, s, 1e15, rmin) / sqrt(rinf * rinf + rmin * rmin) - 1.0) < 1e-9);

    check("monotone_near", coc_lens(f, N, s, 1500.0) < coc_lens(f, N, s, 1200.0));
    check("monotone_far", coc_lens(f, N, s, 2500.0) < coc_lens(f, N, s, 4000.0));

    const char* why;
    check("aperture_ok", aperture_check(50.0, 2.0, 1000.0, &why));
    check("aperture_low_N", !aperture_check(50.0, 0.5, 1000.0, &why));
    check("aperture_too_wide", !aperture_check(400.0, 2.0, 3000.0, &why));
    check("aperture_inside_f", !aperture_check(50.0, 2.0, 30.0, &why));

    /* measured near limit: point on the focal plane renders as a p_diag/2 blob */
    {
        const double W = 36.0, p = W / 160.0;
        double r_min = p * sqrt(2.0) / 2.0;
        Sphere sp = {0.0, s, 0.002 * s};
        double rgeo = sp.r * (f * s / (s - f)) / s;
        double mass = 0.0;
        double rm = measure_blob(0.0, 0.0, rgeo, r_min, f, N, s, r_min, p, &sp, 1, &mass);
        printf("near_limit_meas=%.4f expect=%.4f\n", rm, r_min);
        check("near_limit_measured", mass > 0.0 && fabs(rm / r_min - 1.0) < 0.10);
    }
    /* measured far limit: point at 1e6 mm sits at the d->inf asymptote */
    {
        const double W = 36.0, p = W / 200.0;
        double r_min = p * sqrt(2.0) / 2.0;
        double d = 1e6;
        Sphere sp = {0.0, d, 0.002 * d};
        double rgeo = sp.r * (f * d / (d - f)) / d;
        double re = coc_effective(f, N, s, d, r_min);
        double mass = 0.0;
        double rm = measure_blob(0.0, 0.0, rgeo, re, f, N, s, r_min, p, &sp, 1, &mass);
        printf("far_limit_meas=%.4f expect=%.4f\n", rm, re);
        check("far_limit_measured", mass > 0.0 && fabs(rm / re - 1.0) < 0.10);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc > 1 && strcmp(argv[1], "bokeh") == 0)
        return bokeh(argc > 2 ? argv[2] : "out");
    return run_tests();
}
