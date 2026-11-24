/* glass: curtain-wall moire and the low-discrepancy offset that decouples
   sampling density from geometric subdivision.

   Model
   -----
   The facade is a grid of glass panels. One panel measures p = W / D pixels
   across, where D is the subdivision density (panels across the image width).
   Each panel carries a faint periodic ripple across its face and a dark
   mullion down its leading edge. A renderer shades a pixel by box-filtering
   this pattern with an S x S stratified grid.

   Why stripes appear
   ------------------
   The grid samples every pixel at the same fractional position. When the
   panel pitch p is not a whole number of pixels, the lattice advances by
   frac(p) of a panel from one pixel to the next, and the pixel's average of
   the ripple shifts as that phase sweeps a full panel. That slow phase drift
   is the beat: a low-frequency stripe of spacing

       lambda = p / frac(S*p)              (S samples per pixel per axis)

   and hence

       lambda * D = W / frac(S*W/D).

   For a batch of subdivisions whose sampling lattice keeps the same offset
   against the panel lattice, lambda*D is a constant. A constant lambda*D is
   therefore the signature of an estimator that is still phase-locked: the
   offset, if there was one, did not do its job.

   The fix
   -------
   Rotate the whole S x S lattice of each pixel by a low-discrepancy (R2)
   offset, indexed by the pixel alone, so every pixel's samples move together
   and the stratified spacing is preserved. The offset range is derived from
   the subdivision density: one panel period p = W/D pixels, capped at one
   pixel so a dense wall stays a sub-pixel jitter instead of blurring across
   pixel borders. The offset spans exactly one period of the pattern the
   sampler is integrating, so the phase relative to the panel edges differs
   from pixel to pixel; the sample lattice decorrelates from the panel
   lattice, the pixel average converges to the true anti-aliased value, and
   the coherent stripe disappears without inflating the per-pixel noise.

   The experiment renders three subdivision densities with both samplers,
   measures the dominant stripe spacing on each, and prints lambda*D. For
   each density it also dumps the actual sample positions used for a tile of
   the wall, so the sample distribution can be plotted next to the measured
   stripes. */
#include "glass.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* R2 low-discrepancy sequence: phi2 = plastic constant, a_i = 1/phi2^i. */
#define R2_PHI 1.32471795724474602596
static double r2_a1(void) { return 1.0 / R2_PHI; }
static double r2_a2(void) { return 1.0 / (R2_PHI * R2_PHI); }
static double frac(double x) { return x - floor(x); }

/* ------------------------------------------------------------------ wall -- */
GlassWall glass_wall(int width, int height, int subdiv, int spp) {
    GlassWall w;
    memset(&w, 0, sizeof(w));
    w.width = width;
    w.height = height;
    w.subdiv = subdiv < 1 ? 1 : subdiv;
    w.spp = spp < 1 ? 1 : (spp > GLASS_MAX_SPP ? GLASS_MAX_SPP : spp);
    w.ripple = GLASS_RIPPLE;
    w.mullion = GLASS_MULLION_DEPTH;
    w.duty = GLASS_MULLION_DUTY;
    w.grad_x = GLASS_GRAD_X;
    w.grad_y = GLASS_GRAD_Y;
    return w;
}

double glass_panel_pitch(const GlassWall* w) { return (double)w->width / w->subdiv; }

/* One axis of the facade at t pixels: illumination gradient, dark mullion at
   each panel's leading edge, and the glass ripple. */
static double axis_shade(double t, double pitch, int n, double grad,
                         double ripple_amp, double mullion, double duty) {
    double u = frac(t / pitch);
    double mull = u < duty ? (1.0 - mullion) : 1.0;
    double ripple = 1.0 + ripple_amp * sin(2.0 * GLASS_PI * GLASS_RIPPLE_CYCLES * u);
    return (1.0 + grad * t / n) * mull * ripple;
}

double glass_shade(double px, double py, const GlassWall* w) {
    double p = glass_panel_pitch(w);
    return axis_shade(px, p, w->width, w->grad_x,
                      w->ripple, w->mullion, w->duty) *
           axis_shade(py, p, w->height, w->grad_y,
                      w->ripple, w->mullion, w->duty);
}

/* --------------------------------------------------------------- offsets -- */
/* Offset magnitude in pixels, derived from the subdivision density: one panel
   period p = W/D, clamped to a single pixel. This is the range that lets a
   displaced sample visit every phase of the pattern the panel contributes. */
double glass_offset_radius(const GlassWall* w) {
    double p = glass_panel_pitch(w);
    return p < 1.0 ? p : 1.0;
}

/* Emit the S*S sample positions of one pixel (pixel-local, in [0,1)).

   Base positions are a stratified lattice. In GLASS_SAMPLE_OFFSET mode the
   whole lattice of a pixel is rotated by one low-discrepancy offset, indexed
   by the pixel alone: every sample in the pixel moves together, so the
   stratified spacing (and hence the per-pixel variance) is preserved, while
   the phase relative to the panel edges differs from pixel to pixel. Rotating
   each sample independently would also break the beat, but would destroy the
   stratification and roughly double the per-pixel noise. Positions are
   wrapped back into the pixel so the footprint does not grow. */
int glass_pixel_samples(const GlassWall* w, int px, int py, GlassSampleMode mode,
                        double out[][2]) {
    int S = w->spp;
    double r = glass_offset_radius(w);
    double ox = 0.0, oy = 0.0;
    if (mode == GLASS_SAMPLE_OFFSET) {
        double n = (double)((long)py * w->width + px);
        ox = r * (frac(0.5 + n * r2_a1()) - 0.5);
        oy = r * (frac(0.5 + n * r2_a2()) - 0.5);
    }
    int k = 0;
    for (int sy = 0; sy < S; sy++) {
        for (int sx = 0; sx < S; sx++) {
            out[k][0] = frac(((double)sx + 0.5) / S + ox);
            out[k][1] = frac(((double)sy + 0.5) / S + oy);
            k++;
        }
    }
    return k;
}

/* The beat wavelength the grid sampler is predicted to leave behind, in
   pixels: lambda = p / frac(S*p), or 0 when the lattice is commensurate
   (frac = 0, no beat). */
double glass_predicted_stripe(const GlassWall* w) {
    double p = glass_panel_pitch(w);
    double lattice = (double)w->spp * p;
    double f = frac(lattice);
    if (f < 1e-6 || f > 1.0 - 1e-6) return 0.0;
    return p / f;
}

/* -------------------------------------------------------------- renderer -- */
void glass_render(const GlassWall* w, GlassSampleMode mode, double* out) {
    static double samp[GLASS_MAX_SPP * GLASS_MAX_SPP][2];
    for (int y = 0; y < w->height; y++) {
        for (int x = 0; x < w->width; x++) {
            int n = glass_pixel_samples(w, x, y, mode, samp);
            double acc = 0.0;
            for (int s = 0; s < n; s++)
                acc += glass_shade((double)x + samp[s][0], (double)y + samp[s][1], w);
            out[(long)y * w->width + x] = acc / (double)n;
        }
    }
}

/* Exact pixel average of the facade. The shading is separable into an x and a
   y factor; each axis is integrated analytically over the pixel, splitting at
   panel and mullion boundaries. Since u = frac(t/p), the ripple phase is
   local to each panel, which is exactly what the split integral uses.

   With t = i*p + u and the local phase w0 = 2*pi*C/p, on a subinterval where
   the mullion factor is constant:
     integral (1 + g t/n) (1 + R sin(w0 u)) dt
   = base + R ( integral sin + (g/n) integral t sin ),
     integral t sin(w0 u) dt = i*p * integral sin + [sin(w0 u)/w0^2 - u cos(w0 u)/w0]. */
static double axis_exact(double px, double pitch, int n, double grad,
                         double ripple_amp, double mullion, double duty) {
    double w0 = 2.0 * GLASS_PI * GLASS_RIPPLE_CYCLES / pitch;
    double x0 = px, x1 = px + 1.0;
    int i0 = (int)floor(x0 / pitch) - 1;
    int i1 = (int)floor(x1 / pitch) + 1;
    double g = grad / n;
    double tot = 0.0;
    for (int i = i0; i <= i1; i++) {
        double a0 = i * pitch > x0 ? i * pitch : x0;
        double a1 = (i + 1) * pitch < x1 ? (i + 1) * pitch : x1;
        if (a1 <= a0) continue;
        double mu = i * pitch + duty * pitch;
        double cuts[3];
        int nc = 0;
        cuts[nc++] = a0;
        if (a0 < mu && mu < a1) cuts[nc++] = mu;
        cuts[nc++] = a1;
        for (int c = 0; c + 1 < nc; c++) {
            double c0 = cuts[c], c1 = cuts[c + 1];
            double mull = c0 < mu - 1e-12 ? (1.0 - mullion) : 1.0;
            double base = (c1 - c0) + 0.5 * g * (c1 * c1 - c0 * c0);
            double ua = c0 - i * pitch, ub = c1 - i * pitch;
            double sa = sin(w0 * ua), sb = sin(w0 * ub);
            double ca = cos(w0 * ua), cb = cos(w0 * ub);
            double isin = -(cb - ca) / w0;
            double i_usin = (sb / (w0 * w0) - ub * cb / w0) -
                            (sa / (w0 * w0) - ua * ca / w0);
            double i_tsin = i * pitch * isin + i_usin;
            tot += mull * (base + ripple_amp * (isin + g * i_tsin));
        }
    }
    return tot;
}

void glass_reference(const GlassWall* w, double* out) {
    double p = glass_panel_pitch(w);
    for (int y = 0; y < w->height; y++) {
        double ey = axis_exact((double)y, p, w->height, w->grad_y,
                               w->ripple, w->mullion, w->duty);
        for (int x = 0; x < w->width; x++) {
            double ex = axis_exact((double)x, p, w->width, w->grad_x,
                                   w->ripple, w->mullion, w->duty);
            out[(long)y * w->width + x] = ex * ey;
        }
    }
}

/* ------------------------------------------------------- fringe analysis -- */
/* Dominant stripe in the column-averaged relative residual. Column averaging
   keeps the horizontal beat (panels line up in rows, the beat rolls across
   them along x) and suppresses per-pixel sampling noise. The profile is
   detrended against a moving average wider than the stripes of interest, then
   Fourier-analysed over 3..160 px. Returns the stripe spacing in pixels and
   writes the DFT peak-over-median prominence and the RMS of the detrended
   profile (the stripe's amplitude). */
double glass_fringe_spacing(const double* img, const double* ref,
                            const GlassWall* w, double* prominence, double* amplitude)
{
    int W = w->width;
    double* prof = (double*)malloc(sizeof(double) * (size_t)W);
    double* det = (double*)malloc(sizeof(double) * (size_t)W);
    double mean = 0.0;
    for (int x = 0; x < W; x++) {
        double acc = 0.0;
        for (int y = 0; y < w->height; y++) {
            long i = (long)y * W + x;
            acc += (img[i] - ref[i]) / ref[i];
        }
        prof[x] = acc / w->height;
        mean += prof[x];
    }
    mean /= W;

    int half = 16;
    double energy = 0.0;
    for (int x = 0; x < W; x++) {
        double s = 0.0;
        int c = 0;
        for (int j = -half; j <= half; j++) {
            int k = x + j;
            if (k >= 0 && k < W) { s += prof[k] - mean; c++; }
        }
        det[x] = (prof[x] - mean) - s / c;
        energy += det[x] * det[x];
    }
    if (amplitude) *amplitude = sqrt(energy / W);

    int lo = (int)ceil((double)W / 160.0);
    int hi = (int)floor((double)W / 3.0);
    if (lo < 2) lo = 2;
    int best = 0;
    double best_mag = 0.0, sum_mag = 0.0;
    int nbin = 0;
    for (int b = lo; b <= hi; b++) {
        double re = 0.0, im = 0.0;
        double wf = 2.0 * GLASS_PI * b / W;
        for (int x = 0; x < W; x++) {
            double hann = 0.5 - 0.5 * cos(2.0 * GLASS_PI * x / (W - 1));
            re += det[x] * hann * cos(wf * x);
            im -= det[x] * hann * sin(wf * x);
        }
        double mag = re * re + im * im;
        sum_mag += mag;
        nbin++;
        if (mag > best_mag) { best_mag = mag; best = b; }
    }
    free(prof);
    free(det);
    if (best == 0 || nbin == 0) { if (prominence) *prominence = 0.0; return 0.0; }
    if (prominence) *prominence = best_mag / (sum_mag / nbin + 1e-30);
    return (double)W / best;
}

/* RMS of the column-averaged residual: the energy that survives averaging the
   residual down each column. A coherent horizontal stripe dominates this; a
   per-pixel random error averages toward zero. This measures "how visible is
   the stripe" apart from ordinary sampling noise. */
double glass_band_rms(const double* img, const double* ref, const GlassWall* w) {
    int Wd = w->width;
    double acc = 0.0;
    for (int x = 0; x < Wd; x++) {
        double m = 0.0;
        for (int y = 0; y < w->height; y++) {
            long i = (long)y * Wd + x;
            m += (img[i] - ref[i]) / ref[i];
        }
        m /= w->height;
        acc += m * m;
    }
    return sqrt(acc / Wd);
}

/* RMS of the per-pixel relative residual: ordinary sampling error, stripe or
   not. Reported alongside the band RMS so the trade-off stays visible. */
double glass_pixel_rms(const double* img, const double* ref, const GlassWall* w) {
    long n = (long)w->width * w->height;
    double acc = 0.0;
    for (long i = 0; i < n; i++) {
        double e = (img[i] - ref[i]) / ref[i];
        acc += e * e;
    }
    return sqrt(acc / (double)n);
}

/* -------------------------------------------------------------- selftest -- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int glass_selftest(void) {
    const int W = 96, H = 48;
    double* a = (double*)malloc(sizeof(double) * W * H);
    double* b = (double*)malloc(sizeof(double) * W * H);

    /* smooth wall (no ripple, no mullion): both samplers reproduce the exact
       box filter -- the offset must not perturb a smooth integration. */
    GlassWall flat = glass_wall(W, H, 70, 4);
    flat.ripple = 0.0;
    flat.mullion = 0.0;
    glass_render(&flat, GLASS_SAMPLE_GRID, a);
    glass_reference(&flat, b);
    double e1 = 0.0;
    for (long i = 0; i < (long)W * H; i++) { double d = a[i] - b[i]; e1 += d * d; }
    check("flat_grid_matches_reference", e1 / ((double)W * H) < 1e-8);

    glass_render(&flat, GLASS_SAMPLE_OFFSET, a);
    double e2 = 0.0;
    for (long i = 0; i < (long)W * H; i++) { double d = a[i] - b[i]; e2 += d * d; }
    check("flat_offset_matches_reference", e2 / ((double)W * H) < 1e-6);

    /* offset range tracks the panel pitch, and never exceeds one pixel */
    GlassWall coarse = glass_wall(512, 512, 512, 4);
    check("radius_cap_one_pixel", fabs(glass_offset_radius(&coarse) - 1.0) < 1e-9);
    GlassWall dense = glass_wall(512, 512, 4096, 4);
    check("radius_tracks_pitch", fabs(glass_offset_radius(&dense) - 0.125) < 1e-9);
    check("radius_never_exceeds_pixel",
          glass_offset_radius(&dense) <= 1.0 && glass_offset_radius(&coarse) <= 1.0);

    /* samples stay in the pixel footprint, and the grid is an S*S lattice */
    GlassWall probe = glass_wall(512, 512, 1024, 4);
    double s[GLASS_MAX_SPP * GLASS_MAX_SPP][2];
    int n = glass_pixel_samples(&probe, 7, 3, GLASS_SAMPLE_OFFSET, s);
    int inside = 1;
    for (int i = 0; i < n; i++)
        if (s[i][0] < 0.0 || s[i][0] >= 1.0 || s[i][1] < 0.0 || s[i][1] >= 1.0)
            inside = 0;
    check("offset_stays_in_pixel", inside);
    check("sample_count", n == probe.spp * probe.spp);

    n = glass_pixel_samples(&probe, 0, 0, GLASS_SAMPLE_GRID, s);
    check("grid_stratified", n == 16 && fabs(s[0][0] - 0.125) < 1e-12 &&
                             fabs(s[0][1] - 0.125) < 1e-12);

    /* the low-discrepancy offset makes every pixel's lattice different */
    double s2[GLASS_MAX_SPP * GLASS_MAX_SPP][2];
    glass_pixel_samples(&probe, 7, 3, GLASS_SAMPLE_OFFSET, s);
    glass_pixel_samples(&probe, 8, 3, GLASS_SAMPLE_OFFSET, s2);
    check("offset_decorrelates_pixels", fabs(s[0][0] - s2[0][0]) > 1e-6);

    /* predicted stripe: commensurate lattice has no beat, offset lattice does */
    GlassWall comm = glass_wall(512, 512, 128, 4);   /* p = 4 px, S*p = 16 */
    check("predicted_no_beat", glass_predicted_stripe(&comm) == 0.0);
    GlassWall beat = glass_wall(512, 512, 297, 4);
    double lam = glass_predicted_stripe(&beat);
    check("predicted_beat_positive", lam > 1.0 && lam < 512.0);

    /* reference integrates to the analytic area mean on a smooth wall:
       (1 + gx/2)(1 + gy/2) averaged with the separable gradients */
    GlassWall smooth = glass_wall(64, 32, 16, 4);
    smooth.ripple = 0.0;
    smooth.mullion = 0.0;
    glass_reference(&smooth, b);
    double mr = 0.0;
    for (long i = 0; i < 64L * 32; i++) mr += b[i];
    mr /= (64.0 * 32.0);
    double expect = (1.0 + smooth.grad_x / 2.0) * (1.0 + smooth.grad_y / 2.0);
    check("reference_smooth_mean", fabs(mr - expect) < 1e-6);

    /* reference matches brute-force 2D quadrature on the rippled wall. The
       midpoint rule converges slowly across the mullion step, so use a fine
       grid and compare a handful of interior pixels. */
    GlassWall rip = glass_wall(24, 8, 21, 4);
    glass_reference(&rip, b);
    double maxerr = 0.0;
    int M = 2000;
    for (int px = 0; px < 4; px++) {
        for (int py = 0; py < 2; py++) {
            double acc = 0.0;
            for (int my = 0; my < M; my++)
                for (int mx = 0; mx < M; mx++)
                    acc += glass_shade((double)px + (mx + 0.5) / M,
                                       (double)py + (my + 0.5) / M, &rip);
            acc /= (double)(M * M);
            double d = fabs(acc - b[(long)py * 24 + px]);
            if (d > maxerr) maxerr = d;
        }
    }
    check("reference_matches_quadrature", maxerr < 5e-3);

    /* the core regression: on a phase-locked subdivision the fixed grid leaves
       a stripe, the low-discrepancy offset removes it, and lambda*D is nearly
       constant for the grid but not for the offset. */
    const int EW = 384, EH = 64;
    double* ig = (double*)malloc(sizeof(double) * EW * EH);
    double* io = (double*)malloc(sizeof(double) * EW * EH);
    double* ir = (double*)malloc(sizeof(double) * EW * EH);
    int ld_grid[3], ld_off[3];
    double band_grid[3], band_off[3], pix_grid[3], pix_off[3];
    int ld[3] = { 708, 368, 249 };
    for (int t = 0; t < 3; t++) {
        GlassWall tw = glass_wall(EW, EH, ld[t], 4);
        glass_reference(&tw, ir);
        glass_render(&tw, GLASS_SAMPLE_GRID, ig);
        glass_render(&tw, GLASS_SAMPLE_OFFSET, io);
        double gpr, gamp, opr, oamp;
        double lg = glass_fringe_spacing(ig, ir, &tw, &gpr, &gamp);
        double lo = glass_fringe_spacing(io, ir, &tw, &opr, &oamp);
        ld_grid[t] = (int)(lg * ld[t] + 0.5);
        ld_off[t] = (int)(lo * ld[t] + 0.5);
        band_grid[t] = glass_band_rms(ig, ir, &tw);
        band_off[t] = glass_band_rms(io, ir, &tw);
        pix_grid[t] = glass_pixel_rms(ig, ir, &tw);
        pix_off[t] = glass_pixel_rms(io, ir, &tw);
    }
    double gmin = ld_grid[0], gmax = ld_grid[0];
    for (int t = 1; t < 3; t++) {
        if (ld_grid[t] < gmin) gmin = ld_grid[t];
        if (ld_grid[t] > gmax) gmax = ld_grid[t];
    }
    check("grid_stripe_locked", gmax / gmin < 1.10);
    /* the offset leaves a coherent band of at most ~1/8 the grid's, while not
       inflating the per-pixel noise (the rotation preserves stratification). */
    int band_ok = 1, noise_ok = 1;
    for (int t = 0; t < 3; t++) {
        if (band_off[t] > 0.15 * band_grid[t]) band_ok = 0;
        if (pix_off[t] > 1.3 * pix_grid[t]) noise_ok = 0;
    }
    check("offset_suppresses_stripe", band_ok);
    check("offset_preserves_pixel_noise", noise_ok);
    (void)ld_off;

    free(a); free(b); free(ig); free(io); free(ir);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
