
/* penumbra: analytic penumbra-width verification.

   A rectangular area light (area A, edge length H, across-edge extent L=A/H)
   casts a shadow of a straight-edged occluder onto a receiver line. In the
   cross section perpendicular to the occluder edge the normalized irradiance
   (visible fraction of the emitter) has a closed form:

     umbra edge   x_u = b + (b - a) * d_r / d_l
     penumbra end x_p = b + (b + a) * d_r / d_l      (a = L/2)
     penumbra width W   = x_p - x_u = L * d_r / d_l = (A / H) * d_r / d_l

   With a real umbra f(x) is 0 for |x| <= x_u, (|x| - x_u) / W inside the
   ramp, and 1 past x_p. exact_f() instead uses the general blocked-interval
   overlap in light coordinates, which also covers tangent geometry (x_u = 0)
   and the no-umbra case (x_u < 0, the profile never reaches zero).

   The Monte-Carlo "render" below estimates the same visible fraction, so the
   comparison is like-for-like. The analytic straight-edge profile is C0 with
   a finite slope kink (piecewise linear), not C1: the boundary checks verify
   value continuity (no hard step) and that the measured kink matches the
   analytic ramp slope within sampling noise. A zero-area light degenerates
   to a point light: W must be zero and the profile a binary hard step. */
#include <stdio.h>
#include <math.h>

#define N_SAMPLES 8192          /* MC light samples per receiver point       */
#define N_SCAN_MAX 2048         /* receiver scan points (grid adapts/scene)  */
#define PLATEAU_CELLS 16        /* lit-plateau cells kept past penumbra edge */
#define LEVEL_FRAC_LO 0.2       /* ramp-relative lower crossing level        */
#define LEVEL_FRAC_HI 0.8       /* ramp-relative upper crossing level        */

typedef struct {
    double area;     /* emitter area A = L * H                       */
    double edge_h;   /* emitter edge length H                        */
    double occl_b;   /* occluder half width b                        */
    double d_l;      /* light plane -> occluder distance             */
    double d_r;      /* occluder -> receiver distance                */
} Scene;

static unsigned long long rng_state;
static double randu(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 33) / (double)(1ULL << 31);
}

static double light_width(const Scene* s) { return s->area / s->edge_h; }

/* analytic shadow boundaries at the receiver */
static void edges(const Scene* s, double* xu, double* xp, double* w) {
    double a = light_width(s) * 0.5;
    double b = s->occl_b, dl = s->d_l, dr = s->d_r;
    *xu = b + (b - a) * dr / dl;
    *xp = b + (b + a) * dr / dl;
    *w  = light_width(s) * dr / dl;
}

/* ray from light point (l, 0) to receiver point (xr, D): x at occluder plane */
static int visible(double xr, double l, const Scene* s) {
    double D = s->d_l + s->d_r;
    double x_hit = l + (xr - l) * s->d_l / D;
    return fabs(x_hit) > s->occl_b;
}

/* exact normalized irradiance: fraction of the emitter interval [-a, a]
   whose rays to receiver point xr clear the occluder [-b, b] at z = d_l.
   A ray through light point l hits the occluder plane at
     x_hit = l + (xr - l) * d_l / D,
   so the blocked light points form an interval [lo_b, hi_b] in light
   coordinates. f = 1 - |intersection with [-a, a]| / (2a). The overlap form
   covers umbra, partial penumbra and the no-umbra geometry in one expression. */
static double exact_f(double xr, const Scene* s) {
    double a = light_width(s) * 0.5;
    if (a <= 0.0) return visible(xr, 0.0, s) ? 1.0 : 0.0;
    double D = s->d_l + s->d_r;
    /* invert x_hit = c: l = (c*D - xr*d_l) / d_r; blocked light points are
       those mapping into [-b, b], i.e. l in [l(-b), l(+b)] */
    double lo_b = (-s->occl_b * D - xr * s->d_l) / s->d_r;
    double hi_b = ( s->occl_b * D - xr * s->d_l) / s->d_r;
    double cov = fmax(0.0, fmin(hi_b, a) - fmax(lo_b, -a));
    return fmin(1.0, fmax(0.0, 1.0 - cov / (2.0 * a)));
}

/* Monte-Carlo estimate: stratified samples across the emitter extent.
   A zero-width emitter is sampled as a single point (hard binary result). */
static double mc_f(double xr, const Scene* s) {
    double half = light_width(s) * 0.5;
    if (half <= 0.0) return visible(xr, 0.0, s) ? 1.0 : 0.0;
    int hit = 0;
    for (int j = 0; j < N_SAMPLES; j++) {
        double u = ((double)j + randu()) / (double)N_SAMPLES;
        double l = -half + 2.0 * half * u;
        hit += visible(xr, l, s);
    }
    return (double)hit / (double)N_SAMPLES;
}

/* x where the right-side (rising) ramp crosses level t, by linear
   interpolation inside the bracketing cell. NaN if no crossing exists. */
static double crossing(const double* x, const double* f, int n, double t) {
    for (int i = 0; i + 1 < n; i++) {
        if ((f[i] - t) * (f[i + 1] - t) <= 0.0 && f[i] != f[i + 1])
            return x[i] + (t - f[i]) * (x[i + 1] - x[i]) / (f[i + 1] - f[i]);
    }
    return NAN;
}

/* penumbra width from two level crossings of the right rising ramp. Levels
   are taken relative to the ramp's own range [f_lo, 1] — f_lo is the profile
   value at x = 0 (0 with an umbra, > 0 in no-umbra geometry, where a flat
   bottom plateau would make a fixed 0.2 crossing ill-defined). The ramp
   slope is 1/W, so the level span divided by (t_hi - t_lo) recovers W.
   The profile is symmetric about x = 0, so one ramp suffices. */
static double estimate_width(const double* x, const double* f, int n,
                             double t_lo, double t_hi) {
    return (crossing(x, f, n, t_hi) - crossing(x, f, n, t_lo))
           / (t_hi - t_lo);
}

static void linfit(const double* x, const double* y, int n, double* m, double* c) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 0; i < n; i++) {
        sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i];
    }
    double den = (double)n * sxx - sx * sx;
    *m = ((double)n * sxy - sx * sy) / den;
    *c = (sy - *m * sx) / (double)n;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static void run(const char* tag, const Scene* s) {
    double xu, xp, w;
    edges(s, &xu, &xp, &w);
    /* right-half scan — the profile is symmetric about x = 0. Grid targets
       ~24 cells across the analytic penumbra width and ends several plateau
       cells past the outer edge so no crossing can fall off the grid. */
    double edge = fmax(fabs(xp), fabs(xu));
    double h_ref = w > 0.0 ? w / 24.0 : fmax(edge, 1.0) / 40.0;
    double xmax = edge + PLATEAU_CELLS * h_ref;
    int n = (int)ceil(xmax / h_ref) + 1;
    if (n > N_SCAN_MAX) n = N_SCAN_MAX;
    double h = xmax / (double)(n - 1);

    double xs[N_SCAN_MAX], fe[N_SCAN_MAX], fm[N_SCAN_MAX];
    rng_state = 0x123456789abcdef0ULL;
    for (int i = 0; i < n; i++) {
        xs[i] = h * (double)i;
        fe[i] = exact_f(xs[i], s);
        fm[i] = mc_f(xs[i], s);
    }

    /* ramp-relative levels: [f_lo, 1] mapped to 20% / 80% */
    double f_lo = exact_f(0.0, s);
    double t_lo = f_lo + LEVEL_FRAC_LO * (1.0 - f_lo);
    double t_hi = f_lo + LEVEL_FRAC_HI * (1.0 - f_lo);
    double span = t_hi - t_lo;
    double w_hat = estimate_width(xs, fm, n, t_lo, t_hi);
    double sigma = 1.0 / (2.0 * sqrt((double)N_SAMPLES));      /* max Bernoulli std */
    /* crossing position std <= sigma / slope with slope = 1/W; the width
       combines two crossings (sqrt 2), a 4-sigma noise allowance, plus
       1.5 grid cells of linear-interpolation bracket error */
    double bound_rel = 4.0 * sqrt(2.0) * sigma / span
                     + 1.5 * h / fmax(w, h);
    double bound_abs = bound_rel * fmax(w, h);

    printf("# --- %s ---\n", tag);
    printf("# light_width=%.6f umbra_edge=%.6f penumbra_edge=%.6f\n",
           light_width(s), xu, xp);
    printf("# analytic_width=%.6f measured_width=%.6f abs_err=%.6f rel_err=%.6f\n",
           w, w_hat, fabs(w_hat - w), fabs(w_hat - w) / fmax(w, h));
    printf("# mc_sigma=%.6f err_bound_abs=%.6f err_bound_rel=%.6f grid_h=%.6f n=%d\n",
           sigma, bound_abs, bound_rel, h, n);

    double max_err = 0.0, max_jump = 0.0;
    for (int i = 0; i < n; i++)
        max_err = fmax(max_err, fabs(fm[i] - fe[i]));
    for (int i = 0; i + 1 < n; i++)
        max_jump = fmax(max_jump, fabs(fm[i + 1] - fm[i]));

    if (w > 0.0) {
        check("profile_matches_analytic", max_err <= 4.0 * sigma);
        check("width_within_bound", fabs(w_hat - w) <= bound_abs);
        check("boundary_no_jump",
              max_jump <= h / w + 4.0 * sqrt(2.0) * sigma);

        if (xu > 2.0 * h) {
            /* real umbra: ramp joins the zero plateau at xu and 1 at xp.
               least-squares slopes of ramp interior vs lit plateau */
            double rx[64], ry[64], px[32], py[32];
            int nr = 0, np = 0;
            for (int i = 0; i < n; i++) {
                double r = xs[i];
                if (r > xu + 0.15 * w && r < xp - 0.15 * w && nr < 64) {
                    rx[nr] = r; ry[nr] = fm[i]; nr++;
                }
                if (r > xp + 3.0 * h && r < xp + 12.0 * h && np < 32) {
                    px[np] = r; py[np] = fm[i]; np++;
                }
            }
            double mr, cr, mp, cp;
            linfit(rx, ry, nr, &mr, &cr);
            linfit(px, py, np, &mp, &cp);
            double kink = fabs(mr - mp);                 /* analytic kink = 1/W */
            double join_u = fabs(mr * xu + cr);          /* ramp value at umbra edge ~ 0 */
            double join_p = fabs((mr * xp + cr) - 1.0);  /* ramp value at lit edge ~ 1 */
            printf("# ramp_slope=%.6f plateau_slope=%.6f kink=%.6f analytic_kink=%.6f\n",
                   mr, mp, kink, 1.0 / w);
            check("ramp_fit_has_points", nr >= 4 && np >= 4);
            check("boundary_value_continuous",
                  join_u <= 4.0 * sigma && join_p <= 4.0 * sigma);
            check("boundary_kink_bounded",
                  kink < 3.0 / w && fabs(kink - 1.0 / w) <= 0.15);
        } else if (xu > -2.0 * h) {
            /* tangent geometry: umbra edge sits within two cells of x = 0,
               the exact center value is ~0; no plateau-join check applies */
            double f0 = exact_f(0.0, s);
            printf("# tangent center_exact=%.6f center_mc=%.6f\n", f0, fm[0]);
            check("tangent_center_matches", fabs(fm[0] - f0) <= 4.0 * sigma);
        } else {
            /* no umbra: the center stays lit, the profile bottoms out above 0 */
            double f0 = exact_f(0.0, s);
            printf("# no_umbra center_exact=%.6f center_mc=%.6f\n", f0, fm[0]);
            check("no_umbra_center_lit",
                  f0 > 0.0 && fabs(fm[0] - f0) <= 4.0 * sigma);
        }
    } else {
        /* point light: analytic penumbra is exactly zero; on a finite grid
           both level crossings fall inside the same cell, so the measured
           width can only reach one grid spacing (plus floating-point slack) */
        int binary = 1;
        for (int i = 0; i < n; i++)
            binary = binary && (fm[i] == 0.0 || fm[i] == 1.0);
        check("zero_area_analytic_width_zero", w == 0.0 && xu == xp);
        check("zero_area_measured_width_zero", w_hat <= h * (1.0 + 1e-12));
        check("zero_area_profile_binary", binary);
        check("point_light_hard_edge", max_jump >= 0.9);
    }
}

int main(void) {
    /* area light: A = 4.0, edge H = 4.0 -> across-edge width L = 1.0 */
    Scene area = { 4.0, 4.0, 1.2, 2.0, 3.0 };
    run("area_light", &area);

    /* wide light vs small occluder: L = 2.0 > 2b, umbra edge goes negative,
       i.e. no umbra exists — the overlap formula must still hold */
    Scene no_umbra = { 8.0, 4.0, 0.5, 2.0, 3.0 };
    run("wide_light_no_umbra", &no_umbra);

    /* tangent geometry: umbra edge lands exactly on x = 0, center hits zero
       without an umbra plateau */
    Scene tangent = { 8.0, 4.0, 0.6, 2.0, 3.0 };
    run("tangent_geometry", &tangent);

    /* same geometry, emitter collapsed to a point: penumbra must vanish */
    Scene point = { 0.0, 4.0, 1.2, 2.0, 3.0 };
    run("point_light", &point);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
