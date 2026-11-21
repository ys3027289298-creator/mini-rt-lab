/* dof: thin-lens depth of field — circle of confusion from aperture and
   focal length, lens-area sampling, and checks against analytic depth of
   field and the analytic CoC radius.

   All lengths are millimetres.  The lens is an ideal thin lens with focal
   length f and f-number N; the aperture radius is a = f / (2N).  Focus is
   at object distance s_f, so the sensor sits at the matching image plane.

     image distance   v(s) = f s / (s - f)
     CoC diameter     c(s) = f^2 |s - s_f| / (N s (s_f - f))
     CoC radius       r(s) = a |v(s) - v(s_f)| / v(s)
     hyperfocal       H    = f^2 / (N c_allow) + f
     near limit       s_n  = f^2 s_f / (f^2 + N c_allow (s_f - f))
     far limit        s_r  = f^2 s_f / (f^2 - N c_allow (s_f - f))

   The aperture is sampled uniformly by area: u1,u2 in [0,1] map to
     l = a sqrt(u1) (cos 2pi u2, sin 2pi u2).
   A point at distance s images at (0,0,-v(s)); the ray from lens point l
   through that image crosses the sensor plane z = -v(s_f) at
     l (1 - v(s_f)/v(s)),
   whose length is the sampled CoC radius.  That geometric construction is
   what `lens_sensor_hit` computes, independent of the closed-form c(s).
*/
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846

typedef struct { double f, N, s_f; } Lens;

static int finite_pos(double x) { return isfinite(x) && x > 0.0; }

/* Build a lens.  Rejects zero/negative/non-finite f-number, non-positive
   focal length, focus at or inside the focal length, and an aperture radius
   that overflows.  On failure the struct is left invalid (NaN fields). */
int lens_init(Lens* L, double focal, double fnumber, double focus) {
    if (!L) return -1;
    L->f = L->N = L->s_f = NAN;
    if (!finite_pos(focal) || !finite_pos(fnumber)) return -1;
    if (!isfinite(focus) || !(focus > focal)) return -1;
    double a = focal / (2.0 * fnumber);
    if (!isfinite(a) || !(a > 0.0)) return -1;
    L->f = focal;
    L->N = fnumber;
    L->s_f = focus;
    return 0;
}

int lens_valid(const Lens* L) {
    return L && finite_pos(L->f) && finite_pos(L->N) &&
           isfinite(L->s_f) && L->s_f > L->f;
}

double lens_aperture_radius(const Lens* L) {
    return lens_valid(L) ? L->f / (2.0 * L->N) : NAN;
}

/* Image-side distance for an object at s.  Returns f in the limit s -> inf,
   and NaN where no real image exists (s <= f). */
double lens_image_distance(double f, double s) {
    if (!finite_pos(f)) return NAN;
    if (isinf(s)) return f;
    if (!(s > f)) return NAN;
    return f * s / (s - f);
}

double lens_coc_diameter(const Lens* L, double s) {
    if (!lens_valid(L)) return NAN;
    if (isinf(s))
        return L->f * L->f / (L->N * (L->s_f - L->f));
    if (!(s > L->f)) return NAN;            /* no real image to blur */
    return L->f * L->f * fabs(s - L->s_f) /
           (L->N * s * (L->s_f - L->f));
}

double lens_coc_radius(const Lens* L, double s) {
    return 0.5 * lens_coc_diameter(L, s);
}

double lens_hyperfocal(const Lens* L, double coc_allow) {
    if (!lens_valid(L) || !finite_pos(coc_allow)) return NAN;
    return L->f * L->f / (L->N * coc_allow) + L->f;
}

/* Analytic depth of field for an allowed CoC diameter.  *far is INFINITY
   when the far limit lies beyond the hyperfocal distance. */
int lens_dof(const Lens* L, double coc_allow, double* near, double* far) {
    if (!lens_valid(L) || !finite_pos(coc_allow)) return -1;
    double d = L->s_f - L->f;
    double num = L->f * L->f * L->s_f;
    double den_n = L->f * L->f + L->N * coc_allow * d;
    double den_f = L->f * L->f - L->N * coc_allow * d;
    if (near) *near = num / den_n;
    if (far) *far = (den_f > 0.0) ? num / den_f : INFINITY;
    return 0;
}

/* Uniform-by-area sample of the aperture disk.  u1 = 1 lands on the rim. */
int lens_sample_disk(const Lens* L, double u1, double u2, double* lx, double* ly) {
    if (!lens_valid(L) || !isfinite(u1) || !isfinite(u2)) return -1;
    if (u1 < 0.0) u1 = 0.0;
    if (u1 > 1.0) u1 = 1.0;
    if (u2 < 0.0) u2 = 0.0;
    if (u2 > 1.0) u2 = 1.0;
    double a = lens_aperture_radius(L);
    double r = a * sqrt(u1);
    double t = 2.0 * PI * u2;
    if (lx) *lx = r * cos(t);
    if (ly) *ly = r * sin(t);
    return 0;
}

/* Ray from lens point (lx,ly,0) through the image of a point at distance s,
   intersected with the sensor plane z = -v(s_f).  Returns -1 if the lens or
   the object distance has no real image. */
int lens_sensor_hit(const Lens* L, double lx, double ly, double s,
                    double* hx, double* hy) {
    if (!lens_valid(L)) return -1;
    double v = lens_image_distance(L->f, s);
    double vf = lens_image_distance(L->f, L->s_f);
    if (!isfinite(v) || !isfinite(vf) || !(v > 0.0) || !(vf > 0.0))
        return -1;
    double t = vf / v;
    if (hx) *hx = lx * (1.0 - t);
    if (hy) *hy = ly * (1.0 - t);
    return 0;
}

double lens_spot_radius(const Lens* L, double lx, double ly, double s) {
    double hx, hy;
    if (lens_sensor_hit(L, lx, ly, s, &hx, &hy) != 0) return NAN;
    return hypot(hx, hy);
}

/* ------------------------------------------------------------------ */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static double rim_radius(const Lens* L, double s) {
    double lx, ly;
    lens_sample_disk(L, 1.0, 0.3, &lx, &ly);
    return lens_spot_radius(L, lx, ly, s);
}

int main(void) {
    Lens L;
    /* -- construction and rejection ------------------------------- */
    check("init_ok", lens_init(&L, 50.0, 1.8, 2000.0) == 0 && lens_valid(&L));
    check("aperture_radius", fabs(lens_aperture_radius(&L) - 50.0 / 3.6) < 1e-12);

    Lens bad;
    double bx, by;
    check("reject_zero_aperture",
          lens_init(&bad, 50.0, 0.0, 2000.0) != 0 && !lens_valid(&bad) &&
          isnan(lens_coc_diameter(&bad, 1000.0)) &&
          lens_sample_disk(&bad, 0.5, 0.5, &bx, &by) != 0);
    check("reject_negative_aperture",
          lens_init(&bad, 50.0, -2.0, 2000.0) != 0 && !lens_valid(&bad));
    check("reject_bad_focal",
          lens_init(&bad, 0.0, 2.0, 2000.0) != 0 &&
          lens_init(&bad, -50.0, 2.0, 2000.0) != 0);
    check("reject_focus_inside_focal",
          lens_init(&bad, 50.0, 2.0, 50.0) != 0 &&
          lens_init(&bad, 50.0, 2.0, 25.0) != 0);
    check("reject_nonfinite",
          lens_init(&bad, NAN, 2.0, 2000.0) != 0 &&
          lens_init(&bad, 50.0, INFINITY, 2000.0) != 0 &&
          lens_init(&bad, 50.0, 2.0, INFINITY) != 0);
    check("reject_overflow_aperture",
          lens_init(&bad, 50.0, 1e-308, 2000.0) != 0 && !lens_valid(&bad));

    Lens a8;
    lens_init(&a8, 50.0, 8.0, 2000.0);
    check("coc_known_value",
          fabs(lens_coc_diameter(&a8, 1000.0) - 0.16025641025641025) < 1e-12);
    check("coc_at_infinity",
          fabs(lens_coc_diameter(&a8, INFINITY) - 0.16025641025641025) < 1e-12);

    /* -- focus plane must be exactly sharp ------------------------- */
    int sharp = 1;
    for (int k = 0; k < 16; k++) {
        double lx, ly, hx, hy;
        lens_sample_disk(&L, 1.0, (double)k / 16.0, &lx, &ly);
        if (lens_sensor_hit(&L, lx, ly, L.s_f, &hx, &hy) != 0 ||
            hx != 0.0 || hy != 0.0)
            sharp = 0;
    }
    check("focus_exactly_sharp",
          sharp && lens_coc_diameter(&L, L.s_f) == 0.0);

    /* -- lens disk is sampled uniformly by area -------------------- */
    double a = lens_aperture_radius(&L);
    double sum_r = 0.0, sum_r2 = 0.0, rmax = 0.0;
    int inside = 1, m = 64;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double lx, ly;
            double u1 = ((double)i + 0.5) / m;
            double u2 = ((double)j + 0.5) / m;
            if (lens_sample_disk(&L, u1, u2, &lx, &ly) != 0) {
                inside = 0;
                continue;
            }
            double r = hypot(lx, ly);
            if (r > a) inside = 0;
            sum_r += r;
            sum_r2 += r * r;
            if (r > rmax) rmax = r;
        }
    }
    double ns = (double)m * m;
    check("lens_samples_inside", inside && rmax <= a);
    /* area-uniform disk: E[r] = 2a/3, E[r^2] = a^2/2.  The midpoint rule
       underestimates E[sqrt(U)] by O(1/m^2), so allow 2e-4 at m = 64. */
    check("lens_area_uniform_mean_r",
          fabs(sum_r / ns / a - 2.0 / 3.0) < 2e-4);
    check("lens_area_uniform_mean_r2",
          fabs(sum_r2 / ns / (a * a) - 0.5) < 1e-12);

    double rlx, rly;
    lens_sample_disk(&L, 1.0, 0.25, &rlx, &rly);
    check("lens_rim_exact", fabs(hypot(rlx, rly) - a) < 1e-12);

    /* -- sampled CoC vs analytic, with an explicit error bound ----- */
    const double allow = 0.03;
    const double cfgs[3][3] = {
        { 50.0, 1.8, 2000.0 },
        { 85.0, 2.8, 1200.0 },
        { 24.0, 8.0,  800.0 },
    };
    const double u1s[4] = { 1.0, 0.75, 0.5, 0.25 };
    double err_max = 0.0, bound_max = 0.0, worst = 0.0;
    int sampling_ok = 1, near_coc_ok = 1, far_coc_ok = 1;
    int order_ok = 1, inside_ok = 1, outside_ok = 1;
    int near_sampled_ok = 1, far_sampled_ok = 1;

    for (int ci = 0; ci < 3; ci++) {
        Lens M;
        if (lens_init(&M, cfgs[ci][0], cfgs[ci][1], cfgs[ci][2]) != 0) {
            sampling_ok = 0;
            continue;
        }
        double aM = lens_aperture_radius(&M);
        double svals[] = {
            0.5 * M.s_f, 0.8 * M.s_f, M.s_f, 1.25 * M.s_f,
            2.0 * M.s_f, 10.0 * M.s_f, INFINITY,
        };
        int nsv = (int)(sizeof(svals) / sizeof(svals[0]));
        for (int si = 0; si < nsv; si++) {
            double s = svals[si];
            double rana = lens_coc_radius(&M, s);
            if (!isfinite(rana)) {
                sampling_ok = 0;
                continue;
            }
            for (int k = 0; k < 16; k++) {
                for (int q = 0; q < 4; q++) {
                    double lx, ly;
                    lens_sample_disk(&M, u1s[q],
                                     (double)k / 16.0 + 0.01, &lx, &ly);
                    double lr = hypot(lx, ly);
                    double want = rana * lr / aM;
                    double got = lens_spot_radius(&M, lx, ly, s);
                    double err = fabs(got - want);
                    double bound = 1e-9 + 1e-12 * want;
                    if (!(err <= bound)) sampling_ok = 0;
                    if (err > err_max) err_max = err;
                    if (bound > bound_max) bound_max = bound;
                    if (err / bound > worst) worst = err / bound;
                }
            }
        }

        double nr, fr;
        if (lens_dof(&M, allow, &nr, &fr) != 0) {
            order_ok = 0;
            continue;
        }
        if (!(nr > M.f && nr < M.s_f && fr > M.s_f)) order_ok = 0;
        if (fabs(lens_coc_diameter(&M, nr) - allow) > 1e-9) near_coc_ok = 0;
        if (!isinf(fr) &&
            fabs(lens_coc_diameter(&M, fr) - allow) > 1e-9) far_coc_ok = 0;
        if (lens_coc_diameter(&M, 0.5 * (nr + fr)) > allow + 1e-9) inside_ok = 0;
        if (lens_coc_diameter(&M, nr * 0.99) <= allow) outside_ok = 0;
        if (!isinf(fr) &&
            lens_coc_diameter(&M, fr * 1.01) <= allow) outside_ok = 0;
        if (fabs(rim_radius(&M, nr) - 0.5 * allow) > 1e-9) near_sampled_ok = 0;
        if (!isinf(fr) &&
            fabs(rim_radius(&M, fr) - 0.5 * allow) > 1e-9) far_sampled_ok = 0;
    }

    check("coc_samples_match_analytic", sampling_ok);
    check("coc_err_within_bound", worst <= 1.0);
    printf("coc_err_bound=%.3e\n", bound_max);
    printf("coc_err_max=%.3e\n", err_max);
    check("dof_near_coc", near_coc_ok);
    check("dof_far_coc", far_coc_ok);
    check("dof_order", order_ok);
    check("dof_inside_acceptably_sharp", inside_ok);
    check("dof_outside_blurred", outside_ok);
    check("dof_near_sampled_coc", near_sampled_ok);
    check("dof_far_sampled_coc", far_sampled_ok);

    /* -- aperture dependence --------------------------------------- */
    Lens n18, n8;
    lens_init(&n18, 50.0, 1.8, 2000.0);
    lens_init(&n8, 50.0, 8.0, 2000.0);
    check("stopping_down_reduces_coc",
          lens_coc_diameter(&n8, 1000.0) < lens_coc_diameter(&n18, 1000.0));
    double nr1, fr1, nr8, fr8;
    lens_dof(&n18, allow, &nr1, &fr1);
    lens_dof(&n8, allow, &nr8, &fr8);
    check("stopping_down_deepens_dof",
          (fr8 - nr8) > (fr1 - nr1) && nr8 < nr1 && fr8 > fr1);

    /* -- hyperfocal / infinite far limit --------------------------- */
    double H = lens_hyperfocal(&L, allow);
    check("hyperfocal_beyond_focus", isfinite(H) && H > L.s_f);
    Lens hf;
    lens_init(&hf, 50.0, 1.8, 1.01 * H);
    double hn, hfar;
    check("hyperfocal_far_infinite",
          lens_dof(&hf, allow, &hn, &hfar) == 0 && isinf(hfar));
    Lens hf2;
    lens_init(&hf2, 50.0, 1.8, H);
    check("hyperfocal_coc_at_infinity",
          fabs(lens_coc_diameter(&hf2, INFINITY) - allow) <= 1e-9);

    /* -- other domain checks --------------------------------------- */
    check("reject_zero_coc_for_dof",
          lens_dof(&L, 0.0, NULL, NULL) != 0 &&
          lens_dof(&L, -0.03, NULL, NULL) != 0);
    check("reject_object_inside_focal",
          isnan(lens_coc_diameter(&L, 25.0)) &&
          isnan(lens_spot_radius(&L, 1.0, 0.0, 25.0)));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
