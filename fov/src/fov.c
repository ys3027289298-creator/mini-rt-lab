
/* fov: match resolution to field of view via pixel footprint.
   Per-pixel footprint on the plane at distance d, sampling density
   derived from it, and an upper bound on the geometric error of the
   linear (local) footprint approximation vs the analytic tan-difference. */
#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979323846
/* fov at or beyond this is flagged degenerate: tan(fov/2) grows without
   bound and edge-pixel footprints explode. */
#define FOV_DEGEN (170.0 * PI / 180.0)

typedef struct {
    int ok;          /* fov usable at all */
    int degenerate;  /* usable but footprint model breaks down at edges */
    const char* note;
} FovStatus;

static FovStatus fov_validate(double fov) {
    FovStatus s = { 1, 0, 0 };
    if (!(fov > 0.0)) { /* rejects zero, negative, NaN */
        s.ok = 0; s.note = "fov must be positive";
        return s;
    }
    if (fov >= PI) {
        s.ok = 0; s.note = "fov must be < pi (180 deg)";
        return s;
    }
    if (fov >= FOV_DEGEN) {
        s.degenerate = 1;
        s.note = "fov near 180 deg: tan(fov/2) diverges, edge-pixel footprint "
                 "is unbounded and no uniform sampling density can match it; "
                 "cap the effective fov or reparameterize (fisheye/spherical) "
                 "instead of a flat sensor";
    }
    return s;
}

/* ray angle at fractional position edge in [0,1] within pixel i of h */
static double pixel_angle(double fov, int h, int i, double edge) {
    return (2.0 * ((double)i + edge) / (double)h - 1.0) * fov * 0.5;
}

/* analytic footprint: exact tan-difference across the pixel */
static double footprint_exact(double fov, double dist, int h, int i) {
    double a0 = pixel_angle(fov, h, i, 0.0);
    double a1 = pixel_angle(fov, h, i, 1.0);
    return dist * (tan(a1) - tan(a0));
}

/* linear approximation at the pixel center: d * dphi / cos^2(phi) */
static double footprint_approx(double fov, double dist, int h, int i) {
    double phi = pixel_angle(fov, h, i, 0.5);
    double dphi = fov / (double)h;
    double c = cos(phi);
    return dist * dphi / (c * c);
}

/* upper bound on |exact - approx| for pixel i.
   tan(phi+d) - tan(phi-d) = 2 d sec^2(phi) + (d^3/6)(f'''(x1)+f'''(x2)),
   f = tan, f'''(x) = 2 sec^2(x) (2 tan^2(x) + sec^2(x)) >= 0, growing in |x|,
   so |error| <= d^3/3 * max|f'''| over the pixel's angular interval. */
static double footprint_error_bound(double fov, double dist, int h, int i) {
    double half = fov / (2.0 * (double)h);
    double phi = pixel_angle(fov, h, i, 0.5);
    double lo = fabs(phi - half), hi = fabs(phi + half);
    double xm = lo > hi ? lo : hi;
    double c = cos(xm), sec2 = 1.0 / (c * c), t = tan(xm);
    double f3 = 2.0 * sec2 * (2.0 * t * t + sec2);
    return dist * half * half * half / 3.0 * f3;
}

/* samples per unit area to resolve features at the footprint scale */
static double sampling_density(double footprint) {
    return 1.0 / (footprint * footprint);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    check("reject_zero", fov_validate(0.0).ok == 0);
    check("reject_negative", fov_validate(-0.5).ok == 0);
    check("reject_nan", fov_validate(0.0 / 0.0).ok == 0);
    check("reject_full", fov_validate(PI).ok == 0);

    FovStatus normal = fov_validate(PI / 3.0);
    check("accept_normal", normal.ok == 1 && normal.degenerate == 0);

    FovStatus huge = fov_validate(179.0 * PI / 180.0);
    check("degenerate_huge", huge.ok == 1 && huge.degenerate == 1);
    check("degenerate_note", huge.note != 0);

    double fov = PI / 3.0, dist = 2.0;
    int h = 64;

    /* geometric error of the linear approximation, against analytic */
    double err_center = fabs(footprint_exact(fov, dist, h, 32) -
                             footprint_approx(fov, dist, h, 32));
    check("center_close", err_center < 1e-6);

    int bound_ok = 1;
    for (int i = 0; i < h; i++) {
        double err = fabs(footprint_exact(fov, dist, h, i) -
                          footprint_approx(fov, dist, h, i));
        if (err > footprint_error_bound(fov, dist, h, i)) bound_ok = 0;
    }
    check("bound_holds_all", bound_ok);

    check("bound_shrinks", footprint_error_bound(fov, dist, 128, 64) <
                           footprint_error_bound(fov, dist, 64, 32));

    check("symmetric", fabs(footprint_exact(fov, dist, h, 0) -
                            footprint_exact(fov, dist, h, h - 1)) < 1e-12);

    check("edge_wider", footprint_exact(fov, dist, h, h - 1) >
                        footprint_exact(fov, dist, h, h / 2));

    check("density_inverse_square",
          fabs(sampling_density(0.5) - 4.0 * sampling_density(1.0)) < 1e-12);

    /* matching resolution to fov: at fixed h the footprint grows with
       fov, so samples per pixel must rise as (fp_wide/fp_narrow)^2 to
       keep the same world-space feature resolved */
    double fp_narrow = footprint_exact(PI / 6.0, dist, h, h / 2);
    double fp_wide = footprint_exact(PI / 2.0, dist, h, h / 2);
    check("footprint_grows_with_fov", fp_wide > fp_narrow);
    double ratio = fp_wide / fp_narrow;
    check("samples_scale_quadratic",
          sampling_density(fp_narrow) / sampling_density(fp_wide) >
          0.99 * ratio * ratio);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
