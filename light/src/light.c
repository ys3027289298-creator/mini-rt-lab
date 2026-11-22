
/* light: cone-light sampling with distance solid-angle correction.

   A disk area light of radius R seen from a shading point at distance d
   subtends exactly the same solid angle as a cone of half-angle
   alpha_exact = atan(R/d):

       Omega_disk(d) = 2*pi*(1 - 1/sqrt(1 + (R/d)^2))
                     = Omega_cone(atan(R/d)).

   The cone approximation fits one half-angle alpha_ref at a reference
   distance d_ref and corrects the solid angle with distance as 1/d^2:

       Omega_corr(d) = Omega_cone(alpha_ref) * (d_ref/d)^2.

   The corrected cone is exact at d = d_ref and matches the disk far
   field; its error is measured by q = Omega_corr/Omega_disk and bounded
   for soft-shadow width below.  A zero half-angle degenerates to a
   point (delta) light; angles outside the supported range are reported
   through LightStatus.
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#define LIGHT_PI 3.14159265358979323846
#define LIGHT_TWO_PI (2.0 * LIGHT_PI)
/* Largest half-angle for which the width-error bound stays finite.
   The closed form below loses validity near 40.3 deg; keep a margin. */
#define LIGHT_ALPHA_BOUND_MAX (35.0 * LIGHT_PI / 180.0)
#define LIGHT_ALPHA_HARD_MAX (0.5 * LIGHT_PI)

typedef struct { double x, y, z; } Vec3;

typedef enum {
    LIGHT_OK = 0,
    LIGHT_ERR_NULL,
    LIGHT_ERR_AXIS,
    LIGHT_ERR_RADIUS,
    LIGHT_ERR_DISTANCE,
    LIGHT_ERR_ANGLE_RANGE,
    LIGHT_ERR_ANGLE_BOUND,
    LIGHT_ERR_SOLID_ANGLE,
    LIGHT_ERR_UV
} LightStatus;

typedef struct {
    Vec3 axis;          /* unit direction from shading point to light center */
    double radius;      /* disk radius R (exact area reference) */
    double dist_ref;    /* reference distance the cone is fitted at */
    double alpha_ref;   /* cone half-angle at dist_ref */
    int point;          /* 1 when alpha_ref == 0: degenerate point light */
    LightStatus status;
} ConeLight;

static double dot3(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static double len3(Vec3 a) { return sqrt(dot3(a, a)); }
static Vec3 scale3(Vec3 a, double s) { Vec3 r = { a.x * s, a.y * s, a.z * s }; return r; }
static Vec3 add3(Vec3 a, Vec3 b) { Vec3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }

const char* light_status_str(LightStatus s) {
    switch (s) {
    case LIGHT_OK: return "OK";
    case LIGHT_ERR_NULL: return "ERR_NULL";
    case LIGHT_ERR_AXIS: return "ERR_AXIS";
    case LIGHT_ERR_RADIUS: return "ERR_RADIUS";
    case LIGHT_ERR_DISTANCE: return "ERR_DISTANCE";
    case LIGHT_ERR_ANGLE_RANGE: return "ERR_ANGLE_RANGE";
    case LIGHT_ERR_ANGLE_BOUND: return "ERR_ANGLE_BOUND";
    case LIGHT_ERR_SOLID_ANGLE: return "ERR_SOLID_ANGLE";
    case LIGHT_ERR_UV: return "ERR_UV";
    }
    return "ERR_UNKNOWN";
}

/* Solid angle of a cone with half-angle alpha, zero at alpha = 0.
   A full sphere (alpha = pi) is 4*pi; half-space (alpha = pi/2) is 2*pi. */
static double cone_solid_angle(double alpha) {
    if (alpha <= 0.0) return 0.0;
    if (alpha >= LIGHT_PI) return 2.0 * LIGHT_TWO_PI;
    return LIGHT_TWO_PI * (1.0 - cos(alpha));
}

/* Exact solid angle of a disk of radius R at distance d, seen on axis. */
static double disk_solid_angle(double R, double d) {
    double t = R / d;
    return LIGHT_TWO_PI * (1.0 - 1.0 / sqrt(1.0 + t * t));
}

static double alpha_from_solid_angle(double omega) {
    if (omega <= 0.0) return 0.0;
    if (omega >= LIGHT_TWO_PI) return LIGHT_PI;
    double c = 1.0 - omega / LIGHT_TWO_PI;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    return acos(c);
}

/* F(s) = 1 / (sqrt(1+s) * (sqrt(1+s) + 1)), so that
   Omega_disk = 2*pi*t^2*F(t^2) with t = R/d.  F is decreasing. */
static double f_of_s(double s) {
    double g = sqrt(1.0 + s);
    return 1.0 / (g * (g + 1.0));
}

/* q - 1 for a cone fitted at t_ref and evaluated at t = R/d. */
static double cone_omega_ratio(double t_ref, double t) {
    return f_of_s(t_ref * t_ref) / f_of_s(t * t);
}

static LightStatus cone_light_init(ConeLight* L, Vec3 axis, double radius,
                                   double dist_ref, double alpha_ref) {
    if (!L) return LIGHT_ERR_NULL;
    memset(L, 0, sizeof(*L));
    double al = len3(axis);
    if (!(al > 1e-12)) { L->status = LIGHT_ERR_AXIS; return L->status; }
    L->axis = scale3(axis, 1.0 / al);
    if (radius < 0.0) { L->status = LIGHT_ERR_RADIUS; return L->status; }
    if (!(dist_ref > 0.0)) { L->status = LIGHT_ERR_DISTANCE; return L->status; }
    if (alpha_ref < 0.0 || alpha_ref >= LIGHT_ALPHA_HARD_MAX) {
        L->status = LIGHT_ERR_ANGLE_RANGE; return L->status;
    }
    /* The closed-form distance correction assumes the cone is fitted to the
       disk at the reference distance, i.e. alpha_ref = atan(R/dist_ref). */
    if (fabs(alpha_ref - atan(radius / dist_ref)) > 1e-9) {
        L->status = LIGHT_ERR_ANGLE_RANGE; return L->status;
    }
    L->radius = radius;
    L->dist_ref = dist_ref;
    L->alpha_ref = alpha_ref;
    L->point = (alpha_ref == 0.0);
    L->status = LIGHT_OK;
    return L->status;
}

/* Effective half-angle at distance d after the 1/d^2 solid-angle
   correction.  Reports ERR_SOLID_ANGLE when the corrected solid angle
   leaves (0, 2*pi) (the approximation has broken down). */
static LightStatus cone_light_alpha_eff(const ConeLight* L, double d, double* alpha_eff) {
    if (!L || !alpha_eff) return LIGHT_ERR_NULL;
    if (L->point) { *alpha_eff = 0.0; return LIGHT_OK; }
    if (!(d > 0.0)) return LIGHT_ERR_DISTANCE;
    double omega = cone_solid_angle(L->alpha_ref) *
                   (L->dist_ref * L->dist_ref) / (d * d);
    if (!(omega > 0.0) || omega >= LIGHT_TWO_PI * (1.0 - 1e-15))
        return LIGHT_ERR_SOLID_ANGLE;
    *alpha_eff = alpha_from_solid_angle(omega);
    return LIGHT_OK;
}

/* Orthonormal tangent basis of n. */
static void basis3(Vec3 n, Vec3* t1, Vec3* t2) {
    Vec3 helper = (fabs(n.x) < 0.9) ? (Vec3){ 1, 0, 0 } : (Vec3){ 0, 1, 0 };
    Vec3 a = { n.y * helper.z - n.z * helper.y,
               n.z * helper.x - n.x * helper.z,
               n.x * helper.y - n.y * helper.x };
    double l = len3(a);
    if (l < 1e-12) { a = (Vec3){ 0, 1, 0 }; l = 1.0; }
    *t1 = scale3(a, 1.0 / l);
    Vec3 b = { n.y * t1->z - n.z * t1->y,
               n.z * t1->x - n.x * t1->z,
               n.x * t1->y - n.y * t1->x };
    *t2 = scale3(b, 1.0 / len3(b));
}

/* Uniform sampling inside the effective cone at distance d.
   pdf is per solid angle: 1 / Omega_cone(alpha_eff).  For a point
   light the direction is the axis and pdf is +inf (delta). */
static LightStatus cone_sample(const ConeLight* L, double d, double u1, double u2,
                               Vec3* dir, double* pdf) {
    if (!L || !dir || !pdf) return LIGHT_ERR_NULL;
    if (!(u1 >= 0.0 && u1 <= 1.0 && u2 >= 0.0 && u2 <= 1.0)) return LIGHT_ERR_UV;
    double a;
    LightStatus st = cone_light_alpha_eff(L, d, &a);
    if (st != LIGHT_OK) return st;
    if (L->point) { *dir = L->axis; *pdf = INFINITY; return LIGHT_OK; }

    double cos_a = cos(a);
    double cz = 1.0 - u1 * (1.0 - cos_a);          /* cos(theta) uniform */
    double sz = sqrt(fmax(0.0, 1.0 - cz * cz));
    double phi = LIGHT_TWO_PI * u2;
    Vec3 t1, t2;
    basis3(L->axis, &t1, &t2);
    Vec3 v = add3(scale3(L->axis, cz),
                  add3(scale3(t1, sz * cos(phi)), scale3(t2, sz * sin(phi))));
    *dir = v;
    *pdf = 1.0 / cone_solid_angle(a);
    return LIGHT_OK;
}

/* Compare the distance-corrected cone against the exact disk.
   alpha_eff is the half-angle of the cone whose solid angle equals the
   corrected solid angle, i.e. cos(alpha_eff) = 1 - q*(1 - cos(alpha_exact)).
   The soft-shadow width is the tangent of the half-angle, so the width
   ratio is tan(alpha_eff)/tan(alpha_exact); rel_err is |ratio - 1|. */
static LightStatus cone_width_compare(const ConeLight* L, double d,
                                      double* alpha_exact, double* alpha_eff,
                                      double* ratio, double* rel_err) {
    if (!L) return LIGHT_ERR_NULL;
    if (!(d > 0.0)) return LIGHT_ERR_DISTANCE;
    double a_ex = atan(L->radius / d);
    double t = L->radius / d;
    double u = 1.0 - 1.0 / sqrt(1.0 + t * t);   /* 1 - cos(alpha_exact) */
    double q = (L->point) ? 0.0 : cone_omega_ratio(L->radius / L->dist_ref, t);
    double cos_eff = 1.0 - q * u;
    if (!(cos_eff > -1.0)) return LIGHT_ERR_SOLID_ANGLE;
    double a_ef = (cos_eff >= 1.0) ? 0.0 : acos(cos_eff);
    double r = (t > 0.0) ? tan(a_ef) / t : 1.0;
    if (alpha_exact) *alpha_exact = a_ex;
    if (alpha_eff) *alpha_eff = a_ef;
    if (ratio) *ratio = r;
    if (rel_err) *rel_err = fabs(r - 1.0);
    return LIGHT_OK;
}

/* A-priori bound on |tan(alpha_eff)/tan(alpha_exact) - 1| for every
   t_ref, t in [0, tan(alpha_cap)] with alpha_cap within the supported
   range.  Derived from
       tan^2 ratio - 1 = (q - 1) * A(q,u),
       A(q,u) = (2 - (q+1)u) / ((1 - q u)^2 (2 - u)),
   with |q-1| <= Bq and |A| <= Amax over the domain, giving
       eW <= 1 - sqrt(1 - Bq*Amax).
   Returns NAN and sets *st when the bound is undefined. */
static double cone_width_error_bound(double alpha_cap, LightStatus* st) {
    if (st) *st = LIGHT_OK;
    if (!(alpha_cap >= 0.0) || alpha_cap >= LIGHT_ALPHA_HARD_MAX) {
        if (st) *st = LIGHT_ERR_ANGLE_RANGE;
        return NAN;
    }
    if (alpha_cap > LIGHT_ALPHA_BOUND_MAX) {
        if (st) *st = LIGHT_ERR_ANGLE_BOUND;
        return NAN;
    }
    double tc = tan(alpha_cap);
    double uc = 1.0 - cos(alpha_cap);
    double Ft = f_of_s(tc * tc);
    double qmin = 2.0 * Ft;      /* t_ref = 0,    t = tan(cap) */
    double qmax = 0.5 / Ft;      /* t_ref = cap,  t -> 0      */
    double Bq = fmax(fabs(qmin - 1.0), fabs(qmax - 1.0));
    if (qmax * uc >= 1.0) {
        if (st) *st = LIGHT_ERR_ANGLE_BOUND;
        return NAN;
    }
    /* |A(q,u)| <= (2 + (qmax+1)u) / ((1 - qmax*u)^2 (2 - u)) on the box
       [qmin,qmax] x [0,ucap]; the 2 in the numerator dominates. */
    double num = fmax(2.0, (qmax + 1.0) * uc);
    double Amax = num / ((1.0 - qmax * uc) * (1.0 - qmax * uc) * (2.0 - uc));
    double M = Bq * Amax;
    if (!(M < 1.0)) {
        if (st) *st = LIGHT_ERR_ANGLE_BOUND;
        return NAN;
    }
    return 1.0 - sqrt(1.0 - M);
}

/* ---- deterministic low-discrepancy helpers for the checks ---- */
static double radical_inverse(unsigned i, unsigned base) {
    double f = 1.0, r = 0.0;
    while (i > 0) { f /= (double)base; r += f * (double)(i % base); i /= base; }
    return r;
}

/* Monte Carlo solid angle of the disk by uniform area sampling.
   Shading point at the origin, disk center at axis*d, facing the origin. */
static double disk_solid_angle_area_mc(Vec3 axis, double R, double d, unsigned n) {
    Vec3 t1, t2;
    basis3(axis, &t1, &t2);
    Vec3 c = scale3(axis, d);
    double area = LIGHT_PI * R * R;
    double sum = 0.0;
    for (unsigned i = 0; i < n; i++) {
        double u1 = radical_inverse(i + 1u, 2u);
        double u2 = radical_inverse(i + 1u, 3u);
        double r = R * sqrt(u1);
        double phi = LIGHT_TWO_PI * u2;
        Vec3 p = add3(c, add3(scale3(t1, r * cos(phi)), scale3(t2, r * sin(phi))));
        double dist = len3(p);
        /* dOmega = (n_hat . (o - x)) / |x|^3 dA, o = 0, n_hat = -axis */
        double cos_term = dot3(axis, p);          /* = d > 0 */
        sum += cos_term / (dist * dist * dist);
    }
    return area * sum / (double)n;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    Vec3 z = { 0, 0, 1 };

    /* exact identity: disk solid angle == cone solid angle at atan(R/d) */
    {
        int ok = 1;
        double ts[] = { 0.01, 0.1, 0.5, 1.0, 2.0, 5.0 };
        for (int i = 0; i < 6; i++) {
            double t = ts[i];
            ok &= fabs(disk_solid_angle(t, 1.0) - cone_solid_angle(atan(t))) < 1e-12;
        }
        check("identity_disk_cone", ok);
    }
    check("cone_zero_solid_angle", cone_solid_angle(0.0) == 0.0);
    check("cone_halfpi_solid_angle",
          fabs(cone_solid_angle(LIGHT_PI * 0.5) - LIGHT_TWO_PI) < 1e-12);
    {
        double a = 0.4;
        check("alpha_roundtrip",
              fabs(alpha_from_solid_angle(cone_solid_angle(a)) - a) < 1e-12);
    }

    /* cone sampling: unit directions, inside the cone, correct pdf */
    {
        ConeLight L;
        double a = 25.0 * LIGHT_PI / 180.0;
        double dref = 0.5;
        LightStatus st = cone_light_init(&L, z, dref * tan(a), dref, a);
        check("init_ok", st == LIGHT_OK && L.status == LIGHT_OK && !L.point);
        double cos_a = cos(a);
        int inside = 1, unit = 1;
        double mean_cos = 0.0, pdf_ok = 1.0;
        unsigned n = 4096;
        for (unsigned i = 0; i < n; i++) {
            double u1 = radical_inverse(i + 1u, 2u);
            double u2 = radical_inverse(i + 1u, 3u);
            Vec3 dir; double pdf;
            st = cone_sample(&L, dref, u1, u2, &dir, &pdf);
            if (st != LIGHT_OK) { inside = 0; break; }
            double l = len3(dir);
            if (fabs(l - 1.0) > 1e-12) unit = 0;
            if (dot3(dir, L.axis) < cos_a - 1e-12) inside = 0;
            mean_cos += dot3(dir, L.axis);
            if (fabs(pdf * cone_solid_angle(a) - 1.0) > 1e-12) pdf_ok = 0.0;
        }
        check("sample_unit", unit);
        check("sample_inside_cone", inside);
        check("sample_pdf_normalized", pdf_ok > 0.5);
        double expect = 0.5 * (1.0 + cos_a);
        check("sample_mean_cos", fabs(mean_cos / (double)n - expect) < 2e-3);
    }

    /* distance correction: exact at d_ref, matches F-ratio elsewhere */
    {
        ConeLight L;
        double a = 30.0 * LIGHT_PI / 180.0;
        double dref = 0.5, R = dref * tan(a);
        cone_light_init(&L, z, R, dref, a);
        double ae0;
        LightStatus st = cone_light_alpha_eff(&L, dref, &ae0);
        check("corr_exact_at_ref", st == LIGHT_OK && fabs(ae0 - a) < 1e-12);

        int ok = 1;
        double ds[] = { 0.45, 0.5, 0.8, 1.5, 4.0 };
        for (int i = 0; i < 5; i++) {
            double d = ds[i];
            double ae;
            if (cone_light_alpha_eff(&L, d, &ae) != LIGHT_OK) { ok = 0; continue; }
            double oc = cone_solid_angle(a) * (dref * dref) / (d * d);
            ok &= fabs(cone_solid_angle(ae) - oc) < 1e-12;
            double q = oc / disk_solid_angle(R, d);
            double qf = cone_omega_ratio(R / dref, R / d);
            ok &= fabs(q - qf) < 1e-12;
        }
        check("corr_matches_formula", ok);
    }

    /* exact area sampling agrees with the closed-form disk solid angle */
    {
        int ok = 1;
        double ds[] = { 1.0, 2.0, 5.0 };
        for (int i = 0; i < 3; i++) {
            double d = ds[i], R = 0.5;
            double exact = disk_solid_angle(R, d);
            double mc = disk_solid_angle_area_mc(z, R, d, 32768);
            ok &= fabs(mc - exact) / exact < 0.01;
        }
        check("area_sampling_matches", ok);
    }

    /* width comparison and the a-priori error bound:
       the bound holds for every state with alpha_eff <= alpha_cap. */
    {
        ConeLight L;
        double a = 30.0 * LIGHT_PI / 180.0;
        double dref = 0.5, R = dref * tan(a);
        cone_light_init(&L, z, R, dref, a);
        LightStatus bst;
        double bound = cone_width_error_bound(a, &bst);
        check("bound_defined", bst == LIGHT_OK && bound > 0.0 && bound < 1.0);

        double max_err = 0.0;
        int within = 1, states = 0;
        /* sweep distances so t = R/d covers (0, tan(a)]; t_ref = tan(a). */
        for (int i = 1; i <= 2000; i++) {
            double t = tan(a) * (double)i / 2000.0;
            double d = R / t;
            double a_ex, a_ef, r, err;
            if (cone_width_compare(&L, d, &a_ex, &a_ef, &r, &err) != LIGHT_OK) continue;
            if (a_ef > a + 1e-12) continue;      /* outside the bound's domain */
            states++;
            if (err > max_err) max_err = err;
            if (!(err <= bound + 1e-12)) within = 0;
        }
        check("width_states_exercised", states > 1000);
        check("width_err_within_bound", within);
        printf("width_bound=%.6f\n", bound);
        printf("width_max_err=%.6f\n", max_err);

        /* zero distance correction error at the reference point */
        double a_ex, a_ef, r, err;
        cone_width_compare(&L, L.dist_ref, &a_ex, &a_ef, &r, &err);
        check("width_exact_at_ref", fabs(err) < 1e-12 && fabs(r - 1.0) < 1e-12);

        /* monotonicity: the bound grows with the cap angle */
        double b25 = cone_width_error_bound(25.0 * LIGHT_PI / 180.0, NULL);
        double b35 = cone_width_error_bound(35.0 * LIGHT_PI / 180.0, NULL);
        check("bound_monotone", b25 < bound && b35 > bound);

        /* the bound is global: it holds for every reference half-angle and
           every distance, over all states with alpha_eff <= cap */
        int global_ok = 1, nstates = 0;
        for (int ri = 1; ri <= 6; ri++) {
            double aref = (double)ri * 5.0 * LIGHT_PI / 180.0;
            double dref2 = 1.0, R2 = dref2 * tan(aref);
            ConeLight M;
            if (cone_light_init(&M, z, R2, dref2, aref) != LIGHT_OK) { global_ok = 0; continue; }
            for (int k = 0; k <= 2000; k++) {
                double d = exp(log(0.05) + (log(60.0) - log(0.05)) * (double)k / 2000.0);
                double a_ex, a_ef, r, err;
                if (cone_width_compare(&M, d, &a_ex, &a_ef, &r, &err) != LIGHT_OK) continue;
                if (a_ex > a + 1e-12 || a_ef > a + 1e-12) continue;
                nstates++;
                if (!(err <= bound + 1e-12)) global_ok = 0;
            }
        }
        check("bound_global_across_refs", global_ok && nstates > 3000);
    }

    /* zero half-angle: point (delta) light degenerate */
    {
        ConeLight P;
        LightStatus st = cone_light_init(&P, z, 0.0, 0.5, 0.0);
        check("point_init", st == LIGHT_OK && P.point == 1);
        double ae = -1.0;
        st = cone_light_alpha_eff(&P, 0.7, &ae);
        check("point_alpha_zero", st == LIGHT_OK && ae == 0.0);
        Vec3 dir; double pdf;
        st = cone_sample(&P, 0.7, 0.3, 0.6, &dir, &pdf);
        check("point_sample_axis",
              st == LIGHT_OK && fabs(dir.x) < 1e-15 && fabs(dir.y) < 1e-15 &&
              fabs(dir.z - 1.0) < 1e-15 && isinf(pdf));
        check("point_solid_angle", cone_solid_angle(0.0) == 0.0);
    }

    /* out-of-range angles are reported with distinct statuses */
    {
        ConeLight L;
        LightStatus st;
        double R30 = 0.5 * tan(30.0 * LIGHT_PI / 180.0);

        st = cone_light_init(&L, z, 0.2, 0.5, LIGHT_PI * 0.5);
        check("angle_range_reported", st == LIGHT_ERR_ANGLE_RANGE);

        st = cone_light_init(&L, z, 0.2, 0.5, -0.1);
        check("angle_negative_reported", st == LIGHT_ERR_ANGLE_RANGE);

        /* cone half-angle inconsistent with atan(R/d_ref) is reported */
        st = cone_light_init(&L, z, 0.2, 0.5, 30.0 * LIGHT_PI / 180.0);
        check("angle_inconsistent_reported", st == LIGHT_ERR_ANGLE_RANGE);

        /* physically valid but beyond the finite-bound range */
        double wide = 40.0 * LIGHT_PI / 180.0;
        st = cone_light_init(&L, z, 0.5 * tan(wide), 0.5, wide);
        check("angle_init_ok_wide", st == LIGHT_OK);
        LightStatus bst;
        double b = cone_width_error_bound(wide, &bst);
        check("angle_bound_reported",
              bst == LIGHT_ERR_ANGLE_BOUND && isnan(b));

        /* distance correction overshoot: corrected solid angle >= 2*pi */
        cone_light_init(&L, z, R30, 0.5, 30.0 * LIGHT_PI / 180.0);
        double ae;
        st = cone_light_alpha_eff(&L, 0.05, &ae);
        check("solid_angle_overshoot_reported", st == LIGHT_ERR_SOLID_ANGLE);

        st = cone_light_init(&L, z, -1.0, 0.5, 0.2);
        check("radius_reported", st == LIGHT_ERR_RADIUS);
        st = cone_light_init(&L, z, 0.0, 0.0, 0.0);
        check("distance_reported", st == LIGHT_ERR_DISTANCE);
        Vec3 zero = { 0, 0, 0 };
        st = cone_light_init(&L, zero, 0.2, 0.5, atan(0.2 / 0.5));
        check("axis_reported", st == LIGHT_ERR_AXIS);
        check("status_strings",
              strcmp(light_status_str(LIGHT_ERR_ANGLE_BOUND), "ERR_ANGLE_BOUND") == 0 &&
              strcmp(light_status_str(LIGHT_ERR_SOLID_ANGLE), "ERR_SOLID_ANGLE") == 0 &&
              strcmp(light_status_str(LIGHT_OK), "OK") == 0);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
