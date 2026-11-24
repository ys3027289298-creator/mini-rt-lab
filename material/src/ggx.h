/* ggx.h: GGX microfacet BRDF, visible-normal (VNDF) sampling, and a
   second-order Taylor model of q(t) = pdf * brdf near grazing incidence
   (t = cos(theta_i) -> 0). Double precision throughout: the Taylor
   coefficients are taken by finite differences and the control-variate
   quadrature integrates products of very large pdf/brdf values. */
#ifndef GGX_H
#define GGX_H

#include <math.h>

#define GGX_PI 3.14159265358979323846
#define GGX_ALPHA_MIN 1e-6 /* below this the material is treated as a mirror */

typedef struct { double x, y, z; } V3d;

static double ggx_dot(V3d a, V3d b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3d ggx_add(V3d a, V3d b) { V3d r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static V3d ggx_norm(V3d a) {
    double l = sqrt(ggx_dot(a, a));
    V3d r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}

/* GGX normal distribution, alpha = roughness^2, macro normal = +z */
static double ggx_D(double hz, double alpha) {
    double a2 = alpha * alpha;
    double d = hz * hz * (a2 - 1.0) + 1.0;
    return a2 / (GGX_PI * d * d);
}

/* Smith Lambda for GGX; even in mu (used as a smooth extension for mu<0
   inside the Taylor finite differences, never as a physical masking term). */
static double ggx_lambda(double mu, double alpha) {
    double a2 = alpha * alpha;
    double m2 = mu * mu;
    if (m2 < 1e-18) return 1e9; /* only hit at exactly mu=0, callers avoid */
    return (-1.0 + sqrt(1.0 + a2 * (1.0 - m2) / m2)) * 0.5;
}

static double ggx_G1(double mu, double alpha) {
    return 1.0 / (1.0 + ggx_lambda(mu, alpha));
}

/* height-correlated Smith G2 */
static double ggx_G2(double mui, double muo, double alpha) {
    return 1.0 / (1.0 + ggx_lambda(muo, alpha) + ggx_lambda(mui, alpha));
}

/* G2/mui, stable at mui = 0: t*Lambda_i = (-t + sqrt(t^2 + a^2(1-t^2)))/2
   for t >= 0, so G2/t -> 1/(alpha/2) as t -> 0+. */
static double ggx_G2_over_t(double t, double muo, double alpha) {
    double a2 = alpha * alpha;
    double lam_o = ggx_lambda(muo, alpha);
    double t_lam_i;
    if (t >= 0.0)
        t_lam_i = (-t + sqrt(t * t + a2 * (1.0 - t * t))) * 0.5;
    else
        t_lam_i = (-t - sqrt(t * t + a2 * (1.0 - t * t))) * 0.5;
    return 1.0 / (t * (1.0 + lam_o) + t_lam_i);
}

static double ggx_schlick(double cos_theta, double F0) {
    double m = 1.0 - cos_theta;
    double m2 = m * m;
    return F0 + (1.0 - F0) * m2 * m2 * m;
}

/* full GGX BRDF f_r(wo, wi), macro normal +z, both dirs above horizon */
static double brdf_ggx(V3d wo, V3d wi, double alpha, double F0) {
    V3d h = ggx_norm(ggx_add(wo, wi));
    double D = ggx_D(h.z, alpha);
    double G2 = ggx_G2(wi.z, wo.z, alpha);
    double F = ggx_schlick(ggx_dot(wi, h), F0);
    return F * D * G2 / (4.0 * wi.z * wo.z);
}

/* VNDF sampling (Heitz 2018, "Sampling the GGX Distribution of Visible
   Normals"). Returns the reflected incident direction wi. alpha is clamped
   to GGX_ALPHA_MIN: at alpha = 0 the visible-normal distribution collapses
   to a Dirac delta on the macro normal and sampling degenerates to the
   deterministic mirror reflection wi = 2(wo.z)z - wo. */
static V3d vndf_sample(V3d wo, double alpha, double u1, double u2) {
    if (alpha < GGX_ALPHA_MIN) alpha = GGX_ALPHA_MIN;
    V3d Vh = ggx_norm((V3d){ alpha * wo.x, alpha * wo.y, wo.z });
    double lensq = Vh.x * Vh.x + Vh.y * Vh.y;
    V3d T1 = lensq > 0.0
        ? (V3d){ -Vh.y / sqrt(lensq), Vh.x / sqrt(lensq), 0.0 }
        : (V3d){ 1.0, 0.0, 0.0 };
    V3d T2 = { Vh.y * T1.z - Vh.z * T1.y,
               Vh.z * T1.x - Vh.x * T1.z,
               Vh.x * T1.y - Vh.y * T1.x };
    double r = sqrt(u1);
    double phi = 2.0 * GGX_PI * u2;
    double t1 = r * cos(phi), t2 = r * sin(phi);
    double s = 0.5 * (1.0 + Vh.z);
    t2 = (1.0 - s) * sqrt(1.0 - t1 * t1) + s * t2;
    double z = sqrt(1.0 - t1 * t1 - t2 * t2 > 0.0 ? 1.0 - t1 * t1 - t2 * t2 : 0.0);
    V3d Nh = { t1 * T1.x + t2 * T2.x + z * Vh.x,
               t1 * T1.y + t2 * T2.y + z * Vh.y,
               t1 * T1.z + t2 * T2.z + z * Vh.z };
    V3d Ne = ggx_norm((V3d){ alpha * Nh.x, alpha * Nh.y,
                             Nh.z > 0.0 ? Nh.z : 0.0 });
    double d = ggx_dot(wo, Ne);
    V3d wi = { 2.0 * d * Ne.x - wo.x, 2.0 * d * Ne.y - wo.y, 2.0 * d * Ne.z - wo.z };
    return wi;
}

/* solid-angle pdf of vndf_sample: p(wi) = D(h) G1(wo) / (4 wo.z) */
static double vndf_pdf(V3d wo, V3d wi, double alpha) {
    if (alpha < GGX_ALPHA_MIN) alpha = GGX_ALPHA_MIN;
    V3d h = ggx_norm(ggx_add(wo, wi));
    return ggx_D(h.z, alpha) * ggx_G1(wo.z, alpha) / (4.0 * wo.z);
}

/* q(t, phi) = pdf * brdf along the meridian at azimuth phi, with
   wi = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)), t = cos(theta).
   Smooth in t at t = 0 for alpha > 0 (G2/t has a finite limit), which is
   what makes a low-order grazing expansion possible. */
static double q_product(V3d wo, double t, double phi, double alpha, double F0) {
    if (alpha < GGX_ALPHA_MIN) alpha = GGX_ALPHA_MIN;
    double st = sqrt(1.0 - t * t > 0.0 ? 1.0 - t * t : 0.0);
    V3d wi = { st * cos(phi), st * sin(phi), t };
    V3d h = ggx_norm(ggx_add(wo, wi));
    double D = ggx_D(h.z, alpha);
    double G1o = ggx_G1(wo.z, alpha);
    double p = D * G1o / (4.0 * wo.z);
    double F = ggx_schlick(ggx_dot(wi, h), F0);
    double fr = F * D * ggx_G2_over_t(t, wo.z, alpha) / (4.0 * wo.z);
    return p * fr;
}

/* Second-order Taylor expansion of q(t) at t = 0 (grazing), fixed phi:
   q(t) ~= c0 + c1 t + c2 t^2 / 2.
   c0 is the exact grazing limit; c1, c2 use one-sided second-order finite
   differences (the smooth t>=0 branch only — the even-Lambda extension for
   t<0 is a different branch and must not be mixed into the stencil). */
static void q_taylor(V3d wo, double phi, double alpha, double F0,
                     double *c0, double *c1, double *c2) {
    /* stencil must resolve q's feature scale, which shrinks with alpha */
    double h = 1e-4;
    if (alpha < 0.01) h = 1e-4 * alpha / 0.01;
    if (h < 1e-8) h = 1e-8;
    double q0 = q_product(wo, 0.0, phi, alpha, F0);
    double q1 = q_product(wo, h, phi, alpha, F0);
    double q2 = q_product(wo, 2.0 * h, phi, alpha, F0);
    double q3 = q_product(wo, 3.0 * h, phi, alpha, F0);
    *c0 = q0;
    *c1 = (-3.0 * q0 + 4.0 * q1 - q2) / (2.0 * h);
    *c2 = (2.0 * q0 - 5.0 * q1 + 4.0 * q2 - q3) / (h * h);
}

static double q_taylor_eval(double c0, double c1, double c2, double t) {
    double q = c0 + c1 * t + 0.5 * c2 * t * t;
    return q > 0.0 ? q : 0.0; /* truncated model may dip below 0; clamp */
}

#endif
