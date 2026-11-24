
/* spot: penumbra radius of a spot cone from solid angle x incidence, with
   edge-profile measurement and pixel/aliasing accounting. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

static double clamp01(double x) { return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x); }
static double smoothstep(double x) { x = clamp01(x); return x * x * (3.0 - 2.0 * x); }

typedef struct { double tin, tout, te, band; } Cone;

static Cone cone_make(double tin_deg, double tout_deg) {
    Cone c;
    c.tin = tin_deg * PI / 180.0;
    c.tout = tout_deg * PI / 180.0;
    c.band = c.tout - c.tin;
    c.te = 0.5 * (c.tin + c.tout);
    return c;
}

/* 1 inside the inner cone, 0 outside the outer, smoothstep in between */
static double cone_I(const Cone* c, double th) {
    return 1.0 - smoothstep((th - c->tin) / c->band);
}

static double cone_iso(const Cone* c, double level) {
    double lo = c->tin, hi = c->tout;
    for (int i = 0; i < 200; i++) {
        double m = 0.5 * (lo + hi);
        if (cone_I(c, m) > level) lo = m; else hi = m;
    }
    return 0.5 * (lo + hi);
}

/* solid angle of the 10-90 annulus; exact inverse for the symmetric band */
static double band_domega(const Cone* c, double* dtheta) {
    double t90 = cone_iso(c, 0.9), t10 = cone_iso(c, 0.1);
    double dOm = 2.0 * PI * (cos(t90) - cos(t10));
    *dtheta = 2.0 * asin(dOm / (4.0 * PI * sin(c->te)));
    return dOm;
}

/* receiver plane at axial distance d, its normal at incidence alpha to the
   edge ray. gamma is the normal's tilt from the axis. */
typedef struct { double d, alpha; } Recv;

static double recv_gamma(const Recv* r, const Cone* c) { return c->te - r->alpha; }
static double recv_rho(const Recv* r, const Cone* c) {
    return r->d * cos(recv_gamma(r, c)) / cos(r->alpha);
}
static double recv_cos_alpha(const Recv* r, const Cone* c) { (void)c; return cos(r->alpha); }

/* radius from the two terms: distance rho and incidence 1/cos(alpha) */
static double radius_linear(double rho, double cos_a, double dtheta) {
    return rho * dtheta / cos_a;
}
/* same band mapped exactly onto the receiver plane (tan instead of arc) */
static double radius_exact(double rho, double cos_a, double dtheta) {
    double a = acos(cos_a), h = 0.5 * dtheta;
    return rho * (sin(h) / cos(a + h) + sin(h) / cos(a - h));
}

static void recv_edge_point(const Recv* r, const Cone* c, double* px, double* pz) {
    double rho = recv_rho(r, c);
    *px = rho * sin(c->te);
    *pz = rho * cos(c->te);
}

/* intensity seen by a pixel at surface offset s from the geometric edge */
static double recv_V(const Recv* r, const Cone* c, double s) {
    double px, pz;
    recv_edge_point(r, c, &px, &pz);
    double qx = px + s * cos(recv_gamma(r, c)), qz = pz - s * sin(recv_gamma(r, c));
    return cone_I(c, atan2(qx, qz));
}

/* surface offset at which the edge ray of polar angle th lands */
static double recv_s_of_theta(const Recv* r, const Cone* c, double th) {
    double px, pz;
    recv_edge_point(r, c, &px, &pz);
    double rad = r->d * cos(recv_gamma(r, c)) / cos(th - recv_gamma(r, c));
    double x = rad * sin(th), z = rad * cos(th);
    return cos(recv_gamma(r, c)) * (x - px) - sin(recv_gamma(r, c)) * (z - pz);
}

static double recv_level_s(const Recv* r, const Cone* c, double level) {
    double a = recv_s_of_theta(r, c, c->tin), b = recv_s_of_theta(r, c, c->tout);
    if (recv_V(r, c, a) < recv_V(r, c, b)) { double t = a; a = b; b = t; }
    for (int i = 0; i < 200; i++) {
        double m = 0.5 * (a + b);
        if (recv_V(r, c, m) > level) a = m; else b = m;
    }
    return 0.5 * (a + b);
}

typedef struct { int n_band, n_full; double step_max; } EdgePix;

/* sample the receiver on a grid of pitch h, one pixel boundary on the edge */
static EdgePix recv_pixels(const Recv* r, const Cone* c, double h) {
    EdgePix e; memset(&e, 0, sizeof(e));
    double s_in = recv_s_of_theta(r, c, c->tin), s_out = recv_s_of_theta(r, c, c->tout);
    double s90 = recv_level_s(r, c, 0.9), s10 = recv_level_s(r, c, 0.1);
    double lo = fmin(s_in, s_out), hi = fmax(s_in, s_out);
    int k0 = (int)floor(lo / h) - 2, k1 = (int)ceil(hi / h) + 2;
    double prev = -1.0;
    for (int k = k0; k <= k1; k++) {
        double s = (k + 0.5) * h;
        double v = recv_V(r, c, s);
        if (v > 0.0 && v < 1.0) e.n_full++;
        if ((s - s90) * (s - s10) <= 0.0) e.n_band++;
        if (prev >= 0.0) { double d = fabs(v - prev); if (d > e.step_max) e.step_max = d; }
        prev = v;
    }
    return e;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    const double H = 0.002;   /* pixel pitch, 2 mm */
    Cone c = cone_make(23.5, 26.5);
    double dtheta, dOm = band_domega(&c, &dtheta);
    double t90 = cone_iso(&c, 0.9), t10 = cone_iso(&c, 0.1);

    check("center_half", fabs(cone_I(&c, c.te) - 0.5) < 1e-12);
    check("iso_ordered", t90 < c.te && c.te < t10);
    check("solid_angle_positive", dOm > 0.0 && dOm < 2.0 * PI * (1.0 - cos(c.tout)));
    check("band_inverse_exact", fabs(dtheta - (t10 - t90)) < 1e-12);
    check("edge_offset_zero", fabs(recv_s_of_theta(&(Recv){1.0, 0.3}, &c, c.te)) < 1e-12);

    /* radius must grow with distance at fixed incidence, and with incidence
       at fixed distance */
    {
        Recv near_ = { 2.0, 20.0 * PI / 180.0 }, far_ = { 12.0, 20.0 * PI / 180.0 };
        Recv flat = { 6.0, 20.0 * PI / 180.0 }, graze = { 6.0, 60.0 * PI / 180.0 };
        double cn = recv_cos_alpha(&near_, &c), cf = recv_cos_alpha(&far_, &c);
        double cflat = recv_cos_alpha(&flat, &c), cgraze = recv_cos_alpha(&graze, &c);
        check("radius_distance", radius_linear(recv_rho(&far_, &c), cf, dtheta) >
                                 5.0 * radius_linear(recv_rho(&near_, &c), cn, dtheta));
        check("radius_incidence", radius_linear(recv_rho(&graze, &c), cgraze, dtheta) >
                                  2.0 * radius_linear(recv_rho(&flat, &c), cflat, dtheta));
    }

    printf("dtheta90=%.12f domega90=%.12f band=%.6f deg\n",
           dtheta, dOm, c.band * 180.0 / PI);
    printf("%-6s %6s %7s %8s %10s %10s %10s %7s %8s %8s\n",
           "recv", "d", "alpha", "cos_a", "rho", "R_pred", "W_meas", "N_90", "N_px", "A_max");
    printf("%-6s %6s %7s %8s %10s %10s %10s %7s %8s %8s\n",
           "", "m", "deg", "", "mm", "mm", "mm", "band", "0-1", "per px");

    const char* names[3] = { "A_near", "B_mid", "C_far" };
    /* axial distance d, incidence alpha; near/flat, mid, far/grazing */
    Recv rs[3] = { { 2.0, 5.0 * PI / 180.0 },
                   { 6.0, 40.0 * PI / 180.0 },
                   { 12.0, 65.0 * PI / 180.0 } };
    double Wm[3], Rp[3], Rw[3], cosA[3], rho[3], A[3], Nf[3], Wex[3];
    int Nband[3];
    EdgePix ep[3];

    for (int i = 0; i < 3; i++) {
        double ca = recv_cos_alpha(&rs[i], &c);
        rho[i] = recv_rho(&rs[i], &c);
        cosA[i] = ca;
        Rp[i] = radius_linear(rho[i], ca, dtheta);
        Wex[i] = radius_exact(rho[i], ca, dtheta);
        Wm[i] = fabs(recv_level_s(&rs[i], &c, 0.1) - recv_level_s(&rs[i], &c, 0.9));
        ep[i] = recv_pixels(&rs[i], &c, H);
        Nband[i] = ep[i].n_band;
        A[i] = ep[i].step_max;
        Nf[i] = (double)ep[i].n_full;
        Rw[i] = rho[i] * dtheta * ca;
        printf("%-6s %6.1f %7.2f %8.5f %10.2f %10.2f %10.2f %7d %8.0f %8.5f\n",
               names[i], rs[i].d, rs[i].alpha * 180.0 / PI, ca, rho[i] * 1000.0,
               Rp[i] * 1000.0, Wm[i] * 1000.0, Nband[i], Nf[i], A[i]);
    }

    /* the two-term radius must reproduce the measured 10-90 width */
    double worst = 0.0, second = 0.0;
    for (int i = 0; i < 3; i++) {
        double e = fabs(Wm[i] / Wex[i] - 1.0);
        if (e > worst) worst = e;
        e = Wm[i] / Rp[i] - 1.0;
        if (e > second) second = e;
        char nm[32];
        snprintf(nm, sizeof(nm), "width_%s", names[i]);
        check(nm, fabs(Wm[i] / Wex[i] - 1.0) < 1e-6);
    }
    printf("width_rel_err_exact=%.3e  exceed_linear=%.6f\n", worst, second);
    check("no_exceedance", second < 0.01);

    /* sign audit: if the incidence term were 1/cos instead of cos, the width
       would be rho*dtheta*cos(a); measured/wrong should equal sec^2(a) */
    printf("%-6s %8s %11s %11s %10s %10s %10s\n",
           "recv", "cos_a", "rho*dtheta", "1/cos_a", "R_pred", "W_meas", "W/R_wrong");
    printf("%-6s %8s %11s %11s %10s %10s %10s\n",
           "", "", "mm", "", "mm", "mm", "");
    for (int i = 0; i < 3; i++) {
        printf("%-6s %8.5f %11.2f %11.5f %10.2f %10.2f %10.4f\n",
               names[i], cosA[i], rho[i] * dtheta * 1000.0, 1.0 / cosA[i],
               Rp[i] * 1000.0, Wm[i] * 1000.0, Wm[i] / Rw[i]);
    }
    {
        double sec2 = 1.0 / (cosA[2] * cosA[2]);
        printf("C_far: 1/cos_a=%.4f  1/cos_a^2=%.4f  measured/R_wrong=%.4f\n",
               1.0 / cosA[2], sec2, Wm[2] / Rw[2]);
    }
    check("wrong_sign_rejected",
          fabs(Wm[2] / Rw[2] / (1.0 / (cosA[2] * cosA[2])) - 1.0) < 0.02 &&
          fabs(Wm[2] / Rw[2] / (1.0 / cosA[2]) - 1.0) > 1.0);

    /* radius -> pixel count: more penumbra pixels, smaller per-pixel step */
    check("pixels_track_radius",
          Nband[0] < Nband[1] && Nband[1] < Nband[2] &&
          fabs(Nband[2] / (Rp[2] / H) - 1.0) < 0.05);
    check("step_falls_with_radius", A[0] > A[1] && A[1] > A[2]);
    check("step_times_pixels_flat",
          fabs(A[0] * Nf[0] - 1.5) < 0.1 &&
          fabs(A[1] * Nf[1] - 1.5) < 0.1 &&
          fabs(A[2] * Nf[2] - 1.5) < 0.1);
    printf("A_max*N_full: %.3f %.3f %.3f\n",
           A[0] * Nf[0], A[1] * Nf[1], A[2] * Nf[2]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
