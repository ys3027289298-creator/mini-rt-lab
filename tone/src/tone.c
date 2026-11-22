
/* tone: reversible tone curve.
 *
 * Model
 *   C1 piecewise-cubic Hermite over knots x[0..n-1] inside [0,1].
 *   g[k] is the local gain  d(output)/d(input)  at knot x[k]; g=1 is identity.
 *   Each interval is normalised so f(x[i+1]) - f(x[i]) = x[i+1] - x[i], and the
 *   curve is pinned to f(x[0]) = 0.  Outside [x[0], x[n-1]] it is extended
 *   linearly with the end slopes g[0] / g[n-1], so f is C1 everywhere.
 *
 * Monotonicity (exact, not sampled)
 *   On interval i, u = (t - x[i]) / h, the slope is the quadratic
 *     L(u) = g[i] + (6 - 4 g[i] - 2 g[i+1]) u + (3 g[i] + 3 g[i+1] - 6) u^2
 *   with L(0)=g[i], L(1)=g[i+1].  A quadratic's minimum on [0,1] is at an
 *   endpoint or at its (single) interior vertex, so L<0 can be located
 *   exactly; intervals where L<0 are reported as non-monotone.
 *
 * Inverse
 *   The extension regions invert in closed form; the interior uses bisection
 *   in x, stopped at xTol = 1e-7 * (1 + |x|).  Bisection assumes L >= 0.
 *
 * Round-trip bounds  (derived, see analysis_print)
 *   A: y = f(x); x' = f^-1(y);  |x' - x| <= xTol + yTol/Lmin
 *   B: x = f^-1(y); y' = f(x);  |y' - y| <= Lmax * xTol + yTol
 *   xTol is the bisection tolerance (2e-7 for x in [0,1]),
 *   yTol = 16 ulp (~1.9e-6 for y in [0,1]) covers float evaluation noise.
 *
 * Analytic comparison
 *   For the symmetric single-interval family  x = {0,1}, g[0] = g[1] = g
 *   the curve is the cubic A u^3 + B u^2 + C u with A=2g-2, B=3-3g, C=g,
 *   whose inverse is closed form (Cardano).  The numeric inverse is checked
 *   against it.
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define NMAX 16
#define NONMONO_MAX 8

#define TONE_XTOL_REL 1e-7f          /* bisection stop, relative to 1+|x| */
#define TONE_XTOL_ABS (TONE_XTOL_REL * 2.0)     /* worst case for x in [0,1] */
#define TONE_YTOL_ABS (16.0 * FLT_EPSILON)      /* eval noise for y in [0,1] */

typedef struct { float x[NMAX]; float g[NMAX]; int n; } Curve;

/* ------------------------------------------------------------------ */
/* parameter validation                                                */

typedef enum {
    TONE_OK = 0,
    TONE_E_N, TONE_E_X_ORDER, TONE_E_X_RANGE, TONE_E_X_NAN,
    TONE_E_G_RANGE, TONE_E_G_NAN
} ToneErr;

typedef struct {
    ToneErr err;
    int index;          /* knot index, -1 for n-level errors */
    float value;        /* the offending value */
    const char* field;
} ToneReport;

static const char* tone_err_name(ToneErr e) {
    switch (e) {
    case TONE_OK:        return "ok";
    case TONE_E_N:       return "n_out_of_range";
    case TONE_E_X_ORDER: return "x_not_increasing";
    case TONE_E_X_RANGE: return "x_out_of_unit_interval";
    case TONE_E_X_NAN:   return "x_not_finite";
    case TONE_E_G_RANGE: return "gain_not_positive";
    case TONE_E_G_NAN:   return "gain_not_finite";
    }
    return "?";
}

static ToneReport tone_validate(const Curve* c) {
    ToneReport r;
    memset(&r, 0, sizeof(r));
    r.err = TONE_OK;
    r.index = -1;
    r.field = "n";
    if (c->n < 2 || c->n > NMAX) {
        r.err = TONE_E_N;
        r.value = (float)c->n;
        return r;
    }
    for (int i = 0; i < c->n; i++) {
        if (!isfinite(c->x[i])) {
            r.err = TONE_E_X_NAN; r.index = i; r.value = c->x[i]; r.field = "x";
            return r;
        }
    }
    for (int i = 0; i < c->n; i++) {
        if (c->x[i] < 0.0f || c->x[i] > 1.0f) {
            r.err = TONE_E_X_RANGE; r.index = i; r.value = c->x[i]; r.field = "x";
            return r;
        }
    }
    for (int i = 1; i < c->n; i++) {
        if (!(c->x[i] > c->x[i-1])) {
            r.err = TONE_E_X_ORDER; r.index = i; r.value = c->x[i]; r.field = "x";
            return r;
        }
    }
    for (int i = 0; i < c->n; i++) {
        if (!isfinite(c->g[i])) {
            r.err = TONE_E_G_NAN; r.index = i; r.value = c->g[i]; r.field = "g";
            return r;
        }
    }
    for (int i = 0; i < c->n; i++) {
        if (c->g[i] <= 0.0f) {
            r.err = TONE_E_G_RANGE; r.index = i; r.value = c->g[i]; r.field = "g";
            return r;
        }
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* evaluation                                                          */

/* slope of interval i at local parameter u in [0,1] */
static double seg_slope(const Curve* c, int i, double u) {
    double gi = c->g[i], gj = c->g[i+1];
    return gi + (6.0 - 4.0 * gi - 2.0 * gj) * u
              + (3.0 * gi + 3.0 * gj - 6.0) * u * u;
}

static float tone_eval(const Curve* c, float t) {
    int n = c->n;
    if (t <= c->x[0]) return c->g[0] * (t - c->x[0]);
    if (t >= c->x[n-1])
        return (c->x[n-1] - c->x[0]) + c->g[n-1] * (t - c->x[n-1]);

    int i = 0;
    while (i < n - 2 && t > c->x[i+1]) i++;

    float h = c->x[i+1] - c->x[i];
    float u = (t - c->x[i]) / h;
    float u2 = u * u, u3 = u2 * u;
    float gi = c->g[i], gj = c->g[i+1];
    /* Hermite H(0)=0, H(1)=1, H'(0)=gi, H'(1)=gj */
    float H = gi * (u3 - 2.0f * u2 + u)
            + (-2.0f * u3 + 3.0f * u2)
            + gj * (u3 - u2);
    return (c->x[i] - c->x[0]) + h * H;
}

/* inverse; returns 0 only for a non-finite output value */
static int tone_inv(const Curve* c, float y, float* out) {
    int n = c->n;
    if (!isfinite(y)) return 0;

    float ytop = c->x[n-1] - c->x[0];
    if (y <= 0.0f) {                     /* linear left extension */
        *out = c->x[0] + y / c->g[0];
        return 1;
    }
    if (y >= ytop) {                     /* linear right extension */
        *out = c->x[n-1] + (y - ytop) / c->g[n-1];
        return 1;
    }

    int i = 0;                           /* y[i] == x[i] - x[0] */
    while (i < n - 2 && y > (c->x[i+1] - c->x[0])) i++;

    float lo = c->x[i], hi = c->x[i+1];
    for (int it = 0; it < 100; it++) {
        float mid = 0.5f * (lo + hi);
        if (fabsf(hi - lo) <= TONE_XTOL_REL * (1.0f + fabsf(mid))) break;
        if (mid == lo || mid == hi) break;
        if (tone_eval(c, mid) < y) lo = mid; else hi = mid;
    }
    *out = 0.5f * (lo + hi);
    return 1;
}

/* ------------------------------------------------------------------ */
/* analysis: monotonicity, non-monotone intervals, round-trip bounds   */

typedef struct {
    float lo, hi;        /* x-interval (equal for a single touch point) */
    int kind;            /* 0 = negative slope, 1 = zero-slope touch */
} NonMono;

typedef struct {
    int strict;          /* L > 0 everywhere                          */
    int invertible;      /* L >= 0 everywhere                         */
    double lmin, lmax;   /* slope range over [x0, x_{n-1}]            */
    int nnon;
    NonMono non[NONMONO_MAX];
    double x_tol, y_tol;
    double rtA_bound;    /* |f^-1(f(x)) - x|                          */
    double rtB_bound;    /* |f(f^-1(y)) - y|                          */
} Analysis;

static int quad_roots(double a, double b, double c, double* r1, double* r2) {
    if (fabs(a) < 1e-12) {
        if (fabs(b) < 1e-15) return 0;
        *r1 = *r2 = -c / b;
        return 1;
    }
    double d = b * b - 4.0 * a * c;
    if (d < 0.0) return 0;
    d = sqrt(d);
    double q = (-b - d) / (2.0 * a);
    double w = (-b + d) / (2.0 * a);
    if (q > w) { double t = q; q = w; w = t; }
    *r1 = q; *r2 = w;
    return 1;
}

static void analyse(const Curve* c, Analysis* a) {
    memset(a, 0, sizeof(*a));
    a->strict = 1;
    a->invertible = 1;
    a->lmin = 1e300;
    a->lmax = -1e300;

    for (int i = 0; i < c->n - 1; i++) {
        double gi = c->g[i], gj = c->g[i+1];
        double A = 3.0 * gi + 3.0 * gj - 6.0;
        double B = 6.0 - 4.0 * gi - 2.0 * gj;
        double C = gi;

        double um = -1.0;
        if (fabs(A) > 1e-15) um = -B / (2.0 * A);

        double lo = 1e300, hi = -1e300;
        double cand[3] = { 0.0, 1.0, um };
        int nc = (um > 0.0 && um < 1.0) ? 3 : 2;
        for (int k = 0; k < nc; k++) {
            double L = seg_slope(c, i, cand[k]);
            if (L < lo) lo = L;
            if (L > hi) hi = L;
        }
        if (lo < a->lmin) a->lmin = lo;
        if (hi > a->lmax) a->lmax = hi;

        double h = (double)c->x[i+1] - (double)c->x[i];
        double xlo = (double)c->x[i];
        if (lo < -1e-12) {
            a->strict = 0;
            a->invertible = 0;
            double r1, r2;
            if (quad_roots(A, B, C, &r1, &r2) && a->nnon < NONMONO_MAX) {
                double u0 = r1 > 0.0 ? r1 : 0.0;
                double u1 = r2 < 1.0 ? r2 : 1.0;
                if (u0 < u1) {
                    a->non[a->nnon].lo = (float)(xlo + h * u0);
                    a->non[a->nnon].hi = (float)(xlo + h * u1);
                    a->non[a->nnon].kind = 0;
                    a->nnon++;
                }
            }
        } else if (lo <= 1e-9) {
            /* slope reaches zero at the interior vertex: still increasing,
               but the inverse is not Lipschitz there. */
            a->strict = 0;
            if (um > 0.0 && um < 1.0 && a->nnon < NONMONO_MAX) {
                double xp = xlo + h * um;
                a->non[a->nnon].lo = (float)xp;
                a->non[a->nnon].hi = (float)xp;
                a->non[a->nnon].kind = 1;
                a->nnon++;
            }
        }
    }

    a->x_tol = TONE_XTOL_ABS;
    a->y_tol = TONE_YTOL_ABS;
    if (a->invertible && a->lmin > 0.0) {
        a->rtA_bound = a->x_tol + a->y_tol / a->lmin;
        a->rtB_bound = a->lmax * a->x_tol + a->y_tol;
    } else {
        a->rtA_bound = INFINITY;
        a->rtB_bound = INFINITY;
    }
}

/* ------------------------------------------------------------------ */
/* analytic inverse of the symmetric single-interval cubic            */

/* x such that  A u^3 + B u^2 + C u = y,  A=2g-2, B=3-3g, C=g */
static int cardano_inv(double y, double g, double* xout) {
    if (fabs(g - 1.0) < 1e-12) { *xout = y; return 1; }
    double A = 2.0 * g - 2.0, B = 3.0 - 3.0 * g, C = g;
    double shift = B / (3.0 * A);
    double p = (3.0 * A * C - B * B) / (3.0 * A * A);
    double q = (2.0 * B * B * B - 9.0 * A * B * C + 27.0 * A * A * (-y))
             / (27.0 * A * A * A);
    double disc = (q * 0.5) * (q * 0.5) + (p / 3.0) * (p / 3.0) * (p / 3.0);
    double u;
    if (disc >= 0.0) {
        double s = sqrt(disc);
        double v = cbrt(-0.5 * q + s) + cbrt(-0.5 * q - s);
        u = v - shift;
    } else {
        double r = 2.0 * sqrt(-p / 3.0);
        double th = acos((3.0 * q / (2.0 * p)) * sqrt(-3.0 / p));
        double best = 0.0, bd = 1e300;
        for (int k = 0; k < 3; k++) {
            double vk = r * cos((th - 2.0 * 3.14159265358979323846 * k) / 3.0);
            double uk = vk - shift;
            double d = (uk < 0.0 ? -uk : 0.0) + (uk > 1.0 ? uk - 1.0 : 0.0);
            if (d < bd) { bd = d; best = uk; }
        }
        u = best;
    }
    *xout = u;
    return 1;
}

/* ------------------------------------------------------------------ */
/* self-test                                                           */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static uint32_t rng_state = 0x12345678u;
static float urand(void) {
    rng_state = rng_state * 1664525u + 1013904223u;
    return (float)((rng_state >> 8) * (1.0 / 16777216.0));
}

typedef struct { double rtA, rtB; } RtErr;

static RtErr measure_round_trip(const Curve* c) {
    RtErr e = { 0.0, 0.0 };
    for (int k = 0; k <= 2000; k++) {
        float x = (float)k / 2000.0f;
        float xb;
        tone_inv(c, tone_eval(c, x), &xb);
        double d = fabs((double)xb - (double)x);
        if (d > e.rtA) e.rtA = d;

        float y = (float)k / 2000.0f;
        float xq;
        tone_inv(c, y, &xq);
        double d2 = fabs((double)tone_eval(c, xq) - (double)y);
        if (d2 > e.rtB) e.rtB = d2;
    }
    for (int k = 0; k < 2000; k++) {
        float x = urand();
        float xb;
        tone_inv(c, tone_eval(c, x), &xb);
        double d = fabs((double)xb - (double)x);
        if (d > e.rtA) e.rtA = d;

        float y = urand();
        float xq;
        tone_inv(c, y, &xq);
        double d2 = fabs((double)tone_eval(c, xq) - (double)y);
        if (d2 > e.rtB) e.rtB = d2;
    }
    return e;
}

static void print_analysis(const char* tag, const Analysis* a) {
    printf("%s.strict=%d\n", tag, a->strict);
    printf("%s.lmin=%.9g\n", tag, a->lmin);
    printf("%s.lmax=%.9g\n", tag, a->lmax);
    printf("%s.nonmono_count=%d\n", tag, a->nnon);
    for (int i = 0; i < a->nnon; i++)
        printf("%s.nonmono[%d]=[%.9g,%.9g] kind=%s\n", tag, i,
               a->non[i].lo, a->non[i].hi,
               a->non[i].kind == 0 ? "negative_slope" : "zero_slope_touch");
}

static void show_report(const char* tag, const ToneReport* r) {
    printf("%s.err=%s\n", tag, tone_err_name(r->err));
    printf("%s.field=%s\n", tag, r->field);
    printf("%s.index=%d\n", tag, r->index);
    printf("%s.value=%.9g\n", tag, r->value);
}

int main(void) {
    Curve c;
    Analysis a;

    /* ---- 1. identity curve: f(x) = x exactly ---- */
    memset(&c, 0, sizeof(c));
    c.n = 2; c.x[0] = 0.0f; c.x[1] = 1.0f;
    c.g[0] = 1.0f; c.g[1] = 1.0f;
    ToneReport r = tone_validate(&c);
    check("identity_valid", r.err == TONE_OK);
    float id_err = 0.0f;
    for (int k = 0; k <= 100; k++) {
        float x = (float)k / 100.0f;
        float d = fabsf(tone_eval(&c, x) - x);
        if (d > id_err) id_err = d;
    }
    /* identity: g=1 everywhere, so the Hermite reduces to u up to float
       rounding in the cubic expansion (not bit-exact). */
    check("identity_eval", id_err < 1e-6f);

    analyse(&c, &a);
    print_analysis("identity", &a);
    RtErr ie = measure_round_trip(&c);
    printf("identity.rtA_err=%.9g\n", ie.rtA);
    printf("identity.rtB_err=%.9g\n", ie.rtB);
    check("identity_round_trip", ie.rtA <= a.rtA_bound && ie.rtB <= a.rtB_bound);

    /* ---- 2. S-curve: reversible, strict monotone ---- */
    memset(&c, 0, sizeof(c));
    c.n = 5;
    c.x[0] = 0.00f; c.x[1] = 0.25f; c.x[2] = 0.50f; c.x[3] = 0.75f; c.x[4] = 1.00f;
    c.g[0] = 0.50f; c.g[1] = 0.80f; c.g[2] = 1.50f; c.g[3] = 1.20f; c.g[4] = 0.70f;
    r = tone_validate(&c);
    check("scurve_valid", r.err == TONE_OK);
    analyse(&c, &a);
    print_analysis("scurve", &a);
    check("scurve_strict", a.strict == 1);
    check("scurve_bounds_finite", isfinite(a.rtA_bound) && isfinite(a.rtB_bound));

    RtErr se = measure_round_trip(&c);
    printf("scurve.rtA_bound=%.9g\n", a.rtA_bound);
    printf("scurve.rtA_err=%.9g\n", se.rtA);
    printf("scurve.rtB_bound=%.9g\n", a.rtB_bound);
    printf("scurve.rtB_err=%.9g\n", se.rtB);
    check("scurve_rtA_within_bound", se.rtA <= a.rtA_bound);
    check("scurve_rtB_within_bound", se.rtB <= a.rtB_bound);

    /* expansion/compression actually happened near the knots */
    check("scurve_has_slope_spread", a.lmax > 1.0 && a.lmin < 1.0);

    /* ---- 3. non-monotone curve: report intervals ---- */
    memset(&c, 0, sizeof(c));
    c.n = 3;
    c.x[0] = 0.0f; c.x[1] = 0.5f; c.x[2] = 1.0f;
    c.g[0] = 1.0f; c.g[1] = 6.0f; c.g[2] = 1.0f;
    r = tone_validate(&c);
    check("nonmono_valid_params", r.err == TONE_OK);
    analyse(&c, &a);
    print_analysis("nonmono", &a);
    /* slope 15u^2-10u+1 has roots ~0.1225 and ~0.5441 per half-interval,
       so expect two reversal intervals around x=0.0613..0.2721 and
       0.7279..0.9387 */
    check("nonmono_detected", a.invertible == 0 && a.strict == 0);
    check("nonmono_two_intervals", a.nnon == 2);
    if (a.nnon == 2) {
        check("nonmono_interval_lo",
              fabsf(a.non[0].lo - 0.06125f) < 2e-3f &&
              fabsf(a.non[0].hi - 0.27205f) < 2e-3f);
        check("nonmono_interval_hi",
              fabsf(a.non[1].lo - 0.72795f) < 2e-3f &&
              fabsf(a.non[1].hi - 0.93875f) < 2e-3f);
        check("nonmono_kind", a.non[0].kind == 0 && a.non[1].kind == 0);
    }
    check("nonmono_bounds_infinite", isinf(a.rtA_bound) && isinf(a.rtB_bound));

    /* ---- 4. zero-slope touch: g=3 single interval, monotone but marginal ---- */
    memset(&c, 0, sizeof(c));
    c.n = 2; c.x[0] = 0.0f; c.x[1] = 1.0f;
    c.g[0] = 3.0f; c.g[1] = 3.0f;
    analyse(&c, &a);
    print_analysis("touch", &a);
    check("touch_not_strict", a.strict == 0);
    check("touch_still_invertible", a.invertible == 1);
    check("touch_reported", a.nnon == 1 && a.non[0].kind == 1 &&
                            fabsf(a.non[0].lo - 0.5f) < 1e-3f);

    /* ---- 5. parameter violations report the offending value ---- */
    struct { Curve c; ToneErr want; int idx; float val; const char* tag; } bad[7];
    memset(bad, 0, sizeof(bad));

    bad[0].tag = "viol_n_small"; bad[0].c.n = 1;
    bad[0].want = TONE_E_N; bad[0].idx = -1; bad[0].val = 1.0f;

    bad[1].tag = "viol_x_order";
    bad[1].c.n = 3;
    bad[1].c.x[0] = 0.0f; bad[1].c.x[1] = 0.5f; bad[1].c.x[2] = 0.4f;
    bad[1].c.g[0] = bad[1].c.g[1] = bad[1].c.g[2] = 1.0f;
    bad[1].want = TONE_E_X_ORDER; bad[1].idx = 2; bad[1].val = 0.4f;

    bad[2].tag = "viol_x_range";
    bad[2].c.n = 3;
    bad[2].c.x[0] = 0.0f; bad[2].c.x[1] = 0.5f; bad[2].c.x[2] = 1.2f;
    bad[2].c.g[0] = bad[2].c.g[1] = bad[2].c.g[2] = 1.0f;
    bad[2].want = TONE_E_X_RANGE; bad[2].idx = 2; bad[2].val = 1.2f;

    bad[3].tag = "viol_x_nan";
    bad[3].c.n = 2;
    bad[3].c.x[0] = 0.0f; bad[3].c.x[1] = NAN;
    bad[3].c.g[0] = bad[3].c.g[1] = 1.0f;
    bad[3].want = TONE_E_X_NAN; bad[3].idx = 1;

    bad[4].tag = "viol_g_zero";
    bad[4].c.n = 2;
    bad[4].c.x[0] = 0.0f; bad[4].c.x[1] = 1.0f;
    bad[4].c.g[0] = 1.0f; bad[4].c.g[1] = 0.0f;
    bad[4].want = TONE_E_G_RANGE; bad[4].idx = 1; bad[4].val = 0.0f;

    bad[5].tag = "viol_g_negative";
    bad[5].c.n = 3;
    bad[5].c.x[0] = 0.0f; bad[5].c.x[1] = 0.5f; bad[5].c.x[2] = 1.0f;
    bad[5].c.g[0] = 1.0f; bad[5].c.g[1] = -2.0f; bad[5].c.g[2] = 1.0f;
    bad[5].want = TONE_E_G_RANGE; bad[5].idx = 1; bad[5].val = -2.0f;

    bad[6].tag = "viol_g_nan";
    bad[6].c.n = 2;
    bad[6].c.x[0] = 0.0f; bad[6].c.x[1] = 1.0f;
    bad[6].c.g[0] = INFINITY; bad[6].c.g[1] = 1.0f;
    bad[6].want = TONE_E_G_NAN; bad[6].idx = 0;

    for (int i = 0; i < 7; i++) {
        ToneReport br = tone_validate(&bad[i].c);
        show_report(bad[i].tag, &br);
        int ok = br.err == bad[i].want && br.index == bad[i].idx;
        if (bad[i].want != TONE_E_X_NAN && bad[i].want != TONE_E_G_NAN)
            ok = ok && br.value == bad[i].val;
        check(bad[i].tag, ok);
    }

    /* ---- 6. extrapolation is inverted exactly ---- */
    memset(&c, 0, sizeof(c));
    c.n = 3;
    c.x[0] = 0.2f; c.x[1] = 0.5f; c.x[2] = 0.8f;
    c.g[0] = 0.7f; c.g[1] = 1.3f; c.g[2] = 0.9f;
    analyse(&c, &a);
    int ext_ok = 1;
    float xs[4] = { -0.5f, 0.0f, 1.0f, 1.5f };
    for (int i = 0; i < 4; i++) {
        float xb;
        tone_inv(&c, tone_eval(&c, xs[i]), &xb);
        if (fabsf(xb - xs[i]) > 1e-5f) ext_ok = 0;
    }
    check("extrapolation_inverse", ext_ok);
    printf("ext.endpoints=(%.9g,%.9g)\n",
           (double)tone_eval(&c, c.x[0]), (double)tone_eval(&c, c.x[2]));

    /* ---- 7. numerical inverse vs analytic (Cardano) ---- */
    double gs[3] = { 0.05, 0.5, 2.0 };
    const char* gtags[3] = { "g0p05", "g0p5", "g2" };
    for (int gi = 0; gi < 3; gi++) {
        memset(&c, 0, sizeof(c));
        c.n = 2; c.x[0] = 0.0f; c.x[1] = 1.0f;
        c.g[0] = (float)gs[gi]; c.g[1] = (float)gs[gi];
        analyse(&c, &a);
        double maxe = 0.0;
        for (int k = 0; k <= 2000; k++) {
            float y = (float)k / 2000.0f;
            float xn;
            double xa;
            tone_inv(&c, y, &xn);
            cardano_inv((double)y, gs[gi], &xa);
            double d = fabs((double)xn - xa);
            if (d > maxe) maxe = d;
        }
        printf("analytic_%s_maxerr=%.9g\n", gtags[gi], maxe);
        char nm[64];
        snprintf(nm, sizeof(nm), "analytic_%s_match", gtags[gi]);
        check(nm, maxe < 1e-5);
        check("analytic_monotone", a.strict == 1);
    }

    /* ---- 8. RNG-independent determinism of the inverse ---- */
    float x1, x2;
    tone_inv(&c, 0.25f, &x1);
    tone_inv(&c, 0.25f, &x2);
    check("inverse_deterministic", x1 == x2);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
