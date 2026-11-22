
/* bounce: iterative multi-bounce indirect light for a two-zone enclosure.

   A scene holds direct emission E and a transport matrix T (Lambert albedo
   times form factor); zone 0 is the room corner, zone 1 the walls. Bounce k
   adds T^k E to the indirect accumulator, so after N bounces
       L_N = E + sum_{k=1..N} T^k E.
   With q = ||T||_inf < 1, the geometric-series tail
       r_N = ||E||_inf * q^(N+1) / (1 - q)
   bounds the error against the analytic infinite-bounce solution
       L* = (I - T)^-1 E,
   shrinks by a factor q on every bounce, and drives the stop test. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define ZONES 2
#define BOUNCE_MAX 64

/* return codes */
#define BOUNCE_OK 0
#define BOUNCE_ERR_ZERO -1        /* max_bounces == 0: rejected */
#define BOUNCE_ERR_NONCONTRACT -2 /* ||T||_inf >= 1, energy does not decay */
#define BOUNCE_ERR_PARAM -3
#define BOUNCE_ERR_BUDGET -4      /* tolerance not reached within max_bounces */

typedef struct {
    double T[ZONES][ZONES]; /* transport: albedo * form factor */
    double E[ZONES];        /* direct emission, bounce 0 */
} Scene;

typedef struct {
    double L[ZONES];                       /* E + indirect after k bounces */
    double indirect[ZONES];                /* sum of bounce terms T^k E */
    double term[ZONES];                    /* latest bounce contribution */
    double residual;                       /* tail-energy upper bound */
    double err;                            /* actual error vs analytic L* */
    int bounces;
    int stopped;                           /* residual <= tol */
    double r_hist[BOUNCE_MAX + 1];         /* residual per bounce, r_hist[0] = k0 */
    double e_hist[BOUNCE_MAX + 1];         /* actual error per bounce */
    double term_hist[BOUNCE_MAX + 1];      /* ||bounce term||_inf per bounce */
    int hist;
} BounceResult;

static double v_norm_inf(const double v[ZONES]) {
    double m = 0.0;
    for (int i = 0; i < ZONES; i++) {
        double a = fabs(v[i]);
        if (a > m) m = a;
    }
    return m;
}

static double inf_dist(const double a[ZONES], const double b[ZONES]) {
    double d[ZONES];
    for (int i = 0; i < ZONES; i++) d[i] = a[i] - b[i];
    return v_norm_inf(d);
}

/* contraction factor q = ||T||_inf: largest absolute row sum */
static double scene_q(const Scene* s) {
    double q = 0.0;
    for (int i = 0; i < ZONES; i++) {
        double row = 0.0;
        for (int j = 0; j < ZONES; j++) row += fabs(s->T[i][j]);
        if (row > q) q = row;
    }
    return q;
}

static int valid_scene(const Scene* s) {
    for (int i = 0; i < ZONES; i++) {
        if (s->E[i] < 0.0) return 0;
        for (int j = 0; j < ZONES; j++)
            if (s->T[i][j] < 0.0) return 0;
    }
    return 1;
}

/* analytic infinite-bounce radiance L* = (I - T)^-1 E, explicit 2x2 inverse */
static int bounce_analytic(const Scene* s, double Lstar[ZONES]) {
    double a = 1.0 - s->T[0][0], b = -s->T[0][1];
    double c = -s->T[1][0],       d = 1.0 - s->T[1][1];
    double det = a * d - b * c;
    if (fabs(det) < 1e-15) return -1;
    double inv = 1.0 / det;
    Lstar[0] = inv * (d * s->E[0] - b * s->E[1]);
    Lstar[1] = inv * (-c * s->E[0] + a * s->E[1]);
    return 0;
}

/* smallest k >= 1 for which e0 * q^(k+1) / (1 - q) <= tol */
static int min_bounces(double q, double e0, double tol) {
    if (e0 <= 0.0 || q <= 0.0) return 1;
    if (tol >= e0 * q / (1.0 - q)) return 1; /* tolerance already covers all indirect */
    double k = log(tol * (1.0 - q) / e0) / log(q) - 1.0;
    int n = (int)ceil(k);
    return n < 1 ? 1 : n;
}

static int bounce_iterate(const Scene* s, int max_bounces, double tol, BounceResult* out) {
    memset(out, 0, sizeof(*out));
    double q = scene_q(s);
    double e0 = v_norm_inf(s->E);

    /* zero bounces is never accepted: report the minimum usable count instead */
    if (max_bounces == 0) {
        out->bounces = (tol > 0.0 && q < 1.0 && valid_scene(s))
                           ? min_bounces(q, e0, tol) : 1;
        return BOUNCE_ERR_ZERO;
    }
    if (max_bounces < 0 || tol <= 0.0 || !valid_scene(s)) return BOUNCE_ERR_PARAM;
    if (q >= 1.0) return BOUNCE_ERR_NONCONTRACT;
    if (max_bounces > BOUNCE_MAX) max_bounces = BOUNCE_MAX;

    double Lstar[ZONES];
    bounce_analytic(s, Lstar); /* q < 1 guarantees I - T is invertible */

    for (int i = 0; i < ZONES; i++) {
        out->L[i] = s->E[i];
        out->term[i] = s->E[i];
    }
    double residual = e0 * q / (1.0 - q); /* k = 0 tail: all indirect still missing */
    out->residual = residual;
    out->r_hist[0] = residual;
    out->term_hist[0] = e0;
    out->e_hist[0] = inf_dist(Lstar, s->E);
    out->hist = 1;

    for (int k = 1; k <= max_bounces; k++) {
        double next[ZONES] = { 0.0, 0.0 };
        for (int i = 0; i < ZONES; i++)
            for (int j = 0; j < ZONES; j++)
                next[i] += s->T[i][j] * out->term[j]; /* T * T^(k-1) E = T^k E */

        for (int i = 0; i < ZONES; i++) {
            out->term[i] = next[i];
            out->indirect[i] += next[i];              /* accumulate indirect per bounce */
            out->L[i] = s->E[i] + out->indirect[i];
        }

        residual *= q; /* e0 q^(k+1)/(1-q): geometric tail after one more bounce */
        out->residual = residual;
        out->r_hist[k] = residual;
        out->term_hist[k] = v_norm_inf(next);
        out->e_hist[k] = inf_dist(Lstar, out->L);
        out->hist = k + 1;
        out->bounces = k;

        if (residual <= tol) {
            out->stopped = 1;
            break;
        }
    }
    out->err = out->e_hist[out->hist - 1];
    return out->stopped ? BOUNCE_OK : BOUNCE_ERR_BUDGET;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* corner (zone 0) and walls (zone 1); light mostly reaches the walls */
    Scene room = {
        { { 0.30, 0.25 },
          { 0.20, 0.35 } },
        { 0.20, 1.00 }
    };
    double q = scene_q(&room);
    double e0 = v_norm_inf(room.E);
    check("contractive", fabs(q - 0.55) < 1e-12 && q < 1.0);

    /* zero bounces rejected, minimum usable count reported */
    BounceResult r0;
    check("zero_rejected", bounce_iterate(&room, 0, 1e-6, &r0) == BOUNCE_ERR_ZERO);
    check("min_usable_ge_one", r0.bounces >= 1);

    /* closed-form minimum matches a brute-force search, and is truly minimal */
    int nmin = min_bounces(q, e0, 1e-6);
    int brute = 1;
    while (e0 * pow(q, brute + 1) / (1.0 - q) > 1e-6) brute++;
    check("min_formula_matches_brute", nmin == brute && nmin > 1);
    check("zero_hint_is_min", r0.bounces == nmin);

    BounceResult rn, rm;
    check("min_meets_tol",
          bounce_iterate(&room, nmin, 1e-6, &rn) == BOUNCE_OK &&
          rn.stopped == 1 && rn.bounces == nmin && rn.residual <= 1e-6);
    check("min_minus1_misses",
          bounce_iterate(&room, nmin - 1, 1e-6, &rm) == BOUNCE_ERR_BUDGET &&
          rm.residual > 1e-6);

    /* full run: analytic comparison, error upper bound, monotone residual */
    BounceResult r;
    check("converged", bounce_iterate(&room, BOUNCE_MAX, 1e-9, &r) == BOUNCE_OK);

    double Lstar[ZONES];
    bounce_analytic(&room, Lstar);
    check("analytic_values",
          fabs(Lstar[0] - 0.38 / 0.405) < 1e-10 &&
          fabs(Lstar[1] - 0.74 / 0.405) < 1e-10);

    int bound_ok = 1, mono_ok = 1, contract_ok = 1, split_ok = 1;
    for (int k = 0; k < r.hist; k++) {
        if (r.e_hist[k] > r.r_hist[k] * (1.0 + 1e-12)) bound_ok = 0;       /* error <= tail bound */
        if (k > 0 && !(r.r_hist[k] < r.r_hist[k - 1])) mono_ok = 0;        /* strictly decreasing */
        if (k > 0 && r.term_hist[k] > q * r.term_hist[k - 1] * (1.0 + 1e-12)) contract_ok = 0;
    }
    for (int i = 0; i < ZONES; i++)
        if (fabs(r.L[i] - (room.E[i] + r.indirect[i])) > 1e-15) split_ok = 0;
    check("error_upper_bound", bound_ok);
    check("residual_monotone", mono_ok);
    check("bounce_terms_contract", contract_ok);
    check("total_is_direct_plus_indirect", split_ok);
    check("final_error_within_tol", r.err <= 2e-9);

    /* recurrence residual matches the closed-form geometric tail */
    double tail = e0 * pow(q, r.bounces + 1) / (1.0 - q);
    check("residual_closed_form", fabs(r.residual - tail) <= tail * 1e-12);

    /* uniform enclosure T = rho I: indirect is the scalar series E rho/(1-rho) */
    Scene uni = {
        { { 0.5, 0.0 },
          { 0.0, 0.5 } },
        { 1.0, 2.0 }
    };
    BounceResult u;
    check("uniform_converged", bounce_iterate(&uni, BOUNCE_MAX, 1e-12, &u) == BOUNCE_OK);
    check("uniform_series",
          fabs(u.indirect[0] - 1.0) < 1e-9 && fabs(u.indirect[1] - 2.0) < 1e-9 &&
          fabs(u.L[0] - 2.0) < 1e-9 && fabs(u.L[1] - 4.0) < 1e-9);

    /* no transport: one bounce settles everything, residual hits zero */
    Scene dead = {
        { { 0.0, 0.0 },
          { 0.0, 0.0 } },
        { 1.0, 1.0 }
    };
    BounceResult d;
    check("dead_stops_at_one",
          bounce_iterate(&dead, BOUNCE_MAX, 1e-9, &d) == BOUNCE_OK &&
          d.bounces == 1 && d.stopped == 1 && d.residual == 0.0 &&
          d.indirect[0] == 0.0 && d.indirect[1] == 0.0);
    check("dead_min_is_one", min_bounces(0.0, 1.0, 1e-9) == 1);

    /* parameter validation */
    Scene bad = {
        { { 0.6, 0.5 },
          { 0.4, 0.7 } },
        { 1.0, 1.0 }
    };
    BounceResult tmp;
    check("noncontract_rejected",
          bounce_iterate(&bad, 4, 1e-6, &tmp) == BOUNCE_ERR_NONCONTRACT);
    check("bad_tol_rejected", bounce_iterate(&room, 4, 0.0, &tmp) == BOUNCE_ERR_PARAM);
    check("negative_bounces_rejected",
          bounce_iterate(&room, -1, 1e-6, &tmp) == BOUNCE_ERR_PARAM);

    /* budget exhaustion reports the best estimate instead of silently passing */
    BounceResult lim;
    check("budget_exceeded",
          bounce_iterate(&room, 3, 1e-12, &lim) == BOUNCE_ERR_BUDGET &&
          lim.stopped == 0 && lim.bounces == 3 && lim.residual > 1e-12);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
