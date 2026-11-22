
/* medium: participating medium with emission + absorption + scattering.
   - free-flight distance sampled from extinction coefficient sigma_t
   - direction sampled from the Henyey-Greenstein phase function
   - slab path tracer with next-event estimation, energy-conserving weights
   - Monte Carlo result checked against the analytic single-scattering
     solution, with a Hoeffding high-probability error upper bound. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265358979323846

typedef struct { double x, y, z; } V3;

static V3 v3(double x, double y, double z) { V3 r = { x, y, z }; return r; }
static double dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 norm3(V3 a) {
    double l = sqrt(dot3(a, a));
    return l > 0 ? v3(a.x / l, a.y / l, a.z / l) : v3(0, 0, 0);
}

/* medium coefficients: absorption, scattering, emission */
typedef struct { double sa, ss, se; } Medium;

static double medium_st(const Medium* m) { return m->sa + m->ss; }
static double medium_albedo(const Medium* m) {
    double st = medium_st(m);
    return st > 0 ? m->ss / st : 0.0;
}
static int medium_valid(const Medium* m) {
    return m->sa >= 0 && m->ss >= 0 && m->se >= 0 && medium_st(m) > 0;
}

/* directional beam along +z */
typedef struct { double I0; } Beam;

/* ---------- rng ---------- */
typedef struct { unsigned long long s; } Rng;

static double rng_u01(Rng* r) {
    r->s ^= r->s << 13; r->s ^= r->s >> 7; r->s ^= r->s << 17;
    return (double)(r->s >> 11) / (double)(1ULL << 53);
}

/* ---------- transmittance and free-flight sampling ---------- */
static double transmittance(double st, double d) { return exp(-st * d); }

/* sample distance to the next interaction in a homogeneous medium:
   inverse CDF of p(t) = sigma_t * T(t). */
static double sample_distance(Rng* r, double st) {
    double u = rng_u01(r);
    return -log(1.0 - u) / st;
}

/* returns 1 for an in-slab collision, 0 for escape through the far face */
static int sample_free_flight(Rng* r, double st, double L, double* t) {
    *t = sample_distance(r, st);
    return *t < L;
}

/* ---------- Henyey-Greenstein phase function ---------- */
static double hg_pdf(double cos_t, double g) {
    double d = 1.0 + g * g - 2.0 * g * cos_t;
    return (1.0 - g * g) / (4.0 * PI * sqrt(d) * d);
}

/* sample cos(theta); inverse-CDF of the HG phase function */
static double hg_sample_cos(Rng* r, double g) {
    double u = rng_u01(r);
    if (fabs(g) < 1e-3) return 2.0 * u - 1.0;
    double a = (1.0 - g * g) / (1.0 - g + 2.0 * g * u);
    return (1.0 + g * g - a * a) / (2.0 * g);
}

/* build an orthonormal frame around n and turn a phase sample into a
   world-space direction (xi in n's local frame, z = n) */
static V3 hg_sample_world(Rng* r, V3 n, double g) {
    double cos_t = hg_sample_cos(r, g);
    double sin_t = sqrt(fmax(0.0, 1.0 - cos_t * cos_t));
    double phi = 2.0 * PI * rng_u01(r);
    V3 helper = fabs(n.z) < 0.9 ? v3(0, 0, 1) : v3(1, 0, 0);
    V3 tx = norm3((V3){ helper.y * n.z - helper.z * n.y,
                        helper.z * n.x - helper.x * n.z,
                        helper.x * n.y - helper.y * n.x });
    V3 ty = (V3){ n.y * tx.z - n.z * tx.y,
                  n.z * tx.x - n.x * tx.z,
                  n.x * tx.y - n.y * tx.x };
    V3 d = v3(sin_t * cos(phi) * tx.x + sin_t * sin(phi) * ty.x + cos_t * n.x,
              sin_t * cos(phi) * tx.y + sin_t * sin(phi) * ty.y + cos_t * n.y,
              sin_t * cos(phi) * tx.z + sin_t * sin(phi) * ty.z + cos_t * n.z);
    return norm3(d);
}

/* ---------- analytic single-scattering solution for a slab ----------
   Camera ray enters at z=0 and travels +z through [0, L]; the beam
   travels the same direction, so at any collision cos(theta) = 1.
   L = integral [se + ss * I0 * T(t) * p] * T(t) dt + I0 * T(L). */
static double analytic_emission(const Medium* m, double L) {
    double st = medium_st(m);
    return m->se / st * (1.0 - transmittance(st, L));
}

static double analytic_single_scatter(const Medium* m, const Beam* b,
                                       double L, double g) {
    double st = medium_st(m);
    double p1 = hg_pdf(1.0, g);
    return b->I0 * m->ss * p1 * (1.0 - transmittance(2.0 * st, L)) / (2.0 * st);
}

static double analytic_total(const Medium* m, const Beam* b, double L,
                             double g) {
    double st = medium_st(m);
    return analytic_emission(m, L)
         + analytic_single_scatter(m, b, L, g)
         + b->I0 * transmittance(st, L);
}

/* mean emission depth: where the glow sits along the ray */
static double analytic_emission_centroid(const Medium* m, double L) {
    double st = medium_st(m);
    double TL = transmittance(st, L);
    return (1.0 / st - L * TL / (1.0 - TL));
}

/* ---------- Monte Carlo slab estimator (next-event estimation) ----------
   One path, driven by one uniform u for the free-flight draw.
   View transmittance is implicit in the event probability; the beam's
   own transmittance to the collision point stays explicit. */
typedef struct { double total, emission, depth; int hit; } PathResult;

static PathResult trace_slab(const Medium* m, const Beam* b, double L,
                             double g, double u) {
    PathResult r = { 0, 0, 0, 0 };
    double st = medium_st(m);
    double t = -log(1.0 - u) / st;
    if (t >= L) { r.total = b->I0; return r; }      /* escaped: T implicit */
    r.hit = 1;
    r.depth = t;
    r.emission = m->se / st;                        /* collision estimator */
    double scatter = (m->ss / st) * b->I0 * transmittance(st, t)
                   * hg_pdf(1.0, g);
    r.total = r.emission + scatter;
    return r;
}

/* ---------- running statistics (Welford) ---------- */
typedef struct { double mean, m2, var; int n; } Stat;

static void stat_init(Stat* s) { memset(s, 0, sizeof(*s)); }
static void stat_push(Stat* s, double x) {
    s->n++;
    double d = x - s->mean;
    s->mean += d / s->n;
    s->m2 += d * (x - s->mean);
    s->var = s->n > 1 ? s->m2 / (s->n - 1) : 0.0;
}
static double stat_se(const Stat* s) { return sqrt(s->var / s->n); }

/* Hoeffding bound for iid samples in [0, bound]: with probability
   1-delta the sample mean is within epsilon of the expectation. */
static double hoeffding_eps(double bound, int n, double delta) {
    return bound * sqrt(log(2.0 / delta) / (2.0 * n));
}

/* ---------- photon tracer: explicit energy-budget check ---------- */
typedef struct { double absorbed, transmitted, reflected; } PhotonBudget;

static PhotonBudget photon_trace(const Medium* m, int n, double L, double g,
                                 Rng* r) {
    PhotonBudget pb = { 0, 0, 0 };
    double st = medium_st(m);
    for (int i = 0; i < n; i++) {
        double w = 1.0, z = 0.0;
        V3 d = v3(0, 0, 1);
        for (;;) {
            double t = sample_distance(r, st);
            double nz = z + d.z * t;
            if (nz >= L) { pb.transmitted += w; break; }
            if (nz <= 0.0) { pb.reflected += w; break; }
            z = nz;
            pb.absorbed += w * (m->sa / st);   /* fractions sum to 1 */
            w *= m->ss / st;                   /* albedo: never grows */
            if (w == 0.0) break;               /* fully absorbed */
            d = hg_sample_world(r, d, g);
        }
    }
    return pb;
}

/* ---------- test harness ---------- */
static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* Simpson quadrature of hg_pdf over cos in [-1, 1] */
static double hg_integral(double g) {
    int n = 1000; double h = 2.0 / n, s = hg_pdf(-1, g) + hg_pdf(1, g);
    for (int i = 1; i < n; i++) {
        double x = -1.0 + i * h;
        s += (i % 2 ? 4.0 : 2.0) * hg_pdf(x, g);
    }
    return s * h / 3.0; /* should be 1/(2 pi): normalized over dOmega */
}

int main(void) {
    Rng rng = { 0x123456789abcdefULL };
    Medium m = { 0.4, 0.6, 0.3 };       /* sa, ss, se */
    Beam beam = { 1.0 };
    double L = 2.0, g = 0.0;
    double st = medium_st(&m);

    check("medium_valid", medium_valid(&m));
    Medium bad = { -0.1, 0.6, 0.3 };
    check("medium_rejects_negative", !medium_valid(&bad));
    check("albedo_le_one", medium_albedo(&m) <= 1.0 + 1e-12);

    check("transmittance_zero_dist", fabs(transmittance(st, 0) - 1.0) < 1e-15);
    check("transmittance_decays", transmittance(st, L) < transmittance(st, 1.0));

    /* free-flight draw follows p(t) = st T(t): P(t < L) = 1 - T(L) */
    int nff = 20000, hit = 0;
    for (int i = 0; i < nff; i++) { double t; hit += sample_free_flight(&rng, st, L, &t); }
    double p_hit = (double)hit / nff, p_exact = 1.0 - transmittance(st, L);
    double se_ff = sqrt(p_exact * (1 - p_exact) / nff);
    printf("free_flight_hit: mc=%.5f analytic=%.5f\n", p_hit, p_exact);
    check("free_flight_matches_pdf", fabs(p_hit - p_exact) < 3.0 * se_ff);

    /* phase function: normalization, sampling mean, pdf/histogram agreement */
    check("hg_normalized", fabs(2.0 * PI * hg_integral(0.3) - 1.0) < 1e-9);
    check("hg_isotropic_g0", fabs(hg_pdf(0.7, 0.0) - 1.0 / (4 * PI)) < 1e-15);

    int nph = 8000; double sum_cos = 0;
    int bins = 16; long hist[16] = { 0 };
    for (int i = 0; i < nph; i++) {
        double c = hg_sample_cos(&rng, 0.3);
        sum_cos += c;
        int b = (int)((c + 1.0) * 0.5 * bins);
        if (b == bins) b--;
        hist[b]++;
    }
    check("hg_sample_mean_g", fabs(sum_cos / nph - 0.3) < 0.05);
    int hist_ok = 1;
    for (int b = 0; b < bins; b++) {
        double mid = -1.0 + (b + 0.5) * 2.0 / bins;
        double mass_exact = 2.0 * PI * hg_pdf(mid, 0.3) * (2.0 / bins);
        double mass_mc = (double)hist[b] / nph;
        if (fabs(mass_mc - mass_exact) > 0.02) hist_ok = 0;
    }
    check("hg_histogram_matches_pdf", hist_ok);

    V3 axis = v3(0, 0, 1), mean_d = v3(0, 0, 0);
    int frame_ok = 1;
    for (int i = 0; i < nph; i++) {
        V3 d = hg_sample_world(&rng, axis, 0.3);
        mean_d.x += d.x; mean_d.y += d.y; mean_d.z += d.z;
        if (fabs(dot3(d, d) - 1.0) > 1e-12) frame_ok = 0;
    }
    mean_d.x /= nph; mean_d.y /= nph; mean_d.z /= nph;
    check("world_sample_unit", frame_ok);
    check("world_sample_aligned", fabs(mean_d.z - 0.3) < 0.05
                                  && fabs(mean_d.x) < 0.05 && fabs(mean_d.y) < 0.05);

    /* analytic references */
    double a_em = analytic_emission(&m, L);
    double a_sc = analytic_single_scatter(&m, &beam, L, g);
    double a_tot = analytic_total(&m, &beam, L, g);
    printf("analytic: emission=%.6f single_scatter=%.6f total=%.6f\n",
           a_em, a_sc, a_tot);
    Medium pure_abs = { 1.0, 0.0, 0.0 };
    check("analytic_no_scatter", analytic_single_scatter(&pure_abs, &beam, L, g) == 0.0);

    /* Monte Carlo against the analytic solution */
    int N = 16384;
    Stat st_tot, st_em; stat_init(&st_tot); stat_init(&st_em);
    double depth_num = 0, depth_den = 0, bound = 0.0;
    double safe_bound = beam.I0 + m.se / st
                      + (m.ss / st) * beam.I0 * hg_pdf(1.0, g);
    for (int i = 0; i < N; i++) {
        PathResult pr = trace_slab(&m, &beam, L, g, rng_u01(&rng));
        stat_push(&st_tot, pr.total);
        stat_push(&st_em, pr.emission);
        if (pr.hit) { depth_num += pr.depth; depth_den += 1.0; }
        if (pr.total > bound) bound = pr.total;
    }
    double err = fabs(st_tot.mean - a_tot);
    double eps3 = 3.0 * stat_se(&st_tot);
    double eps_h = hoeffding_eps(safe_bound, N, 0.01);
    printf("mc_total=%.6f err=%.6f se=%.6f 3se=%.6f hoeffding99=%.6f (bound=%.4f)\n",
           st_tot.mean, err, stat_se(&st_tot), eps3, eps_h, safe_bound);
    check("mc_within_3sigma", err < eps3);
    check("mc_within_hoeffding", err < eps_h && bound <= safe_bound + 1e-12);
    check("mc_emission_matches", fabs(st_em.mean - a_em) < 3.0 * stat_se(&st_em));

    /* convergence: stratified free-flight draws, errors must shrink */
    int Ns[] = { 1024, 4096, 16384 };
    double prev_err = 1e9, prev_se = 1e9; int conv_ok = 1;
    for (int k = 0; k < 3; k++) {
        Stat s; stat_init(&s);
        for (int i = 0; i < Ns[k]; i++) {
            double u = ((double)i + 0.5) / Ns[k];
            stat_push(&s, trace_slab(&m, &beam, L, g, u).total);
        }
        double e = fabs(s.mean - a_tot), se = stat_se(&s);
        printf("convergence N=%5d mean=%.6f err=%.2e se=%.2e\n",
               Ns[k], s.mean, e, se);
        if (!(e <= prev_err && se <= prev_se + 1e-15)) conv_ok = 0;
        prev_err = e; prev_se = se;
    }
    check("converges_with_samples", conv_ok);

    /* glow position: mean emission depth vs analytic centroid */
    double mc_centroid = depth_num / depth_den;
    double a_centroid = analytic_emission_centroid(&m, L);
    printf("emission_centroid: mc=%.5f analytic=%.5f\n", mc_centroid, a_centroid);
    check("emission_centroid", fabs(mc_centroid - a_centroid) < 0.02);

    /* zero emission must be reported, not silently rendered as glow */
    Medium m0 = { 0.4, 0.6, 0.0 };
    check("zero_emission_flagged", m0.se == 0.0);
    double a0 = analytic_total(&m0, &beam, L, g);
    Stat s0; stat_init(&s0);
    double em0 = 0.0;
    for (int i = 0; i < N; i++) {
        PathResult pr0 = trace_slab(&m0, &beam, L, g, rng_u01(&rng));
        stat_push(&s0, pr0.total);
        em0 += pr0.emission;
    }
    printf("emission_coeff=0: reported no glow, emission=%.1f (analytic %.1f), "
           "total mc=%.6f analytic=%.6f\n", em0, analytic_emission(&m0, L),
           s0.mean, a0);
    check("zero_emission_contribution", em0 == 0.0 && analytic_emission(&m0, L) == 0.0);
    check("zero_emission_total", fabs(s0.mean - a0) < 3.0 * stat_se(&s0));

    /* vacuum: sigma_t = 0 has no interaction to sample; beam passes through */
    Medium vac = { 0.0, 0.0, 0.0 };
    check("vacuum_reported", !medium_valid(&vac));
    printf("sigma_t=0: no scattering/absorption possible, radiance = I0 = %.1f\n",
           beam.I0);

    /* energy conservation: photon budget partitions exactly */
    PhotonBudget pb = photon_trace(&m, 20000, L, 0.3, &rng);
    double sum = pb.absorbed + pb.transmitted + pb.reflected;
    printf("photon_budget: absorb=%.4f transmit=%.4f reflect=%.4f sum=%.6f\n",
           pb.absorbed / 20000, pb.transmitted / 20000,
           pb.reflected / 20000, sum / 20000);
    check("energy_conservation", fabs(sum / 20000 - 1.0) < 1e-12);

    PhotonBudget pb0 = photon_trace(&m0, 20000, L, 0.3, &rng);
    double sum0 = pb0.absorbed + pb0.transmitted + pb0.reflected;
    check("energy_conservation_no_emission", fabs(sum0 / 20000 - 1.0) < 1e-12);

    PhotonBudget pba = photon_trace(&pure_abs, 20000, L, 0.0, &rng);
    double se_abs = sqrt((1.0 - transmittance(1.0, L)) * transmittance(1.0, L)
                         / 20000.0);
    check("pure_absorption_partition",
          fabs(pba.transmitted / 20000.0 - transmittance(1.0, L)) < 3.0 * se_abs
          && fabs(pba.absorbed / 20000.0 - (1.0 - transmittance(1.0, L))) < 3.0 * se_abs
          && pba.reflected == 0.0);

    Medium mna = { 0.0, 1.0, 0.2 };
    PhotonBudget pbn = photon_trace(&mna, 20000, L, 0.3, &rng);
    check("non_absorbing_no_loss", pbn.absorbed == 0.0
          && fabs((pbn.transmitted + pbn.reflected) / 20000.0 - 1.0) < 1e-12);

    /* anisotropic phase still agrees with its analytic reference */
    double g2 = 0.5;
    double a2 = analytic_total(&m, &beam, L, g2);
    Stat s2; stat_init(&s2);
    for (int i = 0; i < N; i++)
        stat_push(&s2, trace_slab(&m, &beam, L, g2, rng_u01(&rng)).total);
    printf("hg_g0.5: mc=%.6f analytic=%.6f\n", s2.mean, a2);
    check("anisotropic_matches", fabs(s2.mean - a2) < 3.0 * stat_se(&s2));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
