
/* depth: adaptive recursion depth for path tracing.
   Per-bounce continuation is decided from path throughput and surface
   roughness. Three modes are compared against an analytic ground truth:
   fixed depth, adaptive biased (threshold cut with an analytic tail
   bound), adaptive unbiased (Russian roulette). Zero-throughput paths
   terminate early and are accounted. */
#include <stdio.h>
#include <math.h>
#include <string.h>

enum { MODE_FIXED, MODE_BIASED, MODE_UNBIASED };

typedef struct { unsigned s; } Rng;
static float rng_next(Rng* r) {
    r->s = r->s * 1664525u + 1013904223u;
    return (float)(r->s >> 8) / 16777216.0f;
}

/* Bounce chain: each bounce emits `emission` and multiplies throughput
   by a two-point albedo in {0, albedo_hi} with P(hi) = albedo_q.
   Expected radiance is analytic: emission / (1 - q*hi). */
typedef struct {
    float albedo_hi;  /* nonzero albedo support, also the tail-bound multiplier */
    float albedo_q;   /* P(albedo = albedo_hi); zero draw => zero throughput */
    float roughness;  /* per-bounce roughness in [0,1] */
    float emission;   /* per-bounce emitted radiance */
} Scene;

typedef struct {
    int mode;
    int max_depth;     /* hard cap for every mode */
    float threshold;   /* continuation score cutoff (biased cut / RR trigger) */
    float rr_min;      /* unbiased mode: minimum continuation probability */
    float rough_gamma; /* roughness exponent in the continuation score */
} DepthCfg;

typedef struct {
    double radiance;
    int bounces;         /* shading samples used */
    double bias_bound;   /* biased mode: worst-case dropped tail */
    int zero_throughput; /* terminated on exactly zero throughput */
    int rr_killed;       /* unbiased mode: terminated by Russian roulette */
} Path;

static double scene_truth(const Scene* sc) {
    return sc->emission / (1.0 - (double)sc->albedo_q * sc->albedo_hi);
}

static Path trace(const Scene* sc, const DepthCfg* cfg, Rng* rng) {
    Path p; memset(&p, 0, sizeof(p));
    double T = 1.0;
    double keep = pow(1.0 - sc->roughness, cfg->rough_gamma);
    for (int d = 0; d < cfg->max_depth; d++) {
        p.radiance += sc->emission * T;
        p.bounces++;
        double Tn = T * (rng_next(rng) < sc->albedo_q ? sc->albedo_hi : 0.0f);
        if (Tn <= 0.0) { p.zero_throughput = 1; break; } /* terminate + account */
        double score = Tn * keep;
        if (cfg->mode != MODE_FIXED && score < cfg->threshold) {
            if (cfg->mode == MODE_BIASED) {
                /* dropped tail <= Tn * e * (1 + hi + hi^2 + ...) */
                p.bias_bound = Tn * sc->emission / (1.0 - sc->albedo_hi);
                break;
            }
            double pr = score / cfg->threshold; /* in (0,1) */
            if (pr < cfg->rr_min) pr = cfg->rr_min;
            if (rng_next(rng) < pr) Tn /= pr;   /* reweight: stays unbiased */
            else { p.rr_killed = 1; break; }
        }
        T = Tn;
    }
    return p;
}

typedef struct {
    double sum, sum2, bound, bounces;
    long zero_tp, rr_killed;
    long n;
} Stats;

static void stats_add(Stats* s, const Path* p) {
    s->sum += p->radiance;
    s->sum2 += p->radiance * p->radiance;
    s->bound += p->bias_bound;
    s->bounces += p->bounces;
    s->zero_tp += p->zero_throughput;
    s->rr_killed += p->rr_killed;
    s->n++;
}

static double stats_mean(const Stats* s) { return s->sum / s->n; }
static double stats_stderr(const Stats* s) {
    double m = stats_mean(s);
    double v = s->sum2 / s->n - m * m;
    return sqrt(v > 0 ? v / s->n : 0.0);
}

static Stats run(const Scene* sc, const DepthCfg* cfg, long n, unsigned seed) {
    Stats s; memset(&s, 0, sizeof(s));
    Rng rng = { seed };
    for (long i = 0; i < n; i++) {
        Path p = trace(sc, cfg, &rng);
        stats_add(&s, &p);
    }
    return s;
}

static void report(const char* scene, const char* mode, double truth, const Stats* s) {
    printf("scene=%s mode=%s est=%.4f truth=%.4f bias=%.4f bound=%.4f bounces=%.2f zero_tp=%ld rr_kill=%ld\n",
           scene, mode, stats_mean(s), truth, truth - stats_mean(s),
           s->bound / s->n, s->bounces / s->n, s->zero_tp, s->rr_killed);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    /* smooth: deep paths matter (truth=10); rough: energy dies fast (truth=3.33) */
    Scene smooth = { 0.9f / 0.95f, 0.95f, 0.05f, 1.0f };
    Scene rough  = { 0.7f / 0.9f, 0.9f, 0.9f, 1.0f };
    check("energy_bounded", smooth.albedo_hi < 1.0f && rough.albedo_hi < 1.0f);

    /* zero throughput terminates immediately and is accounted */
    Scene dead = { 0.0f, 1.0f, 0.5f, 1.0f };
    DepthCfg probe = { MODE_BIASED, 32, 0.05f, 0.05f, 1.0f };
    Rng r0 = { 7 };
    Path z = trace(&dead, &probe, &r0);
    check("zero_throughput_terminates",
          z.bounces == 1 && z.zero_throughput == 1 && fabs(z.radiance - 1.0) < 1e-9);

    const long n = 40000;
    DepthCfg fixed8 = { MODE_FIXED, 8, 0.05f, 0.05f, 1.0f };
    DepthCfg biased = { MODE_BIASED, 96, 0.05f, 0.05f, 1.0f };
    DepthCfg unbias = { MODE_UNBIASED, 96, 0.05f, 0.05f, 1.0f };

    double ts = scene_truth(&smooth), tr = scene_truth(&rough);

    Stats sf8 = run(&smooth, &fixed8, n, 101);
    Stats sb  = run(&smooth, &biased, n, 102);
    Stats su  = run(&smooth, &unbias, n, 103);
    Stats rf8 = run(&rough, &fixed8, n, 104);
    Stats rb  = run(&rough, &biased, n, 105);
    Stats ru  = run(&rough, &unbias, n, 106);

    report("smooth", "fixed8", ts, &sf8);
    report("smooth", "biased", ts, &sb);
    report("smooth", "unbiased", ts, &su);
    report("rough", "fixed8", tr, &rf8);
    report("rough", "biased", tr, &rb);
    report("rough", "unbiased", tr, &ru);

    /* fixed depth: empirical bias matches analytic truncation e*a^D/(1-a) */
    double a_s = smooth.albedo_q * smooth.albedo_hi;
    double fixed8_bias = smooth.emission * pow(a_s, 8) / (1.0 - a_s);
    check("fixed8_bias_matches_theory",
          fabs(ts - stats_mean(&sf8) - fixed8_bias) < 4.0 * stats_stderr(&sf8) + 1e-6);

    /* biased mode: 0 <= true bias <= reported bound */
    check("biased_bound_smooth",
          ts - stats_mean(&sb) >= -4.0 * stats_stderr(&sb) &&
          ts - stats_mean(&sb) <= sb.bound / sb.n + 4.0 * stats_stderr(&sb));
    check("biased_bound_rough",
          tr - stats_mean(&rb) >= -4.0 * stats_stderr(&rb) &&
          tr - stats_mean(&rb) <= rb.bound / rb.n + 4.0 * stats_stderr(&rb));

    /* unbiased mode: consistent with ground truth within 4 stderr */
    check("unbiased_consistent_smooth", fabs(ts - stats_mean(&su)) <= 4.0 * stats_stderr(&su));
    check("unbiased_consistent_rough", fabs(tr - stats_mean(&ru)) <= 4.0 * stats_stderr(&ru));

    /* termination accounting */
    check("zero_throughput_accounted", sb.zero_tp > 0 && rb.zero_tp > 0);
    check("rr_kills_accounted", su.rr_killed > 0 && ru.rr_killed > 0);

    /* rough paths cut earlier, smooth paths keep depth */
    check("rough_cuts_earlier", rb.bounces / rb.n < sb.bounces / sb.n);
    double reduction = 1.0 - (rb.bounces / rb.n) / (rf8.bounces / rf8.n);
    printf("rough_sample_reduction=%.3f\n", reduction);
    check("rough_saves_samples", reduction > 0.4);

    /* smooth: adaptive is closer to truth than fixed depth 8 */
    check("smooth_gains_accuracy",
          fabs(ts - stats_mean(&sb)) < fabs(ts - stats_mean(&sf8)));

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
