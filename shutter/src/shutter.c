
/* shutter: shutter-interval time sampling for motion blur.
   A hard-edge box object moves along a constant motion vector; exposure is a
   time average over the open shutter [t0, t0+dt). Samples are stratified in
   time and weighted by the motion-vector segment each one covers. Estimates
   are checked against the analytic box-convolution integral, with a
   deterministic error upper bound and convergence as the sample count grows. */
#include <stdio.h>
#include <math.h>

#define N_MAX 256

typedef struct { float t; float p; float w; } TimeSample;

/* hard-edge spatial footprint: 1 within half-width r, 0 outside */
static float box_kernel(float d, float r) { return fabsf(d) <= r ? 1.0f : 0.0f; }

/* stratified midpoint samples on [t0, t0+dt). Each sample represents a
   motion-vector segment of length |v|*dt/N; its weight is that segment over
   the whole path. dt <= 0 (zero shutter) or n <= 0 is rejected. */
static int shutter_samples(float t0, float dt, float p0, float v, int n, TimeSample* out) {
    if (dt <= 0.0f || n <= 0 || n > N_MAX) return -1;
    float h = dt / (float)n;
    float path = fabsf(v) * dt;
    for (int i = 0; i < n; i++) {
        float tc = t0 + ((float)i + 0.5f) * h;
        out[i].t = tc;
        out[i].p = p0 + v * (tc - t0);
        out[i].w = (path > 0.0f) ? (fabsf(v) * h / path) : (1.0f / (float)n);
    }
    return 0;
}

/* jittered variant: same cell weights, one random point inside each stratum */
static unsigned lcg_state;
static float lcg_unit(void) {
    lcg_state = lcg_state * 1103515245u + 12345u;
    return (float)((lcg_state >> 8) & 0xffffffu) / (float)0x1000000;
}
static int shutter_samples_jitter(float t0, float dt, float p0, float v, int n,
                                  unsigned seed, TimeSample* out) {
    if (dt <= 0.0f || n <= 0 || n > N_MAX) return -1;
    lcg_state = seed;
    float h = dt / (float)n;
    float path = fabsf(v) * dt;
    for (int i = 0; i < n; i++) {
        float tc = t0 + ((float)i + lcg_unit()) * h;
        out[i].t = tc;
        out[i].p = p0 + v * (tc - t0);
        out[i].w = (path > 0.0f) ? (fabsf(v) * h / path) : (1.0f / (float)n);
    }
    return 0;
}

/* weighted exposure estimate at pixel x */
static float sample_exposure(const TimeSample* s, int n, float x, float r) {
    double acc = 0.0;
    for (int i = 0; i < n; i++)
        acc += (double)s[i].w * (double)box_kernel(x - s[i].p, r);
    return (float)acc;
}

/* analytic exposure: (1/dt) integral of k(x - p(t)) over the shutter,
   p(t) = p0 + v*(t-t0). With s = v*(t-t0) time becomes path length, so the
   integral is the overlap of [u-r, u+r] with the swept interval, over |disp|. */
static float analytic_exposure(float x, float p0, float r, float disp) {
    float L = fabsf(disp);
    if (L <= 0.0f) return box_kernel(x - p0, r);
    float u = x - p0;
    float a = fminf(0.0f, disp), b = fmaxf(0.0f, disp);
    float ov = fminf(u + r, b) - fmaxf(u - r, a);
    return ov > 0.0f ? ov / L : 0.0f;
}

/* worst |estimate - analytic| over a dense pixel sweep */
static double worst_error(int n, float p0, float v, float r, float dt) {
    TimeSample s[N_MAX];
    if (shutter_samples(0.0f, dt, p0, v, n, s) != 0) return -1.0;
    float disp = v * dt;
    float lo = p0 + fminf(0.0f, disp) - r - 1.0f;
    float hi = p0 + fmaxf(0.0f, disp) + r + 1.0f;
    double worst = 0.0;
    int m = 4001;
    for (int i = 0; i < m; i++) {
        float x = lo + (hi - lo) * (float)i / (float)(m - 1);
        double e = fabs((double)sample_exposure(s, n, x, r) -
                        (double)analytic_exposure(x, p0, r, disp));
        if (e > worst) worst = e;
    }
    return worst;
}

/* RMSE of jittered estimates over seeds and the ramp-zone pixels */
static double jitter_rmse(int n, int seeds, float p0, float v, float r, float dt) {
    const float offsets[4] = { 0.0f, 0.5f, 0.95f, 1.0f };
    double acc = 0.0;
    int cnt = 0;
    for (int sd = 1; sd <= seeds; sd++) {
        TimeSample s[N_MAX];
        shutter_samples_jitter(0.0f, dt, p0, v, n, (unsigned)sd * 7u, s);
        for (int k = 0; k < 4; k++) {
            float x = p0 + v * dt * offsets[k];
            double e = (double)sample_exposure(s, n, x, r) -
                       (double)analytic_exposure(x, p0, r, v * dt);
            acc += e * e;
            cnt++;
        }
    }
    return sqrt(acc / (double)cnt);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    TimeSample s[N_MAX];

    /* zero / negative shutter and bad counts are rejected */
    check("reject_zero_shutter", shutter_samples(0.0f, 0.0f, 0.0f, 1.0f, 8, s) == -1);
    check("reject_neg_shutter", shutter_samples(0.0f, -1.0f, 0.0f, 1.0f, 8, s) == -1);
    check("reject_bad_count", shutter_samples(0.0f, 1.0f, 0.0f, 1.0f, 0, s) == -1);
    check("jitter_reject_zero",
          shutter_samples_jitter(0.0f, 0.0f, 0.0f, 1.0f, 8, 1u, s) == -1);

    /* sample placement and motion-vector weights */
    int n = 16;
    float t0 = 1.0f, dt = 0.5f, p0 = 2.0f, v = 3.0f;
    check("sample_build_ok", shutter_samples(t0, dt, p0, v, n, s) == 0);
    double wsum = 0.0;
    int positive = 1, inside = 1, onvector = 1;
    for (int i = 0; i < n; i++) {
        wsum += s[i].w;
        if (!(s[i].w > 0.0f)) positive = 0;
        if (!(s[i].t > t0 && s[i].t < t0 + dt)) inside = 0;
        float expect_p = p0 + v * (s[i].t - t0);
        if (fabsf(s[i].p - expect_p) > 1e-5f) onvector = 0;
    }
    check("weights_unit_sum", fabs(wsum - 1.0) < 1e-6);
    check("weights_positive", positive);
    check("times_inside_shutter", inside);
    check("positions_on_vector", onvector);

    /* stationary object: zero-length motion vector, equal weights, no blur */
    shutter_samples(0.0f, 1.0f, 5.0f, 0.0f, n, s);
    check("stationary_on", fabsf(sample_exposure(s, n, 5.0f, 0.5f) - 1.0f) < 1e-6f);
    check("stationary_off", sample_exposure(s, n, 9.0f, 0.5f) == 0.0f);
    check("stationary_analytic",
          analytic_exposure(5.0f, 5.0f, 0.5f, 0.0f) == 1.0f &&
          analytic_exposure(9.0f, 5.0f, 0.5f, 0.0f) == 0.0f);

    /* analytic trapezoid: r=1, swept length L=10, plateau 2r/L = 0.2 */
    float r = 1.0f, L = 10.0f;
    check("analytic_zero_outside",
          analytic_exposure(-3.0f, 0.0f, r, L) == 0.0f &&
          analytic_exposure(12.0f, 0.0f, r, L) == 0.0f);
    check("analytic_plateau", fabsf(analytic_exposure(5.0f, 0.0f, r, L) - 0.2f) < 1e-6f);
    check("analytic_ramp", fabsf(analytic_exposure(0.0f, 0.0f, r, L) - 0.1f) < 1e-6f);
    check("analytic_symmetric",
          fabsf(analytic_exposure(3.0f, 0.0f, r, L) -
                analytic_exposure(7.0f, 0.0f, r, L)) < 1e-6f);

    /* energy conservation: integral over x equals the object width 2r */
    double energy = 0.0;
    float step = 1e-3f;
    for (float x = -3.0f; x <= 13.0f; x += step)
        energy += (double)analytic_exposure(x, 0.0f, r, L) * step;
    check("analytic_energy", fabs(energy - 2.0 * r) < 1e-3);

    /* midpoint rule: f(t) in {0,1} has at most two discontinuities, so at
       most two strata err, each by <= its weight 1/N -> |err| <= 2/N */
    double w4 = worst_error(4, 0.0f, 10.0f, r, 1.0f);
    double w16 = worst_error(16, 0.0f, 10.0f, r, 1.0f);
    double w64 = worst_error(64, 0.0f, 10.0f, r, 1.0f);
    check("bound_n4", w4 <= 2.0 / 4 + 1e-6);
    check("bound_n16", w16 <= 2.0 / 16 + 1e-6);
    check("bound_n64", w64 <= 2.0 / 64 + 1e-6);
    check("converges", w64 < w16 && w16 < w4);

    /* the bound is direction-agnostic: negative motion vector */
    check("negative_velocity_bound", worst_error(16, 0.0f, -10.0f, r, 1.0f) <= 2.0 / 16 + 1e-6);

    /* a single time sample reproduces the hard edge: pure 0/1, no soft ramp */
    TimeSample one[1];
    shutter_samples(0.0f, 1.0f, 0.0f, 10.0f, 1, one);
    float hard = sample_exposure(one, 1, 0.0f, r);
    check("single_sample_hard_edge",
          (hard == 0.0f || hard == 1.0f) &&
          fabsf(analytic_exposure(0.0f, 0.0f, r, L) - 0.1f) < 1e-6f);

    /* hand-built unequal cells: weights must follow segment lengths, and the
       same 2-boundary argument gives |err| <= 2 * max(cell weight) */
    float cells[3][2] = { { 0.0f, 0.5f }, { 0.5f, 0.8f }, { 0.8f, 1.0f } };
    TimeSample nu[3];
    double nu_wsum = 0.0;
    float maxw = 0.0f;
    for (int i = 0; i < 3; i++) {
        float tc = 0.5f * (cells[i][0] + cells[i][1]);
        nu[i].t = tc;
        nu[i].p = 10.0f * tc;
        nu[i].w = (cells[i][1] - cells[i][0]); /* |v|*h / (|v|*dt) */
        nu_wsum += nu[i].w;
        if (nu[i].w > maxw) maxw = nu[i].w;
    }
    double nu_worst = 0.0;
    for (int i = 0; i < 2001; i++) {
        float x = -2.0f + 14.0f * (float)i / 2000.0f;
        double e = fabs((double)sample_exposure(nu, 3, x, r) -
                        (double)analytic_exposure(x, 0.0f, r, L));
        if (e > nu_worst) nu_worst = e;
    }
    check("nonuniform_weights", fabs(nu_wsum - 1.0) < 1e-6);
    check("nonuniform_bound", nu_worst <= 2.0 * maxw + 1e-6);

    /* jittered sampling: empirical RMSE drops as N grows */
    double r8 = jitter_rmse(8, 64, 0.0f, 10.0f, r, 1.0f);
    double r256 = jitter_rmse(256, 64, 0.0f, 10.0f, r, 1.0f);
    check("jitter_converges", r256 < r8);
    check("jitter_small", r256 < 0.08);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
