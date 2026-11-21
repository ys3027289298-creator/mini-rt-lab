
/* sample: stratified and low-discrepancy sampling with variance estimate,
   plus environment-map importance sampling over a 2D luminance CDF. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define N_MAX 256

static const double PI = 3.14159265358979323846;

typedef struct { double mean; double var; int n; } Stat;

static double halton2(unsigned i, unsigned base) {
    double f = 1.0, r = 0.0; unsigned k = i;
    while (k > 0) { f /= base; r += f * (k % base); k /= base; }
    return r;
}

/* estimate integral of f over [0,1) with stratified samples */
static Stat stratify(int strata, double (*f)(double)) {
    Stat s; memset(&s, 0, sizeof(s));
    double sum = 0.0, sum2 = 0.0;
    int n = 0;
    for (int i = 0; i < strata; i++) {
        double u = ((double)i + 0.5) / (double)strata;
        double v = f(u);
        sum += v; sum2 += v * v; n++;
    }
    s.n = n;
    s.mean = sum / n;
    double var = (n > 1) ? (sum2 / n - s.mean * s.mean) * n / (n - 1) : 0.0;
    s.var = var > 0 ? var : 0;
    return s;
}

static double f_quad(double u) { return u * u; }
static double f_const(double u) { (void)u; return 1.0; }

/* ---- U[0,1) stream for Monte Carlo runs (xorshift32) ---- */
typedef struct { unsigned s; } Rng;
static double rng_next(Rng* r) {
    unsigned x = r->s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    r->s = x;
    return (double)x / 4294967296.0;
}

/* ---- environment map importance sampling ----
   Radiance is parameterized by theta in [0,pi] (cosine of the polar angle)
   and phi in [0,2pi). Row j spans theta in [j*pi/H,(j+1)*pi/H] and column i
   spans phi in [i*2pi/W,(i+1)*2pi/W). The table is a 2D CDF of the
   luminance-weighted solid angle of each texel:
     w[j][i]  = luminance(midpoint radiance) * solid angle of texel
     row[j+1] = row[j] + sum_i w[j][i]          (theta CDF)
     col[j][i+1] = col[j][i] + w[j][i]          (phi CDF within row j)
   Sampling inverts row-major, then maps the fractional offsets uniformly in
   cos(theta) and phi inside the chosen texel. pdf is therefore the texel's
   luminance times its solid angle over the total, divided by the texel's
   solid angle: piecewise constant over texels, 0 where the table has no
   luminance. */
#define ENV_H 32
#define ENV_W 64
#define ENV_K 40.0

static const double ENV_LOBE[3] = { 12.0, 10.0, 8.0 };   /* sun radiance peak */
static const double ENV_BASE[3] = { 0.02, 0.03, 0.06 };  /* sky base radiance */

typedef struct {
    double w[ENV_H][ENV_W];
    double row[ENV_H + 1];
    double col[ENV_H][ENV_W + 1];
    double total;
} EnvTable;

/* luminance of the radiance field: a diffuse base plus a sharp sun lobe
   exp(-K*(1-cos)) centered on SUN_DIR = (sin pi/3 cos pi/4, sin pi/3 sin pi/4,
   cos pi/3). Its integral over the sphere is closed form (see
   env_luminance_analytic), so the table can be checked against it exactly. */
static const double SUN_DIR[3] = {
    0.6123724356957945, 0.6123724356957945, 0.5
};

static double env_luminance_texel(int j, int i) {
    double theta = ((double)j + 0.5) * PI / ENV_H;
    double phi = ((double)i + 0.5) * 2.0 * PI / ENV_W;
    double st = sin(theta);
    double d[3] = { st * cos(phi), st * sin(phi), cos(theta) };
    double cs = d[0] * SUN_DIR[0] + d[1] * SUN_DIR[1] + d[2] * SUN_DIR[2];
    double g = exp(-ENV_K * (1.0 - cs));
    double lum = 0.0;
    for (int k = 0; k < 3; k++) {
        double c = ENV_BASE[k] + ENV_LOBE[k] * g;
        lum += c * (k == 0 ? 0.2126 : (k == 1 ? 0.7152 : 0.0722));
    }
    return lum;
}

static double env_texel_solid_angle(int j) {
    double t0 = (double)j * PI / ENV_H, t1 = (double)(j + 1) * PI / ENV_H;
    return (2.0 * PI / ENV_W) * (cos(t0) - cos(t1));
}

static void env_build(EnvTable* e) {
    memset(e, 0, sizeof(*e));
    e->row[0] = 0.0;
    for (int j = 0; j < ENV_H; j++) {
        double om = env_texel_solid_angle(j);
        e->col[j][0] = 0.0;
        for (int i = 0; i < ENV_W; i++) {
            e->w[j][i] = env_luminance_texel(j, i) * om;
            e->col[j][i + 1] = e->col[j][i] + e->w[j][i];
        }
        e->row[j + 1] = e->row[j] + e->col[j][ENV_W];
    }
    e->total = e->row[ENV_H];
}

/* pdf of a texel under the table; 0 where it carries no luminance */
static double env_pdf_texel(const EnvTable* e, int j, int i) {
    if (e->w[j][i] <= 0.0) return 0.0;
    return e->w[j][i] / e->total / env_texel_solid_angle(j);
}

/* pdf of a direction under the table, or 0 outside the support */
static double env_pdf(const EnvTable* e, const double d[3]) {
    if (e->total <= 0.0) return 1.0 / (4.0 * PI); /* mirrored: uniform */
    double ct = d[2] < -1.0 ? -1.0 : (d[2] > 1.0 ? 1.0 : d[2]);
    double theta = acos(ct);
    double phi = atan2(d[1], d[0]);
    if (phi < 0.0) phi += 2.0 * PI;
    int j = (int)(theta / (PI / ENV_H));
    int i = (int)(phi / (2.0 * PI / ENV_W));
    if (j > ENV_H - 1) j = ENV_H - 1;
    if (i > ENV_W - 1) i = ENV_W - 1;
    return env_pdf_texel(e, j, i);
}

/* draw a direction from the table; returns its pdf through *pdf_out.
   When the table has no luminance the distribution degenerates to uniform
   over the sphere and the pdf is 1/(4pi). */
static void env_sample(const EnvTable* e, double u1, double u2, double d[3], double* pdf_out) {
    if (e->total <= 0.0) {
        double ct = 1.0 - 2.0 * u1;
        double st = sqrt(1.0 - ct * ct > 0.0 ? 1.0 - ct * ct : 0.0);
        double phi = 2.0 * PI * u2;
        d[0] = st * cos(phi); d[1] = st * sin(phi); d[2] = ct;
        *pdf_out = 1.0 / (4.0 * PI);
        return;
    }
    double t = u1 * e->total;
    int j = 0;
    while (j < ENV_H - 1 && e->row[j + 1] <= t) j++;
    double rw = e->row[j + 1] - e->row[j];
    double f = rw > 0.0 ? (t - e->row[j]) / rw : 0.0;
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;

    double t0 = (double)j * PI / ENV_H, t1 = (double)(j + 1) * PI / ENV_H;
    double c0 = cos(t0), c1 = cos(t1);
    double spine = c0 + f * (c1 - c0);   /* sampled cos(theta) */
    double st = sqrt(1.0 - spine * spine > 0.0 ? 1.0 - spine * spine : 0.0);

    double s = u2 * rw;
    int i = 0;
    while (i < ENV_W - 1 && e->col[j][i + 1] <= s) i++;
    double cw = e->col[j][i + 1] - e->col[j][i];
    double fc = cw > 0.0 ? (s - e->col[j][i]) / cw : 0.0;
    if (fc < 0.0) fc = 0.0;
    if (fc > 1.0) fc = 1.0;
    double phi = ((double)i + fc) * 2.0 * PI / ENV_W;

    d[0] = st * cos(phi); d[1] = st * sin(phi); d[2] = spine;
    *pdf_out = e->w[j][i] / e->total / env_texel_solid_angle(j);
}

/* closed-form luminance integral over the sphere for the above field:
   base*4pi + sum(lobe chrominance)*2pi*(1-exp(-2K))/K */
static double env_luminance_analytic(void) {
    double lum_base = 0.2126 * ENV_BASE[0] + 0.7152 * ENV_BASE[1] + 0.0722 * ENV_BASE[2];
    double lum_lobe = 0.2126 * ENV_LOBE[0] + 0.7152 * ENV_LOBE[1] + 0.0722 * ENV_LOBE[2];
    return lum_base * 4.0 * PI + lum_lobe * 2.0 * PI * (1.0 - exp(-2.0 * ENV_K)) / ENV_K;
}

/* luminance at a direction, same field as the table */
static double env_luminance_dir(const double d[3]) {    double cs = d[0] * SUN_DIR[0] + d[1] * SUN_DIR[1] + d[2] * SUN_DIR[2];
    double g = exp(-ENV_K * (1.0 - cs));
    return 0.2126 * (ENV_BASE[0] + ENV_LOBE[0] * g) + 0.7152 * (ENV_BASE[1] + ENV_LOBE[1] * g) +
           0.0722 * (ENV_BASE[2] + ENV_LOBE[2] * g);
}

/* narrow GGX-like reflection lobe around the sun: (n.w)^power, w = sun dir */
#define SPEC_POWER 64.0
static double env_spec_lobe(const double d[3]) {
    double cs = d[0] * SUN_DIR[0] + d[1] * SUN_DIR[1] + d[2] * SUN_DIR[2];
    return cs > 0.0 ? pow(cs, SPEC_POWER) : 0.0;
}

typedef struct { double mean; double uniform; double importance; } VarSet;

/* Reference mean and variances of the integrand L*dot over the sphere, with
   the same piecewise-constant pdf the table produces. Texels are subdivided
   QUAD_SUB x QUAD_SUB to keep the quadrature error well below the ratio. */
#define QUAD_SUB 4
static VarSet env_reference_variance(const EnvTable* e) {
    VarSet v;
    double sum = 0.0, sum_u2 = 0.0, sum_i2 = 0.0;
    for (int j = 0; j < ENV_H; j++) {
        double t0 = (double)j * PI / ENV_H, t1 = (double)(j + 1) * PI / ENV_H;
        double dct = (cos(t0) - cos(t1)) / QUAD_SUB;
        for (int i = 0; i < ENV_W; i++) {
            double p0 = (double)i * 2.0 * PI / ENV_W;
            double p1 = (double)(i + 1) * 2.0 * PI / ENV_W;
            double dph = (p1 - p0) / QUAD_SUB;
            double pdf = e->w[j][i] > 0.0 ? e->w[j][i] / e->total / env_texel_solid_angle(j) : 0.0;
            for (int a = 0; a < QUAD_SUB; a++) {
                double c = cos(t0) - ((double)a + 0.5) * dct;
                double th = acos(c), st = sin(th);
                for (int b = 0; b < QUAD_SUB; b++) {
                    double ph = p0 + ((double)b + 0.5) * dph;
                    double d[3] = { st * cos(ph), st * sin(ph), c };
                    double g = env_luminance_dir(d) * env_spec_lobe(d);
                    double om = dct * dph;
                    sum += g * om;
                    sum_u2 += (4.0 * PI * g) * (4.0 * PI * g) * om / (4.0 * PI);
                    if (pdf > 0.0) sum_i2 += g * g / pdf * om;
                }
            }
        }
    }
    v.mean = sum;
    v.uniform = sum_u2 - sum * sum;
    v.importance = sum_i2 - sum * sum;
    return v;
}

typedef struct {
    double mean_u, mean_i;
    double var_u, var_i;
    int n;
} EnvStats;

/* Monte Carlo estimate of the same integral under both samplers */
static EnvStats env_estimate(const EnvTable* e, int n, unsigned seed) {
    EnvStats s;
    memset(&s, 0, sizeof(s));
    Rng rng = { seed };
    double su = 0.0, su2 = 0.0, si = 0.0, si2 = 0.0;
    for (int k = 0; k < n; k++) {
        double u1 = rng_next(&rng), u2 = rng_next(&rng);

        /* uniform over the sphere: solid angle 4pi */
        double ct = 1.0 - 2.0 * u1;
        double st = sqrt(1.0 - ct * ct > 0.0 ? 1.0 - ct * ct : 0.0);
        double ph = 2.0 * PI * u2;
        double du[3] = { st * cos(ph), st * sin(ph), ct };
        double eu = 4.0 * PI * env_luminance_dir(du) * env_spec_lobe(du);
        su += eu; su2 += eu * eu;

        /* environment importance sampling */
        double di[3], pdf;
        env_sample(e, u1, u2, di, &pdf);
        double ei = pdf > 0.0 ? env_luminance_dir(di) * env_spec_lobe(di) / pdf : 0.0;
        si += ei; si2 += ei * ei;
    }
    s.n = n;
    s.mean_u = su / n; s.var_u = (su2 / n - s.mean_u * s.mean_u) * n / (n - 1);
    s.mean_i = si / n; s.var_i = (si2 / n - s.mean_i * s.mean_i) * n / (n - 1);
    return s;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    check("halton_in_unit", halton2(1, 2) >= 0.0 && halton2(1, 2) < 1.0);
    check("halton_first", fabs(halton2(0, 2)) < 1e-12);

    Stat a = stratify(64, f_quad);
    check("quad_mean_close", fabs(a.mean - 1.0 / 3.0) < 1e-3);

    Stat b = stratify(16, f_quad);
    Stat c = stratify(64, f_quad);
    check("variance_drops", c.var <= b.var + 1e-6);

    Stat d = stratify(32, f_const);
    check("const_exact", fabs(d.mean - 1.0) < 1e-12 && d.var < 1e-12);

    check("halton_distinct", halton2(1, 2) != halton2(2, 2));
    check("halton_monotone_base", halton2(1, 2) > halton2(0, 2));

    /* ---- environment importance sampling ---- */
    EnvTable env;
    env_build(&env);

    /* table mass vs analytic luminance integral */
    double analytic = env_luminance_analytic();
    printf("env_table_mass=%.9f env_analytic=%.9f env_rel_err=%.2e\n",
           env.total, analytic, fabs(env.total - analytic) / analytic);
    check("env_table_matches_analytic", fabs(env.total - analytic) / analytic < 1e-3);

    /* pdf times solid angle must sum to 1 over the sphere */
    double mass = 0.0;
    for (int j = 0; j < ENV_H; j++)
        for (int i = 0; i < ENV_W; i++)
            if (env.w[j][i] > 0.0) mass += env_pdf_texel(&env, j, i) * env_texel_solid_angle(j);
    check("env_pdf_normalized", fabs(mass - 1.0) < 1e-9);

    /* the direction returned by env_sample must be the one its own pdf
       describes: re-deriving the texel from (theta,phi) must give the same pdf */
    int inside = 0;
    Rng rng_draw = { 12345u };   /* independent stream per test */
    for (int k = 0; k < 4096; k++) {
        double d[3], pdf;
        env_sample(&env, rng_next(&rng_draw), rng_next(&rng_draw), d, &pdf);
        double pa = env_pdf(&env, d);
        if (fabs(pdf - pa) <= 1e-12 * (pa > 0.0 ? pa : 1.0)) inside++;
    }
    check("env_sample_pdf_consistent", inside == 4096);

    /* chi-square: grouped sample counts follow the table's luminance mass */
#define CHI_J 8
#define CHI_I 16
    int counts[CHI_J * CHI_I] = { 0 };
    double expect[CHI_J * CHI_I] = { 0.0 };
    int nsamp = 1 << 19;
    for (int j = 0; j < ENV_H; j++)
        for (int i = 0; i < ENV_W; i++)
            expect[(j * CHI_J / ENV_H) * CHI_I + (i * CHI_I / ENV_W)] += nsamp * env.w[j][i] / env.total;
    Rng rng_chi = { 777u };
    for (int k = 0; k < nsamp; k++) {
        double d[3], pdf;
        env_sample(&env, rng_next(&rng_chi), rng_next(&rng_chi), d, &pdf);
        double th = acos(d[2] < -1.0 ? -1.0 : (d[2] > 1.0 ? 1.0 : d[2]));
        double ph = atan2(d[1], d[0]);
        if (ph < 0.0) ph += 2.0 * PI;
        int gj = (int)(th / (PI / CHI_J)); if (gj > CHI_J - 1) gj = CHI_J - 1;
        int gi = (int)(ph / (2.0 * PI / CHI_I)); if (gi > CHI_I - 1) gi = CHI_I - 1;
        counts[gj * CHI_I + gi]++;
    }
    double chi2 = 0.0;
    for (int k = 0; k < CHI_J * CHI_I; k++)
        chi2 += (counts[k] - expect[k]) * (counts[k] - expect[k]) / expect[k];
    /* 128 bins whose expected counts follow the luminance; the smoothed
       table leaves ~100 for a correct sampler, so 200 is a loose reject */
    check("env_chi2", chi2 < 200.0);

    /* variance versus uniform sampling; ratio against the reference too */
    VarSet ref = env_reference_variance(&env);
    EnvStats st = env_estimate(&env, 1 << 19, 2463534242u);
    double ref_ratio = ref.uniform / ref.importance;
    double mc_ratio = st.var_u / st.var_i;
    printf("env_ref_mean=%.6f env_ref_ratio=%.2f env_mc_mean_u=%.6f env_mc_mean_i=%.6f\n",
           ref.mean, ref_ratio, st.mean_u, st.mean_i);
    printf("env_var_u=%.4f env_var_i=%.4f env_var_ratio=%.2f (sigma_u=%.4f sigma_i=%.4f)\n",
           st.var_u, st.var_i, mc_ratio, sqrt(st.var_u), sqrt(st.var_i));
    check("env_importance_matches_ref", fabs(mc_ratio - ref_ratio) / ref_ratio < 0.05);
    check("env_importance_better", mc_ratio > 10.0);
    check("env_mc_means_agree", fabs(st.mean_u / st.mean_i - 1.0) < 0.05);

    /* all-zero luminance: both samplers fall back to uniform, exactly */
    EnvTable dark;
    memset(&dark, 0, sizeof(dark));
    Rng rng_dark = { 987654321u };
    int dark_ok = 1;
    for (int k = 0; k < 256; k++) {
        double d[3], pdf;
        env_sample(&dark, rng_next(&rng_dark), rng_next(&rng_dark), d, &pdf);
        if (fabs(pdf - 1.0 / (4.0 * PI)) > 1e-12) dark_ok = 0;
        if (fabs(env_pdf(&dark, d) - 1.0 / (4.0 * PI)) > 1e-12) dark_ok = 0;
        if (fabs(sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) - 1.0) > 1e-9) dark_ok = 0;
    }
    check("env_dark_fallback", dark_ok);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
