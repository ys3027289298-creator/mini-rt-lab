/* hdr: encoding error model for a float format with E exponent bits and M
 * fraction (mantissa) bits.
 *
 * DERIVATION
 *   A normal value is  v = s * 1.f * 2^k  with  2^k <= |v| < 2^(k+1).
 *   Each binade [2^k, 2^(k+1)) is a uniform grid of step
 *       ulp(k) = 2^(k-M).
 *   Storing v rounds it (nearest-even) to that grid, so
 *       |err|      <= ulp(k)/2 = 2^(k-M-1)                 (absolute, per binade)
 *       |err|/|v|  <= 2^(k-M-1) / 2^k = 2^-(M+1)           (relative, normals)
 *   The relative bound depends only on M, not on E.  E fixes the dynamic range
 *   [2^emin, 2^emax] and therefore the absolute stair height; it does not change
 *   relative precision.  Inverting the bound:
 *       2^-(M+1) <= target   <=>   M >= -log2(target) - 1.
 *   The absolute error doubles with every binade, so |err|(k) is a staircase
 *   that rises by 2x per octave of magnitude while |err|/|v| stays flat.
 *
 * ENDS (absolute, not relative)
 *   Overflow:  |v| >= qmax + 2^(emax-M-1)  rounds to infinity;
 *              qmax = 2^emax * (2 - 2^-M) is the largest finite value.
 *   Underflow: |v| <  2^(emin-M-1)         rounds to zero; that threshold is
 *              exactly half the smallest subnormal 2^(emin-M).  At the tie
 *              itself nearest-even picks zero.
 *   In the subnormal range [2^(emin-M), 2^emin) the grid is fixed at
 *   2^(emin-M), so the relative error grows without bound (up to 1 at the
 *   smallest value); only the absolute bound 2^(emin-M-1) survives.
 *
 * E = 0
 *   With no exponent field the only representable value is zero, so every
 *   nonzero value is unencodable.  encode() reports that as status -2.
 *
 * This program prints the derivation, verifies the model reproduces float32
 * bit-for-bit, measures the staircase, prints both-end limits as concrete
 * numbers, checks the measured relative error against the bound across several
 * formats, and reports the E = 0 case.  Output is key=value lines;
 * all_ok=1 means every check passed. */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

/* ---------- format arithmetic ---------- */

static int emax_of(int E) { return (1 << (E - 1)) - 1; }   /* largest normal exponent */
static int emin_of(int E) { return 2 - (1 << (E - 1)); }   /* smallest normal exponent */

/* largest finite value 2^emax * (2 - 2^-M) */
static double qmax_of(int E, int M) {
    return ldexp((double)(((int64_t)1 << (M + 1)) - 1), emax_of(E) - M);
}
static double smallest_subnormal(int E, int M) { return ldexp(1.0, emin_of(E) - M); }
static double rel_bound(int M) { return ldexp(1.0, -(M + 1)); }

/* smallest M with 2^-(M+1) <= target */
static int mantissa_bits(double target) {
    double x = -log2(target) - 1.0;
    int n = (int)ceil(x - 1e-12);
    return n < 0 ? 0 : n;
}

/* Round finite x to the (E,M) format with nearest-even.
 * out gets the nearest representable value as a float (exact for M <= 23).
 * Returns  0 representable
 *          1 overflow  (|x| >= qmax + half an ulp at emax; result is inf)
 *         -1 underflow (|x| < half the smallest subnormal; result is zero)
 *         -2 format cannot encode a nonzero value (E < 1) */
static int encode(double x, int E, int M, float *out) {
    if (x == 0.0) { *out = signbit(x) ? -0.0f : 0.0f; return 0; }
    if (E < 1 || M < 0) return -2;      /* no exponent field: only zero fits */
    double s = signbit(x) ? -1.0 : 1.0;
    double ax = fabs(x);
    int emax = emax_of(E), emin = emin_of(E);

    if (!(ax < qmax_of(E, M) + ldexp(1.0, emax - M - 1))) return 1;
    if (ax < ldexp(1.0, emin - M - 1)) { *out = s * 0.0f; return -1; }

    int e;
    (void)frexp(ax, &e);              /* ax = f * 2^e, f in [0.5,1) */
    int eb = e - 1;                   /* exponent of the leading bit */
    int es = eb < emin ? emin : eb;   /* fixed grid below the normal range */
    double ulp = ldexp(1.0, es - M);
    double q = nearbyint(ax / ulp);   /* nearest-even; carry into next binade is exact */
    double v = q * ulp;
    if (v == 0.0) { *out = s * 0.0f; return -1; }
    *out = (float)(s * v);
    return 0;
}

/* ---------- harness ---------- */

static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static float bits_to_float(uint32_t b) { float f; memcpy(&f, &b, sizeof f); return f; }

/* (1) the model must reproduce float32 exactly */
static void roundtrip_sweep(void) {
    uint64_t mism = 0, n = 0;
    double worst = 0.0;
    for (uint64_t b = 0; b < (1ull << 32); b += 4093) {
        uint32_t bits = (uint32_t)b;
        if ((bits & 0x7F800000u) == 0x7F800000u) continue;   /* skip inf / nan */
        float f = bits_to_float(bits), g;
        int st = encode((double)f, 8, 23, &g);
        if (st != 0 || memcmp(&f, &g, sizeof f) != 0) mism++;
        if (f != 0.0f) {
            double rel = fabs((double)g - (double)f) / fabs((double)f);
            if (rel > worst) worst = rel;
        }
        n++;
    }
    /* binade transition points: 2^k, one ulp above, one ulp below 2^(k+1) */
    for (int k = -125; k <= 124; k++) {
        double lo = ldexp(1.0, k), half = ldexp(0.5, k - 23);
        double pts[3] = { lo, lo + half, ldexp(2.0, k + 1) - half };
        for (int j = 0; j < 3; j++) {
            float f = (float)pts[j], g;
            int st = encode(pts[j], 8, 23, &g);
            if (st != 0 || f != g) mism++;
            n++;
        }
    }
    printf("roundtrip_samples=%llu roundtrip_mismatches=%llu roundtrip_worst_rel=%g\n",
           (unsigned long long)n, (unsigned long long)mism, worst);
    check("roundtrip_exact", mism == 0);
}

/* (2) staircase: per binade, worst abs error and worst rel error */
static void staircase(int E, int M) {
    double bound = rel_bound(M);
    double prev_abs = 0.0, max_rel = 0.0;
    int have_prev = 0, flat_ok = 1, abs_ok = 1, step_ok = 1;
    for (int k = -3; k <= 3; k++) {
        double lo = ldexp(1.0, k), step = ldexp(1.0, k) / 4096.0;
        double mabs = 0.0, mrel = 0.0;
        for (int j = 0; j < 4096; j++) {
            double x = lo + step * j;                   /* spans [2^k, 2^(k+1)) */
            float g;
            if (encode(x, E, M, &g) != 0) continue;
            double err = fabs((double)g - x);
            if (err > mabs) mabs = err;
            if (err / x > mrel) mrel = err / x;
        }
        if (mrel > bound * (1.0 + 1e-12)) flat_ok = 0;
        if (mabs > ldexp(1.0, k - M - 1) * (1.0 + 1e-12)) abs_ok = 0;
        if (have_prev) {
            double ratio = mabs / prev_abs;
            if (ratio < 1.9 || ratio > 2.1) step_ok = 0;   /* abs doubles per octave */
        }
        printf("binade=%d max_abs=%g max_rel=%g rel_bound=%g\n", k, mabs, mrel, bound);
        if (mrel > max_rel) max_rel = mrel;
        have_prev = 1; prev_abs = mabs;
    }
    printf("staircase E=%d M=%d measured_worst_rel=%g rel_bound=%g\n",
           E, M, max_rel, bound);
    check("staircase_rel_within_bound", flat_ok);
    check("staircase_abs_within_half_ulp", abs_ok);
    check("staircase_abs_doubles_per_octave", step_ok);
}

/* (3) measured relative error across a format's whole normal range */
static void format_sweep(int E, int M, const char *name) {
    double bound = rel_bound(M);
    double worst = 0.0;
    long tested = 0, overflow = 0;
    for (int k = emin_of(E); k < emax_of(E); k++) {
        for (int j = 0; j < 257; j++) {
            double x = ldexp(1.0 + j / 257.0, k);
            float g;
            int st = encode(x, E, M, &g);
            if (st == 1) { overflow++; continue; }
            if (st != 0) continue;
            double rel = fabs((double)g - x) / x;
            if (rel > worst) worst = rel;
            tested++;
        }
    }
    printf("format=%s E=%d M=%d rel_bound=%g measured_worst_rel=%g within=%d samples=%ld overflows=%ld\n",
           name, E, M, bound, worst, worst <= bound * (1.0 + 1e-12), tested, overflow);
    char chk[64];
    snprintf(chk, sizeof chk, "%s_measured_within_bound", name);
    check(chk, worst <= bound * (1.0 + 1e-12));
}

/* (4) both ends as concrete numbers */
static void endpoints(int E, int M, const char *name) {
    double qmax = qmax_of(E, M);
    double ovf = qmax + ldexp(1.0, emax_of(E) - M - 1);     /* first value rounding to inf */
    double ss = smallest_subnormal(E, M);
    double tie = ldexp(1.0, emin_of(E) - M - 1);            /* half the smallest subnormal */
    double mn = ldexp(1.0, emin_of(E));
    float g_max, g_ovf, g_ss, g_tie, g_tie_hi;

    int st_max = encode(qmax, E, M, &g_max);
    int st_ovf = encode(ovf, E, M, &g_ovf);
    int st_ss = encode(ss, E, M, &g_ss);
    int st_tie = encode(tie, E, M, &g_tie);
    int st_tie_hi = encode(tie + ldexp(1.0, emin_of(E) - M - 40), E, M, &g_tie_hi);

    printf("endpoints %s E=%d M=%d range=[2^%d,2^%d] max_finite=%.9g (%a) overflow_at=%.9g (%a) "
           "min_normal=%.9g min_subnormal=%.9g (%a) underflow_tie=%.9g (%a)\n",
           name, E, M, emin_of(E), emax_of(E), qmax, qmax, ovf, ovf,
           mn, ss, ss, tie, tie);
    printf("endpoints_hex %s max_finite=%a overflow_at=%a min_normal=%a min_subnormal=%a underflow_tie=%a\n",
           name, qmax, ovf, mn, ss, tie);

    char chk[64];
    snprintf(chk, sizeof chk, "%s_max_finite_exact", name);
    check(chk, st_max == 0 && (double)g_max == qmax);
    snprintf(chk, sizeof chk, "%s_overflow_rounds_to_inf", name);
    check(chk, st_ovf == 1);
    snprintf(chk, sizeof chk, "%s_min_subnormal_exact", name);
    check(chk, st_ss == 0 && (double)g_ss == ss);
    snprintf(chk, sizeof chk, "%s_underflow_tie_to_zero", name);
    check(chk, st_tie == -1 && g_tie == 0.0f);
    snprintf(chk, sizeof chk, "%s_just_above_tie_rounds_up", name);
    check(chk, st_tie_hi == 0 && (double)g_tie_hi == ss);

    /* subnormal range: absolute bound holds, relative bound does not */
    float g;
    double sub_abs = 0.0, sub_rel = 0.0, sub_ulp = ldexp(1.0, emin_of(E) - M);
    for (int j = 0; j < 4096; j++) {
        double x = sub_ulp * (j + 0.5);
        if (encode(x, E, M, &g) != 0) continue;
        double err = fabs((double)g - x);
        if (err > sub_abs) sub_abs = err;
        if (err / x > sub_rel) sub_rel = err / x;
    }
    printf("subnormal %s max_abs=%.9g abs_bound=%.9g max_rel=%g rel_bound=%g\n",
           name, sub_abs, ldexp(1.0, emin_of(E) - M - 1), sub_rel, rel_bound(M));
    snprintf(chk, sizeof chk, "%s_subnormal_abs_within_bound", name);
    check(chk, sub_abs <= ldexp(1.0, emin_of(E) - M - 1) * (1.0 + 1e-12));
    snprintf(chk, sizeof chk, "%s_subnormal_rel_exceeds_normal_bound", name);
    check(chk, sub_rel > rel_bound(M));
}

/* (5) HDR radiance through the float32 model: identity, no colour drift */
static void hdr_identity(void) {
    const double vals[] = { 0.1, 1.0 / 3.0, 1e-3, 1.0, 1000.0, 65504.0, 1e6, 3.402823466e38 };
    const int nv = (int)(sizeof vals / sizeof vals[0]);
    int bad = 0;
    for (int i = 0; i < nv; i++) {
        float g, expect = (float)vals[i];
        int st = encode(vals[i], 8, 23, &g);
        double rel = st == 0 ? fabs((double)g - vals[i]) / vals[i] : -1.0;
        if (st != 0 || g != expect) bad++;
        printf("hdr_float32 in=%.9g out=%.9g status=%d rel=%g\n", vals[i], (double)g, st, rel);
    }
    check("hdr_float32_readback_identity", bad == 0);
}

int main(void) {
    printf("model rel_bound=2^-(M+1) (independent of E); abs_bound_k=2^(k-M-1)\n");
    const double targets[] = { 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 };
    for (int i = 0; i < 7; i++) {
        int M = mantissa_bits(targets[i]);
        printf("derive target=%g min_M=%d bound=%g fits=%d\n",
               targets[i], M, rel_bound(M), rel_bound(M) <= targets[i]);
    }
    check("derive_minM_for_1e-4_is_13", mantissa_bits(1e-4) == 13);
    check("derive_minM_for_1e-7_is_23", mantissa_bits(1e-7) == 23);

    roundtrip_sweep();
    staircase(8, 10);

    /* relative bound holds across formats; E changes range, not relative error */
    format_sweep(4, 10, "e4m10");
    format_sweep(5, 10, "half");
    format_sweep(8, 7, "bfloat16");
    format_sweep(4, 23, "e4m23");
    format_sweep(5, 23, "e5m23");
    format_sweep(8, 23, "float32");

    endpoints(4, 23, "e4m23");
    endpoints(5, 10, "half");
    endpoints(5, 23, "e5m23");
    endpoints(8, 23, "float32");

    /* exponent width zero: no nonzero value can be encoded */
    float z;
    int st0 = encode(1.0, 0, 23, &z);
    int st0z = encode(0.0, 0, 23, &z);
    printf("exp_bits=0 nonzero_status=%d zero_status=%d encodable_nonzero=0\n", st0, st0z);
    check("zero_exp_bits_unencodable", st0 == -2);
    check("zero_exp_bits_zero_only", st0z == 0 && z == 0.0f);

    hdr_identity();

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
