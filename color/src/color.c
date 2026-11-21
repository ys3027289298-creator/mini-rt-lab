
/* color: linear-light render pipeline — transfer functions, linear accumulation,
   output encoding, round-trip bounds, and strict parameter validation. */
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>

typedef struct { double r, g, b; } RGB;
typedef struct { RGB rgb; double weight; } Sample;

typedef enum { E_OK = 0, E_ARG = 1, E_PARAM = 2 } Status;

typedef enum { TF_LINEAR = 0, TF_GAMMA = 1, TF_SRGB = 2 } TransferKind;

/* Transfer description. Nonlinear kinds carry no defaults: every coefficient
   must be supplied and validated. A zeroed/heap-copied struct without the
   proper *_set flag is rejected as "missing parameters". */
typedef struct {
    TransferKind kind;

    /* TF_GAMMA */
    int    gamma_set;
    double gamma;

    /* TF_SRGB: piecewise slope/power curve, all five coefficients required */
    int    curve_set;
    double cut;     /* encoded-domain breakpoint T               */
    double slope;   /* linear-segment multiplier m               */
    double scale;   /* power-segment multiplier a                */
    double offset;  /* power-segment offset o                    */
    double power;   /* power-segment exponent g (2.4 for sRGB)   */
} Transfer;

static const char* status_str(Status s) {
    switch (s) {
    case E_OK:    return "ok";
    case E_ARG:   return "bad_argument";
    case E_PARAM: return "bad_parameter";
    }
    return "unknown";
}

/* ---------- validation ---------- */

static Status transfer_validate(const Transfer* tf) {
    if (!tf) return E_ARG;
    switch (tf->kind) {
    case TF_LINEAR:
        return E_OK;
    case TF_GAMMA:
        if (!tf->gamma_set) return E_PARAM;
        if (!isfinite(tf->gamma) || !(tf->gamma > 0.0)) return E_PARAM;
        return E_OK;
    case TF_SRGB:
        if (!tf->curve_set) return E_PARAM;
        if (!isfinite(tf->cut) || !isfinite(tf->slope) ||
            !isfinite(tf->scale) || !isfinite(tf->offset) ||
            !isfinite(tf->power))
            return E_PARAM;
        if (!(tf->cut > 0.0) || !(tf->slope > 0.0) ||
            !(tf->scale > 0.0) || !(tf->offset >= 0.0) ||
            !(tf->power > 0.0))
            return E_PARAM;
        return E_OK;
    }
    return E_PARAM;
}

/* ---------- transfer functions ---------- */

static double linear_cut(const Transfer* tf) { return tf->cut / tf->slope; }

static double power_encode(const Transfer* tf, double linear) {
    return tf->scale * pow(linear, 1.0 / tf->power) - tf->offset;
}

static double power_decode(const Transfer* tf, double encoded) {
    return pow((encoded + tf->offset) / tf->scale, tf->power);
}

static Status encode1(const Transfer* tf, double linear, double* out) {
    Status st = transfer_validate(tf);
    if (st != E_OK) return st;
    if (!out || !isfinite(linear) || linear < 0.0) return E_ARG;

    switch (tf->kind) {
    case TF_LINEAR:
        *out = linear;
        return E_OK;
    case TF_GAMMA:
        *out = pow(linear, 1.0 / tf->gamma);
        return E_OK;
    case TF_SRGB:
        *out = (linear <= linear_cut(tf))
             ? tf->slope * linear
             : power_encode(tf, linear);
        return E_OK;
    }
    return E_PARAM;
}

static Status decode1(const Transfer* tf, double encoded, double* out) {
    Status st = transfer_validate(tf);
    if (st != E_OK) return st;
    if (!out || !isfinite(encoded) || encoded < 0.0) return E_ARG;

    switch (tf->kind) {
    case TF_LINEAR:
        *out = encoded;
        return E_OK;
    case TF_GAMMA:
        *out = pow(encoded, tf->gamma);
        return E_OK;
    case TF_SRGB:
        *out = (encoded <= tf->cut)
             ? encoded / tf->slope
             : power_decode(tf, encoded);
        return E_OK;
    }
    return E_PARAM;
}

#define TRY(expr) do { Status st_ = (expr); if (st_ != E_OK) return st_; } while (0)

static Status decode_rgb(const Transfer* tf, RGB in, RGB* out) {
    RGB v;
    TRY(decode1(tf, in.r, &v.r));
    TRY(decode1(tf, in.g, &v.g));
    TRY(decode1(tf, in.b, &v.b));
    if (!out) return E_ARG;
    *out = v;
    return E_OK;
}

static Status encode_rgb(const Transfer* tf, RGB in, RGB* out) {
    RGB v;
    TRY(encode1(tf, in.r, &v.r));
    TRY(encode1(tf, in.g, &v.g));
    TRY(encode1(tf, in.b, &v.b));
    if (!out) return E_ARG;
    *out = v;
    return E_OK;
}

/* ---------- disciplined pipeline: decode -> accumulate linearly -> encode ---------- */

typedef struct {
    Transfer tf;
    double r, g, b; /* linear light accumulator */
    double weight;
    int count;
} Pipeline;

static Status pipeline_init(Pipeline* p, const Transfer* tf) {
    if (!p) return E_ARG;
    Status st = transfer_validate(tf);
    if (st != E_OK) return st;
    memset(p, 0, sizeof(*p));
    p->tf = *tf;
    return E_OK;
}

/* Each encoded sample is converted to linear on entry; accumulation is linear. */
static Status pipeline_add(Pipeline* p, Sample s) {
    if (!p) return E_ARG;
    if (!isfinite(s.weight) || !(s.weight >= 0.0)) return E_ARG;
    RGB lin;
    TRY(decode_rgb(&p->tf, s.rgb, &lin));
    p->r += lin.r * s.weight;
    p->g += lin.g * s.weight;
    p->b += lin.b * s.weight;
    p->weight += s.weight;
    p->count++;
    return E_OK;
}

/* One conversion back to encoded space at the end, after normalization. */
static Status pipeline_output(const Pipeline* p, RGB* out) {
    if (!p || !out) return E_ARG;
    if (!(p->weight > 0.0)) return E_ARG;
    RGB avg = { p->r / p->weight, p->g / p->weight, p->b / p->weight };
    return encode_rgb(&p->tf, avg, out);
}

/* Reference: convert each sample exactly (per pixel), average in linear light,
   encode once. The pipeline above must agree with this to rounding. */
static Status pixelwise_reference(const Transfer* tf, const Sample* s, int n, RGB* out) {
    if (!tf || !s || !out || n <= 0) return E_ARG;
    double r = 0.0, g = 0.0, b = 0.0, w = 0.0;
    for (int i = 0; i < n; i++) {
        if (!isfinite(s[i].weight) || !(s[i].weight >= 0.0)) return E_ARG;
        RGB lin;
        TRY(decode_rgb(tf, s[i].rgb, &lin));
        r += lin.r * s[i].weight;
        g += lin.g * s[i].weight;
        b += lin.b * s[i].weight;
        w += s[i].weight;
    }
    if (!(w > 0.0)) return E_ARG;
    RGB avg = { r / w, g / w, b / w };
    return encode_rgb(tf, avg, out);
}

/* ---------- round-trip error and bounds ---------- */

/* Rounding floor for a pair of pow() calls, in absolute units on [0,1]. */
#define ROUND_FLOOR (128.0 * DBL_EPSILON)

typedef struct {
    double encoded; /* max |encode(decode(e)) - e| for display codes e in [0,1] */
    double linear;  /* max |decode(encode(l)) - l| on linear values in [0,1]  */
    double delta_e; /* upper bound on CIE76 dE between a linear colour and its
                       decode(encode(.)) round trip                            */
} RoundtripBound;

static double max2(double a, double b) { return a > b ? a : b; }

/* Lab sensitivity: f'(t) <= 1/(3*delta^2) for the sRGB/D65 observer, so a
   per-channel linear-RGB perturbation of eps changes Lab by at most K*eps. */
#define LAB_XN 0.95047
#define LAB_YN 1.0
#define LAB_ZN 1.08883
#define LAB_DELTA (6.0 / 29.0)

static double lab_lipschitz(void) {
    const double fp = 1.0 / (3.0 * LAB_DELTA * LAB_DELTA);
    const double rx = 0.4124564 + 0.3575761 + 0.1804375;
    const double ry = 0.2126729 + 0.7151522 + 0.0721750;
    const double rz = 0.0193339 + 0.1191920 + 0.9503041;
    const double dL = 116.0 * fp * ry;
    const double da = 500.0 * fp * (rx / LAB_XN + ry / LAB_YN);
    const double db = 200.0 * fp * (ry / LAB_YN + rz / LAB_ZN);
    return sqrt(dL * dL + da * da + db * db);
}

/* Analytic upper bound on the round-trip error.
   For the piecewise curve, let Lb = T/m be the linear end of the slope branch
   and Pb = a*Lb^(1/g) - o the power branch evaluated there. The two branches
   meet only if Pb == T. With delta = T - Pb:
     delta > 0 : power branch dips below the seam. encode(decode(.)) is exact;
                 linear values in (Lb, D(T)) come back off by at most
                 max(delta/m, D(T) - Lb).
     delta < 0 : power branch overshoots. decode(encode(.)) is exact; encoded
                 codes in (T, Pb) are unreachable and come back off by at most
                 max(|delta|, T - m*D(T)).
   Both directions carry a pow-rounding floor. */
static RoundtripBound roundtrip_bound(const Transfer* tf) {
    /* Rounding terms: pow/eval roundoff in both directions, plus a margin for
       the literal mismatch between the two piecewise branches. */
    const double seam_floor = 64.0 * DBL_EPSILON;
    RoundtripBound b = { 0.0, 0.0, 0.0 };
    if (transfer_validate(tf) != E_OK) {
        b.encoded = b.linear = b.delta_e = INFINITY;
        return b;
    }
    if (tf->kind == TF_LINEAR) return b;

    if (tf->kind == TF_GAMMA) {
        b.encoded = ROUND_FLOOR;
        b.linear = ROUND_FLOOR;
    } else {
        double lb = linear_cut(tf);
        double pb = power_encode(tf, lb);
        double dt = power_decode(tf, tf->cut);
        double delta = tf->cut - pb;
        if (delta >= 0.0) {
            b.encoded = seam_floor;
            b.linear = max2(delta / tf->slope, dt - lb) + seam_floor;
        } else {
            /* Power branch overshoots the seam: encoded codes in (T, Pb) are
               unreachable; both the gap itself and the slope-branch image of
               D(T) bound the encoded-side mismatch. */
            double d = -delta;
            b.encoded = max2(d, fabs(tf->cut - tf->slope * dt)) + seam_floor;
            b.linear = seam_floor;
        }
    }
    b.delta_e = lab_lipschitz() * b.linear;
    return b;
}

/* ---------- measured round trip (seam-aware) ---------- */

typedef struct { double encoded; double linear; } RoundtripMeasured;

static void measure_encoded_at(const Transfer* tf, double e, double* worst) {
    double lin, back;
    if (decode1(tf, e, &lin) != E_OK) return;
    if (encode1(tf, lin, &back) != E_OK) return;
    double err = fabs(back - e);
    if (err > *worst) *worst = err;
}

static void measure_linear_at(const Transfer* tf, double l, double* worst) {
    double enc, back;
    if (encode1(tf, l, &enc) != E_OK) return;
    if (decode1(tf, enc, &back) != E_OK) return;
    double err = fabs(back - l);
    if (err > *worst) *worst = err;
}

/* Uniform sweep plus the seam points, which a uniform grid can step over:
   the branch mismatch lives in an interval only a few 1e-9 wide. */
static RoundtripMeasured measure_roundtrip(const Transfer* tf, int steps) {
    RoundtripMeasured m = { 0.0, 0.0 };
    if (transfer_validate(tf) != E_OK || steps <= 0) {
        m.encoded = m.linear = INFINITY;
        return m;
    }
    for (int i = 0; i <= steps; i++) {
        double x = (double)i / (double)steps;
        measure_encoded_at(tf, x, &m.encoded);
        measure_linear_at(tf, x, &m.linear);
    }
    if (tf->kind == TF_SRGB) {
        double lb = linear_cut(tf);
        double dt = power_decode(tf, tf->cut);
        double seam[6] = { tf->cut, nextafter(tf->cut, 0.0), nextafter(tf->cut, 1.0),
                           lb, dt, 0.5 * (lb + dt) };
        for (int i = 0; i < 6; i++) {
            measure_encoded_at(tf, seam[i], &m.encoded);
            measure_linear_at(tf, seam[i], &m.linear);
        }
    }
    return m;
}

/* ---------- colour difference ---------- */

static double lab_f(double t) {
    return (t > LAB_DELTA * LAB_DELTA * LAB_DELTA)
         ? cbrt(t)
         : t / (3.0 * LAB_DELTA * LAB_DELTA) + 4.0 / 29.0;
}

static void rgb_to_lab(RGB c, double* L, double* A, double* B) {
    double X = 0.4124564 * c.r + 0.3575761 * c.g + 0.1804375 * c.b;
    double Y = 0.2126729 * c.r + 0.7151522 * c.g + 0.0721750 * c.b;
    double Z = 0.0193339 * c.r + 0.1191920 * c.g + 0.9503041 * c.b;
    double fx = lab_f(X / LAB_XN), fy = lab_f(Y / LAB_YN), fz = lab_f(Z / LAB_ZN);
    *L = 116.0 * fy - 16.0;
    *A = 500.0 * (fx - fy);
    *B = 200.0 * (fy - fz);
}

static double delta_e76(RGB a, RGB b) {
    double l1, a1, b1, l2, a2, b2;
    rgb_to_lab(a, &l1, &a1, &b1);
    rgb_to_lab(b, &l2, &a2, &b2);
    double dl = l1 - l2, da = a1 - a2, db = b1 - b2;
    return sqrt(dl * dl + da * da + db * db);
}

/* Max measured dE between a linear grey and its round trip, seam-aware. */
static double measure_roundtrip_delta_e(const Transfer* tf, int steps) {
    double worst = 0.0;
    if (transfer_validate(tf) != E_OK || steps <= 0) return INFINITY;
    for (int i = 0; i <= steps; i++) {
        double l = (double)i / (double)steps;
        double enc, back;
        if (encode1(tf, l, &enc) != E_OK || decode1(tf, enc, &back) != E_OK) continue;
        RGB a = { l, l, l }, b = { back, back, back };
        double de = delta_e76(a, b);
        if (de > worst) worst = de;
    }
    if (tf->kind == TF_SRGB) {
        double lb = linear_cut(tf);
        double dt = power_decode(tf, tf->cut);
        double seam[5] = { lb, dt, 0.5 * (lb + dt), nextafter(lb, 1.0), nextafter(dt, 0.0) };
        for (int i = 0; i < 5; i++) {
            double l = seam[i], enc, back;
            if (encode1(tf, l, &enc) != E_OK || decode1(tf, enc, &back) != E_OK) continue;
            RGB a = { l, l, l }, b = { back, back, back };
            double de = delta_e76(a, b);
            if (de > worst) worst = de;
        }
    }
    return worst;
}

/* ---------- self-test ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    double v;

    /* --- nonlinear parameters are required, never defaulted --- */
    /* A struct with no transfer kind set is not silently treated as linear. */
    Transfer missing = { .kind = (TransferKind)99 };
    check("unknown_kind_error", decode1(&missing, 0.5, &v) == E_PARAM);

    Transfer srgb_missing = { .kind = TF_SRGB };
    check("missing_srgb_error", decode1(&srgb_missing, 0.5, &v) == E_PARAM);
    check("missing_srgb_encode_error", encode1(&srgb_missing, 0.5, &v) == E_PARAM);

    Transfer gamma_missing = { .kind = TF_GAMMA };
    check("missing_gamma_error", decode1(&gamma_missing, 0.5, &v) == E_PARAM);
    check("zero_gamma_error",
          decode1(&(Transfer){ .kind = TF_GAMMA, .gamma_set = 1, .gamma = 0.0 }, 0.5, &v) == E_PARAM);
    check("negative_gamma_error",
          decode1(&(Transfer){ .kind = TF_GAMMA, .gamma_set = 1, .gamma = -2.2 }, 0.5, &v) == E_PARAM);
    check("nan_gamma_error",
          decode1(&(Transfer){ .kind = TF_GAMMA, .gamma_set = 1, .gamma = NAN }, 0.5, &v) == E_PARAM);
    check("inf_gamma_error",
          decode1(&(Transfer){ .kind = TF_GAMMA, .gamma_set = 1, .gamma = INFINITY }, 0.5, &v) == E_PARAM);

    /* a partial piecewise curve is rejected as well */
    Transfer partial = { .kind = TF_SRGB, .curve_set = 1, .cut = 0.04045,
                         .slope = 12.92, .scale = 1.055 };
    check("partial_srgb_error", decode1(&partial, 0.5, &v) == E_PARAM);
    check("zero_threshold_error",
          decode1(&(Transfer){ .kind = TF_SRGB, .curve_set = 1, .cut = 0.0,
                               .slope = 12.92, .scale = 1.055, .offset = 0.055,
                               .power = 2.4 }, 0.5, &v) == E_PARAM);

    Transfer lin = { .kind = TF_LINEAR };
    check("linear_ok", decode1(&lin, 0.5, &v) == E_OK && fabs(v - 0.5) < 1e-15);

    /* --- sRGB reference points (D65 / IEC 61966-2-1 constants) --- */
    Transfer srgb = {
        .kind = TF_SRGB, .curve_set = 1,
        .cut = 0.04045, .slope = 12.92, .scale = 1.055,
        .offset = 0.055, .power = 2.4
    };
    check("srgb_decode_half", decode1(&srgb, 0.5, &v) == E_OK &&
                              fabs(v - 0.21404114048223255) < 1e-12);
    check("srgb_encode_half", encode1(&srgb, 0.5, &v) == E_OK &&
                              fabs(v - 0.7353569830524495) < 1e-12);
    check("srgb_white", decode1(&srgb, 1.0, &v) == E_OK && fabs(v - 1.0) < 1e-12);

    /* invalid channels propagate as errors */
    check("negative_channel", decode_rgb(&srgb, (RGB){ 0.5, -0.1, 0.5 }, &(RGB){0,0,0}) == E_ARG);
    check("nan_channel", decode_rgb(&srgb, (RGB){ 0.5, 0.5, NAN }, &(RGB){0,0,0}) == E_ARG);
    check("negative_linear", encode1(&srgb, -0.01, &v) == E_ARG);

    /* --- encoded-space blending is too bright --- */
    {
        double l2, l3, correct;
        decode1(&srgb, 0.2, &l2);
        decode1(&srgb, 0.3, &l3);
        encode1(&srgb, l2 + l3, &correct);
        check("encoded_sum_too_bright", fabs(0.5 - correct) > 0.1);
        printf("note: encoded 0.2+0.3=0.5, linear-correct=%.12g\n", correct);
    }

    /* --- pipeline matches per-pixel exact conversion --- */
    const Sample samples[6] = {
        { { 0.05, 0.05, 0.05 }, 1.0 },
        { { 0.10, 0.10, 0.10 }, 1.0 },
        { { 0.20, 0.20, 0.20 }, 1.0 },
        { { 0.40, 0.40, 0.40 }, 1.0 },
        { { 0.60, 0.60, 0.60 }, 1.0 },
        { { 0.90, 0.90, 0.90 }, 1.0 },
    };
    Pipeline pipe;
    int ok_init = (pipeline_init(&pipe, &srgb) == E_OK);
    check("pipeline_init", ok_init);
    int ok_add = 1;
    for (int i = 0; i < 6; i++) ok_add &= (pipeline_add(&pipe, samples[i]) == E_OK);
    check("pipeline_add", ok_add);
    RGB out, ref;
    check("pipeline_output", pipeline_output(&pipe, &out) == E_OK);
    check("pixelwise_ref_ok", pixelwise_reference(&srgb, samples, 6, &ref) == E_OK);
    check("pipeline_matches_pixelwise",
          fabs(out.r - ref.r) < 1e-12 && fabs(out.g - ref.g) < 1e-12 &&
          fabs(out.b - ref.b) < 1e-12);
    check("pipeline_is_linear_mean", fabs(out.r - 0.5002959330885246) < 1e-12);
    check("pipeline_not_encoded_mean", fabs(out.r - 0.375) > 0.1);
    printf("note: 6-sample output=%.12g (encoded mean would be 0.375)\n", out.r);

    /* weighted samples */
    const Sample weighted[2] = {
        { { 0.1, 0.1, 0.1 }, 1.0 },
        { { 0.9, 0.9, 0.9 }, 3.0 },
    };
    Pipeline pw;
    pipeline_init(&pw, &srgb);
    pipeline_add(&pw, weighted[0]);
    pipeline_add(&pw, weighted[1]);
    RGB wo, wr;
    pipeline_output(&pw, &wo);
    pixelwise_reference(&srgb, weighted, 2, &wr);
    check("weighted_matches_pixelwise", fabs(wo.r - wr.r) < 1e-12);

    /* missing parameters fail at every pipeline entry point */
    Pipeline badpipe;
    check("pipeline_init_missing", pipeline_init(&badpipe, &srgb_missing) == E_PARAM);
    pipeline_init(&badpipe, &srgb);
    {
        Transfer broken = srgb;
        broken.curve_set = 0;
        Pipeline bp2;
        pipeline_init(&bp2, &srgb);
        bp2.tf = broken;
        check("pipeline_add_missing", pipeline_add(&bp2, samples[0]) == E_PARAM);
    }
    check("pixelwise_missing", pixelwise_reference(&srgb_missing, samples, 6, &ref) == E_PARAM);
    check("bound_missing_params", isinf(roundtrip_bound(&srgb_missing).encoded));

    /* --- round-trip error: measured vs analytic upper bound --- */
    RoundtripMeasured meas = measure_roundtrip(&srgb, 2000000);
    RoundtripBound bound = roundtrip_bound(&srgb);
    double meas_de = measure_roundtrip_delta_e(&srgb, 2000000);
    printf("measurement: encoded=%.12g linear=%.12g dE=%.12g\n",
           meas.encoded, meas.linear, meas_de);
    printf("bound:       encoded=%.12g linear=%.12g dE=%.12g\n",
           bound.encoded, bound.linear, bound.delta_e);

    check("roundtrip_encoded_within_bound", meas.encoded <= bound.encoded);
    check("roundtrip_linear_within_bound", meas.linear <= bound.linear);
    check("roundtrip_delta_e_within_bound", meas_de <= bound.delta_e);
    check("roundtrip_encoded_tiny", meas.encoded < 1e-12);
    check("roundtrip_linear_tiny", meas.linear < 1e-8);
    check("roundtrip_delta_e_tiny", bound.delta_e < 1e-3);
    check("roundtrip_encoded_below_8bit", meas.encoded * 255.0 < 0.5);
    check("roundtrip_linear_below_8bit", meas.linear * 255.0 < 0.5);

    /* the bound is in linear units; check an off-grey colour whose channels the
       seam actually perturbs (channels far from the breakpoint round-trip
       exactly, so use values whose encoded forms sit just above the seam). */
    {
        double lb = linear_cut(&srgb);
        RGB a = { lb, 0.25, 0.6 }, ea, eb;
        encode_rgb(&srgb, a, &ea);
        decode_rgb(&srgb, ea, &eb);
        double worst = 0.0;
        const double d[3] = { eb.r - a.r, eb.g - a.g, eb.b - a.b };
        for (int i = 0; i < 3; i++) if (fabs(d[i]) > worst) worst = fabs(d[i]);
        check("colour_roundtrip_linear_within_bound", worst <= bound.linear);
        printf("note: colour round-trip per-channel=%.12g\n", worst);
    }

    /* gamma transfer round-trip (single branch, no seam) */
    Transfer g22 = { .kind = TF_GAMMA, .gamma_set = 1, .gamma = 2.2 };
    RoundtripMeasured gm = measure_roundtrip(&g22, 500000);
    RoundtripBound gb = roundtrip_bound(&g22);
    check("gamma_roundtrip_within_bound",
          gm.encoded <= gb.encoded && gm.linear <= gb.linear);
    check("gamma_roundtrip_tiny", gm.encoded < 1e-12 && gm.linear < 1e-12);

    (void)status_str;
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
