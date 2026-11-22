/* texture: wrap-aware screen-space derivatives (ddx/ddy).  Neighboring
   pixel coordinates are unwrapped per addressing mode so that crossing a
   wrap seam does not turn into a huge derivative (and the wrong mip
   level); clamp differentiates the clamped coordinate instead.  The
   result is checked against the analytic derivative with an error upper
   bound, a seam-pixel count that must be zero, and an invalid mode
   reports the offending value. */
#include <stdio.h>
#include <math.h>
#include <string.h>

enum { WRAP_REPEAT = 0, WRAP_CLAMP = 1 };

#define DERIV_ERR_MAX 1e-5f  /* max |wrap-aware - analytic| per pixel */
#define SEAM_SLOPE    0.5f   /* |dcoord/dpixel| beyond half the texture is a seam */

static float clamp01(float u) {
    return u < 0.0f ? 0.0f : (u > 1.0f ? 1.0f : u);
}

/* coordinate as the sampler sees it under this addressing mode */
static float wrap_coord(int mode, float u) {
    if (mode == WRAP_REPEAT) return u - floorf(u);
    return clamp01(u);
}

/* returns -1 and writes the offending numeric mode for anything unsupported */
static int wrap_mode_error(int mode, char* err, size_t errn) {
    if (mode == WRAP_REPEAT || mode == WRAP_CLAMP) return 0;
    snprintf(err, errn, "invalid wrap mode: %d", mode);
    return -1;
}

/* Screen-space derivative between two adjacent pixels.
   repeat: take the shorter way around the seam, so 0.97 -> 0.03 is +0.06
           instead of -0.94.
   clamp:  differentiate the clamped coordinate, which is zero outside
           [0,1] and partial on the single pixel that crosses the border.
   c0/c1 may be raw or already wrapped coordinates; both are idempotent. */
static int wrap_deriv(int mode, float c0, float c1, float* out, char* err, size_t errn) {
    if (wrap_mode_error(mode, err, errn) != 0) return -1;
    float d = wrap_coord(mode, c1) - wrap_coord(mode, c0);
    if (mode == WRAP_REPEAT) d -= rintf(d);  /* wrap difference into [-0.5, 0.5] */
    *out = d;
    return 0;
}

/* the derivative a shader would get without wrap correction, for counting
   the seams the fix removes */
static float naive_deriv(int mode, float c0, float c1) {
    return wrap_coord(mode, c1) - wrap_coord(mode, c0);
}

typedef struct { float max_err; int seams; int naive_seams; } Scan;

/* Walk one row/column of a linear coordinate field base + slope*i and
   compare the wrap-aware derivative against the analytic value `exact`
   (constant over the whole interval). */
static Scan scan(int mode, float slope, float base, int n, float exact) {
    Scan s = { 0.0f, 0, 0 };
    char err[64] = { 0 };
    for (int i = 0; i < n; i++) {
        float u0 = base + slope * (float)i;
        float u1 = base + slope * (float)(i + 1);
        float d = 0.0f;
        if (wrap_deriv(mode, u0, u1, &d, err, sizeof err) != 0) {
            s.max_err = INFINITY;
            break;
        }
        float e = fabsf(d - exact);
        if (e > s.max_err) s.max_err = e;
        if (fabsf(d) > SEAM_SLOPE) s.seams++;
        if (fabsf(naive_deriv(mode, u0, u1)) > SEAM_SLOPE) s.naive_seams++;
    }
    return s;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

int main(void) {
    check("repeat_maps_unit", wrap_coord(WRAP_REPEAT, 3.25f) >= 0.0f &&
                              wrap_coord(WRAP_REPEAT, 3.25f) < 1.0f);
    check("clamp_bounds", wrap_coord(WRAP_CLAMP, -2.0f) == 0.0f &&
                          wrap_coord(WRAP_CLAMP, 5.0f) == 1.0f);

    /* repeat: a 1/32-per-pixel UV ramp crosses 8 seams over 256 pixels.
       The wrap-aware derivative must recover the slope everywhere and
       flag no seam, while the naive derivative does show the artifact. */
    Scan rr = scan(WRAP_REPEAT, 1.0f / 32.0f, 0.917f, 256, 1.0f / 32.0f);
    printf("repeat_max_err=%.9g repeat_seams=%d repeat_naive_seams=%d\n",
           (double)rr.max_err, rr.seams, rr.naive_seams);
    check("repeat_err_within_bound", rr.max_err <= DERIV_ERR_MAX);
    check("repeat_seams_zero", rr.seams == 0);
    check("repeat_naive_seams_seen", rr.naive_seams > 0);

    /* negative slope through negative, multi-period coordinates */
    Scan rn = scan(WRAP_REPEAT, -1.0f / 24.0f, -3.3f, 240, -1.0f / 24.0f);
    check("repeat_negative_err_within_bound", rn.max_err <= DERIV_ERR_MAX);
    check("repeat_negative_seams_zero", rn.seams == 0);

    /* the same correction on the v axis (ddy) */
    Scan dv = scan(WRAP_REPEAT, 1.0f / 20.0f, 0.43f, 200, 1.0f / 20.0f);
    check("ddy_repeat_err_within_bound", dv.max_err <= DERIV_ERR_MAX);

    /* clamp: interior pixels differentiate like an unwrapped coordinate */
    Scan ci = scan(WRAP_CLAMP, 0.01f, 0.2f, 64, 0.01f);
    printf("clamp_interior_max_err=%.9g clamp_interior_seams=%d\n",
           (double)ci.max_err, ci.seams);
    check("clamp_interior_err_within_bound", ci.max_err <= DERIV_ERR_MAX);
    check("clamp_interior_seams_zero", ci.seams == 0);

    /* clamp: fully outside the texture on either side -> zero derivative */
    Scan cl = scan(WRAP_CLAMP, 0.01f, -0.6f, 40, 0.0f);
    Scan cr = scan(WRAP_CLAMP, 0.01f, 1.2f, 40, 0.0f);
    check("clamp_left_zero", cl.max_err == 0.0f && cl.seams == 0);
    check("clamp_right_zero", cr.max_err == 0.0f && cr.seams == 0);

    /* clamp: the one pixel spanning a border gets a partial, never a spiked,
       derivative */
    char err[64];
    float d = 0.0f;
    int ok_hi = wrap_deriv(WRAP_CLAMP, 0.99f, 1.01f, &d, err, sizeof err) == 0 &&
                d >= 0.0f && d <= 0.02f;
    int ok_lo = wrap_deriv(WRAP_CLAMP, -0.01f, 0.01f, &d, err, sizeof err) == 0 &&
                d >= 0.0f && d <= 0.02f;
    check("clamp_border_partial_high", ok_hi);
    check("clamp_border_partial_low", ok_lo);

    /* invalid modes report the value they were given and leave *out alone */
    d = 123.0f;
    err[0] = 0;
    int bad = wrap_deriv(7, 0.1f, 0.2f, &d, err, sizeof err);
    printf("bad_mode_msg=\"%s\"\n", err);
    check("bad_mode_fails", bad == -1);
    check("bad_mode_reports_value", strstr(err, "7") != NULL);
    check("bad_mode_leaves_output", d == 123.0f);

    err[0] = 0;
    bad = wrap_deriv(-3, 0.1f, 0.2f, &d, err, sizeof err);
    check("negative_mode_reports_value", bad == -1 && strstr(err, "-3") != NULL);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
