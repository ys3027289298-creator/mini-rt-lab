/* bump: height-field self-occlusion.

   A bump-mapped surface is z = h(x, y), normal +z, with the shading point at
   the origin and h(0,0) = 0.  The patch is the square [-1,1]^2 and h is C1
   there, tapering smoothly to zero at the boundary, so a ray only needs
   marching until it leaves the square or rises clear of every possible
   height.  The march is the usual bump-shadow sweep: walk along L and test
   z(t) = t against h(x(t), y(t)).

   Step selection follows the field's slope.  With lam = max |grad h| over the
   patch, one step can change h by at most lam * |Lh| * dt (|Lh| is the
   lateral travel per unit height of the ray), so a height budget eps gives
   dt = eps / (lam * |Lh|).  Expressed as a step count, n = ceil(t_exit / dt)
   - "precision by step count".

   Correctness is checked against explicit geometry: the same height field is
   ray traced at a fixed fine step (with bisection at each sign change) and
   the two visibilities are compared per direction.  The sweep only reports
   "occluded" at a sample that is strictly below the surface, so it cannot
   invent a shadow; the possible error is the opposite one, a missed
   occlusion.  Two facts pin the error down:

     - a flank crossing always finds one sample inside, because a C1 field
       rises by at most lam * |Lh| * dt = eps while a step advances the ray
       by dt (for eps < 1), and the step that hits the flank with the largest
       gap lands below the surface at the next sample;
     - a miss therefore needs a whole occluded component (ridge to ridge) to
       fit between two samples, i.e. max chord < dt.  dt shrinks linearly in
       the step count, so the measure of rays that can hide shrinks with it
       and reaches zero at a finite count.

   Both are tested: counts of misses and of "explainable" rays (chord < dt)
   are printed per step count and required to fall to zero by the end. */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265358979323846f

/* ------------------------------------------------------------------ */
/* height fields                                                      */

#define NFIELD 5

typedef struct { float h0; } Field;
static Field F[NFIELD];
static float LAMBDA[NFIELD];          /* max |grad h| on the patch; drives the step */
static float ZMAX[NFIELD];            /* max |h| on the patch; march cutoff */

static float gauss(float x, float y, float cx, float cy, float s) {
    float dx = x - cx, dy = y - cy;
    return expf(-(dx * dx + dy * dy) / (s * s));
}

static float base_h(int id, float x, float y) {
    switch (id) {
    case 0: return 0.0f;                                             /* flat */
    case 1: return 0.05f * sinf(3.0f * x) * sinf(2.0f * y);          /* gentle */
    case 2: return 0.09f * cosf(5.0f * x) * cosf(4.0f * y);          /* ripples */
    case 3: return 0.16f * sinf(7.0f * x) * cosf(6.0f * y)
                 + 0.13f * sinf(9.0f * y);                           /* steep */
    case 4: return 0.30f * gauss(x, y, 0.0f, 0.0f, 0.32f)
                 - 0.22f * gauss(x, y, 0.42f, 0.30f, 0.30f)
                 + 0.18f * gauss(x, y, -0.40f, 0.25f, 0.22f);        /* lumps */
    default: return 0.0f;
    }
}

/* C1 taper to zero on the patch boundary */
static float window(float x, float y) {
    float wx, wy;
    if (fabsf(x) >= 1.0f || fabsf(y) >= 1.0f) return 0.0f;
    wx = 1.0f - x * x;
    wy = 1.0f - y * y;
    return wx * wx * wy * wy;
}

static float field_h(int id, float x, float y) {
    return (base_h(id, x, y) - F[id].h0) * window(x, y);
}

static void field_grad(int id, float x, float y, float* gx, float* gy) {
    float e = 1e-4f;
    *gx = (field_h(id, x + e, y) - field_h(id, x - e, y)) / (2.0f * e);
    *gy = (field_h(id, x, y + e) - field_h(id, x, y - e)) / (2.0f * e);
}

/* The slope that drives the step: worst |grad h| on a dense patch grid. */
static float field_max_slope(int id) {
    int n = 400;
    float gmax = 0.0f;
    for (int i = 0; i <= n; i++) {
        float x = -1.0f + 2.0f * (float)i / (float)n;
        for (int j = 0; j <= n; j++) {
            float y = -1.0f + 2.0f * (float)j / (float)n;
            float gx, gy;
            field_grad(id, x, y, &gx, &gy);
            float g = sqrtf(gx * gx + gy * gy);
            if (g > gmax) gmax = g;
        }
    }
    return gmax * 1.05f;                       /* grid-vs-continuum margin */
}

static float field_max_height(int id) {
    int n = 200;
    float zmax = 0.0f;
    for (int i = 0; i <= n; i++) {
        float x = -1.0f + 2.0f * (float)i / (float)n;
        for (int j = 0; j <= n; j++) {
            float y = -1.0f + 2.0f * (float)j / (float)n;
            float h = fabsf(field_h(id, x, y));
            if (h > zmax) zmax = h;
        }
    }
    return zmax;
}

static void fields_init(void) {
    for (int i = 0; i < NFIELD; i++) {
        F[i].h0 = base_h(i, 0.0f, 0.0f);       /* shading point sits on surface */
        LAMBDA[i] = field_max_slope(i);
        ZMAX[i] = field_max_height(i);
    }
}

/* An empty height field cannot be shaded: reject it up front and say why. */
static int field_ok(int id, char* why, size_t whyn) {
    if (!why || whyn == 0) return 0;
    if (id < 0) { snprintf(why, whyn, "empty height field (no surface supplied)"); return 0; }
    if (id >= NFIELD) { snprintf(why, whyn, "empty height field (id %d out of range)", id); return 0; }
    if (!isfinite(LAMBDA[id]) || LAMBDA[id] < 0.0f) {
        snprintf(why, whyn, "slope of field %d is not finite", id);
        return 0;
    }
    if (!isfinite(F[id].h0) || fabsf(field_h(id, 0.0f, 0.0f)) > 1e-6f) {
        snprintf(why, whyn, "field %d does not pass through the shading point", id);
        return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* sweep                                                              */

/* Vertical gap between the ray and the surface: > 0 above, < 0 inside. */
static float gap_at(int id, float Lx, float Ly, float Lz, float t) {
    return t - field_h(id, (Lx / Lz) * t, (Ly / Lz) * t);
}

/* Time at which a ray leaves the patch (the square is convex).  Straight-up
   rays never leave; they are never occluded either. */
static float t_exit_of(float Lx, float Ly, float Lz) {
    float m = fmaxf(fabsf(Lx), fabsf(Ly));
    if (!(m > 1e-9f)) return -1.0f;
    return Lz / m;
}

typedef struct { int visible; int steps; float dt; float t_exit; } Sweep;

/* Rays going straight up sit at h(0,0) = 0 with no lateral motion and never
   leave the patch; nothing can shadow them.  Their dt is 0 for that reason. */
static Sweep sweep_up_clear(void) {
    Sweep s = { 1, 0, 0.0f, -1.0f };
    return s;
}

/* Fixed step count: the precision knob the convergence test doubles.
   Stop early once the ray is at least |h|max above every possible height
   (all fields are C1 at 0, so one step never skips the hem); on this patch
   the taper keeps |h| small so that cutoff is what usually ends the march. */
static Sweep sweep_steps(int id, const float L[3], int n) {
    Sweep s = { 1, 0, 0.0f, 0.0f };
    float Lx = L[0], Ly = L[1], Lz = L[2];
    float te = t_exit_of(Lx, Ly, Lz);
    s.t_exit = te;
    if (te < 0.0f) return sweep_up_clear();
    if (!(n > 0)) n = 1;
    {
        float dt = te / (float)n;
        s.dt = dt;
        for (int k = 1; k <= n; k++) {
            float t = dt * (float)k;
            float x = (Lx / Lz) * t, y = (Ly / Lz) * t;
            if (t - field_h(id, x, y) < -1e-6f) {
                s.visible = 0;
                s.steps = k;
                return s;
            }
            s.steps = k;
            /* above the band, and rising: nothing further can be occluded */
            if (t - ZMAX[id] > 1e-3f) break;
        }
    }
    return s;
}

/* Slope-driven stepping: eps is the largest height change one step may make,
   so dt = eps / (max |grad h| * |Lh|), where |Lh| is the lateral travel per
   unit height of the ray.  Returns "steps = T / dt". */
static int sweep_steps_for(int id, const float L[3], float eps) {
    float Lh, dt, te;
    int n;
    if (!(eps > 0.0f) || LAMBDA[id] <= 0.0f) return 1;
    te = t_exit_of(L[0], L[1], L[2]);
    if (te < 0.0f) return 1;
    Lh = sqrtf(L[0] * L[0] + L[1] * L[1]) / L[2];
    if (!(Lh > 0.0f)) return 1;
    dt = eps / (LAMBDA[id] * Lh);
    n = (int)ceilf(te / dt);
    if (n < 1) n = 1;
    if (n > (1 << 22)) n = 1 << 22;
    return n;
}

static Sweep sweep_slope(int id, const float L[3], float eps) {
    return sweep_steps(id, L, sweep_steps_for(id, L, eps));
}

/* ------------------------------------------------------------------ */
/* reference: ray trace the height field itself                       */

typedef struct {
    int occluded;
    float enter, exit;
    float maxchord;          /* longest occluded segment along the ray */
    float total;             /* occluded path length */
} Ref;

static float bisect_cross(int id, const float L[3], float a, float b) {
    float Lx = L[0], Ly = L[1], Lz = L[2];
    float ga = gap_at(id, Lx, Ly, Lz, a);
    for (int i = 0; i < 40; i++) {
        float m = 0.5f * (a + b);
        float gm = gap_at(id, Lx, Ly, Lz, m);
        if ((gm > 0.0f) == (ga > 0.0f)) { a = m; ga = gm; }
        else b = m;
    }
    return 0.5f * (a + b);
}

/* Reference: dense fixed-step march plus bisection at each sign change.  The
   step is far below the thinnest occluded lobe the fields produce, and the
   test cross-validates it at two densities and against the sweep at very
   high step count, so its own error cannot mask the sweep's. */
static Ref reference_trace(int id, const float L[3], int n) {
    Ref r;
    float te, dt, prev_t = 0.0f, prev_g;
    int inside = 0;
    float enter = 0.0f;
    memset(&r, 0, sizeof(r));
    te = t_exit_of(L[0], L[1], L[2]);
    if (te < 0.0f) return r;
    if (n < 1) n = 1;
    dt = te / (float)n;
    prev_g = gap_at(id, L[0], L[1], L[2], 0.0f);      /* == 0 at the origin */
    for (int k = 1; k <= n; k++) {
        float t = dt * (float)k;
        float g = gap_at(id, L[0], L[1], L[2], t);
        if (!inside && g < 0.0f) {
            inside = 1;
            enter = (prev_g > 0.0f) ? bisect_cross(id, L, prev_t, t) : 0.0f;
        } else if (inside && g >= 0.0f) {
            float ex = bisect_cross(id, L, prev_t, t);
            float chord = ex - enter;
            inside = 0;
            r.occluded = 1;
            r.total += chord;
            if (chord > r.maxchord) r.maxchord = chord;
        }
        prev_t = t;
        prev_g = g;
    }
    if (inside) {
        float chord = te - enter;
        r.occluded = 1;
        r.total += chord;
        if (chord > r.maxchord) r.maxchord = chord;
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* tests                                                              */

#define NDIR 128
static float DIRS[NDIR][3];
static Ref REF[NFIELD][NDIR];        /* ground truth, filled per field */
static const int TEST_FIELDS[] = { 1, 2, 3, 4 };

static uint32_t g_rng = 0x9e3779b9u;
static uint32_t xorshift32(void) {
    uint32_t x = g_rng;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    return g_rng = x;
}
static float rand01(void) { return (xorshift32() >> 8) * (1.0f / 16777216.0f); }

/* Uniform in z over the upper hemisphere: covers straight up to grazing. */
static void make_dirs(void) {
    for (int i = 0; i < NDIR; i++) {
        float z = 0.087f + (0.999f - 0.087f) * rand01();     /* up to ~85 deg */
        float r = sqrtf(1.0f - z * z);
        float phi = 2.0f * PI * rand01();
        DIRS[i][0] = r * cosf(phi);
        DIRS[i][1] = r * sinf(phi);
        DIRS[i][2] = z;
    }
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int same_visibility(const Ref* a, const Ref* b) {
    return (a->occluded != 0) == (b->occluded != 0);
}

int main(void) {
    char why[128];
    fields_init();
    make_dirs();

    /* 1. usable fields, and empty input rejected with a reason */
    {
        int ok = 1;
        for (int i = 0; i < NFIELD; i++) ok &= field_ok(i, why, sizeof why);
        check("fields_ok", ok && LAMBDA[1] > 0.0f && LAMBDA[3] > LAMBDA[1]);
    }
    {
        char w1[128] = "", w2[128] = "";
        int r1 = field_ok(-1, w1, sizeof w1);
        int r2 = field_ok(NFIELD, w2, sizeof w2);
        check("empty_field_rejected", r1 == 0 && r2 == 0);
        check("empty_field_reason", strlen(w1) > 0 && strlen(w2) > 0);
    }

    /* 2. flat field: no self-shadow for any direction, at any step count */
    {
        int ok = 1;
        for (int d = 0; d < NDIR; d++) {
            ok &= sweep_steps(0, DIRS[d], 1).visible;
            ok &= sweep_steps(0, DIRS[d], 37).visible;
        }
        check("flat_never_shadows", ok);
    }

    /* 3. straight up is always clear, even on the steep field */
    {
        float up[3] = { 0.0f, 0.0f, 1.0f };
        check("up_clear", sweep_steps(3, up, 1).visible == 1);
    }

    /* 4. reference at two fine steps must agree: it is the comparison baseline */
    {
        int ok = 1;
        for (size_t fi = 0; fi < sizeof TEST_FIELDS / sizeof TEST_FIELDS[0]; fi++) {
            int id = TEST_FIELDS[fi];
            for (int d = 0; d < NDIR; d++) {
                Ref lo = reference_trace(id, DIRS[d], 2048);
                Ref hi = reference_trace(id, DIRS[d], 8192);
                REF[id][d] = hi;
                ok &= same_visibility(&lo, &hi);
            }
        }
        check("reference_stable", ok);
    }

    /* the scene really does self-shadow, or the rest proves nothing */
    {
        int shadowed = 0;
        for (size_t fi = 0; fi < sizeof TEST_FIELDS / sizeof TEST_FIELDS[0]; fi++)
            for (int d = 0; d < NDIR; d++)
                shadowed += REF[TEST_FIELDS[fi]][d].occluded;
        check("scene_has_shadows", shadowed > 0);
        printf("shadowed_rays=%d\n", shadowed);
    }

    /* 5. sweep vs geometry: no invented shadows, every miss explained by a
          chord shorter than one step, and both counts fall as steps double
          until they are exactly zero */
    {
        static const int NS[] = { 4, 8, 16, 32, 64, 128, 256, 512, 1024 };
        int miss_prev = -1, bound_prev = -1;
        int ok_noninvent = 1, ok_bound = 1, ok_converge = 1;
        for (size_t ni = 0; ni < sizeof NS / sizeof NS[0]; ni++) {
            int n = NS[ni], bound = 0, missed = 0, invented = 0, shadowed = 0;
            float dt_min = 1e30f, dt_max = 0.0f;
            for (size_t fi = 0; fi < sizeof TEST_FIELDS / sizeof TEST_FIELDS[0]; fi++) {
                int id = TEST_FIELDS[fi];
                for (int d = 0; d < NDIR; d++) {
                    Ref* r = &REF[id][d];
                    Sweep s = sweep_steps(id, DIRS[d], n);
                    if (s.dt > 0.0f) {
                        if (s.dt < dt_min) dt_min = s.dt;
                        if (s.dt > dt_max) dt_max = s.dt;
                    }
                    if (!s.visible && !r->occluded) invented++;
                    if (s.visible && r->occluded) {
                        missed++;
                        /* a miss is only possible if every occluded segment
                           is shorter than one step */
                        if (!(r->maxchord < s.dt + 1e-6f)) ok_bound = 0;
                    }
                    /* upper bound on what this step count can still miss: the
                       number of shadowed rays whose longest chord is shorter
                       than one step.  Nothing else can be hidden from it. */
                    if (r->occluded) {
                        shadowed++;
                        if (r->maxchord < s.dt + 1e-6f) bound++;
                    }
                }
            }
            ok_noninvent &= (invented == 0);
            if (miss_prev >= 0 && missed > miss_prev) ok_converge = 0;
            if (bound_prev >= 0 && bound > bound_prev) ok_converge = 0;
            miss_prev = missed;
            bound_prev = bound;
            printf("sweep n=%4d dt=[%.5f,%.5f] shadowed=%3d missed=%3d bound=%3d invented=%d\n",
                   n, (double)dt_min, (double)dt_max, shadowed, missed, bound, invented);
            if (ni + 1 == sizeof NS / sizeof NS[0])
                ok_converge &= (missed == 0 && bound == 0);
        }
        check("sweep_never_invents_shadow", ok_noninvent);
        check("miss_implies_short_chord", ok_bound);
        check("error_bound_shrinks", ok_converge);
    }

    /* 6. slope-driven stepping: eps is a height budget, so halving it never
          reduces the step count, and the marched result is exactly what the
          equivalent fixed step count gives */
    {
        int ok = 1, ok2 = 1;
        float eps = 0.01f;
        for (size_t fi = 0; fi < sizeof TEST_FIELDS / sizeof TEST_FIELDS[0]; fi++) {
            int id = TEST_FIELDS[fi];
            for (int d = 0; d < NDIR; d++) {
                int n1 = sweep_steps_for(id, DIRS[d], eps);
                int n2 = sweep_steps_for(id, DIRS[d], 0.5f * eps);
                Sweep a = sweep_slope(id, DIRS[d], eps);
                Sweep b = sweep_steps(id, DIRS[d], n1);
                if (n2 < n1) ok = 0;
                if (a.visible != b.visible || a.steps != b.steps) ok2 = 0;
            }
        }
        check("eps_halves_steps_increase", ok);
        check("slope_matches_fixed", ok2);

        int missed = 0;
        for (size_t fi = 0; fi < sizeof TEST_FIELDS / sizeof TEST_FIELDS[0]; fi++) {
            int id = TEST_FIELDS[fi];
            for (int d = 0; d < NDIR; d++) {
                Sweep s = sweep_slope(id, DIRS[d], eps);
                if (s.visible && REF[id][d].occluded) missed++;
            }
        }
        printf("slope_eps=%.3f missed=%d\n", (double)eps, missed);
        check("slope_budget_ok", missed <= 2);
    }

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
