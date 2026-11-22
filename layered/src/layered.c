
/* layered: multi-layer dielectric stacks.
   Light crosses a planar stack media[0] | media[1] ... | media[n-1], where
   media[0] and media[n-1] are semi-infinite half-spaces and each interior
   medium has a physical thickness. Per RGB channel the ray is handled with
   that channel's IOR and absorption coefficient.

   Two solvers are provided:
   1. layered_trace: deterministic path-split transport. At every interface the
      ray refracts/reflects per Snell's law, energy is split by the unpolarized
      Fresnel factor, and each interior crossing picks up Beer-Lambert
      attenuation exp(-sigma * d / cos_theta). Splitting stops at max_depth and
      the weight of unresolved paths is returned as a certified upper bound on
      the truncation error versus the infinite-bounce answer.
   2. layered_transport: analytic incoherent (net-radiation) recurrence that
      sums every multiple-reflection bounce, used as the reference answer.

   Energy accounting: absorption is tallied explicitly at every attenuation
   event, so the tracer obeys T + R + A + unresolved = 1 exactly, and both
   solvers satisfy T + R <= 1 (no more energy leaves than entered). A stack with
   fewer than two media (zero layers) is rejected with LAYERED_EMPTY. */
#include <stdio.h>
#include <math.h>
#include <string.h>

#define LAYERED_MAX_MEDIA 16

enum { LAYERED_OK = 0, LAYERED_EMPTY = -1, LAYERED_CAP = -2, LAYERED_PARAM = -3 };

typedef struct { double eta[3]; double sigma[3]; } Medium;
typedef struct {
    Medium m[LAYERED_MAX_MEDIA];
    double d[LAYERED_MAX_MEDIA]; /* thickness of interior medium j; 0/last unused */
    int n;
} Stack;
typedef struct { double T[3], R[3], A[3]; } Transport;
typedef struct { double T[3], R[3], A[3], unresolved[3]; } TraceResult;

/* ---- Fresnel (unpolarized intensity, non-magnetic dielectrics) ---- */
static double fresnel(double cos_i, double eta_i, double eta_t, int *tir) {
    double sin2_t = (eta_i * eta_i) / (eta_t * eta_t) * (1.0 - cos_i * cos_i);
    if (sin2_t >= 1.0) { *tir = 1; return 1.0; } /* total internal reflection */
    *tir = 0;
    double cos_t = sqrt(1.0 - sin2_t);
    double rs = (eta_i * cos_i - eta_t * cos_t) / (eta_i * cos_i + eta_t * cos_t);
    double rp = (eta_i * cos_t - eta_t * cos_i) / (eta_i * cos_t + eta_t * cos_i);
    return 0.5 * (rs * rs + rp * rp);
}

static int valid_medium(const Medium *m) {
    for (int c = 0; c < 3; c++)
        if (!(m->eta[c] > 0.0) || m->sigma[c] < 0.0) return 0;
    return 1;
}

int layered_from_media(Stack *s, const Medium *media, const double *thick, int nmedia) {
    if (nmedia < 2) return LAYERED_EMPTY;
    if (nmedia > LAYERED_MAX_MEDIA) return LAYERED_CAP;
    for (int i = 0; i < nmedia; i++) {
        if (!valid_medium(&media[i])) return LAYERED_PARAM;
        if (i > 0 && i < nmedia - 1 && thick[i] < 0.0) return LAYERED_PARAM;
    }
    memset(s, 0, sizeof(*s));
    s->n = nmedia;
    for (int i = 0; i < nmedia; i++) {
        s->m[i] = media[i];
        s->d[i] = (i > 0 && i < nmedia - 1) ? thick[i] : 0.0;
    }
    return LAYERED_OK;
}

/* Convenience builder for nslab solid slabs separated by a gap medium:
   ambient | slab0 | gap | slab1 | ... | slab(n-1) | ambient.
   nslab == 0 is a zero-layer request and is reported, not silently accepted. */
int layered_from_slabs(Stack *s, const Medium *ambient, const Medium *gap, double gap_d,
                       const Medium *slabs, const double *slab_d, int nslab) {
    if (nslab <= 0) return LAYERED_EMPTY;
    if (2 * nslab + 1 > LAYERED_MAX_MEDIA) return LAYERED_CAP;
    if (gap_d < 0.0 || !valid_medium(ambient) || !valid_medium(gap)) return LAYERED_PARAM;
    for (int i = 0; i < nslab; i++) {
        if (slab_d[i] < 0.0 || !valid_medium(&slabs[i])) return LAYERED_PARAM;
    }
    memset(s, 0, sizeof(*s));
    s->m[0] = *ambient;
    for (int k = 0; k < nslab; k++) {
        s->m[2 * k + 1] = slabs[k];
        s->d[2 * k + 1] = slab_d[k];
        s->m[2 * k + 2] = *gap;
        s->d[2 * k + 2] = gap_d;
    }
    s->m[2 * nslab] = *ambient; /* overwrite the trailing gap slot */
    s->n = 2 * nslab + 1;
    return LAYERED_OK;
}

/* Per-channel geometry: cosine in every medium, forward/backward Fresnel at
   every interface, and one-way Beer-Lambert attenuation of every medium. */
static int setup_channel(const Stack *s, int c, double cos0,
                         double *ca, double *ff, double *fb, double *tau) {
    if (!(cos0 > 0.0) || cos0 > 1.0 + 1e-15) return LAYERED_PARAM;
    int M = s->n;
    ca[0] = cos0;
    int tir_at = -1;
    for (int k = 0; k < M - 1; k++) {
        int tir;
        double ei = s->m[k].eta[c], et = s->m[k + 1].eta[c];
        double F = fresnel(ca[k], ei, et, &tir);
        if (tir) { tir_at = k; ff[k] = fb[k] = 1.0; ca[k + 1] = 0.0; break; }
        ff[k] = F;
        double sin2_t = (ei * ei) / (et * et) * (1.0 - ca[k] * ca[k]);
        ca[k + 1] = sqrt(1.0 - sin2_t);
        /* back side: ray inside medium k+1 at the refracted angle; reciprocity
           makes the unpolarized Fresnel intensity factor identical both ways */
        fb[k] = fresnel(ca[k + 1], et, ei, &tir);
    }
    if (tir_at >= 0) { /* unreachable interfaces beyond the TIR point */
        for (int k = tir_at + 1; k < M - 1; k++) { ff[k] = 1.0; fb[k] = 1.0; }
    }
    for (int j = 0; j < M; j++) {
        if (j == 0 || j == M - 1 || (tir_at >= 0 && j > tir_at)) { tau[j] = 1.0; continue; }
        tau[j] = exp(-s->m[j].sigma[c] * s->d[j] / ca[j]);
    }
    return LAYERED_OK;
}

/* Analytic incoherent stack: effective (R, T) of the sub-stack behind each
   interface, accumulated from the last interface backward. One round trip
   through medium k+1 contributes tau^2 and a (1-F)^2 interface factor. */
int layered_transport(const Stack *s, double cos0, Transport *out) {
    if (s->n < 2) return LAYERED_EMPTY;
    int M = s->n;
    for (int c = 0; c < 3; c++) {
        double ca[LAYERED_MAX_MEDIA], ff[LAYERED_MAX_MEDIA - 1];
        double fb[LAYERED_MAX_MEDIA - 1], tau[LAYERED_MAX_MEDIA];
        int rc = setup_channel(s, c, cos0, ca, ff, fb, tau);
        if (rc != LAYERED_OK) return rc;

        double Re = ff[M - 2];             /* last interface: half-space behind */
        double Te = 1.0 - ff[M - 2];
        for (int k = M - 3; k >= 0; k--) {
            double F = ff[k], t = tau[k + 1], t2 = t * t;
            /* round trip behind interface k: reflect here (F), cross medium
               k+1 twice (t^2), reflect off the sub-stack (Re) */
            double den = 1.0 - F * t2 * Re;
            if (den < 1e-15) { Re = 1.0; Te = 0.0; continue; } /* degenerate 0/0 */
            double Rn = F + (1.0 - F) * (1.0 - F) * t2 * Re / den;
            double Tn = (1.0 - F) * t * Te / den;
            Re = Rn; Te = Tn;
        }
        if (Re < 0.0 && Re > -1e-12) Re = 0.0;
        if (Te < 0.0 && Te > -1e-12) Te = 0.0;
        out->R[c] = Re;
        out->T[c] = Te;
        out->A[c] = 1.0 - Re - Te;
    }
    return LAYERED_OK;
}

/* ---- deterministic path-split tracer ----
   State machine on media, not interfaces: fwd(j)/bck(j) hold a ray at one end
   of medium j about to traverse it. Crossing an interior medium applies
   Beer-Lambert attenuation and books the lost weight as absorption immediately,
   so T + R + A + unresolved = 1 at every step. */
typedef struct {
    const Stack *s;
    int depth_limit;
    double ff[LAYERED_MAX_MEDIA - 1];
    double fb[LAYERED_MAX_MEDIA - 1], tau[LAYERED_MAX_MEDIA];
    double T, R, A, unresolved;
} Ctx;

static void bck(Ctx *x, int j, double w, int depth);

/* traverse medium j forward; j == M-1 exits into the right half-space */
static void fwd(Ctx *x, int j, double w, int depth) {
    if (depth > x->depth_limit) { x->unresolved += w; return; }
    if (j == x->s->n - 1) { x->T += w; return; }
    double t = x->tau[j];
    x->A += w * (1.0 - t);
    double w1 = w * t;
    double F = x->ff[j];
    bck(x, j, w1 * F, depth + 1);                 /* reflect: stay in medium j, head back */
    fwd(x, j + 1, w1 * (1.0 - F), depth + 1);     /* transmit: enter medium j+1 */
}

/* traverse medium j backward; j == 0 exits into the left half-space */
static void bck(Ctx *x, int j, double w, int depth) {
    if (depth > x->depth_limit) { x->unresolved += w; return; }
    if (j == 0) { x->R += w; return; }
    double t = x->tau[j];
    x->A += w * (1.0 - t);
    double w1 = w * t;
    double F = x->fb[j - 1];                       /* interface on the left of medium j */
    fwd(x, j, w1 * F, depth + 1);                  /* reflect: stay in medium j, head forward */
    bck(x, j - 1, w1 * (1.0 - F), depth + 1);      /* transmit: enter medium j-1 */
}

int layered_trace(const Stack *s, double cos0, int max_depth, TraceResult *out) {
    if (s->n < 2) return LAYERED_EMPTY;
    if (max_depth < 0) return LAYERED_PARAM;
    for (int c = 0; c < 3; c++) {
        Ctx x;
        memset(&x, 0, sizeof(x));
        x.s = s; x.depth_limit = max_depth;
        double ca[LAYERED_MAX_MEDIA];
        int rc = setup_channel(s, c, cos0, ca, x.ff, x.fb, x.tau);
        if (rc != LAYERED_OK) return rc;
        fwd(&x, 0, 1.0, 0);
        out->T[c] = x.T; out->R[c] = x.R;
        out->A[c] = x.A; out->unresolved[c] = x.unresolved;
    }
    return LAYERED_OK;
}

/* ---------------- self tests ---------------- */
static int failures = 0;
static void check(const char *name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static Medium mk(double eta, double sigma) {
    Medium m;
    for (int c = 0; c < 3; c++) { m.eta[c] = eta; m.sigma[c] = sigma; }
    return m;
}
static Medium mk3(double eta, double s0, double s1, double s2) {
    Medium m;
    m.eta[0] = m.eta[1] = m.eta[2] = eta;
    m.sigma[0] = s0; m.sigma[1] = s1; m.sigma[2] = s2;
    return m;
}

static unsigned long long rng_state = 0x123456789abcdef0ULL;
static double randu(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return (double)(rng_state >> 11) / 9007199254740992.0;
}

int main(void) {
    /* 1. zero layers must be reported, not treated as vacuum */
    Stack empty;
    memset(&empty, 0, sizeof(empty));
    Medium amb = mk(1.0, 0.0), air = mk(1.0, 0.0);
    Medium one_slab = mk(1.5, 0.0);
    double one_d = 0.5;
    check("zero_slab_reported",
          layered_from_slabs(&empty, &amb, &air, 0.01, &one_slab, &one_d, 0) == LAYERED_EMPTY);
    Transport tr;
    TraceResult ttr;
    check("empty_transport_reported", layered_transport(&empty, 1.0, &tr) == LAYERED_EMPTY);
    check("empty_trace_reported", layered_trace(&empty, 1.0, 8, &ttr) == LAYERED_EMPTY);
    Medium solo = mk(1.3, 0.0);
    double solo_d = 0.0;
    check("single_medium_reported", layered_from_media(&empty, &solo, &solo_d, 1) == LAYERED_EMPTY);

    /* 2. parameter validation */
    Medium bad = mk(0.0, 0.0);
    check("bad_ior_reported", layered_from_media(&empty, (Medium[]){ bad, air },
                                                  (double[]){ 0, 0 }, 2) == LAYERED_PARAM);
    Stack good;
    Medium pair[2] = { mk(1.0, 0.0), mk(1.5, 0.0) };
    layered_from_media(&good, pair, (double[]){ 0, 0 }, 2);
    check("bad_angle_reported", layered_transport(&good, 0.0, &tr) == LAYERED_PARAM);

    /* 3. Fresnel golden values: normal incidence air->glass F = 0.04 */
    int tir;
    double F = fresnel(1.0, 1.0, 1.5, &tir);
    check("fresnel_normal", !tir && fabs(F - 0.04) < 1e-12);
    /* oblique reciprocity: 30 deg air->glass, backward ray at the refracted angle */
    double cos_i = cos(30.0 * 3.14159265358979323846 / 180.0);
    double sin2_t = (1.0 / 2.25) * (1.0 - cos_i * cos_i);
    double cos_t = sqrt(1.0 - sin2_t);
    /* backward ray travels glass -> air, meeting the interface at cos_t */
    double Fback = fresnel(cos_t, 1.5, 1.0, &tir);
    check("fresnel_reciprocity", !tir && fabs(Fback - fresnel(cos_i, 1.0, 1.5, &tir)) < 1e-12);

    /* 4. bare interface: T = 1 - F, energy exactly conserved */
    layered_transport(&good, 1.0, &tr);
    check("bare_interface_T", fabs(tr.T[0] - 0.96) < 1e-10);
    check("bare_interface_R", fabs(tr.R[0] - 0.04) < 1e-10);
    check("bare_interface_energy", fabs(tr.T[0] + tr.R[0] - 1.0) < 1e-10);
    layered_trace(&good, 1.0, 12, &ttr);
    check("bare_trace_identity",
          fabs(ttr.T[0] + ttr.R[0] + ttr.A[0] + ttr.unresolved[0] - 1.0) < 1e-12);
    check("bare_trace_bound", fabs(ttr.T[0] - tr.T[0]) <= ttr.unresolved[0] + 1e-12);

    /* 5. single non-absorbing slab, incoherent golden value T = (1-F)/(1+F) */
    Medium slab_med[3] = { mk(1.0, 0.0), mk(1.5, 0.0), mk(1.0, 0.0) };
    Stack slab;
    layered_from_media(&slab, slab_med, (double[]){ 0, 0.5, 0 }, 3);
    layered_transport(&slab, 1.0, &tr);
    double gold_T = 0.96 / 1.04, gold_R = 0.08 / 1.04;
    check("slab_T_golden", fabs(tr.T[0] - gold_T) < 1e-10);
    check("slab_R_golden", fabs(tr.R[0] - gold_R) < 1e-10);
    check("slab_no_absorption", fabs(tr.T[0] + tr.R[0] - 1.0) < 1e-10);
    layered_trace(&slab, 1.0, 16, &ttr);
    check("slab_trace_identity",
          fabs(ttr.T[0] + ttr.R[0] + ttr.A[0] + ttr.unresolved[0] - 1.0) < 1e-12);
    check("slab_trace_bound", fabs(ttr.T[0] - tr.T[0]) <= ttr.unresolved[0] + 1e-12);
    check("slab_trace_one_sided", tr.T[0] - ttr.T[0] >= -1e-12 &&
                                  tr.T[0] - ttr.T[0] <= ttr.unresolved[0] + 1e-12);

    /* 6. absorbing slab: energy leaves as heat, thicker slab transmits less */
    Medium abs_med[3] = { mk(1.0, 0.0), mk(1.5, 2.0), mk(1.0, 0.0) };
    Stack abs1, abs2;
    layered_from_media(&abs1, abs_med, (double[]){ 0, 1.0, 0 }, 3);
    layered_from_media(&abs2, abs_med, (double[]){ 0, 2.0, 0 }, 3);
    Transport tr2;
    layered_transport(&abs1, 1.0, &tr);
    layered_transport(&abs2, 1.0, &tr2);
    check("absorption_loses_energy", tr.T[0] + tr.R[0] < 1.0 && tr.A[0] > 0.0);
    check("thicker_darker", tr2.T[0] < tr.T[0]);
    layered_trace(&abs1, 1.0, 16, &ttr);
    check("abs_trace_identity",
          fabs(ttr.T[0] + ttr.R[0] + ttr.A[0] + ttr.unresolved[0] - 1.0) < 1e-12);
    check("abs_trace_bound", fabs(ttr.T[0] - tr.T[0]) <= ttr.unresolved[0] + 1e-12);
    check("abs_trace_heat_match", fabs(ttr.A[0] - tr.A[0]) <= ttr.unresolved[0] + 1e-12);
    /* 7. total internal reflection: eta 2 -> 1 at cos 0.5, nothing transmits */
    Medium tir_med[2] = { mk(2.0, 0.0), mk(1.0, 0.0) };
    Stack tir_stack;
    layered_from_media(&tir_stack, tir_med, (double[]){ 0, 0 }, 2);
    layered_transport(&tir_stack, 0.5, &tr);
    check("tir_blocks_transmission", tr.T[0] == 0.0 && fabs(tr.R[0] - 1.0) < 1e-12);
    layered_trace(&tir_stack, 0.5, 12, &ttr);
    check("tir_trace", ttr.T[0] == 0.0 && fabs(ttr.R[0] - 1.0) < 1e-12);

    /* 8. tinted three-pane glass: wavelength-dependent sigma, oblique incidence */
    Medium glass = mk3(1.5, 0.05, 0.25, 0.60);
    Medium panes[3] = { glass, glass, glass };
    double pane_d[3] = { 0.8, 0.8, 0.8 };
    Stack windows;
    int rc = layered_from_slabs(&windows, &amb, &air, 0.05, panes, pane_d, 3);
    check("three_pane_build", rc == LAYERED_OK && windows.n == 7);
    layered_transport(&windows, 0.9, &tr);
    layered_trace(&windows, 0.9, 16, &ttr);
    int energy_ok = 1, bound_ok = 1, identity_ok = 1;
    for (int c = 0; c < 3; c++) {
        printf("info_pane_c%d T=%.6f trace=%.6f err=%.2e bound=%.2e R=%.6f A=%.6f\n",
               c, tr.T[c], ttr.T[c], fabs(ttr.T[c] - tr.T[c]), ttr.unresolved[c],
               tr.R[c], tr.A[c]);
        if (!(tr.T[c] >= 0.0 && tr.R[c] >= 0.0 && tr.T[c] + tr.R[c] <= 1.0 + 1e-12))
            energy_ok = 0;
        if (fabs(ttr.T[c] - tr.T[c]) > ttr.unresolved[c] + 1e-12) bound_ok = 0;
        if (fabs(ttr.T[c] + ttr.R[c] + ttr.A[c] + ttr.unresolved[c] - 1.0) > 1e-10)
            identity_ok = 0;
    }
    check("tint_order_red_to_blue", tr.T[0] > tr.T[1] && tr.T[1] > tr.T[2]);
    check("three_pane_energy", energy_ok);
    check("three_pane_bound", bound_ok);
    check("three_pane_identity", identity_ok);

    /* 9. certified bound tightness: on a lossless slab the unresolved weight is
       the geometric tail of the round-trip series (ratio F^2), so it decays
       geometrically with depth and dominates the actual error at every depth */
    int conv = 1;
    double prev_un = 1.0;
    for (int depth = 2; depth <= 14; depth += 2) {
        TraceResult tt;
        layered_trace(&slab, 1.0, depth, &tt);
        double err = fabs(tt.T[0] - gold_T);
        if (err > tt.unresolved[0] + 1e-12) conv = 0;            /* bound holds at every depth */
        if (tt.unresolved[0] >= prev_un) conv = 0;               /* bound shrinks each step */
        if (depth > 2 && tt.unresolved[0] > prev_un * 0.25) conv = 0; /* ratio <= F^2 + slack */
        prev_un = tt.unresolved[0];
    }
    check("bound_converges_geometric", conv);
    TraceResult deep;
    layered_trace(&slab, 1.0, 40, &deep);
    check("deep_trace_matches_analytic", fabs(deep.T[0] - gold_T) < 1e-10 &&
                                         deep.unresolved[0] < 1e-10);
    TraceResult shallow;
    layered_trace(&slab, 1.0, 0, &shallow);
    check("zero_depth_all_unresolved",
          fabs(shallow.unresolved[0] - 1.0) < 1e-12 && shallow.T[0] == 0.0);

    /* 10. randomized fuzz over stack geometry, materials and angles */
    int fuzz_range = 1, fuzz_energy = 1, fuzz_bound = 1, fuzz_identity = 1, fuzz_cons = 1;
    for (int trial = 0; trial < 30; trial++) {
        int M = 2 + (int)(randu() * 7.0);
        Medium media[LAYERED_MAX_MEDIA];
        double thick[LAYERED_MAX_MEDIA];
        int absorbing = (int)(randu() * 2.0);
        for (int i = 0; i < M; i++) {
            double base_eta = 0.9 + randu() * 1.3;
            for (int c = 0; c < 3; c++) {
                media[i].eta[c] = base_eta * (1.0 + 0.02 * (randu() - 0.5));
                media[i].sigma[c] = absorbing ? randu() * 1.5 : 0.0;
            }
            thick[i] = (i > 0 && i < M - 1) ? randu() * 1.2 : 0.0;
        }
        Stack st;
        if (layered_from_media(&st, media, thick, M) != LAYERED_OK) { fuzz_range = 0; continue; }
        double cos0 = 0.15 + 0.85 * randu();
        Transport ta;
        TraceResult tt;
        if (layered_transport(&st, cos0, &ta) != LAYERED_OK ||
            layered_trace(&st, cos0, 12, &tt) != LAYERED_OK) { fuzz_range = 0; continue; }
        for (int c = 0; c < 3; c++) {
            if (!(ta.T[c] >= -1e-12 && ta.R[c] >= -1e-12 &&
                  ta.T[c] <= 1.0 + 1e-12 && ta.R[c] <= 1.0 + 1e-12 &&
                  ta.T[c] + ta.R[c] <= 1.0 + 1e-10)) fuzz_energy = 0;
            if (fabs(tt.T[c] + tt.R[c] + tt.A[c] + tt.unresolved[c] - 1.0) > 1e-9)
                fuzz_identity = 0;
            double dT = fabs(tt.T[c] - ta.T[c]);
            if (dT > tt.unresolved[c] + 1e-10) fuzz_bound = 0;
            if (ta.T[c] - tt.T[c] < -1e-10 || ta.T[c] - tt.T[c] > tt.unresolved[c] + 1e-10)
                fuzz_bound = 0; /* one-sided: truncation can only lose transmitted paths */
            if (!absorbing && fabs(ta.T[c] + ta.R[c] - 1.0) > 1e-9) fuzz_cons = 0;
        }
    }
    check("fuzz_range", fuzz_range);
    check("fuzz_energy_within_incident", fuzz_energy);
    check("fuzz_trace_identity", fuzz_identity);
    check("fuzz_error_bound", fuzz_bound);
    check("fuzz_lossless_conservation", fuzz_cons);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
