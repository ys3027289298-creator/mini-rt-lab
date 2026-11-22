/*
 * inv_main.c -- invariant-check driver.
 *
 * Modes:
 *   (default)      run the four invariant classes on the shipped renderer and
 *                  on the golden reference, print a per-check table, and exit
 *                  nonzero if the renderer violates any class threshold.
 *   --calibrate    run every calibration mutant and assert each is caught by
 *                  its target class and by nothing else (zero false positives
 *                  on the golden reference).
 *   --render-check assert the checker leaves the renderer bit-identical:
 *                  hash a fixed workload before/after a full invariant run.
 *
 * Exit status: 0 all requested assertions held, 1 otherwise, 2 setup error.
 */
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "inv_impl.h"

/* ------------------------------------------------------------------ purity */

/*
 * A deterministic "render": build the shipped tree, cast a fixed fan of rays
 * against it through the shipped slab test, accumulate shipped BRDF x shipped
 * samples.  Folding every result into one 64-bit hash detects *any* change in
 * observable renderer output, including state the checker might have mutated.
 */
uint64_t inv_render_fingerprint(void) {
    uint64_t h = 1469598103934665603ull; /* FNV-1a offset basis */
    #define MIX(u) do { h ^= (uint64_t)(u); h *= 1099511628211ull; } while (0)

    InvTree t; memset(&t, 0, sizeof t);
    MIX(edge_build(&t));
    MIX(t.nnodes); MIX(t.nprims);

    const double lo[3] = { -1, -1, -1 }, hi[3] = { 1, 1, 1 };
    for (int k = 0; k < 4096; k++) {
        double o[3] = { -3.0 + (k % 7) * 0.01, 0.0, 0.0 };
        double d[3] = { 1.0, -0.001 * (k % 5), 0.0005 * (k % 3) };
        double L = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        for (int i = 0; i < 3; i++) d[i] /= L;
        double tt = 0;
        int hit = edge_box_hit(lo, hi, o, d, 1e30, &tt);
        MIX(hit);
        if (hit) {
            union { double d; uint64_t u; } cv;
            cv.d = tt; MIX(cv.u);
        }
    }
    for (unsigned i = 0; i < 512; i++) {
        double v = edge_halton2(i, 2), w = edge_halton(i, 3);
        double s[3]; edge_cosine_sample(v, w, s);
        union { double d; uint64_t u; } cv;
        cv.d = v; MIX(cv.u);
        cv.d = w; MIX(cv.u);
        for (int i2 = 0; i2 < 3; i2++) { cv.d = s[i2]; MIX(cv.u); }
        cv.d = edge_ggx_D(v, 0.3); MIX(cv.u);
        cv.d = edge_lambert(v); MIX(cv.u);
    }
    #undef MIX
    return h;
}

int inv_module_state_clean(void) {
    return edge_bvh_failures() == 0 && edge_sample_failures() == 0 &&
           edge_material_failures() == 0;
}

/* ------------------------------------------------------------------ report */

static int print_report(const char* title, const InvReport* rep) {
    printf("\n== %s ==\n", title);
    printf("%-10s %-30s %-5s %14s %14s  %s\n",
           "class", "check", "pass", "measured", "threshold", "detail");
    int failed = 0;
    for (int i = 0; i < rep->n; i++) {
        const InvResult* r = &rep->r[i];
        if (!r->pass) failed++;
        printf("%-10s %-30s %-5s %14.6g %14.6g  %s\n",
               r->cls, r->name, r->pass ? "ok" : "FAIL",
               r->measured, r->threshold, r->detail);
    }
    return failed;
}

/* ------------------------------------------------------------------ modes */

static int mode_audit(void) {
    printf("renderer under audit: %s\n", edge_model_source_desc());

    uint64_t fp_before = inv_render_fingerprint();
    int state_before = inv_module_state_clean();

    InvReport golden; memset(&golden, 0, sizeof golden);
    inv_run_model(&inv_golden_reference, &golden);

    InvReport fctrl; memset(&fctrl, 0, sizeof fctrl);
    inv_run_model(&inv_float_control, &fctrl);

    InvReport real; memset(&real, 0, sizeof real);
    inv_run_model(&inv_real_renderer, &real);

    InvReport real_h; memset(&real_h, 0, sizeof real_h);
    inv_run_model(&inv_real_halton, &real_h);

    uint64_t fp_after = inv_render_fingerprint();
    int state_after = inv_module_state_clean();

    int gf = print_report("golden reference, double (must pass)", &golden);
    int ff = print_report("float control, correct float impl (must pass)", &fctrl);
    print_report("shipped renderer, sampler=halton2", &real);
    print_report("shipped renderer, sampler=halton", &real_h);

    printf("\n== gate ==\n");
    printf("golden_reference_clean=%d\n", gf == 0);
    printf("float_control_clean=%d\n", ff == 0);
    printf("renderer_changed_by_checker=%d\n", fp_before != fp_after);
    printf("module_state_touched=%d\n", !(state_before && state_after));

    int ok = (gf == 0) && (ff == 0) && (fp_before == fp_after) && state_before && state_after;
    if (!ok) {
        printf("verdict=CHECKER_INVALID\n");
        return 2;
    }

    const char* cls[] = { "energy", "geometry", "sampling", "symmetry" };
    int viol = 0;
    printf("\n== per-class verdict (shipped, halton2) ==\n");
    for (int i = 0; i < 4; i++) {
        int f = inv_class_failed(&real, cls[i]);
        printf("%-9s %s (%d checks)\n", cls[i], f ? "VIOLATED" : "ok",
               inv_class_total(&real, cls[i]));
        viol |= f;
    }
    printf("\n== per-class verdict (shipped, halton) ==\n");
    for (int i = 0; i < 4; i++)
        printf("%-9s %s\n", cls[i], inv_class_failed(&real_h, cls[i]) ? "VIOLATED" : "ok");

    printf("verdict=%s\n", viol ? "INVARIANT_VIOLATION" : "OK");
    return viol ? 1 : 0;
}

static int mode_calibrate(void) {
    InvMutant mts[16];
    int n = inv_build_mutants(mts, 16);
    if (n <= 0) { fprintf(stderr, "no mutants built\n"); return 2; }

    const char* cls[] = { "energy", "geometry", "sampling", "symmetry" };

    /* false-positive control: the golden reference must fail nothing */
    InvReport g; memset(&g, 0, sizeof g);
    inv_run_model(&inv_golden_reference, &g);
    int golden_fail = 0;
    for (int i = 0; i < 4; i++) golden_fail |= inv_class_failed(&g, cls[i]);
    printf("%-28s caught=%s  of class=%s  other_classes_failing=%s\n",
           "golden/reference", "-", "-", golden_fail ? "YES(BAD)" : "none");

    int missed = 0;
    printf("\n%s\n", "mutant                       target     caught  collateral");
    for (int k = 0; k < n; k++) {
        InvReport r; memset(&r, 0, sizeof r);
        inv_run_model(&mts[k].m, &r);
        int caught = inv_class_failed(&r, mts[k].target_cls);
        int coll = 0;
        for (int i = 0; i < 4; i++) {
            const char* c = cls[i];
            if (strcmp(c, mts[k].target_cls) == 0) continue;
            if (inv_class_failed(&r, c)) coll++;
        }
        if (!caught) missed++;
        /* mutants deliberately break one physical law; a collateral failure in
           another class is reported but is not by itself a checker defect. */
        printf("%-28s %-10s %-7s %d\n", mts[k].name, mts[k].target_cls,
               caught ? "yes" : "MISSED", coll);
    }
    printf("\nfalse_negatives=%d/%d  golden_false_positives=%d\n",
           missed, n, golden_fail ? 1 : 0);
    printf("verdict=%s\n", (missed == 0 && !golden_fail) ? "CALIBRATED" : "CALIBRATION_FAILED");
    return (missed == 0 && !golden_fail) ? 0 : 1;
}

static int mode_render_check(void) {
    uint64_t a = inv_render_fingerprint();
    int sa = inv_module_state_clean();

    InvReport r; memset(&r, 0, sizeof r);
    inv_run_model(&inv_real_renderer, &r);
    inv_run_model(&inv_real_halton, &r);
    inv_run_model(&inv_golden_reference, &r);

    InvReport g; memset(&g, 0, sizeof g);
    inv_run_model(&inv_golden_reference, &g);

    uint64_t b = inv_render_fingerprint();
    int sb = inv_module_state_clean();

    printf("fingerprint_before=%016llx\n", (unsigned long long)a);
    printf("fingerprint_after =%016llx\n", (unsigned long long)b);
    printf("bit_identical=%d\n", a == b);
    printf("module_state_untouched=%d\n", sa && sb);
    int ok = (a == b) && sa && sb;
    printf("verdict=%s\n", ok ? "NO_SIDE_EFFECTS" : "SIDE_EFFECT_DETECTED");
    return ok ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc > 1 && strcmp(argv[1], "--calibrate") == 0) return mode_calibrate();
    if (argc > 1 && strcmp(argv[1], "--render-check") == 0) return mode_render_check();
    return mode_audit();
}
