
/* budget: fixed total sample budget B = pixels * spp trade-off.
   Error model E(P, s) = a / P + v / s, where a/P is the discretization
   bias of a P-pixel image and v/s the Monte Carlo variance with s
   samples per pixel. Enumerates every feasible (resolution, spp) combo
   on discrete grids, finds the exhaustive discrete optimum and a pick
   derived from the continuous relaxation, and certifies an upper bound
   on the pick's deviation from the exhaustive optimum. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define MAX_GRID 64
#define MAX_COMBO 256
#define EPS 1e-12

typedef struct { int w, h; long long px; } Res;
typedef struct {
    int w, h, spp;
    long long px, used;
    double err;
} Combo;
typedef struct {
    Combo c[MAX_COMBO];
    int n;
    int exh_i, pick_i;   /* indices into c[], -1 when infeasible */
    double pc, sc, ec;   /* continuous relaxation optimum */
    int feasible;
} Analysis;

static double combo_err(double a, double v, const Combo* q) {
    return a / (double)q->px + v / (double)q->spp;
}

/* Returns the number of feasible combos; 0 means no feasible solution. */
static int analyze(long long B, double a, double v,
                   const Res* res, int nres, const int* spp, int nspp,
                   Analysis* z) {
    memset(z, 0, sizeof(*z));
    z->exh_i = z->pick_i = -1;

    for (int i = 0; i < nres; i++) {
        for (int j = 0; j < nspp; j++) {
            long long used = res[i].px * (long long)spp[j];
            if (used > B) continue;                 /* over budget */
            if (z->n >= MAX_COMBO) return -1;
            Combo* q = &z->c[z->n++];
            q->w = res[i].w; q->h = res[i].h; q->px = res[i].px;
            q->spp = spp[j]; q->used = used;
            q->err = combo_err(a, v, q);
        }
    }
    z->feasible = z->n > 0;

    /* continuous relaxation: any P in (0, B], s = B/P */
    z->pc = sqrt(a * (double)B / v);
    z->sc = (double)B / z->pc;
    z->ec = 2.0 * sqrt(a * v / (double)B);
    if (!z->feasible) return 0;

    /* exhaustive discrete optimum (first minimum wins, deterministic) */
    z->exh_i = 0;
    for (int i = 1; i < z->n; i++)
        if (z->c[i].err < z->c[z->exh_i].err) z->exh_i = i;

    /* fast pick: feasible resolution nearest to P* on the log scale,
       ties resolved to the smaller P; then largest feasible spp. */
    long long pick_px = -1;
    double best_d = 0.0;
    for (int i = 0; i < nres; i++) {
        int any = 0;
        for (int j = 0; j < nspp; j++) {
            if (res[i].px * (long long)spp[j] <= B) { any = 1; break; }
        }
        if (!any) continue;
        double d = fabs(log((double)res[i].px / z->pc));
        if (pick_px < 0 || d < best_d - 1e-15 ||
            (fabs(d - best_d) <= 1e-15 && res[i].px < pick_px)) {
            best_d = d;
            pick_px = res[i].px;
        }
    }
    /* combos are stored P-ascending then spp-ascending, so the last
       combo with this P carries the largest feasible spp */
    for (int i = 0; i < z->n; i++)
        if (z->c[i].px == pick_px) z->pick_i = i;
    return z->n;
}

/* Structural a-priori guarantee applies only when the picked P is within
   a factor of 2 of the continuous optimum and the picked spp wastes at
   most half of the per-pixel budget. */
static int structural_applies(long long B, const Analysis* z) {
    if (!z->feasible) return 0;
    const Combo* p = &z->c[z->pick_i];
    double r = (double)p->px / z->pc;
    double s_cont = (double)B / (double)p->px;
    return r >= 0.5 - EPS && r <= 2.0 + EPS &&
           (double)p->spp >= 0.5 * s_cont - EPS;
}

static void print_combo(const char* tag, const Combo* q) {
    printf("%s: w=%d h=%d pixels=%lld spp=%d used=%lld err=%.9g\n",
           tag, q->w, q->h, q->px, q->spp, q->used, q->err);
}

static void print_report(long long B, double a, double v,
                         const Res* res, int nres, const int* spp, int nspp) {
    Analysis z;
    int n = analyze(B, a, v, res, nres, spp, nspp, &z);
    long long min_cost = res[0].px * (long long)spp[0];

    printf("budget B=%lld a=%g v=%g\n", B, a, v);
    printf("res_grid=%d spp_grid=%d min_cost=%lld feasible=%d\n",
           nres, nspp, min_cost, n);

    if (n == 0) {
        printf("no feasible solution: min cost %lld exceeds budget %lld\n",
               min_cost, B);
        return;
    }

    printf("curve:\n");
    printf("w h pixels spp used err\n");
    for (int i = 0; i < n; i++) {
        const Combo* q = &z.c[i];
        printf("%d %d %lld %d %lld %.9g\n",
               q->w, q->h, q->px, q->spp, q->used, q->err);
    }

    printf("continuous_opt: P*=%.6g s*=%.6g E*=%.9g\n", z.pc, z.sc, z.ec);
    print_combo("exhaustive_opt", &z.c[z.exh_i]);
    print_combo("fast_pick", &z.c[z.pick_i]);

    double exh = z.c[z.exh_i].err, pick = z.c[z.pick_i].err;
    double gap = pick - exh;            /* actual deviation */
    double bound = pick - z.ec;         /* certified instance upper bound */
    printf("gap=%.9g bound=%.9g bound_holds=%d\n",
           gap, bound, gap <= bound + EPS * fmax(1.0, fabs(bound)));
    printf("rel_gap=%.6g rel_bound=%.6g\n",
           gap / exh, bound / z.ec);

    if (structural_applies(B, &z)) {
        double r = (double)z.c[z.pick_i].px / z.pc;
        printf("structural_bound: applicable=1 r=%.6g "
               "pick/E*<=2.25 rel_gap<=1.25 (observed %.6g, %.6g)\n",
               r, pick / z.ec, gap / exh);
    } else {
        printf("structural_bound: applicable=0 (instance bound still valid)\n");
    }
}

static void print_csv(long long B, double a, double v,
                      const Res* res, int nres, const int* spp, int nspp) {
    Analysis z;
    int n = analyze(B, a, v, res, nres, spp, nspp, &z);
    printf("# feasible=%d B=%lld a=%g v=%g pc=%.10g sc=%.10g ec=%.10g",
           n > 0, B, a, v, z.pc, z.sc, z.ec);
    if (n > 0)
        printf(" exh_err=%.10g pick_err=%.10g structural=%d",
               z.c[z.exh_i].err, z.c[z.pick_i].err,
               structural_applies(B, &z));
    printf("\n");
    if (n == 0) return;
    printf("w,h,pixels,spp,used,err,role\n");
    for (int i = 0; i < n; i++) {
        const Combo* q = &z.c[i];
        const char* role = (i == z.exh_i && i == z.pick_i) ? "both"
                         : (i == z.exh_i) ? "exhaustive"
                         : (i == z.pick_i) ? "pick" : "-";
        printf("%d,%d,%lld,%d,%lld,%.10g,%s\n",
               q->w, q->h, q->px, q->spp, q->used, q->err, role);
    }
}

/* ---------- self tests ---------- */

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* hand-computed 2x2 grid: P in {4,16}, spp in {1,2}, B=16, a=v=1 */
static void test_hand_grid(void) {
    Res res[2] = {{2, 2, 4}, {4, 4, 16}};
    int spp[2] = {1, 2};
    Analysis z;
    int n = analyze(16, 1.0, 1.0, res, 2, spp, 2, &z);

    check("hand_feasible_count", n == 3);          /* (16,2) costs 32 */
    check("hand_budget_respected",
          z.c[0].used == 4 && z.c[1].used == 8 && z.c[2].used == 16);
    check("hand_err_values",
          fabs(z.c[0].err - 1.25) < 1e-12 &&
          fabs(z.c[1].err - 0.75) < 1e-12 &&
          fabs(z.c[2].err - 1.0625) < 1e-12);
    check("hand_exhaustive", z.exh_i == 1 &&
          fabs(z.c[z.exh_i].err - 0.75) < 1e-12);
    check("hand_continuous",
          fabs(z.pc - 4.0) < 1e-12 && fabs(z.ec - 0.5) < 1e-12);
    check("hand_pick", z.pick_i == 1);             /* P=4, largest spp 2 */

    double gap = z.c[z.pick_i].err - z.c[z.exh_i].err;
    double bound = z.c[z.pick_i].err - z.ec;
    check("hand_gap_zero", fabs(gap) < 1e-12);
    check("hand_bound_holds", gap <= bound + EPS);
    check("hand_structural", structural_applies(16, &z) &&
          z.c[z.pick_i].err <= 2.25 * z.ec + EPS);
}

static void test_empty(const Res* res, int nres, const int* spp, int nspp) {
    Analysis z;
    long long min_cost = res[0].px * (long long)spp[0];
    long long B = min_cost - 1;                   /* always infeasible */
    int n = analyze(B, 1.0, 1.0, res, nres, spp, nspp, &z);
    check("empty_no_feasible", n == 0 && z.feasible == 0);
    check("empty_sentinels", z.exh_i == -1 && z.pick_i == -1);
    check("empty_continuous_still_set", z.ec > 0.0 && z.pc > 0.0);
}

static void test_grids(long long B, double a, double v,
                       const Res* res, int nres, const int* spp, int nspp) {
    Analysis z;
    int n = analyze(B, a, v, res, nres, spp, nspp, &z);
    if (n == 0) return;

    /* every combo respects the budget and error is monotone in spp
       for a fixed resolution */
    int budget_ok = 1, mono_ok = 1;
    for (int i = 0; i < n; i++) {
        if (z.c[i].used > B) budget_ok = 0;
        for (int j = i + 1; j < n; j++) {
            if (z.c[j].px == z.c[i].px && z.c[j].spp > z.c[i].spp &&
                z.c[j].err >= z.c[i].err) mono_ok = 0;
        }
    }
    check("sweep_budget", budget_ok);
    check("sweep_monotone_spp", mono_ok);

    double exh = z.c[z.exh_i].err, pick = z.c[z.pick_i].err;
    double gap = pick - exh;
    double bound = pick - z.ec;
    double tol = 1e-9 * fmax(1.0, fabs(bound));

    /* exhaustive optimum must beat every combo */
    int exh_ok = 1;
    for (int i = 0; i < n; i++)
        if (z.c[i].err < exh - EPS) exh_ok = 0;
    check("sweep_exhaustive_min", exh_ok);

    /* universal instance bound, no grid-spacing assumption needed:
       E_pick >= E_exh >= E*_c  =>  gap <= E_pick - E*_c = bound */
    check("sweep_pick_above_relaxation", pick >= z.ec - EPS);
    check("sweep_bound_holds", gap >= -EPS && gap <= bound + tol);

    if (structural_applies(B, &z)) {
        check("sweep_structural_225", pick <= 2.25 * z.ec + EPS);
        check("sweep_structural_rel125",
              gap / exh <= 1.25 + 1e-12);
    }
}

/* optimum lands exactly on a grid point: deviation must be zero */
static void test_exact_match(const Res* res, int nres,
                             const int* spp, int nspp) {
    long long B = 1048576;            /* P*=1024 (32x32), s*=1024 */
    Analysis z;
    int n = analyze(B, 1.0, 1.0, res, nres, spp, nspp, &z);
    check("exact_feasible", n > 0);
    check("exact_pick_pixels", z.c[z.pick_i].px == 1024);
    check("exact_pick_spp", z.c[z.pick_i].spp == 1024);
    double gap = z.c[z.pick_i].err - z.c[z.exh_i].err;
    check("exact_gap_zero", fabs(gap) < 1e-12);
    check("exact_matches_relaxation",
          fabs(z.c[z.pick_i].err - z.ec) < 1e-10);
}

static void selftest(void) {
    test_hand_grid();

    int sides[] = {8, 16, 32, 64, 128, 256};
    int nres = (int)(sizeof(sides) / sizeof(sides[0]));
    Res res[MAX_GRID];
    for (int i = 0; i < nres; i++)
        res[i] = (Res){sides[i], sides[i], (long long)sides[i] * sides[i]};
    int spp[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    int nspp = (int)(sizeof(spp) / sizeof(spp[0]));

    test_empty(res, nres, spp, nspp);
    test_exact_match(res, nres, spp, nspp);

    double av[] = {0.25, 1.0, 4.0, 16.0};
    long long Bv[] = {256, 4096, 1LL << 16, 1LL << 22};
    int structural_seen = 0;
    for (size_t i = 0; i < sizeof(av) / sizeof(av[0]); i++) {
        for (size_t j = 0; j < sizeof(Bv) / sizeof(Bv[0]); j++) {
            Analysis z;
            int n = analyze(Bv[j], av[i], 1.0, res, nres, spp, nspp, &z);
            if (n > 0 && structural_applies(Bv[j], &z)) structural_seen = 1;
            test_grids(Bv[j], av[i], 1.0, res, nres, spp, nspp);
        }
    }
    check("structural_case_observed", structural_seen);

    printf("all_ok=%d\n", failures == 0);
}

int main(int argc, char** argv) {
    int sides[] = {8, 16, 32, 64, 128, 256};
    int nres = (int)(sizeof(sides) / sizeof(sides[0]));
    Res res[MAX_GRID];
    for (int i = 0; i < nres; i++)
        res[i] = (Res){sides[i], sides[i], (long long)sides[i] * sides[i]};
    int spp[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024};
    int nspp = (int)(sizeof(spp) / sizeof(spp[0]));

    if (argc == 1 || strcmp(argv[1], "selftest") == 0) {
        selftest();
        return failures == 0 ? 0 : 1;
    }

    if (strcmp(argv[1], "empty") == 0) {
        print_report(32, 1.0, 1.0, res, nres, spp, nspp);
        return 0;
    }

    if (argc >= 2 && (strcmp(argv[1], "report") == 0 ||
                      strcmp(argv[1], "csv") == 0)) {
        long long B = (argc > 2) ? atoll(argv[2]) : 1048576;
        double a = (argc > 3) ? atof(argv[3]) : 1.0;
        double v = (argc > 4) ? atof(argv[4]) : 1.0;
        if (B <= 0 || !(a > 0.0) || !(v > 0.0)) {
            fprintf(stderr, "usage: %s [report|csv] B a v  (all positive)\n",
                    argv[0]);
            return 2;
        }
        if (strcmp(argv[1], "csv") == 0)
            print_csv(B, a, v, res, nres, spp, nspp);
        else
            print_report(B, a, v, res, nres, spp, nspp);
        return 0;
    }

    fprintf(stderr,
            "usage: %s [selftest|report B a v|csv B a v|empty]\n", argv[0]);
    return 2;
}
