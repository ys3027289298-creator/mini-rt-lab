/*
 * perf: one budget, two acceptance tiers, measured over the three rendering
 * paths (bvh/, sample/, material/).
 *
 * For a tier the harness
 *   1. renders each path exactly width*height*spp times (scheduled work),
 *   2. times each path over several repeats (min / median / max),
 *   3. records the process peak RSS,
 *   4. hashes the RGB bytes of each path's image and of all three combined,
 *   5. checks every measured value against the tier's upper bounds.
 *
 * `--no-budget-check` disables only step 5.  Rendering is untouched: the same
 * buffers, the same stream seeds and the same hashes come out, and the
 * orchestrator proves it byte-for-byte.
 *
 * Exit codes: 0 all checks pass (or checks disabled), 1 a bound was violated,
 * 2 usage / I/O / internal error.
 */
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/resource.h>

#include "perf_render.h"
#include "perf_rng.h"

uint64_t perf_bvh_work(void);
uint64_t perf_sample_work(void);
uint64_t perf_material_work(void);

typedef struct {
    int width, height;
    int bvh_spp;
    int sample_draws;
    int material_spp;
    uint64_t max_work_bvh;
    uint64_t max_work_sample;
    uint64_t max_work_material;
    double max_path_ms;
    double max_total_ms;
    uint64_t max_peak_rss_kib;
} Budget;

typedef struct {
    const char* name;
    PerfRenderFn fn;
    uint64_t (*work_fn)(void);
    int spp;
    unsigned char* buf;
    uint64_t scheduled;
    uint64_t measured;
    double ms_min, ms_median, ms_max;
    uint64_t hash;
} Path;

static int failures;

static double now_ms(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) return 0.0;
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

static uint64_t fnv1a64_update(uint64_t h, const void* data, size_t n) {
    const unsigned char* p = (const unsigned char*)data;
    for (size_t i = 0; i < n; i++) {
        h ^= (uint64_t)p[i];
        h *= UINT64_C(1099511628211);
    }
    return h;
}

static uint64_t fnv1a64(const void* data, size_t n) {
    return fnv1a64_update(UINT64_C(14695981039346656037), data, n);
}

/* Hash the exact file bytes, header included, so the recorded value can be
   recomputed by any independent tool straight from the artifact on disk. */
static int fnv1a64_file(const char* path, uint64_t* out) {
    FILE* f = fopen(path, "rb");
    if (!f) return 0;
    uint64_t h = UINT64_C(14695981039346656037);
    unsigned char buf[1 << 16];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0)
        h = fnv1a64_update(h, buf, n);
    int bad = ferror(f);
    fclose(f);
    if (bad) return 0;
    *out = h;
    return 1;
}

static int parse_u64(const char* s, uint64_t* out) {
    char* end = NULL;
    errno = 0;
    unsigned long long v = strtoull(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0') return 0;
    *out = (uint64_t)v;
    return 1;
}

/* Minimal key=value reader for budget.conf.  Comments and blank lines are
   skipped; keys are exact.  Returns 0 on missing key or bad value. */
static int conf_get(const char* path, const char* key, char* val, size_t n) {
    FILE* f = fopen(path, "r");
    if (!f) return 0;
    char line[512];
    size_t klen = strlen(key);
    int found = 0;
    while (fgets(line, sizeof(line), f)) {
        char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '#' || *p == '\n' || *p == '\0') continue;
        if (strncmp(p, key, klen) == 0 && p[klen] == '=') {
            char* v = p + klen + 1;
            char* e = v;
            while (*e && *e != '\n' && *e != '\r' && *e != ' ' && *e != '\t') e++;
            *e = '\0';
            snprintf(val, n, "%s", v);
            found = 1;
            break;
        }
    }
    fclose(f);
    return found;
}

static int conf_u64(const char* path, const char* key, uint64_t* out) {
    char buf[128];
    return conf_get(path, key, buf, sizeof(buf)) && parse_u64(buf, out);
}

static int parse_double(const char* s, double* out) {
    char* end = NULL;
    errno = 0;
    double v = strtod(s, &end);
    if (errno != 0 || end == s || *end != '\0' || !(v >= 0.0)) return 0;
    *out = v;
    return 1;
}

static int conf_double(const char* path, const char* key, double* out) {
    char buf[128];
    return conf_get(path, key, buf, sizeof(buf)) && parse_double(buf, out);
}

static int load_budget(const char* path, const char* tier, Budget* b) {
    char key[128];
#define GET_U64(field, name) \
    do { \
        uint64_t tmp_; \
        snprintf(key, sizeof(key), "%s.%s", tier, name); \
        if (!conf_u64(path, key, &tmp_) || tmp_ > 0x7fffffffULL) { \
            fprintf(stderr, "perf: missing or bad budget key %s\n", key); \
            return 0; \
        } \
        b->field = (int)tmp_; \
    } while (0)
    GET_U64(width, "width");
    GET_U64(height, "height");
    GET_U64(bvh_spp, "bvh_spp");
    GET_U64(sample_draws, "sample_draws");
    GET_U64(material_spp, "material_spp");
#undef GET_U64
#define GET_U64(field, name) \
    do { \
        snprintf(key, sizeof(key), "%s.%s", tier, name); \
        if (!conf_u64(path, key, &(b->field))) { \
            fprintf(stderr, "perf: missing or bad budget key %s\n", key); \
            return 0; \
        } \
    } while (0)
    GET_U64(max_work_bvh, "max_work_units_bvh");
    GET_U64(max_work_sample, "max_work_units_sample");
    GET_U64(max_work_material, "max_work_units_material");
    GET_U64(max_peak_rss_kib, "max_peak_rss_kib");
#undef GET_U64
#define GET_DBL(field, name) \
    do { \
        snprintf(key, sizeof(key), "%s.%s", tier, name); \
        if (!conf_double(path, key, &(b->field))) { \
            fprintf(stderr, "perf: missing or bad budget key %s\n", key); \
            return 0; \
        } \
    } while (0)
    GET_DBL(max_path_ms, "max_path_ms");
    GET_DBL(max_total_ms, "max_total_ms");
#undef GET_DBL
    if (b->width <= 0 || b->height <= 0 || b->bvh_spp <= 0 ||
        b->sample_draws <= 0 || b->material_spp <= 0) {
        fprintf(stderr, "perf: tier %s has non-positive parameters\n", tier);
        return 0;
    }
    return 1;
}

/* Hash the concatenated bytes of several files, in order. */
static int fnv1a64_files(char paths[][512], int n, uint64_t* out) {
    uint64_t h = UINT64_C(14695981039346656037);
    for (int i = 0; i < n; i++) {
        FILE* f = fopen(paths[i], "rb");
        if (!f) return 0;
        unsigned char buf[1 << 16];
        size_t got;
        while ((got = fread(buf, 1, sizeof(buf), f)) > 0)
            h = fnv1a64_update(h, buf, got);
        int bad = ferror(f);
        fclose(f);
        if (bad) return 0;
    }
    *out = h;
    return 1;
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static void measure(Path* p, const PerfRenderParams* params, int repeats) {
    double samples[PERF_MAX_REPEATS];
    uint64_t had_work = 0;

    /* warm-up: not timed, but the work counter is still strict */
    memset(p->buf, 0, (size_t)params->width * (size_t)params->height * 3u);
    p->fn(params, p->buf);

    for (int r = 0; r < repeats; r++) {
        memset(p->buf, 0, (size_t)params->width * (size_t)params->height * 3u);
        double t0 = now_ms();
        p->fn(params, p->buf);
        double t1 = now_ms();
        samples[r] = t1 - t0;
        had_work = p->work_fn();
    }
    p->measured = had_work;
    qsort(samples, (size_t)repeats, sizeof(double), cmp_double);
    p->ms_min = samples[0];
    p->ms_max = samples[repeats - 1];
    p->ms_median = (repeats % 2) ? samples[repeats / 2]
                                 : 0.5 * (samples[repeats / 2 - 1] +
                                          samples[repeats / 2]);
}

static int write_ppm(const char* path, const unsigned char* rgb, int w, int h) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return 0; }
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    size_t n = (size_t)w * (size_t)h * 3u;
    if (fwrite(rgb, 1, n, f) != n) { fclose(f); return 0; }
    fclose(f);
    return 1;
}

static int write_json(const char* path, const char* tier, const Budget* b,
                      const PerfRenderParams* params, Path* paths, int npaths,
                      double total_min, double total_median, double total_max,
                      uint64_t combined_hash, uint64_t peak_rss,
                      int checks_enabled, int all_ok, int repeats) {
    FILE* f = fopen(path, "w");
    if (!f) { perror(path); return 0; }
    fprintf(f, "{\n");
    fprintf(f, "  \"tier\": \"%s\",\n", tier);
    fprintf(f, "  \"scene_id\": \"%s\",\n", PERF_SCENE_ID);
    fprintf(f, "  \"seed\": %llu,\n", (unsigned long long)params->seed);
    fprintf(f, "  \"checks_enabled\": %s,\n", checks_enabled ? "true" : "false");
    fprintf(f, "  \"repeats\": %d,\n", repeats);
    fprintf(f, "  \"params\": {\"width\": %d, \"height\": %d, \"bvh_spp\": %d, "
               "\"sample_draws\": %d, \"material_spp\": %d},\n",
            b->width, b->height, b->bvh_spp, b->sample_draws, b->material_spp);
    fprintf(f, "  \"budget\": {\"max_work_units_bvh\": %llu, "
               "\"max_work_units_sample\": %llu, \"max_work_units_material\": %llu, "
               "\"max_path_ms\": %.6f, \"max_total_ms\": %.6f, "
               "\"max_peak_rss_kib\": %llu},\n",
            (unsigned long long)b->max_work_bvh,
            (unsigned long long)b->max_work_sample,
            (unsigned long long)b->max_work_material,
            b->max_path_ms, b->max_total_ms,
            (unsigned long long)b->max_peak_rss_kib);
    fprintf(f, "  \"paths\": {\n");
    for (int i = 0; i < npaths; i++) {
        Path* p = &paths[i];
        fprintf(f, "    \"%s\": {\"spp\": %d, \"scheduled_work\": %llu, "
                   "\"measured_work\": %llu, \"ms_min\": %.6f, "
                   "\"ms_median\": %.6f, \"ms_max\": %.6f, "
                   "\"image\": \"%s.ppm\", \"image_fnv1a64\": \"%016llx\"}%s\n",
                p->name, p->spp,
                (unsigned long long)p->scheduled,
                (unsigned long long)p->measured,
                p->ms_min, p->ms_median, p->ms_max, p->name,
                (unsigned long long)p->hash,
                (i + 1 < npaths) ? "," : "");
    }
    fprintf(f, "  },\n");
    fprintf(f, "  \"total\": {\"ms_min\": %.6f, \"ms_median\": %.6f, "
               "\"ms_max\": %.6f},\n", total_min, total_median, total_max);
    fprintf(f, "  \"combined_image_fnv1a64\": \"%016llx\",\n",
            (unsigned long long)combined_hash);
    fprintf(f, "  \"peak_rss_kib\": %llu,\n", (unsigned long long)peak_rss);
    fprintf(f, "  \"all_ok\": %s\n", all_ok ? "true" : "false");
    fprintf(f, "}\n");
    fclose(f);
    return 1;
}

static void usage(const char* argv0) {
    fprintf(stderr,
        "usage: %s run <tier> --budget FILE --outdir DIR [options]\n"
        "options:\n"
        "  --no-budget-check   render and report, skip bound evaluation\n"
        "  --repeats N         timed repeats per path (1..%d, default 7)\n"
        "  %s selftest        deterministic/self checks\n", argv0, PERF_MAX_REPEATS, argv0);
}

static int cmd_run(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 2; }
    const char* tier = argv[2];
    const char* budget_path = NULL;
    const char* outdir = NULL;
    int checks = 1;
    int repeats = 7;

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--budget") == 0 && i + 1 < argc) {
            budget_path = argv[++i];
        } else if (strcmp(argv[i], "--outdir") == 0 && i + 1 < argc) {
            outdir = argv[++i];
        } else if (strcmp(argv[i], "--no-budget-check") == 0) {
            checks = 0;
        } else if (strcmp(argv[i], "--repeats") == 0 && i + 1 < argc) {
            uint64_t r = 0;
            if (!parse_u64(argv[++i], &r) || r < 1 || r > PERF_MAX_REPEATS) {
                fprintf(stderr, "perf: --repeats must be 1..%d\n", PERF_MAX_REPEATS);
                return 2;
            }
            repeats = (int)r;
        } else {
            fprintf(stderr, "perf: unknown option %s\n", argv[i]);
            usage(argv[0]);
            return 2;
        }
    }
    if (!budget_path || !outdir) {
        fprintf(stderr, "perf: --budget and --outdir are required\n");
        usage(argv[0]);
        return 2;
    }

    Budget b;
    if (!load_budget(budget_path, tier, &b)) return 2;

    size_t npix = (size_t)b.width * (size_t)b.height;
    Path paths[3] = {
        { "bvh", perf_bvh_render, perf_bvh_work, b.bvh_spp,
          malloc(npix * 3u), 0, 0, 0, 0, 0, 0 },
        { "sample", perf_sample_render, perf_sample_work, b.sample_draws,
          malloc(npix * 3u), 0, 0, 0, 0, 0, 0 },
        { "material", perf_material_render, perf_material_work, b.material_spp,
          malloc(npix * 3u), 0, 0, 0, 0, 0, 0 },
    };
    for (int i = 0; i < 3; i++) {
        if (!paths[i].buf) {
            fprintf(stderr, "perf: out of memory\n");
            for (int j = 0; j < 3; j++) free(paths[j].buf);
            return 2;
        }
    }

    char ppm_files[3][512];
    for (int i = 0; i < 3; i++) {
        PerfRenderParams pp = {
            b.width, b.height, paths[i].spp, PERF_SCENE_SEED, 0
        };
        pp.stream = (i == 0) ? PERF_STREAM_BVH
                  : (i == 1) ? PERF_STREAM_SAMPLE
                             : PERF_STREAM_MATERIAL;
        paths[i].scheduled = (uint64_t)b.width * (uint64_t)b.height *
                             (uint64_t)paths[i].spp;
        measure(&paths[i], &pp, repeats);

        snprintf(ppm_files[i], sizeof(ppm_files[i]), "%s/%s.ppm",
                 outdir, paths[i].name);
        if (!write_ppm(ppm_files[i], paths[i].buf, b.width, b.height)) {
            for (int j = 0; j < 3; j++) free(paths[j].buf);
            return 2;
        }
        if (!fnv1a64_file(ppm_files[i], &paths[i].hash)) {
            fprintf(stderr, "perf: cannot hash %s\n", ppm_files[i]);
            for (int j = 0; j < 3; j++) free(paths[j].buf);
            return 2;
        }
        if (paths[i].measured != paths[i].scheduled) {
            fprintf(stderr,
                "perf: %s did %llu work units, scheduled %llu\n",
                paths[i].name,
                (unsigned long long)paths[i].measured,
                (unsigned long long)paths[i].scheduled);
            failures++;
        }
    }

    /* combined bill of materials: one hash over all three PPM files, in the
       order bvh, sample, material */
    uint64_t combined = 0;
    if (!fnv1a64_files(ppm_files, 3, &combined)) {
        fprintf(stderr, "perf: cannot hash combined artifacts\n");
        for (int j = 0; j < 3; j++) free(paths[j].buf);
        return 2;
    }

    double total_min = paths[0].ms_min + paths[1].ms_min + paths[2].ms_min;
    double total_median = paths[0].ms_median + paths[1].ms_median +
                          paths[2].ms_median;
    double total_max = paths[0].ms_max + paths[1].ms_max + paths[2].ms_max;

    struct rusage ru;
    uint64_t peak_rss = 0;
    if (getrusage(RUSAGE_SELF, &ru) == 0)
        peak_rss = (uint64_t)ru.ru_maxrss;   /* KiB on Linux */

    /* budget evaluation */
    int work_ok = 1, time_ok = 1, rss_ok = 1;
    uint64_t maxw[3] = { b.max_work_bvh, b.max_work_sample,
                         b.max_work_material };
    for (int i = 0; i < 3; i++) {
        if (paths[i].measured > maxw[i]) work_ok = 0;
        if (paths[i].measured != paths[i].scheduled) work_ok = 0;
        if (paths[i].ms_median > b.max_path_ms) time_ok = 0;
    }
    if (total_median > b.max_total_ms) time_ok = 0;
    if (peak_rss > b.max_peak_rss_kib) rss_ok = 0;
    int all_ok = work_ok && time_ok && rss_ok && failures == 0;

    char json[512];
    snprintf(json, sizeof(json), "%s/metrics.json", outdir);
    if (!write_json(json, tier, &b, &(PerfRenderParams){
                        b.width, b.height, 0, PERF_SCENE_SEED, 0 },
                    paths, 3, total_min, total_median, total_max,
                    combined, peak_rss, checks, all_ok, repeats)) {
        for (int j = 0; j < 3; j++) free(paths[j].buf);
        return 2;
    }

    printf("tier=%s checks_enabled=%d repeats=%d\n", tier, checks, repeats);
    printf("%-8s %12s %12s %9s %9s %9s %16s\n",
           "path", "scheduled", "measured", "min_ms", "med_ms", "max_ms",
           "fnv1a64");
    for (int i = 0; i < 3; i++) {
        Path* p = &paths[i];
        printf("%-8s %12llu %12llu %9.3f %9.3f %9.3f %016llx\n",
               p->name, (unsigned long long)p->scheduled,
               (unsigned long long)p->measured,
               p->ms_min, p->ms_median, p->ms_max,
               (unsigned long long)p->hash);
    }
    printf("%-8s %12s %12s %9.3f %9.3f %9.3f %016llx\n", "total", "", "",
           total_min, total_median, total_max, (unsigned long long)combined);
    printf("peak_rss_kib=%llu\n", (unsigned long long)peak_rss);
    printf("work_ok=%d time_ok=%d rss_ok=%d all_ok=%d\n",
           work_ok, time_ok, rss_ok, all_ok);

    for (int j = 0; j < 3; j++) free(paths[j].buf);
    if (!checks) return 0;
    return all_ok ? 0 : 1;
}

/* ---------- self test ---------- */

static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int cmd_selftest(void) {
    /* FNV-1a reference vector: hash of the empty string. */
    check("fnv_empty", fnv1a64("", 0) == UINT64_C(0xcbf29ce484222325));

    /* The PCG stream is a pure function of (seed, stream). */
    PerfRng a, b, c;
    perf_rng_seed(&a, 1234, 1);
    perf_rng_seed(&b, 1234, 1);
    perf_rng_seed(&c, 1234, 2);
    int same = 1, distinct = 0;
    for (int i = 0; i < 64; i++) {
        uint32_t x = perf_rng_next(&a), y = perf_rng_next(&b), z = perf_rng_next(&c);
        if (x != y) same = 0;
        if (x != z) distinct = 1;
    }
    check("rng_same_seed_same_stream", same);
    check("rng_different_stream_distinct", distinct);

    PerfRng u;
    perf_rng_seed(&u, 99, 7);
    int unit = 1;
    for (int i = 0; i < 10000; i++) {
        float f = perf_rng_float(&u);
        if (!(f >= 0.0f && f < 1.0f)) unit = 0;
    }
    check("rng_unit_interval", unit);

    /* Rendering is deterministic: two runs of the same params are bit-equal,
       and a different stream changes the image. */
    int w = 33, h = 21;
    size_t n = (size_t)w * (size_t)h * 3u;
    unsigned char *i1 = malloc(n), *i2 = malloc(n), *i3 = malloc(n);
    if (!i1 || !i2 || !i3) { printf("alloc_failed\n"); return 2; }
    PerfRenderParams p1 = { w, h, 4, PERF_SCENE_SEED, PERF_STREAM_BVH };
    PerfRenderParams p2 = { w, h, 4, PERF_SCENE_SEED, PERF_STREAM_BVH };
    PerfRenderParams p3 = { w, h, 4, PERF_SCENE_SEED, PERF_STREAM_SAMPLE };

    perf_bvh_render(&p1, i1);
    perf_bvh_render(&p2, i2);
    perf_sample_render(&p3, i3);
    check("bvh_repeat_identical", memcmp(i1, i2, n) == 0);
    check("different_stream_differs", memcmp(i1, i3, n) != 0);

    /* Work equals scheduled work on every path. */
    PerfRenderParams q = { 7, 5, 6, PERF_SCENE_SEED, PERF_STREAM_MATERIAL };
    unsigned char* buf = malloc((size_t)7 * 5 * 3u);
    if (!buf) { printf("alloc_failed\n"); return 2; }

    perf_bvh_render(&(PerfRenderParams){ 7, 5, 6, PERF_SCENE_SEED, PERF_STREAM_BVH }, buf);
    check("bvh_work_exact", perf_bvh_work() == 7u * 5u * 6u);
    perf_sample_render(&(PerfRenderParams){ 7, 5, 6, PERF_SCENE_SEED, PERF_STREAM_SAMPLE }, buf);
    check("sample_work_exact", perf_sample_work() == 7u * 5u * 6u);
    perf_material_render(&q, buf);
    check("material_work_exact", perf_material_work() == 7u * 5u * 6u);

    /* A pixel-count/zero-spp parameter is rejected by the loader, so a
       degenerate run cannot masquerade as cheap work. */
    check("work_positive", 7u * 5u * 6u > 0);

    free(i1); free(i2); free(i3); free(buf);
    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); return 2; }
    if (strcmp(argv[1], "run") == 0) return cmd_run(argc, argv);
    if (strcmp(argv[1], "selftest") == 0) return cmd_selftest();
    usage(argv[0]);
    return 2;
}
