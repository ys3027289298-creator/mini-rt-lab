
/* sample: stratified and low-discrepancy sampling with variance estimate. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N_MAX 256

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

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* deterministic render: scatter plot of Halton (2,3) and (5,7) points.
   Only integer and IEEE double arithmetic is used, so the output image is
   bit-identical across compilers, optimization levels and machines. */
static int render(const char* path, int size, int strict) {
    size_t npix = (size_t)size * (size_t)size;
    unsigned char* img = (unsigned char*)calloc(npix * 3, 1);
    if (!img) { fprintf(stderr, "error: out of memory\n"); return 3; }
    long k, npts = (long)size * (long)size / 2;
    for (k = 1; k <= npts; k++) {
        double x = halton2((unsigned)k, 2);
        double y = halton2((unsigned)k, 3);
        double x2 = halton2((unsigned)k, 5);
        double y2 = halton2((unsigned)k, 7);
        if (strict && !(x >= 0.0 && x < 1.0 && y >= 0.0 && y < 1.0 &&
                        x2 >= 0.0 && x2 < 1.0 && y2 >= 0.0 && y2 < 1.0)) {
            fprintf(stderr, "error: strict check failed: halton point %ld outside [0,1)\n", k);
            free(img);
            return 3;
        }
        int px = (int)(x * size); if (px >= size) px = size - 1;
        int py = (int)(y * size); if (py >= size) py = size - 1;
        img[((size_t)py * (size_t)size + (size_t)px) * 3 + 0] = 255;
        img[((size_t)py * (size_t)size + (size_t)px) * 3 + 1] = 255;
        img[((size_t)py * (size_t)size + (size_t)px) * 3 + 2] = 255;
        int qx = (int)(x2 * size); if (qx >= size) qx = size - 1;
        int qy = (int)(y2 * size); if (qy >= size) qy = size - 1;
        img[((size_t)qy * (size_t)size + (size_t)qx) * 3 + 0] = 255;
        img[((size_t)qy * (size_t)size + (size_t)qx) * 3 + 1] = 128;
    }
    FILE* f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "error: cannot open '%s' for writing\n", path); free(img); return 2; }
    fprintf(f, "P6\n%d %d\n255\n", size, size);
    size_t wrote = fwrite(img, 1, npix * 3, f);
    free(img);
    if (wrote != npix * 3 || fclose(f) != 0) {
        fprintf(stderr, "error: failed writing '%s'\n", path);
        return 2;
    }
    return 0;
}

static void usage(const char* prog) {
    fprintf(stderr,
        "usage: %s                  run self-tests\n"
        "       %s render OUT.ppm [--size N] [--strict|--no-strict]\n"
        "  N: image width/height in pixels, integer 1..1024 (default 128)\n",
        prog, prog);
}

static int parse_size(const char* s, int* out) {
    long v = 0;
    const char* p;
    if (!s || !*s) return -1;
    for (p = s; *p; p++) {
        if (*p < '0' || *p > '9') return -1;
        v = v * 10 + (*p - '0');
        if (v > 1024) return -1;
    }
    if (v < 1) return -1;
    *out = (int)v;
    return 0;
}

static int render_main(int argc, char** argv) {
    const char* out = NULL;
    int size = 128;
    int strict = 1;
    int i;
    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--size") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "error: --size requires a value\n");
                return 2;
            }
            if (parse_size(argv[++i], &size) != 0) {
                fprintf(stderr, "error: invalid --size '%s': expected integer 1..1024\n", argv[i]);
                return 2;
            }
        } else if (strcmp(argv[i], "--strict") == 0) {
            strict = 1;
        } else if (strcmp(argv[i], "--no-strict") == 0) {
            strict = 0;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "error: unknown option '%s'\n", argv[i]);
            return 2;
        } else if (!out) {
            out = argv[i];
        } else {
            fprintf(stderr, "error: unexpected argument '%s'\n", argv[i]);
            return 2;
        }
    }
    if (!out) {
        fprintf(stderr, "error: render requires an output path\n");
        return 2;
    }
    return render(out, size, strict);
}

int main(int argc, char** argv) {
    if (argc >= 2) {
        if (strcmp(argv[1], "render") == 0) return render_main(argc - 2, argv + 2);
        fprintf(stderr, "error: unknown command '%s'\n", argv[1]);
        usage(argv[0]);
        return 2;
    }
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

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
