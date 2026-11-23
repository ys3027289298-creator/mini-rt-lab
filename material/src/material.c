
/* material: Lambert + GGX BRDF evaluation and cosine-weighted sampling,
   with a PDF check against the analytic cosine distribution. */
#include <stdio.h>
#include <math.h>
#include <string.h>

typedef struct { float x, y, z; } V3;

static float dot3(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 norm3(V3 a) {
    float l = sqrtf(dot3(a, a));
    V3 r = { 0, 0, 0 };
    if (l > 0) { r.x = a.x / l; r.y = a.y / l; r.z = a.z / l; }
    return r;
}

/* Lambert: albedo / pi */
static float brdf_lambert(float albedo) { return albedo / 3.14159265f; }

/* GGX normal distribution */
static float ggx_ndf(float nh, float rough) {
    float a = rough * rough;
    float a2 = a * a;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / (3.14159265f * d * d);
}

/* cosine-weighted hemisphere: pdf = cos/pi */
static float cosine_pdf(float cos_theta) { return cos_theta / 3.14159265f; }

static V3 cosine_sample(float u1, float u2) {
    float r = sqrtf(u1);
    float phi = 2.0f * 3.14159265f * u2;
    V3 v = { r * cosf(phi), r * sinf(phi), sqrtf(1.0f - u1 > 0 ? 1.0f - u1 : 0) };
    return norm3(v);
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

/* deterministic render: unit sphere at origin, orthographic camera,
   Lambert + GGX specular. Only +,-,*,/ and sqrtf (IEEE correctly rounded),
   so the image is bit-identical across compilers, optimization levels
   and machines. */
static int render(const char* path, int size, int strict) {
    V3 l = { 0.5f, 0.7f, 1.0f };
    V3 view = { 0.0f, 0.0f, 1.0f };
    V3 h;
    FILE* f;
    int px, py;
    l = norm3(l);
    h = norm3((V3){ l.x + view.x, l.y + view.y, l.z + view.z });
    f = fopen(path, "wb");
    if (!f) { fprintf(stderr, "error: cannot open '%s' for writing\n", path); return 2; }
    fprintf(f, "P6\n%d %d\n255\n", size, size);
    for (py = 0; py < size; py++) {
        for (px = 0; px < size; px++) {
            float u = ((float)px + 0.5f) / (float)size * 2.4f - 1.2f;
            float v = 1.2f - ((float)py + 0.5f) / (float)size * 2.4f;
            unsigned char rgb[3] = { 12, 12, 16 };
            float r2 = u * u + v * v;
            if (r2 <= 1.0f) {
                V3 n = { u, v, sqrtf(1.0f - r2) };
                float ndl = dot3(n, l); if (ndl < 0.0f) ndl = 0.0f;
                float nh = dot3(n, h); if (nh < 0.0f) nh = 0.0f;
                float dif = brdf_lambert(0.8f) * ndl;
                float sp = ggx_ndf(nh, 0.35f) * ndl * 0.08f;
                float cr = dif * 1.00f + sp;
                float cg = dif * 0.85f + sp;
                float cb = dif * 0.70f + sp;
                if (strict && !(cr >= 0.0f && cg >= 0.0f && cb >= 0.0f &&
                                cr < 10.0f && cg < 10.0f && cb < 10.0f)) {
                    fprintf(stderr, "error: strict check failed at pixel (%d,%d): color out of range\n", px, py);
                    fclose(f);
                    return 3;
                }
                if (cr > 1.0f) cr = 1.0f;
                if (cg > 1.0f) cg = 1.0f;
                if (cb > 1.0f) cb = 1.0f;
                rgb[0] = (unsigned char)(cr * 255.0f);
                rgb[1] = (unsigned char)(cg * 255.0f);
                rgb[2] = (unsigned char)(cb * 255.0f);
            }
            fwrite(rgb, 1, 3, f);
        }
    }
    if (fclose(f) != 0) { fprintf(stderr, "error: failed writing '%s'\n", path); return 2; }
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
    check("lambert_unit", fabsf(brdf_lambert(3.14159265f) - 1.0f) < 1e-4f);
    check("lambert_zero", brdf_lambert(0.0f) == 0.0f);

    check("ggx_peak", ggx_ndf(1.0f, 0.5f) > ggx_ndf(0.2f, 0.5f));
    check("ggx_positive", ggx_ndf(0.5f, 0.3f) > 0.0f);

    check("pdf_half", fabsf(cosine_pdf(0.5f) - 0.5f / 3.14159265f) < 1e-6f);

    V3 s = cosine_sample(0.5f, 0.25f);
    check("sample_unit", fabsf(sqrtf(dot3(s, s)) - 1.0f) < 1e-4f);
    check("sample_z_nonneg", s.z >= 0.0f);

    V3 up = { 0, 0, 1 };
    check("sample_hemisphere", dot3(s, up) >= -1e-6f);

    /* mean of cos over samples ~ 2/3 for cosine sampling */
    double acc = 0.0; int n = 512;
    for (int i = 0; i < n; i++) {
        float u1 = (float)((i % 16) + 0.5f) / 16.0f;
        float u2 = (float)((i / 16) + 0.5f) / 32.0f;
        V3 v = cosine_sample(u1, u2);
        acc += v.z;
    }
    check("mean_cos_close", fabs(acc / n - 2.0 / 3.0) < 0.02);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
