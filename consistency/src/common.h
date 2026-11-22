/*
 * common.h: shared scene definition, camera math and framebuffer IO
 * for the rasterizer / ray-tracer consistency check.
 *
 * Both renderers consume the exact same scene file and shading formula.
 * Deliberate precision split: rasterizer works in float32 end to end,
 * tracer intersects in float64 and only quantizes at framebuffer write.
 */
#ifndef CONSISTENCY_COMMON_H
#define CONSISTENCY_COMMON_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_TRI 64
#define FB_CH 5 /* r, g, b, linear view depth, primitive id */
#define FB_BG_DEPTH 1e30f
#define FB_BG_ID -1.0f
#define FB_MAGIC 0x46425546u /* 'FBUF' */
#define PI_F 3.14159265358979323846

typedef struct {
    double a[3], b[3], c[3];
    float col[3];
} Tri;

typedef struct {
    double fovx;          /* horizontal field of view, degrees */
    double light[3];      /* point light position, camera space */
    double ambient;
    Tri tri[MAX_TRI];
    int ntri;
} Scene;

static int scene_load(Scene* s, const char* path) {
    FILE* f = fopen(path, "r");
    if (!f) { perror(path); return -1; }
    s->ntri = 0;
    s->fovx = 60.0;
    s->ambient = 0.15;
    s->light[0] = 1.5; s->light[1] = 2.0; s->light[2] = 0.5;

    char line[512];
    int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        /* strip an inline/full-line comment, then tokenize what remains */
        char* hash = strchr(line, '#');
        if (hash) *hash = '\0';
        char* kw = strtok(line, " \t\r\n");
        if (!kw) continue;

        if (!strcmp(kw, "fovx")) {
            char* v = strtok(NULL, " \t\r\n");
            if (!v || sscanf(v, "%lf", &s->fovx) != 1) goto bad;
        } else if (!strcmp(kw, "light")) {
            char* a = strtok(NULL, " \t\r\n");
            char* b = strtok(NULL, " \t\r\n");
            char* c = strtok(NULL, " \t\r\n");
            if (!a || !b || !c ||
                sscanf(a, "%lf", &s->light[0]) != 1 ||
                sscanf(b, "%lf", &s->light[1]) != 1 ||
                sscanf(c, "%lf", &s->light[2]) != 1) goto bad;
        } else if (!strcmp(kw, "ambient")) {
            char* v = strtok(NULL, " \t\r\n");
            if (!v || sscanf(v, "%lf", &s->ambient) != 1) goto bad;
        } else if (!strcmp(kw, "tri")) {
            if (s->ntri >= MAX_TRI) { fprintf(stderr, "scene too large\n"); fclose(f); return -1; }
            Tri* t = &s->tri[s->ntri++];
            double* dv[9] = { &t->a[0], &t->a[1], &t->a[2],
                              &t->b[0], &t->b[1], &t->b[2],
                              &t->c[0], &t->c[1], &t->c[2] };
            for (int i = 0; i < 9; i++) {
                char* tok = strtok(NULL, " \t\r\n");
                if (!tok || sscanf(tok, "%lf", dv[i]) != 1) goto bad;
            }
            for (int i = 0; i < 3; i++) {
                char* tok = strtok(NULL, " \t\r\n");
                if (!tok || sscanf(tok, "%f", &t->col[i]) != 1) goto bad;
            }
            if (strtok(NULL, " \t\r\n")) goto bad; /* reject trailing junk */
        } else {
            fprintf(stderr, "%s:%d: unknown keyword: %s\n", path, lineno, kw);
            fclose(f);
            return -1;
        }
        continue;
    bad:
        fprintf(stderr, "%s:%d: malformed line starting with %s\n", path, lineno, kw);
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}

/* framebuffer: 4 x uint32 header (magic, w, h, channels), then float32 rows */
static int fb_write(const char* path, int w, int h, const float* pix) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return -1; }
    uint32_t hdr[4] = { FB_MAGIC, (uint32_t)w, (uint32_t)h, FB_CH };
    int ok = fwrite(hdr, sizeof(hdr), 1, f) == 1 &&
             fwrite(pix, sizeof(float) * FB_CH * (size_t)w * h, 1, f) == 1;
    if (fclose(f) != 0) return -1;
    return ok ? 0 : -1;
}

static double focal_px(int w, double fovx_deg) {
    return (w * 0.5) / tan(fovx_deg * PI_F / 360.0);
}

/* float shading path used by the rasterizer */
static inline void shade_f(const float p[3], const float n[3], const float lp[3],
                           float ambient, const float albedo[3], float out[3]) {
    float lx = lp[0] - p[0], ly = lp[1] - p[1], lz = lp[2] - p[2];
    float ll = sqrtf(lx * lx + ly * ly + lz * lz);
    lx /= ll; ly /= ll; lz /= ll;
    float ndl = n[0] * lx + n[1] * ly + n[2] * lz;
    if (ndl < 0.0f) ndl = 0.0f;
    float k = ambient + (1.0f - ambient) * ndl;
    for (int i = 0; i < 3; i++) {
        float v = albedo[i] * k;
        out[i] = v > 1.0f ? 1.0f : (v < 0.0f ? 0.0f : v);
    }
}

/* double shading path used by the ray tracer */
static inline void shade_d(const double p[3], const double n[3], const double lp[3],
                           double ambient, const float albedo[3], double out[3]) {
    double lx = lp[0] - p[0], ly = lp[1] - p[1], lz = lp[2] - p[2];
    double ll = sqrt(lx * lx + ly * ly + lz * lz);
    lx /= ll; ly /= ll; lz /= ll;
    double ndl = n[0] * lx + n[1] * ly + n[2] * lz;
    if (ndl < 0.0) ndl = 0.0;
    double k = ambient + (1.0 - ambient) * ndl;
    for (int i = 0; i < 3; i++) {
        double v = (double)albedo[i] * k;
        out[i] = v > 1.0 ? 1.0 : (v < 0.0 ? 0.0 : v);
    }
}

#endif
