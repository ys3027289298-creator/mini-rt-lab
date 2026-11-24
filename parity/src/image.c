#include "image.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

HdrImage image_alloc(int w, int h) {
    HdrImage img;
    img.w = w; img.h = h;
    img.px = (V3*)calloc((size_t)w * h, sizeof(V3));
    return img;
}

void image_free(HdrImage* img) {
    free(img->px);
    img->px = NULL;
    img->w = img->h = 0;
}

/* PFM on disk is bottom-up; the in-memory HdrImage is top-left, so rows are
   flipped on write. image_read_pfm flips them back, making IO an identity. */
int image_write_pfm(const HdrImage* img, const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    fprintf(f, "PF\n%d %d\n-1.0\n", img->w, img->h);
    size_t row = 3 * (size_t)img->w;
    int ok = 1;
    for (int y = img->h - 1; y >= 0; y--) {
        if (fwrite(img->px + (size_t)y * img->w, sizeof(float), row, f) != row) ok = 0;
    }
    fclose(f);
    return ok;
}

/* PFM tokens: magic, width, height, scale (sign = endianness). */
static int pfm_header(FILE* f, int* w, int* h, int* little) {
    char magic[8];
    if (fscanf(f, "%7s", magic) != 1 || strcmp(magic, "PF") != 0) return 0;
    int c;
    do { c = fgetc(f); } while (c == '#');           /* skip start of line / comments */
    ungetc(c, f);
    float scale;
    if (fscanf(f, "%d %d %f", w, h, &scale) != 3) return 0;
    fgetc(f);                                       /* single whitespace after header */
    *little = scale < 0.0f;
    return 1;
}

int image_read_pfm(HdrImage* img, const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return 0;
    int w, h, little;
    if (!pfm_header(f, &w, &h, &little)) { fclose(f); return 0; }
    *img = image_alloc(w, h);
    size_t n = 3 * (size_t)w * h;
    if (fread(img->px, sizeof(float), n, f) != n) { fclose(f); return 0; }
    fclose(f);
    if (!little) {
        unsigned char* b = (unsigned char*)img->px;
        for (size_t i = 0; i < n * sizeof(float); i += 4) {
            unsigned char t = b[i]; b[i] = b[i + 3]; b[i + 3] = t;
            t = b[i + 1]; b[i + 1] = b[i + 2]; b[i + 2] = t;
        }
    }
    /* PFM stores rows bottom-up; flip to our top-left convention. */
    for (int y = 0; y < h / 2; y++) {
        V3* a = img->px + (size_t)y * w;
        V3* b2 = img->px + (size_t)(h - 1 - y) * w;
        for (int x = 0; x < w; x++) { V3 t = a[x]; a[x] = b2[x]; b2[x] = t; }
    }
    return 1;
}

int image_write_ppm(const unsigned char* rgb, int w, int h, const char* path) {
    FILE* f = fopen(path, "wb");
    if (!f) return 0;
    fprintf(f, "P6\n%d %d\n255\n", w, h);
    if (fwrite(rgb, 1, 3 * (size_t)w * h, f) != 3 * (size_t)w * h) {
        fclose(f); return 0;
    }
    fclose(f);
    return 1;
}

/* Read a header token; on return the single delimiter whitespace has been
   consumed, so the caller must NOT skip another byte before binary data. */
static int ppm_token(FILE* f, char* tok, int cap) {
    int c;
    for (;;) {
        do { c = fgetc(f); } while (c != EOF && (c == ' ' || c == '\t' ||
                                                 c == '\n' || c == '\r'));
        if (c == '#') {                       /* comment: skip to end of line */
            do { c = fgetc(f); } while (c != EOF && c != '\n');
            continue;
        }
        break;
    }
    if (c == EOF) return 0;
    int n = 0;
    while (c != EOF && c != ' ' && c != '\t' && c != '\n' && c != '\r') {
        if (n < cap - 1) tok[n++] = (char)c;
        c = fgetc(f);
    }
    tok[n] = 0;
    return n > 0;
}

int image_read_ppm(unsigned char** rgb, int* w, int* h, const char* path) {
    FILE* f = fopen(path, "rb");
    if (!f) return 0;
    char tok[16];
    if (!ppm_token(f, tok, sizeof tok) || strcmp(tok, "P6") != 0) { fclose(f); return 0; }
    if (!ppm_token(f, tok, sizeof tok)) { fclose(f); return 0; }
    int iw = atoi(tok);
    if (!ppm_token(f, tok, sizeof tok)) { fclose(f); return 0; }
    int ih = atoi(tok);
    if (!ppm_token(f, tok, sizeof tok)) { fclose(f); return 0; }
    int maxv = atoi(tok);
    if (iw <= 0 || ih <= 0 || maxv != 255) { fclose(f); return 0; }
    /* ppm_token above consumed exactly the one whitespace after maxval. */
    unsigned char* buf = (unsigned char*)malloc(3 * (size_t)iw * ih);
    if (!buf || fread(buf, 1, 3 * (size_t)iw * ih, f) != 3 * (size_t)iw * ih) {
        free(buf); fclose(f); return 0;
    }
    fclose(f);
    *rgb = buf; *w = iw; *h = ih;
    return 1;
}
