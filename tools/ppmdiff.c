/* ppmdiff: compare two binary PPM (P6) images pixel by pixel.
   Prints differing-pixel count and max per-channel delta.
   Exit 0 if identical, 1 if pixels differ, 2 on usage/IO error. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned char* read_ppm(const char* path, int* w, int* h) {
    FILE* f = fopen(path, "rb");
    char magic[3] = { 0, 0, 0 };
    int maxv = 0;
    size_t n;
    unsigned char* buf;
    if (!f) { fprintf(stderr, "ppmdiff: cannot open %s\n", path); exit(2); }
    if (fscanf(f, "%2s", magic) != 1 || strcmp(magic, "P6") != 0) {
        fprintf(stderr, "ppmdiff: %s: not a P6 PPM\n", path); exit(2);
    }
    if (fscanf(f, "%d %d %d", w, h, &maxv) != 3 || maxv != 255) {
        fprintf(stderr, "ppmdiff: %s: bad header\n", path); exit(2);
    }
    fgetc(f); /* single whitespace after maxval */
    n = (size_t)(*w) * (size_t)(*h) * 3;
    buf = (unsigned char*)malloc(n);
    if (!buf || fread(buf, 1, n, f) != n) {
        fprintf(stderr, "ppmdiff: %s: short read\n", path); exit(2);
    }
    fclose(f);
    return buf;
}

int main(int argc, char** argv) {
    int wa, ha, wb, hb;
    unsigned char *a, *b;
    long i, npix, diff_pixels = 0;
    int max_delta = 0;
    if (argc != 3) { fprintf(stderr, "usage: ppmdiff A.ppm B.ppm\n"); return 2; }
    a = read_ppm(argv[1], &wa, &ha);
    b = read_ppm(argv[2], &wb, &hb);
    if (wa != wb || ha != hb) {
        printf("size_mismatch: %dx%d vs %dx%d\n", wa, ha, wb, hb);
        return 1;
    }
    npix = (long)wa * (long)ha;
    for (i = 0; i < npix; i++) {
        int d = 0, c;
        for (c = 0; c < 3; c++) {
            int delta = a[i * 3 + c] - b[i * 3 + c];
            if (delta < 0) delta = -delta;
            if (delta > d) d = delta;
        }
        if (d > 0) diff_pixels++;
        if (d > max_delta) max_delta = d;
    }
    printf("%s vs %s: %dx%d pixels, diff_pixels=%ld, max_channel_delta=%d\n",
           argv[1], argv[2], wa, ha, diff_pixels, max_delta);
    free(a);
    free(b);
    return diff_pixels == 0 ? 0 : 1;
}
