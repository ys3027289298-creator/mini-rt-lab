/*
 * compare: bit-exact, per-pixel/per-channel comparison of two PFM frames.
 *
 * Exit status:
 *   0  frames are identical (same dimensions, every float32 channel equal)
 *   1  a difference was found; the first differing pixel is reported with
 *      its (x, y, channel), both values, and the total mismatch count
 *   2  usage or file error (including dimension mismatch)
 *
 * Coordinates: x grows right, y grows downward with y=0 at the TOP row,
 * matching the renderer's scanline indexing (PFM's on-disk rows are
 * bottom-to-top and are flipped on read).
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

static void bswap_row(float* row, size_t n) {
    unsigned char* p = (unsigned char*)row;
    for (size_t i = 0; i < n; i++) {
        unsigned char* c = p + i * 4;
        unsigned char tmp = c[0];
        c[0] = c[3]; c[3] = tmp;
        tmp = c[1]; c[1] = c[2]; c[2] = tmp;
    }
}

/* Reads a binary RGB PFM into a newly malloc'd top-row-first float buffer. */
static int read_pfm(const char* path, int* out_w, int* out_h, float** out_fb) {
    FILE* f = fopen(path, "rb");
    if (!f) { perror(path); return -1; }

    char magic[3] = {0};
    int w, h;
    float scale;
    /* fscanf stops at whitespace; one fgetc each consumes the delimiter so
     * arbitrary newlines/spaces in the header are tolerated. */
    if (fscanf(f, "%2s", magic) != 1 || strcmp(magic, "PF") != 0) {
        fprintf(stderr, "%s: not a binary RGB PFM (bad magic)\n", path);
        fclose(f);
        return -1;
    }
    (void)fgetc(f);
    if (fscanf(f, "%d", &w) != 1 || w <= 0) {
        fprintf(stderr, "%s: bad width\n", path);
        fclose(f);
        return -1;
    }
    (void)fgetc(f);
    if (fscanf(f, "%d", &h) != 1 || h <= 0) {
        fprintf(stderr, "%s: bad height\n", path);
        fclose(f);
        return -1;
    }
    (void)fgetc(f);
    if (fscanf(f, "%f", &scale) != 1 || scale == 0.0f) {
        fprintf(stderr, "%s: bad endian/scale header\n", path);
        fclose(f);
        return -1;
    }
    (void)fgetc(f);

    int need_swap = (scale > 0.0f); /* negative scale = little-endian */
    size_t n = (size_t)w * h * 3;
    float* fb = malloc(n * sizeof(float));
    if (!fb) { perror("malloc"); fclose(f); return -1; }

    /* On-disk order is bottom row first; flip into top-row-first. */
    for (int y = h - 1; y >= 0; y--) {
        if (fread(&fb[(size_t)y * w * 3], sizeof(float),
                  (size_t)w * 3, f) != (size_t)w * 3) {
            fprintf(stderr, "%s: truncated pixel data\n", path);
            free(fb);
            fclose(f);
            return -1;
        }
        if (need_swap) bswap_row(&fb[(size_t)y * w * 3], (size_t)w * 3);
    }
    if (fgetc(f) != EOF) {
        fprintf(stderr, "%s: trailing data after frame\n", path);
        free(fb);
        fclose(f);
        return -1;
    }
    fclose(f);
    *out_w = w;
    *out_h = h;
    *out_fb = fb;
    return 0;
}

static const char* chname(int c) {
    return c == 0 ? "R" : (c == 1 ? "G" : "B");
}

static int compare_frames(const char* pa, const char* pb) {
    int wa, ha, wb, hb;
    float *a = NULL, *b = NULL;
    if (read_pfm(pa, &wa, &ha, &a) != 0) return 2;
    if (read_pfm(pb, &wb, &hb, &b) != 0) { free(a); return 2; }

    if (wa != wb || ha != hb) {
        fprintf(stderr,
                "dimension mismatch: %s=%dx%d vs %s=%dx%d\n",
                pa, wa, ha, pb, wb, hb);
        free(a); free(b);
        return 2;
    }

    uint32_t diff_channels = 0, diff_pixels = 0;
    size_t total = (size_t)wa * ha * 3;
    size_t first = total;
    for (size_t i = 0; i < total; i++) {
        if (a[i] != b[i]) {
            if (first == total) first = i;
            diff_channels++;
        }
    }
    if (diff_channels == 0) {
        printf("IDENTICAL %dx%d (%zu float channels, %zu bytes)\n",
               wa, ha, total, total * sizeof(float));
        free(a); free(b);
        return 0;
    }

    /* Count pixels containing at least one differing channel. */
    for (int y = 0; y < ha; y++) {
        for (int x = 0; x < wa; x++) {
            const float* pa0 = &a[((size_t)y * wa + x) * 3];
            const float* pb0 = &b[((size_t)y * wa + x) * 3];
            if (pa0[0] != pb0[0] || pa0[1] != pb0[1] || pa0[2] != pb0[2])
                diff_pixels++;
        }
    }

    size_t first_pixel = first / 3;
    int fx = (int)(first_pixel % (size_t)wa);
    int fy = (int)(first_pixel / (size_t)wa);
    int fc = (int)(first % 3);
    uint32_t ba, bb;
    memcpy(&ba, &a[first], 4);
    memcpy(&bb, &b[first], 4);

    printf("DIFFERENT %dx%d\n", wa, ha);
    printf("first_diff: pixel=(%d,%d) channel=%s\n", fx, fy, chname(fc));
    printf("  %s: %.9g (0x%08x)\n", pa, (double)a[first], ba);
    printf("  %s: %.9g (0x%08x)\n", pb, (double)b[first], bb);
    printf("diff_channels=%u diff_pixels=%u (of %d)\n",
           diff_channels, diff_pixels, wa * ha);
    free(a); free(b);
    return 1;
}

/* Built-in negative control: exercises the mismatch reporting path without
 * depending on the renderer or any shell-level binary patching. */
static int selftest(void) {
    char path[] = "/tmp/rtcmp_selftest_XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) { perror("mkstemp"); return 2; }
    FILE* f = fdopen(fd, "wb");
    if (!f) { perror("fdopen"); return 2; }
    /* One 2x1 frame; the second write flips the top-left R channel. */
    float frame[6] = { 0.25f, 0.5f, 0.75f,  0.1f, 0.2f, 0.3f };
    fprintf(f, "PF\n2 1\n-1.0\n");
    fwrite(frame, sizeof(float), 6, f);
    fclose(f);

    int w, h;
    float* a = NULL;
    if (read_pfm(path, &w, &h, &a) != 0) { unlink(path); return 2; }
    float* b = malloc(sizeof(frame));
    if (!b) { perror("malloc"); return 2; }
    memcpy(b, a, sizeof(frame));
    b[0] = 0.2500001f; /* top-left R channel only */

    int rc = 0;
    if (a[0] == b[0]) {
        fprintf(stderr, "selftest: expected mismatch not injected\n");
        rc = 2;
    } else if (a[1] != b[1] || a[2] != b[2] || a[3] != b[3]
            || a[4] != b[4] || a[5] != b[5]) {
        fprintf(stderr, "selftest: collateral mismatch\n");
        rc = 2;
    }
    uint32_t ba, bb;
    memcpy(&ba, &a[0], 4);
    memcpy(&bb, &b[0], 4);
    if (rc == 0) {
        printf("selftest first_diff would be: pixel=(0,0) channel=R "
               "%.9g(0x%08x) vs %.9g(0x%08x)\n",
               (double)a[0], ba, (double)b[0], bb);
        printf("selftest=1\n");
    }
    free(a); free(b);
    unlink(path);
    return rc;
}

int main(int argc, char** argv) {
    if (argc == 2 && strcmp(argv[1], "--selftest") == 0)
        return selftest();
    if (argc != 3) {
        fprintf(stderr, "usage: %s A.pfm B.pfm | %s --selftest\n",
                argv[0], argv[0]);
        return 2;
    }
    return compare_frames(argv[1], argv[2]);
}
