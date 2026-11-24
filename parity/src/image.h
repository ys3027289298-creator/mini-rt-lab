/* image: linear HDR float framebuffer plus PFM/PPM IO.
   Both back ends produce the same HDR container, so diffing is format-neutral. */
#ifndef PARITY_IMAGE_H
#define PARITY_IMAGE_H

#include "math3.h"

typedef struct {
    int w, h;
    V3* px;          /* linear radiance, row-major, origin top-left */
} HdrImage;

HdrImage image_alloc(int w, int h);
void image_free(HdrImage* img);

/* PFM: little-endian RGB float, "-1.0" scale flag. */
int image_write_pfm(const HdrImage* img, const char* path);
int image_read_pfm(HdrImage* img, const char* path);

/* PPM: 8-bit sRGB bytes already produced by the shared tone mapper. */
int image_write_ppm(const unsigned char* rgb, int w, int h, const char* path);
int image_read_ppm(unsigned char** rgb, int* w, int* h, const char* path);

#endif
