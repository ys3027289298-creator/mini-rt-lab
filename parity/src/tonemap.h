/* tonemap: the single display transform for both back ends.
   linear HDR radiance -> ACES(AP1) filmic curve -> sRGB encode -> 8-bit. */
#ifndef PARITY_TONEMAP_H
#define PARITY_TONEMAP_H

#include "math3.h"

float tonemap_curve(float x);
float srgb_encode(float c);
unsigned char tonemap_channel(float c);

/* writes 3 bytes per pixel into dst (size >= w*h*3) */
void tonemap_image(unsigned char* dst, const V3* hdr, int w, int h);

#endif
