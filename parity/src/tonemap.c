#include "tonemap.h"
#include <stddef.h>

/* Narkowicz ACES filmic approximation; identical curve in both pipelines. */
float tonemap_curve(float x) {
    const float a = 2.51f, b = 0.03f, c = 2.43f, d = 0.59f, e = 0.14f;
    x = clampf(x, 0.0f, 1e6f);
    return clampf((x * (a * x + b)) / (x * (c * x + d) + e), 0.0f, 1.0f);
}

float srgb_encode(float c) {
    c = clampf(c, 0.0f, 1.0f);
    return c <= 0.0031308f ? 12.92f * c
                           : 1.055f * powf(c, 1.0f / 2.4f) - 0.055f;
}

unsigned char tonemap_channel(float c) {
    return (unsigned char)clampf(srgb_encode(tonemap_curve(c)) * 255.0f + 0.5f, 0.0f, 255.0f);
}

void tonemap_image(unsigned char* dst, const V3* hdr, int w, int h) {
    size_t n = (size_t)w * h;
    for (size_t i = 0; i < n; i++) {
        dst[i * 3 + 0] = tonemap_channel(hdr[i].x);
        dst[i * 3 + 1] = tonemap_channel(hdr[i].y);
        dst[i * 3 + 2] = tonemap_channel(hdr[i].z);
    }
}
