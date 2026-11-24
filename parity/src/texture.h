/* texture: one procedural RGBA8 mip-level plus a single sampler contract.
   Both the path tracer and the rasterizer MUST go through tex_sample;
   wrap = repeat, filter = bilinear 4-tap. No backend-specific sampling. */
#ifndef PARITY_TEXTURE_H
#define PARITY_TEXTURE_H

#include "math3.h"

#define TEX_SIZE 128

typedef struct {
    unsigned char rgba[TEX_SIZE * TEX_SIZE * 4];
    int w, h;
} Texture;

void tex_checker(Texture* t, int cells, V3 a, V3 b);
V3 tex_rgb(const Texture* t, int x, int y);

/* u,v in any range: repeat wrapping then bilinear filtering. Linear radiance. */
V3 tex_sample(const Texture* t, float u, float v);

#endif
