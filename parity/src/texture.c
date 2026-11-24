#include "texture.h"

static float srgb_to_linear(float c) {
    c = clampf(c, 0.0f, 1.0f);
    return c <= 0.04045f ? c / 12.92f : powf((c + 0.055f) / 1.055f, 2.4f);
}

/* Textures are authored in sRGB but stored as linear radiance, so every
   downstream consumer works in one color space. */
void tex_checker(Texture* t, int cells, V3 a, V3 b) {
    t->w = t->h = TEX_SIZE;
    for (int y = 0; y < TEX_SIZE; y++) {
        for (int x = 0; x < TEX_SIZE; x++) {
            int cx = (int)((float)x / TEX_SIZE * cells);
            int cy = (int)((float)y / TEX_SIZE * cells);
            V3 s = ((cx + cy) & 1) ? b : a;
            unsigned char* p = &t->rgba[(y * TEX_SIZE + x) * 4];
            p[0] = (unsigned char)clampf(srgb_to_linear(s.x) * 255.0f + 0.5f, 0, 255);
            p[1] = (unsigned char)clampf(srgb_to_linear(s.y) * 255.0f + 0.5f, 0, 255);
            p[2] = (unsigned char)clampf(srgb_to_linear(s.z) * 255.0f + 0.5f, 0, 255);
            p[3] = 255;
        }
    }
}

V3 tex_rgb(const Texture* t, int x, int y) {
    const unsigned char* p = &t->rgba[(y * t->w + x) * 4];
    return v3(p[0] / 255.0f, p[1] / 255.0f, p[2] / 255.0f);
}

V3 tex_sample(const Texture* t, float u, float v) {
    /* repeat wrap, including for negative coordinates */
    u = u - floorf(u);
    v = v - floorf(v);
    float fx = u * t->w;
    float fy = v * t->h;
    int x0 = (int)floorf(fx), y0 = (int)floorf(fy);
    float tx = fx - x0, ty = fy - y0;
    int x1 = (x0 + 1) % t->w;
    int y1 = (y0 + 1) % t->h;
    x0 %= t->w; y0 %= t->h;

    V3 c00 = tex_rgb(t, x0, y0), c10 = tex_rgb(t, x1, y0);
    V3 c01 = tex_rgb(t, x0, y1), c11 = tex_rgb(t, x1, y1);
    V3 top = vlerp(c00, c10, tx);
    V3 bot = vlerp(c01, c11, tx);
    return vlerp(top, bot, ty);
}
