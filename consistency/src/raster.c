/*
 * raster.c: float32 rasterizer.
 * Edge-function fill (top-left rule), perspective-correct interpolation,
 * depth test on linear view-space depth. First triangle wins depth ties.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"

typedef struct { float x, y, z, invw; } Vf;

int main(int argc, char** argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: %s scene.scene width height out.fbuf\n", argv[0]);
        return 2;
    }
    Scene s;
    if (scene_load(&s, argv[1]) != 0) return 2;
    int W = atoi(argv[2]), H = atoi(argv[3]);
    if (W <= 0 || H <= 0) { fprintf(stderr, "bad dimensions\n"); return 2; }

    float* fb = malloc((size_t)W * H * FB_CH * sizeof(float));
    float* depth = malloc((size_t)W * H * sizeof(float));
    if (!fb || !depth) { fprintf(stderr, "oom\n"); return 1; }
    for (int i = 0; i < W * H; i++) {
        depth[i] = FB_BG_DEPTH;
        float* o = &fb[(size_t)i * FB_CH];
        o[0] = o[1] = o[2] = 0.0f;
        o[3] = FB_BG_DEPTH;
        o[4] = FB_BG_ID;
    }

    float focal = (float)focal_px(W, s.fovx);

    for (int ti = 0; ti < s.ntri; ti++) {
        const Tri* T = &s.tri[ti];
        float va[3] = { (float)T->a[0], (float)T->a[1], (float)T->a[2] };
        float vb[3] = { (float)T->b[0], (float)T->b[1], (float)T->b[2] };
        float vc[3] = { (float)T->c[0], (float)T->c[1], (float)T->c[2] };
        const float* v3[3] = { va, vb, vc };

        /* Normal is evaluated in double from the original scene vertices
           so shading matches the tracer's double path; the fill/interp
           path below stays float32 on purpose. */
        double nd[3];
        {
            double e1[3] = { T->b[0] - T->a[0], T->b[1] - T->a[1], T->b[2] - T->a[2] };
            double e2[3] = { T->c[0] - T->a[0], T->c[1] - T->a[1], T->c[2] - T->a[2] };
            nd[0] = e1[1] * e2[2] - e1[2] * e2[1];
            nd[1] = e1[2] * e2[0] - e1[0] * e2[2];
            nd[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double nl = sqrt(nd[0] * nd[0] + nd[1] * nd[1] + nd[2] * nd[2]);
            for (int k = 0; k < 3; k++) nd[k] /= nl;
        }

        Vf sv[3];
        for (int k = 0; k < 3; k++) {
            float z = v3[k][2];
            float invw = (z > 1e-4f) ? 1.0f / z : 1e30f;
            sv[k].x = (W * 0.5f) + focal * v3[k][0] * invw;
            sv[k].y = (H * 0.5f) - focal * v3[k][1] * invw;
            sv[k].z = z;
            sv[k].invw = invw;
        }

        float area = (sv[1].x - sv[0].x) * (sv[2].y - sv[0].y)
                   - (sv[1].y - sv[0].y) * (sv[2].x - sv[0].x);
        if (fabsf(area) <= 1e-7f) continue; /* degenerate only; both windings render */
        if (area < 0.0f) {
            /* back-facing screen winding: reorder for positive fill math.
               The shading normal is oriented separately below using the
               view ray, so it must NOT be flipped here. */
            Vf tv = sv[1]; sv[1] = sv[2]; sv[2] = tv;
            const float* vv = v3[1]; v3[1] = v3[2]; v3[2] = vv;
            area = -area;
        }
        /* two-sided shading: orient the geometric normal toward the camera.
           Primary rays run along +z, so a +z normal faces away and is
           flipped; this is the identical convention used by the tracer. */
        if (nd[2] > 0.0) {
            nd[0] = -nd[0]; nd[1] = -nd[1]; nd[2] = -nd[2];
        }

        int x0 = W, y0 = H, x1 = -1, y1 = -1;
        for (int k = 0; k < 3; k++) {
            int bx = (int)floorf(sv[k].x), by = (int)floorf(sv[k].y);
            if (bx < x0) x0 = bx;
            if (by < y0) y0 = by;
            if (bx > x1) x1 = bx;
            if (by > y1) y1 = by;
        }
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 > W - 1) x1 = W - 1;
        if (y1 > H - 1) y1 = H - 1;
        if (x1 < x0 || y1 < y0) continue;

        /* top-left fill rule: edge is counted when its top/left side is inside */
        float bias[3];
        const Vf* ev[3] = { &sv[0], &sv[1], &sv[2] };
        for (int e = 0; e < 3; e++) {
            const Vf* p = ev[e];
            const Vf* q = ev[(e + 1) % 3];
            float ex = q->x - p->x, ey = q->y - p->y;
            bias[e] = ((ey > 0.0f) || (ey == 0.0f && ex < 0.0f)) ? 0.0f : -1e-6f;
        }

        for (int y = y0; y <= y1; y++) {
            for (int x = x0; x <= x1; x++) {
                float px = (float)x + 0.5f, py = (float)y + 0.5f;
                float w0 = (sv[1].x - sv[0].x) * (py - sv[0].y)
                         - (sv[1].y - sv[0].y) * (px - sv[0].x);
                float w1 = (sv[2].x - sv[1].x) * (py - sv[1].y)
                         - (sv[2].y - sv[1].y) * (px - sv[1].x);
                float w2 = (sv[0].x - sv[2].x) * (py - sv[2].y)
                         - (sv[0].y - sv[2].y) * (px - sv[2].x);
                if (w0 + bias[0] < 0.0f || w1 + bias[1] < 0.0f || w2 + bias[2] < 0.0f)
                    continue;
                float l0 = w0 / area, l1 = w1 / area, l2 = w2 / area;
                float wsum = l0 * sv[0].invw + l1 * sv[1].invw + l2 * sv[2].invw;
                float c0 = l0 * sv[0].invw / wsum;
                float c1 = l1 * sv[1].invw / wsum;
                float c2 = l2 * sv[2].invw / wsum;
                float z = c0 * sv[0].z + c1 * sv[1].z + c2 * sv[2].z;

                int idx = y * W + x;
                if (z >= depth[idx]) continue; /* strict near wins; first triangle keeps ties */
                depth[idx] = z;

                /* shading inputs promoted to double; the barycentric weights
                   c0..c2 still come from the float32 edge-function path, which
                   is the intended precision-class residual */
                double p3d[3];
                for (int k = 0; k < 3; k++)
                    p3d[k] = (double)(c0 * v3[0][k] + c1 * v3[1][k] + c2 * v3[2][k]);
                double cold[3];
                shade_d(p3d, nd, s.light, s.ambient, T->col, cold);
                float col[3] = { (float)cold[0], (float)cold[1], (float)cold[2] };
                float* o = &fb[(size_t)idx * FB_CH];
                o[0] = col[0]; o[1] = col[1]; o[2] = col[2];
                o[3] = z;
                o[4] = (float)ti;
            }
        }
    }

    int rc = fb_write(argv[4], W, H, fb);
    free(fb);
    free(depth);
    return rc == 0 ? 0 : 1;
}
