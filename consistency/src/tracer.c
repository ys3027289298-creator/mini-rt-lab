/*
 * tracer.c: float64 ray tracer.
 * One primary ray through each pixel center, Moller-Trumbore intersection,
 * nearest hit. Depth stored is the linear view-space z (same convention as
 * raster.c) so depth deltas are directly comparable.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "common.h"

int main(int argc, char** argv) {
    if (argc != 5) {
        fprintf(stderr, "usage: %s scene.scene width height out.fbuf\n", argv[0]);
        return 2;
    }
    Scene s;
    if (scene_load(&s, argv[1]) != 0) return 2;
    int W = atoi(argv[2]), H = atoi(argv[3]);
    if (W <= 0 || H <= 0) { fprintf(stderr, "bad dimensions\n"); return 2; }

    float* fb = calloc((size_t)W * H * FB_CH, sizeof(float));
    if (!fb) { fprintf(stderr, "oom\n"); return 1; }

    double focal = focal_px(W, s.fovx);
    double aspect = (double)H / (double)W;

    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            double ndcx = ((double)x + 0.5) / (double)W - 0.5;
            double ndcy = 0.5 - ((double)y + 0.5) / (double)H;
            double dx = ndcx;
            double dy = ndcy * aspect;
            double dz = focal / W; /* == 1/(2 tan(fovx/2)): same framing as raster */
            double dl = sqrt(dx * dx + dy * dy + dz * dz);
            dx /= dl; dy /= dl; dz /= dl;

            double best_t = 1e30;
            int best_id = -1;
            double best_u = 0.0, best_v = 0.0;
            for (int ti = 0; ti < s.ntri; ti++) {
                const Tri* T = &s.tri[ti];
                double e1[3], e2[3], p[3], q[3], tvec[3];
                for (int i = 0; i < 3; i++) {
                    e1[i] = T->b[i] - T->a[i];
                    e2[i] = T->c[i] - T->a[i];
                }
                p[0] = dy * e2[2] - dz * e2[1];
                p[1] = dz * e2[0] - dx * e2[2];
                p[2] = dx * e2[1] - dy * e2[0];
                double det = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];
                if (fabs(det) <= 1e-12) continue; /* degenerate; both windings render */
                double inv = 1.0 / det;
                for (int i = 0; i < 3; i++) tvec[i] = -T->a[i];
                double u = (tvec[0] * p[0] + tvec[1] * p[1] + tvec[2] * p[2]) * inv;
                if (u < 0.0 || u > 1.0) continue;
                q[0] = tvec[1] * e1[2] - tvec[2] * e1[1];
                q[1] = tvec[2] * e1[0] - tvec[0] * e1[2];
                q[2] = tvec[0] * e1[1] - tvec[1] * e1[0];
                double v = (dx * q[0] + dy * q[1] + dz * q[2]) * inv;
                if (v < 0.0 || u + v > 1.0) continue;
                double t = (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]) * inv;
                if (t > 1e-6 && t < best_t) {
                    best_t = t; best_id = ti; best_u = u; best_v = v;
                }
            }

            float* o = &fb[((size_t)y * W + x) * FB_CH];
            if (best_id < 0) {
                o[3] = FB_BG_DEPTH;
                o[4] = FB_BG_ID;
                continue;
            }
            const Tri* T = &s.tri[best_id];
            double hit[3];
            for (int i = 0; i < 3; i++)
                hit[i] = T->a[i] + best_u * (T->b[i] - T->a[i])
                               + best_v * (T->c[i] - T->a[i]);
            double e1[3], e2[3], n[3];
            for (int i = 0; i < 3; i++) {
                e1[i] = T->b[i] - T->a[i];
                e2[i] = T->c[i] - T->a[i];
            }
            n[0] = e1[1] * e2[2] - e1[2] * e2[1];
            n[1] = e1[2] * e2[0] - e1[0] * e2[2];
            n[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            for (int i = 0; i < 3; i++) n[i] /= nl;
            /* two-sided shading: face the camera, matching the raster normal */
            if (n[0] * dx + n[1] * dy + n[2] * dz > 0.0) {
                n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
            }

            double col[3];
            shade_d(hit, n, s.light, s.ambient, T->col, col);
            o[0] = (float)col[0];
            o[1] = (float)col[1];
            o[2] = (float)col[2];
            o[3] = (float)hit[2]; /* linear view depth, identical convention to raster */
            o[4] = (float)best_id;
        }
    }

    int rc = fb_write(argv[4], W, H, fb);
    free(fb);
    return rc == 0 ? 0 : 1;
}
