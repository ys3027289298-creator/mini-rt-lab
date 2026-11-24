/* rasterizer: online-preview back end (software stand-in for a GPU pipeline).
   Edge-function triangle setup, perspective-correct interpolation of
   world position / normal / uv, depth test, then the SAME shade_direct +
   ambient_term + texture sampler + tone mapper as the path tracer.
   The only transport difference is structural: no indirect bounces. */
#include "math3.h"
#include "scene.h"
#include "tri.h"
#include "camera.h"
#include "radiometry.h"
#include "image.h"
#include "tonemap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SPHERE_RINGS 24
#define SPHERE_SECTORS 48

/* optional debug: dump the per-pixel shading hit point as an RGB32F PFM */
static void dump_points(const V3* p, int covered, int w, int h, const char* path) {
    HdrImage img = image_alloc(w, h);
    for (int i = 0; i < w * h; i++) img.px[i] = covered ? p[i] : v3(0, 0, 0);
    image_write_pfm(&img, path);
    image_free(&img);
}

typedef struct {
    float p[3];       /* world position */
    float n[3];
    float u, v;
    float inv_z;
} Vary;

int main(int argc, char** argv) {
    int w = 240, h = 180;
    const char* hdr_path = "out/rast_hdr.pfm";
    const char* ldr_path = "out/rast.ppm";
    const char* points_path = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-size") && i + 2 < argc) { w = atoi(argv[++i]); h = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-hdr") && i + 1 < argc) hdr_path = argv[++i];
        else if (!strcmp(argv[i], "-out") && i + 1 < argc) ldr_path = argv[++i];
        else if (!strcmp(argv[i], "-points") && i + 1 < argc) points_path = argv[++i];
    }

    Scene sc; scene_build(&sc);
    /* rasterizer representation of the analytic spheres: tessellation */
    for (int i = 0; i < sc.nsph; i++) {
        if (scene_tessellate(&sc, i, SPHERE_RINGS, SPHERE_SECTORS) != 0) {
            fprintf(stderr, "sphere tessellation exceeded MAX_TRI\n");
            return 1;
        }
    }

    Bvh bvh; bvh_build(&bvh, &sc);
    Camera cam = camera_make(v3(0.0f, 2.2f, 6.2f), v3(0.0f, 0.6f, 0.0f),
                             v3(0, 1, 0), 50.0f, w, h);

    HdrImage img = image_alloc(w, h);
    float* depth = (float*)malloc((size_t)w * h * sizeof(float));
    V3* points = points_path ? (V3*)calloc((size_t)w * h, sizeof(V3)) : NULL;
    for (size_t i = 0; i < (size_t)w * h; i++) depth[i] = 1e30f;

    for (int ti = 0; ti < sc.ntri; ti++) {
        const Tri* t = &sc.tri[ti];
        float sx[3], sy[3], sz[3];
        Vary var[3];
        int visible = 1;
        for (int k = 0; k < 3; k++) {
            sx[k] = sy[k] = sz[k] = 0.0f;
            if (!camera_project(&cam, t->v[k], &sx[k], &sy[k], &sz[k])) visible = 0;
            var[k].p[0] = t->v[k].x; var[k].p[1] = t->v[k].y; var[k].p[2] = t->v[k].z;
            var[k].n[0] = t->n[k].x; var[k].n[1] = t->n[k].y; var[k].n[2] = t->n[k].z;
            var[k].u = t->uv[k][0]; var[k].v = t->uv[k][1];
            var[k].inv_z = sz[k] > 1e-4f ? 1.0f / sz[k] : 0.0f;
        }
        if (!visible) continue;

        float area = (sx[1] - sx[0]) * (sy[2] - sy[0]) -
                     (sx[2] - sx[0]) * (sy[1] - sy[0]);
        if (fabsf(area) < 1e-8f) continue;

        int x0 = (int)floorf(fminf(sx[0], fminf(sx[1], sx[2])));
        int x1 = (int)ceilf(fmaxf(sx[0], fmaxf(sx[1], sx[2])));
        int y0 = (int)floorf(fminf(sy[0], fminf(sy[1], sy[2])));
        int y1 = (int)ceilf(fmaxf(sy[0], fmaxf(sy[1], sy[2])));
        x0 = fmaxf(0, x0); y0 = fmaxf(0, y0);
        x1 = fminf(w - 1, x1); y1 = fminf(h - 1, y1);

        for (int py = y0; py <= y1; py++) {
            for (int px = x0; px <= x1; px++) {
                float cx = px + 0.5f, cy = py + 0.5f;
                float w0 = ((sx[1] - cx) * (sy[2] - cy) - (sx[2] - cx) * (sy[1] - cy)) / area;
                float w1 = ((sx[2] - cx) * (sy[0] - cy) - (sx[0] - cx) * (sy[2] - cy)) / area;
                float w2 = 1.0f - w0 - w1;
                if (w0 < 0 || w1 < 0 || w2 < 0) continue;

                /* perspective-correct interpolation */
                float iw = w0 * var[0].inv_z + w1 * var[1].inv_z + w2 * var[2].inv_z;
                float c0 = w0 * var[0].inv_z / iw;
                float c1 = w1 * var[1].inv_z / iw;
                float c2 = w2 * var[2].inv_z / iw;
                float z = 1.0f / iw;
                size_t idx = (size_t)py * w + px;
                if (z >= depth[idx]) continue;
                depth[idx] = z;

                V3 wp = v3(c0 * var[0].p[0] + c1 * var[1].p[0] + c2 * var[2].p[0],
                           c0 * var[0].p[1] + c1 * var[1].p[1] + c2 * var[2].p[1],
                           c0 * var[0].p[2] + c1 * var[1].p[2] + c2 * var[2].p[2]);
                V3 nn = vnorm(v3(c0 * var[0].n[0] + c1 * var[1].n[0] + c2 * var[2].n[0],
                                 c0 * var[0].n[1] + c1 * var[1].n[1] + c2 * var[2].n[1],
                                 c0 * var[0].n[2] + c1 * var[1].n[2] + c2 * var[2].n[2]));
                Hit hit;
                hit.p = wp; hit.n = nn; hit.t = z; hit.mat = t->mat;
                hit.u = c0 * var[0].u + c1 * var[1].u + c2 * var[2].u;
                hit.v = c0 * var[0].v + c1 * var[1].v + c2 * var[2].v;

                /* view vector straight from the shared camera basis */
                V3 wo = vnorm(vsub(cam.eye, wp));
                if (vdot(nn, wo) < 0.0f) hit.n = vscl(nn, -1.0f);

                V3 L = ambient_term(&sc, &hit, wo);
                for (int k = 0; k < sc.nlight; k++)
                    L = vadd(L, shade_direct(&sc, &bvh, &hit, wo, &sc.light[k]));
                img.px[idx] = L;
                if (points) points[idx] = wp;
            }
        }
    }

    /* background pixels: shared environment radiance through the same
       tone curve rather than a hard-coded clear color */
    for (size_t i = 0; i < (size_t)w * h; i++)
        if (depth[i] >= 1e29f) img.px[i] = sc.env;

    unsigned char* ldr = (unsigned char*)malloc(3 * (size_t)w * h);
    tonemap_image(ldr, img.px, w, h);
    int ok = image_write_pfm(&img, hdr_path) && image_write_ppm(ldr, w, h, ldr_path);
    if (points_path && points) dump_points(points, 1, w, h, points_path);
    printf("rasterizer: size=%dx%d tris=%d hdr=%s ldr=%s ok=%d\n",
           w, h, sc.ntri, hdr_path, ldr_path, ok);

    free(depth); free(ldr); free(points); image_free(&img);
    return ok ? 0 : 1;
}
