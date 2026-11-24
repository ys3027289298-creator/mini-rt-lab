/* pathtracer: offline back end.
   Unidirectional path tracing, direct light sampled explicitly at every
   bounce (delta point light), next-event estimation, flat environment
   on misses, up to MAX_BOUNCE diffuse/specular bounces.
   All shading calls land in the shared radiometry module. */
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

#ifdef _OPENMP
#include <omp.h>
#endif

#define MAX_BOUNCE 3
#define DEFAULT_SPP 4096

/* van der Corput / Halton style low-discrepancy pixel samples */
static float halton(unsigned i, unsigned base) {
    float f = 1.0f, r = 0.0f;
    while (i > 0) { f /= (float)base; r += f * (float)(i % base); i /= base; }
    return r;
}

static V3 trace_path(const Scene* sc, const Bvh* bvh, Ray ray, unsigned* rng, int depth) {
    Hit hit;
    if (!scene_intersect(sc, bvh, ray, 1e6f, &hit)) return sc->env;

    V3 wo = vscl(ray.d, -1.0f);
    const Material* m = &sc->mat[hit.mat];
    V3 albedo = surface_albedo(sc, &hit);

    V3 L = v3(0, 0, 0);
    for (int k = 0; k < sc->nlight; k++)
        L = vadd(L, shade_direct(sc, bvh, &hit, wo, &sc->light[k]));
    L = vadd(L, ambient_term(sc, &hit, wo));
    if (depth >= MAX_BOUNCE) return L;

    /* one indirect bounce sample */
    V3 wi; float pdf = 1.0f; V3 f;
    if (m->type == MAT_LAMBERT) {
        wi = sample_cosine_hemisphere(hit.n, rng_next(rng), rng_next(rng), &pdf);
        f = brdf_lambert(albedo);
    } else {
        V3 h = sample_ggx(hit.n, m->rough, rng_next(rng), rng_next(rng));
        wi = reflect_v(vscl(wo, -1.0f), h);
        if (vdot(wi, hit.n) <= 0.0f) return L;
        pdf = pdf_ggx(hit.n, wo, wi, m->rough);
        f = brdf_ggx(albedo, m->rough, hit.n, wo, wi);
    }
    if (pdf <= 0.0f) return L;
    float nl = fabsf(vdot(wi, hit.n));
    V3 throughput = vscl(vmul(f, v3(nl, nl, nl)), 1.0f / pdf);

    Ray next = { vadd(hit.p, vscl(hit.n, 2e-2f)), wi };
    V3 Li = trace_path(sc, bvh, next, rng, depth + 1);
    return vadd(L, vmul(throughput, Li));
}

int main(int argc, char** argv) {
    int spp = DEFAULT_SPP;
    int w = 240, h = 180;
    const char* hdr_path = "out/rt_hdr.pfm";
    const char* ldr_path = "out/rt.ppm";
    const char* points_path = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-spp") && i + 1 < argc) spp = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-size") && i + 2 < argc) { w = atoi(argv[++i]); h = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "-hdr") && i + 1 < argc) hdr_path = argv[++i];
        else if (!strcmp(argv[i], "-out") && i + 1 < argc) ldr_path = argv[++i];
        else if (!strcmp(argv[i], "-points") && i + 1 < argc) points_path = argv[++i];
    }

    Scene sc; scene_build(&sc);
    Bvh bvh; bvh_build(&bvh, &sc);
    Camera cam = camera_make(v3(0.0f, 2.2f, 6.2f), v3(0.0f, 0.6f, 0.0f),
                             v3(0, 1, 0), 50.0f, w, h);

    HdrImage img = image_alloc(w, h);
    double* acc = (double*)calloc((size_t)w * h * 3, sizeof(double));
    if (!acc) { fprintf(stderr, "oom\n"); return 1; }

    #pragma omp parallel for schedule(dynamic, 4)
    for (int y = 0; y < h; y++) {
        unsigned rng = 0x9e3779b9u ^ ((unsigned)y * 2654435761u);
        for (int x = 0; x < w; x++) {
            V3 sum = v3(0, 0, 0);
            for (int s = 0; s < spp; s++) {
                /* Halton(2,3) sub-pixel jitter, decorrelated per pixel by rng offset */
                float jx = halton((unsigned)s + (rng & 1023u), 2);
                float jy = halton((unsigned)s + (rng & 1023u), 3);
                Ray ray = camera_ray(&cam, ((float)x + jx) / w, ((float)y + jy) / h);
                V3 L = trace_path(&sc, &bvh, ray, &rng, 0);
                sum = vadd(sum, L);
            }
            V3 mean = vscl(sum, 1.0f / spp);
            img.px[(size_t)y * w + x] = mean;
            acc[((size_t)y * w + x) * 3 + 0] = mean.x;
            acc[((size_t)y * w + x) * 3 + 1] = mean.y;
            acc[((size_t)y * w + x) * 3 + 2] = mean.z;
        }
    }

    unsigned char* ldr = (unsigned char*)malloc(3 * (size_t)w * h);
    tonemap_image(ldr, img.px, w, h);

    /* deterministic centered-ray primary hit points, for geometry parity checks */
    if (points_path) {
        HdrImage pts = image_alloc(w, h);
        for (int y = 0; y < h; y++) for (int x = 0; x < w; x++) {
            Ray ray = camera_ray(&cam, ((float)x + 0.5f) / w, ((float)y + 0.5f) / h);
            Hit hit;
            if (scene_intersect(&sc, &bvh, ray, 1e6f, &hit)) pts.px[(size_t)y * w + x] = hit.p;
        }
        image_write_pfm(&pts, points_path);
        image_free(&pts);
    }

    int ok = image_write_pfm(&img, hdr_path) && image_write_ppm(ldr, w, h, ldr_path);
    printf("pathtracer: spp=%d size=%dx%d hdr=%s ldr=%s ok=%d\n",
           spp, w, h, hdr_path, ldr_path, ok);

    free(acc); free(ldr); image_free(&img);
    return ok ? 0 : 1;
}
