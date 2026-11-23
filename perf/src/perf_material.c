/*
 * material rendering path for the performance harness.
 *
 * The module source is compiled in UNMODIFIED (only its self-test entry point
 * is renamed).  Each pixel draws exactly `spp` cosine-weighted hemisphere
 * directions from the module's own `cosine_sample()` (two PCG floats per
 * direction) and evaluates the module's Lambert BRDF and GGX NDF on them.
 *
 * Work unit = one material sample.  The counter must equal width*height*spp.
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define main material_module_selftest_main
#include "../../material/src/material.c"
#undef main

#include "../include/perf_render.h"
#include "../include/perf_rng.h"

static uint64_t g_material_work;

void perf_material_render(const PerfRenderParams* p, unsigned char* rgb) {
    g_material_work = 0;
    PerfRng rng;
    perf_rng_seed(&rng, p->seed, p->stream);

    const float invw = 1.0f / (float)p->width;
    const float invh = 1.0f / (float)p->height;

    for (int py = 0; py < p->height; py++) {
        for (int px = 0; px < p->width; px++) {
            /* View angle and roughness vary across the frame. */
            float rough = 0.15f + 0.70f * (float)px * invw;
            float albedo = 0.35f + 0.55f * (float)py * invh;
            float lambert = brdf_lambert(albedo);

            float acc = 0.0f;
            for (int s = 0; s < p->spp; s++) {
                float u1 = perf_rng_float(&rng);
                float u2 = perf_rng_float(&rng);
                V3 wi = cosine_sample(u1, u2);

                float nh = wi.z;
                if (nh < 0.0f) nh = 0.0f;
                float d = ggx_ndf(nh, rough);
                /* cosine-weighted estimator folded with the NDF response */
                acc += lambert * (0.25f + 0.75f * d / (1.0f + d));
                g_material_work++;
            }

            float v = acc / (float)p->spp;
            float lum = 255.0f * (0.20f + 0.80f * (v > 1.0f ? 1.0f : v));
            if (!(lum > 0.0f)) lum = 0.0f;

            size_t off = ((size_t)py * (size_t)p->width + (size_t)px) * 3u;
            rgb[off + 0] = (unsigned char)(lum + 0.5f);
            rgb[off + 1] = (unsigned char)(lum * 0.82f + 0.5f);
            rgb[off + 2] = (unsigned char)(lum * 0.64f + 0.5f);
        }
    }
}

uint64_t perf_material_work(void) { return g_material_work; }
