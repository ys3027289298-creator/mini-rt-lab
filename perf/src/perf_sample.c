/*
 * sample rendering path for the performance harness.
 *
 * The module source is compiled in UNMODIFIED (only its self-test entry point
 * is renamed).  Each pixel performs exactly `spp` stratified draws through the
 * module's own `stratify()` / `f_quad()`; a spatial term from the module's
 * Halton sequence makes the image vary across the frame, with the scene seed
 * choosing the sequence offset.
 *
 * Work unit = one stratified draw.  The counter must equal width*height*spp;
 * anything else is a budget violation even on a fast machine.
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define main sample_module_selftest_main
#include "../../sample/src/sample.c"
#undef main

#include "../include/perf_render.h"

static uint64_t g_sample_work;

void perf_sample_render(const PerfRenderParams* p, unsigned char* rgb) {
    g_sample_work = 0;
    const unsigned stride = (unsigned)p->width;
    unsigned base = (unsigned)((p->seed ^ p->stream) & 0xFFFFu);

    for (int py = 0; py < p->height; py++) {
        for (int px = 0; px < p->width; px++) {
            /* Stratified estimate of integral u^2 over [0,1): 1/3. */
            Stat s = stratify(p->spp, f_quad);
            g_sample_work += (uint64_t)p->spp;

            unsigned idx = base + (unsigned)py * stride + (unsigned)px;
            double h = halton2(idx, 2);
            double h3 = halton2(idx, 3);

            /* Center the estimate around a mid grey and modulate by the
               low-discrepancy sequence so neighbouring pixels differ. */
            double v = 0.5 * s.mean + 0.25 * h + 0.25 * h3;
            double lum = 235.0 * v * (0.6 + 0.4 * (double)px / (double)p->width);

            size_t off = ((size_t)py * (size_t)p->width + (size_t)px) * 3u;
            rgb[off + 0] = (unsigned char)(lum > 255.0 ? 255 : (lum < 0 ? 0 : lum));
            rgb[off + 1] = (unsigned char)((double)rgb[off + 0] * 0.8);
            rgb[off + 2] = (unsigned char)((double)rgb[off + 0] * 0.6);
        }
    }
}

uint64_t perf_sample_work(void) { return g_sample_work; }
