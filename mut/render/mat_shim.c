/*
 * mat_shim.c -- compile the real material/src/material.c into the renderer.
 */
#include "iface.h"

#ifndef MAT_SRC
#error "MAT_SRC must be defined to the path of material.c"
#endif

#define main material_module_selftest_main
#include MAT_SRC
#undef main

float rt_brdf_lambert(float albedo) { return brdf_lambert(albedo); }
float rt_ggx_ndf(float nh, float rough) { return ggx_ndf(nh, rough); }
float rt_cosine_pdf(float cos_theta) { return cosine_pdf(cos_theta); }

void rt_cosine_sample(float u1, float u2, float out[3]) {
    V3 v = cosine_sample(u1, u2);
    out[0] = v.x; out[1] = v.y; out[2] = v.z;
}

void rt_norm3(const float in[3], float out[3]) {
    V3 a = { in[0], in[1], in[2] };
    V3 r = norm3(a);
    out[0] = r.x; out[1] = r.y; out[2] = r.z;
}
