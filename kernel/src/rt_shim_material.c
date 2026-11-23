/* rt_shim_material.c -- exposes material/src/material.c:cosine_sample()
   without modifying or linking the module's main(). */
#define main material_module_selftest_main
#include "../../material/src/material.c"
#undef main

void rt_shim_cosine_sample(float u1, float u2, float out[3]) {
    V3 v = cosine_sample(u1, u2);
    out[0] = v.x; out[1] = v.y; out[2] = v.z;
}
