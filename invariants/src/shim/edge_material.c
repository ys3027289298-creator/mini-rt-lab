/*
 * edge_material.c -- observer shim for material/src/material.c (compiled unmodified).
 * Pure forwarders to the module's static functions.
 */
#define main material_module_main
#define failures material_module_failures
#include "../../material/src/material.c"
#undef main
#undef failures

#include "inv_impl.h"

double edge_ggx_D(double nh, double rough) { return (double)ggx_ndf((float)nh, (float)rough); }
double edge_lambert(double albedo)         { return (double)brdf_lambert((float)albedo); }
double edge_cosine_pdf(double cos_theta)   { return (double)cosine_pdf((float)cos_theta); }

void edge_cosine_sample(double u1, double u2, double out[3]) {
    V3 v = cosine_sample((float)u1, (float)u2);
    out[0] = v.x; out[1] = v.y; out[2] = v.z;
}

int edge_material_failures(void) { return material_module_failures; }
