/* rt_shim.h -- read-only access to the existing module functions.
 *
 * The kernel calls the modules' own code; nothing in bvh/, sample/ or
 * material/ is edited.  Each shim #includes one module source with its entry
 * point renamed, so the module's static helpers stay visible to a thin
 * wrapper in the same translation unit.
 */
#ifndef RT_SHIM_H
#define RT_SHIM_H

/* material/src/material.c: cosine_sample(); the returned V3 is split into
   components because V3 is local to that source file. */
void rt_shim_cosine_sample(float u1, float u2, float out[3]);

/* sample/src/sample.c: f_quad(). */
double rt_shim_f_quad(double u);

#endif /* RT_SHIM_H */
