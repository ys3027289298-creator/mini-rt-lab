#ifndef RT_IFACE_H
#define RT_IFACE_H

/*
 * Interface between the renderer and the modules under test.
 *
 * bvh_shim.c / mat_shim.c `#include` the real source files (bvh/src/bvh.c,
 * material/src/material.c) with their self-test `main` renamed, then export
 * only the entry points below.  The renderer therefore exercises the actual
 * code that the mutation harness patches -- there is no copy to drift.
 */

/* ---- intersection (bvh module) ---- */
/* boxes: n*6 floats, minx,miny,minz,maxx,maxy,maxz per box. 1 = ok. */
int rt_bvh_build(const float* boxes, int n);
/* nearest box hit; returns 1 and writes box index + entry parameter t. */
int rt_bvh_hit(const float o[3], const float d[3], float tmax, int* box, float* t);
int rt_bvh_prim_count(void);
/* Unaccelerated reference path: test every primitive box.  Used by the
   validation tool to check that traversal loses no hits. */
int rt_bvh_hit_brute(const float o[3], const float d[3], float tmax, int* box, float* t);

/* ---- material module ---- */
float rt_brdf_lambert(float albedo);
float rt_ggx_ndf(float nh, float rough);
float rt_cosine_pdf(float cos_theta);
void rt_cosine_sample(float u1, float u2, float out[3]);
void rt_norm3(const float in[3], float out[3]);

#endif
