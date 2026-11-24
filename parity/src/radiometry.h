/* radiometry: the shared shading contract for both renderers.
   Every function here is the ONLY place lighting math is defined.
   Units: radiance [W/(sr*m^2)] everywhere; point lights carry radiant
   intensity [W/sr] and fall off with the inverse-square law. */
#ifndef PARITY_RADIOMETRY_H
#define PARITY_RADIOMETRY_H

#include "scene.h"
#include "tri.h"

/* surface albedo: texture lookup through the shared sampler, or base color */
V3 surface_albedo(const Scene* sc, const Hit* hit);

/* Lambert BRDF: albedo / pi [1/sr] */
V3 brdf_lambert(V3 albedo);

/* GGX microfacet: D (GGX), G (Smith), F (Schlick) */
float ggx_d(float nh, float rough);
float ggx_g1(float nv, float rough);
float ggx_g(float nv, float nl, float rough);
V3 fresnel_schlick(V3 f0, float u);
V3 brdf_ggx(V3 albedo, float rough, V3 n, V3 wo, V3 wi);

/* direct lighting from one delta point light, including the shared BVH
   shadow-ray visibility test. Returns radiance; 0 if occluded. */
V3 shade_direct(const Scene* sc, const Bvh* bvh, const Hit* hit, V3 wo,
                const PointLight* light);

/* closed-form environment contribution; same term in both pipelines.
   The path tracer also samples the environment on ray misses, which
   converges to exactly this for flat environments. */
V3 ambient_term(const Scene* sc, const Hit* hit, V3 wo);

/* importance sampling used only by the path tracer bounce steps */
V3 sample_cosine_hemisphere(V3 n, float u1, float u2, float* pdf);
V3 sample_ggx(V3 n, float rough, float u1, float u2);
float pdf_ggx(V3 n, V3 wo, V3 wi, float rough);

#endif
