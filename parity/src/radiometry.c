#include "radiometry.h"

V3 surface_albedo(const Scene* sc, const Hit* hit) {
    const Material* m = &sc->mat[hit->mat];
    if (m->tex >= 0) return vmul(m->base, tex_sample(&sc->tex[m->tex], hit->u, hit->v));
    return m->base;
}

V3 brdf_lambert(V3 albedo) { return vscl(albedo, INV_PI); }

float ggx_d(float nh, float rough) {
    float a2 = rough * rough; a2 *= a2;
    float d = nh * nh * (a2 - 1.0f) + 1.0f;
    return a2 / (PI * d * d);
}

float ggx_g1(float nv, float rough) {
    float a2 = rough * rough; a2 *= a2;
    return 2.0f * nv / (nv + sqrtf(a2 + (1.0f - a2) * nv * nv));
}

float ggx_g(float nv, float nl, float rough) {
    return ggx_g1(nv, rough) * ggx_g1(nl, rough);
}

V3 fresnel_schlick(V3 f0, float u) {
    return vadd(f0, vmul(vsub(v3(1, 1, 1), f0), v3(powf(1.0f - u, 5.0f),
                                                      powf(1.0f - u, 5.0f),
                                                      powf(1.0f - u, 5.0f))));
}

V3 brdf_ggx(V3 albedo, float rough, V3 n, V3 wo, V3 wi) {
    V3 h = vnorm(vadd(wo, wi));
    float nh = clampf(vdot(n, h), 0.0f, 1.0f);
    float nv = clampf(vdot(n, wo), 1e-4f, 1.0f);
    float nl = clampf(vdot(n, wi), 0.0f, 1.0f);
    if (nl <= 0.0f) return v3(0, 0, 0);
    float d = ggx_d(nh, rough);
    float g = ggx_g(nv, nl, rough);
    /* half-metal approximation: tint the dielectric 4% base with albedo */
    V3 f0 = vlerp(v3(0.04f, 0.04f, 0.04f), albedo, 0.5f);
    V3 f = fresnel_schlick(f0, nh);
    float denom = 4.0f * nv * nl + 1e-5f;
    V3 spec = vscl(f, d * g / denom);
    /* energy-preserving diffuse: (1-F) * Lambert, weighted by metalness 0 */
    V3 kd = vsub(v3(1, 1, 1), f);
    V3 diff = vscl(vmul(kd, albedo), INV_PI);
    return vadd(spec, diff);
}

V3 shade_direct(const Scene* sc, const Bvh* bvh, const Hit* hit, V3 wo,
                const PointLight* light) {
    V3 to_l = vsub(light->pos, hit->p);
    float dist = vlen(to_l);
    V3 wi = vscl(to_l, 1.0f / dist);
    float nl = vdot(hit->n, wi);
    if (nl <= 0.0f) return v3(0, 0, 0);

    /* Shadow-ray bias. A tiny fixed offset is smaller than the chord error of
       a tessellated curved surface, so the ray re-hits a neighboring triangle
       of the SAME surface (shadow acne) and kills real highlights. The offset
       is scene-scale relative and shared verbatim by both back ends. */
    const float shadow_bias = 2e-2f;
    V3 origin = vadd(hit->p, vscl(hit->n, shadow_bias));
    Ray shadow = { origin, wi };
    float t_end = dist - shadow_bias;
    if (t_end <= 0.0f) return v3(0, 0, 0);
    if (scene_occluded(sc, bvh, shadow, t_end)) return v3(0, 0, 0);

    /* irradiance E = I / r^2 * cos(theta); reflected radiance L = f * E */
    V3 irradiance = vscl(light->intensity, nl / (dist * dist));
    const Material* m = &sc->mat[hit->mat];
    V3 albedo = surface_albedo(sc, hit);
    V3 f = (m->type == MAT_LAMBERT)
               ? brdf_lambert(albedo)
               : brdf_ggx(albedo, m->rough, hit->n, wo, wi);
    return vmul(f, irradiance);
}

V3 ambient_term(const Scene* sc, const Hit* hit, V3 wo) {
    const Material* m = &sc->mat[hit->mat];
    V3 albedo = surface_albedo(sc, hit);
    float nv = clampf(vdot(hit->n, wo), 0.0f, 1.0f);
    if (m->type == MAT_LAMBERT) {
        /* integral of Lambert over the hemisphere = albedo, times Le */
        return vmul(albedo, sc->env);
    }
    V3 f0 = vlerp(v3(0.04f, 0.04f, 0.04f), albedo, 0.5f);
    V3 f = fresnel_schlick(f0, nv);
    V3 kd = vsub(v3(1, 1, 1), f);
    /* diffuse hemisphere integral + mirror-ish environment reflection */
    return vadd(vmul(vmul(kd, albedo), sc->env), vmul(f, sc->env));
}

static V3 tangent_basis(V3 n) {
    V3 helper = fabsf(n.x) > 0.9f ? v3(0, 1, 0) : v3(1, 0, 0);
    return vnorm(vcross(n, helper));
}

V3 sample_cosine_hemisphere(V3 n, float u1, float u2, float* pdf) {
    float r = sqrtf(u1);
    float phi = 2.0f * PI * u2;
    V3 local = v3(r * cosf(phi), sqrtf(fmaxf(0.0f, 1.0f - u1)), r * sinf(phi));
    V3 t = tangent_basis(n);
    V3 b = vcross(n, t);
    V3 wi = vnorm(vadd(vadd(vscl(t, local.x), vscl(n, local.y)), vscl(b, local.z)));
    if (pdf) *pdf = fmaxf(vdot(n, wi), 0.0f) * INV_PI;
    return wi;
}

V3 sample_ggx(V3 n, float rough, float u1, float u2) {
    float a2 = rough * rough; a2 *= a2;
    float phi = 2.0f * PI * u2;
    float cos2 = (1.0f - u1) / (1.0f + (a2 - 1.0f) * u1);
    float cos_t = sqrtf(fmaxf(0.0f, cos2));
    float sin_t = sqrtf(fmaxf(0.0f, 1.0f - cos2));
    V3 local_h = v3(sin_t * cosf(phi), cos_t, sin_t * sinf(phi));
    V3 t = tangent_basis(n);
    V3 b = vcross(n, t);
    return vnorm(vadd(vadd(vscl(t, local_h.x), vscl(n, local_h.y)), vscl(b, local_h.z)));
}

float pdf_ggx(V3 n, V3 wo, V3 wi, float rough) {
    if (vdot(n, wi) <= 0.0f) return 0.0f;
    V3 h = vnorm(vadd(wo, wi));
    float nh = clampf(vdot(n, h), 0.0f, 1.0f);
    float wh = clampf(vdot(wo, h), 0.0f, 1.0f);
    return ggx_d(nh, rough) * nh / (4.0f * wh + 1e-6f);
}
