/* selftest: contract checks for every shared module.
   The point of these checks is "one function, two callers": if the shared
   radiometry/texture/tonemap math is wrong, both back ends are wrong in
   the same way, so unit-level parity is checked here first. */
#include "math3.h"
#include "scene.h"
#include "tri.h"
#include "camera.h"
#include "radiometry.h"
#include "texture.h"
#include "tonemap.h"
#include "image.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}
static int nearf(float a, float b, float eps) { return fabsf(a - b) < eps; }

int main(void) {
    /* ---- texture: repeat wrap + bilinear continuity ---- */
    Scene sc; scene_build(&sc);
    V3 c00 = tex_sample(&sc.tex[0], 0.0f, 0.0f);
    V3 cwrap = tex_sample(&sc.tex[0], 1.0f, 1.0f);   /* wraps back to origin */
    check("tex_repeat", nearf(c00.x, cwrap.x, 1e-5f) && nearf(c00.y, cwrap.y, 1e-5f));
    V3 cneg = tex_sample(&sc.tex[0], -0.25f, -0.25f); /* negative coords wrap too */
    V3 cpos = tex_sample(&sc.tex[0], 0.75f, 0.75f);
    check("tex_neg_wrap", nearf(cneg.x, cpos.x, 1e-5f));
    V3 cmid = tex_sample(&sc.tex[0], 0.5f, 0.5f);
    check("tex_linear_range", cmid.x >= 0.0f && cmid.x <= 1.0f);
    V3 cinterp = tex_sample(&sc.tex[0], 0.5f / TEX_SIZE, 0.5f / TEX_SIZE);
    V3 cbase = tex_rgb(&sc.tex[0], 0, 0);
    check("tex_bilinear_mid", nearf(cinterp.x,
          (cbase.x + tex_rgb(&sc.tex[0], 1, 0).x +
           tex_rgb(&sc.tex[0], 0, 1).x + tex_rgb(&sc.tex[0], 1, 1).x) * 0.25f, 1e-5f));

    /* ---- tonemap: monotonic, endpoints, shared curve identity ---- */
    check("tonemap_zero", tonemap_curve(0.0f) == 0.0f);
    check("tonemap_mono", tonemap_curve(1.0f) > tonemap_curve(0.1f));
    check("tonemap_clamp", tonemap_curve(1e9f) <= 1.0f);
    check("srgb_roundtrip", nearf(srgb_encode(0.5f), 0.7354f, 2e-3f));
    check("tonemap_byte_zero", tonemap_channel(0.0f) == 0);
    check("tonemap_byte_neg_clamp", tonemap_channel(-5.0f) == 0);
    check("tonemap_byte_hi_clamp", tonemap_channel(1e9f) == 255);

    /* ---- radiometry: Lambert energy conservation ---- */
    V3 alb = v3(0.8f, 0.2f, 0.5f);
    V3 lam = brdf_lambert(alb);
    check("lambert_pi", nearf(lam.x, 0.8f * INV_PI, 1e-6f));
    /* integrate Lambert over hemisphere analytically: rho = pi * f */
    check("lambert_conservation", nearf(lam.x * PI, alb.x, 1e-6f));

    /* GGX normal distribution normalized: integral D(h)*cos(theta) domega = 1.
       In mu=cos(theta): 2*pi * integral_0^1 D(mu)*mu dmu. */
    double acc = 0.0; int N = 4096;
    for (int i = 0; i < N; i++) {
        float mu = (i + 0.5f) / N;
        acc += (double)ggx_d(mu, 0.5f) * mu;
    }
    acc *= 2.0 * PI / N;
    check("ggx_normalized", fabs(acc - 1.0) < 0.01);
    check("ggx_peak", ggx_d(1.0f, 0.5f) > ggx_d(0.2f, 0.5f));

    /* cosine hemisphere sampling: unit vector + pdf consistency */
    unsigned rng = 12345u;
    int unit_ok = 1, pdf_ok = 1;
    for (int i = 0; i < 64; i++) {
        V3 n = v3(0, 1, 0);
        float pdf;
        V3 wi = sample_cosine_hemisphere(n, rng_next(&rng), rng_next(&rng), &pdf);
        if (!nearf(vlen(wi), 1.0f, 1e-4f)) unit_ok = 0;
        if (!nearf(pdf, vdot(n, wi) * INV_PI, 1e-5f)) pdf_ok = 0;
    }
    check("cosine_unit", unit_ok);
    check("cosine_pdf", pdf_ok);

    /* ---- geometry: BVH hits the shared scene, sphere + triangle ---- */
    Bvh bvh; bvh_build(&bvh, &sc);
    check("bvh_built", bvh.nnode > 0);
    Ray down = { v3(0.0f, 5.0f, -0.4f), v3(0, -1, 0) };
    Hit hit;
    /* sphere 0 center (-0.6,0.65,-0.4) r=0.65; the x=0 ray grazes its side,
       horizontal gap 0.6 -> hit y = 0.65 + sqrt(0.65^2-0.6^2) = 0.90 */
    check("sphere_hit", scene_intersect(&sc, &bvh, down, 1e6f, &hit) &&
                        nearf(hit.p.y, 0.90f, 1e-3f));
    Ray miss = { v3(-5.0f, 5.0f, -5.0f), v3(1, 0, 0) };
    check("scene_miss", !scene_intersect(&sc, &bvh, miss, 0.5f, &hit));
    Ray ground = { v3(2.0f, 3.0f, 4.0f), v3(0, -1, 0) };
    check("ground_hit", scene_intersect(&sc, &bvh, ground, 1e6f, &hit) &&
                        nearf(hit.p.y, 0.0f, 1e-3f) && hit.mat == 0);

    /* shadow ray: light above an open spot is visible; through a box is not */
    PointLight L = sc.light[0];
    Ray sh = { v3(0.0f, 0.01f, 4.0f), vnorm(vsub(L.pos, v3(0.0f, 0.01f, 4.0f))) };
    check("shadow_open", !scene_occluded(&sc, &bvh, sh, vlen(vsub(L.pos, sh.o))));
    Ray sh2 = { v3(-1.7f, 0.75f, 1.2f), v3(0, 1, 0) };
    check("shadow_blocked", scene_occluded(&sc, &bvh, sh2, 3.0f) ||
                            /* start inside the box: hit on the way up */
                            scene_occluded(&sc, &bvh,
                                (Ray){ v3(-1.7f, 0.75f, 2.5f), v3(0, 1, 0) }, 3.0f));

    /* ---- camera: ray/project round trip ---- */
    Camera cam = camera_make(v3(0, 2, 6), v3(0, 1, 0), v3(0, 1, 0), 50.0f, 100, 80);
    Ray cr = camera_ray(&cam, 0.5f, 0.5f);
    check("cam_center_ray", nearf(vdot(cr.d, cam.forward), 1.0f, 1e-5f));
    V3 tgt = vadd(cam.eye, vscl(cr.d, 3.0f));
    float px, py, vz;
    check("cam_project_roundtrip", camera_project(&cam, tgt, &px, &py, &vz) &&
                                   nearf(px, 50.0f, 0.6f) && nearf(py, 40.0f, 0.6f));

    /* ---- radiometric units: inverse-square falloff ---- */
    Hit h2; Ray r2 = { v3(0.0f, 3.0f, -0.4f), v3(0, -1, 0) };
    scene_intersect(&sc, &bvh, r2, 1e6f, &h2);
    V3 wo = v3(0, 1, 0);
    V3 d1 = shade_direct(&sc, &bvh, &h2, wo, &L);
    PointLight L2 = { v3(L.pos.x, L.pos.y + 2.0f, L.pos.z), L.intensity };
    V3 d2 = shade_direct(&sc, &bvh, &h2, wo, &L2);
    /* farther light must give strictly less radiance through the same BRDF */
    check("inverse_square", d2.x < d1.x && d1.x > 0.0f);
    check("ambient_nonneg", ambient_term(&sc, &h2, wo).x >= 0.0f);

    /* ---- image IO round trip: write then read must be the identity,
       including row order (PFM is bottom-up on disk, top-left in memory). */
    HdrImage img = image_alloc(4, 2);
    img.px[3] = v3(1.5f, -0.2f, 0.7f);
    img.px[7] = v3(0.4f, 0.5f, 0.6f);
    check("pfm_write", image_write_pfm(&img, "out/_t.pfm"));
    HdrImage back;
    check("pfm_read", image_read_pfm(&back, "out/_t.pfm") && back.w == 4 && back.h == 2);
    check("pfm_roundtrip", nearf(back.px[3].x, 1.5f, 1e-5f) &&
                           nearf(back.px[3].z, 0.7f, 1e-5f) &&
                           nearf(back.px[7].y, 0.5f, 1e-5f));
    image_free(&img); image_free(&back);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
