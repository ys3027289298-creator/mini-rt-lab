#include "scene.h"
#include <string.h>

static int add_material(Scene* s, Material m) {
    int id = s->nmat++;
    s->mat[id] = m;
    return id;
}

static void add_quad(Scene* s, V3 a, V3 b, V3 c, V3 d, V3 n, int mat) {
    float (*u0)[2] = s->tri[s->ntri].uv; (void)u0;
    Tri t1 = {
        .v = { a, b, c }, .n = { n, n, n }, .mat = mat,
        .uv = { {0, 0}, {1, 0}, {1, 1} }
    };
    Tri t2 = {
        .v = { a, c, d }, .n = { n, n, n }, .mat = mat,
        .uv = { {0, 0}, {1, 1}, {0, 1} }
    };
    s->tri[s->ntri++] = t1;
    s->tri[s->ntri++] = t2;
}

/* axis-aligned box centered at c with half-size h, uniform material */
static void add_box(Scene* s, V3 c, V3 h, int mat) {
    V3 p0 = vsub(c, h), p1 = vadd(c, h);
    V3 x0 = v3(p0.x, 0, 0), x1 = v3(p1.x, 0, 0);
    V3 y0 = v3(0, p0.y, 0), y1 = v3(0, p1.y, 0);
    V3 z0 = v3(0, 0, p0.z), z1 = v3(0, 0, p1.z);
    V3 v000 = vadd(vadd(x0, y0), z0), v100 = vadd(vadd(x1, y0), z0);
    V3 v110 = vadd(vadd(x1, y1), z0), v010 = vadd(vadd(x0, y1), z0);
    V3 v001 = vadd(vadd(x0, y0), z1), v101 = vadd(vadd(x1, y0), z1);
    V3 v111 = vadd(vadd(x1, y1), z1), v011 = vadd(vadd(x0, y1), z1);
    add_quad(s, v000, v100, v110, v010, v3(0, 0, -1), mat); /* front -z */
    add_quad(s, v101, v001, v011, v111, v3(0, 0, 1), mat);  /* back +z */
    add_quad(s, v001, v000, v010, v011, v3(-1, 0, 0), mat); /* left -x */
    add_quad(s, v100, v101, v111, v110, v3(1, 0, 0), mat);  /* right +x */
    add_quad(s, v010, v110, v111, v011, v3(0, 1, 0), mat);  /* top +y */
    add_quad(s, v001, v101, v100, v000, v3(0, -1, 0), mat); /* bottom -y */
}

int scene_tessellate(Scene* s, int sph_idx, int rings, int sectors) {
    Sphere sp = s->sph[sph_idx];
    for (int r = 0; r < rings; r++) {
        float theta0 = PI * (float)r / rings;
        float theta1 = PI * (float)(r + 1) / rings;
        for (int q = 0; q < sectors; q++) {
            float phi0 = 2.0f * PI * (float)q / sectors;
            float phi1 = 2.0f * PI * (float)(q + 1) / sectors;
            V3 p00 = vadd(sp.c, vscl(v3(sinf(theta0) * cosf(phi0), cosf(theta0), sinf(theta0) * sinf(phi0)), sp.r));
            V3 p01 = vadd(sp.c, vscl(v3(sinf(theta0) * cosf(phi1), cosf(theta0), sinf(theta0) * sinf(phi1)), sp.r));
            V3 p10 = vadd(sp.c, vscl(v3(sinf(theta1) * cosf(phi0), cosf(theta1), sinf(theta1) * sinf(phi0)), sp.r));
            V3 p11 = vadd(sp.c, vscl(v3(sinf(theta1) * cosf(phi1), cosf(theta1), sinf(theta1) * sinf(phi1)), sp.r));
            V3 n00 = vnorm(vsub(p00, sp.c)), n01 = vnorm(vsub(p01, sp.c));
            V3 n10 = vnorm(vsub(p10, sp.c)), n11 = vnorm(vsub(p11, sp.c));
            float v0 = (float)r / rings, v1 = (float)(r + 1) / rings;
            float u0 = (float)q / sectors, u1 = (float)(q + 1) / sectors;
            Tri t1 = {
                .v = { p00, p10, p11 }, .n = { n00, n10, n11 }, .mat = sp.mat,
                .uv = { {u0, v0}, {u0, v1}, {u1, v1} }
            };
            Tri t2 = {
                .v = { p00, p11, p01 }, .n = { n00, n11, n01 }, .mat = sp.mat,
                .uv = { {u0, v0}, {u1, v1}, {u1, v0} }
            };
            if (s->ntri + 2 > MAX_TRI) return -1;
            s->tri[s->ntri++] = t1;
            s->tri[s->ntri++] = t2;
        }
    }
    return 0;
}

void scene_build(Scene* s) {
    memset(s, 0, sizeof(*s));

    /* textures authored in sRGB, converted to linear inside tex_checker */
    tex_checker(&s->tex[0], 8, v3(0.85f, 0.85f, 0.80f), v3(0.25f, 0.28f, 0.32f));
    s->ntex = 1;

    Material floor_m = { MAT_LAMBERT, v3(1, 1, 1), 1.0f, 0 };
    Material cube_l  = { MAT_LAMBERT, v3(0.75f, 0.25f, 0.18f), 1.0f, -1 };
    Material cube_r  = { MAT_GGX,     v3(0.80f, 0.78f, 0.70f), 0.35f, -1 };
    Material sph_l   = { MAT_GGX,     v3(0.18f, 0.45f, 0.80f), 0.15f, -1 };
    Material sph_m   = { MAT_GGX,     v3(0.85f, 0.65f, 0.20f), 0.55f, -1 };
    Material sph_r   = { MAT_LAMBERT, v3(0.30f, 0.70f, 0.35f), 1.0f, -1 };
    add_material(s, floor_m);   /* 0 */
    add_material(s, cube_l);    /* 1 */
    add_material(s, cube_r);    /* 2 */
    add_material(s, sph_l);     /* 3 */
    add_material(s, sph_m);     /* 4 */
    add_material(s, sph_r);     /* 5 */

    /* ground: 12x12 plane at y=0, uv 0..1 per 2 world units */
    V3 gn = v3(0, 1, 0);
    /* ground kept fully in front of the camera (eye z=6.2 looking toward -z),
       so no ground triangle straddles the near plane. */
    V3 ga = v3(-6, 0, -3), gb = v3(6, 0, -3), gc = v3(6, 0, 5.5f), gd = v3(-6, 0, 5.5f);
    add_quad(s, ga, gb, gc, gd, gn, 0);
    /* remap ground uv to world scale after add_quad's 0..1 defaults */
    for (int i = 0; i < 2; i++) {
        Tri* t = &s->tri[i];
        for (int k = 0; k < 3; k++) {
            t->uv[k][0] = (t->v[k].x + 6.0f) * 0.5f;
            t->uv[k][1] = (t->v[k].z + 3.0f) * 0.5f;
        }
    }

    add_box(s, v3(-1.7f, 0.75f, 1.2f), v3(0.75f, 0.75f, 0.75f), 1);
    add_box(s, v3(1.9f, 0.55f, 0.6f), v3(0.55f, 0.55f, 0.55f), 2);

    s->sph[s->nsph++] = (Sphere){ v3(-0.6f, 0.65f, -0.4f), 0.65f, 3 };
    s->sph[s->nsph++] = (Sphere){ v3(0.9f, 0.5f, 2.1f), 0.5f, 4 };
    s->sph[s->nsph++] = (Sphere){ v3(0.4f, 0.35f, -1.4f), 0.35f, 5 };

    /* point light: radiant intensity W/sr. One delta light, so both
       back ends use direct light sampling only (no MIS weight needed). */
    s->light[s->nlight++] = (PointLight){ v3(2.2f, 4.2f, 1.5f),
                                          v3(60.0f, 56.0f, 48.0f) };

    /* environment radiance: flat ambient term, same value for both */
    s->env = v3(0.12f, 0.14f, 0.18f);
}
