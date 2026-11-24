/* scene: one shared scene description consumed by both renderers.
   All light transport values are radiometric:
     - point light: radiant intensity I [W/sr], color-tinted
     - surfaces:    BRDF [1/sr], linear albedo
     - environment: radiance Le [W/(sr*m^2)]
   No arbitrary "intensity" knobs live anywhere else. */
#ifndef PARITY_SCENE_H
#define PARITY_SCENE_H

#include "math3.h"
#include "texture.h"

#define MAX_TRI 8192
#define MAX_SPHERE 8
#define MAX_LIGHT 4
#define MAX_MAT 16
#define MAX_TEX 4

typedef enum { MAT_LAMBERT, MAT_GGX } MatType;

typedef struct {
    MatType type;
    V3 base;            /* linear albedo */
    float rough;        /* GGX linear roughness alpha */
    int tex;            /* texture index, -1 = use base only */
} Material;

typedef struct {
    V3 v[3], n[3];
    float uv[3][2];
    int mat;
} Tri;

typedef struct { V3 c; float r; int mat; } Sphere;
typedef struct { V3 pos; V3 intensity; } PointLight;

typedef struct {
    Tri tri[MAX_TRI];
    Sphere sph[MAX_SPHERE];
    Material mat[MAX_MAT];
    Texture tex[MAX_TEX];
    PointLight light[MAX_LIGHT];
    V3 env;
    int ntri, nsph, nmat, ntex, nlight;
} Scene;

typedef struct {
    float t;
    V3 p;               /* hit point */
    V3 n;               /* shading normal, oriented against the ray */
    float u, v;         /* surface uv */
    int mat;
} Hit;

void scene_build(Scene* s);

/* tessellate an analytic sphere into triangles (rasterizer representation
   of the same primitive; the path tracer intersects the sphere directly). */
int scene_tessellate(Scene* s, int sph_idx, int rings, int sectors);

#endif
