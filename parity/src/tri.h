/* tri: median-split BVH over scene triangles, stack traversal, sphere hits.
   Shared by the path tracer (primary + shadow + bounce rays) and the
   rasterizer (shadow rays), so occlusion decisions cannot diverge. */
#ifndef PARITY_TRI_H
#define PARITY_TRI_H

#include "scene.h"

#define BVH_MAX_NODE 16384
#define BVH_STACK 64

typedef struct { float min[3], max[3]; } AABB;
typedef struct { AABB box; int left, right; int first, count; } BvhNode;

typedef struct {
    const Scene* sc;
    BvhNode node[BVH_MAX_NODE];
    int order[MAX_TRI];
    int nnode;
} Bvh;

void bvh_build(Bvh* b, const Scene* sc);
int bvh_intersect(const Bvh* b, Ray ray, float tmax, Hit* hit);
int bvh_occluded(const Bvh* b, Ray ray, float tmax);

/* analytic spheres, checked alongside triangles */
int scene_intersect(const Scene* sc, const Bvh* b, Ray ray, float tmax, Hit* hit);
int scene_occluded(const Scene* sc, const Bvh* b, Ray ray, float tmax);

#endif
