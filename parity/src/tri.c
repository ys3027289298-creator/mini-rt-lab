#include "tri.h"
#include <string.h>
#include <float.h>

static void aabb_empty(AABB* b) {
    for (int i = 0; i < 3; i++) { b->min[i] = FLT_MAX; b->max[i] = -FLT_MAX; }
}
static void aabb_grow_p(AABB* b, V3 p) {
    if (p.x < b->min[0]) b->min[0] = p.x;
    if (p.x > b->max[0]) b->max[0] = p.x;
    if (p.y < b->min[1]) b->min[1] = p.y;
    if (p.y > b->max[1]) b->max[1] = p.y;
    if (p.z < b->min[2]) b->min[2] = p.z;
    if (p.z > b->max[2]) b->max[2] = p.z;
}
static void aabb_grow(AABB* b, const AABB* o) {
    for (int i = 0; i < 3; i++) {
        if (o->min[i] < b->min[i]) b->min[i] = o->min[i];
        if (o->max[i] > b->max[i]) b->max[i] = o->max[i];
    }
}
/* Robust slab test. A zero direction component means the ray is parallel
   to that slab: it misses iff its origin is outside. Computing 1/d up front
   turns this case into inf arithmetic that wrongly collapses t1 to 0 for
   zero-thickness boxes (the ground plane), so direction is tested directly. */
static int aabb_hit(const AABB* b, Ray ray, float tmax) {
    float t0 = 0.0f, t1 = tmax;
    float o[3] = { ray.o.x, ray.o.y, ray.o.z };
    float d[3] = { ray.d.x, ray.d.y, ray.d.z };
    for (int i = 0; i < 3; i++) {
        float lo = b->min[i], hi = b->max[i];
        if (fabsf(d[i]) < 1e-10f) {
            if (o[i] < lo || o[i] > hi) return 0;
        } else {
            float inv = 1.0f / d[i];
            float ta = (lo - o[i]) * inv, tb = (hi - o[i]) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            if (ta > t0) t0 = ta;
            if (tb < t1) t1 = tb;
            if (t0 > t1) return 0;
        }
    }
    return 1;
}

static void tri_bounds(const Scene* sc, int ti, AABB* out) {
    aabb_empty(out);
    const Tri* t = &sc->tri[ti];
    for (int k = 0; k < 3; k++) aabb_grow_p(out, t->v[k]);
}

static V3 tri_centroid(const Tri* t) {
    return vscl(vadd(vadd(t->v[0], t->v[1]), t->v[2]), 1.0f / 3.0f);
}

static int build_rec(Bvh* b, int first, int count) {
    int ni = b->nnode++;
    BvhNode* n = &b->node[ni];
    aabb_empty(&n->box);
    for (int i = first; i < first + count; i++) {
        AABB tb; tri_bounds(b->sc, b->order[i], &tb);
        aabb_grow(&n->box, &tb);
    }
    if (count <= 2) {
        n->left = n->right = -1;
        n->first = first; n->count = count;
        return ni;
    }
    float w[3] = { n->box.max[0] - n->box.min[0],
                   n->box.max[1] - n->box.min[1],
                   n->box.max[2] - n->box.min[2] };
    int axis = (w[0] > w[1]) ? (w[0] > w[2] ? 0 : 2) : (w[1] > w[2] ? 1 : 2);
    float mid = 0.5f * (n->box.min[axis] + n->box.max[axis]);
    int p = first;
    for (int i = first; i < first + count; i++) {
        V3 c = tri_centroid(&b->sc->tri[b->order[i]]);
        float cv = axis == 0 ? c.x : (axis == 1 ? c.y : c.z);
        if (cv < mid) {
            int tmp = b->order[i]; b->order[i] = b->order[p]; b->order[p] = tmp;
            p++;
        }
    }
    if (p == first || p == first + count) p = first + count / 2;
    n->first = -1; n->count = 0;
    n->left = build_rec(b, first, p - first);
    n->right = build_rec(b, p, first + count - p);
    return ni;
}

void bvh_build(Bvh* b, const Scene* sc) {
    b->sc = sc; b->nnode = 0;
    for (int i = 0; i < sc->ntri; i++) b->order[i] = i;
    if (sc->ntri > 0) build_rec(b, 0, sc->ntri);
}

/* Moller-Trumbore; returns t > eps or 0 on miss. Fills shading data on hit. */
static int tri_hit(const Scene* sc, int ti, Ray ray, float tmax, Hit* hit) {
    const Tri* t = &sc->tri[ti];
    V3 e1 = vsub(t->v[1], t->v[0]);
    V3 e2 = vsub(t->v[2], t->v[0]);
    V3 h = vcross(ray.d, e2);
    float a = vdot(e1, h);
    if (fabsf(a) < 1e-9f) return 0;
    float f = 1.0f / a;
    V3 s = vsub(ray.o, t->v[0]);
    float u = f * vdot(s, h);
    if (u < 0.0f || u > 1.0f) return 0;
    V3 q = vcross(s, e1);
    float v = f * vdot(ray.d, q);
    if (v < 0.0f || u + v > 1.0f) return 0;
    float tt = f * vdot(e2, q);
    if (tt <= 1e-4f || tt >= tmax) return 0;
    if (hit) {
        float w = 1.0f - u - v;
        hit->t = tt;
        hit->p = vadd(ray.o, vscl(ray.d, tt));
        V3 gn = vnorm(vcross(e1, e2));
        V3 sn = vnorm(vadd(vadd(vscl(t->n[0], w), vscl(t->n[1], u)), vscl(t->n[2], v)));
        hit->n = vlen2(sn) > 1e-12f ? sn : gn;
        if (vdot(hit->n, ray.d) > 0.0f) hit->n = vscl(hit->n, -1.0f);
        hit->u = w * t->uv[0][0] + u * t->uv[1][0] + v * t->uv[2][0];
        hit->v = w * t->uv[0][1] + u * t->uv[1][1] + v * t->uv[2][1];
        hit->mat = t->mat;
    }
    return 1;
}

int bvh_intersect(const Bvh* b, Ray ray, float tmax, Hit* hit) {
    int stack[BVH_STACK]; int sp = 0;
    stack[sp++] = 0;
    float best = tmax; int found = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        const BvhNode* n = &b->node[ni];
        if (!aabb_hit(&n->box, ray, best)) continue;
        if (n->count > 0) {
            for (int i = n->first; i < n->first + n->count; i++) {
                Hit cand;
                if (tri_hit(b->sc, b->order[i], ray, best, &cand)) {
                    best = cand.t; *hit = cand; found = 1;
                }
            }
        } else {
            stack[sp++] = n->left;
            stack[sp++] = n->right;
        }
    }
    return found;
}

int bvh_occluded(const Bvh* b, Ray ray, float tmax) {
    int stack[BVH_STACK]; int sp = 0;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        const BvhNode* n = &b->node[ni];
        if (!aabb_hit(&n->box, ray, tmax)) continue;
        if (n->count > 0) {
            for (int i = n->first; i < n->first + n->count; i++)
                if (tri_hit(b->sc, b->order[i], ray, tmax, NULL)) return 1;
        } else {
            stack[sp++] = n->left;
            stack[sp++] = n->right;
        }
    }
    return 0;
}

static int sphere_hit(const Sphere* s, Ray ray, float tmax, Hit* hit) {
    V3 oc = vsub(ray.o, s->c);
    float b = vdot(oc, ray.d);
    float c = vdot(oc, oc) - s->r * s->r;
    float disc = b * b - c;
    if (disc < 0.0f) return 0;
    float sq = sqrtf(disc);
    float t = -b - sq;
    if (t <= 1e-4f || t >= tmax) {
        t = -b + sq;
        if (t <= 1e-4f || t >= tmax) return 0;
    }
    if (hit) {
        hit->t = t;
        hit->p = vadd(ray.o, vscl(ray.d, t));
        hit->n = vnorm(vscl(vsub(hit->p, s->c), 1.0f / s->r));
        if (vdot(hit->n, ray.d) > 0.0f) hit->n = vscl(hit->n, -1.0f);
        hit->mat = s->mat;
        hit->u = 0.5f + atan2f(hit->n.z, hit->n.x) / (2.0f * PI);
        hit->v = 0.5f - asinf(clampf(hit->n.y, -1.0f, 1.0f)) / PI;
    }
    return 1;
}

int scene_intersect(const Scene* sc, const Bvh* b, Ray ray, float tmax, Hit* hit) {
    Hit best; float bt = tmax; int found = 0;
    if (sc->ntri > 0 && bvh_intersect(b, ray, bt, &best)) { bt = best.t; found = 1; }
    for (int i = 0; i < sc->nsph; i++) {
        Hit cand;
        if (sphere_hit(&sc->sph[i], ray, bt, &cand)) { bt = cand.t; best = cand; found = 1; }
    }
    if (found) *hit = best;
    return found;
}

int scene_occluded(const Scene* sc, const Bvh* b, Ray ray, float tmax) {
    if (sc->ntri > 0 && bvh_occluded(b, ray, tmax)) return 1;
    for (int i = 0; i < sc->nsph; i++)
        if (sphere_hit(&sc->sph[i], ray, tmax, NULL)) return 1;
    return 0;
}
