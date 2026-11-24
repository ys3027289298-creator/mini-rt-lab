
/* nested: path tracer for a small opaque sphere nested inside a glass sphere.
 *
 * Two tracing modes:
 *   buggy - reproduce the legacy behavior: at the outer dielectric surface
 *           the refracted ray performs exactly ONE more intersection query;
 *           if that query leaves the outer medium the ray stops, so paths
 *           that re-enter after total internal reflection are missed and the
 *           silhouette of the inner sphere breaks into blocky chunks.
 *   walk  - continuous ray walking: every interface (entry, exit, TIR bounce,
 *           Lambertian bounce on the inner sphere) becomes the next segment of
 *           the same path until the throughput falls below an energy threshold
 *           or the penetration cap is reached. Both limits are derived from
 *           the inner sphere's average reflectance rho (see derive_limits).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define EPS 1e-4f
#define EVENT_CAP 64          /* hard safety cap for glass-only TIR loops */
#define LSB (1.0f / 255.0f)   /* one quantization step of the 8-bit output */

typedef struct { float x, y, z; } V3;

static V3 vadd(V3 a, V3 b) { V3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static V3 vsub(V3 a, V3 b) { V3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static V3 vscale(V3 a, float s) { V3 r = { a.x * s, a.y * s, a.z * s }; return r; }
static V3 vmul(V3 a, V3 b) { V3 r = { a.x * b.x, a.y * b.y, a.z * b.z }; return r; }
static float vdot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static float vlen(V3 a) { return sqrtf(vdot(a, a)); }
static V3 vnorm(V3 a) {
    float l = vlen(a);
    if (l < 1e-12f) { V3 z = { 0, 0, 0 }; return z; }
    return vscale(a, 1.0f / l);
}
static V3 vcross(V3 a, V3 b) {
    V3 r = { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return r;
}
static float vmaxc(V3 a) {
    float m = a.x > a.y ? a.x : a.y;
    return m > a.z ? m : a.z;
}
static V3 vlerp(V3 a, V3 b, float t) { return vadd(a, vscale(vsub(b, a), t)); }

typedef struct { V3 c; float r; float ior; float rho; int inner; } Sphere;
typedef struct { float t; V3 n; V3 p; int id; } Hit;

/* Scene: glass shell (id 0) containing an opaque Lambertian sphere (id 1). */
static const float GLASS_IOR = 1.5f;
static V3 inner_center = { 0.10f, 0.06f, 0.08f };
static const float INNER_R = 0.55f;
static float inner_rho = 0.25f;

static Sphere outer_sphere(void) {
    V3 zero = { 0, 0, 0 };
    Sphere s = { zero, 1.0f, GLASS_IOR, 0.0f, 0 };
    return s;
}
static Sphere inner_sphere(void) {
    Sphere s = { inner_center, INNER_R, 1.0f, inner_rho, 1 };
    return s;
}

static int sphere_hit(const Sphere* s, V3 o, V3 d, float tmax, Hit* h) {
    V3 oc = vsub(o, s->c);
    float b = vdot(oc, d);
    float c = vdot(oc, oc) - s->r * s->r;
    float disc = b * b - c;
    if (disc < 0.0f) return 0;
    float sq = sqrtf(disc);
    float t = -b - sq;
    if (t < EPS) t = -b + sq;
    if (t < EPS || t >= tmax) return 0;
    h->t = t;
    h->p = vadd(o, vscale(d, t));
    h->n = vscale(vsub(h->p, s->c), 1.0f / s->r);
    h->id = s->inner;
    return 1;
}

/* Nearest of the two spheres; returns -1 on a miss. */
static int nearest(V3 o, V3 d, Hit* h) {
    Sphere so = outer_sphere(), si = inner_sphere();
    Hit h0, h1;
    int a = sphere_hit(&so, o, d, 1e30f, &h0);
    int b = sphere_hit(&si, o, d, 1e30f, &h1);
    if (a && b) { *h = (h0.t < h1.t) ? h0 : h1; }
    else if (a) { *h = h0; }
    else if (b) { *h = h1; }
    else return -1;
    return h->id;
}

static V3 env_color(V3 d) {
    V3 sky = { 0.72f, 0.80f, 0.90f };
    V3 ground = { 0.20f, 0.22f, 0.26f };
    float t = 0.5f * d.y + 0.5f;
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    return vlerp(ground, sky, t);
}

/* R = Fresnel reflectance. eta = n1 / n2 with n1 the incident medium. */
static float fresnel(V3 d, V3 n, float eta) {
    float cosi = -vdot(d, n);
    if (cosi < 0.0f) { cosi = -cosi; eta = 1.0f / eta; }
    float sin2i = 1.0f - cosi * cosi;
    float sin2t = eta * eta * sin2i;
    if (sin2t >= 1.0f) return 1.0f;                 /* total internal reflection */
    float cost = sqrtf(1.0f - sin2t);
    float rs = (eta * cosi - cost) / (eta * cosi + cost);
    float rp = (cosi - eta * cost) / (cosi + eta * cost);
    return 0.5f * (rs * rs + rp * rp);
}

static V3 reflect3(V3 d, V3 n) {
    return vsub(d, vscale(n, 2.0f * vdot(d, n)));
}

/* Refracted direction; returns 0 on TIR. eta = n1 / n2. */
static int refract3(V3 d, V3 n, float eta, V3* td) {
    V3 nf = vdot(d, n) > 0.0f ? vscale(n, -1.0f) : n;
    float cosi = -vdot(d, nf);
    float k = 1.0f - eta * eta * (1.0f - cosi * cosi);
    if (k < 0.0f) return 0;
    *td = vnorm(vadd(vscale(d, eta), vscale(nf, eta * cosi - sqrtf(k))));
    return 1;
}

typedef struct { uint32_t s; } Rng;

static float rnd(Rng* r) {
    r->s = r->s * 1664525u + 1013904223u;
    return (float)(r->s >> 8) * (1.0f / 16777216.0f);
}

/* Cosine-weighted hemisphere sample around normal n. */
static V3 lambert_sample(V3 n, Rng* r) {
    V3 ax = fabsf(n.x) > 0.9f ? (V3){ 0, 1, 0 } : (V3){ 1, 0, 0 };
    V3 t = vnorm(vcross(ax, n));
    V3 b = vcross(n, t);
    float u1 = rnd(r), u2 = rnd(r);
    float rad = sqrtf(u1), phi = 2.0f * (float)M_PI * u2;
    V3 v = vadd(vadd(vscale(t, rad * cosf(phi)), vscale(b, rad * sinf(phi))),
                vscale(n, sqrtf(1.0f - u1)));
    return vnorm(v);
}

/*
 * Termination limits derived from the inner sphere's average reflectance rho.
 *
 * Each penetration of the inner sphere multiplies throughput by rho, so after
 * K penetrations the surviving energy is rho^K. The display cannot resolve a
 * contribution below one 8-bit quantization step (LSB = 1/255). The smallest
 * K satisfying rho^K <= LSB therefore bounds the truncation error below what
 * one output pixel can represent, and the matching energy threshold is the
 * residual energy at that cap. Nothing here is tuned by hand: both numbers
 * are functions of rho and the output bit depth.
 */
static void derive_limits(float rho, int* K, float* tau) {
    int k = 0;
    float e = rho;                        /* energy after the first penetration */
    while (e > LSB) { e *= rho; k++; }
    *K = k + 1;                           /* smallest K with rho^K <= LSB */
    *tau = e;
}

static V3 black(void) { V3 z = { 0, 0, 0 }; return z; }

/* Termination counters: why paths died (evidence for the threshold analysis). */
static long n_capkill = 0, n_taukill = 0, n_inner_events = 0, n_tir = 0;

/* Continuous path walking. *hit_inner flags any path that reached the inner
 * sphere (used to build the silhouette mask for the edge statistics). */
static V3 trace_walk(V3 o, V3 d, Rng* r, int cap, float tau, int* hit_inner) {
    V3 T = { 1, 1, 1 };
    int inner_hits = 0, events = 0;
    for (;;) {
        if (++events > EVENT_CAP) return black();
        Hit h;
        int id = nearest(o, d, &h);
        if (id < 0) return vmul(T, env_color(d));

        if (id == 1) {                                  /* inner Lambertian sphere */
            *hit_inner = 1;
            n_inner_events++;
            if (++inner_hits > cap) { n_capkill++; return black(); } /* penetration cap */
            T = vscale(T, inner_rho);
            if (vmaxc(T) < tau) { n_taukill++; return black(); }     /* energy threshold */
            V3 nd = lambert_sample(h.n, r);
            o = vadd(h.p, vscale(h.n, EPS));
            d = nd;
            continue;
        }

        int entering = vdot(d, h.n) < 0.0f;
        float eta = entering ? 1.0f / GLASS_IOR : GLASS_IOR;
        float F = fresnel(d, h.n, eta);
        V3 nd;
        if (rnd(r) < F) {
            /* single-sample Fresnel: P(reflect)=F, so the carried weight F/F=1 */
            nd = reflect3(d, h.n);
        } else {
            V3 td;
            if (!refract3(d, h.n, eta, &td)) {          /* TIR: deterministic reflect */
                n_tir++;
                nd = reflect3(d, h.n);
            } else {
                /* P(refract)=1-F, weight (1-F)/(1-F)=1; glass keeps throughput */
                nd = td;
            }
        }
        V3 nside = vdot(nd, h.n) > 0.0f ? h.n : vscale(h.n, -1.0f);
        o = vadd(h.p, vscale(nside, EPS));
        d = nd;
        if (vmaxc(T) < tau) return black();
    }
}

/* Legacy behavior: one intersection per refracted ray, stop on exiting glass. */
static V3 trace_buggy(V3 o, V3 d, Rng* r, int* hit_inner) {
    Hit h;
    int id = nearest(o, d, &h);
    if (id < 0) return env_color(d);
    if (id == 1) {
        *hit_inner = 1;
        V3 nd = lambert_sample(h.n, r);
        return vscale(env_color(nd), inner_rho);
    }

    /* First hit is the outer shell entered from vacuum: refract once and run
     * exactly one more intersection query. */
    V3 td;
    if (!refract3(d, h.n, 1.0f / GLASS_IOR, &td)) return env_color(d);
    V3 no = vadd(h.p, vscale(td, EPS));
    Hit h2;
    int id2 = nearest(no, td, &h2);
    if (id2 < 0) return env_color(d);                   /* ray left the scene */
    if (id2 == 0) return env_color(d);                  /* exited the shell: stop */
    *hit_inner = 1;                                     /* direct hit on inner sphere */
    V3 nd = lambert_sample(h2.n, r);
    return vscale(env_color(nd), inner_rho);
}

static uint8_t encode(float c) {
    if (c < 0.0f) c = 0.0f;
    if (c > 1.0f) c = 1.0f;
    c = powf(c, 1.0f / 2.2f);
    int v = (int)(c * 255.0f + 0.5f);
    if (v < 0) v = 0;
    if (v > 255) v = 255;
    return (uint8_t)v;
}

static int write_ppm(const char* path, const float* img, const uint8_t* mask,
                     int W, int H) {
    FILE* f = fopen(path, "wb");
    if (!f) { perror(path); return 1; }
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    for (int i = 0; i < W * H; i++) {
        uint8_t px[3] = { encode(img[3 * i]), encode(img[3 * i + 1]), encode(img[3 * i + 2]) };
        fwrite(px, 1, 3, f);
    }
    fclose(f);

    size_t L = strlen(path);
    char* mpath = malloc(L + 16);
    if (L >= 4 && strcmp(path + L - 4, ".ppm") == 0) {
        memcpy(mpath, path, L - 4);
        strcpy(mpath + L - 4, ".mask.ppm");
    } else {
        strcpy(mpath, path);
        strcat(mpath, ".mask.ppm");
    }
    f = fopen(mpath, "wb");
    if (!f) { perror(mpath); free(mpath); return 1; }
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    for (int i = 0; i < W * H; i++) {
        uint8_t v = mask[i] ? 255 : 0;
        uint8_t px[3] = { v, v, v };
        fwrite(px, 1, 3, f);
    }
    fclose(f);
    free(mpath);
    return 0;
}

static int render(int W, int H, int spp, int buggy, int cap, float tau,
                  uint32_t seed, const char* out) {
    float* img = calloc((size_t)W * H * 3, sizeof(float));
    uint8_t* mask = calloc((size_t)W * H, sizeof(uint8_t));
    if (!img || !mask) { fprintf(stderr, "oom\n"); return 1; }

    n_capkill = n_taukill = n_inner_events = n_tir = 0;
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            V3 acc = { 0, 0, 0 };
            int seen = 0;
            for (int s = 0; s < spp; s++) {
                Rng r = { seed + (uint32_t)((y * W + x) * spp + s) * 2654435761u };
                /* jitter of half a pixel keeps stratified samples inside the pixel */
                float jx = spp > 1 ? rnd(&r) - 0.5f : 0.0f;
                float jy = spp > 1 ? rnd(&r) - 0.5f : 0.0f;
                float fx = ((x + 0.5f + jx) / W - 0.5f) * 2.3f;
                float fy = (0.5f - (y + 0.5f + jy) / H) * 2.3f;
                V3 o = { fx, fy, 3.0f };
                V3 d = { 0, 0, -1 };
                int hit_inner = 0;
                V3 c = buggy ? trace_buggy(o, d, &r, &hit_inner)
                             : trace_walk(o, d, &r, cap, tau, &hit_inner);
                acc = vadd(acc, c);
                if (hit_inner) seen = 1;
            }
            size_t idx = (size_t)(y * W + x) * 3;
            img[idx + 0] = acc.x / spp;
            img[idx + 1] = acc.y / spp;
            img[idx + 2] = acc.z / spp;
            mask[y * W + x] = seen ? 1 : 0;
        }
    }
    printf("events: inner=%ld tir=%ld capkill=%ld taukill=%ld\n",
           n_inner_events, n_tir, n_capkill, n_taukill);
    int rc = write_ppm(out, img, mask, W, H);
    free(img);
    free(mask);
    return rc;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int selftest(void) {
    /* ray through sphere center: entry at t = d - R, exit at t = d + R */
    Sphere s = outer_sphere();
    V3 o = { 0, 0, 3 }, d = { 0, 0, -1 };
    Hit h;
    check("sphere_entry", sphere_hit(&s, o, d, 1e30f, &h) &&
          fabsf(h.t - 2.0f) < 1e-5f && fabsf(h.n.z - 1.0f) < 1e-5f);

    V3 om = { 5, 0, 3 };
    check("sphere_miss", !sphere_hit(&s, om, d, 1e30f, &h));

    /* Fresnel at normal incidence for air->glass: ((eta-1)/(eta+1))^2 */
    V3 dn = { 0, 0, -1 }, nn = { 0, 0, 1 };
    float F0 = fresnel(dn, nn, 1.0f / GLASS_IOR);
    float expect = (GLASS_IOR - 1.0f) / (GLASS_IOR + 1.0f);
    expect *= expect;
    check("fresnel_normal", fabsf(F0 - expect) < 1e-4f);
    V3 dg = vnorm((V3){ 1.0f, 0, -0.005f }); /* ~89.7 deg from the normal */
    check("fresnel_grazing", fresnel(dg, nn, 1.0f / GLASS_IOR) > 0.95f);

    /* TIR for glass->air beyond the critical angle asin(1/1.5) ~ 41.8 deg */
    V3 inside = { 0, 0, 1 };                 /* ray heading outward, normal +z */
    V3 tilt = vnorm((V3){ 0.8f, 0, 0.6f });  /* 53 deg from normal */
    V3 td;
    check("tir_reflects", !refract3(tilt, inside, GLASS_IOR, &td));
    V3 shallow = vnorm((V3){ 0.3f, 0, 0.95f });
    check("subcritical_refracts", refract3(shallow, inside, GLASS_IOR, &td));

    /* derived limits: rho^K <= LSB and rho^(K-1) > LSB */
    int K; float tau;
    derive_limits(0.25f, &K, &tau);
    check("cap_rho25", K == 4 && tau <= LSB &&
          powf(0.25f, (float)(K - 1)) > LSB);
    derive_limits(0.75f, &K, &tau);
    check("cap_rho75", K == 20 && tau <= LSB &&
          powf(0.75f, (float)(K - 1)) > LSB);

    /* RNG determinism and rough uniformity */
    Rng a = { 1234 }, b = { 1234 };
    int same = 1;
    double mean = 0.0;
    int n = 4096;
    for (int i = 0; i < n; i++) {
        float x = rnd(&a), y = rnd(&b);
        if (x != y) same = 0;
        mean += x;
    }
    check("rng_deterministic", same);
    check("rng_uniform", fabs(mean / n - 0.5) < 0.02);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
    if (argc == 2 && strcmp(argv[1], "selftest") == 0) return selftest();
    if (argc != 10) {
        fprintf(stderr,
            "usage: %s W H SPP rho buggy|walk K(-1=derive) tau(-1=derive) seed out.ppm\n",
            argv[0]);
        return 2;
    }
    int W = atoi(argv[1]), H = atoi(argv[2]), spp = atoi(argv[3]);
    inner_rho = (float)atof(argv[4]);
    int buggy = strcmp(argv[5], "buggy") == 0;
    int K = atoi(argv[6]);
    float tau = (float)atof(argv[7]);
    uint32_t seed = (uint32_t)strtoul(argv[8], NULL, 10);
    const char* out = argv[9];

    if (K < 0 || tau < 0.0f) {
        int dK; float dtau;
        derive_limits(inner_rho, &dK, &dtau);
        if (K < 0) K = dK;
        if (tau < 0.0f) tau = dtau;
    }
    printf("rho=%.4f mode=%s K=%d tau=%.6f LSB=%.6f spp=%d seed=%u out=%s\n",
           inner_rho, buggy ? "buggy" : "walk", K, tau, LSB, spp, seed, out);
    return render(W, H, spp, buggy, K, tau, seed, out);
}
