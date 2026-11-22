
/* tile: deterministic tiled rendering.
   Each pixel's RNG stream is derived from (base_seed, x, y) alone, never
   from tile id, tile size, or completion order — so the tiled render is
   bit-identical to the serial one for any partition. Tiles are staged
   separately and merged into the framebuffer in fixed row-major order. */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#define MAX_W 64
#define MAX_H 64
/* staging pool: tiles cover the image, edge tiles included, so the worst
   case is bounded by rounding both dims up one full tile */
#define MAX_STAGE ((2 * MAX_W) * (2 * MAX_H))

typedef struct {
    int width, height;  /* image size in pixels */
    int tile_w, tile_h; /* tile size in pixels */
    int spp;            /* samples per pixel */
    unsigned long long base_seed;
} RenderCfg;

typedef struct { float r, g, b; } Pixel;

static void set_err(char* err, size_t errn, const char* fmt, ...) {
    if (!err || errn == 0) return;
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(err, errn, fmt, ap);
    va_end(ap);
}

/* splitmix64-style finalizer */
static unsigned long long mix64(unsigned long long z) {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* Per-pixel seed: function of (base, x, y) only. This is what makes the
   result independent of tile size and scheduling. */
static unsigned long long pixel_seed(unsigned long long base, int x, int y) {
    return mix64(base ^ mix64((unsigned long long)(unsigned)x * 0x9e3779b97f4a7c15ULL
                            + (unsigned long long)(unsigned)y * 0xbf58476d1ce4e5b9ULL
                            + 0x94d049bb133111ebULL));
}

/* Tile seed: for tile-local bookkeeping (logs, work stealing) only.
   Must never feed pixel samples, or tile size would change the image. */
static unsigned long long tile_seed(unsigned long long base, int tx, int ty) {
    return mix64((base + 0x632be59bd9b4e019ULL)
               ^ mix64((unsigned long long)(unsigned)tx * 0x9e3779b97f4a7c15ULL
                     + (unsigned long long)(unsigned)ty * 0xbf58476d1ce4e5b9ULL));
}

/* Validated seed derivation: illegal parameters are reported, not used. */
static int seed_derive(unsigned long long base, int x, int y, int w, int h,
                       unsigned long long* out, char* err, size_t errn) {
    if (w <= 0 || h <= 0) {
        set_err(err, errn, "seed_derive: image size must be positive, got %dx%d", w, h);
        return 0;
    }
    if (x < 0 || x >= w || y < 0 || y >= h) {
        set_err(err, errn, "seed_derive: pixel (%d,%d) outside %dx%d image", x, y, w, h);
        return 0;
    }
    if (!out) {
        set_err(err, errn, "seed_derive: null output pointer");
        return 0;
    }
    *out = pixel_seed(base, x, y);
    return 1;
}

static int cfg_check(const RenderCfg* c, char* err, size_t errn) {
    if (!c) { set_err(err, errn, "render: null config"); return 0; }
    if (c->width <= 0 || c->height <= 0) {
        set_err(err, errn, "render: image size must be positive, got %dx%d",
                c->width, c->height);
        return 0;
    }
    if (c->width > MAX_W || c->height > MAX_H) {
        set_err(err, errn, "render: image %dx%d exceeds max %dx%d",
                c->width, c->height, MAX_W, MAX_H);
        return 0;
    }
    if (c->tile_w <= 0 || c->tile_h <= 0) {
        set_err(err, errn, "render: tile size must be positive, got %dx%d",
                c->tile_w, c->tile_h);
        return 0;
    }
    if (c->spp <= 0) {
        set_err(err, errn, "render: spp must be positive, got %d", c->spp);
        return 0;
    }
    return 1;
}

typedef struct { unsigned long long s; } Rng;

static unsigned long long rng_next(Rng* r) {
    r->s += 0x9e3779b97f4a7c15ULL;
    return mix64(r->s);
}
static double rng_double(Rng* r) {
    return (double)(rng_next(r) >> 11) * (1.0 / 9007199254740992.0);
}

/* Per-pixel Monte Carlo estimate; sample order within a pixel is fixed. */
static Pixel shade(const RenderCfg* c, int x, int y) {
    Rng rng = { pixel_seed(c->base_seed, x, y) };
    double u0 = (double)x / (double)c->width;
    double v0 = (double)y / (double)c->height;
    double acc = 0.0, acc2 = 0.0;
    for (int s = 0; s < c->spp; s++) {
        double u = rng_double(&rng);
        double v = rng_double(&rng);
        double fu = u0 + u / (double)c->width;
        double fv = v0 + v / (double)c->height;
        acc += fu * fu + fv;
        acc2 += fu * fv;
    }
    Pixel p;
    p.r = (float)(acc / c->spp);
    p.g = (float)(acc2 / c->spp);
    p.b = (float)(u0 * v0);
    return p;
}

static int render_serial(const RenderCfg* c, Pixel* fb, char* err, size_t errn) {
    if (!cfg_check(c, err, errn)) return 0;
    for (int y = 0; y < c->height; y++)
        for (int x = 0; x < c->width; x++)
            fb[y * c->width + x] = shade(c, x, y);
    return 1;
}

/* Tiled render. Tiles are computed in arbitrary completion order
   (reverse when `reverse` is set, simulating out-of-order workers),
   staged per tile, then merged in fixed row-major tile order. */
static int render_tiled(const RenderCfg* c, Pixel* fb, int reverse, char* err, size_t errn) {
    if (!cfg_check(c, err, errn)) return 0;
    int ntx = (c->width + c->tile_w - 1) / c->tile_w;
    int nty = (c->height + c->tile_h - 1) / c->tile_h;
    int ntiles = ntx * nty;
    size_t tile_px = (size_t)c->tile_w * (size_t)c->tile_h;
    if ((size_t)ntiles * tile_px > (size_t)MAX_STAGE) {
        set_err(err, errn, "render: staging overflow for %d tiles of %dx%d",
                ntiles, c->tile_w, c->tile_h);
        return 0;
    }
    static Pixel stage[MAX_STAGE];

    /* parallel phase: completion order may be anything */
    for (int k = 0; k < ntiles; k++) {
        int t = reverse ? ntiles - 1 - k : k;
        int tx = t % ntx, ty = t / ntx;
        Pixel* buf = stage + (size_t)t * tile_px;
        for (int j = 0; j < c->tile_h; j++)
            for (int i = 0; i < c->tile_w; i++) {
                int x = tx * c->tile_w + i, y = ty * c->tile_h + j;
                if (x < c->width && y < c->height)
                    buf[j * c->tile_w + i] = shade(c, x, y);
            }
    }
    /* merge phase: fixed row-major tile order, fixed pixel positions */
    for (int t = 0; t < ntiles; t++) {
        int tx = t % ntx, ty = t / ntx;
        const Pixel* buf = stage + (size_t)t * tile_px;
        for (int j = 0; j < c->tile_h; j++)
            for (int i = 0; i < c->tile_w; i++) {
                int x = tx * c->tile_w + i, y = ty * c->tile_h + j;
                if (x < c->width && y < c->height)
                    fb[y * c->width + x] = buf[j * c->tile_w + i];
            }
    }
    return 1;
}

static int failures = 0;
static void check(const char* name, int ok) {
    printf("%s=%d\n", name, ok ? 1 : 0);
    if (!ok) failures++;
}

static int fb_equal(const Pixel* a, const Pixel* b, int n) {
    return memcmp(a, b, (size_t)n * sizeof(Pixel)) == 0;
}

int main(void) {
    char err[128];
    RenderCfg c = { 37, 29, 8, 8, 16, 12345 }; /* non-multiple of tile size */
    static Pixel serial[MAX_W * MAX_H], tiled[MAX_W * MAX_H], other[MAX_W * MAX_H];
    int npix = c.width * c.height;

    check("serial_ok", render_serial(&c, serial, err, sizeof(err)));
    check("tiled_ok", render_tiled(&c, tiled, 0, err, sizeof(err)));
    check("serial_tiled_identical", fb_equal(serial, tiled, npix));

    /* out-of-order completion still merges to the same image */
    check("tiled_reverse_ok", render_tiled(&c, other, 1, err, sizeof(err)));
    check("order_independent", fb_equal(tiled, other, npix));

    /* tile size must not change a single pixel */
    int sizes[][2] = { {1, 1}, {3, 5}, {7, 16}, {13, 2}, {64, 64} };
    int all_same = 1;
    for (unsigned i = 0; i < sizeof(sizes) / sizeof(sizes[0]); i++) {
        RenderCfg c2 = c;
        c2.tile_w = sizes[i][0];
        c2.tile_h = sizes[i][1];
        if (!render_tiled(&c2, other, i % 2, err, sizeof(err)) ||
            !fb_equal(serial, other, npix))
            all_same = 0;
    }
    check("tile_size_invariant", all_same);

    /* seed derivation: deterministic, per-pixel distinct, tile-independent */
    unsigned long long s1 = pixel_seed(7, 3, 4), s2 = pixel_seed(7, 3, 4);
    check("seed_deterministic", s1 == s2);
    check("seed_distinct_xy", pixel_seed(7, 3, 4) != pixel_seed(7, 4, 3) &&
                              pixel_seed(7, 3, 4) != pixel_seed(7, 3, 5));
    check("seed_distinct_base", pixel_seed(7, 3, 4) != pixel_seed(8, 3, 4));
    check("tile_seed_distinct", tile_seed(7, 0, 0) != tile_seed(7, 1, 0) &&
                                tile_seed(7, 0, 0) != pixel_seed(7, 0, 0));

    /* illegal seed-derivation parameters are reported */
    unsigned long long out;
    check("seed_reject_neg", !seed_derive(1, -1, 0, 8, 8, &out, err, sizeof(err)) && err[0]);
    check("seed_reject_oob", !seed_derive(1, 8, 0, 8, 8, &out, err, sizeof(err)) && err[0]);
    check("seed_reject_dim", !seed_derive(1, 0, 0, 0, 8, &out, err, sizeof(err)) && err[0]);
    check("seed_accept", seed_derive(1, 7, 7, 8, 8, &out, err, sizeof(err)) &&
                         out == pixel_seed(1, 7, 7));

    /* illegal render configs are reported */
    RenderCfg bad = c;
    bad.width = 0;
    check("cfg_reject_width", !render_serial(&bad, other, err, sizeof(err)) && err[0]);
    bad = c; bad.tile_w = -2;
    check("cfg_reject_tile", !render_tiled(&bad, other, 0, err, sizeof(err)) && err[0]);
    bad = c; bad.spp = 0;
    check("cfg_reject_spp", !render_serial(&bad, other, err, sizeof(err)) && err[0]);
    check("cfg_reject_null", !render_serial(NULL, other, err, sizeof(err)) && err[0]);

    printf("all_ok=%d\n", failures == 0);
    return failures == 0 ? 0 : 1;
}
