/* glass.h: facade moire model, sample-position generator, fringe measurement. */
#ifndef GLASS_H
#define GLASS_H

#define GLASS_PI 3.14159265358979323846

/* grid: fixed stratified lattice, identical in every pixel.
   offset: same lattice displaced by a low-discrepancy R2 offset. */
typedef enum { GLASS_SAMPLE_GRID = 0, GLASS_SAMPLE_OFFSET = 1 } GlassSampleMode;

/* Facade tuning. Kept as macros so the wall and the analysis cannot drift. */
#define GLASS_RIPPLE          0.28   /* ripple amplitude, fraction of brightness */
#define GLASS_RIPPLE_CYCLES   1.0    /* ripple cycles across one panel */
#define GLASS_MULLION_DEPTH   0.80   /* mullion darkening */
#define GLASS_MULLION_DUTY    0.08   /* mullion width, fraction of panel pitch */
#define GLASS_GRAD_X          0.18   /* horizontal illumination gradient */
#define GLASS_GRAD_Y          0.10   /* vertical illumination gradient */

#define GLASS_MAX_SPP 8

typedef struct {
    int width, height;   /* image size in pixels */
    int subdiv;          /* D: panels across the image width */
    int spp;             /* S: stratified samples per pixel per axis (S*S total) */
    double ripple;       /* ripple amplitude (defaults to GLASS_RIPPLE) */
    double mullion;      /* mullion darkening (defaults to GLASS_MULLION_DEPTH) */
    double duty;         /* mullion width fraction (defaults to GLASS_MULLION_DUTY) */
    double grad_x, grad_y;
} GlassWall;

GlassWall glass_wall(int width, int height, int subdiv, int spp);

/* p = W/D, the panel pitch in pixels. */
double glass_panel_pitch(const GlassWall* w);

/* Shading at a continuous pixel-space position. */
double glass_shade(double px, double py, const GlassWall* w);

/* The offset range derived from the subdivision density, in pixels. */
double glass_offset_radius(const GlassWall* w);

/* Predicted grid-sampler beat wavelength lambda = p / frac(S*p) in pixels,
   or 0 when the lattice is commensurate with the panels (no beat). */
double glass_predicted_stripe(const GlassWall* w);

/* Fill out[][2] with S*S pixel-local sample positions; returns the count. */
int glass_pixel_samples(const GlassWall* w, int px, int py, GlassSampleMode mode,
                        double out[][2]);

/* Render the whole facade (w*h doubles, row-major) with the chosen sampler. */
void glass_render(const GlassWall* w, GlassSampleMode mode, double* out);

/* Exact (analytically integrated) box-filtered facade: the reference. */
void glass_reference(const GlassWall* w, double* out);

/* Dominant stripe spacing (px) between a render and the reference; writes
   DFT peak-over-median prominence and the detrended-profile RMS amplitude.
   Returns 0 if no stripe is measurable. */
double glass_fringe_spacing(const double* img, const double* ref,
                            const GlassWall* w, double* prominence,
                            double* amplitude);

/* residual RMS after averaging each column: how strong the coherent stripe is */
double glass_band_rms(const double* img, const double* ref, const GlassWall* w);

/* residual RMS per pixel: ordinary sampling error, stripe or not */
double glass_pixel_rms(const double* img, const double* ref, const GlassWall* w);

int glass_selftest(void);

#endif
