# shadowmap — where a point-light shadow map spends its texels

A point light at height `H` above a wall, one hard-edged occluder field, and a
polar shadow map whose *radial row budget* is re-allocated with distance from
the light. The question: given the "square patch ring near the light" symptom,
which allocation exponent is right, and is it larger or smaller than the one a
solid-angle derivation predicts?

## Build and run

```sh
make            # builds ./shadowmap and runs the self-test
make run        # writes out/*.csv and the wall captures = out/*.ppm
make figures    # out/*.csv -> figs/*.png   (needs python3 + matplotlib)
```

Everything numerical comes from `src/shadowmap.c` (plain C99, no deps).
`./shadowmap --selftest` checks the allocation law, the truth intersection and
the map round-trip before any measurement is trusted. Output files land next to
the CSV path you pass, so `./shadowmap out/points.csv 0` writes `out/summary.csv`
etc. A third argument of `1` also dumps wall and map images as PPM.

## Setup

* wall `y = 0`, albedo 1; light `(0, 1, 0)`; orthographic camera from `+y`,
  so one image pixel is one fixed wall area (no perspective term in the metric).
* 14 small blockers hang between light and wall at `y ≈ 0.55`; their umbrae
  land at target radii spread over `0.2 … 2.2 H`. Visibility is analytic
  (ray–sphere), so the map is compared against exact truth, not a reference map.
* the layout is re-drawn from a fixed seed for each of 12 trials; every
  per-radius number below is a mean over those trials.
* shadow map: `NH = 256` rows × `NP = 512` columns, azimuth always uniform.
  `NH·NP` is held fixed for every allocation, so runs differ only in *where the
  rows go*. A lookup reads the nearest texel's binary visibility — the error a
  real map has before filtering.

Twenty sample radii per run, geometric from 0.06 H to 2.40 H.

## The derivation

Solid angle of a wall element at in-plane radius `r`, light distance
`r_L² = r² + H²`:

```
dΩ = dA · cosθ_L / r_L²,   cosθ_L = H / r_L
dΩ/dr = 2π H r / (r² + H²)^{3/2}
```

Equal solid angle per texel needs `n(r) = dN/dr ∝ dΩ/dr`, so the *wall area one
texel covers* is

```
A_t(r) = dA / (n dr) ∝ (r² + H²)^{3/2} = r_L³
```

and the local log-log exponent is

```
α(r) = d ln A_t / d ln r = 3r² / (r² + H²)     ∈ [0, 3)
```

`α → 0` at the light, `1.5` at `r = H`, `→ 3` far out. **There is no single
constant exponent**; fitting `r_L³` over the map's whole band gives
`α = 0.36` (LS), which is a band average, not the law. The allocation itself is
exact and does not need a constant — the constant only enters when the warp is
forced into a one-parameter family.

Note also what a uniform-area (planar / atlas) map is: `α = 0` on `A_t`, i.e.
*every texel covers the same wall area*, so each texel subtends the largest
solid angle closest to the light — exactly the compression the symptom
describes. `α > 0` is the requested near-dense re-allocation.

## Result

Leak = wall area where a lookup says *lit* although the analytic shadow is
dark. Notch = the reverse. Both in H², mean over 12 trials.

| allocation | near leak (r<0.5H) | far leak (r≥0.5H) | total |
|---|---|---|---|
| uniform, α = 0 | 0.00443 | 0.02397 | 0.02840 |
| constant α = 0.36 (band fit) | 0.00225 | 0.02305 | 0.02530 |
| **constant α = 1.25** | **0.00089** | **0.02216** | **0.02305** |
| constant α = 2.0 | 0.00091 | 0.02874 | 0.02965 |
| derived, exact warp | 0.00125 | 0.02201 | 0.02326 |
| reversed, α = −1 | 0.02366 | 0.03041 | 0.05407 |

Sweep over α (fig 3, fig 6):

* near field (r < 0.5H): leak falls steeply from α = −1 to a **minimum at
  α ≈ 1.75**, then rises; α = 3 is 5.5× worse than α = 2.
* far field (r ≥ 0.5H): **minimum at α ≈ 1.25**; above that the far wall runs
  out of texels and the curve turns up steeply.
* whole wall: minimum at α ≈ 1.25.

So: **the derived near-dense direction is right** — every α from 0.2 to 2 beats
the uniform map on the near field, and α ≈ 1.25–2.0 cuts the near-light leak by
4.8–5.8×. **The exponent is under-stated, not over-stated.** The band fit (0.36)
is far too small for the near field; the near-field optimum is ≈ 1.75. But a
single constant that large is bad for the far wall, which is why the exact,
non-constant warp (equal solid angle everywhere) is the right target if the
engine can express it — it is near-best in both bands without being the
optimum of either.

The "is the density relation backwards?" branch does **not** fire: the
sign-reversed map is worse everywhere (5.3× the uniform near leak). Fig 2
overlays the measured density on the derived curve to show this is not a
measurement artifact; the measured points track `r/(r²+H²)^{3/2}` and the
uniform line is flat.

### Why the far-field minimum is at 1.25, not 3

The wall is sampled only to `r = 2.4H`, so the ring beyond `r = 1.3H` is 70% of
the wall's area but carries only the largest, sparsest blockers. Pushing texels
outward past α ≈ 1.25 starves the middle band (`r ≈ 1H`, where `α(r) ≈ 1.5`)
faster than it helps the rim. The optimum therefore sits between the near and
far demands, as a single exponent must; this is also why the near and far optima
differ and why the exact warp, which follows `α(r)` rather than a constant, is
competitive on both.

## Files

| path | contents |
|---|---|
| `src/shadowmap.c` | simulator + self-test (`--selftest`) |
| `figs/plot.py` | figures 1–6 from `out/*.csv` |
| `out/points.csv` | per config × intensity × sample radius: coverage, wall area/texel, leak, notch, visible leak (raw and per unit area) |
| `out/sweep.csv` | α sweep at three intensities, split near/far |
| `out/summary.csv` | the table above, plus per-trial s.d. and area-normalised visible leak |
| `out/density.csv` | analytic derived density and local exponent vs r |
| `out/density_obs.csv` | measured density of the derived and uniform maps |
| `out/wall_*.ppm` | wall captures, three intensities × allocations × two layouts |

Figures: `fig1_exponent` (α(r) vs constants), `fig2_density` (predicted vs
measured density), `fig3_sweep` (near / far / whole-wall leak vs α),
`fig4_intensity` (visible leak density vs radius at I = 0.35 / 1 / 3),
`fig5_summary` (the table as bars), `fig6_pareto` (near vs far trade-off).

## Caveats

* Visibility is binary, no depth bias or PCF: this isolates resolution error.
  A real renderer adds bias/acne, which changes the magnitude but not the
  direction of the allocation result.
* The shadow map here is polar, so the exponent maps one-to-one onto row
  spacing. A perspective map needs the same exponent applied along both image
  axes, which is only equivalent for the near-field ring.
* The visible-patch threshold (`I·cosθ/r² > 0.10`, after a mild tonemap) only
  decides which leaked texels are counted as "visible"; it does not change
  which allocation is best. Measured geometry is intensity-independent, as the
  per-intensity tables show.
