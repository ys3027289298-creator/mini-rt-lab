# mut — image-oracle mutation testing for `bvh/` and `material/`

Answers one question: *if someone broke the intersector or the BRDF code, would
the picture change visibly?*  Each mutant is a real source patch, compiled
unmodified into a small path tracer; frames are compared against pristine ones.

No third-party libraries.  Python stdlib + `gcc`; PNG writing uses `zlib`.
The module's own self-test is still compiled and run for every mutant, as a
baseline for "what `make test` would have caught".

## Layout

```
mut/
  mutants.json          the catalogue: find/replace patches, kind, expectation
  render/
    rt.c                path tracer (scenes, camera, RNG, output format)
    iface.h             the only surface the renderer sees
    bvh_shim.c          #includes bvh/src/bvh.c, exports rt_bvh_*
    mat_shim.c          #includes material/src/material.c, exports rt_*
  tools/
    mutate.py           --check / --list / apply one patch
    validate.py         determinism, build stability, tree==brute, noise floor
    run_mutations.py    the suite -> results/results.json + per-mutant evidence
    sensitivity.py      smallest visible photometric error -> results/sensitivity.json
    raw2png.py          RT1 raw frame -> PNG (or amplified diff PNG)
    report.py           results/*.json -> REPORT.md
    compare_rev.py      "did my edit change the picture?" -- rev vs worktree
  results/              generated: results.json, validate.json, sensitivity.json,
                        frames/, mutants/<id>/{patched source, frames, diffs}
  REPORT.md             generated: tables, miss list, hashes, reproduce commands
```

## Run

```sh
cd mut
python3 -I tools/mutate.py --check
python3 -I tools/validate.py
python3 -I tools/run_mutations.py --out results     # full catalogue
python3 -I tools/sensitivity.py
python3 -I tools/report.py
```

One mutant, higher resolution:

```sh
python3 -I tools/run_mutations.py --only mat-pdf-no-pi --out /tmp/one \
        --width 640 --height 480 --spp 64
python3 -I tools/raw2png.py /tmp/one/mutants/mat-pdf-no-pi/grid.raw /tmp/grid.png
```

## Did my last edit change the picture?

```sh
python3 -I tools/compare_rev.py                 # HEAD vs working tree
python3 -I tools/compare_rev.py --rev de9c470  # any revision
```

Prints per-scene visible-pixel and max-channel deltas, writes
`results/compare/compare.json` plus frames and amplified diff images.  Exit 0
means byte-identical; exit 1 means the edit is visible.  Measured sensitivity
is ~2% relative photometric error, so "same" here means "below what 8-bit
output can show", not "mathematically identical".

## How a mutant is judged

- **Caught** iff some scene has ≥8 pixels differing by ≥3 levels in some RGB
  channel.  Applies to every mutant identically; scenes and threshold are not
  tuned per mutant.
- Each catalogue entry declares `expected: caught|survive`.  A mismatch is
  reported as a violation rather than silently tolerated.
  - `equivalent`: same rendered function, different code (rationale in the
    entry's `equivalence_note`).  Surviving is correct.
  - `defect-subthreshold`: a real defect measured to be below the 8-bit
    sensitivity floor.  Excluded from the headline rate and quantified by
    `sensitivity.py`.

The comparison is paired (same camera/resolution/spp/seed/RNG).  Monte-Carlo
noise cancels exactly, which `validate.py` confirms; unpaired renders differ
in ~30–50% of pixels, so any unpaired comparison would be meaningless.

## Two things that make the oracle trustworthy

1. **The modules are not copied.** `bvh_shim.c` / `mat_shim.c` `#include` the
   real source file and only rename its `main`.  A mutation cannot be missed
   because a stale copy was compiled.
2. **Tree shape must not matter.** The BVH build permutes the primitive array
   and boxes can be hit at exactly the same `t`; the shim resolves ties by
   lower primitive id.  `validate.py` proves the BVH path matches an
   unaccelerated all-boxes path on every scene, so a shape-only change is
   invisible by construction, not by luck.

## Adding a mutant

Append an entry to `mutants.json` with a `find` string that occurs exactly once
in the pristine module.  `python3 -I tools/mutate.py --check` verifies that
before anything is built.  Mutants must not need new files or link flags; the
harness compiles `rt.c`, the two shims, and `-lm` only.

## Requirements / limits

- Frames are only comparable within one build and one seed.  The oracle cannot
  compare two renderers that cannot share a sample stream.
- Sensitivity floor (measured): ~2% relative photometric error for any scene,
  ~5% to be visible in all five.  Smaller defects are not observable at this
  exposure.
- Five scenes, one camera each.  Geometry/camera configurations outside them
  are not covered.
