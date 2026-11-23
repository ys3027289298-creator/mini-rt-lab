# kernel -- one random source behind three rendering sites

The three rendering modules (`bvh/`, `sample/`, `material/`) and their public
self-tests are **untouched**. This directory adds a small kernel that drives
three stochastic sites -- attenuation, material sampling and stratified
sampling -- and proves that, with the unified random source enabled, each
site's output depends only on the work item, never on thread, visit order or
the machine that built it.

## Why there are two configurations

The checkout at `de9c470` has no random source, no threads and no attenuation
code at all; the three modules are deterministic and expose nothing but
`main`. So "before" is defined here, once, and frozen:

- `src/rt_legacy.c` -- the pre-unification shape: one private generator per
  site (`LCG` for attenuation, a second `LCG` for material, `xorshift64` for
  sampling), each with its own mutable state and draw order. This is the
  common pattern that makes a scene change when work moves between threads or
  machines.
- `src/rt_unified.c` -- all three sites draw from `src/rt_rng.c`, a
  counter-based source with **no generator state**: every draw is a pure
  function of `(site, item, stream)`.

Both implement the same frozen header `include/rt_kernel.h`.

Provenance, stated plainly: this checkout contained no random source to take a
"before" from, so `src/rt_legacy.c` **is** the before -- written first, its
output frozen into `baseline/legacy/`, and only then was `src/rt_unified.c`
written. The OFF-switch check therefore proves "the new code can reproduce the
behaviour that was frozen", not "an unknown pre-existing build was recovered".
`baseline/legacy/` is newer than `src/rt_legacy.c` and older than
`src/rt_unified.c`; `make test` re-derives it from the current OFF build every
run and fails if it drifts.

## Run it

```sh
cd kernel
make test          # build + full acceptance suite, prints PASS/FAIL per check
make portability   # same results across gcc/clang and -O0..-O3
./scripts/report.sh        # regenerate evidence/REPORT.md
```

From the repository root the same thing is `make test` (module self-tests
first, then the kernel suite).

Nothing is installed and no tracked file is modified; build products go to
`kernel/build/`, evidence to `kernel/evidence/`.

## Layout

```
kernel/
├── include/rt_kernel.h     frozen public API (hash-checked)
├── include/rt_rng.h        the unified counter-based source
├── include/rt_shim.h       read-only access to the modules' own functions
├── src/rt_rng.c            splitmix64 mixer + domain-separated keys
├── src/rt_unified.c        sites -> unified source (default)
├── src/rt_legacy.c         sites -> per-site generators (OFF switch)
├── src/rt_shim_*.c         #include the module source, rename main only
├── tools/site_probe.c      per-site output, N threads, seq/strided order
├── tools/site_stats.c      bias / standard error / RMSE vs analytic target
├── tools/scene_hash.c      all three sites -> one 64-bit scene fingerprint
├── tools/rng_selftest.c    purity, domain separation, uniformity
├── tools/sigcheck.c        compile-time signature assertions
├── baseline/               frozen before-unification goldens
├── scripts/check.sh        the acceptance suite
├── scripts/portability.sh  compiler / optimisation matrix
└── scripts/report.sh       evidence/REPORT.md
```

## What "verified" means here

Each of these is a `cmp`/hash comparison of bytes or doubles, not a visual
judgement:

| Requirement | Where | What is compared |
|---|---|---|
| Public signatures unchanged | `check.sh` §2, `sigcheck.c` | header sha256; every function assigned to a fixed-signature pointer; exported symbol tables of both builds identical |
| OFF switch reproduces current output byte-for-byte | `check.sh` §4 | `cmp` of 256-item stdout against `baseline/legacy/` |
| Thread independence | `check.sh` §5 | identical per-site output at 1,2,3,4,5,7,8 threads and on repeat |
| Schedule independence | `check.sh` §5b, §5c | same items in a different, **thread-free** visit order; whole-scene FNV hash unchanged |
| Error before/after | `check.sh` §6 | 20 000-item mean vs analytic target for each site, within 5 standard errors; 80 000-item convergence line |
| No hidden state | `check.sh` §7 | `nm` shows no writable data in the unified objects |
| Another machine | `make portability` | identical scene hash from 8 gcc/clang × -O0..-O3 builds; OFF build still matches the baseline under both compilers |
| Memory / data races | `check.sh` §8 | ASan+UBSan on both builds, TSan on the unified 8-thread probe |

The legacy path is deliberately kept and shown to **fail** the schedule check
(§5b) and to differ from the unified scene hash (§5c) -- that contrast is what
demonstrates the check is sensitive to the defect rather than vacuous.

## Sites

| Site | Estimator | Analytic target |
|---|---|---|
| attenuation | `t = -ln(1-u)/sigma`, weight `exp(-sigma t)` | `E[weight] = 1/2` for every `sigma > 0` |
| material | mean `cos` of `n` cosine-weighted hemisphere samples, using `material/`'s own `cosine_sample()` | `2/3` |
| sampling | jittered stratified estimate of `∫u² du`, using `sample/`'s own `f_quad()` | `1/3` |

The material and sampling sites call the existing modules through the shims,
so the modules' code is what is being tested; the kernels only supply the
random draw.

## A caveat about "another machine"

The checkout here is one WSL2 x86_64 box. The portability check is the
strongest available proxy: it varies compiler and optimisation level, which is
where unordered evaluation and floating-point contraction would show up. (That
check caught a real one: two `material_u()` calls in one argument list bound
`u1`/`u2` differently under gcc and clang. Fixed by sequencing the draws.) It
does not prove reproducibility on a different CPU architecture; the design
argument for that is that the unified source is integer-only with no floating
point in the generator and no dependence on draw order.
