# Random-source unification -- acceptance evidence

Generated 2026-11-23T03:30:44Z from `kernel/evidence/`.
Regenerate with `make -C kernel test && kernel/scripts/report.sh`.

## 1. Public API frozen

Header hash (baseline vs now):

```
baseline: f667dfd3286d95537fb1068bab72f4ecceaf43fb133fa72e99570836667e9bcb  include/rt_kernel.h
now:      f667dfd3286d95537fb1068bab72f4ecceaf43fb133fa72e99570836667e9bcb  include/rt_kernel.h
```

Compile-time signature assertions:

```
sig_atten=1
sig_mean_cos=1
sig_strat=1
sig_configure=1
record_fields_t=1
record_ptr_t=1
all_ok=1
```

Exported symbols (identical between the OFF build and the unified build):

```
T material_module_selftest_main
T rt_atten_sample
T rt_material_mean_cos
T rt_rng_configure
T rt_rng_key
T rt_rng_mix64
T rt_rng_u64
T rt_rng_unit
T rt_shim_cosine_sample
T rt_shim_f_quad
T rt_strat_quad
T sample_module_selftest_main
```

## 2. Existing modules untouched

Self-test stdout byte-compared against the frozen baseline:

```
bvh       baseline ddaa0f9f631fb76976932495820c40d20f7e96ba97353785385e4996c81283a7
bvh       now      ddaa0f9f631fb76976932495820c40d20f7e96ba97353785385e4996c81283a7
sample    baseline e8b7368660e532f28a1334cd521d49a6b6b72b79038eff66cab848855266b79c
sample    now      e8b7368660e532f28a1334cd521d49a6b6b72b79038eff66cab848855266b79c
material  baseline 40e65e5cfe92eeb741bb179dd2f164829cfe82e7bb2f9e63d75f5c1deb64dd5d
material  now      40e65e5cfe92eeb741bb179dd2f164829cfe82e7bb2f9e63d75f5c1deb64dd5d
```

Tracked-file diff over bvh/, sample/, material/ (empty means untouched):

```
```

## 3. Random source property tests

```
purity_same_key_same_draw=1
domain_separation_distinct_keys=1
uniform_in_half_open_unit=1
uniform_mean=1
bit_balance=1
no_short_cycle=1
draw_order_independent=1
all_ok=1
```

## 4. OFF switch: byte-identical to the frozen pre-change output

The OFF build uses the pre-unification per-site generators. For each site,
its 256-item stdout was compared with `cmp` against the baseline captured
before unification:

```
site0 baseline 99c293a8de7078ce7d3155a0e75afce4093bd1f8abf81e3faaa566e156c1cf70
site0 current  99c293a8de7078ce7d3155a0e75afce4093bd1f8abf81e3faaa566e156c1cf70
site1 baseline e5e0b7d3a5886aea90bd47328dbb74bbf9fcf14f984f936aef8452403159cfde
site1 current  e5e0b7d3a5886aea90bd47328dbb74bbf9fcf14f984f936aef8452403159cfde
site2 baseline 58b02525e8f39e0a2cd659bc88b8fd1bc73b6f60dbc7333ca85dfd4079455c4d
site2 current  58b02525e8f39e0a2cd659bc88b8fd1bc73b6f60dbc7333ca85dfd4079455c4d
```

## 5. Thread and schedule independence

Scene fingerprint, all three sites folded into one FNV-1a hash, across
thread counts and visit orders:

```
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=2 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=3 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=4 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=7 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=8 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=2 order=strided
scene_hash 95ce887092e4a8a2 items=256 threads=3 order=strided
scene_hash 95ce887092e4a8a2 items=256 threads=8 order=strided
-- legacy path, for contrast
scene_hash c4541511fc650094 items=256 threads=1 order=seq
scene_hash f150cd99fd482f08 items=256 threads=3 order=strided
```

Per-site hashes at each thread count are compared in `check.sh` section 5.

## 6. Error before vs after

Each site estimates an analytic target; 20000 items per run. The verdict
requires the mean to be within 5 standard errors of the target.

```
### attenuation (site 0, param 1.7, 20000 items)
-- before: legacy per-site RNG
site=0 param=1.7 items=20000
mean=0.49785201111733912 target=0.5 bias=-0.00214799 stderr=0.0020406 rmse=0.00296275
verdict=PASS
-- after: unified counter-based source
site=0 param=1.7 items=20000
mean=0.50013417578329611 target=0.5 bias=0.000134176 stderr=0.00204481 rmse=0.0020492
verdict=PASS

### material (site 1, param 8, 20000 items)
-- before: legacy per-site RNG
site=1 param=8 items=20000
mean=0.66768852949690916 target=0.66666666666666663 bias=0.00102186 stderr=0.00059258 rmse=0.00118125
verdict=PASS
-- after: unified counter-based source
site=1 param=8 items=20000
mean=0.6665053926559551 target=0.66666666666666663 bias=-0.000161274 stderr=0.000585491 rmse=0.000607297
verdict=PASS

### sampling (site 2, param 16, 20000 items)
-- before: legacy per-site RNG
site=2 param=16 items=20000
mean=0.33328651655279684 target=0.33333333333333331 bias=-4.68168e-05 stderr=3.6918e-05 rmse=5.96217e-05
verdict=PASS
-- after: unified counter-based source
site=2 param=16 items=20000
mean=0.33331002719604597 target=0.33333333333333331 bias=-2.33061e-05 stderr=3.66538e-05 rmse=4.34359e-05
verdict=PASS

### convergence (4x items)
80k items: site=0 param=1.7 items=80000 mean=0.50169693435466356 target=0.5 bias=0.00169693 stderr=0.00102201 rmse=0.00198093 verdict=PASS 
80k items: site=1 param=8 items=80000 mean=0.66642449329052178 target=0.66666666666666663 bias=-0.000242173 stderr=0.000294721 rmse=0.000381456 verdict=PASS 
80k items: site=2 param=16 items=80000 mean=0.3333285499896973 target=0.33333333333333331 bias=-4.78334e-06 stderr=1.84002e-05 rmse=1.90118e-05 verdict=PASS 

```

## 7. Mutable generator state

```
--- writable data symbols per object ---
rt_impl_l.o    s_atten_state(d) s_material_state(d) s_sampling_state(d) 
rt_impl_u.o    
rt_rng_l.o     
rt_rng_u.o     
```

## 7b. Negative controls (the checks are not vacuous)

Each control injects one defect into a throwaway copy and requires the
corresponding check to fail.  A suite that passes mutants too would
measure nothing; these must all report the defect as detected.

```
PASS  control stateful_order_dependence: defect detected
PASS  control weak_mixer: defect detected
PASS  control bad_signature: defect detected
PASS  control biased_estimator: defect detected

RESULT: PASS (every injected defect was detected)
```

## 8. Cross-compiler / cross-optimisation

```
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
scene_hash 95ce887092e4a8a2 items=256 threads=1 order=seq
```
