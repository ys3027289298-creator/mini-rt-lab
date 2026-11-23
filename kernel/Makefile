# kernel/Makefile
#
#   make            # build every probe and checker (both configurations)
#   make test       # build + run the full acceptance suite (scripts/check.sh)
#   make portability# cross-compiler / cross-optimisation reproducibility
#   make off        # the same build with the unified source disabled
#   make clean
#
# Both configurations share the same public header and the same exported
# symbols; only src/rt_unified.c versus src/rt_legacy.c differs.
#
# -ffp-contract=off and no -ffast-math: results must reproduce bit-for-bit
# across runs, builds and machines, not merely agree closely.
CC      ?= gcc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -ffp-contract=off -fno-fast-math
CPPFLAGS = -Iinclude
LDLIBS  = -lm -pthread

RNG_UNIFIED_SRC = src/rt_rng.c src/rt_unified.c
RNG_LEGACY_SRC  = src/rt_rng.c src/rt_legacy.c
SHIM_SRC        = src/rt_shim_material.c src/rt_shim_sample.c

ITEMS  ?= 256
SIGMA  ?= 1.7
N      ?= 8
STRATA ?= 16

SITE_FLAGS = -DPROBE_SIGMA=$(SIGMA) -DPROBE_N=$(N) -DPROBE_STRATA=$(STRATA) \
             -DPROBE_ITEMS=$(ITEMS)

all: all-probes

# The switch: which implementation object the kernel links against.
#   unified (default) -> src/rt_unified.c : one counter-based source, no state
#   off               -> src/rt_legacy.c  : pre-unification per-site generators
# Both expose the identical header and symbols; only the random source differs.
unified: all-probes
	@echo "unified random source -> build/unified_scene (default configuration)"

off: all-probes
	@echo "unified random source OFF -> build/legacy_scene (pre-unification behaviour)"

all-probes: build/sigcheck build/rng_selftest build/stats_legacy build/stats_unified \
            build/legacy_site0 build/legacy_site1 build/legacy_site2 \
            build/unified_site0 build/unified_site1 build/unified_site2 \
            build/legacy_api.o build/unified_api.o \
            build/symbols_legacy.txt build/symbols_unified.txt \
            build/scene_hash_legacy build/scene_hash_unified

build:
	@mkdir -p build

# ---- probes: one binary per site per configuration -----------------------

build/legacy_site%: tools/site_probe.c $(RNG_LEGACY_SRC) $(SHIM_SRC) \
                    include/rt_kernel.h include/rt_rng.h include/rt_shim.h | build
	$(CC) $(CFLAGS) -DPROBE_SITE=$* $(SITE_FLAGS) $(CPPFLAGS) \
	      tools/site_probe.c $(RNG_LEGACY_SRC) $(SHIM_SRC) -o $@ $(LDLIBS)

build/unified_site%: tools/site_probe.c $(RNG_UNIFIED_SRC) $(SHIM_SRC) \
                     include/rt_kernel.h include/rt_rng.h include/rt_shim.h | build
	$(CC) $(CFLAGS) -DPROBE_SITE=$* $(SITE_FLAGS) $(CPPFLAGS) \
	      tools/site_probe.c $(RNG_UNIFIED_SRC) $(SHIM_SRC) -o $@ $(LDLIBS)

build/scene_hash_unified: tools/scene_hash.c $(RNG_UNIFIED_SRC) $(SHIM_SRC) \
                          include/rt_kernel.h include/rt_rng.h include/rt_shim.h | build
	$(CC) $(CFLAGS) $(CPPFLAGS) tools/scene_hash.c $(RNG_UNIFIED_SRC) $(SHIM_SRC) -o $@ $(LDLIBS)
	cp $@ build/unified_scene

build/scene_hash_legacy: tools/scene_hash.c $(RNG_LEGACY_SRC) $(SHIM_SRC) \
                         include/rt_kernel.h include/rt_rng.h include/rt_shim.h | build
	$(CC) $(CFLAGS) $(CPPFLAGS) tools/scene_hash.c $(RNG_LEGACY_SRC) $(SHIM_SRC) -o $@ $(LDLIBS)
	cp $@ build/legacy_scene

# ---- accuracy statistics -------------------------------------------------

build/stats_legacy: tools/site_stats.c $(RNG_LEGACY_SRC) $(SHIM_SRC) \
                    include/rt_kernel.h include/rt_rng.h include/rt_shim.h | build
	$(CC) $(CFLAGS) $(CPPFLAGS) tools/site_stats.c $(RNG_LEGACY_SRC) $(SHIM_SRC) -o $@ $(LDLIBS)

build/stats_unified: tools/site_stats.c $(RNG_UNIFIED_SRC) $(SHIM_SRC) \
                     include/rt_kernel.h include/rt_rng.h include/rt_shim.h | build
	$(CC) $(CFLAGS) $(CPPFLAGS) tools/site_stats.c $(RNG_UNIFIED_SRC) $(SHIM_SRC) -o $@ $(LDLIBS)

# ---- API checks ----------------------------------------------------------

# -Werror=incompatible-pointer-types: a changed signature must be a hard build
# error, not a warning that scrolls past.
build/sigcheck: tools/sigcheck.c include/rt_kernel.h $(RNG_UNIFIED_SRC) $(SHIM_SRC) \
                include/rt_rng.h include/rt_shim.h | build
	$(CC) $(CFLAGS) -Werror=incompatible-pointer-types $(CPPFLAGS) \
	      tools/sigcheck.c $(RNG_UNIFIED_SRC) $(SHIM_SRC) -o $@ $(LDLIBS)

build/rng_selftest: tools/rng_selftest.c src/rt_rng.c include/rt_rng.h | build
	$(CC) $(CFLAGS) $(CPPFLAGS) tools/rng_selftest.c src/rt_rng.c -o $@ $(LDLIBS)

# Partial links, so the exported symbol tables can be compared directly.
build/legacy_api.o: $(RNG_LEGACY_SRC) $(SHIM_SRC) | build
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_rng.c -o build/rt_rng_l.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_legacy.c -o build/rt_impl_l.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_shim_material.c -o build/rt_shim_m_l.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_shim_sample.c -o build/rt_shim_s_l.o
	ld -r build/rt_rng_l.o build/rt_impl_l.o build/rt_shim_m_l.o build/rt_shim_s_l.o -o $@

build/unified_api.o: $(RNG_UNIFIED_SRC) $(SHIM_SRC) | build
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_rng.c -o build/rt_rng_u.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_unified.c -o build/rt_impl_u.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_shim_material.c -o build/rt_shim_m_u.o
	$(CC) $(CFLAGS) $(CPPFLAGS) -c src/rt_shim_sample.c -o build/rt_shim_s_u.o
	ld -r build/rt_rng_u.o build/rt_impl_u.o build/rt_shim_m_u.o build/rt_shim_s_u.o -o $@

build/symbols_legacy.txt: build/legacy_api.o
	nm -g --defined-only $< | awk '{print $$2, $$3}' | sort > $@

build/symbols_unified.txt: build/unified_api.o
	nm -g --defined-only $< | awk '{print $$2, $$3}' | sort > $@

# ---- entry points --------------------------------------------------------

test: all-probes
	./scripts/check.sh

portability: all-probes
	./scripts/portability.sh

negative: all-probes
	./scripts/negative_controls.sh

clean:
	rm -rf build

.PHONY: all unified off all-probes test portability negative clean
