# Top-level orchestrator: unified RNG builds, three acceptance artifacts,
# legacy byte-identity check and the cross-artifact Python verification.
CC ?= gcc
MODULES := bvh sample material
ARTDIR := artifacts/acceptance
GOLDEN := artifacts/golden

.PHONY: all verify legacy-check cross-check clean

all: verify

# Build the three verifiers (unified source).
$(ARTDIR):
	mkdir -p $(ARTDIR)

verify: | $(ARTDIR)
	$(MAKE) -C sample verify_hist CC=$(CC)
	$(MAKE) -C material verify_variance CC=$(CC)
	$(MAKE) -C bvh verify_image CC=$(CC)
	sample/verify_hist $(ARTDIR)
	material/verify_variance $(ARTDIR)
	bvh/verify_image $(ARTDIR)
	python3 -I verify_artifacts.py $(ARTDIR)

# Off-switch proof: legacy builds must reproduce the captured golden stdout
# byte-for-byte in every module.
legacy-check:
	@for m in $(MODULES); do \
	  $(MAKE) -C $$m clean >/dev/null; \
	  $(MAKE) -C $$m $$m-legacy CC=$(CC) >/dev/null; \
	  $$m/$$m-legacy > /tmp/$$m-legacy.out; \
	  cmp /tmp/$$m-legacy.out $(GOLDEN)/$$m.out \
	    && echo "$$m legacy: byte-identical to golden"; \
	done

# Determinism proof: rebuild with clang and diff every artifact byte-for-byte.
cross-check: | $(ARTDIR)
	mkdir -p artifacts/clang
	clang -O2 -std=c11 -DRT_UNIFIED_RNG -Ishared -Dmain=sample_selftest_main \
	  -Wno-unused-function sample/verify_hist.c shared/rng.c \
	  -o artifacts/tmp/vh_clang -lm
	clang -O2 -std=c11 -DRT_UNIFIED_RNG -Ishared -Dmain=material_selftest_main \
	  -Wno-unused-function material/verify_variance.c shared/rng.c \
	  -o artifacts/tmp/vv_clang -lm
	clang -O2 -std=c11 -DRT_UNIFIED_RNG -Ishared -Dmain=bvh_selftest_main \
	  -Wno-unused-function bvh/verify_image.c shared/rng.c \
	  -o artifacts/tmp/vi_clang -lm
	artifacts/tmp/vh_clang artifacts/clang >/dev/null
	artifacts/tmp/vv_clang artifacts/clang >/dev/null
	artifacts/tmp/vi_clang artifacts/clang >/dev/null
	@for f in histogram.csv variance.json image.ppm image_hash.json; do \
	  cmp $(ARTDIR)/$$f artifacts/clang/$$f \
	    && echo "$$f: gcc == clang byte-identical"; \
	done

clean:
	@for m in $(MODULES); do $(MAKE) -C $$m clean >/dev/null; done
	rm -rf artifacts/acceptance artifacts/runB artifacts/clang artifacts/tmp
