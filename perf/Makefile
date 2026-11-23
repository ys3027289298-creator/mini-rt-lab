# perf: performance-budget acceptance harness over the three rendering paths.
#
# Each path is compiled with the module source included unmodified; the
# module's own `main` is renamed so the harness supplies the entry point.
CC      ?= gcc
# -Wno-unused-function: the module sources are included unmodified and each
# carries helpers its self-test does not call (sample.c's halton, for one).
# Our own code in perf_main.c still builds under -Werror.
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra -Werror -Wno-unused-function -ffp-contract=off
CPPFLAGS = -D_POSIX_C_SOURCE=200809L -Iinclude
LDLIBS   = -lm
BIN      = perf
ARTDIR   = ../artifacts/perf

SRCS = src/perf_main.c src/perf_bvh.c src/perf_sample.c src/perf_material.c

$(BIN): $(SRCS) include/perf_render.h include/perf_rng.h budget.conf
	$(CC) $(CFLAGS) $(CPPFLAGS) $(SRCS) -o $@ $(LDLIBS)

all: $(BIN)

test: $(BIN)
	./$(BIN) selftest

low: $(BIN)
	mkdir -p $(ARTDIR)/low
	./$(BIN) run low --budget budget.conf --outdir $(ARTDIR)/low

high: $(BIN)
	mkdir -p $(ARTDIR)/high
	./$(BIN) run high --budget budget.conf --outdir $(ARTDIR)/high

# Off-switch proof: with budget checks disabled the rendered pixels must be
# unchanged, including the PPMs on disk.  Two no-check runs of each tier also
# prove the measurement itself is deterministic.
no-check: $(BIN)
	mkdir -p $(ARTDIR)/nocheck-low $(ARTDIR)/nocheck-high
	./$(BIN) run low --budget budget.conf --outdir $(ARTDIR)/nocheck-low --no-budget-check
	./$(BIN) run high --budget budget.conf --outdir $(ARTDIR)/nocheck-high --no-budget-check
	python3 -I tools/verify_tiers.py --dir $(ARTDIR) \
	    --nocheck-low $(ARTDIR)/nocheck-low \
	    --nocheck-high $(ARTDIR)/nocheck-high

# Full gate: self-test, both tiers, verification, then the off-switch proof.
accept: $(BIN)
	./$(BIN) selftest
	mkdir -p $(ARTDIR)
	{ echo "host: $$(uname -srm)"; \
	  echo "cc: $$($(CC) --version | head -1)"; \
	  echo "nproc: $$(nproc)"; \
	  echo "date: $$(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > $(ARTDIR)/env.txt
	mkdir -p $(ARTDIR)/low $(ARTDIR)/high
	./$(BIN) run low --budget budget.conf --outdir $(ARTDIR)/low
	./$(BIN) run high --budget budget.conf --outdir $(ARTDIR)/high
	python3 -I tools/verify_tiers.py --dir $(ARTDIR)
	mkdir -p $(ARTDIR)/nocheck-low $(ARTDIR)/nocheck-high
	./$(BIN) run low --budget budget.conf --outdir $(ARTDIR)/nocheck-low --no-budget-check
	./$(BIN) run high --budget budget.conf --outdir $(ARTDIR)/nocheck-high --no-budget-check
	python3 -I tools/verify_tiers.py --dir $(ARTDIR) \
	    --nocheck-low $(ARTDIR)/nocheck-low \
	    --nocheck-high $(ARTDIR)/nocheck-high

clean:
	rm -f $(BIN)

.PHONY: all test low high accept no-check clean
