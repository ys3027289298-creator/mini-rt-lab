PYTHON ?= python3
CACHE := $(CURDIR)/.pycache

# Isolated interpreter mode (-B) keeps the module from writing .pyc next to the
# sources; PYTHONDONTWRITEBYTECODE is the belt-and-braces version.
export PYTHONDONTWRITEBYTECODE := 1
RUN := $(PYTHON) -B

FRAME_ARGS := --width 320 --height 240 --spp 64 --block 16 --bits 16 --rank-max 16

all: test

test:
	$(RUN) src/selftest.py

run:
	$(RUN) src/run_compression.py $(FRAME_ARGS) --outdir out

# smaller run, for a quick check on a slow machine
run-small:
	$(RUN) src/run_compression.py --width 160 --height 112 --spp 16 \
		--block 16 --bits 16 --rank-max 16 --outdir out

clean:
	rm -rf out

.PHONY: all test run run-small clean
