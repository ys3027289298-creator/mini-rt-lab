CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
LDLIBS = -lm -lpthread

all: render compare

render: src/render.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

compare: src/compare.c
	$(CC) $(CFLAGS) $< -o $@ -lm

# Temporal consistency: same parameters must render bit-identically,
# across runs and across thread counts.
test: render compare
	./scripts/consistency_test.sh

clean:
	rm -f render compare

.PHONY: all test clean
