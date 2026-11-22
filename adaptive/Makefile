CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
CSV ?= artifacts/adaptive_pixels.csv

all: adaptive

adaptive: src/adaptive.c
	$(CC) $(CFLAGS) $< -o $@ -lm

# Per-pixel details go to artifacts/adaptive_pixels.csv; override with CSV=...
test: adaptive
	mkdir -p $(dir $(CSV))
	./adaptive $(CSV)

clean:
	rm -f adaptive
