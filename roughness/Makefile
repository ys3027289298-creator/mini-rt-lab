CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: roughness

roughness: src/roughness.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: roughness
	./roughness

plot: out/mapping_curve.csv
	python3 tools/plot_curves.py

out/mapping_curve.csv: roughness
	./roughness

clean:
	rm -f roughness
	rm -rf out
