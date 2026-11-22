CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: spectrum

spectrum: src/spectrum.c src/cie1931.h
	$(CC) $(CFLAGS) src/spectrum.c -o $@ -lm

test: spectrum
	./spectrum

# spectral-grating strip rendered from the module (P6 PPM, 256x64)
image: spectrum
	./spectrum grating.ppm

clean:
	rm -f spectrum grating.ppm
