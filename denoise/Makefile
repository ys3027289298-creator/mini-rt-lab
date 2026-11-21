CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: denoise

denoise: src/denoise.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: denoise
	./denoise

clean:
	rm -f denoise
