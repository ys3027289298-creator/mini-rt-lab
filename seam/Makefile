CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: seam

seam: src/seam.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: seam
	./seam

clean:
	rm -f seam
