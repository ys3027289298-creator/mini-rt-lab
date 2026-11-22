CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: hybrid

hybrid: src/hybrid.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: hybrid
	./hybrid

images: hybrid
	mkdir -p out
	./hybrid --dump out

sweep: hybrid
	./hybrid --sweep

clean:
	rm -f hybrid
	rm -rf out
