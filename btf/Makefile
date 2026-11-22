CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: btf

btf: src/btf.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: btf
	./btf

clean:
	rm -f btf
