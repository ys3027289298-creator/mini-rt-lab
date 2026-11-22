CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: upsample

upsample: src/upsample.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: upsample
	./upsample

clean:
	rm -f upsample
