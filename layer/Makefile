CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: layer

layer: src/layer.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: layer
	./layer

clean:
	rm -f layer
