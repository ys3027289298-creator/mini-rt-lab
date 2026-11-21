CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: texture

texture: src/texture.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: texture
	./texture

clean:
	rm -f texture
