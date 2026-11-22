CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: fov

fov: src/fov.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: fov
	./fov

clean:
	rm -f fov
