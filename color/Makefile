CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: color

color: src/color.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: color
	./color

clean:
	rm -f color
