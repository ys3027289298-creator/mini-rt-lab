CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: light

light: src/light.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: light
	./light

clean:
	rm -f light
