CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: shutter

shutter: src/shutter.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: shutter
	./shutter

clean:
	rm -f shutter
