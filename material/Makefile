CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: material

material: src/material.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: material
	./material

clean:
	rm -f material
