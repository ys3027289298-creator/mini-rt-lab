CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Werror -ffp-contract=off

all: material

material: src/material.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: material
	./material

render: material
	./material render material.ppm

clean:
	rm -f material material.ppm
