CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: material albedo_ec

material: src/material.c
	$(CC) $(CFLAGS) $< -o $@ -lm

albedo_ec: src/albedo_ec.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: material albedo_ec
	./material
	./albedo_ec

clean:
	rm -f material albedo_ec
