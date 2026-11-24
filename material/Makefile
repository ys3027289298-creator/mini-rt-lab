CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: material

material: src/material.c src/ggx.h
	$(CC) $(CFLAGS) src/material.c -o $@ -lm

highlight: src/highlight.c src/ggx.h
	$(CC) $(CFLAGS) src/highlight.c -o $@ -lm

test: material
	./material

exp: highlight
	mkdir -p ../artifacts
	./highlight ../artifacts
	python3 tools/compose.py ../artifacts

clean:
	rm -f material highlight
