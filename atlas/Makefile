CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: atlas

atlas: src/atlas.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: atlas
	./atlas

clean:
	rm -f atlas *.ppm
