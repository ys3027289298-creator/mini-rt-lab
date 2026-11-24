CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: spot

spot: src/spot.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: spot
	./spot

clean:
	rm -f spot
