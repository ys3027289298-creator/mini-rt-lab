CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: differential

differential: src/differential.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: differential
	./differential

clean:
	rm -f differential
