CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: reflect

reflect: src/reflect.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: reflect
	./reflect

clean:
	rm -f reflect
