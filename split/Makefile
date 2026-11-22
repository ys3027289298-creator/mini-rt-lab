CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: splitting

splitting: src/splitting.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: splitting
	./splitting

clean:
	rm -f splitting
