CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: mis

mis: src/mis.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: mis
	./mis

clean:
	rm -f mis
