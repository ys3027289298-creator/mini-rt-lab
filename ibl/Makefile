CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: ibl

ibl: src/ibl.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: ibl
	./ibl

clean:
	rm -f ibl
