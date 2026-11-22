CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: decorr

decorr: src/decorr.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: decorr
	./decorr

clean:
	rm -f decorr
