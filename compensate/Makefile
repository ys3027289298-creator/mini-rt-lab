CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: compensate

compensate: src/compensate.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: compensate
	./compensate

clean:
	rm -f compensate
