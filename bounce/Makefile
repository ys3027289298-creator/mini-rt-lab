CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: bounce

bounce: src/bounce.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: bounce
	./bounce

clean:
	rm -f bounce
