CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: ao

ao: src/ao.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: ao
	./ao

clean:
	rm -f ao
