CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: oit

oit: src/oit.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: oit
	./oit

clean:
	rm -f oit
