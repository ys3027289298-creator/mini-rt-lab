CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: allocate

allocate: src/allocate.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: allocate
	./allocate

clean:
	rm -f allocate
