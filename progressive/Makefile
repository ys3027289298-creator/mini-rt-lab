CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: progressive

progressive: src/progressive.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: progressive
	./progressive

clean:
	rm -f progressive
