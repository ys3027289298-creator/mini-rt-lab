CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: temporal

temporal: src/temporal.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: temporal
	./temporal

clean:
	rm -f temporal
