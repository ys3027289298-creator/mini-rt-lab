CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: sample

sample: src/sample.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sample
	./sample

clean:
	rm -f sample
