CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Werror -ffp-contract=off

all: sample

sample: src/sample.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sample
	./sample

render: sample
	./sample render sample.ppm

clean:
	rm -f sample sample.ppm
