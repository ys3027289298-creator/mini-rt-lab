CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: sched

sched: src/sched.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: sched
	./sched

clean:
	rm -f sched
