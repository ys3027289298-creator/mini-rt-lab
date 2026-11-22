CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: viscache

viscache: src/viscache.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: viscache
	./viscache

clean:
	rm -f viscache
