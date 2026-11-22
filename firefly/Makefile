CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: firefly

firefly: src/firefly.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: firefly
	./firefly

clean:
	rm -f firefly
