CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: shadow

shadow: src/shadow.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: shadow
	./shadow

clean:
	rm -f shadow
