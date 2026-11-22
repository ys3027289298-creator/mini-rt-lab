CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: cylinder

cylinder: src/cylinder.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: cylinder
	./cylinder

clean:
	rm -f cylinder
