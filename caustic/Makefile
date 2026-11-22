CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: caustic

caustic: src/caustic.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: caustic
	./caustic

clean:
	rm -f caustic
