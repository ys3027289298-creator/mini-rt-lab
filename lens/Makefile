CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: lens

lens: src/lens.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: lens
	./lens

clean:
	rm -f lens
