CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: select

select: src/select.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: select
	./select

clean:
	rm -f select
