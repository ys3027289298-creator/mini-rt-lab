CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: dof

dof: src/dof.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: dof
	./dof

clean:
	rm -f dof
