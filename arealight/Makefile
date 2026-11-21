CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: arealight

arealight: src/arealight.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: arealight
	./arealight

clean:
	rm -f arealight
