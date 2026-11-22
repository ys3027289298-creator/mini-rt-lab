CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
USER_SITE := $(shell python3 -c 'import site; print(site.getusersitepackages())')

BIN = raster tracer

all: $(BIN)

raster: src/raster.c src/common.h
	$(CC) $(CFLAGS) src/raster.c -o $@ -lm

tracer: src/tracer.c src/common.h
	$(CC) $(CFLAGS) src/tracer.c -o $@ -lm

test: all
	PYTHONPATH="$(USER_SITE)" python3 -E tools/run_check.py

clean:
	rm -f $(BIN)
	rm -rf out
