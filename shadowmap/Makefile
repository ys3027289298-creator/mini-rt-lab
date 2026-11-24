CC      ?= cc
CFLAGS  ?= -O2 -std=c99 -Wall -Wextra
LDLIBS  ?= -lm
BIN     := shadowmap

all: $(BIN) test

$(BIN): src/shadowmap.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

test: $(BIN)
	./$(BIN) --selftest

run: $(BIN)
	./$(BIN) out/points.csv 1

figures: run
	python3 figs/plot.py

clean:
	rm -f $(BIN)
	rm -f out/*.csv out/*.ppm out/*.png

.PHONY: all test run figures clean
