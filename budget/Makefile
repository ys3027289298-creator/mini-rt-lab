CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: budget

budget: src/budget.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: budget
	./budget

report: budget
	./budget report 1048576 4 1

csv: budget
	./budget csv 1048576 4 1

plot: budget
	python3 tools/plot_budget.py 1048576 4 1 ../artifacts

empty: budget
	./budget empty

clean:
	rm -f budget
