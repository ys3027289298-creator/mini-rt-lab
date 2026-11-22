CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra

all: budget

budget: src/budget.c
	$(CC) $(CFLAGS) $< -o $@ -lm

test: budget
	./budget

clean:
	rm -f budget error_heatmap_uniform.ppm error_heatmap_opt.ppm
