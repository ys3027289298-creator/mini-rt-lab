# Minimal in-tree test entry. No third-party test library: the test C
# files #include the module source directly and rename its demo main.
CC ?= gcc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra
REPO := ..

all: bvh_test material_test

bvh_test: bvh_test.c $(REPO)/bvh/src/bvh.c
	$(CC) $(CFLAGS) -DBVH_SRC='"$(REPO)/bvh/src/bvh.c"' bvh_test.c -o $@ -lm

material_test: material_test.c $(REPO)/material/src/material.c
	$(CC) $(CFLAGS) -DMATERIAL_SRC='"$(REPO)/material/src/material.c"' material_test.c -o $@ -lm

# The demo programs' success-path output must stay byte-identical.
check-baseline:
	@$(MAKE) -s -C $(REPO)/bvh
	@$(MAKE) -s -C $(REPO)/material
	@$(MAKE) -s -C $(REPO)/sample
	@$(REPO)/bvh/bvh > /tmp/bvh.out 2>&1
	@$(REPO)/material/material > /tmp/material.out 2>&1
	@$(REPO)/sample/sample > /tmp/sample.out 2>&1
	@cd /tmp && sha256sum -c $(CURDIR)/baseline/bvh.sha256 \
	  $(CURDIR)/baseline/material.sha256 $(CURDIR)/baseline/sample.sha256

test: bvh_test material_test
	@./bvh_test
	@./material_test

mutation:
	@python3 -I mutate_test.py

clean:
	rm -f bvh_test material_test

.PHONY: all test check-baseline mutation clean
