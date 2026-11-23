# Top-level build for mini-rt-lab.
#
#   make test         # all three module self-tests + kernel acceptance suite
#   make modules      # just bvh/, sample/, material/ self-tests
#   make kernel       # just the kernel acceptance suite
#   make portability  # kernel results across compilers and -O levels
#   make clean        # remove build products (tracked sources are untouched)
#
# Each module also has its own Makefile: `make -C bvh test`.

MODULES = bvh sample material

.PHONY: test modules kernel portability clean

test: modules
	$(MAKE) -C kernel test

modules:
	@for m in $(MODULES); do \
	    echo "=== $$m ==="; \
	    $(MAKE) -s -C $$m clean >/dev/null 2>&1 || true; \
	    $(MAKE) -s -C $$m test || exit 1; \
	done

kernel:
	$(MAKE) -C kernel test

portability:
	$(MAKE) -C kernel portability

clean:
	@for m in $(MODULES); do $(MAKE) -s -C $$m clean >/dev/null 2>&1 || true; done
	$(MAKE) -s -C kernel clean >/dev/null 2>&1 || true
