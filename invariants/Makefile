CC      ?= gcc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra
INCDIR   = include
SRCDIR   = src
BIN      = invcheck

# The three module sources are compiled UNMODIFIED into the observer shims.
# Editing them is never necessary and would defeat the point of the exercise.
RT_SRCS  = $(SRCDIR)/shim/edge_bvh.c \
           $(SRCDIR)/shim/edge_sample.c \
           $(SRCDIR)/shim/edge_material.c
OWN_SRCS = $(SRCDIR)/golden.c $(SRCDIR)/mutants.c $(SRCDIR)/inv_checks.c \
           $(SRCDIR)/renderer_edge.c $(SRCDIR)/inv_main.c

.PHONY: all audit validate clean

all: $(BIN)

$(BIN): $(RT_SRCS) $(OWN_SRCS) $(INCDIR)/invariants.h $(INCDIR)/inv_impl.h
	$(CC) $(CFLAGS) -I$(INCDIR) $(RT_SRCS) $(OWN_SRCS) -o $@ -lm

# The two halves, reported separately on purpose:
#   validate -> the checker is trustworthy (no side effects, no false
#               positives on two known-good controls, every planted violation
#               caught).  Exit 0 required before trusting any audit result.
#   audit    -> run the invariant classes on the shipped renderer.  Exit 1 here
#               is the *finding* (a violated invariant), not a checker fault.
validate: $(BIN)
	./$(BIN) --render-check
	./$(BIN) --calibrate

audit: $(BIN)
	./$(BIN)

# Convenience: validate first; only audit if the checker itself is sound.
.PHONY: check
check: validate
	@echo "--- checker validated; running audit (exit 1 here = renderer violation) ---"
	-./$(BIN)

.PHONY: clean
clean:
	rm -f $(BIN)
