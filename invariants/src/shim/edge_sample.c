/*
 * edge_sample.c -- observer shim for sample/src/sample.c (compiled unmodified).
 * Pure forwarders to the module's static functions.
 */
#define main sample_module_main
#define failures sample_module_failures
#include "../../sample/src/sample.c"
#undef main
#undef failures

#include "inv_impl.h"

double edge_halton(unsigned i, unsigned base)  { return halton(i, base); }
double edge_halton2(unsigned i, unsigned base) { return halton2(i, base); }

int edge_sample_failures(void) { return sample_module_failures; }
