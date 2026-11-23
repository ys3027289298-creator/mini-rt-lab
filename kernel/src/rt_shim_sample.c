/* rt_shim_sample.c -- exposes sample/src/sample.c:f_quad()
   without modifying or linking the module's main(). */
#define main sample_module_selftest_main
#include "../../sample/src/sample.c"
#undef main

double rt_shim_f_quad(double u) { return f_quad(u); }
