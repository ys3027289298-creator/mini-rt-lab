/*
 * renderer_edge.c -- the single translation unit that adapts the three shipped
 * modules to the checker's InvModel surface.
 *
 * Each module source is included *unmodified* exactly once, with only its
 * entry point renamed.  Their static functions are therefore visible to the
 * small forwarders in src/shim/ (each shim .c) only inside those TUs; the forwarders are
 * declared in inv_impl.h and collected into edge_model_* below.
 */
#include <string.h>

#include "inv_impl.h"

/* real_renderer / real_halton are the public aliases declared in
 * invariants.h; they are the same objects as edge_model_*, kept separate so the
 * header does not have to mention the adapter layer. */
const InvModel inv_real_renderer = {
    "shipped/halton2", edge_box_hit, edge_halton2, edge_ggx_D,
    edge_lambert, edge_cosine_pdf, edge_cosine_sample, edge_build
};

const InvModel inv_real_halton = {
    "shipped/halton", edge_box_hit, edge_halton, edge_ggx_D,
    edge_lambert, edge_cosine_pdf, edge_cosine_sample, edge_build
};

const char* edge_model_source_desc(void) {
    return "bvh/src/bvh.c, sample/src/sample.c, material/src/material.c (unmodified)";
}
