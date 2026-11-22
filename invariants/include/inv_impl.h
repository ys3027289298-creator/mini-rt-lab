/*
 * inv_impl.h -- internal interface between the golden reference, the mutants
 * and the shipped-module adapters.  Not part of the public checker API.
 */
#ifndef INV_IMPL_H
#define INV_IMPL_H

#include "invariants.h"

int    golden_box_hit(const double lo[3], const double hi[3],
                      const double o[3], const double d[3],
                      double tmax, double* t);
double golden_halton(unsigned i, unsigned base);
double golden_ggx_D(double nh, double rough);
double golden_lambert(double albedo);
double golden_cosine_pdf(double cos_theta);
void   golden_cosine_sample(double u1, double u2, double out[3]);
int    golden_build(InvTree* out);

/* adapter: the shipped modules, compiled unmodified into renderer_edge.c */
int    edge_box_hit(const double lo[3], const double hi[3],
                    const double o[3], const double d[3],
                    double tmax, double* t);
int    edge_build(InvTree* out);
double edge_halton(unsigned i, unsigned base);
double edge_halton2(unsigned i, unsigned base);
double edge_ggx_D(double nh, double rough);
double edge_lambert(double albedo);
double edge_cosine_pdf(double cos_theta);
void   edge_cosine_sample(double u1, double u2, double out[3]);

extern const InvModel edge_model_halton2;
extern const InvModel edge_model_halton;
extern const char* edge_model_source_desc(void);

/* Purity state of the three modules' own self-test counters. */
int edge_bvh_failures(void);
int edge_sample_failures(void);
int edge_material_failures(void);

#endif
