#ifndef GTBETA_H_
#define GTBETA_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Get joint likelihoods */ 
void gtbeta_22_jl(double *jl, const unsigned long *counts, size_t samples, double we, double me, size_t at);
/* Zap genotype */
void gtbeta_22_zap(double *jl, size_t samples, size_t si, size_t gi);
/* Get genotypic likelihoods */
void gtbeta_22_gl(double *gl, const double *jl, size_t samples);
/* Get genotype */
const char *gtbeta_22_gt(const double *gl);

#ifdef __cplusplus
}
#endif

#endif /* GTBETA_H_ */
