#ifndef JACOBI_RULE_H_
#define JACOBI_RULE_H_

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Gauss-Jacobi quarature */
void jacobi_rule(double *x, double *y, size_t n, double a, double b);

#ifdef __cplusplus
}
#endif

#endif /* JACOBI_RULE_H_ */
