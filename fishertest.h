#ifndef FISHERTEST_H_
#define FISHERTEST_H_

#ifdef __cplusplus
extern "C" {
#endif

/* Performs a Fisher's exact test on a 2-by-2 contingency table */
double fishertest(const unsigned long *counts, int tail);
/* As above, returns logarithmic p-value */
double fishertest_log(const unsigned long *counts, int tail);

#ifdef __cplusplus
}
#endif

#endif /* FISHERTEST_H_ */
