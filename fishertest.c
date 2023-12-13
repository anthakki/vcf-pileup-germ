
/* This file is subject to Mozilla Public License *
 * Copyright (c) 2016 Antti Hakkinen              */

#include "fishertest.h"
#include <assert.h>
#include <math.h> /* exp(), log() */
#include <stddef.h>

static
double
lnnchoosek(unsigned long n, unsigned long k)
{
	double y;
	unsigned long j;

	/*
	 * Computes the logarithm of a binomial coefficient n choose k
	 */

	/* Factor in the values */
	y = 0.;
	for (j = 1; j++ < n-k;)
		y -= log((double)j);
	for (j = k; j++ < n;)
		y += log((double)j);

	return y;
}

static
double
lnhygepdf(unsigned long x, unsigned long m, unsigned long n, unsigned long k)
{
	/* Compute logarithm of a hypergeometric pmf */
	return lnnchoosek(m,x) + lnnchoosek(n,k-x) - lnnchoosek(m+n,k);
}

static
double
lnhygepdfdn(double y, unsigned long x, unsigned long m, unsigned long n, unsigned long k)
{
	/* Compute lndhyper(x,m,n,k) given y := lndhyper(x+1,m,n,k) */
	return y + log(( (double)(x+1) / (m-x) * (n-k+x+1) / (k-x) ));
}

static
double
lnhygepdfup(double y, unsigned long x, unsigned long m, unsigned long n, unsigned long k)
{
	/* Compute lndhyper(x,m,n,k) given y := lndhyper(x-1,m,n,k) */
	return y + log(( (double)(m-x+1) / x * (k-x+1) / (n-k+x) ));
}

double
fishertest(const unsigned long *counts, int tail)
{
	assert(counts != NULL);

	/* We compute using log-probabilities anyway.. */
	return exp(fishertest_log(counts, tail));
}

static
double
log1pexp(double x)
{ return log1p(exp(x)); }

static
double
logspace_add(double x, double y)
{
	if (!(x < y))
		return x + log1pexp( y-x );
	else
		return y + log1pexp( x-y );
}

double
fishertest_log(const unsigned long *counts, int tail)
{
	static const double lnrelerr = 1e-7;

	unsigned long m, n, k, x11, lo, hi, x;
	double log_s, y;

	assert(counts != NULL);

	/* Get partial sums */  
	m = counts[0] + counts[1];
	n = counts[2] + counts[3];
	k = counts[0] + counts[2];
	x11 = counts[0];

	/* Get support for the pmf */
	lo = 0;
	if (k > n)
		lo = k-n;
	hi = k;
	if (k > m)
		hi = m;

	/* Compute */
	log_s = -HUGE_VAL;
	switch ((tail > 0) - (tail < 0))
	{
		case -1:
			/* Sum left tail */
			y = 0.;
			for (x = x11; x-- > lo;)
			{
				y = lnhygepdfdn(y, x, m, n, k);
				log_s = logspace_add(log_s, y);
			}
			break;

		case +1:
			/* Sum right tail */
			y = 0.;
			for (x = x11; x++ < hi;)
			{
				y = lnhygepdfup(y, x, m, n, k);
				log_s = logspace_add(log_s, y);
			}
			break;

		default:
			/* Sum left tail */
			y = 0.;
			for (x = x11; x-- > lo;)
			{
				y = lnhygepdfdn(y, x, m, n, k);
				if (y <= lnrelerr)
					log_s = logspace_add(log_s, y);
			}

			/* Sum right tail */
			y = 0.;
			for (x = x11; x++ < hi;)
			{
				y = lnhygepdfup(y, x, m, n, k);
				if (y <= lnrelerr)
					log_s = logspace_add(log_s, y);
			}

			break;

	} /* switch (tail) */

	/* Scale the result */
	y = lnhygepdf(x11, m, n, k);
	return y + log1pexp(log_s);
}
