
#include "gtbeta.h"
#include "jacobi_rule.h"
#include <alloca.h>
#include <assert.h>
#include <math.h>

static
void
unif_rule(double *x, double *y, size_t n)
{
	size_t i;

	/* Get rule */
	jacobi_rule(x, y, n, 0., 0.);

	/* Transform to (0, 1) */
	for (i = 0; i < n; ++i)
		x[i] = .5 + .5 * x[i];

#if 1   /* NB. scale is necessary for the scale of integrals (0/1 and 1/1) to match non-int (0/0) */
	/* Scale */
	for (i = 0; i < n; ++i)
		y[i] = .5 * y[i];
#endif
}

static
__attribute__((unused))
double
lbeta(double x, double y)
{ return lgamma(x) + lgamma(y) - lgamma(x+y); }

void
beta_rule(double *x, double *y, size_t n, double a, double b)
{
	size_t i;

	/* Get rule */
	jacobi_rule(x, y, n, b-1, a-1);

	/* Transform to (0, 1) */
	for (i = 0; i < n; ++i)
		x[i] = .5 + .5 * x[i];

#ifdef HAVE_SCALE
	/* Scale */
	{ double ys;
	ys = .5 * exp( -(a+b-2)*log(2) - lbeta(a, b) );
	for (i = 0; i < n; ++i)
		y[i] = ys * y[i];
	}
#endif
}

static
__attribute__((unused))
double
lfac(double k)
{ return lgamma(k+1); }

static
double
bino_logl(double a, double b, double p)
{
#ifdef HAVE_SCALE
	return lfac(a+b) - lfac(a) - lfac(b) + a * log(p) + b * log1p(-p);
#else
	return a * log(p) + b * log1p(-p);
#endif
}

static
double
logspace_add(double x, double y)
{
	if (!(x < y))
		return x + log1p(exp(y-x));
	else
		return y + log1p(exp(x-y));
}

static
void
dec2base(size_t *d, size_t x, size_t b, size_t n)
{
	size_t k;
	
	/* Convert to digits.. */
	for (k = 0; k < n; ++k)
	{
		d[k] = x % b;
		x /= b;
	}
}

void
gtbeta_22_jl(double *jl, const unsigned long *counts, size_t samples, double we, double me)
{
	size_t combs, j, eord, pord, k, i, c, *inds;
	unsigned long tot_dp, a, b;
	double *ex, *ey, *px, *py, *v, *v1, u, u1;

	assert(jl != NULL || samples < 1);
	assert(counts != NULL || samples < 1);
	assert(we > 0.);
	assert(me > 0.);

	/* Get number of genotype combinations */
	combs = samples > 0;
	for (j = 0; j < samples; ++j)
		combs *= 2*(2+1)/2;

	/* Get total depth */
	tot_dp = 0;
	for (j = 0; j < samples; ++j)
	{
		/* Get counts */
		b = ( &counts[2*j] )[0];
		a = ( &counts[2*j] )[1];

		/* Aggregate */
		tot_dp += a + b;
	}

	/* TODO: optimize quadrature order - can be a bit lower*/

	/* Make error quadrature */
	eord = tot_dp + 1;
	ex = (double *)alloca( eord*2 * sizeof(*ex) );
	ey = &ex[eord];
	beta_rule(ex, ey, eord, we*me, we*(1.-me));

	/* Allocate space */
	px = (double *)alloca( eord*2 * sizeof(*px) );
	py = &px[eord];
	v = (double *)alloca( eord*(2*(2+1)/2)*samples * sizeof(*v) );
	inds = (size_t *)alloca( samples * sizeof(*inds) );

	/* Integrate VAF distributions */
	for (j = 0; j < samples; ++j)
	{
		/* Get counts */
		b = ( &counts[2*j] )[0];
		a = ( &counts[2*j] )[1];

		/* Make VAF quadrature */
		pord = a + b + 1;
		unif_rule(px, py, pord);

		/* Get 0/0 model */
		{
			v1 = &( &( &v[eord*3*j] )[eord*0] )[0];
		for (i = 0; i < eord; ++i)
			v1[i] = bino_logl( a, b, ex[i] );
		}

		/* Get 0/1 model */
		{
			v1 = &( &( &v[eord*3*j] )[eord*1] )[0];
		{
			for (i = 0; i < eord; ++i)
				v1[i] = -HUGE_VAL;
		}
		for (k = 0; k < pord; ++k)
		{
			double px_k, dp, log_py_k;

			/* Precompute constants */
			px_k = .5 * px[k];  /* NB. 0/1 variants */
			dp = ( 1 - px_k ) - px_k;
			log_py_k = log( py[k] );

			/* Get model */
			for (i = 0; i < eord; ++i)
				v1[i] = logspace_add( v1[i], bino_logl( a, b, px_k + ex[i] * dp ) + log_py_k );
		}
		}

		/* Get 1/1 model */
		{
			v1 = &( &( &v[eord*3*j] )[eord*2] )[0];
		{
			for (i = 0; i < eord; ++i)
				v1[i] = -HUGE_VAL;
		}
		for (k = 0; k < pord; ++k)
		{
			double px_k, dp, log_py_k;

			/* Precompute constants */
			px_k = px[k];  /* NB. 1/1 variants */
			dp = ( 1 - px_k ) - px_k;
			log_py_k = log( py[k] );

			/* Get model */
			for (i = 0; i < eord; ++i)
				v1[i] = logspace_add( v1[i], bino_logl( a, b, px_k + ex[i] * dp ) + log_py_k );
		}
		}
	}

	/* Integrate the sample models together */
	for (c = 0; c < combs; ++c)
	{
		/* Get genotypes of this combination */
		dec2base(inds, c, 2*(2+1)/2, samples);

		/* Merge models */
		u = -HUGE_VAL;
		for (i = 0; i < eord; ++i)
		{
			/* Factor models */
			u1 = log( ey[i] );
			for (j = 0; j < samples; ++j)
				u1 += ( &( &v[eord*3*j] )[eord*inds[j]] )[i];

			/* Integrate over error */
			u = logspace_add( u, u1 );
		}

		/* Stash result */
		jl[c] = u;
	}
}

void
gtbeta_22_gl(double *gl, const double *jl, size_t samples)
{
	size_t combs, c, j, *inds, i;
	double s;

	assert(gl != NULL || samples < 1);
	assert(jl != NULL || samples < 1);

	/* Get number of genotype combinations */
	combs = samples > 0;
	for (j = 0; j < samples; ++j)
		combs *= 2*(2+1)/2;

	/* Allocate space */
	inds = (size_t *)alloca( samples * sizeof(*inds) );

	/* Zero accumulators */
	for (j = 0; j < samples; ++j)
		for (i = 0; i < 2*(2+1)/2; ++i)
			( &gl[2*(2+1)/2*j] )[i] = -HUGE_VAL;

	/* Scatter each probability */
	for (c = 0; c < combs; ++c)
	{
		/* Get genotypes for this combination */
		dec2base(inds, c, 2*(2+1)/2, samples);

		/* Scatter */
		for (j = 0; j < samples; ++j)
			( &gl[2*(2+1)/2*j] )[inds[j]] = logspace_add( ( &gl[2*(2+1)/2*j] )[inds[j]], jl[c] );
	}

	/* Normalize */
	for (j = 0; j < samples; ++j)
	{
			s = -HUGE_VAL;
		for (i = 0; i < 2*(2+1)/2; ++i)
			s = logspace_add( s, ( &gl[2*(2+1)/2*j] )[i] );
		for (i = 0; i < 2*(2+1)/2; ++i)
			( &gl[2*(2+1)/2*j] )[i] -= s;
	}
}

#if 0
void
gtbeta_22_gl_RR(double *gl, const double *jl, const int *samples)
{
	gtbeta_22_gl(gl, jl, *samples);
}
#endif

const char *
gtbeta_22_gt(const double *gl)
{
	static const char *const gt[] = { "0/0", "0/1", "1/1" };
	size_t p, i;

	assert(gl != NULL);

	/* Just find the best */ 
	p = 0;
	for (i = 1; i < 2*(2+1)/2; ++i)
		if (gl[i] > gl[p])
			p = i;

	return gt[p];
}
