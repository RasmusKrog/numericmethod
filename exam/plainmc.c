#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <gsl/gsl_rng.h>

void pseudo_random_vector(int dim, double *a, double *b, double *x, gsl_rng * r);

int plainmc(int d, double a[], double b[], double func(double*), int N, double *result, double *error){
	
	int i;
	double x[d];
	double sum1 = 0, sum2 = 0, f;
	double vol = 1;
for (int i = 0; i < d; i++)
	vol *= (b[i] - a[i]);
gsl_rng *r;
#pragma omp master
{
	r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, omp_get_num_threads());

for (i = 0; i < N; i++) {
	pseudo_random_vector(d, a, b, x, r);
	f = func(x);
	sum1 += f;
	sum2 += f * f;
	}
	gsl_rng_free(r);
	}
	double avr = sum1 / N, var = sum2 / N - avr * avr;
	*result = avr * vol;
	*error = sqrt(var / N) * vol;
return 0;
}
