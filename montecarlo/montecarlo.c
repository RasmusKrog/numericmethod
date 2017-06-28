#include<stdlib.h>
#include<math.h>
#include "gsl/gsl_vector.h"
#define RND ((double)rand()/RAND_MAX)

void randomx(gsl_vector* a, gsl_vector* b, gsl_vector* x);
void plainmc(gsl_vector* a, gsl_vector* b, double f(gsl_vector* x), int N, double* result, double* err);

void randomx(gsl_vector* a, gsl_vector* b, gsl_vector* x){
	double xi;
	for (int i=0; i < a->size; i++){
	xi = gsl_vector_get(a,i) + RND * (gsl_vector_get(b,i) - gsl_vector_get(a,i));
	gsl_vector_set(x,i,xi);
	}
}

// This section will calculate the integral using simple montecarlo method

void plainmc(gsl_vector* a, gsl_vector* b, double f(gsl_vector* x), int N, double* result, double* err){
	gsl_vector* x = gsl_vector_alloc(a->size);
	double V=1;
	int dim = a->size;
	for (int i=0; i<dim; i++){
	V *= (gsl_vector_get(b,i) - gsl_vector_get(a,i));
	}
	double sum1=0, sum2=0, fx;
	for (int i=0; i<N; i++){
	randomx(a,b,x);
	fx=f(x);
	sum1 += fx; //fx
	sum2 += fx*fx; //fx^2
	}
	double avr = sum1/N, var = sum2/N - avr*avr;
	*result = avr * V;
	*err = sqrt(var / N) * V;
	gsl_vector_free(x);
}
