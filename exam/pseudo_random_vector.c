#include <stdlib.h>
#include <gsl/gsl_rng.h>
#define RND ((double)rand()/RAND_MAX)

//Pseudo-random sampling creates regions with high density of points and other regions with low density.

void pseudo_random_vector(int dim,double* a,double* b, double* x, gsl_rng *r){
	
	for(int i=0;i<dim;i++) x[i]=a[i]+gsl_rng_uniform(r)*(b[i]-a[i]);
}

// gsl_rng_uniform returns a double precision floating point number uniformly distributed in the range from 0 to 1
