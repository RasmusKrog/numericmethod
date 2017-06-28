#include<stdlib.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include<math.h>
#include<assert.h>
#include "linearequation.h"

int main(int argc, char const *argv[]){
	gsl_vector* x = gsl_vector_alloc(10);
	gsl_vector* y = gsl_vector_alloc(10);
	gsl_vector* dy = gsl_vector_alloc(10);

	gsl_vector_set(x,0,0.100);
	gsl_vector_set(x,1,0.145);
	gsl_vector_set(x,2,0.211);
	gsl_vector_set(x,3,0.307);
	gsl_vector_set(x,4,0.447);
	gsl_vector_set(x,5,0.649);
	gsl_vector_set(x,6,0.944);
	gsl_vector_set(x,7,1.372);
	gsl_vector_set(x,8,1.995);
	gsl_vector_set(x,9,2.900);

	gsl_vector_set(y,0,12.644);
	gsl_vector_set(y,1,9.235);
	gsl_vector_set(y,2,7.377);
	gsl_vector_set(y,3,6.460);
	gsl_vector_set(y,4,5.555);
	gsl_vector_set(y,5,5.896);
	gsl_vector_set(y,6,5.673);
	gsl_vector_set(y,7,6.964);
	gsl_vector_set(y,8,8.896);
	gsl_vector_set(y,9,11.355);

	gsl_vector_set(dy,0,0.858);
	gsl_vector_set(dy,1,0.359);
	gsl_vector_set(dy,2,0.505);
	gsl_vector_set(dy,3,0.403);
	gsl_vector_set(dy,4,0.683);
	gsl_vector_set(dy,5,0.605);
	gsl_vector_set(dy,6,0.856);
	gsl_vector_set(dy,7,0.351);
	gsl_vector_set(dy,8,1.083);
	gsl_vector_set(dy,9,1.002);

	printf("# x y dy\n");
	for (int i=0; i<x->size; i++){
	double xi=gsl_vector_get(x,i);
	double yi=gsl_vector_get(y,i);
	double dyi=gsl_vector_get(dy,i);	
	printf("%g %g %g \n",xi,yi,dyi );
	}
	printf("\n\n");

	int n = x->size;
	int m = 3;
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_matrix* S = gsl_matrix_alloc(m,m);
	leastsquares(b, S, y, x, dy);
	//printf("c=\n");
	//vector_print(b);
	//printf("S=\n");
	//matrix_print(S);
	double c0 = gsl_vector_get(b,0);
	double c1 = gsl_vector_get(b,1);
	double c2 = gsl_vector_get(b,2);
	printf("\n\n");

	gsl_vector* ai = gsl_vector_alloc(b->size);
	gsl_vector* SA = gsl_vector_alloc(b->size);
	printf("# x f df\n");
	for (double xi = 0.09; xi < 3; xi+=0.03){
	double f0 = fitfunction( 0, xi);
	double f1 = fitfunction( 1,xi);
	double f2 = fitfunction( 2, xi);
	double f = c0*f0+c1*f1+c2*f2;
	gsl_vector_set(ai,0,1/xi);
	gsl_vector_set(ai,1,1);
	gsl_vector_set(ai,2,xi);
	gsl_vector_set(SA,0,1/xi);
	gsl_vector_set(SA,1,1);
	gsl_vector_set(SA,2,xi);
	matrix_vector_prod(S,SA);
	double df = 0;
	for (int i = 0; i < ai->size; i++){
	df+= gsl_vector_get(SA,i)*gsl_vector_get(ai,i);
	} 
	double fplus = sqrt(df) + f;
	double fminus = f- sqrt(df);
	printf("%g %g %g %g %g \n",xi, f, sqrt(df), fplus, fminus);
}

// part B

	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_vector_free(dy);
	gsl_vector_free(b);
	gsl_matrix_free(S);
return 0;
}

