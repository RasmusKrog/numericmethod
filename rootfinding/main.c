#include<stdlib.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include<math.h>
#include<assert.h>

int newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart, double dx, double epsilon);

void funA(gsl_vector* x, gsl_vector* fx){
	int A = 10000;
	gsl_vector_set(fx,0,A*gsl_vector_get(x,0)*gsl_vector_get(x,1)-1);
	gsl_vector_set(fx,1,exp(-gsl_vector_get(x,0))*exp(-gsl_vector_get(x,1))-1-1/A);	
}

int main(){
	int n = 2;
	double epsilon = 0.00001;
	double dx = 0.001;
	gsl_vector* xstart = gsl_vector_alloc(n);
	gsl_vector_set(xstart,0,0);
	gsl_vector_set(xstart,1,9);

	int steps = newton(funA, xstart, dx,epsilon);
	printf("x0=%g\n",gsl_vector_get(xstart,0));
	printf("x1=%g\n",gsl_vector_get(xstart,1));
	printf("steps = %i\n", steps);
return 0;
}
