#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include "optimization.h"

int main(){
	gsl_vector* xstart = gsl_vector_alloc(2);
	double eps = 0.0001;
	int steps;
	steps = newtonMinimization(rosenbrock, dfrosenbrock, Hrosenbrock, xstart, eps);
	printf("The minimum in the Rosenbrock function is:\n");
	printf("x=%g\n", gsl_vector_get(xstart,0));
	printf("y=%g\n", gsl_vector_get(xstart,1));
	printf("The number of steps is: %i\n", steps);
	gsl_vector_set(xstart,0,10);
	gsl_vector_set(xstart,1,-5);
	steps =newtonMinimization(himmelblau, dfhimmelblau, Hhimmelblau, xstart, eps);
	printf("The minimum in the Himmelblau function is:\n");
	printf("x=%g\n", gsl_vector_get(xstart,0));
	printf("y=%g\n", gsl_vector_get(xstart,1));
	printf("The number of steps is: %i\n", steps);

return 0;
}
