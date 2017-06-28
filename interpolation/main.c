#include "spline.h"
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<assert.h>

#define RND (double)rand()/RAND_MAX
/*
struct qspline{int n; double *x, *y, *b, *c;};
struct qspline * qspline_alloc(int n, double *x,double* y);
double qspline_evaluate(struct qspline* s, double z);
void qspline_free(struct qspline* s);
*/
int main(int argc, char const *argv[]){
	int n = atof(argv[1]);
	//int n=6;
	gsl_vector *y = gsl_vector_alloc(n);
	gsl_vector *x = gsl_vector_alloc(n);

	/*make data with n data points*/
	fprintf(stderr, "The original data points:\n");
	for (int i=0; i < x->size; i++){
		gsl_vector_set(x,i,M_PI/n*2*i);
		gsl_vector_set(y,i,sin(gsl_vector_get(x,i)));
		printf("%g %g\n", gsl_vector_get(x,i),gsl_vector_get(y,i));
	}
	printf("\n\n");
	double dz=0.001;

	qspline *qs = qspline_alloc(n,x,y);
	cspline *cs = cspline_alloc(n,x,y);
	for (double z=gsl_vector_get(x,0); z<gsl_vector_get(x,n-1); z+=dz){
		double linearData = lspline(n,x,y,z);
		double linearInt = lspline_integ(n,x,y,z);		
		double squareData = qspline_evaluate(qs,z);
		double squareInt = qspline_integ(qs,z);
		double squareDer = qspline_deriv(qs,z);
		double cubicData = cspline_evaluate(cs,z);
		double cubicInt = cspline_integ(cs,z);
		double cubicDer = cspline_deriv(cs,z);
		printf("%g %g %g %g %g %g %g %g %g\n", z,linearData,linearInt,squareData,squareInt,squareDer,cubicData,cubicInt,cubicDer);
	}
	fprintf(stderr, "after data generating loop\n");

	gsl_vector_free(x);
	gsl_vector_free(y);
	qspline_free(qs);
	cspline_free(cs);
	return 0;
}





