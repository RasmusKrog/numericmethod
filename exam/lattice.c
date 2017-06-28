#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

void quasi_random_vector (int dim,double* a,double* b,double* x, int sec);

int lattice(
	int d, double a[], double b[], double (*func)(double*),
	int N, double* result, double* error)
{
	double x[d];
	double vol=1; for(int i=0;i<d;i++) vol*=(b[i]-a[i]);
	double sum1=0, sum2=0, f; int i;
	N /= 2;
	quasi_random_vector(d,a,b,x,0);
	quasi_random_vector(d,a,b,x,1);
#pragma omp parallel private(i,x,f) //Specifies that each thread should have its own instance of a variable.
{
#pragma omp sections //The omp sections directive distributes work among threads bound to a defined parallel region
	{
#pragma omp section
	for(i=0;i<N;i++){
	quasi_random_vector(d,a,b,x,0); f=(*func)(x); sum1+=f;
	}
#pragma omp section
	for(i=0;i<N;i++){
	quasi_random_vector(d,a,b,x,1); f=(*func)(x); sum2+=f;
	}
}}
	*result=(sum1+sum2)/2/N*vol; *error=fabs(sum1-sum2)/2/N*vol;
return 0;
}
