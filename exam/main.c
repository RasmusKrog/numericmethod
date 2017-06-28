#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#define PI 3.1415926535

void quasi_random_vector(int dim, double* a, double* b, double* x);
void pseudo_random_vector(int dim, double* a, double* b, double* x, gsl_rng *r);

int plainmc(int d, double a[], double b[], double (*func)(double*), int N, double* result, double* error);

int lattice(int d, double a[], double b[], double (*func)(double*), int N, double* result, double* error);

int quasimc(int d, double a[], double b[], double (*func)(double*), int N, double* result, double* error);

double fun(double x[]);

// Just to check if main.c works. If not removed please ignore the following 7 lines
//int main(int argc, char** argv){
//      double x[] = {0.1, 0.2, 0.3};
//	double r;
//	r = fun(x);
//	printf("r%f",r);
//	return 0;
//}

double fun(double x[]){
	double r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	if(r<1) 
return r;
	else
return 0;
}

int main(int argc, char** argv){
	int N;
double exact=PI/8;
int d=3; double a[] ={0,0,0}; double b[]={PI,PI,PI};
double result,error;
N=(argc>1?atoi(argv[1]):1000000);

//Results of tests

plainmc(d,a,b,&fun,N,&result,&error);
	printf("--- pseudo-random: ---\n");
	printf("N\t= %d\n",N);
	printf("result\t= %g\n",result);
	printf("exact\t= %g\n",exact);
	printf("estimated error\t= %g\n",error);
	printf("actual error\t= %g\n",fabs(result-exact));

lattice(d,a,b,&fun,N,&result,&error);
	printf("--- lattice: ---\n");
	printf("N\t= %d\n",N);
	printf("result\t= %g\n",result);
	printf("exact\t= %g\n",exact);
	printf("estimated error\t= %g\n",error);
	printf("actual error\t= %g\n",fabs(result-exact));

quasimc(d,a,b,&fun,N,&result,&error);
	printf("--- halton: ---\n");
	printf("N\t= %d\n",N);
	printf("result\t= %g\n",result);
	printf("exact\t= %g\n",exact);
	printf("estimated error\t= %g\n",error);
	printf("actual error\t= %g\n",fabs(result-exact));

return 0;
}
