#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

double adaptint(double f(double), double a, double b, double acc, double eps, double *err, double QTrue);
double adaptint24(double f(double), double a, double b, double acc, double eps, double f2, double f3, int nrec,int inf_inf,double a0, double b0);
double f_inf_inf(double f(double),double x,double a0, double b0, int inf_inf);


double f_inf_inf(double f(double),double x, double a0, double b0, int inf_inf){
	if (inf_inf == 1){
	double y = x/(1-pow(x,2));
	return f(y)*(1 + 1/pow(x,2))*pow(y,2);
	}
	if(inf_inf == 2){
	double y = x/(1-x);
	return f(y+a0)*pow(y,2)/pow(x,2);
	}
	if (inf_inf == -2) {
	double y = x/(1+x);
	return f(y+b0)*pow(y,2)/pow(x,2);
	}
	if (inf_inf == 0){
	return f(x);
	}
	else{
	return 0;
	}
}

double adaptint(double f(double), double a, double b, double acc, double eps, double *err, double QTrue){
	double a0=a, b0=b;
	int inf_inf = 0;
	if ((a==-INFINITY && b==INFINITY))	{
	a=-1;
	b=1;
	inf_inf = 1;
	}
	if ((a!=-INFINITY && b==INFINITY))	{
	a=0;
	b=1;
	inf_inf = 2;
	}
	if ((a==-INFINITY && b!=INFINITY))	{
	a=-1;
	b=0;
	inf_inf = -2;
	fprintf(stderr, "a=%g,b=%g,a0=%g,b0=%g, index = %i\n", a,b,a0,b0,inf_inf);
	}
	printf("inf %i\n",inf_inf );
	double f2 = f_inf_inf(f,a + 2*(b-a)/6,a0,b0,inf_inf);
	double f3 = f_inf_inf(f,a + 4*(b-a)/6,a0,b0,inf_inf); 
	int nrec = 0;
	double Q = adaptint24(f,a,b,acc,eps,f2,f3,nrec,inf_inf,a0,b0);
	*err= (fabs(Q - QTrue));
	return Q;
}

double adaptint24(double f(double), double a, double b, double acc, double eps, double f2, double f3, int nrec, int inf_inf, double a0, double b0){
	assert(nrec<1000000);

	double f1 = f_inf_inf(f,a + (b-a)/6,a0,b0,inf_inf);
	double f4 = f_inf_inf(f,a + 5*(b-a)/6,a0,b0,inf_inf);
	double Q = (2*f1 + f2 + f3 + 2*f4)/6*(b-a);
	double q = (f1 + f4 + f3 + f2)/4*(b-a);
	double tol = acc + eps*fabs(Q);
	double error = fabs(Q-q);
	if (error < tol) return Q;
	else{
	double Q1 = adaptint24(f,a,(a+b)/2,acc/sqrt(2.0),eps,f1,f2,nrec+1, inf_inf,a0,b0);
	double Q2 = adaptint24(f,(a+b)/2,b,acc/sqrt(2.0),eps,f3,f4,nrec+1, inf_inf,a0,b0);
	return Q1+Q2;
	}
}
	
