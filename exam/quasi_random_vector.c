#include <stdlib.h>
#include <tgmath.h>
#define fracl(x) ((x)-floorl(x))
#define E (long double)2.7182818284
#define PI (long double)3.1415926535

//Quasi-random sampling distribute points in a highly correlated manner with a specific requirement of low discrepancy.

void quasi_random_vector(int dim, double *a, double *b, double *x, int sec) {
	static long int n1=0, n2=0; //global variable
	static long double *alpha, *beta;
	
	if(n1==0){
		alpha=(long double*)malloc(dim*sizeof(long double));
		for(int i=0;i<dim;i++) alpha[i]=fracl(cbrtl(11*i+E)); //cbrtl is a long double arg
	}
					
	if(n2==0){
		beta =(long double*)malloc(dim*sizeof(long double));
		for(int i=0;i<dim;i++) beta[i]=fracl(cbrtl(13*i+PI));
	}
						
	if(sec == 0){
		n1++;
		for(int i=0;i<dim;i++) x[i]=a[i]+fracl(n1*alpha[i])*(b[i]-a[i]);
	}
						
	else{
		n2++;
		for(int i=0;i<dim;i++) x[i]=a[i]+fracl(n2*beta [i])*(b[i]-a[i]);
	}
}
