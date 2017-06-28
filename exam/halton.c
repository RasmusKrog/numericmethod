#include<stdio.h>
#include<assert.h>
#include<math.h>

// A Halton sequence is a generalizaatin of the van der Corput sequence to d-dimensional spaces where a van der Corput sequence is a low-discrepancy sequence over the unit interval.

double corput(int n, int base){
	double q=0, bk=(double)1/base;
	while(n>0) { q += (n % base)*bk; n /= base; bk /= base; }
return q;
}

void halton(int n, int d, double *a, double *b, double *x){
	int dmax=20; assert(d <= dmax);
	int base[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79}; //base consisting of prime numbers
for(int i=0;i<d;i++)x[i]=a[i]+corput(n+1,base[i])*(b[i]-a[i]);
	}

int quasimc(int d, double a[], double b[], double f(double*),
	int N, double* result, double* error){
	
	double x[d];
	double vol=1; for(int i=0;i<d;i++) vol*=(b[i]-a[i]);
	double sum1=0, sum2=0;
	for(int i=0;i<N/2;i++){
	halton(i,d,a,b,x); sum1+=f(x);
	}
	
	for(int i=N/2;i<N;i++){
	halton(i,d,a,b,x); sum2+=f(x);
	}
	
	*result=(sum1+sum2)/N*vol; *error=fabs(sum1-sum2)/N*vol;
return 0;
}

