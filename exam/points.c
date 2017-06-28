#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>

void quasi_random_vector(int n,double *a,double *b,double *x,int sec);
void pseudo_random_vector(int d,double a[],double b[],double x[],gsl_rng *r);
void halton(int n, int d, double *a, double *b, double *x);

int main(int argc, char** argv){

	int N=1000;
	int d=2; double a[] ={0,0}; double b[]={1,1}; double x[]={0,0};

	printf("# m=0 S=4\n");

if(strcmp(argv[1],"quasi")==0){ //compares the string pointed to
	for(int i=0;i<N;i++){
	quasi_random_vector(d,a,b,x,0);
	printf("%g %g\n",x[0],x[1]);
	}
}

else if(strcmp(argv[1],"halton")==0){
	for(int i=0;i<N;i++){
	halton(i,d,a,b,x);
	printf("%g %g\n",x[0],x[1]);
	}
}
	
else if(strcmp(argv[1],"pseudo")==0){
	gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
	for(int i=0;i<N;i++){
	pseudo_random_vector(d,a,b,x,r);
	printf("%g %g\n",x[0],x[1]);
	}
}

else printf("usage: %s quasi[pseudo]\n",argv[0]);

return 0;
}
