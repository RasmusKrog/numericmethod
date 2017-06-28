#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "gsl/gsl_vector.h"
#define pi 4*atan(1)

void plainmc(gsl_vector* a, gsl_vector* b, double f(gsl_vector* x), int N, double* result, double* err);

// Here are some funny integrals to calculate :-)

double f(gsl_vector* x){
	return gsl_vector_get(x,0)*gsl_vector_get(x,1)*pow(gsl_vector_get(x,2),2);
}

double g(gsl_vector* x){
	return sin(gsl_vector_get(x,0)*gsl_vector_get(x,1));
}

double h(gsl_vector* x){
	return 1/(pow(M_PI,3)*(1-cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2))));
}

double errfun(gsl_vector* x){
	return gsl_vector_get(x,0)*gsl_vector_get(x,1);
} 

int main(){

// Exercise A: calculating some funnt integrals

	printf("# A, Calculate some interesting integrals:\n");
	int N=1e5, dim=3; double err, result;

	gsl_vector* a = gsl_vector_alloc(dim);
	gsl_vector* b = gsl_vector_alloc(dim);
	for (int i = 0; i < dim; ++i){
	gsl_vector_set(a,i,0);
	gsl_vector_set(b,i,1);
	}

	plainmc(a, b, f, N, &result, &err);
		
	printf("# Integral of f(x,y,z)=x*y*z^2 on [0,1]x[0,1]x[0,1]:\n");
	printf("# N\t= %d\n",N);
	printf("# result\t= %g\n",result);
	printf("# exact\t= %g\n",1/12.);
	printf("# estimated error\t= %g\n",err);
	gsl_vector_free(a);
	gsl_vector_free(b);

	dim=2;
	gsl_vector* a2 = gsl_vector_alloc(dim);
	gsl_vector* b2 = gsl_vector_alloc(dim);
	for (int i = 0; i < dim; ++i){
	gsl_vector_set(a2,i,0);
	gsl_vector_set(b2,i,M_PI);
	}

	plainmc(a2, b2, g, N, &result, &err);
		
	printf("\n # Integral of f(x,y)=sin(xy) on [0,pi]x[0,pi]:\n");
	printf("# N\t= %d\n",N);
	printf("# result\t= %g\n",result);
	printf("# exact\t= %g\n",2.90068);
	printf("# estimated error\t= %g\n",err);
	gsl_vector_free(a2);
	gsl_vector_free(b2);

	printf("# A, Calculate the integral:\n");
	dim=3;
	gsl_vector* a3 = gsl_vector_alloc(dim);
	gsl_vector* b3 = gsl_vector_alloc(dim);
	for (int i = 0; i < dim; ++i){
	gsl_vector_set(a3,i,0);
	gsl_vector_set(b3,i,M_PI);
	}
	plainmc(a3,b3,h,N,&result,&err);

	printf("\n #  Integral of f(x,y,z)=1/(pi^3*(1-cos(x)cos(y)cos(z))) on [0,pi]x[0,pi]x[0,pi]:\n");
	printf("#  N\t= %d\n",N);
	printf("# result\t= %g\n",result);
	printf("# exact\t= %g\n",1.3932039296856768591842462603255);
	printf("# estimated error\t= %g\n",err);
	gsl_vector_free(a3);
	gsl_vector_free(b3);

// Exercise B:

	printf("\n\n");
	printf("# Exercise B, data:\n");
	int M;
	dim=2;
	gsl_vector* a4 = gsl_vector_alloc(dim);
	gsl_vector* b4 = gsl_vector_alloc(dim);
	for (int i = 0; i < dim; ++i){
	gsl_vector_set(a4,i,0);
	gsl_vector_set(b4,i,1);
	}

	printf("# m=0, S=4\n");

	for(int i=100; i<10000	; i=i+100){
	M=i+100;
	plainmc(a4,b4,errfun,M,&result,&err);
	printf("%d %g\n",M,err);
	}

	printf("\n\n");
	printf("# m=1, S=0\n");
	for(int i=100; i<10000; i++){
	printf("%d %g\n", i, 0.22*1/sqrt(i));
	}
	gsl_vector_free(a4);
	gsl_vector_free(b4);

return 0;
}	
