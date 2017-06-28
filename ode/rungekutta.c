#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

void rkstep12(void f(int n, double x, gsl_vector*yx, gsl_vector*dydx),
	int n, double x, gsl_vector* yx, double h, gsl_vector* yh, gsl_vector* dy){
	gsl_vector* k0 = gsl_vector_alloc(n);
	gsl_vector* yt = gsl_vector_alloc(n);
	gsl_vector* k12 = gsl_vector_alloc(n);
	f(n,x,yx,k0);  
	for(int i=0;i<n;i++){
	gsl_vector_set(yt,i,gsl_vector_get(yx,i)+ gsl_vector_get(k0,i)*h/2);
	}
	f(n,x+h/2,yt,k12); 
	for(int i=0;i<n;i++){
	gsl_vector_set(yh,i,gsl_vector_get(yx,i)+gsl_vector_get(k12,i)*h);
	}
	for(int i=0;i<n;i++){
	gsl_vector_set(dy,i,(gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*h/2);
	}
	gsl_vector_free(k0);
	gsl_vector_free(yt);
	gsl_vector_free(k12);
			        
}
