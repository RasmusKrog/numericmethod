#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_vector.h"

void rkstep12(void f(int n, double x, gsl_vector*yx, gsl_vector*dydx),
	int n, double x, gsl_vector* yx, double h, gsl_vector* yh, gsl_vector* dy);

int ode_driver(void f(int n,double x,gsl_vector*y,gsl_vector*dydx),
	int n,gsl_vector* xlist, gsl_matrix* ylist,
	double b,double h,double acc,double eps,int max){
int i,k=0; 
	double x,s,err,normy,tol;
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* yh = gsl_vector_alloc(n);
	gsl_vector* dy = gsl_vector_alloc(n);
	double a = gsl_vector_get(xlist,0);

while(gsl_vector_get(xlist,k)<b){
	x = gsl_vector_get(xlist,k);
	for (int i = 0; i < n; ++i)    {
	gsl_vector_set(y,i,gsl_matrix_get(ylist,k,i)); 
}
	if(x+h>b) h=b-x;
	ode_stepper(f,n,x,y,h,yh,dy);
	s=0; 
	for(i=0;i<n;i++){ 
	s+= pow(gsl_vector_get(dy,i),2); 
}
	err = sqrt(s);    
	s=0; 
	for(i=0;i<n;i++){
	s+= pow(gsl_vector_get(yh,i),2); 
}
	normy = sqrt(s);  
	tol = (normy*eps + acc) * sqrt(h/(b-a));
	if(err<tol){ 
       	k++; 
	if(k>max-1) return -k; 
	gsl_vector_set(xlist,k,x+h); 
for (int i = 0; i < n; ++i)      {
	gsl_matrix_set(ylist,k,i,gsl_vector_get(yh,i));
}
}
	if(err>0) h*=pow(tol/err,0.25)*0.95; else h*=2;
	} /* end while here */

	gsl_vector_free(y);
	gsl_vector_free(yh);
	gsl_vector_free(dy);  
return k+1; 
} /* This will return the amount of entries in the x- and ylist */
