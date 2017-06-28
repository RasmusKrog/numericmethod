#include<stdlib.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include<math.h>
#include "linearequation.h"

void inv_R(gsl_matrix* E, gsl_matrix* R){
	gsl_vector* ei= gsl_vector_alloc(E->size1);
	for(int i=0; i< E->size2; i++){
	gsl_matrix_set(E,i,i,1);
	matrix_coloumn_get(E,ei,i);
	qr_gs_solve(E, R, ei);
	matrix_coloumn_set(E,ei,i); // E --> inv(R)
	}
}

double fitfunction(int i,double x){
	switch(i){
	case 0: return 1.0/x; break;
	case 1: return 1.0; break;
	case 2: return x; break;
	default: {fprintf(stderr, "fitfunction: wrong i: %d\n",i); return NAN;}
	}
}

void leastsquares(gsl_vector* b, gsl_matrix* R, gsl_vector* y, gsl_vector* x, gsl_vector* dy){
	int n = x->size;
	int m = 3;
	// gsl_vector* b = gsl_vector_alloc(n);
	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_matrix* E = gsl_matrix_alloc(m,m);
	gsl_matrix* ET = gsl_matrix_alloc(m,m);
	//gsl_matrix* R = gsl_matrix_alloc(m,m);
	for (int i = 0; i < n; i++){
	gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
	for (int k = 0; k < m; k++){
	double f = fitfunction(k,gsl_vector_get(x,i));
	gsl_matrix_set(A,i,k,f/gsl_vector_get(dy,i));
	}
}
	qr_gs_decomp(A, R);;
	qr_gs_solve(A, R, b);
	inv_R(E, R);
	gsl_matrix_transpose_memcpy(ET,E);
	matrix_prod(R,E,ET);
}
