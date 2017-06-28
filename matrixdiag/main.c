#include "stdlib.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl_shortcut.h"
#include "math.h"
#define RND (double)rand()/RAND_MAX

int main(){
	int n=4;
	gsl_matrix* A = gsl_matrix_alloc(n,n);
	gsl_matrix* C = gsl_matrix_alloc(n,n);
	gsl_matrix* CT = gsl_matrix_alloc(n,n);
	gsl_matrix* V = gsl_matrix_alloc(n,n);
	gsl_matrix* D = gsl_matrix_alloc(n,n);
	gsl_vector* e = gsl_vector_alloc(n);
	gsl_matrix* AVD = gsl_matrix_alloc(n,n);
	for (int i = 0; i < n; ++i){
	for (int j = 0; j < n; ++j){
	gsl_matrix_set(C,i,j,RND);
	}
}

	gsl_matrix_transpose_memcpy(CT,C);
	matrix_prod(A,CT,C);
	printf("A original=\n");
	matrix_print(A);
	int h = jacobi(A, e, V);
	printf("A=\n");
	matrix_print(A);
	printf("V =\n");
	matrix_print(V);
	printf("diag(D) =\n");
	vector_print(e);
	printf("h=%i\n", h);

/*	for (int i = 0; i < D->size1; ++i){
	gsl_matrix_set(D,i,i,gsl_vector_get(e,i));
	} */

	for (int i = 0; i < n; ++i){
	for (int j = i; j < n; ++j){
	gsl_matrix_set(A,i,j,gsl_matrix_get(A,j,i));
	}
	}
	printf("A restored=\n");
	matrix_print(A);
	gsl_matrix_transpose_memcpy(D,V);
	matrix_prod(AVD,D,A);
	matrix_prod(D,AVD,V);
	printf("D=\n");
	matrix_print(D);
	gsl_matrix_free(A);
	gsl_matrix_free(V);
	gsl_vector_free(e);
return 0;
}
