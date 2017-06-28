#include<stdlib.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include<math.h>
#include "gsl_shortcut.h"

double dot_prod(gsl_matrix* A, int i, int j){
	double Akl=0;
	for(int k=0; k < A->size1; k++){
	Akl += gsl_matrix_get(A,k,i)*gsl_matrix_get(A,k,j);
	}
	return Akl;
}

void matrix_coloumn_get(gsl_matrix* A, gsl_vector* b, int i){
		
	int k;
	for(k=0; k<A->size1; k++){
	gsl_vector_set(b,k,gsl_matrix_get(A,k,i));
	}
}

void matrix_coloumn_set(gsl_matrix* A, gsl_vector* b, int i){
	int k;
	for(k=0; k<A->size1; k++){
	gsl_matrix_set(A,k,i,gsl_vector_get(b,k));
	}
}

void matrix_print(gsl_matrix* A){
	for(int i=0;i < A->size1; i++){
	for(int j=0; j< A->size2; j++){
	printf("%g \t ", gsl_matrix_get(A,i,j));
	}
	printf("\n");
	}
	printf("\n");
}

void vector_print(gsl_vector* b){
	for(int i=0;i < b->size; i++){
	printf("%g \t ", gsl_vector_get(b,i));
	}
	printf("\n\n");
}

void matrix_prod(gsl_matrix* AB, gsl_matrix* A, gsl_matrix* B){
	for (int k = 0; k < A->size1; k++){
	for (int l = 0; l < B->size2; l++){
	double P=0;
	for(int i=0;i < A->size2; i++){
	P+=gsl_matrix_get(A,k,i)*gsl_matrix_get(B,i,l);
	}
	gsl_matrix_set(AB,k,l,P);
	}
	}
}

void matrix_vector_prod(gsl_matrix* A, gsl_vector* b){
	gsl_vector* c = gsl_vector_alloc(A->size1);
	for (int i = 0; i < A->size1; i++){
	double P=0;
	for(int j=0;j < b->size; j++){
	P+=gsl_matrix_get(A,i,j)*gsl_vector_get(b,j);
	}
	gsl_vector_set(c,i,P);
	}
	for (int i = 0; i < c->size; ++i){
	gsl_vector_set(b,i,gsl_vector_get(c,i));
	}
}
