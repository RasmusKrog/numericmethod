#include<stdlib.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include<math.h>
#include "lineq.h"

double dot_prod(gsl_matrix* A, int i, int j){
	double Akl=0;
	//fprintf(stderr, "No of rows=%zu\n", A->size1);
	for(int k=0; k < A->size1; k++){
		Akl += gsl_matrix_get(A,k,i)*gsl_matrix_get(A,k,j);
	//	fprintf(stderr, "prod %g for row %i\n",Akl,k);
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
		printf("%g ", gsl_matrix_get(A,i,j));
		}
	printf("\n");
	}
}

void vector_print(gsl_vector* b){
	for(int i=0;i < b->size; i++){
	printf("%g, ", gsl_vector_get(b,i));
	}
	printf("\n");
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

void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R){
	double Rii, Rij;
	gsl_vector* ai = gsl_vector_alloc(A->size1);
	//	gsl_vector* aj = gsl_vector_alloc(A->size1);
	int i, j;
	for(i=0; i< A->size2; i++){
	gsl_vector* aj = gsl_vector_alloc(A->size1);
	Rii = sqrt(dot_prod(A,i,i)); 
	gsl_matrix_set(R,i,i,Rii); //fill the R matrix
	matrix_coloumn_get(A,ai,i); //extract the vector ai from A
	gsl_vector_scale(ai,1/Rii); // scale ai --> ai/Rii = qi
	matrix_coloumn_set(A,ai,i);  // A(i) --> Q(i)
	for (j = i+1; j < A->size2; j++){
	matrix_coloumn_get(A,ai,i);	// reset ai again		
	Rij = dot_prod(A,i,j); // qi^T*aj
	gsl_matrix_set(R,i,j,Rij); // fill R matrix
	matrix_coloumn_get(A,aj,j); // get aj from A
	gsl_vector_scale(ai,Rij); // ai --> ai*Rij - but not saving into A
	gsl_vector_sub(aj,ai);   // aj --> aj - ai
	matrix_coloumn_set(A,aj,j); // A(j) --> Q(j)
}

gsl_vector_free(aj);
	}
	gsl_vector_free(ai);
//	gsl_vector_free(aj);
}

void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b){
	gsl_matrix* QT = gsl_matrix_alloc(Q->size2, Q->size1);
	gsl_matrix_transpose_memcpy(QT,Q);
	matrix_vector_prod(QT,b);
	for (int i = b->size-1; i >= 0; i--){
	double s = gsl_vector_get(b,i);
	for (int k = i+1; k < b->size; ++k){
	s -= gsl_matrix_get(R,i,k)*gsl_vector_get(b,k);
}
	gsl_vector_set(b,i,s/gsl_matrix_get(R,i,i));
}
}
