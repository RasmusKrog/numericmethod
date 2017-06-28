#include<stdlib.h>
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include<math.h>
#include<assert.h>
#include "lineq.h"

#define RND (double)rand()/RAND_MAX

int main(int argc, char const *argv[]){
	int m=5;
	int n=m;
	gsl_matrix* B = gsl_matrix_alloc(m,n);
	gsl_matrix* R = gsl_matrix_alloc(n,n);
	gsl_matrix* BT = gsl_matrix_alloc(n,m);
	gsl_matrix* I = gsl_matrix_alloc(n,n);
	gsl_matrix* Res = gsl_matrix_alloc(m,n);
	gsl_vector* b = gsl_vector_alloc(B->size1);
	gsl_matrix* E = gsl_matrix_alloc(B->size1,B->size2);
	gsl_vector* ei= gsl_vector_alloc(B->size1);
	printf("B*x=b, b is originally=\n" );
	for (int i=0; i<m; i++){
	for (int j = 0; j < n; j++){
	double Bij = RND;
	gsl_matrix_set(B,i,j,Bij);
	}
	gsl_vector_set(b,i,RND);
	}
	vector_print(b);
	printf("B is originally =\n");
	matrix_print(B);

// 1A
	fprintf(stderr, "A1\n");
	qr_gs_decomp(B, R);  //B-->Q
	printf("B=QR\n");
	printf("Upper triangular matrix R=\n"); // check R
	matrix_print(R);
	printf("Ortogonal matrix Q=\n"); //check Q
	matrix_print(B);
	double f = dot_prod(B,0,1);
	double ff = dot_prod(B,1,2);
	double fff = dot_prod(B,0,2);
	printf("check Q(0)*Q(1) = %g, should be 0\n", f);	
	printf("check Q(1)*Q(2) = %g, should be 0\n", ff);	
	printf("check Q(0)*Q(2) = %g, should be 0\n", fff);	
											gsl_matrix_transpose_memcpy(BT,B);
	matrix_prod(I,BT,B);
	printf("Check Q on normality Q^TQ =\n");
	matrix_print(I);

	printf("check QR = \n" );
	matrix_prod(Res,B,R);
	matrix_print(Res);

// 1a - part 2

	fprintf(stderr, "A2\n" );
	qr_gs_solve(B, R, b);
	printf("Solution to Bx=b, x=");
	vector_print(b);
	matrix_vector_prod(Res,b);
	printf("check Bx=\n");
	vector_print(b);
// 1B
	fprintf(stderr, "B\n" );
	for(int i=0; i< E->size2; i++){
	gsl_matrix_set(E,i,i,1);
	matrix_coloumn_get(E,ei,i);
	qr_gs_solve(B, R, ei);
	matrix_coloumn_set(E,ei,i); // E --> inv(A)
}
	printf("Ainv=\n");
	matrix_print(E);
	matrix_prod(I,B,R);
	matrix_prod(BT,E,I);
	printf("check Ainv*A=\n" );
	matrix_print(BT);
	gsl_vector_free(ei);
	gsl_matrix_free(E);
	gsl_matrix_free(I);
	gsl_matrix_free(BT);
	gsl_matrix_free(B);
	gsl_matrix_free(R);
	gsl_vector_free(b);
	gsl_matrix_free(Res);
										return 0;
}






