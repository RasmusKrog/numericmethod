double dot_prod(gsl_matrix* A, int i, int j);
void matrix_coloumn_get(gsl_matrix* A, gsl_vector* b,int i);
void matrix_coloumn_set(gsl_matrix* A, gsl_vector* b,int i);
void matrix_print(gsl_matrix* A);
void matrix_prod(gsl_matrix* AB, gsl_matrix* A, gsl_matrix* B);
void matrix_vector_prod(gsl_matrix* A, gsl_vector* b);
void vector_print(gsl_vector* b);
int jacobi(gsl_matrix *A, gsl_vector* e, gsl_matrix* V);
int jacobi1(gsl_matrix *A, gsl_vector* e, gsl_matrix* V);
