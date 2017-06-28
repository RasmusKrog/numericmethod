void qr_gs_decomp(gsl_matrix* A, gsl_matrix* R);
double dot_prod(gsl_matrix* A, int i, int j);
void matrix_coloumn_get(gsl_matrix* A, gsl_vector* b,int i);
void matrix_coloumn_set(gsl_matrix* A, gsl_vector* b,int i);
void matrix_print(gsl_matrix* A);
void matrix_prod(gsl_matrix* AB, gsl_matrix* A, gsl_matrix* B);
void qr_gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b);
void matrix_vector_prod(gsl_matrix* A, gsl_vector* b);
void vector_print(gsl_vector* b);
void inv_R(gsl_matrix* E, gsl_matrix* R);
double fitfunction(int i,double x);
void leastsquares(gsl_vector* b, gsl_matrix* R, gsl_vector* y, gsl_vector* x, gsl_vector* dy);

