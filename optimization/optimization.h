#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
 
double rosenbrock(gsl_vector* x);
double himmelblau(gsl_vector* x);
void dfrosenbrock(gsl_vector* x, gsl_vector* df);
void dfhimmelblau(gsl_vector* x, gsl_vector* df);
void Hrosenbrock(gsl_vector* x, gsl_matrix* H);
void Hhimmelblau(gsl_vector* x, gsl_matrix* H);
double vector_norm(gsl_vector* x);
void matrix_vector_prod(gsl_matrix* A, gsl_vector* b, gsl_vector* out);
void givens_qr(gsl_matrix* A);
void givens_qr_QTvec(gsl_matrix* QR, gsl_vector* v);
void givens_qr_solve(gsl_matrix* QR, gsl_vector* b, gsl_vector* x);
int newtonMinimization(double f(gsl_vector* x), void gradient(gsl_vector *x, gsl_vector* df), void hessian(gsl_vector* x,gsl_matrix* H), gsl_vector* xstart, double eps);

