#include "gsl/gsl_vector.h"

double lspline(int n, gsl_vector *x, gsl_vector *y, double z);
double lspline_integ(int n, gsl_vector *x, gsl_vector *y, double z);
typedef struct{int n; double *x, *y, *b, *c;} qspline;
qspline * qspline_alloc(int n, gsl_vector *x,gsl_vector* y);
double qspline_evaluate(qspline* s, double z);
double qspline_deriv(qspline *s, double z);
double qspline_integ(qspline *s, double z);
void qspline_free(qspline* s);

typedef struct{int n; double *x, *y, *b, *c, *d;} cspline;
cspline * cspline_alloc(int n, gsl_vector *x,gsl_vector* y);
double cspline_evaluate(cspline* s, double z);
double cspline_deriv(cspline *s, double z);
double cspline_integ(cspline *s, double z);
void cspline_free(cspline* s);
