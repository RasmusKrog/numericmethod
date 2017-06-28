#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

double adaptint(double f(double), double a, double b, double acc, double eps, double *err, double QTrue);

int main(int argc, char const *argv[]){
	int calls = 0;
	double a = 0,b = 1, err = 0, acc = 0.001, eps = 0.001;
	double Q;

	double sqrt_x(double x){
	calls++;
	//True Q-value = 2/3;
	return sqrt(x);
	}

double
	inv_sqrt_x(double x){
	calls++;
	//True Q-value = 2;
	return 1/sqrt(x);
	}

double
	ln_sqrt_x(double x){
	calls++;
	//True Q-value= = -4;
	return log(x)/sqrt(x);
	}

double
	f_A2(double x){
	calls++;
	//True Q-value = -pi;
	return 4*sqrt(1-(1-x)*(1-x));
	}

double 
	f_B(double x){
	calls++;
	//True Q-value = sqrt(pi), double y=x/(1-pow(x,2));
	return exp(-pow(x,2));
	}

// Exercise A

	fprintf(stderr,"Exercise A, part 1: \n The code is tested on following test functions:\n");

	double Qt = 2.0/3.0;
	Q = adaptint(sqrt_x,a,b,acc,eps,&err, Qt);
	printf("Int_theo(sqrt(x)) = %g, Q = %lg, error = %lg, calls = %i\n", Qt,Q,err,calls);

	Q = adaptint(inv_sqrt_x,a,b,acc,eps,&err, 2);
	printf("Int_theo(1/sqrt(x)) = %g, Q = %lg, error = %lg, calls = %i\n", 2.0,Q,err,calls);

	Q = adaptint(ln_sqrt_x,a,b,acc,eps,&err, -4);
	printf("Int_theo(ln(x)/sqrt(x)) = %g, Q = %lg, error = %lg, calls = %i\n", -4.0,Q,err,calls);

	printf("Exercise A, part 2: \n Following function is integrated:\n");

	Qt = M_PI;
	printf("%g\n", Qt);

	Q = adaptint(f_A2,a,b,acc,eps,&err, Qt);
	printf("4*sqrt(1-(1-x)*(1-x)) = %g, Q = %lg, error = %lg, and the number of calls = %i\n", Qt,Q,err,calls);

// Exercise B

	fprintf(stderr,"Exercise B\n");
	a = -INFINITY;
	b = INFINITY;

	Qt = sqrt(M_PI);
	Q = adaptint(f_B,a,b,acc,eps,&err, Qt);
	printf("exp(-x^2) = %g, Q = %lg, error = %lg, and the number of calls = %i\n", Qt,Q,err,calls);
				
	Qt = sqrt(M_PI)/2;
	Q = adaptint(f_B2,a,b,acc,eps,&err, Qt);
	printf("x^2*exp(-x^2) = %g, Q = %lg, error = %lg, and the number of calls = %i\n", Qt,Q,err,calls);

	a = -INFINITY;
	b = 0;
									
	Qt = sqrt(M_PI)/2;
	Q = adaptint(f_B,a,b,acc,eps,&err, Qt);
	printf("exp(-x²) = %g, Q = %lg, error = %lg, and the number of calls = %i\n", Qt,Q,err,calls);

	a = 0;
	b = INFINITY;
	
	Qt = sqrt(M_PI)/2;
	Q = adaptint(f_B,a,b,acc,eps,&err, Qt);
	printf("exp(-x²) = %g, Q = %lg, error = %lg, and the number of calls = %i\n", Qt,Q,err,calls);

return 0;
}
