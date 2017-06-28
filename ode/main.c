Last login: Mon May 22 15:19:13 on ttys000
Rasmus-Krogs-MacBook-Pro:rootfinding rasmuspedersen$ cd 
Rasmus-Krogs-MacBook-Pro:~ rasmuspedersen$ locate laura_feuerhake

WARNING: The locate database (/var/db/locate.database) does not exist.
To create the database, run the following command:

  sudo launchctl load -w /System/Library/LaunchDaemons/com.apple.locate.plist

  Please be aware that the database can take some time to generate; once
  the database has been created, this message will no longer appear.

  Rasmus-Krogs-MacBook-Pro:~ rasmuspedersen$ cd Downloads/laura_feuerhake-numeric-cc88470d6ae9
  Rasmus-Krogs-MacBook-Pro:laura_feuerhake-numeric-cc88470d6ae9 rasmuspedersen$ ls
  ODE		hello		leastsquaresfit	matrixdiag	roots
  f		interpolation	linearequations	optim
  Rasmus-Krogs-MacBook-Pro:laura_feuerhake-numeric-cc88470d6ae9 rasmuspedersen$ cd ODE/
  Rasmus-Krogs-MacBook-Pro:ODE rasmuspedersen$ ls
  Makefile	main.c		ode_driver.c	plot.png
  data		main.o		ode_driver.o	rkstep12.c
  main		ode.h		plot.gpi	rkstep12.o
  Rasmus-Krogs-MacBook-Pro:ODE rasmuspedersen$ cd Makefile 
  -bash: cd: Makefile: Not a directory
  Rasmus-Krogs-MacBook-Pro:ODE rasmuspedersen$ vim Makefile 
  Rasmus-Krogs-MacBook-Pro:ODE rasmuspedersen$ vim main.c




#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

  int ode_driver(
	void f(int n, double x, gsl_vector*y, gsl_vector*dydx),
	int n, gsl_vector* xlist, gsl_matrix* ylist,
	double b, double h, double acc, double eps, int max);

void f(int n, double x, gsl_vector* y, gsl_vector* dydx){
	gsl_vector_set(dydx,0, gsl_vector_get(y,1));
	gsl_vector_set(dydx,1, -1*gsl_vector_get(y,0));  // sinus fkt
return;
i}

int main(){
	int n=2;
	int max=1000;
	// Allocate xlist and ylist: xlist is containing times in each loop, while ylist containg the solution and derivative of the solution y dydx for each of the loops.

	gsl_vector* xlist = gsl_vector_alloc(max);
	gsl_matrix* ylist = gsl_matrix_alloc(max,n);
	double pi=atan(1.)*4;
	double a=0, b=8*pi, acc=0.01, eps=0.01;
	double h=0.1;
	gsl_vector_set(xlist,0,a);
	gsl_matrix_set(ylist,0,0,0);
	gsl_matric_set(ylist,0,1,1);
	int k = ode_driver(f,n,xlist,ylist,b,h,acc,eps,max);

	if(k<0)
		printf("max steps reached in ode_drives is \n");
	printf("Solution to exercise A: x=%g, y=%g and dydx=%g\n",gsl_vector_set(xlist,k-1),gsl_matrix_set(ylist,k-1,0),gsl_matric_get(ylist,k-1,1));
	printf("Theoretical solution to the problem is sin(%g)=%g\n",gsl_vector_get(xlist,k-1),sin(gsl_vector_get(xlist,k-1)));

	printf("The solution to exercise B is: \n");
	prinf("\n\n");

for(int i=0;i<k;i++)
	printf("%g %g\n",gsl_vector_get(xlist,i),gsl_matrix_get(ylist,i,0));
	printf("\n\n");
	printf("# m=1, S=0\n");

for(int i=0; i<k; i++)
	printf("%g %g\n",gsl_vector_get(xlist,i),sin(gsl_vector_get(xlist,i)));

	gsl_vector_free(xlist);
	gsl_matrix_free(ylist);

return 0;
}
