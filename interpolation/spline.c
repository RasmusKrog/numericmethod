#include "spline.h"
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<assert.h>
#include<math.h>
#define RND (double)rand()/RAND_MAX

double lspline(int n, gsl_vector*x,gsl_vector*y, double z){
	//fprintf(stderr, "n=%g z=%g xmin=%g xmax=%g\n",n,z,gsl_vector_get(x,0),gsl_vector_get(x,n-1));
	assert(n>1 && z>=gsl_vector_get(x,0) && z<=gsl_vector_get(x,n-1));
	int i=0, j=n-1;
	while (j-i>1){
		int m=(i+j)/2;
		if (z>gsl_vector_get(x,m)) i=m; else j=m;
	}
	double dx = gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
	double dy = gsl_vector_get(y,i+1)-gsl_vector_get(y,i);
	return gsl_vector_get(y,i)+dy/dx*(z-gsl_vector_get(x,i));
}

double lspline_integ(int n, gsl_vector*x,gsl_vector*y, double z){
	//fprintf(stderr, "n=%g z=%g xmin=%g xmax=%g\n",n,z,gsl_vector_get(x,0),gsl_vector_get(x,n-1));
	double x0 = gsl_vector_get(x,0);
	assert(n>1 && z>=x0 && z<=gsl_vector_get(x,n-1));
	int i=0, j=n-1;
	while (j-i>1){
		int m=(i+j)/2;
		if (z>gsl_vector_get(x,m)) i=m; else j=m;
	}
	double dx = gsl_vector_get(x,i+1)-gsl_vector_get(x,i);
	double dy = gsl_vector_get(y,i+1)-gsl_vector_get(y,i);
	double I0 = 0;
	for (int l=0; l<i; l++){
		double dxl = gsl_vector_get(x,l+1)-gsl_vector_get(x,l);
		double dyl = gsl_vector_get(y,l+1)-gsl_vector_get(y,l);
		I0+= gsl_vector_get(y,l)*dxl+0.5*dyl/dxl*dxl*dxl; 
	}
	return I0+gsl_vector_get(y,i)*(z-gsl_vector_get(x,i))+0.5*dy/dx*(z-gsl_vector_get(x,i))*(z-gsl_vector_get(x,i));
}

/* quadratic */

qspline *qspline_alloc(int n, gsl_vector *x, gsl_vector *y){
	qspline *s = (qspline*) malloc(sizeof(qspline));
	s->n = n;
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->b = (double*) malloc((n-1)*sizeof(double));
	s->c = (double*) malloc((n-1)*sizeof(double));

	for (int i = 0; i < n; ++i){
		s->x[i] = gsl_vector_get(x,i);
		s->y[i] = gsl_vector_get(y,i);
	}
	int i; double p[n-1], dx[n-1], dy[n-1];
	for (int i = 0; i < n-1; i++){
		dx[i] = gsl_vector_get(x,i+1) - gsl_vector_get(x,i);
		dy[i] = gsl_vector_get(y,i+1) - gsl_vector_get(y,i);
		p[i] = dy[i]/dx[i];
	 } 
	s->c[0]=0;
	//up:
	for(i=0; i<=n-3; i++){
		s->c[i+1] = (p[i+1] - p[i] - s->c[i]*dx[i]) /dx[i+1];
	}
	//down:
	for(i=n-2; i>=0; i--){
		s->c[i] = (p[i+1] - p[i] - s->c[i+1]*dx[i+1]) /dx[i];
	}
	for (int i = 0; i <= n-2; ++i){
		s->b[i] = p[i] - s->c[i] * dx[i];
	}
	return s;
}

double qspline_evaluate(qspline* s, double z){
	//fprintf(stderr, "z=%g, xmin=%g, xmax =%g\n",z,s->x[0],s->x[s->n-1]);
	assert(z >= s->x[0] && z<= s->x[s->n-1]);
	int i=0, j= s->n -1;
	while(j-i>1){
		int m=(i+j)/2;
	//	fprintf(stderr, "i=%i,j=%i,m=%i\n",i,j,m);
		if (z > s->x[m]) i=m;
		else j=m;
	}
	double h = z - s->x[i];
	//fprintf(stderr, "h = %g\n", h);
	return s->y[i] + s->b[i]*h + s->c[i] * h*h;
}

double qspline_deriv(qspline* s, double z){
	//fprintf(stderr, "z=%g, xmin=%g, xmax =%g\n",z,s->x[0],s->x[s->n-1]);
	assert(z >= s->x[0] && z<= s->x[s->n-1]);
	int i=0, j= s->n -1;
	while(j-i>1){
		int m=(i+j)/2;
	//	fprintf(stderr, "i=%i,j=%i,m=%i\n",i,j,m);
		if (z > s->x[m]) i=m;
		else j=m;
	}
	double h = z - s->x[i];
	//fprintf(stderr, "h = %g\n", h);
	return s->b[i] + 2*s->c[i] * h;
}

double qspline_integ(qspline* s, double z){
	//fprintf(stderr, "z=%g, xmin=%g, xmax =%g\n",z,s->x[0],s->x[s->n-1]);
	assert(z >= s->x[0] && z<= s->x[s->n-1]);
	int i=0, j= s->n -1;
	while(j-i>1){
		int m=(i+j)/2;
	//	fprintf(stderr, "i=%i,j=%i,m=%i\n",i,j,m);
		if (z > s->x[m]) i=m;
		else j=m;
	}
	double h = z - s->x[i];
	//fprintf(stderr, "h = %g\n", h);
	double I0 = 0;
	for (int l=0; l<i; l++){
		double hl = s->x[l+1] - s->x[l];
		I0+= s->y[l]*hl + 0.5*s->b[l]*hl*hl + 1/3*s->c[l]*hl*hl*hl; 
	}
	return I0 + s->y[i]*h + 0.5*s->b[i]*h*h + 1/3*s->c[i]*h*h*h; 
}

void qspline_free(qspline* s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s);
}

/* Cubic spline */

cspline *cspline_alloc(int n, gsl_vector *x, gsl_vector *y){
	cspline *s = (cspline*) malloc(sizeof(cspline));
	s->n = n;
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->b = (double*) malloc(n*sizeof(double));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->d = (double*) malloc((n-1)*sizeof(double));

	for (int i = 0; i < n; ++i){
		s->x[i] = gsl_vector_get(x,i);
		s->y[i] = gsl_vector_get(y,i);
	}
	int i; double p[n-1], h[n-1], dy[n-1];
	for (int i = 0; i < n-1; i++){
		h[i] = gsl_vector_get(x,i+1) - gsl_vector_get(x,i);
		assert(h[i]>0);
		dy[i] = gsl_vector_get(y,i+1) - gsl_vector_get(y,i);
		p[i] = dy[i]/h[i];
	}
	//tridiagonal system:
	double D[n], Q[n-1], B[n];
	D[0]=2; D[n-2]=2; Q[0]=1;
	B[0]=3*p[0]; B[n-1]=3*p[n-2];
	for (int i = 0; i < n-2; ++i){
		D[i+1]=2*h[i]/h[i+1]+2;
	 	Q[i+1]=h[i]/h[i+1];
	 	B[i+1]=3*(p[i]+p[i+1]*h[i]/h[i+1]);
	}
	//gauss elimination:
	for(int i=1;i<n;i++){
		D[i]-= Q[i-1]/D[i-1];
		B[i]-= B[i-1]/D[i-1];
	}
	s->b[n-1] = B[n-1]/D[n-1];
	//back substitution:
	for (i=n-2;i>=0;i--){
		s->b[i] = (B[i]-Q[i]*s->b[i+1])/D[i];
	}
	for(int i=0;i<n-1;i++){
		s->c[i] = (-2*s->b[i] - s->b[i+1] +3*p[i])/h[i];
		s->d[i] = (s->b[i] + s->b[i+1] -2*p[i])/h[i]/h[i];
	}
	return s;
}

double cspline_evaluate(cspline* s, double z){
	//fprintf(stderr, "z=%g, xmin=%g, xmax =%g\n",z,s->x[0],s->x[s->n-1]);
	assert(z >= s->x[0] && z<= s->x[s->n-1]);
	int i=0, j= s->n -1;
	while(j-i>1){
		int m=(i+j)/2;
	//	fprintf(stderr, "i=%i,j=%i,m=%i\n",i,j,m);
		if (z > s->x[m]) i=m;
		else j=m;
	}
	double h = z - s->x[i];
	//fprintf(stderr, "h = %g\n", h);
	return s->y[i] + s->b[i]*h + s->c[i] * h*h + s->d[i]*h*h*h;
}

double cspline_deriv(cspline* s, double z){
	//fprintf(stderr, "z=%g, xmin=%g, xmax =%g\n",z,s->x[0],s->x[s->n-1]);
	assert(z >= s->x[0] && z<= s->x[s->n-1]);
	int i=0, j= s->n -1;
	while(j-i>1){
		int m=(i+j)/2;
	//	fprintf(stderr, "i=%i,j=%i,m=%i\n",i,j,m);
		if (z > s->x[m]) i=m;
		else j=m;
	}
	double h = z - s->x[i];
	//fprintf(stderr, "h = %g\n", h);
	return s->b[i] + 2*s->c[i] * h + 3*s->d[i]*h*h;
}

double cspline_integ(cspline* s, double z){
	//fprintf(stderr, "z=%g, xmin=%g, xmax =%g\n",z,s->x[0],s->x[s->n-1]);
	assert(z >= s->x[0] && z<= s->x[s->n-1]);
	int i=0, j= s->n -1;
	while(j-i>1){
		int m=(i+j)/2;
	//	fprintf(stderr, "i=%i,j=%i,m=%i\n",i,j,m);
		if (z > s->x[m]) i=m;
		else j=m;
	}
	double h = z - s->x[i];
	//fprintf(stderr, "h = %g\n", h);
	double I0 = 0;
	for (int l=0; l<i; l++){
		double hl = s->x[l+1] - s->x[l];
		I0+= s->y[l]*hl + 0.5*s->b[l]*hl*hl + 1/3*s->c[l]*hl*hl*hl + 1/4*s->d[l]*hl*hl*hl*hl; 
	}
	return I0 + s->y[i]*h + 0.5*s->b[i]*h*h + 1/3*s->c[i]*h*h*h + 1/4*s->d[i]*h*h*h*h; 
}

void cspline_free(cspline* s){
	free(s->x);
	fprintf(stderr, "x\n" );
	free(s->y);
	fprintf(stderr, "y\n" );
	free(s->b);
	fprintf(stderr, "b\n" );
	free(s->c);
	fprintf(stderr, "c\n" );
	free(s->d);
	fprintf(stderr, "d\n" );
	free(s);
	fprintf(stderr, "s\n" );
}



