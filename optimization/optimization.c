#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>
#include<stdlib.h>
#include "optimization.h"

double vector_norm(gsl_vector* x){
	int n=x->size;
	double sum=0;
	for(int i=0;i<n;i++){
	sum+=pow(gsl_vector_get(x,i),2);
	}
return pow(sum,0.5);
}

void matrix_vector_prod(gsl_matrix* A, gsl_vector* b, gsl_vector* out){
	gsl_vector* c = gsl_vector_alloc(A->size1);
	for (int i = 0; i < A->size1; i++){
	double P=0;
	for(int j=0;j < b->size; j++){
	P+=gsl_matrix_get(A,i,j)*gsl_vector_get(b,j);
	}
	gsl_vector_set(c,i,P);
	}
	for (int i = 0; i < c->size; ++i){
	gsl_vector_set(out,i,gsl_vector_get(c,i));
	}
}

double rosenbrock(gsl_vector* x){
	double x0 = gsl_vector_get(x,0);
	double x1 = gsl_vector_get(x,1);
	double f = pow((1-x0),2) + 100*pow((x1-x0*x0),2);
return f;
}

double himmelblau(gsl_vector* x){
	double x0 = gsl_vector_get(x,0);
	double x1 = gsl_vector_get(x,1);
	double f = pow((x1+x0*x0-11),2) + pow((x0+x1*x1-7),2);
return f;
}

void dfrosenbrock(gsl_vector* x, gsl_vector* df){
	double x0 = gsl_vector_get(x,0);
	double x1 = gsl_vector_get(x,1);
	double f0 = 2*(x0-1) - 400*(x1-pow(x0,2))*x0;
	double f1 = 200*(x1 - pow(x0,2));
	gsl_vector_set(df,0,f0);
	gsl_vector_set(df,1,f1);	
}

void dfhimmelblau(gsl_vector* x, gsl_vector* df){
	double x0 = gsl_vector_get(x,0);
	double x1 = gsl_vector_get(x,1);
	double f0 = 4*x0*(pow(x0,2)+x1-11) + 2*(x0+pow(x1,2)-7);
	double f1 = 2*(pow(x0,2)+x1-11) + 4*x1*(x0+pow(x1,2)-7);
	gsl_vector_set(df,0,f0);
	gsl_vector_set(df,1,f1);	
}

void Hrosenbrock(gsl_vector* x, gsl_matrix* H){
	double x0=gsl_vector_get(x,0);
	double x1=gsl_vector_get(x,1);
	double H00, H11, H01, H10;
	H00 = 2-400*(x1-3*pow(x0,2));
	H01 = -400*x0;
	H10 = -400*x0;
	H11 = 200;
	gsl_matrix_set(H,0,0,H00);
	gsl_matrix_set(H,0,1,H01);
	gsl_matrix_set(H,1,0,H10);
	gsl_matrix_set(H,1,1,H11);
}

void Hhimmelblau(gsl_vector* x, gsl_matrix* H){
	double x0=gsl_vector_get(x,0);
	double x1=gsl_vector_get(x,1);
	double H00, H11, H01, H10;
	H00 = 12*pow(x0,2) + 4*x1 - 42;
	H01 = 4*(x0 + x1);
	H10 = 4*(x0 + x1);
	H11 = 12*pow(x1,2) + 4*x0 - 26;
	gsl_matrix_set(H,0,0,H00);
	gsl_matrix_set(H,0,1,H01);
	gsl_matrix_set(H,1,0,H10);
	gsl_matrix_set(H,1,1,H11);
}

void givens_qr(gsl_matrix* A){
	for(int q=0;q<A->size2;q++){
	for(int p=q+1;p<A->size1;p++){
	double theta=atan2(gsl_matrix_get(A,p,q),gsl_matrix_get(A,q,q));
	for(int k=q;k<A->size2;k++){
	double xq=gsl_matrix_get(A,q,k);
	double xp=gsl_matrix_get(A,p,k);
	gsl_matrix_set(A,q,k,xq*cos(theta)+xp*sin(theta));
	gsl_matrix_set(A,p,k,-xq*sin(theta)+xp*cos(theta));
	}
	gsl_matrix_set(A,p,q,theta); 
	}
	}
}

void givens_qr_QTvec(gsl_matrix* QR, gsl_vector* v){
	for(int q=0; q<QR->size2; q++){
	for(int p=q+1;p<QR->size1;p++){
	double theta=gsl_matrix_get(QR,p,q);
	double vq=gsl_vector_get(v,q);
	double vp=gsl_vector_get(v,p);
	gsl_vector_set(v,q,vq*cos(theta)+vp*sin(theta));
	gsl_vector_set(v,p,-vq*sin(theta)+vp*cos(theta));
	}
	}
}

void givens_qr_solve(gsl_matrix* QR, gsl_vector* b, gsl_vector* x){
	givens_qr_QTvec(QR,b); 
	for (int i=QR->size2-1; i>=0; i--){ //back-substitution
	double s=0;
	for(int k=i+1; k<QR->size2; k++) s+=gsl_matrix_get(QR,i,k)*gsl_vector_get(x,k);
	gsl_vector_set(x,i,(gsl_vector_get(b,i)-s)/gsl_matrix_get(QR,i,i));
	}
}

int newtonMinimization(double f(gsl_vector* x), void gradient(gsl_vector *x, gsl_vector* df), void hessian(gsl_vector* x,gsl_matrix* H), gsl_vector* xstart, double eps){
	int n=xstart->size;
	gsl_vector* Dx  = gsl_vector_alloc(n);
	gsl_vector* y   = gsl_vector_alloc(n);
	gsl_matrix* H   = gsl_matrix_alloc(n,n);
	gsl_vector* HDx = gsl_vector_alloc(n);
	double fx, Dxd;
	gsl_vector* dfxxx = gsl_vector_alloc(n);
	gsl_vector* dfx   = gsl_vector_alloc(n);
	gsl_vector* df    = gsl_vector_alloc(n);
	double dx = 0.001;
	int steps = 0;

	do{
	steps++;
	fx = f(xstart);
	gradient(xstart,dfx);
	hessian(xstart,H);
	givens_qr(H);
	gsl_vector_scale(dfx,-1.0);
	givens_qr_solve(H,dfx,Dx);
	gsl_vector_scale(dfx,-1.0);
	matrix_vector_prod(H,Dx,HDx);
	for (int i = 0; i < n; i++){
	gsl_vector_set(dfxxx,i,gsl_vector_get(dfx,i)+gsl_vector_get(HDx,i));
	}
	double lambda=2;
	do{
	lambda/=2;
	gsl_vector_memcpy(y,Dx);
	gsl_vector_scale(y,lambda);
	gsl_vector_add(y,xstart);
	Dxd = 0;
	for (int i = 0; i < n; ++i){
	Dxd+=gsl_vector_get(Dx,i)*gsl_vector_get(dfx,i);
	}
	}while(fabs(f(y)) > fabs(fx) + 0.0001*lambda*(Dxd)  && lambda>0.02);
	gsl_vector_memcpy(xstart,y);	
	}while(vector_norm(Dx)>dx && vector_norm(dfx)>eps);

	gsl_vector_free(Dx);
	gsl_vector_free(y);
	gsl_vector_free(df);
	gsl_matrix_free(H);
	gsl_vector_free(HDx);
	gsl_vector_free(dfxxx);
	return steps;
}

