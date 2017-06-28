#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<math.h>

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

double vector_norm(gsl_vector* x){
	int n=x->size;
	double sum=0;
	for(int i=0;i<n;i++){
	sum+=pow(gsl_vector_get(x,i),2);
	}
return pow(sum,0.5);
}

int newton(void f(gsl_vector* x, gsl_vector* fx), gsl_vector* xstart, double dx, double epsilon){
	int n=xstart->size;
	gsl_matrix* J=gsl_matrix_alloc(n,n);
	gsl_vector* Dx=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_alloc(n);
	gsl_vector* fy=gsl_vector_alloc(n);
	gsl_vector* fx=gsl_vector_alloc(n);
	gsl_vector* df=gsl_vector_alloc(n);
	int steps=0;
	do{
		steps++;
		f(xstart,fx);
		for(int j=0; j<n;j++){
		gsl_vector_set(xstart,j,gsl_vector_get(xstart,j)+dx);
		f(xstart,df);
		gsl_vector_sub(df,fx);
		for(int i=0;i<n;i++){
		gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/dx);
		}
		gsl_vector_set(xstart,j,gsl_vector_get(xstart,j)-dx);
		}
		givens_qr(J);
		gsl_vector_scale(fx,-1.0);
		givens_qr_solve(J,fx,Dx);
		gsl_vector_scale(fx,-1.0);
		double lambda=2;

		do{
		lambda/=2;
		gsl_vector_memcpy(y,Dx);
		gsl_vector_scale(y,lambda);
		gsl_vector_add(y,xstart);
		f(y,fy);
		}while(vector_norm(fy)>(1-lambda/2)*vector_norm(fx) && lambda>0.02);
		gsl_vector_memcpy(xstart,y);		
		gsl_vector_memcpy(fx,fy);
		}while(vector_norm(Dx)>dx && vector_norm(fx)>epsilon);

	gsl_matrix_free(J);
	gsl_vector_free(Dx);
	gsl_vector_free(y);
	gsl_vector_free(fy);
	gsl_vector_free(fx);
	gsl_vector_free(df);
return steps;
}
