#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x);
void Ravn(double *x, int n, int l, int r);
void Cheb(double *x, int n, int l, int r);
double P(double x, double *a, int n);
double Inter(double *x, double y, int n);
void out(double *x, double *a, FILE *fo, int n);
int SolveSystem(int n, double* a, double* b, double* x);

double f(double x){
	return exp(x);
}

void Ravn(double *x, int n, int l, int r){
    int i;
	for(i=0; i<n; i++)
		x[i] = l + (double)i*(r - l)/(n - 1);
}

void Cheb(double *x, int n, int l, int r){
    int i;
	for(i = 1; i<n+1; i++)
		x[i-1] = (cos(M_PI*(2*i-1)/2/n)*(r - l) + (r + l))/2;
}

void out(double *x, double *a, FILE *output, int n){
    int i;
	double *xnew = malloc((3*n - 2)*sizeof(double));
	double *h = malloc((n - 1)*sizeof(double));

	for (i=0; i<n-1; i++){
		h[i] = (x[i+1] - x[i])/3;
	}
	for (i=0; i<3*n-3; i++){
		xnew[i] = x[i/3] + h[i/3]*(i%3);
	}
	xnew[3*n-3] = x[n-1];
	for(i=0; i<3*n-2; i++){
		fprintf(output, "x= %19.17lf  f(x)= %19.17lf  P1= %19.17lf P2= %19.17lf  diff1= %19.17lf  diff2= %19.17lf \r\n", xnew[i], f(xnew[i]), P(xnew[i],a,n), Inter(x,xnew[i],n), f(xnew[i])-P(xnew[i],a,n),f(xnew[i])-Inter(x,xnew[i],n));
	}
}

double P(double x,double *a, int n){
    int i;
	double p=0;
	for(i=0; i<n; i++) p+= a[i]*pow(x,i);
	return p;
}

double Inter(double *x, double y, int n){
    int i,j;
	double P=0;
	double *Phi = malloc(n*sizeof(double));
	for (i=0; i<n; i++) Phi[i]=1;
	for (i=0; i<n; i++){
		for (j=0; j<i; j++)
			Phi[i]*=((y - x[j])/(x[i] - x[j]));
		for (j=i+1; j<n; j++)
			Phi[i]*=((y - x[j])/(x[i] - x[j]));
	}
	for(i=0; i<n; i++){
		P+= f(x[i])*Phi[i];
	}
	return P;
}

int SolveSystem(int n, double* a, double* b, double* x){
	int i,j,k;
	double r;
	double tmp1, tmp2;
	double cosPhi, sinPhi;

	for (i=0; i<n; i++){
		for (j=i+1; j<n; j++){
			tmp1 = a[i*n + i];
			tmp2 = a[j*n + i];

			r = sqrt(tmp1*tmp1 + tmp2*tmp2);

			if (r<1e-100)
				return -1;

			cosPhi =  tmp1 / r;
			sinPhi = -tmp2 / r;

			a[i * n + i] = r;
			a[j * n + i] = 0.0;

			for (k=i+1; k<n; k++){
				tmp1 = a[i*n + k];
				tmp2 = a[j*n + k];

				a[i*n + k] = tmp1*cosPhi - tmp2*sinPhi;
				a[j*n + k] = tmp1*sinPhi + tmp2*cosPhi;
			}

			tmp1 = b[i];
			tmp2 = b[j];

			b[i] = tmp1*cosPhi - tmp2*sinPhi;
			b[j] = tmp1*sinPhi + tmp2*cosPhi;
		}
	}

	for (i=n-1; i>=0; i--){
		tmp1 = b[i];
		for (j=i+1; j<n; j++)
			tmp1-= a[i*n + j]*x[j];
		x[i] = tmp1/a[i*n + i];
	}

	return 0;
}

int main (void){
    int n,b,l,r;
    double *x;
    double *y;
    double *A;
    double *a;
    FILE *output;

printf("Enter quantity of points: ");
scanf("%d", &n);

x = (double*)malloc(n*sizeof(double));
y = (double*)malloc(n*sizeof(double));
A = (double*)malloc(n*n*sizeof(double));
a = (double*)malloc(n*sizeof(double));

printf("Choose <0> - Cheb or <1> - Ravn: ");
scanf("%d", &b);
printf("Interval:\nFrom ");
scanf("%d", &l);
printf("To ");
scanf("%d", &r);


if (b==0) {Cheb(x,n,l,r);}
else {Ravn(x,n,l,r);}


output = fopen("output.txt", "w");
for(int i=0; i<n; i++) {
	y[i] = f(x[i]);
}
for (int i=0; i<n; i++)
	for (int j=0; j<n; j++)
		A[i*n+j] = pow(x[i],j);
SolveSystem(n,A,y,a);

out(x,a,output,n);

fclose(output);
free(x);
return 0;
}
