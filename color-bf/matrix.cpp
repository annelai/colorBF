#include "matrix.h"

// 3x3 matrix inverse
void Matrix33::inv()
{	
	double a22a33_a23a32 = a22*a33 - a23*a32;
	double a13a32_a12a33 = a13*a32 - a12*a33;
	double a12a23_a13a22 = a12*a23 - a13*a22;
	double a23a31_a21a33 = a23*a31 - a21*a33;
	double a11a33_a13a31 = a11*a33 - a13*a31;
	double a13a21_a11a23 = a13*a21 - a11*a23;
	double a21a32_a22a31 = a21*a32 - a22*a31;
	double a12a31_a11a32 = a12*a31 - a11*a32;
	double a11a22_a12a21 = a11*a22 - a12*a21;
	// inv(det(A))
	double detA = a11*a22a33_a23a32 + a12*a23a31_a21a33 + a13*a21a32_a22a31;
	detA = 1.f/detA;
	
	a11 = a22a33_a23a32*detA;
	a12 = a13a32_a12a33*detA;
	a13 = a12a23_a13a22*detA;
	a21 = a23a31_a21a33*detA;
	a22 = a11a33_a13a31*detA;
	a23 = a13a21_a11a23*detA;
	a31 = a21a32_a22a31*detA;
	a32 = a12a31_a11a32*detA;
	a33 = a11a22_a12a21*detA;
}

void CGSolver(double* x, const double* A, const double* b, const int n)
{
	const double thr_stopping = 0.00001;

	double* res = new double[n];
	double* p = new double[n];
	double* Ap = new double[n];
	double rs_old, rs_new, alpha, beta;
	
	for(int i = 0; i < n; ++i) {
		x[i] = b[i];
	}

	rs_old = 0;
	for(int i = 0; i < n; ++i) {
		// A*x;
		Ap[i] = 0;
		for(int j = 0; j < n; ++j) {
			Ap[i] += A[i*n+j]*x[j];
		}
		// r = b - A*x;
		res[i] = b[i] - Ap[i];
		// p = r;
		p[i] = res[i];
		// rs_old = r'*r;
		rs_old += res[i]*res[i];
	}

	for(int iter = 0; iter < n; ++iter) {
		alpha = 0;
		for(int i = 0; i < n; ++i) {
			// A*p
			Ap[i] = 0;
			for(int j = 0; j < n; ++j) {
				Ap[i] += A[i*n+j]*p[j];
			}
			// p'*A*p
			alpha += p[i]*Ap[i];
		}
		alpha = rs_old/alpha;
		rs_new = 0;
		for(int i = 0; i < n; ++i) {
			// x = x + alpha*p;
			x[i] += alpha*p[i];
			// r = r - alpha*Ap;
			res[i] -= alpha*Ap[i];
			// rs_new = r'*r;
			rs_new += res[i]*res[i];
		}
		if(rs_new < thr_stopping) {
			break;
		}
		beta = rs_new/rs_old;
		for(int i = 0; i < n; ++i) {
			p[i] = res[i] + beta*p[i];
		}
		rs_old = rs_new;
	}
}
