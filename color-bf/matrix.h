#ifndef MATRIX_H
#define MATRIX_H

template <class T>
class vec3 {
public:
	vec3() {}
	~vec3() {}
	T x;
	T y;
	T z;
};

// A 3x3 Matrix
class Matrix33 {
public:
	Matrix33() {}
	~Matrix33() {}
	// inverse
	void inv();
// defined as public data members for fast access
	double a11;
	double a12;
	double a13;
	double a21;
	double a22;
	double a23;
	double a31;
	double a32;
	double a33;
};

/*
 Conjugate gradient
 http://en.wikipedia.org/wiki/Conjugate_gradient_method
 Solve Ax = b
 A: n*n elements array, representing a n-by-n PSD matrix
 b: n elements array, representing a n-by-1 vector
 x: n elements array, representing a n-by-1 vector
*/
void CGSolver(double* x, const double* A, const double* b, const int n);

#endif