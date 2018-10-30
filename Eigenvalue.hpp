#ifndef Eigenvalue_hpp
#define Eigenvalue_hpp

#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "Matrix.hpp"

class Eigenvalue
{
private:

public:
	int dimension;
	Matrix * original_A; 
	Matrix * eigenvalues;
	Matrix * eigenvectors;
	bool preliminary_reduction;

	Eigenvalue(Matrix & matrix_A, int m, int n) {
		if (m != n) {cout << "row and column dimensions not matched!!!"; }
		else {dimension = m; }

		original_A = new Matrix(matrix_A);
		eigenvalues = NULL;
		eigenvectors = new Matrix(dimension, dimension);
		
		preliminary_reduction = false;
    }

    void power_iteration();
    void QR_iteration();
    void Preliminary_Reduction();
    ~Eigenvalue() {
        delete original_A;
        delete eigenvalues;
        delete eigenvectors;
    }

};

void Eigenvalue::QR_iteration() {
	if (preliminary_reduction == true) {

	}
	else {

	}

}



#endif