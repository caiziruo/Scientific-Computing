#ifndef Eigenvalue_hpp
#define Eigenvalue_hpp

#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>
#include <string>

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
		eigenvalues = new Matrix(dimension, 1);
		eigenvectors = new Matrix(dimension, dimension);

		preliminary_reduction = false;
	}

	void Get_eigenvalues(string method);
	void power_iteration(int iterations);
	void QR_iteration(int iterations);
	void Preliminary_Reduction();

	void Print_eigenvalues() const;
	void Print_eigenvectors() const;

	~Eigenvalue() {
		delete original_A;
		delete eigenvalues;
		delete eigenvectors;
	}

};

void Eigenvalue::Get_eigenvalues(string method) {
	if (method == "QR_iteration") {
		QR_iteration(1000);
	}
}

void Eigenvalue::Print_eigenvalues() const {
	for (int i = 0; i < dimension; ++i) {
		cout << eigenvalues -> matrix[i][0] << '\n';
	}
} 

void Eigenvalue::Print_eigenvectors() const {
	for (int i = 0; i < dimension; ++i) {
		for (int j = 0; j < dimension; ++j) {
			cout << setw(14) <<  eigenvectors -> matrix[i][j];
		}
		cout << '\n';
	}
}

void Eigenvalue::QR_iteration(int iterations) {
	Matrix Ak(dimension, dimension); 
	Matrix eigenvectors_matrix(dimension, dimension);
	eigenvectors_matrix.Generate_Identity();

	if (preliminary_reduction == true) {

	}
	else {
		Ak = (* original_A); 
	}

	Matrix Q(dimension, dimension);
	Matrix R(dimension, dimension);

	for (int i = 0; i < iterations; ++i) {
		if (i % 10 == 0) {Ak.Matrix_approximation(1e-10); }

		Ak.QR_factorization();
		Q = Ak.QR_factorization_Q();
		R = Ak.QR_factorization_R();
		Ak = R * Q;
		eigenvectors_matrix = eigenvectors_matrix * Q; 
	}

	for (int i = 0; i < dimension; i++) {
		eigenvalues -> matrix[i][0] = Ak.matrix[i][i];
	}
	*eigenvectors = eigenvectors_matrix;
}



#endif