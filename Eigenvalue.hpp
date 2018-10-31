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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Eigenvalue
{
private:
	void power_iteration(int iterations);
	void QR_iteration(int iterations);
	void Preliminary_Reduction(string method);

public:
	int dimension;  // row_dimension = column_dimension
	Matrix * original_A;   // original_A record the original matrix
	Matrix * eigenvalues;  // eigenvalues matrix, n x 1.
	Matrix * eigenvectors;  // eigenvectors matrix, n x n, each column is an eigenvector.

	Matrix * reduced_Hessenberg; // H is the Hessenberg matrix preliminarily reduced from original_A
	Matrix * reduced_Q; // preliminary reduction: QT * A * Q = H 
	bool preliminary_reduction;  

	Eigenvalue(Matrix & matrix_A) {
		if (matrix_A.row_dimension != matrix_A.column_dimension) {
			cout << "row and column dimensions not matched!!!"; 
		}
		else {dimension = matrix_A.row_dimension; }

		original_A = new Matrix(matrix_A);

		eigenvalues = NULL;
		eigenvectors = NULL;
		reduced_Hessenberg = NULL;
		reduced_Q = NULL;

		preliminary_reduction = false;
	}

	void Get_eigenvalues(string method, bool use_preliminary_reduction, string reduction_method = "Arnoldi");

	void Print_eigenvalues() const {eigenvalues -> Print_matrix(); }
	void Print_eigenvectors() const {eigenvectors -> Print_matrix(); }

	~Eigenvalue() {
		delete original_A;
		delete eigenvalues;
		delete eigenvectors;
	}

};

void Eigenvalue::Get_eigenvalues(string method, bool use_preliminary_reduction, string reduction_method) {
	// reduction method: Arnoldi, Lanczos.

	if (use_preliminary_reduction == true) {Preliminary_Reduction(reduction_method); }
	else {
		preliminary_reduction = false;
	}

	if (method == "QR_iteration") {
		QR_iteration(1000);
	}
	else {
		cout << "Failed to get eigenvalues...";
		return;
	}

	cout << "Successfully get eigenvalues!\n";
}

void Eigenvalue::QR_iteration(int iterations) {
	Matrix Ak(dimension, dimension); 
	Matrix eigenvectors_matrix(dimension, dimension);

	if (preliminary_reduction == true) {
		eigenvectors_matrix = (* reduced_Q);
		Ak = (* reduced_Hessenberg);
	}
	else {
		eigenvectors_matrix.Generate_Identity();
		Ak = (* original_A);
	} 

	Matrix Q(Ak);
	Matrix R(Ak);

	Matrix A_record(Ak);
	int eigen_dimension = Ak.row_dimension;

	// QR iteration starts.
	for (int i = 0; i < iterations; ++i) {
		if (i % 5 == 4) {
			Ak.Matrix_approximation(1e-10); 
			if ((A_record - Ak).Frobenius_norm() < 1e-5 * eigen_dimension) {
				cout << "QR iterations:" << i << '\n';
				break;
			}
			else {A_record = Ak; }
		}

		Ak.QR_factorization();
		Q = Ak.QR_factorization_Q();
		R = Ak.QR_factorization_R();
		Ak = R * Q;
		eigenvectors_matrix = eigenvectors_matrix * Q;
	}

	
	eigenvalues = new Matrix(eigen_dimension, 1);
	eigenvectors = new Matrix(dimension, eigen_dimension);

	for (int i = 0; i < eigen_dimension; i++) {
		eigenvalues -> matrix[i][0] = Ak.matrix[i][i];
	}
	*eigenvectors = eigenvectors_matrix;

	eigenvalues -> Matrix_approximation();
	eigenvectors -> Matrix_approximation();
	// Ak.Print_matrix();
}

void Eigenvalue::Preliminary_Reduction(string method) {
	if (method == "Arnoldi") { // P62 chapter 04, Scientiﬁc Computing: An Introductory Survey.pdf
		preliminary_reduction = true;

		Matrix x0(dimension, 1); // x0 is an arbitrary nonzero starting vector.
		x0.Generate_Identity();
		// x0.Generate_Random();

		Matrix Q(dimension, dimension);
		Matrix H(dimension, dimension);

		Q.Set_column(0, x0 * (1.0 / x0.Frobenius_norm()) );
		int k = 0;
		for (k = 0; k < dimension; k++) {
			Matrix uk(dimension, 1);
			uk = (*original_A) * (Q.column_matrix(k));

			for (int j = 0; j <= k; j++) {
				H.matrix[j][k] = (Q.column_matrix(j).Transpose() * uk).matrix[0][0];
				uk = uk - Q.column_matrix(j) * H.matrix[j][k];
			}

			if (k < dimension - 1) {
				H.matrix[k + 1][k] = uk.Frobenius_norm(); 
			}
			else {break; }

			if (uk.Frobenius_norm() < 1e-8) {break; }

			Q.Set_column(k + 1, uk * (1.0 / H.matrix[k + 1][k]) );
		}

		reduced_Q = new Matrix(dimension, k + 1);
		reduced_Hessenberg = new Matrix(k + 1, k + 1);
		for (int i = 0; i < dimension; ++i) {
			for (int j = 0; j < k + 1; ++j) {
				reduced_Q -> matrix[i][j] = Q.matrix[i][j];
			}
		}
		for (int i = 0; i < k + 1; ++i) {
			for (int j = 0; j < k + 1; ++j) {
				reduced_Hessenberg -> matrix[i][j] = H.matrix[i][j];
			}
		}
	}
	else if (method == "Lanczos") {
		preliminary_reduction = true;


		Matrix x0(dimension, 1); // x0 is an arbitrary nonzero starting vector.
		x0.Generate_Identity();
		// x0.Generate_Random();

		Matrix Q(dimension, dimension);
		Matrix H(dimension, dimension);

		Q.Set_column(0, x0 * (1.0 / x0.Frobenius_norm()) );

		int k = 0;
		Matrix uk(dimension, 1);

		uk = (*original_A) * (Q.column_matrix(k));
		H.matrix[k][k] = (Q.column_matrix(k).Transpose() * uk).matrix[0][0];
		uk = uk - Q.column_matrix(k) * H.matrix[k][k];
		H.matrix[k + 1][k] = uk.Frobenius_norm();
		Q.Set_column(k + 1, uk * (1.0 / H.matrix[k + 1][k]) );

		for (k = 1; k < dimension; k++) {
			uk = (*original_A) * (Q.column_matrix(k));
			H.matrix[k][k] = (Q.column_matrix(k).Transpose() * uk).matrix[0][0];
			H.matrix[k - 1][k] = H.matrix[k][k - 1];
			uk = uk - Q.column_matrix(k) * H.matrix[k][k] - Q.column_matrix(k - 1) * H.matrix[k - 1][k];

			if (k < dimension - 1) {
				H.matrix[k + 1][k] = uk.Frobenius_norm(); 
			}
			else {break; }

			if (uk.Frobenius_norm() < 1e-8) {break; }

			Q.Set_column(k + 1, uk * (1.0 / H.matrix[k + 1][k]) );
		}

		reduced_Q = new Matrix(dimension, k + 1);
		reduced_Hessenberg = new Matrix(k + 1, k + 1);
		for (int i = 0; i < dimension; ++i) {
			for (int j = 0; j < k + 1; ++j) {
				reduced_Q -> matrix[i][j] = Q.matrix[i][j];
			}
		}
		for (int i = 0; i < k + 1; ++i) {
			for (int j = 0; j < k + 1; ++j) {
				reduced_Hessenberg -> matrix[i][j] = H.matrix[i][j];
			}
		}

	}
}

#endif