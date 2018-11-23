#ifndef Interpolation_hpp
#define Interpolation_hpp

#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "Matrix.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Runge_function(double x) {
	return 1.0 / (1.0 + 25.0 * x * x);
}

class Interpolation
{
private:
	double Newton_basis_function(double t, int j) { // j = 1, 2, ... n, n = dimension_data
		double x = 1.0;
		for (int i = 0; i < j - 1; i++) {
			x *= (t - data_t -> matrix[i][0]); 
		}
		return x;
	}
	double Newton_basis_interpolation_polynomial(double t, Matrix parameters) {
		double output = 0;

		for (int i = 0; i < dimension_data; ++i) {
			output += parameters.matrix[i][0] * Newton_basis_function(t, i + 1);
		}
		return output;
	}
	// double B_spline_function(double t, int i, int k) {  // t_{i} = 2.0 * (i - 1) / (dimension_data - 1.0) - 1.0 
		// double ti = 2.0 * (i - 1.0) / (dimension_data - 1.0) - 1.0;
		// double h = 2.0 / (dimension_data - 1.0);

		// if (k == 0) {
		// 	if ( (t >= ti) && (t < ti + h) ) {return 1.0; }
		// 	else {return 0.0;}
		// }
		// else {
		// 	return (t - ti) / (k * h) * B_spline_function(t, i, k - 1) 
		// 	+ (ti + (k + 1) * h - t) / (k * h) * B_spline_function(t, i + 1, k - 1);
		// }
	// }
	double B_spline_function(double t, int i, int k) {
		double ti = 2.0 * (i - 1.0) / (dimension_data - 1.0) - 1.0;
		double h = 2.0 / (dimension_data - 1.0);

		if (k == 1) {
			if ((t > ti) && (t < ti+h)) {
				return (t - ti) / h;
			}
			else if ((t > ti + h) && (t < ti + 2.0 * h)) {
				return (ti + 2 * h - t) / h;
			}
			else if (t == ti + h) {return 1;}
			else {return 0;}
		}
		else {
			return (t - ti) / (k * h) * B_spline_function(t, i, k - 1) 
			+ (ti + (k + 1) * h - t) / (k * h) * B_spline_function(t, i + 1, k - 1);
		}
	}

public:
	int dimension_data;
	Matrix * data_t;
	Matrix * data_y;
	Matrix * Interpolation_t;
	Matrix * Interpolation_y;

	Interpolation(int n) { // n > 2
		dimension_data = n;

		data_t = new Matrix(dimension_data, 1);
		data_y = new Matrix(dimension_data, 1);

		// initialize data points t_{i} and y_{i}
		for (int i = 0; i < dimension_data; i++) {
			data_t -> matrix[i][0] = 2.0 * i / (dimension_data - 1.0) - 1.0;
			data_y -> matrix[i][0] = Runge_function(data_t -> matrix[i][0]);
		}

	}

	// void Check_B_spline() {
	// 	for (int i = 0; i < 101; i++) {
	// 		double tj = 2.0 * i / (101.0 - 1.0) - 1.0;
	// 		// cout << 2.0 * i / (101.0 - 1.0) - 1.0 << "\n";
	// 		// cout << B_spline_function(2.0 * i / (101.0 - 1.0), 2, 1) << " ";
	// 		// for (int j = 0; j < dimension_data + 2; j++) {
	// 			// cout << B_spline_function(data_t -> matrix[i][0], j - 2, 3) << " "; 
	// 		cout << tj << " " << B_spline_function(tj, 1, 3) << " ";
	// 		// }
	// 		cout << '\n';
	// 	}
	// }
	void Newton_Interpolation(int m); // Print m points on [-1, 1] of the interpolation function.
	void Shooting_method(int m);
	void B_spline(int m);
	void Print_Interpolation_t() {Interpolation_t -> Print_matrix();}
	void Print_Interpolation_y() {Interpolation_y -> Print_matrix();}

	~Interpolation() {
		delete data_t;
		delete data_y;
	}
	
};

void Interpolation::Newton_Interpolation(int m) {
	Matrix parameters_linear_system(dimension_data, dimension_data);
	Matrix parameters(dimension_data, 1);

	for (int i = 0; i < dimension_data; ++i) {
		for (int j = 0; j < dimension_data; ++j) {
			parameters_linear_system.matrix[i][j] = Newton_basis_function(data_t -> matrix[i][0], j + 1);
		}
	}
	
	parameters = parameters_linear_system.Inverse() * (*data_y);

	Interpolation_t = new Matrix(m, 1);
	Interpolation_y = new Matrix(m, 1);
	for (int i = 0; i < m; ++i) {
		double t_i = 2.0 * i / (m - 1.0) - 1.0;
		Interpolation_t -> matrix[i][0] = t_i;
		Interpolation_y -> matrix[i][0] = Newton_basis_interpolation_polynomial(t_i, parameters);
	}

}

void Interpolation::Shooting_method(int m) {
	Matrix * parameters_linear_system [dimension_data - 1];
	Matrix * parameters [dimension_data - 1];
	Matrix * b [dimension_data - 1];

	for (int i = 0; i < dimension_data - 1; ++i) { // initialize 4(n-1) parameters and (n-1) equations.
		parameters_linear_system[i] = new Matrix(4, 4);
		parameters[i] = new Matrix(4, 1);
		b[i] = new Matrix(4, 1);

		parameters_linear_system[i] -> matrix[0][0] = 1.0;
		parameters_linear_system[i] -> matrix[0][1] = data_t -> matrix[i][0];
		parameters_linear_system[i] -> matrix[0][2] = pow(data_t -> matrix[i][0], 2);
		parameters_linear_system[i] -> matrix[0][3] = pow(data_t -> matrix[i][0], 3);
		parameters_linear_system[i] -> matrix[1][0] = 1.0;
		parameters_linear_system[i] -> matrix[1][1] = data_t -> matrix[i+1][0];
		parameters_linear_system[i] -> matrix[1][2] = pow(data_t -> matrix[i+1][0], 2);
		parameters_linear_system[i] -> matrix[1][3] = pow(data_t -> matrix[i+1][0], 3);
		parameters_linear_system[i] -> matrix[2][0] = 0.0;
		parameters_linear_system[i] -> matrix[2][1] = 1.0;
		parameters_linear_system[i] -> matrix[2][2] = 2.0 * data_t -> matrix[i][0];
		parameters_linear_system[i] -> matrix[2][3] = 3.0 * pow(data_t -> matrix[i][0], 2);
		parameters_linear_system[i] -> matrix[3][0] = 0.0;
		parameters_linear_system[i] -> matrix[3][1] = 0.0;
		parameters_linear_system[i] -> matrix[3][2] = 2.0;
		parameters_linear_system[i] -> matrix[3][3] = 6.0 * data_t -> matrix[i][0];
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// initialize f'(t_{1})
	double f_first_derivative_t_1 = 1.0;
	double f_second_derivative_t_n_1;
	double f_first_derivative_t_i = f_first_derivative_t_1;
	double f_second_derivative_t_i = 0.0;

	// compute f_second_derivative_t_n_1
	for (int i = 0; i < dimension_data - 1; ++i) {
		b[i] -> matrix[0][0] = data_y -> matrix[i][0];
		b[i] -> matrix[1][0] = data_y -> matrix[i+1][0];
		b[i] -> matrix[2][0] = f_first_derivative_t_i;
		b[i] -> matrix[3][0] = f_second_derivative_t_i;

		*parameters[i] = parameters_linear_system[i] -> Inverse() * (*b[i]);

		f_first_derivative_t_i = 
		1.0 * parameters[i] -> matrix[1][0] + 
		2.0 * parameters[i] -> matrix[2][0] * data_t -> matrix[i+1][0] + 
		3.0 * parameters[i] -> matrix[3][0] * pow(data_t -> matrix[i+1][0], 2);
		f_second_derivative_t_i = 2.0 * parameters[i] -> matrix[2][0] + 6.0 * parameters[i] -> matrix[3][0] * data_t -> matrix[i+1][0];
	}
	f_second_derivative_t_n_1 = f_second_derivative_t_i;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	f_first_derivative_t_1 = 0.0;
	double f_second_derivative_t_n_2;
	f_first_derivative_t_i = f_first_derivative_t_1;
	f_second_derivative_t_i = 0.0;
	// compute f_second_derivative_t_n_2
	for (int i = 0; i < dimension_data - 1; ++i) {
		b[i] -> matrix[0][0] = data_y -> matrix[i][0];
		b[i] -> matrix[1][0] = data_y -> matrix[i+1][0];
		b[i] -> matrix[2][0] = f_first_derivative_t_i;
		b[i] -> matrix[3][0] = f_second_derivative_t_i;

		*parameters[i] = parameters_linear_system[i] -> Inverse() * (*b[i]);

		f_first_derivative_t_i = 
		1.0 * parameters[i] -> matrix[1][0] + 
		2.0 * parameters[i] -> matrix[2][0] * data_t -> matrix[i+1][0] + 
		3.0 * parameters[i] -> matrix[3][0] * pow(data_t -> matrix[i+1][0], 2);
		f_second_derivative_t_i = 2.0 * parameters[i] -> matrix[2][0] + 6.0 * parameters[i] -> matrix[3][0] * data_t -> matrix[i+1][0];
	}
	f_second_derivative_t_n_2 = f_second_derivative_t_i;
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// compute true f'(t1)
	f_first_derivative_t_i = - f_second_derivative_t_n_2 / (f_second_derivative_t_n_1 - f_second_derivative_t_n_2);
	f_second_derivative_t_i = 0.0;
	// compute true parameters
	for (int i = 0; i < dimension_data - 1; ++i) {
		b[i] -> matrix[0][0] = data_y -> matrix[i][0];
		b[i] -> matrix[1][0] = data_y -> matrix[i+1][0];
		b[i] -> matrix[2][0] = f_first_derivative_t_i;
		b[i] -> matrix[3][0] = f_second_derivative_t_i;

		*parameters[i] = parameters_linear_system[i] -> Inverse() * (*b[i]);

		f_first_derivative_t_i = 
		1.0 * parameters[i] -> matrix[1][0] + 
		2.0 * parameters[i] -> matrix[2][0] * data_t -> matrix[i+1][0] + 
		3.0 * parameters[i] -> matrix[3][0] * pow(data_t -> matrix[i+1][0], 2);
		f_second_derivative_t_i = 2.0 * parameters[i] -> matrix[2][0] + 6.0 * parameters[i] -> matrix[3][0] * data_t -> matrix[i+1][0];
	}

	// compute 
	Interpolation_t = new Matrix(m, 1);
	Interpolation_y = new Matrix(m, 1);
	int i = 0;
	for (int j = 0; j < m; j++) {
		double t_j = 2.0 * j / (m - 1.0) - 1.0;
		Interpolation_t -> matrix[j][0] = t_j;

		if (t_j > data_t -> matrix[i+1][0]) i++;
		Interpolation_y -> matrix[j][0] = 
		parameters[i] -> matrix[0][0] * 1.0 + 
		parameters[i] -> matrix[1][0] * t_j + 
		parameters[i] -> matrix[2][0] * pow(t_j, 2) + 
		parameters[i] -> matrix[3][0] * pow(t_j, 3);
	}

}

void Interpolation::B_spline(int m) {
	Matrix parameters_linear_system(dimension_data + 2, dimension_data + 2);
	Matrix parameters(dimension_data + 2, 1);
	Matrix b(dimension_data + 2, 1);

	for (int i = 0; i < dimension_data; i++) {
		for (int j = 0; j < dimension_data + 2; j++) {
			parameters_linear_system.matrix[i][j] = B_spline_function(data_t -> matrix[i][0], j - 2, 3);
		}
		b.matrix[i][0] = data_y -> matrix[i][0];
	}

	b.matrix[dimension_data][0] = 0;
	b.matrix[dimension_data + 1][0] = 0;
	double h = 0.000001;
	// finite difference to compute second derivative
	parameters_linear_system.matrix[dimension_data][0] = (
		B_spline_function(data_t -> matrix[0][0] - h, -2, 3) - 
		B_spline_function(data_t -> matrix[0][0], -2, 3) * 2 + 
		B_spline_function(data_t -> matrix[0][0] + h, -2, 3)
		) / pow(h, 2);
	parameters_linear_system.matrix[dimension_data][1] = (
		B_spline_function(data_t -> matrix[0][0] - h, -1, 3) - 
		B_spline_function(data_t -> matrix[0][0], -1, 3) * 2 + 
		B_spline_function(data_t -> matrix[0][0] + h, -1, 3)
		) / pow(h, 2);
	parameters_linear_system.matrix[dimension_data][2] = (
		B_spline_function(data_t -> matrix[0][0] - h, 0, 3) - 
		B_spline_function(data_t -> matrix[0][0], 0, 3) * 2 + 
		B_spline_function(data_t -> matrix[0][0] + h, 0, 3)
		) / pow(h, 2);
	parameters_linear_system.matrix[dimension_data + 1][dimension_data - 1] = (
		B_spline_function(data_t -> matrix[dimension_data - 1][0] - h, dimension_data - 3, 3) - 
		B_spline_function(data_t -> matrix[dimension_data - 1][0], dimension_data - 3, 3) * 2 + 
		B_spline_function(data_t -> matrix[dimension_data - 1][0] + h, dimension_data - 3, 3)
		) / pow(h, 2);
	parameters_linear_system.matrix[dimension_data + 1][dimension_data] = (
		B_spline_function(data_t -> matrix[dimension_data - 1][0] - h, dimension_data - 2, 3) - 
		B_spline_function(data_t -> matrix[dimension_data - 1][0], dimension_data - 2, 3) * 2 + 
		B_spline_function(data_t -> matrix[dimension_data - 1][0] + h, dimension_data - 2, 3)
		) / pow(h, 2);
	parameters_linear_system.matrix[dimension_data + 1][dimension_data + 1] = (
		B_spline_function(data_t -> matrix[dimension_data - 1][0] - h, dimension_data - 1, 3) - 
		B_spline_function(data_t -> matrix[dimension_data - 1][0], dimension_data - 1, 3) * 2 + 
		B_spline_function(data_t -> matrix[dimension_data - 1][0] + h, dimension_data - 1, 3)
		) / pow(h, 2);

	parameters_linear_system.Matrix_approximation();
	// parameters_linear_system.Print_matrix();
	// b.Print_matrix();

	parameters = parameters_linear_system.Inverse() * b;

	Interpolation_t = new Matrix(m, 1);
	Interpolation_y = new Matrix(m, 1);
	Interpolation_y -> Generate_Zero();
	for (int j = 0; j < m; j++) {
		double t_j = 2.0 * j / (m - 1.0) - 1.0;
		Interpolation_t -> matrix[j][0] = t_j;

		for (int k = 0; k < dimension_data + 2; ++k) {
			Interpolation_y -> matrix[j][0] += parameters.matrix[k][0] * B_spline_function(t_j, k - 2, 3);
		}
	}


}

#endif