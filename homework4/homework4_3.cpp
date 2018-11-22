#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>
#include <string>

using namespace std;

#include "../Matrix.hpp"

Matrix function_f(Matrix x) {
	int m = 2;
	Matrix f(m, 1);

	f.matrix[0][0] = x.matrix[0][0] + 2 * x.matrix[1][0] - 2;
	f.matrix[1][0] = pow(x.matrix[0][0], 2) + 4 * pow(x.matrix[1][0], 2) - 4;

	// ...
	// f[m - 1][0] = 

	return f;
}

Matrix Jacobian_matrix(Matrix x) {
	int m = 2;
	Matrix J(m, x.row_dimension);

	J.matrix[0][0] = 1;
	J.matrix[0][1] = 2;
	J.matrix[1][0] = 2 * x.matrix[0][0];
	J.matrix[1][1] = 8 * x.matrix[1][0];

	return J;
}

class Non_Linear_System
{
public:
	// x = (x1, x2, ... , xn).T
	// f(x) = (f1(x), f2(x), ... , fm(x)).T
	int solution_dimension;  //  n = solution_dimension
    Matrix * solution; // x
    
	Non_Linear_System(int n) {
		solution_dimension = n;

		solution = new Matrix(solution_dimension, 1);
		solution -> Generate_Zero();
	}

	void Newton_method();
	void Broyden_method();
	void Print_solution() const {
		for (int i = 0; i < solution_dimension; ++i) {
			cout << setw(14) << solution -> matrix[i][0] << '\n';
		}
	}
	
	~Non_Linear_System() {
		delete solution;
	}
};

void Non_Linear_System::Newton_method() {
	Matrix x(solution_dimension, 1);
	Matrix x_next(solution_dimension, 1);
	x.Generate_Random();

	for (int i = 0; i < 1000; ++i) {
		Matrix J(Jacobian_matrix(x));
		Matrix f(function_f(x));

		// J.Print_matrix();
		// cout << '\n';
		J = J.Inverse();
		// J.Print_matrix();

		// f.Print_matrix();
		// cout << '\n';

		x_next = x - J * f;
		if ((x_next - x).Frobenius_norm() < 1e-10 * pow(solution_dimension, 0.5) ) break;
		x = x_next;
	}
	*solution = x;
}

void Non_Linear_System::Broyden_method() {
	Matrix x(solution_dimension, 1);
	Matrix x_next(solution_dimension, 1);
	Matrix B(solution_dimension, solution_dimension);
	Matrix B_next(solution_dimension, solution_dimension);
	x.Generate_Random();
	B = Jacobian_matrix(x);

	Matrix f(function_f(x));
	Matrix s(solution_dimension, 1);
	Matrix y(solution_dimension, 1);
	for (int i = 0; i < 1000; ++i) {
		f = function_f(x);
		s = B.Inverse() * f.Opposite();
		x_next = x + s;
		if ((x_next - x).Frobenius_norm() < 1e-10 * pow(solution_dimension, 0.5) ) break;

		y = function_f(x_next) - function_f(x);
		B_next = B + ((y - B * s) * s.Transpose()) * pow((s.Transpose() * s).matrix[0][0], -1);

		x = x_next;
		B = B_next;
	}
	*solution = x;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
	Non_Linear_System non_linear_equation(2);
	non_linear_equation.Newton_method();
	cout << "Newton's method solution:\n";
	non_linear_equation.Print_solution();

	non_linear_equation.Broyden_method();
	cout << "Broyden's method solution:\n";
	non_linear_equation.Print_solution();

	return 0;
}
