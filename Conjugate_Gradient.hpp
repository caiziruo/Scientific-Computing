#ifndef Conjugate_Gradient_hpp
#define Conjugate_Gradient_hpp


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

// function_f is defined by user.
// double function_f(Matrix x) {
// 	double fx;
// 	f = 0.5 * pow(x.matrix[0][0], 2) + 2.5 * pow(x.matrix[1][0], 2);

// 	return f;
// }

class Conjugate_Gradient
{
public:
	// x = (x1, x2, ... , xn).T
	// conjugate gradient method for solving the minimization problem of f(x) = 0.5 * x.T A x - b.T x + c

	int solution_dimension;  //  n = solution_dimension
    Matrix * solution; // x
    Matrix * A;
    Matrix * b;
    
	Conjugate_Gradient(int n, Matrix AA, Matrix bb, Matrix start) {
		solution_dimension = n;

		A = new Matrix(AA);
		b = new Matrix(bb);
		solution = new Matrix(start);
	}

	void compute_minimization();

	void Print_solution() const {
		for (int i = 0; i < solution_dimension; ++i) {
			cout << setw(14) << solution -> matrix[i][0] << '\n';
		}
	}

	~Conjugate_Gradient() {
		delete solution;
	}
};

void Conjugate_Gradient::compute_minimization() {
	double alpha; // line search parameter
	double beta;
	Matrix x(solution_dimension, 1);  // initial x_{0} and then x_{k}
	Matrix x_next(solution_dimension, 1);  // x_{k+1}
	Matrix g(solution_dimension, 1);  // g_{0} = \nabla f(x_{0})
	Matrix g_next(solution_dimension, 1);
	Matrix s(solution_dimension, 1);  
	Matrix s_next(solution_dimension, 1);
	 
	// initialization
	x = (*solution);
	// x.Generate_Random();
	g = (*A) * x - (*b);
	if (g.Frobenius_norm() < 1e-10 * pow(solution_dimension, 0.5)) {
		(*solution) = x;
		return;
	}
	s = g.Opposite();

	for (int k = 0; k < 1000; k++) {
		// compute line search parameter alpha
		alpha = (s.Transpose() * ((*b) - (*A) * x) * 
			pow((s.Transpose() * (*A) * s).matrix[0][0], -1) 
			).matrix[0][0];

		x_next = x + s * alpha;
		if ((x_next - x).Frobenius_norm() < 1e-10 * pow(solution_dimension, 0.5)) {break; }

		g_next = (*A) * x_next - (*b);
		if (g_next.Frobenius_norm() < 1e-10 * pow(solution_dimension, 0.5)) {x = x_next; break; }
		beta = ((g_next.Transpose() * g_next) * pow((g.Transpose() * g).matrix[0][0], -1)).matrix[0][0];
		s_next = s * beta - g_next;

		x = x_next;
		g = g_next;
		s = s_next;
	}

	(*solution) = x;


}



#endif