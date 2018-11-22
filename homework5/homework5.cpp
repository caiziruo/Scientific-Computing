#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "../Matrix.hpp"
#include "../Conjugate_Gradient.hpp"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
	int n = 5;
	Matrix A(n, n);
	Matrix b(n, 1);
	Matrix x0(n, 1);
	A.Generate_Positive_Definite();
	// A.Generate_Diagonally_Dominant();
	b.Generate_Random();
	// x0.Generate_Random();
	x0.Generate_Zero();

	cout << "A:\n";
	A.Print_matrix();
	cout << "\nb:\n";
	b.Print_matrix();
	cout << "\nstart from:\n";
	x0.Print_matrix();


	Conjugate_Gradient cg_method(n, A, b, x0);
	cg_method.compute_minimization();

	cout << "\nsolution:\n";
	cg_method.Print_solution();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	n = 2;
	Matrix AA(n, n);
	Matrix bb(n, 1);
	Matrix xx(n, 1);
	AA.matrix[0][0] = 0.5;
	AA.matrix[1][1] = 2.5;
	xx.matrix[0][0] = 5;
	xx.matrix[1][0] = 1;
	bb.Generate_Zero();

	cout << "\nAA:\n";
	AA.Print_matrix();
	cout << "\nbb:\n";
	bb.Print_matrix();
	cout << "\nstart from:\n";
	xx.Print_matrix();

	Conjugate_Gradient cg_method_2(n, AA, bb, xx);
	cg_method_2.compute_minimization();

	cout << "\nsolution:\n";
	cg_method_2.Print_solution();





	return 0;
}