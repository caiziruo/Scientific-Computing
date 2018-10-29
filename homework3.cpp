#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "Matrix.hpp"
#include "Linear_System.hpp"


int main() {
    int m = 21;
    int n = 12;
    
    Linear_system linear_equation(m, n); // Initialize a mxn linear system for least square problem.
    for (int i = 0; i < m; i++) { // Set A
        for (int j = 0; j < n; j++) {
            linear_equation.Set_A_element(i, j, pow(i / (m - 1.0), j));
        }
    }

    Matrix initial_solution(n, 1);
    for (int i = 0; i < n; i++) {
        initial_solution.Set_Element(i, 0, 1);
    }

    for (int i = 0; i < m; i++) { // Set b
        double y_i = 0;
        for (int j = 0; j < n; j++) {
            y_i += linear_equation.Matrix_A_element(i, j);
        }
        linear_equation.Set_b_element(i, y_i);
    }
    // linear_equation.Show_system();
    // cout << '\n';

    double perturbation_epsilon = pow(10, -10);
    for (int i = 0; i < m; i++) {
        double tmp = (double)rand()/ RAND_MAX;
        linear_equation.Set_b_element(i, linear_equation.Vector_b_element(i) + (2 * tmp - 1) * perturbation_epsilon);
    }


    Linear_system normal_equation_method(n, n);
    Matrix A(m, n);
    Matrix AT(n, m);
    Matrix ATA(n, n);
    Matrix ATb(n, 1);
    Matrix perturbated_b(m, 1);

    A = linear_equation.Matrix_A();
    AT = linear_equation.Matrix_A().Transpose();
    ATA = AT * AT.Transpose();                 //  define matrix A.T * A
    normal_equation_method.Set_A(ATA);
    perturbated_b = linear_equation.Vector_b();

    ATb = AT * perturbated_b;
    normal_equation_method.Set_b(ATb);

    // normal_equation_method.Show_system();
    normal_equation_method.Solve_From_Cholesky();
    cout << "Normal solution method(Cholesky):\n";
    normal_equation_method.Show_Solution();
    cout << '\n';

    linear_equation.Solve_From_QR();
    cout << "Solution using QR factorization:" << "\n";
    linear_equation.Show_Solution();
    // cout << '\n';
    // linear_equation.Show_system();

    return 0;

}