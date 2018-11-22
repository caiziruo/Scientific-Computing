#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "../Matrix.hpp"
#include "../Linear_System.hpp"


int main() {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int n = 1000;
    Linear_system linear_equation(n, n); // Initialize a nxn linear system.
    
    // linear_equation.Generate_Positive_Definite();    // Generate a positive definite matrix A.
    // linear_equation.Generate_Diagonally_Dominant();   // Generate a diagonally dominant matrix A.
    linear_equation.Generate_Nonsingular();           // Generate a nonsingular matrix A.
    
    // cout << "\nLinear system:\n";
    // linear_equation.Show_system();      // Show A and b.
    
    Matrix matrix_A(n, n);
    Matrix vector_b(n, 0);
    Matrix solution_without_pivoting(n, 0);
    Matrix solution_with_pivoting(n, 0);
    Matrix residual_matrix(n, n);
    
    matrix_A = linear_equation.Matrix_A();
    vector_b = linear_equation.Vector_b();
    
    linear_equation.Solve();
    solution_without_pivoting = linear_equation.Solution_Matrix();  // Compute the solution without pivoting
    // cout << "\nSolution:\n";
    // solution_without_pivoting.Print_matrix();
    
    linear_equation.LU_Use_Pivoting(true);
    linear_equation.Solve();            // Compute the solution with pivoting.
    solution_with_pivoting = linear_equation.Solution_Matrix();
    // cout << "\nSolution(use pivoting):\n";
    // solution_with_pivoting.Print_matrix();
    
    residual_matrix = matrix_A * solution_without_pivoting - vector_b;
    double norm_of_residual_without_pivoting = residual_matrix.Frobenius_norm();
    residual_matrix = matrix_A * solution_with_pivoting - vector_b;
    double norm_of_residual_with_pivoting = residual_matrix.Frobenius_norm();
    
    cout <<
    "Norm of residual(without pivoting): "  << norm_of_residual_without_pivoting    <<  ",\n" <<
    "Norm of residual(with pivoting): "     << norm_of_residual_with_pivoting        << '\n';
    
}
