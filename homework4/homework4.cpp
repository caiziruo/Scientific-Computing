#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "../Matrix.hpp"
#include "../Eigenvalue.hpp"


int main() {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int n = 5;
    Matrix A(n, n);
    // A.Generate_Nonsingular();
    A.Generate_Positive_Definite();
    // A.Generate_Diagonally_Dominant();

    A.Print_matrix();
    cout << '\n'; 

    Eigenvalue eigenvalues_computation(A);

    eigenvalues_computation.Get_eigenvalues("QR_iteration", false);
    // cout << "Eigenvalues:\n";
    // eigenvalues_computation.Print_eigenvalues();
    // cout << "\nEigenvectors:\n";
    // eigenvalues_computation.Print_eigenvectors();
    // cout << '\n';

    eigenvalues_computation.Get_eigenvalues("QR_iteration", true, "Arnoldi");
    // eigenvalues_computation.Print_reduced_matrix();
    // cout << "Eigenvalues:\n";
    // eigenvalues_computation.Print_eigenvalues();
    // cout << "\nEigenvectors:\n";
    // eigenvalues_computation.Print_eigenvectors();
    // cout << '\n';

    eigenvalues_computation.Get_eigenvalues("QR_iteration", true, "Lanczos");
    eigenvalues_computation.Print_reduced_matrix();
    cout << "Eigenvalues:\n";
    eigenvalues_computation.Print_eigenvalues();
    cout << "\nEigenvectors:\n";
    eigenvalues_computation.Print_eigenvectors();
    cout << '\n';
}
