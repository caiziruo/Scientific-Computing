#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "Matrix.hpp"
#include "Eigenvalue.hpp"


int main() {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int n = 9;
    Matrix A(n, n);
    // A.Generate_Nonsingular();
    // A.Generate_Positive_Definite();
    A.Generate_Diagonally_Dominant();



    Eigenvalue eigenvalues_computation(A, 9, 9);
    eigenvalues_computation.QR_iteration(1000);

    cout << "Eigenvalues:\n";
    eigenvalues_computation.Print_eigenvalues();
    cout << "\nEigenvectors:\n";
    eigenvalues_computation.Print_eigenvectors();

}
