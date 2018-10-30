#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "Matrix.hpp"


int main() {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int n = 9;
    Matrix A(n, n);
    // A.Generate_Nonsingular();
    A.Generate_Positive_Definite();
    // A.Generate_Diagonally_Dominant();

    Matrix eigenvectors(n, n);
    eigenvectors.Generate_Identity();
    Matrix Q(n, n);
    Matrix R(n, n);
    

    // A.Print_matrix();
    Matrix Ak(A);
    cout << '\n';
    for (int i = 0; i < 1000; ++i) {
        if (i % 10 == 0) {
            Ak.Matrix_approximation(1e-10);
        }
        Ak.QR_factorization(true);
        Q = Ak.QR_factorization_Q();
        R = Ak.QR_factorization_R();
        Ak = R * Q;
        eigenvectors = eigenvectors * Q;
    }

    Ak.Print_matrix();
    // cout << '\n';
    // A = (eigenvectors.Transpose() * A) * eigenvectors;
    // A.Matrix_approximation(1e-10);
    // A.Print_matrix();

}
