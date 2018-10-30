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
    int m = 9;
    int n = 5;
    Matrix A(m, n);
    A.Generate_Nonsingular();
    // A.Generate_Positive_Definite();
    // A.Generate_Diagonally_Dominant();
    

    Matrix Q(m, n);
    Matrix R(n, n);
    
    A.QR_factorization(true);
    
    cout << "\nMatrix Q:\n";
    Q = A.QR_factorization_Q();
    Q.Print_matrix();
    
    cout << "\nMatrix R:\n";
    R = A.QR_factorization_R();
    R.Print_matrix();
    
    cout << "\nMatrix Q * R:\n";
    Matrix D(n, n);
    D = Q * R;
    D.Print_matrix();
    
    cout << "\nMatrix A:\n";
    A.Print_matrix();

    cout << "\nQR - A:\n";
    A = D - A;
    A.Print_matrix();

}
