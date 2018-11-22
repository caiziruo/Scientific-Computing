#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "../Matrix.hpp"


int main() {
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int n = 9;
    Matrix A(n, n);
    // A.Generate_Nonsingular();
    A.Generate_Positive_Definite();
    // A.Generate_Diagonally_Dominant();
    
    Matrix B(n, n);
    Matrix C(n, n);
    C = A;
    
    C.LU_factorization(true);
    
    cout << "Matrix L:\n";
    A = C.LU_factorization_L();
    A.Print_matrix();
    
    cout << "Matrix U:\n";
    B = C.LU_factorization_U();
    B.Print_matrix();
    
    cout << "Matrix L * U:\n";
    Matrix D(n, n);
    D = A * B;
    D.Print_matrix();
    
    cout << "Matrix A:\n";
    C.Print_matrix();

}