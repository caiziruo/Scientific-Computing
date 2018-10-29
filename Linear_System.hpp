#ifndef Linear_System_hpp
#define Linear_System_hpp

#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;

#include "Matrix.hpp"

class Linear_system  //Solving Ax = b.
{
private:
    
public:
    int row_dimension;
    int column_dimension;
    Matrix * matrix_A; // A
    Matrix * solution; // x
    Matrix * vector_b; // b
    Matrix * factorization_L; // L
    Matrix * factorization_U; // U
    bool use_pivoting;
    
    Linear_system(int m, int n) {
        row_dimension = m;
        column_dimension = n;
        
        // Initialize matrix A, x, b.
        matrix_A = new Matrix(row_dimension, column_dimension);
        vector_b = new Matrix(row_dimension, 1);
        
        for (int i = 0; i < row_dimension; i++) { // Initialize b
            vector_b -> matrix[i][0] = 1.0; // * i;
        }
        
        solution = new Matrix(column_dimension, 1);
        factorization_L = NULL;
        factorization_U = NULL;
        use_pivoting = false;
    }
    
    void Set_A_element(int i, int j, double x) {matrix_A -> matrix[i][j] = x;}
    void Set_b_element(int i, double x) {vector_b -> matrix[i][0] = x;}
    void Set_A(Matrix & tmp);
    void Set_b(Matrix & tmp);
    
    void LU_Use_Pivoting(bool use_or_not) {use_pivoting = use_or_not;}
    void Generate_Nonsingular() {matrix_A -> Generate_Nonsingular();}
    void Generate_Diagonally_Dominant() {matrix_A -> Generate_Diagonally_Dominant();}
    void Generate_Positive_Definite() {matrix_A -> Generate_Positive_Definite();}
    void Solve_From_Cholesky();
    void Solve_From_LU();
    void Solve_From_QR(); // A = mxn, m > n
    void Solve(); // Use Gaussian elimination to get a solution;
    void Show_Solution() const {solution -> Print_matrix();}
    double Matrix_A_element(int i, int j) const {return matrix_A -> matrix[i][j];}
    double Vector_b_element(int i) const {return vector_b -> matrix[i][0];}
    Matrix Matrix_A() const;
    Matrix Vector_b() const;
    Matrix Solution_Matrix() const;
    void Show_system() {
        cout << "A:\n";
        matrix_A -> Print_matrix();
        cout << "b:\n";
        vector_b -> Print_matrix();
    }
    
    ~Linear_system() {
        delete matrix_A;
        delete vector_b;
        delete solution;
    };
};



void Linear_system::Solve_From_LU() {
    matrix_A -> LU_factorization(use_pivoting);
    
    if (factorization_L != NULL) {
        delete factorization_L;
    }
    if (factorization_U != NULL) {
        delete factorization_U;
    }
    
    factorization_L = new Matrix(matrix_A -> LU_factorization_L());
    factorization_U = new Matrix(matrix_A -> LU_factorization_U());
    Matrix * tmp_b = new Matrix(*vector_b);
    
    Matrix * solution_y = new Matrix(row_dimension, 1);
    
    queue<int> permutation_index_queue;
    if (true) {
        for (int i = 0; i < matrix_A -> permutation_index.size(); i++) {
            permutation_index_queue.push(matrix_A -> permutation_index[i]);
        }
        double tmp;
        int i;
        int j;
        while (permutation_index_queue.size() > 0) {
            i = permutation_index_queue.front();
            permutation_index_queue.pop();
            j = permutation_index_queue.front();
            permutation_index_queue.pop();
            
            // Permutation of vector b
            tmp = tmp_b -> matrix[i][0];
            tmp_b -> matrix[i][0] = tmp_b -> matrix[j][0];
            tmp_b -> matrix[j][0] = tmp;
        }
    }
    
    for (int i = 0; i < row_dimension; i++) { //solve Ly = b
        double tmp = 0;
        for (int j = 0; j < i; j++) {
            tmp += factorization_L -> matrix[i][j] * solution_y -> matrix[j][0];
        }
        solution_y -> matrix[i][0] = (tmp_b -> matrix[i][0] - tmp) / factorization_L -> matrix[i][i];
    }
    
    for (int i = row_dimension - 1; i >= 0; i--) {
        double tmp = 0;
        for (int j = row_dimension - 1; j > i; j--) {
            tmp += factorization_U -> matrix[i][j] * solution -> matrix[j][0];
        }
        solution -> matrix[i][0] = (solution_y -> matrix[i][0] - tmp) / factorization_U -> matrix[i][i];
    }
}

void Linear_system::Solve_From_Cholesky() {
    Matrix copy_A(*matrix_A);
    Matrix L(*matrix_A);
    Matrix copy_b(*vector_b);
    Matrix solution_y(*vector_b);

    L.Generate_Identity();
    
    for (int i = 0; i < row_dimension; i++) {  // Compute L
        for (int j = i; j < row_dimension; j++) {
            double tmp = 0;
            for (int k = 0; k < i; k++) {
                tmp += L.matrix[i][k] * L.matrix[j][k];
            }
            if (j == i) {
                L.matrix[i][i] = pow(copy_A.matrix[j][i] - tmp, 0.5);
            }
            else {
                L.matrix[j][i] = (copy_A.matrix[j][i] - tmp) / L.matrix[i][i];
            }
        }
    }
    
    // solve x from Ly = b and L.T x = y
    for (int i = 0; i < row_dimension; i++) {
        double tmp = 0;
        for (int j = 0; j < i; j++) {
            tmp += L.matrix[i][j] * solution_y.matrix[j][0];
        }
        solution_y.matrix[i][0] = (copy_b.matrix[i][0] - tmp) / L.matrix[i][i];
    }

    for (int i = row_dimension - 1; i >= 0; i--) {
        double tmp = 0;
        for (int j = i + 1; j < column_dimension; j++) {
            tmp += L.matrix[j][i] * copy_b.matrix[j][0];
        }
        copy_b.matrix[i][0] = (solution_y.matrix[i][0] - tmp) / L.matrix[i][i];
    }
    
    if (solution != NULL) {
        delete solution;
    }
    solution = new Matrix(copy_b);
}

void Linear_system::Solve() {
    Matrix copy_A(*matrix_A);
    Matrix copy_b(*vector_b);
    
    for (int i = 0; i < row_dimension - 1; i++) {  // LU factorization. L and U are stored in copy_A.
        int biggest_pivot_index = i;
        if (use_pivoting) {
            double biggest_pivot = abs(copy_A.matrix[i][i]);
            for (int j = i + 1; j < row_dimension; j++) { // find the biggest pivot
                if (abs(copy_A.matrix[j][i]) > biggest_pivot) {
                    biggest_pivot_index = j;
                    biggest_pivot = abs(copy_A.matrix[j][i]);
                }
            }
        }
        else {
            int biggest_pivot_index = i;
            if (copy_A.matrix[i][i] == 0) {
                for (int j = i + 1; j < row_dimension; j++) {
                    if (copy_A.matrix[j][i] != 0){
                        biggest_pivot_index = j;
                        break;
                    }
                }
            }
        }
        
        //  interchange the rows of copy_A and copy_b if needed.
        if (biggest_pivot_index != i) {
            double tmp;
            for (int k = i; k < column_dimension; k++) {
                tmp = copy_A.matrix[i][k]; // interchange copy_A
                copy_A.matrix[i][k] = copy_A.matrix[biggest_pivot_index][k];
                copy_A.matrix[biggest_pivot_index][k] = tmp;
            }
            tmp = copy_b.matrix[i][0];  // interchange copy_b
            copy_b.matrix[i][0] = copy_b.matrix[biggest_pivot_index][0];
            copy_b.matrix[biggest_pivot_index][0] = tmp;
        }
        
        // Gaussian elimination to copy_A and copy_b and record L .
        for (int l = i + 1; l < row_dimension; l++) {
            double m_li = copy_A.matrix[l][i] / copy_A.matrix[i][i];
            copy_A.matrix[l][i] = m_li; // record L to lower triangular area of A.
            
            // Do the elimination to copy_A and copy_b.
            for (int k = i + 1; k < column_dimension; k++) {
                copy_A.matrix[l][k] -= m_li * copy_A.matrix[i][k];
            }
            copy_b.matrix[l][0] -= m_li * copy_b.matrix[i][0];
        }
        
    }
    
    // solve from U
    for (int i = row_dimension - 1; i >= 0; i--) {
        for (int j = i + 1; j < column_dimension; j++) {
            copy_b.matrix[i][0] -= copy_A.matrix[i][j] * copy_b.matrix[j][0];
        }
        copy_b.matrix[i][0] /= copy_A.matrix[i][i];
        copy_A.matrix[i][i] = 1;
    }
    
    
    if (solution != NULL) {
        delete solution;
    }
    solution = new Matrix(copy_b);
}

void Linear_system::Solve_From_QR() {
    Matrix copy_A(*matrix_A);  //  R
    Matrix copy_b(*vector_b);
    Matrix copy_Q(copy_A); //  Q  
    copy_Q.Generate_Identity();

    for (int i = 0; i < column_dimension; i++) {  // QR factorization. Q and R are stored in copy_A.
        int biggest_pivot_index = i;
        if (use_pivoting) {
            double biggest_pivot = abs(copy_A.matrix[i][i]);
            for (int j = i + 1; j < column_dimension; j++) { // find the biggest pivot
                if (abs(copy_A.matrix[j][i]) > biggest_pivot) {
                    biggest_pivot_index = j;
                    biggest_pivot = abs(copy_A.matrix[j][i]);
                }
            }
        }
        else {
            int biggest_pivot_index = i;
            if (copy_A.matrix[i][i] == 0) {
                for (int j = i + 1; j < column_dimension; j++) {
                    if (copy_A.matrix[j][i] != 0){
                        biggest_pivot_index = j;
                        break;
                    }
                }
            }
        }
        
        //  interchange the rows of copy_A and copy_b if needed.
        if (biggest_pivot_index != i) {
            double tmp;
            for (int k = i; k < column_dimension; k++) {
                tmp = copy_A.matrix[i][k]; // interchange copy_A
                copy_A.matrix[i][k] = copy_A.matrix[biggest_pivot_index][k];
                copy_A.matrix[biggest_pivot_index][k] = tmp;
            }
            tmp = copy_b.matrix[i][0];  // interchange copy_b
            copy_b.matrix[i][0] = copy_b.matrix[biggest_pivot_index][0];
            copy_b.matrix[biggest_pivot_index][0] = tmp;
        }
        
        // Q R factorization and record Q and R
        double alpha = 0;
        for (int j = i; j < row_dimension; j++) {
            alpha += pow(copy_A.matrix[j][i], 2);
        }
        alpha = pow(alpha, 0.5);  

        if (copy_A.matrix[i][i] > 0) {alpha = -alpha;} // to avoid cancellation


        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //   v_i = copy_A[i:m][i] - alpha * e_i
        Matrix v_i(row_dimension, 1);

        copy_Q.matrix[i][i] = copy_A.matrix[i][i] - alpha;
        v_i.matrix[i][0] = copy_Q.matrix[i][i];
        for (int j = 0; j < i; j++) {
            v_i.matrix[j][0] = 0;
        }
        for (int j = i + 1; j < row_dimension; j++) { // record v_{i} and Q[][i]
            copy_Q.matrix[j][i] = copy_A.matrix[j][i];
            v_i.matrix[j][0] = copy_A.matrix[j][i];
        }

        double v_i_norm  = pow(copy_Q.matrix[i][i], 2);  // compute the norm of v_{i}
        for (int j = i + 1; j < row_dimension; j++) {      // compute the norm of v_{i}
            v_i_norm += pow(copy_Q.matrix[j][i], 2);
        }
        // v_i_norm = pow(v_i_norm, 0.5);
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // compute H{i} a_{i}
        copy_A.matrix[i][i] = alpha;
        for (int j = i + 1; j < row_dimension; j++) {
            copy_A.matrix[j][i] = 0;
        }

        // Hi = I - 2  v_{i} * v_{i}^{T}  / || v_{i} ||_{2}^{2}
        Matrix old_a_i(row_dimension, 1); // old_a_i to record the old a_{k}
        Matrix new_a_i(row_dimension, 1);
        for (int k = i + 1; k < column_dimension; k++) { //  compute  H_{i} a_{k} 
            for (int j = 0; j < row_dimension; j++) {
                old_a_i.matrix[j][0] = copy_A.matrix[j][k];
            }

            // Hi a_{k} = a_{k} - 2  v_{i} * v_{i}^{T} * a_{k}  / || v_{i} ||_{2}^{2}
            new_a_i = old_a_i - v_i * v_i.Transpose() * old_a_i * 2 * (1.0 / v_i_norm);
            for (int j = 0; j < row_dimension; j++) {
                copy_A.matrix[j][k] = new_a_i.matrix[j][0];
            }
        }

        // compute H_{i} b
        new_a_i = copy_b - v_i * v_i.Transpose() * copy_b * 2 * (1.0 / v_i_norm);
        for (int j = 0; j < row_dimension; j++) {
            copy_b.matrix[j][0] = new_a_i.matrix[j][0];
        }
        
    }

    // solve from R
    for (int i = column_dimension - 1; i >= 0; i--) {
        for (int j = i + 1; j < column_dimension; j++) {
            copy_b.matrix[i][0] -= copy_A.matrix[i][j] / copy_A.matrix[j][j] * copy_b.matrix[j][0]; 
        } 
        // copy_b.matrix[i][0] /= copy_A.matrix[i][i]; 
        // copy_A.matrix[i][i] = 1; 
    }
    for (int i = 0; i < column_dimension; i++) {
        solution -> matrix[i][0] = copy_b.matrix[i][0] / copy_A.matrix[i][i];
        // solution -> matrix[i][0] = copy_b.matrix[i][0];
    }
}

void Linear_system::Set_A(Matrix & tmp) {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            matrix_A -> matrix[i][j] = tmp.matrix[i][j];
        }
    }
}

void Linear_system::Set_b(Matrix & tmp) {
    for (int i = 0; i < row_dimension; i++) {
        vector_b -> matrix[i][0] = tmp.matrix[i][0];
    }
}

Matrix Linear_system::Matrix_A() const {
    Matrix tmp(row_dimension, column_dimension);
    
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < row_dimension; j++) {
            tmp.matrix[i][j] = matrix_A -> matrix[i][j];
        }
    }
    
    return tmp;
}

Matrix Linear_system::Vector_b() const {
    Matrix tmp(row_dimension, 1);
    
    for (int i = 0; i < row_dimension; i++) {
        tmp.matrix[i][0] = vector_b -> matrix[i][0];
    }
    
    return tmp;
}


Matrix Linear_system::Solution_Matrix() const {
    Matrix tmp(row_dimension, 1);
    
    for (int i = 0; i < row_dimension; i++) {
        tmp.matrix[i][0] = solution -> matrix[i][0];
    }
    
    return tmp;
}


#endif