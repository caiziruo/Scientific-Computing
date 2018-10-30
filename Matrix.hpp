#ifndef Matrix_hpp
#define Matrix_hpp

#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>
#include <vector>
#include <stack>
#include <queue>

using namespace std;


// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Matrix
{
private:
    
public:
    int row_dimension;
    int column_dimension;
    double ** matrix;
    double ** factorization_L;
    double ** factorization_U;
    double ** factorization_Q;
    double ** factorization_R;

    vector<int> permutation_index;
    
    Matrix(int m, int n) { // Construction function, construct a mxn matrix.
        row_dimension = m;
        column_dimension = n;
        
        // Initialize matrix
        matrix = new double * [row_dimension];
        for (int i = 0; i < row_dimension; i++) {
            matrix[i] = new double[column_dimension];
        }
        
        factorization_L = NULL;
        factorization_U = NULL;
        factorization_Q = NULL;
        factorization_R = NULL;
    }
    
    Matrix(Matrix const & B) { // Construction function, construct from the existed matrix B.
        row_dimension = B.row_dimension;
        column_dimension = B.column_dimension;
        
        // Initialize matrix and copy B
        matrix = new double * [row_dimension];
        for (int i = 0; i < row_dimension; i++) {
            matrix[i] = new double[column_dimension];
            for (int j = 0; j < column_dimension; j++) {
                matrix[i][j] = B.matrix[i][j];
            }
        }
        
        factorization_L = NULL;
        factorization_U = NULL;
        factorization_Q = NULL;
        factorization_R = NULL;
    }

    Matrix(double ** & pointer_A, int m, int n) { // Construction function, constructed from array pointer_A.
        row_dimension = m;
        column_dimension = n;

        // Initialize matrix and copy B
        matrix = new double * [row_dimension];
        for (int i = 0; i < row_dimension; i++) {
            matrix[i] = new double[column_dimension];
            for (int j = 0; j < column_dimension; j++) {
                matrix[i][j] = pointer_A[i][j];
            }
        }
        
        factorization_L = NULL;
        factorization_U = NULL;
        factorization_Q = NULL;
        factorization_R = NULL;
    }
    
    void Set_Element(int i, int j, double x) {matrix[i][j] = x;}
    // double Matrix_element(int i, int j) const {return matrix[i][j]; }
    // int Matrix_row_dimension() const {return row_dimension;}
    // int Matrix_column_dimension() const {return column_dimension;}
    
    void Interchange_row(int m, int n);
    
    void Matrix_approximation(double epsilon);
    void Generate_Identity();
    void Generate_Random();
    void Generate_Nonsingular();
    void Generate_Positive_Definite();
    void Generate_Diagonally_Dominant();

    double Frobenius_norm() const;
    
    void Print_Permutation() const;
    void Print_matrix() const;   // Print all the elements of the matrix
    void Print_L() const;
    void Print_U() const;

    void LU_factorization(bool use_pivoting = true);
    void QR_factorization(bool use_pivoting = true);
    Matrix Transpose() const;
    Matrix inverse_upper_triangular(); 
    Matrix LU_factorization_L();
    Matrix LU_factorization_U();
    Matrix QR_factorization_Q();
    Matrix QR_factorization_R();
    
    // ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Matrix operator+(const Matrix & B) {
        Matrix tmp(row_dimension, column_dimension);
        
        for (int i = 0; i < row_dimension; i++) {
            for (int j = 0; j < column_dimension; j++) {
                tmp.matrix[i][j] = matrix[i][j] + B.matrix[i][j];
            }
        }
        return tmp;
    }
    
    Matrix operator-(const Matrix & B) {
        Matrix tmp(row_dimension, column_dimension);
        
        for (int i = 0; i < row_dimension; i++) {
            for (int j = 0; j < column_dimension; j++) {
                tmp.matrix[i][j] = matrix[i][j] - B.matrix[i][j];
            }
        }
        return tmp;
    }
    
    Matrix operator*(const Matrix & B) {
        Matrix tmp(row_dimension, B.column_dimension);
        
        for (int i = 0; i < row_dimension; i++) {
            for (int j = 0; j < B.column_dimension; j++) {
                double sum = 0;
                for (int k = 0; k < column_dimension; k++) {
                    sum += matrix[i][k] * B.matrix[k][j];
                }
                tmp.matrix[i][j] = sum;
            }
        }
        return tmp;
    }
    
    Matrix operator*(double scalar) { // matrix * scalar
        Matrix tmp(row_dimension, column_dimension);
        
        for (int i = 0; i < row_dimension; i++) {
            for (int j = 0; j < column_dimension; j++) {
                tmp.matrix[i][j] = matrix[i][j] * scalar;
            }
        }
        return tmp;
    }
    
    Matrix & operator=(const Matrix & B) {
        if (this == & B) return * this;// = itself
        
        for (int i = 0; i < row_dimension; i++) {
            delete [] matrix[i];
        }
        delete [] matrix;
        
        //initialize a new matrix
        row_dimension = B.row_dimension;
        column_dimension = B.column_dimension;
        matrix = new double * [row_dimension];
        for (int i = 0; i < row_dimension; i++) {matrix[i] = new double[column_dimension];}
        
        for (int i = 0; i < row_dimension; i++) {
            for (int j = 0; j < column_dimension; j++) {
                matrix[i][j] = B.matrix[i][j];
            }
        }
        
        return * this;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ~Matrix() {
        for (int i = 0; i < row_dimension; ++i) {
            delete [] matrix[i];
        }
        delete [] matrix;
    }
};

Matrix Matrix::Transpose() const {
    Matrix tmp(column_dimension, row_dimension);
    
    for (int i = 0; i < column_dimension; i++) {
        for (int j = 0; j < row_dimension; j++) {
            tmp.matrix[i][j] = matrix[j][i];
        }
    }
    
    return tmp;
}

double Matrix::Frobenius_norm() const {
    double norm = 0;
    
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            norm +=  pow(abs(matrix[i][j]), 2);
        }
    }
    
    norm = pow(norm, 0.5);
    return norm;
}

void Matrix::Interchange_row(int m, int n) {
    if (m == n) return;
    double tmp;
    
    for (int i = 0; i < row_dimension; i++) {
        tmp = matrix[m][i];
        matrix[m][i] = matrix[n][i];
        matrix[n][i] = tmp;
    }
}

void Matrix::Print_Permutation() const {
    for (int i = 0; i < permutation_index.size(); i++) {
        cout << permutation_index[i] << ' ';
        if (i % 2 == 1) {
            cout << ";";
        }
    }
}


void Matrix::Print_matrix() const {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            cout << setw(14) <<  matrix[i][j];
        }
        cout << '\n';
    }
}

void Matrix::Print_L() const {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            cout << setw(14) <<  factorization_L[i][j];
        }
        cout << '\n';
    }
}

void Matrix::Print_U() const {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            cout << setw(14) <<  factorization_U[i][j];
        }
        cout << '\n';
    }
}

void Matrix::Matrix_approximation(double epsilon = 1e-10) {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            if (abs(matrix[i][j]) < epsilon) {matrix[i][j] = 0; }
        }
    }
}

void Matrix::Generate_Identity() {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            if (i == j) matrix[i][j] = 1;
            else matrix[i][j] = 0;
        }
    }
}

void Matrix::Generate_Random() {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            matrix[i][j] = (double)rand()/ RAND_MAX;
        }
    }
}

void Matrix::Generate_Nonsingular() {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            matrix[i][j] = (double)rand()/ RAND_MAX;
        }
    }
}

void Matrix::Generate_Positive_Definite() {
    Matrix random_matrix(row_dimension, row_dimension);
    random_matrix.Generate_Random();
    
    random_matrix = random_matrix.Transpose() * random_matrix;
    
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            matrix[i][j] = random_matrix.matrix[i][j];
        }
    }
}

void Matrix::Generate_Diagonally_Dominant() {
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            if (i == j) matrix[i][j] = row_dimension;
            else matrix[i][j] = 1;
        }
    }
}

Matrix Matrix::inverse_upper_triangular() {
    Matrix tmp(column_dimension, column_dimension);
    tmp.Generate_Identity();

    for (int i = column_dimension - 2; i >= 0; i--) {
        for (int j = i + 1; j < column_dimension; j++) {
            double aj = matrix[i][j] / matrix[j][j];
            for (int k = j; k < column_dimension; k++) {
                tmp.matrix[i][k] -= aj * tmp.matrix[j][k];
            }
        }
    }
    for (int i = 0; i < column_dimension; i++) {
        for (int j = i; j < column_dimension; ++j) {
            tmp.matrix[i][j] = tmp.matrix[i][j] / matrix[i][i];
        }
    }

    return tmp;
}

Matrix Matrix::LU_factorization_L() {
    Matrix tmp(row_dimension, column_dimension);
    
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            tmp.matrix[i][j] = factorization_L[i][j];
        }
    }
    
    return tmp;
}

Matrix Matrix::LU_factorization_U() {
    Matrix tmp(row_dimension, column_dimension);
    
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            tmp.matrix[i][j] = factorization_U[i][j];
        }
    }
    
    return tmp;
}

Matrix Matrix::QR_factorization_Q() {
    Matrix tmp(row_dimension, column_dimension);
    
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            tmp.matrix[i][j] = factorization_Q[i][j];
        }
    }
    
    return tmp;
}

Matrix Matrix::QR_factorization_R() {
    Matrix tmp(column_dimension, column_dimension);
    
    for (int i = 0; i < column_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            tmp.matrix[i][j] = factorization_R[i][j];
        }
    }
    
    return tmp;
}

void Matrix::LU_factorization(bool use_pivoting) {
    if (factorization_L == NULL) {
        factorization_L = new double * [row_dimension];
        for (int i = 0; i < row_dimension; i++) {
            factorization_L[i] = new double [row_dimension];
        }
    }
    if (factorization_U == NULL) {
        factorization_U = new double * [row_dimension];
        for (int i = 0; i < row_dimension; i++) {
            factorization_U[i] = new double [row_dimension];
        }
    }
    while (permutation_index.size() > 0) {
        permutation_index.pop_back();
    }
    
    // copy the original matrix to U and initialize L with identity matrix
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < row_dimension; j++) {
            factorization_U[i][j] = matrix[i][j];
            if (i == j) factorization_L[i][j] = 1;
            else factorization_L[i][j] = 0;
        }
    }
    
    // compute the factorization of A.
    for (int i = 0; i < column_dimension - 1; i++) {
        if (use_pivoting) {
            int biggest_pivot_index = i;
            double biggest_pivot = abs(matrix[i][i]);
            
            for (int j = i + 1; j < row_dimension; j++) {
                if (abs(matrix[j][i]) > biggest_pivot) {
                    biggest_pivot_index = j;
                    biggest_pivot = abs(matrix[j][i]);
                }
            }
            
            if (biggest_pivot_index != i) { // record which two rows are interchanged.
                permutation_index.push_back(i);
                permutation_index.push_back(biggest_pivot_index);
            }
            
            double tmp;  // interchange i, biggest_pivot_index-th row of U
            for (int k = i; k < column_dimension; k++) {
                tmp = factorization_U[biggest_pivot_index][k];
                factorization_U[biggest_pivot_index][k] = factorization_U[i][k];
                factorization_U[i][k] = tmp;
            }
        }
        else if (factorization_U[i][i] == 0) {
            int j = i;
            for (j = i; j < row_dimension; j++) {   //find the first nonzero pivot
                if (matrix[j][i] != 0) {break;}
            }
            
            permutation_index.push_back(i);
            permutation_index.push_back(j);
            
            double tmp;  // interchange i, j-th row of U
            for (int k = i; k < column_dimension; k++) {
                tmp = factorization_U[j][k];
                factorization_U[j][k] = factorization_U[i][k];
                factorization_U[i][k] = tmp;
            }
        }
        
        //////////////////////     Elimination
        for (int l = i + 1; l < row_dimension; l++) {
            double m_li = factorization_U[l][i] / factorization_U[i][i];
            factorization_L[l][i] = m_li;
            
            // Do the elimination
            for (int k = i + 1; k < column_dimension; k++) {
                factorization_U[l][k] -= m_li * factorization_U[i][k];
            }
            factorization_U[l][i] = 0;
        }
    }
}

void Matrix::QR_factorization(bool use_pivoting) {
    Matrix copy_A(matrix, row_dimension, column_dimension);  //  R
    Matrix copy_Q(copy_A); //  Q 
    copy_Q.Generate_Identity();

    for (int i = 0; i < column_dimension; i++) {  // QR factorization. Q and R are stored in copy_A.
        // int biggest_pivot_index = i;
        // if (use_pivoting) {
        //     double biggest_pivot = abs(copy_A.matrix[i][i]);
        //     for (int j = i + 1; j < column_dimension; j++) { // find the biggest pivot
        //         if (abs(copy_A.matrix[j][i]) > biggest_pivot) {
        //             biggest_pivot_index = j;
        //             biggest_pivot = abs(copy_A.matrix[j][i]);
        //         }
        //     }
        // }
        // else {
        //     int biggest_pivot_index = i;
        //     if (copy_A.matrix[i][i] == 0) {
        //         for (int j = i + 1; j < column_dimension; j++) {
        //             if (copy_A.matrix[j][i] != 0){
        //                 biggest_pivot_index = j;
        //                 break;
        //             }
        //         }
        //     }
        // }
        
        // //  interchange the rows of copy_A and copy_b if needed.
        // if (biggest_pivot_index != i) {
        //     double tmp;
        //     for (int k = i; k < column_dimension; k++) {
        //         tmp = copy_A.matrix[i][k]; // interchange copy_A
        //         copy_A.matrix[i][k] = copy_A.matrix[biggest_pivot_index][k];
        //         copy_A.matrix[biggest_pivot_index][k] = tmp;
        //     }
        // }
        
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
    }

    Matrix reduced_Q(row_dimension, column_dimension);  // Q: mxn
    Matrix reduced_R(column_dimension, column_dimension); // R: nxn
    Matrix inverse_R(column_dimension, column_dimension);
    inverse_R.Generate_Identity();
    for (int i = 0; i < column_dimension; ++i) {
        for (int j = 0; j < column_dimension; ++j) {
            reduced_R.matrix[i][j] = copy_A.matrix[i][j];
        }
    }

    inverse_R = reduced_R.inverse_upper_triangular();

    // compute Q = A * R^{-1}
    Matrix original_A(matrix, row_dimension, column_dimension);
    reduced_Q = original_A * inverse_R;


    if (factorization_Q == NULL) {
        factorization_Q = new double * [row_dimension]; 
        for (int i = 0; i < row_dimension; i++) {
            factorization_Q[i] = new double [column_dimension]; 
        }
    }
    for (int i = 0; i < row_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            factorization_Q[i][j] = reduced_Q.matrix[i][j];
        }
    }

    if (factorization_R == NULL) {
        factorization_R = new double * [column_dimension];
        for (int i = 0; i < column_dimension; i++) {
            factorization_R[i] = new double [column_dimension];
        }
    }
    for (int i = 0; i < column_dimension; i++) {
        for (int j = 0; j < column_dimension; j++) {
            factorization_R[i][j] = reduced_R.matrix[i][j];
        }
    }

}

#endif