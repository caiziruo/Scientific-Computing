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
    }
    
    Matrix(Matrix const & B) { // Construction function, construct a existed matrix.
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
    }
    
    void Set_Element(int i, int j, double x) {matrix[i][j] = x;}
    // double Matrix_element(int i, int j) const {return matrix[i][j]; }
    // int Matrix_row_dimension() const {return row_dimension;}
    // int Matrix_column_dimension() const {return column_dimension;}
    // void copy_matrix(Matrix const B); // copy matrix B to this matrix
    
    void Interchange_row(int m, int n);
    
    void Generate_Identity();
    void Generate_Random();
    void Generate_Nonsingular();
    void Generate_Positive_Definite();
    void Generate_Diagonally_Dominant();
    
    void LU_factorization(bool use_pivoting);
    double Norm(int p) const;
    double Frobenius_norm() const;
    
    void Print_Permutation() const;
    void Print_matrix() const;   // Print all the elements of the matrix
    void Print_L() const;
    void Print_U() const;
    Matrix Transpose() const;
    Matrix LU_factorization_L();
    Matrix LU_factorization_U();
    
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
    
    Matrix operator*(double scalar) {
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

double Matrix::Norm(int p) const {
    double norm = 0;
    // to be continued....
    
    return norm;
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


///////////////////////////////////////////////////////////////////////////////////////////




// int main() {
    
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // int m = 21;
    // int n = 12;
    
    // Linear_system linear_equation(m, n); // Initialize a mxn linear system for least square problem.
    // for (int i = 0; i < m; i++) { // Set A
    //     for (int j = 0; j < n; j++) {
    //         linear_equation.Set_A_element(i, j, pow(i / (m - 1.0), j));
    //     }
    // }

    // Matrix initial_solution(n, 1);
    // for (int i = 0; i < n; i++) {
    //     initial_solution.Set_Element(i, 0, 1);
    // }

    // for (int i = 0; i < m; i++) { // Set b
    //     double y_i = 0;
    //     for (int j = 0; j < n; j++) {
    //         y_i += linear_equation.Matrix_A_element(i, j);
    //     }
    //     linear_equation.Set_b_element(i, y_i);
    // }
    // // linear_equation.Show_system();
    // // cout << '\n';

    // double perturbation_epsilon = pow(10, -10);
    // for (int i = 0; i < m; i++) {
    //     double tmp = (double)rand()/ RAND_MAX;
    //     linear_equation.Set_b_element(i, linear_equation.Vector_b_element(i) + (2 * tmp - 1) * perturbation_epsilon);
    // }


    // Linear_system normal_equation_method(n, n);
    // Matrix A(m, n);
    // Matrix AT(n, m);
    // Matrix ATA(n, n);
    // Matrix ATb(n, 1);
    // Matrix perturbated_b(m, 1);

    // A = linear_equation.Matrix_A();
    // AT = linear_equation.Matrix_A().Transpose();
    // ATA = AT * AT.Transpose();                 //  define matrix A.T * A
    // normal_equation_method.Set_A(ATA);
    // perturbated_b = linear_equation.Vector_b();

    // ATb = AT * perturbated_b;
    // normal_equation_method.Set_b(ATb);

    // // normal_equation_method.Show_system();
    // normal_equation_method.Solve_From_Cholesky();
    // cout << "Normal solution method(Cholesky):\n";
    // normal_equation_method.Show_Solution();
    // cout << '\n';

    // linear_equation.Solve_From_QR();
    // cout << "Solution using QR factorization:" << "\n";
    // linear_equation.Show_Solution();
    // // cout << '\n';
    // // linear_equation.Show_system();
    

    
//     return 0;
// }

#endif