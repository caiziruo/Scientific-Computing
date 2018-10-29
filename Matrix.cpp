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



int main() {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // int n = 9;
    // Matrix A(n, n);
    // // A.Generate_Nonsingular();
    // A.Generate_Positive_Definite();
    // // A.Generate_Diagonally_Dominant();
    
    // Matrix B(n, n);
    // Matrix C(n, n);
    // C = A;
    
    // C.LU_factorization(true);
    
    // cout << "Matrix L:\n";
    // A = C.LU_factorization_L();
    // A.Print_matrix();
    
    // cout << "Matrix U:\n";
    // B = C.LU_factorization_U();
    // B.Print_matrix();
    
    // cout << "Matrix L * U:\n";
    // Matrix D(n, n);
    // D = A * B;
    // D.Print_matrix();
    
    // cout << "Matrix A:\n";
    // C.Print_matrix();
    
    // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // int n = 1000;
    // Linear_system linear_equation(n, n); // Initialize a nxn linear system.
    
    // // linear_equation.Generate_Positive_Definite();    // Generate a positive definite matrix A.
    // // linear_equation.Generate_Diagonally_Dominant();   // Generate a diagonally dominant matrix A.
    // linear_equation.Generate_Nonsingular();           // Generate a nonsingular matrix A.
    
    // // cout << "\nLinear system:\n";
    // // linear_equation.Show_system();      // Show A and b.
    
    // Matrix matrix_A(n, n);
    // Matrix vector_b(n, 0);
    // Matrix solution_without_pivoting(n, 0);
    // Matrix solution_with_pivoting(n, 0);
    // Matrix residual_matrix(n, n);
    
    // matrix_A = linear_equation.Matrix_A();
    // vector_b = linear_equation.Vector_b();
    
    // linear_equation.Solve();
    // solution_without_pivoting = linear_equation.Solution_Matrix();  // Compute the solution without pivoting
    // // cout << "\nSolution:\n";
    // // solution_without_pivoting.Print_matrix();
    
    // linear_equation.LU_Use_Pivoting(true);
    // linear_equation.Solve();            // Compute the solution with pivoting.
    // solution_with_pivoting = linear_equation.Solution_Matrix();
    // // cout << "\nSolution(use pivoting):\n";
    // // solution_with_pivoting.Print_matrix();
    
    // residual_matrix = matrix_A * solution_without_pivoting - vector_b;
    // double norm_of_residual_without_pivoting = residual_matrix.Frobenius_norm();
    // residual_matrix = matrix_A * solution_with_pivoting - vector_b;
    // double norm_of_residual_with_pivoting = residual_matrix.Frobenius_norm();
    
    // cout <<
    // "Norm of residual(without pivoting): "  << norm_of_residual_without_pivoting    <<  ",\n" <<
    // "Norm of residual(with pivoting): "     << norm_of_residual_with_pivoting        << '\n';
    
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
    

    
    return 0;
}

