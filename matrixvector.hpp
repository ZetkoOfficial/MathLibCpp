#ifndef MATRIX_VECTOR_LIB
#define MATRIX_VECTOR_LIB

#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <string>
#include <array>
#include <vector>

namespace matrixvector {

template<int dimensions>
struct vector_f {
    std::array<double, dimensions> data;
    
    vector_f() { data.fill(0); };
    vector_f(double x, double y, double z = 0) {
        if(dimensions == 2) { data[0] = x; data[1] = y; }
        if(dimensions >= 3) { data[0] = x; data[1] = y; data[2] = z; }
    }
    vector_f(const std::vector<double>& data) { for(int i = 0; i < dimensions; i++) this->data[i] = data[i]; }

    constexpr int size() const { return dimensions; }

    double& operator[](int i) { return data[i]; }
    double operator[](int i) const { return data[i]; }
};

/*Adds vector_fs entry-wise.*/
template <int dimensions>
vector_f<dimensions> operator+ (const vector_f<dimensions>& v1, const vector_f<dimensions>& v2) {
    vector_f<dimensions> res; for(int i = 0; i < res.size(); i++) res[i] = v1[i] + v2[i]; return res;
}
/*Subtracts vector_fs entry-wise.*/
template <int dimensions>
vector_f<dimensions> operator- (const vector_f<dimensions>& v1, const vector_f<dimensions>& v2){
    vector_f<dimensions> res; for(int i = 0; i < res.size(); i++) res[i] = v1[i] - v2[i]; return res;
}
/*Returns vector_f with every entry scaled by a constant.*/
template <int dimensions>
vector_f<dimensions> operator* (const double& d, const vector_f<dimensions>& vec) {
    vector_f<dimensions> res; for(int i = 0; i < res.size(); i++) res[i] = vec[i] * d; return res;
}
/*Returns vector_f with every entry scaled by a constant.*/
template <int dimensions>
vector_f<dimensions> operator* (const vector_f<dimensions>& vec, const double& d) {
    vector_f<dimensions> res; for(int i = 0; i < res.size(); i++) res[i] = vec[i] * d; return res;
}
/*Returns vector_f with every entry divided by a constant.*/
template <int dimensions>
vector_f<dimensions> operator/ (const vector_f<dimensions>& vec, const double& d) {
    vector_f<dimensions> res; for(int i = 0; i < res.size(); i++) res[i] = vec[i] / d; return res;
}

/*
Returns the dot product of vector_fs.
*/
template <int dimensions>
double operator* (const vector_f<dimensions>& v1, const vector_f<dimensions>& v2) {
    double res = 0;
    for(int i = 0; i < v1.size(); i++) res += v1[i] * v2[i];

    return res;
}

/*
Writes the vector_f to a stream. 
The vector_f entries are of fixed length 5.
*/
template <int dimensions>
std::ostream& operator<<(std::ostream& stream, const vector_f<dimensions>& vec) {
    stream << "Vector " << vec.size() << "x1 \n";
    
    stream << "| ";
    for(int i = 0; i < vec.size(); i++) stream << std::fixed << std::setprecision(5) << std::setw(10) << vec[i] << " ";
    stream << " |\n";

    return stream;
}

/*Compares if two vector_fs are equal*/
template <int dimensions>
bool operator==(const vector_f<dimensions>& v1, const vector_f<dimensions>& v2) {
    return v1.data == v2.data;
}

/*
Applies a function which determines an entry i of the result vector_f, from entry i. 
*/
template <int dimensions>
static vector_f<dimensions> apply_function(const vector_f<dimensions>& vec, std::function<double(double)> f){
    vector_f<dimensions> res;
    for(int i = 0; i < res.size(); i++) res[i] = f(vec[i]);

    return res;
}
/*
Applies a function which determines an entry i of the result vector_f, from entries i of two vector_fs. 
*/
template <int dimensions>
vector_f<dimensions> apply_function(const vector_f<dimensions>& v1, const vector_f<dimensions>& v2, const std::function<double(double, double)> f){
    vector_f<dimensions> res;
    for(int i = 0; i < res.size(); i++) res[i] = f(v1[i], v2[i]);

    return res;
}

/*Returns length of the vector_f*/
template <int dimensions>
double abs(const vector_f<dimensions>& vec){ return sqrt(vec * vec); }
    
/*Creates a unit vector_f with the same direction as vec*/
template <int dimensions>
vector_f<dimensions> normalize(const vector_f<dimensions>& vec){ 
    if(vec[0] == 0 && vec[1] == 0) return {0,0};
    double len = abs(vec);
    return vec / len;
}

typedef vector_f<2> vector;
typedef vector_f<2> vector2;
typedef vector_f<3> vector3;

template <int SIZE>
struct matrix {
    std::array<matrixvector::vector_f<SIZE>, SIZE> data;

    matrix() { data.fill(matrixvector::vector_f<SIZE>()); }
    matrix(std::vector<std::vector<double>> data) {
        for(int j = 0; j < SIZE; j++){
            for(int i = 0; i < SIZE; i++) this->data[j][i] = data[j][i];
        }
    }
    int size() const { return SIZE; } 

    /*Returns reference to the row vector*/
    matrixvector::vector_f<SIZE>& operator[](int j) { return data[j]; }
    /*Returns copy of the row vector*/
    matrixvector::vector_f<SIZE> operator[](int j) const { return data[j]; }

    /*Creates the column vector */
    matrixvector::vector_f<SIZE> column(int i) const {
        matrixvector::vector_f<SIZE> column;
        for(int j = 0; j < SIZE; j++) column[j] = data[j][i];

        return column;
    }

    matrix<SIZE-1> create_minor(int j, int i) const {
        matrix<SIZE-1> result;
        int jc = 0, ic = 0;
        for(int jt = 0; jt < SIZE; jt++) {
            if(jt == j) continue;
            for(int it = 0; it < SIZE; it++){
                if(it == i) continue;
                result[jc][ic] = data[jt][it];
                ic++;
            }
            jc++;
            ic = 0;
        }

        return result;
    }

    /*Creates the transpose of this matrix*/
    matrix<SIZE> T() const {
        matrix<SIZE> result;
        for(int j = 0; j < SIZE; j++) result[j] = column(j);
        return result;
    }

    /*Creates an identity matrix*/
    static matrix<SIZE> create_identity() {
        matrix<SIZE> res;
        for(int i = 0; i < SIZE; i++) res[i][i] = 1;
        
        return res;
    }
};

/*Adds matrices entry-wise.*/
template <int SIZE>
matrix<SIZE> operator+ (const matrix<SIZE>& m1, const matrix<SIZE>& m2){
    matrix<SIZE> res;
    for(int j = 0; j < SIZE; j++) res[j] = m1[j] + m2[j];

    return res;
}

/*Subtracts matrices entry-wise.*/
template <int SIZE>
matrix<SIZE> operator- (const matrix<SIZE>& m1, const matrix<SIZE>& m2){
    matrix<SIZE> res;
    for(int j = 0; j < SIZE; j++) res[j] = m1[j] - m2[j];

    return res;
}

/*Returns matrix with every entry scaled by a constant*/
template <int SIZE>
matrix<SIZE> operator* (const matrix<SIZE>& m1, double d){
    matrix<SIZE> res;
    for(int j = 0; j < SIZE; j++) res[j] = m1[j] * d;

    return res;
}

/*Returns matrix with every entry scaled by a constant*/
template <int SIZE>
matrix<SIZE> operator* (double d, const matrix<SIZE>& m2){ return m2 * d; }

/*Returns matrix with every entry divided by a constant*/
template <int SIZE>
matrix<SIZE> operator/ (const matrix<SIZE>& m2, double d){ return m2 * (1.0/d); }

/*Compares matrices entry wise*/
template <int SIZE1, int SIZE2>
bool operator== (const matrix<SIZE1>& m1, const matrix<SIZE2>& m2){
    if(SIZE1 != SIZE2) return false;
    for(int j = 0; j < SIZE1; j++){
        if(!(m1[j] == m2[j])) return false;
    }

    return true;
}

/*
Returns matrix which is the result of matrix multiplication.
*/ 
template <int SIZE>
matrix<SIZE> operator* (const matrix<SIZE>& m1, const matrix<SIZE>& m2) {
    matrix<SIZE> res;
    
    for(int i = 0; i < SIZE; i++){
        for(int j = 0; j < SIZE; j++) {
            double cell_value = 0;
            for(int r = 0; r < SIZE; r++) cell_value += m1[i][r] * m2[r][j];
            res[i][j] = cell_value; 
        }
    }

    return res;
}

/*
Writes the matrix to a stream. 
The matrix entries are of fixed length 5.
*/
template <int SIZE>
std::ostream& operator<<(std::ostream& stream, const matrix<SIZE>& mat) {
    stream << "Matrix " << SIZE << "x" << SIZE << "\n";
    for(int j = 0; j < SIZE; j++){
        stream << "| ";
        for(int i = 0; i < SIZE; i++) stream << std::fixed << std::setprecision(5) << std::setw(10) << mat[j][i] << " ";
        stream << " |\n";
    }
    return stream;
}

/*
Returns vector which is the result of matrix-vector multiplication.
*/
template <int SIZE>
matrixvector::vector_f<SIZE> operator* (const matrix<SIZE>& mat, const matrixvector::vector_f<SIZE>& vec) {
    matrixvector::vector_f<SIZE> res;
    for(int j = 0; j < SIZE; j++) res = res + vec[j] * mat.column(j);

    return res;
}


/*
Uses binary exponentiation to quickly calculate power of square matrix.
*/
template <int SIZE>
matrix<SIZE> pow(const matrix<SIZE>& mat, int power) {
    if(power == 0) return matrix<SIZE>::create_identity();

    if(power % 2 == 0) return pow(mat, power/2) * pow(mat, power/2);
    return mat * pow(mat, power - 1);
}


double abs(const matrix<1>& mat){
    return mat[0][0];
}

template <int SIZE>
double abs(const matrix<SIZE>& mat){
    double res = 0;
    for(int i = 0; i < SIZE; i++){
        res += mat[0][i] * (i % 2 == 0 ? 1 : -1 ) * abs(mat.create_minor(0, i)); 
    }
    return res;
}

template <int SIZE>
matrix<SIZE> invert(const matrix<SIZE>& mat) {
    matrix<SIZE> res;
    for(int j = 0; j < SIZE; j++){
        for(int i = 0; i < SIZE; i++){
            res[j][i] = ((i+j) % 2 == 0 ? 1 : -1) * abs(mat.create_minor(j, i));
        }
    }
    return res.T() / abs(mat);
}

template <int SIZE>
vector_f<SIZE> solve_gauss(matrix<SIZE> mat, vector_f<SIZE> b) {
    for(int i = 0; i< SIZE-1; i++) {
        b[i] = b[i] / mat[i][i]; mat[i] = mat[i] / mat[i][i]; 
        for(int j = i+1; j < SIZE; j++){
            b[j] = b[j] - mat[j][i] * b[i]; mat[j] = mat[j] - mat[j][i] * mat[i];;
        }
    }
    for(int i = SIZE-1; i >= 0; i--) {
        b[i] = b[i] / mat[i][i]; mat[i] = mat[i] / mat[i][i]; 
        for(int j = i-1; j >= 0; j--){
            b[j] = b[j] - mat[j][i] * b[i]; mat[j] = mat[j] - mat[j][i] * mat[i]; 
        }
    }
    
    return b;
}

}

typedef std::vector<std::vector<double>> m_def;
#endif