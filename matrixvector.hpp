#ifndef MATRIX_VECTOR_LIB
#define MATRIX_VECTOR_LIB

#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <vector>
#include <string>

typedef std::vector<std::vector<double>> m_def;

namespace matrixvector {
typedef std::vector<double> vector;

class matrix {
    std::vector<std::vector<double>> matrix_arr;
    public:
    int width, height;

    matrix() = default;
    matrix(int height, int width) {
        this->height = height;
        this->width = width;
        this->matrix_arr = std::vector<std::vector<double>> (height, std::vector<double>(width));
    }
    matrix(std::vector<std::vector<double>> matrix_arr) {
        this->height = matrix_arr.size();
        this->width = matrix_arr[0].size();
        this->matrix_arr = matrix_arr;
    }

    //Reference of value at row i and column j
    double& operator() (int i, int j) { return matrix_arr[i][j]; }

    //Value at row i and column j
    double operator() (int i, int j) const { return matrix_arr[i][j]; }

    //Returns vector representing row i.
    vector row(int i) const { return matrix_arr[i]; }

    //Returns vector representing column j.
    vector column(int j) const { 
        vector res (height);
        for(int i = 0; i < height; i++) res[i] = matrix_arr[i][j];

        return res; 
    }

    /*
    Returns the transpose of this matrix.

    Time complexity: O(width*height).
    */
    matrix T() {
        matrix res (width, height);
        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++) res(j, i) = matrix_arr[i][j];
        }

        return res;
    }

    static matrix create_identity(int size) {
        matrix res (size, size);
        for(int i = 0; i < size; i++) res(i, i) = 1;

        return res;
    }
};
/*
Applies a function which determines an entry ij of the result matrix, from entry ij. 
*/
matrix apply_function(const matrix& mat, std::function<double(double)> f){
    matrix res (mat.height, mat.width);
    for(int i = 0; i < mat.height; i++){
        for(int j = 0; j < mat.width; j++) res(i, j) = f(mat(i, j));
    }

    return res;
}

/*
Applies a function which determines an entry ij of the result matrix, from entries ij of two matrices. 
*/
matrix apply_function(const matrix& m1, const matrix& m2, std::function<double(double, double)> f){
    matrix res (m1.height, m1.width);
    for(int i = 0; i < m1.height; i++){
        for(int j = 0; j < m1.width; j++) res(i, j) = f(m1(i, j), m2(i, j));
    }

    return res;
}

/*
Applies a function which determines an entry i of the result vector, from entry i. 
*/
vector apply_function(const vector& vec, std::function<double(double)> f){
    vector res (vec.size());
    for(int i = 0; i < vec.size(); i++) res[i] = f(vec[i]);

    return res;
}

/*
Applies a function which determines an entry i of the result vector, from entries i of two vectors. 
*/
vector apply_function(const vector& v1, const vector& v2, std::function<double(double, double)> f){
    vector res (v1.size());
    for(int i = 0; i < v1.size(); i++) res[i] = f(v1[i], v2[i]);

    return res;
}

/*
Returns the length of the vector.
*/
double abs(const vector& vec) {
    double res = 0;
    for(int i = 0; i < vec.size(); i++) res += vec[i] * vec[i];

    return sqrt(res);
}

/*Adds matrices entry-wise.*/
matrix operator+ (const matrix& m1, const matrix& m2){
    return apply_function(m1, m2, [](double x, double y){return x + y;});
}
/*Adds vectors entry-wise.*/
vector operator+ (const vector& v1, const vector& v2){
    return apply_function(v1, v2, [](double x, double y){return x + y;});
}


/*Subtracts matrices entry-wise.*/
matrix operator- (const matrix& m1, const matrix& m2){
    return apply_function(m1, m2, [](double x, double y){return x - y;});
}
/*Subtracts vectors entry-wise.*/
vector operator- (const vector& v1, const vector& v2){
    return apply_function(v1, v2, [](double x, double y){return x - y;});
}


/*Returns matrix with every entry scaled by a constant.*/
matrix operator* (const double& d, const matrix& mat) {
    return apply_function(mat, [d](double i) {return i * d;});
}
/*Returns vector with every entry scaled by a constant.*/
vector operator* (const double& d, const vector& vec) {
    return apply_function(vec, [d](double i) {return i * d;});
}


/*Returns matrix with every entry scaled by a constant.*/
matrix operator* (const matrix& mat, const double& d) {
    return apply_function(mat, [d](double i) {return i * d;});
}
/*Returns vector with every entry scaled by a constant.*/
vector operator* (const vector& vec, const double& d) {
    return apply_function(vec, [d](double i) {return i * d;});
}

/*
Returns matrix which is the result of matrix multiplication.
*/ 
matrix operator* (const matrix& m1, const matrix& m2) {
    matrix res (m1.height, m2.width);
    
    for(int i = 0; i < m1.height; i++){
        for(int j = 0; j < m2.width; j++) {
            double cell_value = 0;
            for(int r = 0; r < m1.width; r++) cell_value += m1(i, r) * m2(r, j);
            res(i, j) = cell_value; 
        }
    }

    return res;
}

/*
Returns vector which is the result of matrix-vector multiplication.
*/
vector operator* (const matrix& mat, const vector& vec) {
    vector res (vec.size());
    for(int j = 0; j < vec.size(); j++) res = res + vec[j] * mat.column(j);

    return res;
}

/*
Returns the dot product of vectors.
*/
double operator* (const vector& v1, const vector& v2) {
    double res = 0;
    for(int i = 0; i < v1.size(); i++) res += v1[i] * v2[i];

    return res;
}

/*
Writes the matrix to a stream. 
The matrix entries are of fixed length 5.
*/
std::ostream& operator<<(std::ostream& stream, const matrix& mat) {
    stream << "Matrix " << mat.height << "x" << mat.width << "\n";

    for(int i = 0; i < mat.height; i++){
        stream << "| ";
        for(int j = 0; j < mat.width; j++) stream << std::fixed << std::setprecision(5) << std::setw(10) << mat(i, j) << " ";
        stream << " |\n";
    }

    return stream;
}
/*
Writes the vector to a stream. 
The vector entries are of fixed length 5.
*/
std::ostream& operator<<(std::ostream& stream, const vector& vec) {
    stream << "Vector " << vec.size() << "x1 \n";
    
    stream << "| ";
    for(int i = 0; i < vec.size(); i++) stream << std::fixed << std::setprecision(5) << std::setw(10) << vec[i] << " ";
    stream << " |\n";

    return stream;
}

/*
Uses binary exponentiation to quickly calculate power of square matrix.
*/
matrix pow(const matrix& mat, int power) {
    if(power == 0) return matrix::create_identity(mat.height);

    if(power % 2 == 0) return pow(mat, power/2) * pow(mat, power/2);
    return mat * pow(mat, power - 1);
}

}

/*
Allows usage of operators without namespace.
*/
using matrixvector::operator*;
using matrixvector::operator+;
using matrixvector::operator-;
using matrixvector::operator<<;

#endif