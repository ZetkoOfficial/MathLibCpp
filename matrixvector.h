#include <iostream>
#include <iomanip>
#include <functional>
#include <vector>
#include <string>

using namespace std;

typedef vector<vector<double>> m_def;
namespace matrixvector {

class matrix {
    vector<vector<double>> matrix_arr;
    public:
    int width, height;

    matrix() = default;
    matrix(int height, int width) {
        this->height = height;
        this->width = width;
        this->matrix_arr = vector<vector<double>> (height, vector<double>(width));
    }
    matrix(vector<vector<double>> matrix_arr) {
        this->height = matrix_arr.size();
        this->width = matrix_arr[0].size();
        this->matrix_arr = matrix_arr;
    }

    //Value at row i and column j
    double& operator() (int i, int j) {
        return matrix_arr[i][j];
    }

    double operator() (int i, int j) const {
        return matrix_arr[i][j];
    }

    //Transpose of matrix
    matrix T() {
        matrix res (width, height);
        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++) res(j, i) = matrix_arr[i][j];
        }

        return res;
    }
};
//Applies a function to the matrix
matrix apply_function(const matrix& mat, double function(double)){
    matrix res (mat.height, mat.width);
    for(int i = 0; i < mat.height; i++){
        for(int j = 0; j < mat.width; j++) res(i, j) = function(mat(i, j));
    }

    return res;
}

//Applies a function to the matrix
matrix apply_function(const matrix& m1, const matrix& m2, double function(double, double)){
    matrix res (m1.height, m1.width);
    for(int i = 0; i < m1.height; i++){
        for(int j = 0; j < m1.width; j++) res(i, j) = function(m1(i, j), m2(i, j));
    }

    return res;
}

//Matrix addtion
matrix operator+ (const matrix& m1, const matrix& m2){
    return apply_function(m1, m2, [](double x, double y){return x + y;});
}

//Matrix subtraction
matrix operator- (const matrix& m1, const matrix& m2){
    return apply_function(m1, m2, [](double x, double y){return x - y;});
}

//Matrix multiplication
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

matrix operator* (const double& d, const matrix& mat) {
    matrix res (mat.height, mat.width);
    for(int i = 0; i < mat.height; i++){
        for(int j = 0; j < mat.width; j++) res(i, j) = mat(i, j) * d;
    }

    return res;
}

//Prints matrix to stream
ostream& operator<<(ostream& stream, const matrix& mat){
    for(int i = 0; i < mat.height; i++){
        for(int j = 0; j < mat.width; j++) stream << setw(5) << mat(i, j) << " ";
        cout << "\n";
    }

    return stream;
}

}