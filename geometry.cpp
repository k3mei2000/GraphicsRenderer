#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include "geometry.h"

template <> template <> Vec3<int>::Vec3(const Vec3<float> &v) : x(int(v.x+.5)), y(int(v.y+.5)), z(int(v.z+.5)) {

}

template <> template <> Vec3<float>::Vec3(const Vec3<int> &v) : x(v.x), y(v.y), z(v.z) {

}

Matrix::Matrix(int r, int c) : m(std::vector<std::vector<float> >(r, std::vector<float>(c, 0.f))), rows(r), cols(c) {

}

// int Matrix::nrows() {
//     return rows;
// }

// int Matrix::ncols() {
//     return cols;
// }

Matrix Matrix::identity(int dim) {
    Matrix I(dim, dim);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            I[i][j] = (i==j ? 1.f : 0.f);
        }
    }
    return I;
}

std::vector<float>& Matrix::operator[](const int i) {
    assert(i>=0 && i<rows);
    return m[i];
}

Matrix Matrix::operator*(const Matrix& a) {
    assert(cols == a.rows);
    Matrix result(rows, a.cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.m[i][j] = 0.f;
            for (int k = 0; k < cols; k++) {
                result.m[i][j] += m[i][k]*a.m[k][j];
            }
        }
    }

    return result;
}

Matrix Matrix::transpose() {
    Matrix result(cols,rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[j][i] = m[i][j];
        }
    }
    return result;
}

Matrix Matrix::inverse() {
    assert(rows==cols);

    Matrix result(rows, cols*2);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = m[i][j];
        }
    }
    for (int i = 0; i < rows; i++) {
        result[i][i+cols] = 1;
    }

    for (int i = 0; i < rows-1; i++) {
        for (int j = result.cols-1; j >= 0; j--) {
            result[i][j] /= result[i][i];
        }
        for (int k = i+1; k < rows; k++ ) {
            float coeff = result[k][i];
            for (int j=0; j<result.cols; j++) {
                result[k][j] -= result[i][j]*coeff;
            }
        }
    }

    for (int j = result.cols-1; j >= rows - 1; j--) {
        result[rows-1][j] /= result[rows-1][rows-1];
    }

    Matrix truncate(rows,cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            truncate[i][j] = result[i][j+cols];
        }
    }
    return truncate;
}

std::ostream& operator<<(std::ostream& s, Matrix& m) {
    for (int i = 0; i < m.nrows(); i++) {
        for (int j = 0; j < m.ncols(); j++) {
            s << m[i][j];
            if (j < m.ncols() - 1) s << "\t";
        }
        s << "\n";
    }
    return s;
}