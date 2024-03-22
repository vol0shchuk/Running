// Running.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class MathUtils {
public:
    static double sine(double x) {
        return sin(x);
    }

    static double cosine(double x) {
        return cos(x);
    }

    static double tangent(double x) {
        return tan(x);
    }

    static double cotangent(double x) {
        return 1.0 / tan(x);
    }

    static double power(double base, double exponent) {
        return pow(base, exponent);
    }
};

class Matrix {
private:
    vector<vector<double>> data;
    size_t rows;
    size_t cols;

public:
    Matrix(size_t rows, size_t cols) : rows(rows), cols(cols) {
        data.resize(rows, vector<double>(cols, 0.0));
    }

    void setElement(size_t row, size_t col, double value) {
        if (row >= rows || col >= cols) {
            throw out_of_range("Index out of range");
        }
        data[row][col] = value;
    }

    double getElement(size_t row, size_t col) const {
        if (row >= rows || col >= cols) {
            throw out_of_range("Index out of range");
        }
        return data[row][col];
    }

    size_t getRows() const {
        return rows;
    }

    size_t getCols() const {
        return cols;
    }

    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            throw invalid_argument("Matrix dimensions mismatch");
        }
        Matrix result(rows, cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result.setElement(i, j, data[i][j] + other.getElement(i, j));
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw invalid_argument("Matrix dimensions mismatch");
        }
        Matrix result(rows, other.cols);
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < other.cols; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < cols; ++k) {
                    sum += data[i][k] * other.getElement(k, j);
                }
                result.setElement(i, j, sum);
            }
        }
        return result;
    }

    friend ostream& operator<<(ostream& os, const Matrix& matrix) {
        for (size_t i = 0; i < matrix.rows; ++i) {
            for (size_t j = 0; j < matrix.cols; ++j) {
                os << matrix.getElement(i, j) << " ";
            }
            os << endl;
        }
        return os;
    }
};

int main() {
    cout << "sin(1.57) = " << MathUtils::sine(1.57) << endl;
    cout << "cos(0) = " << MathUtils::cosine(0) << endl;
    cout << "tan(0.785) = " << MathUtils::tangent(0.785) << endl;
    cout << "cot(1.047) = " << MathUtils::cotangent(1.047) << endl;
    cout << "2^3 = " << MathUtils::power(2, 3) << endl;

    Matrix A(2, 2);
    A.setElement(0, 0, 1);
    A.setElement(0, 1, 2);
    A.setElement(1, 0, 3);
    A.setElement(1, 1, 4);

    Matrix B(2, 2);
    B.setElement(0, 0, 5);
    B.setElement(0, 1, 6);
    B.setElement(1, 0, 7);
    B.setElement(1, 1, 8);

    cout << "Matrix A:\n" << A << endl;
    cout << "Matrix B:\n" << B << endl;

    Matrix C = A + B;
    cout << "Matrix A + B:\n" << C << endl;

    Matrix D = A * B;
    cout << "Matrix A * B:\n" << D << endl;

    return 0;
}



// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
