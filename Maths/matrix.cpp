////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: matrix.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "maths.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix::Matrix(Field* field, vector<vector<Element>> elements) {
    #ifdef DEBUG
        for (int i = 1; i < elements.size(); i ++) {
            if (elements[i].size() != elements[0].size()) {
                throw "Maths\\Matrix\\Matrix(Field*, vector<vector<Element>>)\\elements\\size";
            }
        }
        for (int i = 0; i < elements.size(); i ++) {
            for (int j = 0; j < elements[0].size(); j ++) {
                if (*elements[i][j].field != *field) {
                    throw "Maths\\Matrix\\Matrix(Field*, vector<vector<Element>>)\\elements\\field";
                }
            }
        }
    #endif
    this->field = field;
    this->elements = elements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool Matrix::operator==(const Matrix& matrix) const {
    if (this->field == matrix.field && this->elements == matrix.elements) return true;
    return *this->field == *matrix.field && this->elements == matrix.elements;
}

bool Matrix::operator!=(const Matrix& matrix) const {
    return !(*this == matrix);
}

ostream& operator<<(ostream& output, const Matrix& matrix) {
    output << "{";
    for (int i = 0; i < matrix.row(); i ++) {
        output << "{";
        for (int j = 0; j < matrix.column(); j ++) {
            output << matrix.elements[i][j];
            if (j != matrix.column() - 1) {
                output << ", ";
            }
        }
        output << "}";
        if (i != matrix.row() - 1) {
            output << ", ";
        }
    }
    output << "}";
    return output;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int Matrix::row() const {
    return this->elements.size();
}

int Matrix::column() const {
    return this->elements[0].size();
}

Matrix Matrix::resize(int row, int column) const {
    #ifdef DEBUG
        if (row <= 0) {
            throw "Maths\\Matrix\\resize(int, int)\\row";
        } else if (column <= 0) {
            throw "Maths\\Matrix\\resize(int, int)\\column";
        }
    #endif
    vector<vector<Element>> new_elements(row, vector<Element>(column, this->field->zero_element()));
    for (int i = 0; i < min(row, this->row()); i ++) {
        for (int j = 0; j < min(column, this->column()); j ++) {
            new_elements[i][j] = this->elements[i][j];
        }
    }
    return Matrix(this->field, new_elements);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix Matrix::operator~() const {
    #ifdef DEBUG
        if (this->row() != this->column() || this->determinant() == this->field->zero_element()) {
            throw "Maths\\Matrix\\operator~()";
        }
    #endif
    vector<vector<Element>> new_elements(this->row(), vector<Element>(this->column(), this->field->zero_element()));
    for (int i = 0; i < this->row(); i ++) {
        new_elements[i][i] = this->field->unit_element();
    }
    vector<vector<Element>> temp_0_elements = this->elements;
    for (int i = 0; i < this->row(); i ++) {
        if (temp_0_elements[i][i] == this->field->zero_element()) {
            for (int j = i + 1; j < this->row(); j ++) {
                if (temp_0_elements[j][i] != this->field->zero_element()) {
                    swap(temp_0_elements[i], temp_0_elements[j]);
                    swap(new_elements[i], new_elements[j]);
                    break;
                }
                if (j == this->row() - 1) throw "ERROR 8712";
            }
        }
        Element temp_1_diagonal = temp_0_elements[i][i];
        for (int j = 0; j < this->row(); j ++) {
            temp_0_elements[i][j] = temp_0_elements[i][j] / temp_1_diagonal;
            new_elements[i][j] = new_elements[i][j] / temp_1_diagonal;
        }
        for (int j = 0; j < this->row(); j ++) {
            if (j != i) {
                Element temp_2_factor = temp_0_elements[j][i];
                for (int k = 0; k < this->row(); k ++) {
                    temp_0_elements[j][k] = temp_0_elements[j][k] - (temp_2_factor * temp_0_elements[i][k]);
                    new_elements[j][k] = new_elements[j][k] - (temp_2_factor * new_elements[i][k]);
                }
            }
        }
    }
    return Matrix(this->field, new_elements);
}

Matrix Matrix::operator+(const Matrix& matrix) const {
    #ifdef DEBUG
        if (*matrix.field != *this->field) {
            throw "Maths\\Matrix\\operator+(const Matrix&)\\matrix\\field";
        } else if (matrix.row() != this->row()) {
            throw "Maths\\Matrix\\operator+(const Matrix&)\\matrix\\row";
        } else if (matrix.column() != this->column()) {
            throw "Maths\\Matrix\\operator+(const Matrix&)\\matrix\\column";
        }
    #endif
    int new_row = this->row();
    int new_column = this->column();
    vector<vector<Element>> new_elements(new_row, vector<Element>(new_column, this->field->zero_element()));
    for (int i = 0; i < new_row; i ++) {
        for (int j = 0; j < new_column; j ++) {
            new_elements[i][j] = this->elements[i][j] + matrix.elements[i][j];
        }
    }
    return Matrix(this->field, new_elements);
}

Matrix Matrix::operator*(const Matrix& matrix) const {
    #ifdef DEBUG
        if (*matrix.field != *this->field) {
            throw "Maths\\Matrix\\operator*(const Matrix&)\\matrix\\field";
        } else if (matrix.row() != this->column()) {
            throw "Maths\\Matrix\\operator*(const Matrix&)\\matrix\\row";
        }
    #endif
    int new_row = this->row();
    int new_column = matrix.column();
    vector<vector<Element>> new_elements(new_row, vector<Element>(new_column, this->field->zero_element()));
    for (int i = 0; i < new_row; i ++) {
        for (int j = 0; j < new_column; j ++) {
            for (int k = 0; k < this->column(); k ++) {
                new_elements[i][j] = new_elements[i][j] + (this->elements[i][k] * matrix.elements[k][j]);
            }
        }
    }
    return Matrix(this->field, new_elements);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix Matrix::operator-() const {
    int new_row = this->row();
    int new_column = this->column();
    vector<vector<Element>> new_elements(new_row, vector<Element>(new_column, this->field->zero_element()));
    for (int i = 0; i < new_row; i ++) {
        for (int j = 0; j < new_column; j ++) {
            new_elements[i][j] = -this->elements[i][j];
        }
    }
    return Matrix(this->field, new_elements);
}

Matrix Matrix::operator-(const Matrix& matrix) const {
    #ifdef DEBUG
        if (*matrix.field != *this->field) {
            throw "Maths\\Matrix\\operator-(const Matrix&)\\matrix\\field";
        } else if (matrix.row() != this->row()) {
            throw "Maths\\Matrix\\operator-(const Matrix&)\\matrix\\row";
        } else if (matrix.column() != this->column()) {
            throw "Maths\\Matrix\\operator-(const Matrix&)\\matrix\\column";
        }
    #endif
    return (*this) + (-matrix);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Matrix Matrix::operator*(const Element& scalar) const {
    #ifdef DEBUG
        if (scalar.field != this->field) {
            throw "Maths\\Matrix\\operator*(const Element&)\\scalar\\field";
        }
    #endif
    int new_row = this->row();
    int new_column = this->column();
    vector<vector<Element>> new_elements(new_row, vector<Element>(new_column, this->field->zero_element()));
    for (int i = 0; i < new_row; i ++) {
        for (int j = 0; j < new_column; j ++) {
            new_elements[i][j] = this->elements[i][j] * scalar;
        }
    }
    return Matrix(this->field, new_elements);
}

Matrix Matrix::transpose() const {
    int new_row = this->column();
    int new_column = this->row();
    vector<vector<Element>> new_elements(new_row, vector<Element>(new_column, this->field->zero_element()));
    for (int i = 0; i < new_row; i ++) {
        for (int j = 0; j < new_column; j ++) {
            new_elements[i][j] = this->elements[j][i];
        }
    }
    return Matrix(this->field, new_elements);
}

Element Matrix::determinant() const {
    #ifdef DEBUG
        if (this->row() != this->column()) {
            throw "Maths\\Matrix\\determinant()";
        }
    #endif
    vector<int> temp_0_permutation(this->row());
    for (int i = 0; i < this->row(); i ++) {
        temp_0_permutation[i] = i;
    }
    Element new_element = this->field->zero_element();
    bool temp_1_sign = false;
    do {
        Element temp_2_element = this->field->unit_element();
        for (int i = 0; i < this->row(); i ++) {
            temp_2_element = temp_2_element * this->elements[i][temp_0_permutation[i]];
        }
        temp_1_sign = !temp_1_sign;
        if (temp_1_sign) {
            new_element = new_element + temp_2_element;
        } else {
            new_element = new_element - temp_2_element;
        }
    } while (next_permutation(temp_0_permutation.begin(), temp_0_permutation.end()));
    return new_element;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////