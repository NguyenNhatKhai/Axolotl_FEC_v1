////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: polynomial.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "maths.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial::Polynomial() {
    this->field = &fields::default_field;
    this->coefficients = vector<Element>();
}

Polynomial::Polynomial(Field* field, vector<Element> coefficients) {
    #ifdef DEBUG
        if (coefficients.size() == 0) {
            throw "Maths\\Polynomial\\Polynomial(Field*, vector<Element>)\\coefficients\\size";
        }
        for (int i = 0; i < coefficients.size(); i ++) {
            if (*coefficients[i].field != *field) {
                throw "Maths\\Polynomial\\Polynomial(Field*, vector<Element>)\\coefficients\\field";
            }
        }
    #endif
    this->field = field;
    this->coefficients = coefficients;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool Polynomial::operator==(const Polynomial& polynomial) const {
    if (this == &polynomials::default_polynomial && &polynomial != &polynomials::default_polynomial || this != &polynomials::default_polynomial && &polynomial == &polynomials::default_polynomial) return false;
    if (this->field == polynomial.field && this->align().coefficients == polynomial.align().coefficients) return true;
    return *this->field == *polynomial.field && this->align().coefficients == polynomial.align().coefficients;
}

bool Polynomial::operator!=(const Polynomial& polynomial) const {
    return !(*this == polynomial);
}

ostream& operator<<(ostream& output, const Polynomial& polynomial) {
    output << "{";
    for (int i = 0; i < polynomial.coefficients.size(); i ++) {
        output << polynomial.coefficients[i];
        if (i != polynomial.coefficients.size() - 1) {
            output << ", ";
        }
    }
    output << "}";
    return output;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int Polynomial::degree() const {
    for (int i = 0; i < this->coefficients.size(); i ++) {
        if (this->coefficients[this->coefficients.size() - i - 1] != this->field->zero_element()) {
            return this->coefficients.size() - i - 1;
        }
    }
    return 0;
}

Polynomial Polynomial::redegree(int degree) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\redegree(int)";
        } else if (degree < 0) {
            throw "Maths\\Polynomial\\redegree(int)\\degree";
        }
    #endif
    vector<Element> new_coefficients(degree + 1, this->field->zero_element());
    for (int i = 0; i <= min(degree, this->degree()); i ++) {
        new_coefficients[i] = this->coefficients[i];
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::align() const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\align()";
        }
    #endif
    return this->redegree(this->degree());
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial Polynomial::operator+(const Polynomial& polynomial) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator+(const Polynomial&)";
        } else if (polynomial == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator+(const Polynomial&)\\polynomial";
        } else if (*polynomial.field != *this->field) {
            throw "Maths\\Polynomial\\operator+(const Polynomial&)\\polynomial\\field";
        }
    #endif
    int new_degree = max(this->degree(), polynomial.degree());
    Polynomial temp_0_polynomial = this->redegree(new_degree);
    Polynomial temp_1_polynomial = polynomial.redegree(new_degree);
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[i] = temp_0_polynomial.coefficients[i] + temp_1_polynomial.coefficients[i];
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::operator*(const Polynomial& polynomial) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator*(const Polynomial&)";
        } else if (polynomial == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator*(const Polynomial&)\\polynomial";
        } else if (polynomial.field != this->field) {
            throw "Maths\\Polynomial\\operator*(const Polynomial&)\\polynomial\\field";
        }
    #endif
    int new_degree = this->degree() + polynomial.degree();
    Polynomial temp_0_polynomial = this->align();
    Polynomial temp_1_polynomial = polynomial.align();
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= temp_0_polynomial.degree(); i ++) {
        for (int j = 0; j <= temp_1_polynomial.degree(); j ++) {
            new_coefficients[i + j] = new_coefficients[i + j] + (temp_0_polynomial.coefficients[i] * temp_1_polynomial.coefficients[j]);
        }
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::operator/(const Polynomial& polynomial) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator/(const Polynomial&)";
        } else if (polynomial == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator/(const Polynomial&)\\polynomial";
        } else if (polynomial.field != this->field) {
            throw "Maths\\Polynomial\\operator/(const Polynomial&)\\polynomial\\field";
        }
    #endif
    int new_degree = this->degree() - polynomial.degree();
    if (new_degree < 0) return Polynomial(this->field, {this->field->zero_element()});
    Polynomial temp_0_polynomial = this->align();
    Polynomial temp_1_polynomial = polynomial.align();
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[new_degree - i] = temp_0_polynomial.coefficients[this->degree() - i] / temp_1_polynomial.coefficients[temp_1_polynomial.degree()];
        for (int j = 0; j < polynomial.degree(); j ++) {
            temp_0_polynomial.coefficients[new_degree - i + j] = temp_0_polynomial.coefficients[new_degree - i + j] - (new_coefficients[new_degree - i] * temp_1_polynomial.coefficients[j]);
        }
    }
    return Polynomial(this->field, new_coefficients);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial Polynomial::operator-() const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator-()";
        }
    #endif
    int new_degree = this->degree();
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[i] = -this->coefficients[i];
    }
    return Polynomial(this->field, new_coefficients);
}

Polynomial Polynomial::operator-(const Polynomial& polynomial) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator-(const Polynomial&)";
        } else if (polynomial == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator-(const Polynomial&)\\polynomial";
        } else if (polynomial.field != this->field) {
            throw "Maths\\Polynomial\\operator-(const Polynomial&)\\polynomial\\field";
        }
    #endif
    return (*this) + (-polynomial);
}

Polynomial Polynomial::operator%(const Polynomial& polynomial) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator%(const Polynomial&)";
        } else if (polynomial == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator%(const Polynomial&)\\polynomial";
        } else if (polynomial.field != this->field) {
            throw "Maths\\Polynomial\\operator%(const Polynomial&)\\polynomial\\field";
        }
    #endif
    int new_degree = polynomial.degree() - 1;
    if (new_degree < 0) return Polynomial(this->field, {this->field->zero_element()});
    return (*this - (polynomial * (*this / polynomial))).redegree(new_degree);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial Polynomial::operator*(const Element& scalar) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\operator*(const Element&)";
        } else if (scalar.field != this->field) {
            throw "Maths\\Polynomial\\operator*(const Element&)\\scalar\\field";
        }
    #endif
    int new_degree = this->degree();
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[i] = this->coefficients[i] * scalar;
    }
    return Polynomial(this->field, new_coefficients);
}

Element Polynomial::evaluate(const Element& argument) const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\evaluate()";
        } else if (argument.field != this->field) {
            throw "Maths\\Polynomial\\evaluate()\\argument\\field";
        }
    #endif
    Element new_element = this->field->zero_element();
    for (int i = 0; i <= this->degree(); i ++) {
        new_element = new_element * argument + this->coefficients[this->degree() - i];
    }
    return new_element;
}

Polynomial Polynomial::derivative() const {
    #ifdef DEBUG
        if (*this == polynomials::default_polynomial) {
            throw "Maths\\Polynomial\\derivative()";
        }
    #endif
    int new_degree = this->degree() - 1;
    if (new_degree < 0) return Polynomial(this->field, {this->field->zero_element()});
    vector<Element> new_coefficients(new_degree + 1, this->field->zero_element());
    for (int i = 0; i <= new_degree; i ++) {
        new_coefficients[i] = (this->coefficients[i + 1] * (i + 1));
    }
    return Polynomial(this->field, new_coefficients);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////