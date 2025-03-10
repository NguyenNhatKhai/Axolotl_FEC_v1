////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: ffa.h
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _FFA_H_
#define _FFA_H_

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class Polynomial;
class Field;
class Element;

////////////////////////////////////////////////////////////////////////////////////////////////////

class Field {
    public:
    Polynomial* primitive_polynomial;
    vector<Element> general_elements;

    public:
    Field();
    Field(Polynomial* primitive_polynomial);
    ~Field() = default;

    public:
    bool operator==(const Field& field) const;
    bool operator!=(const Field& field) const;

    public:
    int size() const;
    Element zero_element() const;
    Element unit_element() const;
    Element primitive_element() const;
};

namespace fields {
    static Field default_field;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

class Element {
    public:
    Field* field;
    bool value;
    vector<Element> values;

    public:
    Element() = delete;
    Element(bool value);
    Element(Field* field, vector<Element> values);
    ~Element() = default;

    public:
    bool operator==(const Element& element) const;
    bool operator!=(const Element& element) const;

    public:
    Element operator-() const;
    Element operator~() const;
    Element operator*(int times) const;
    Element operator^(int power) const;
    Element operator+(const Element& element) const;
    Element operator-(const Element& element) const;
    Element operator*(const Element& element) const;
    Element operator/(const Element& element) const;

    public:
    int size() const;
};

ostream& operator<<(ostream& output, const Element& element);

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../Maths/maths.h"
#include "element.cpp"
#include "field.cpp"

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////