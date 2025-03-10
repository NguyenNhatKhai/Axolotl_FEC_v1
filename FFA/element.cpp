////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: element.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "ffa.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Element::Element(bool value) {
    this->field = &fields::default_field;
    this->value = value;
    this->values = vector<Element>();
}

Element::Element(Field* field, vector<Element> values) {
    #ifdef DEBUG
        if (values.size() != field->primitive_polynomial->degree()) {
            throw "FFA\\Element\\Element(Field*, vector<Element>)\\values\\size";
        }
    #endif
    this->field = field;
    this->value = false;
    this->values = values;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool Element::operator==(const Element& element) const {
    if (this->field == element.field && this->value == element.value && this->values == element.values) return true;
    return *this->field == *element.field && this->value == element.value && this->values == element.values;
}

bool Element::operator!=(const Element& element) const {
    return !(*this == element);
}

ostream& operator<<(ostream& output, const Element& element) {
    if (*element.field == fields::default_field) {
        output << element.value;
    } else {
        output << "{";
        for (int i = 0; i < element.values.size(); i ++) {
            output << element.values[i];
            if (i != element.values.size() - 1) {
                output << ", ";
            }
        }
        output << "}";
    }
    return output;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int Element::size() const {
    return this->field->primitive_polynomial->degree();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Element Element::operator+(const Element& element) const {
    #ifdef DEBUG
        if (*element.field != *this->field) {
            throw "FFA\\Element\\operator+(const Element&)\\element\\field";
        }
    #endif
    if (*this->field == fields::default_field) {
        bool new_value = this->value ^ element.value;
        return Element(new_value);
    } else {
        Polynomial temp_0_polynomial(this->field->primitive_polynomial->field, this->values);
        Polynomial temp_1_polynomial(element.field->primitive_polynomial->field, element.values);
        vector<Element> new_values = (temp_0_polynomial + temp_1_polynomial).redegree(this->size() - 1).coefficients;
        return Element(this->field, new_values);
    }
}

Element Element::operator*(const Element& element) const {
    #ifdef DEBUG
        if (*element.field != *this->field) {
            throw "FFA\\Element\\operator*(const Element&)\\element\\field";
        }
    #endif
    if (*this->field == fields::default_field) {
        bool new_value = this->value & element.value;
        return Element(new_value);
    } else {
        Polynomial temp_0_polynomial(this->field->primitive_polynomial->field, this->values);
        Polynomial temp_1_polynomial(element.field->primitive_polynomial->field, element.values);
        vector<Element> new_values = (temp_0_polynomial * temp_1_polynomial % *this->field->primitive_polynomial).redegree(this->size() - 1).coefficients;
        return Element(this->field, new_values);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Element Element::operator-() const {
    for (int i = 0; i < this->field->general_elements.size(); i ++) {
        if (*this + this->field->general_elements[i] == this->field->zero_element()) return this->field->general_elements[i];
    }
    throw "FFA\\Element\\operator-()";
}

Element Element::operator~() const {
    for (int i = 0; i < this->field->general_elements.size(); i ++) {
        if (*this * this->field->general_elements[i] == this->field->unit_element()) return this->field->general_elements[i];
    }
    throw "FFA\\Element\\operator~()";
}

Element Element::operator*(int times) const {
    #ifdef DEBUG
        if (times < 0) {
            throw "FFA\\Element\\operator*(int)\\times";
        }
    #endif
    Element new_element = this->field->zero_element();
    for (int i = 0; i < times; i ++) {
        new_element = new_element + *this;
    }
    return new_element;
}

Element Element::operator^(int power) const {
    #ifdef DEBUG
        if (power < 0) {
            throw "FFA\\Element\\operator^(int)\\power";
        }
    #endif
    Element new_element = this->field->unit_element();
    for (int i = 0; i < power; i ++) {
        new_element = new_element * *this;
    }
    return new_element;
}

Element Element::operator-(const Element& element) const {
    #ifdef DEBUG
        if (*element.field != *this->field) {
            throw "FFA\\Element\\operator-(const Element&)\\element\\field";
        }
    #endif
    return (*this) + (-element);
}

Element Element::operator/(const Element& element) const {
    #ifdef DEBUG
        if (*element.field != *this->field) {
            throw "FFA\\Element\\operator~(const Element&)\\element\\field";
        }
    #endif
    return (*this) * (~element);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////