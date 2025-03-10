////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: field.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "ffa.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Field::Field() {
    this->primitive_polynomial = &polynomials::default_polynomial;
    this->general_elements = {Element(true), Element(false)};
}

Field::Field(Polynomial* primitive_polynomial) {
    this->primitive_polynomial = primitive_polynomial;
    int temp_0_size = pow(primitive_polynomial->field->size(), primitive_polynomial->degree());
    this->general_elements = vector<Element>(temp_0_size, Element(this, vector<Element>(primitive_polynomial->degree(), primitive_polynomial->field->zero_element())));
    vector<Element> temp_1_elements = vector<Element>(temp_0_size, Element(this, vector<Element>(primitive_polynomial->degree(), primitive_polynomial->field->zero_element())));
    for (int i = 0; i < temp_0_size; i ++) {
        int temp_2_index = i;
        for (int j = primitive_polynomial->degree() - 1; j >= 0; j --) {
            temp_1_elements[i].values[j] = primitive_polynomial->field->general_elements[temp_2_index % primitive_polynomial->field->size()];
            temp_2_index = temp_2_index / primitive_polynomial->field->size();
        }
    }
    Element temp_3_element(this, vector<Element>(primitive_polynomial->degree(), primitive_polynomial->field->zero_element()));
    for (int i = 0; i < primitive_polynomial->degree(); i ++) {
        temp_3_element.values[i] = primitive_polynomial->coefficients[i];
    }
    for (int i = 0; i < temp_0_size; i ++) {
        if (temp_1_elements[i] + temp_1_elements[i] == temp_1_elements[i]) {
            this->general_elements[temp_0_size - 1] = temp_1_elements[i];
        }
    }
    for (int i = 0; i < temp_0_size; i ++) {
        if (temp_1_elements[i] * temp_1_elements[i] == temp_1_elements[i] && temp_1_elements[i] != this->zero_element()) {
            this->general_elements[0] = temp_1_elements[i];
        }
    }
    for (int i = 0; i < temp_0_size; i ++) {
        if ((temp_1_elements[i] ^ primitive_polynomial->degree()) == temp_3_element) {
            this->general_elements[1] = temp_1_elements[i];
        }
    }
    for (int i = 2; i < temp_0_size - 1; i ++) {
        this->general_elements[i] = (this->primitive_element() ^ i);
    }
    for (int i = 0; i < temp_0_size; i ++) {
        for (int j = i + 1; j < temp_0_size; j ++) {
            if (this->general_elements[i] == this->general_elements[j]) {
                throw "FFA\\Field\\Field(Polynomial*)";
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool Field::operator==(const Field& field) const {
    if (this == &fields::default_field && &field != &fields::default_field || this != &fields::default_field && &field == &fields::default_field) return false;
    if (this->primitive_polynomial == field.primitive_polynomial) return true;
    return *this->primitive_polynomial == *field.primitive_polynomial;
}

bool Field::operator!=(const Field& field) const {
    return !(*this == field);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int Field::size() const {
    return this->general_elements.size();
}

Element Field::zero_element() const {
    return this->general_elements[this->size() - 1];
}

Element Field::unit_element() const {
    return this->general_elements[0];
}

Element Field::primitive_element() const {
    return this->general_elements[1];
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////