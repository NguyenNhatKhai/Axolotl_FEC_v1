////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: rs.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "fec.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

RS::RS(Field* symbol_field, int correction_capability) {
    #ifdef DEBUG
        if (correction_capability <= 0 || 2 * correction_capability >= symbol_field->size() - 1) {
            throw "FEC\\RS\\RS(Field*, int)\\correction_capability";
        }
    #endif
    this->symbol_field = symbol_field;
    this->generator_polynomial = Polynomial(this->symbol_field, {this->symbol_field->unit_element()});
    for (int i = 0; i < 2 * correction_capability; i ++) {
        this->generator_polynomial = this->generator_polynomial * Polynomial(this->symbol_field, {this->symbol_field->general_elements[i + 1], this->symbol_field->unit_element()});
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int RS::codeword_length() const {
    return this->symbol_field->size() - 1;
}

int RS::message_length() const {
    return this->symbol_field->size() - this->generator_polynomial.degree() - 1;
}

int RS::parity_length() const {
    return this->codeword_length() - this->message_length();
}

int RS::symbol_size() const {
    return this->symbol_field->primitive_polynomial->degree();
}

int RS::detection_capability() const {
    return this->parity_length();
}

int RS::correction_capability() const {
    return this->parity_length() / 2;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::systematic_encode(const Polynomial& message) const {
    #ifdef DEBUG
        if (message.coefficients.size() != this->message_length()) {
            throw "FEC\\RS\\systematic_encode(const Polynomial&, const string&)\\message\\coefficients\\size";
        } else if (*message.field != *this->symbol_field) {
            throw "FEC\\RS\\systematic_encode(const Polynomial&, const string&)\\message\\field";
        }
    #endif
    Polynomial temp_0_polynomial(this->symbol_field, vector<Element>(this->parity_length() + 1, this->symbol_field->zero_element()));
    temp_0_polynomial.coefficients[this->parity_length()] = this->symbol_field->unit_element();
    Polynomial temp_1_polynomial = message * temp_0_polynomial;
    temp_0_polynomial = (temp_1_polynomial + (temp_1_polynomial % this->generator_polynomial)).redegree(this->codeword_length() - 1);
    return temp_0_polynomial;
}

Polynomial RS::nonsystematic_encode(const Polynomial& message) const {
    #ifdef DEBUG
        if (message.coefficients.size() != this->message_length()) {
            throw "FEC\\RS\\nonsystematic_encode(const Polynomial&, const string&)\\message\\coefficients\\size";
        } else if (*message.field != *this->symbol_field) {
            throw "FEC\\RS\\nonsystematic_encode(const Polynomial&, const string&)\\message\\field";
        }
    #endif
    Polynomial temp_0_polynomial = (message * this->generator_polynomial).redegree(this->codeword_length() - 1);
    return temp_0_polynomial;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::add_error(const Polynomial& codeword, const Polynomial& error) const {
    #ifdef DEBUG
        if (codeword.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\codeword\\coefficients\\size";
        } else if (error.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\error\\coefficients\\size";
        } else if (*codeword.field != *this->symbol_field) {
            throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\codeword\\field";
        } else if (*error.field != *this->symbol_field) {
            throw "FEC\\RS\\add_error(const Polynomial&, const Polynomial&, const string&)\\error\\field";
        }
    #endif
    Polynomial temp_0_polynomial = (codeword + error).redegree(this->codeword_length() - 1);
    return temp_0_polynomial;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::pgz_decode(const Polynomial& received) const {
    #ifdef DEBUG
        if (received.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\pgz_decode(const Polynomial&, const string&)\\received\\coefficients\\size";
        } else if (*received.field != *this->symbol_field) {
            throw "FEC\\RS\\pgz_decode(const Polynomial&, const string&)\\received\\field";
        }
    #endif
    vector<Element> temp_0_syndrome = this->syndrome(received);
    Polynomial temp_1_syndrome = this->syndrome(temp_0_syndrome);
    Polynomial temp_2_error_location = this->pgz_error_location(temp_0_syndrome);
    Polynomial temp_3_error_evaluator = this->error_evaluator(temp_1_syndrome, temp_2_error_location);
    vector<Element> temp_4_error_numbers = this->chien_search(temp_2_error_location);
    vector<Element> temp_5_error_values = this->forney_algorithm(temp_2_error_location, temp_3_error_evaluator, temp_4_error_numbers);
    Polynomial temp_6_estimated_error = this->estimated_error(temp_4_error_numbers, temp_5_error_values);
    Polynomial temp_7_estimated_codeword = this->estimated_codeword(received, temp_6_estimated_error);
    Polynomial temp_8_estimated_message = this->estimated_message(temp_7_estimated_codeword);
    return temp_8_estimated_message;
}

Polynomial RS::bm_decode(const Polynomial& received) const {
    #ifdef DEBUG
        if (received.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\bm_decode(const Polynomial&, const string&)\\received\\coefficients\\size";
        } else if (*received.field != *this->symbol_field) {
            throw "FEC\\RS\\bm_decode(const Polynomial&, const string&)\\received\\field";
        }
    #endif
    vector<Element> temp_0_syndrome = this->syndrome(received);
    Polynomial temp_1_syndrome = this->syndrome(temp_0_syndrome);
    Polynomial temp_2_error_location = this->bm_error_location(temp_0_syndrome);
    Polynomial temp_3_error_evaluator = this->error_evaluator(temp_1_syndrome, temp_2_error_location);
    vector<Element> temp_4_error_numbers = this->chien_search(temp_2_error_location);
    vector<Element> temp_5_error_values = this->forney_algorithm(temp_2_error_location, temp_3_error_evaluator, temp_4_error_numbers);
    Polynomial temp_6_estimated_error = this->estimated_error(temp_4_error_numbers, temp_5_error_values);
    Polynomial temp_7_estimated_codeword = this->estimated_codeword(received, temp_6_estimated_error);
    Polynomial temp_8_estimated_message = this->estimated_message(temp_7_estimated_codeword);
    return temp_8_estimated_message;
}

Polynomial RS::euclidean_decode(const Polynomial& received) const {
    #ifdef DEBUG
        if (received.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\euclidean_decode(const Polynomial&, const string&)\\received\\coefficients\\size";
        } else if (*received.field != *this->symbol_field) {
            throw "FEC\\RS\\euclidean_decode(const Polynomial&, const string&)\\received\\field";
        }
    #endif
    vector<Element> temp_0_syndrome = this->syndrome(received);
    Polynomial temp_1_syndrome = this->syndrome(temp_0_syndrome);
    Polynomial temp_2_error_location = this->euclidean_error_location(temp_1_syndrome);
    Polynomial temp_3_error_evaluator = this->error_evaluator(temp_1_syndrome, temp_2_error_location);
    vector<Element> temp_4_error_numbers = this->chien_search(temp_2_error_location);
    vector<Element> temp_5_error_values = this->forney_algorithm(temp_2_error_location, temp_3_error_evaluator, temp_4_error_numbers);
    Polynomial temp_6_estimated_error = this->estimated_error(temp_4_error_numbers, temp_5_error_values);
    Polynomial temp_7_estimated_codeword = this->estimated_codeword(received, temp_6_estimated_error);
    Polynomial temp_8_estimated_message = this->estimated_message(temp_7_estimated_codeword);
    return temp_8_estimated_message;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

vector<Element> RS::syndrome(const Polynomial& received) const {
    #ifdef DEBUG
        if (received.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\syndrome(const Polynomial&)\\received\\coefficients\\size";
        } else if (*received.field != *this->symbol_field) {
            throw "FEC\\RS\\syndrome(const Polynomial&)\\received\\field";
        }
    #endif
    vector<Element> temp_0_elements(this->parity_length(), this->symbol_field->zero_element());
    for (int i = 0; i < temp_0_elements.size(); i ++) {
        temp_0_elements[i] = received.evaluate(this->symbol_field->general_elements[i + 1]);
    }
    return temp_0_elements;
}

Polynomial RS::syndrome(const vector<Element>& syndrome) const {
    #ifdef DEBUG
        if (syndrome.size() != 2 * this->correction_capability()) {
            throw "FEC\\RS\\syndrome(const vector<Element>&)\\syndrome\\size";
        } else if (true) {
            for (int i = 0; i < syndrome.size(); i ++) {
                if (*syndrome[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\syndrome(const vector<Element>&)\\syndrome\\field";
                }
            }
        }
    #endif
    return Polynomial(this->symbol_field, syndrome);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::pgz_error_location(const vector<Element>& syndrome) const {
    #ifdef DEBUG
        if (syndrome.size() != 2 * this->correction_capability()) {
            throw "FEC\\RS\\pgz_error_location(const vector<Element>&)\\syndrome\\size";
        } else if (true) {
            for (int i = 0; i < syndrome.size(); i ++) {
                if (*syndrome[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\pgz_error_location(const vector<Element>&)\\syndrome\\field";
                }
            }
        }
    #endif
    vector<Element> temp_0_coefficients(1, this->symbol_field->unit_element());
    for (int i = this->correction_capability(); i > 0; i --) {
        vector<vector<Element>> temp_1_elements(i, vector<Element>(i, this->symbol_field->zero_element()));
        for (int j = 0; j < i; j ++) {
            for (int k = 0; k < i; k ++) {
                temp_1_elements[j][k] = syndrome[j + k];
            }
        }
        Matrix temp_2_matrix(this->symbol_field, temp_1_elements);
        if (temp_2_matrix.determinant() != this->symbol_field->zero_element()) {
            vector<vector<Element>> temp_3_elements(i, vector<Element>(1, this->symbol_field->zero_element()));
            for (int j = 0; j < i; j ++) {
                temp_3_elements[j][0] = -syndrome[i + j];
            }
            Matrix temp_4_matrix = (~temp_2_matrix) * Matrix(this->symbol_field, temp_3_elements);
            for (int j = 1; j <= i; j ++) {
                temp_0_coefficients.push_back(temp_4_matrix.elements[i - j][0]);
            }
            break;
        }
    }
    return Polynomial(this->symbol_field, temp_0_coefficients).redegree(this->correction_capability());
}

Polynomial RS::bm_error_location(const vector<Element>& syndrome) const {
    #ifdef DEBUG
        if (syndrome.size() != 2 * this->correction_capability()) {
            throw "FEC\\RS\\bm_error_location(const vector<Element>&)\\syndrome\\size";
        } else if (true) {
            for (int i = 0; i < syndrome.size(); i ++) {
                if (*syndrome[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\bm_error_location(const vector<Element>&)\\syndrome\\field";
                }
            }
        }
    #endif
    vector<Polynomial> error_locations(2 * this->correction_capability() + 2, Polynomial(this->symbol_field, {this->symbol_field->zero_element()}));
    vector<Element> discrepancies(2 * this->correction_capability() + 2, this->symbol_field->zero_element());
    error_locations[0].coefficients[0] = this->symbol_field->unit_element();
    discrepancies[0] = this->symbol_field->unit_element();
    error_locations[1].coefficients[0] = this->symbol_field->unit_element();
    for (int i = 1; i < 2 * this->correction_capability() + 1; i ++) {
        for (int j = 0; j <= error_locations[i].degree(); j ++) {
            discrepancies[i] = discrepancies[i] + (error_locations[i].coefficients[j] * syndrome[i - j - 1]);
        }
        if (discrepancies[i] == this->symbol_field->zero_element()) {
            error_locations[i + 1] = error_locations[i];
        } else {
            int temp_0_int = 0;
            int temp_1_min = i + error_locations[0].degree();
            for (int j = 1; j < i; j ++) {
                if (discrepancies[j] != this->symbol_field->zero_element() && i - j + error_locations[j].degree() < temp_1_min) {
                    temp_0_int = j;
                    temp_1_min = i - j + error_locations[j].degree();
                }
            }
            Polynomial temp_2_polynomial(this->symbol_field, vector<Element>(i - temp_0_int + 1, this->symbol_field->zero_element()));
            temp_2_polynomial.coefficients[i - temp_0_int] = discrepancies[i] * (~discrepancies[temp_0_int]);
            error_locations[i + 1] = error_locations[i] + temp_2_polynomial * error_locations[temp_0_int];
        }
    }
    return error_locations[2 * this->correction_capability() + 1].redegree(this->correction_capability());
}

Polynomial RS::bma_error_location(const vector<Element>& syndrome) const {
    #ifdef DEBUG
        if (syndrome.size() != 2 * this->correction_capability()) {
            throw "FEC\\RS\\bma_error_location(const vector<Element>&)\\syndrome\\size";
        } else if (true) {
            for (int i = 0; i < syndrome.size(); i ++) {
                if (*syndrome[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\bma_error_location(const vector<Element>&)\\syndrome\\field";
                }
            }
        }
    #endif
    vector<Element> discrepancies;
    vector<Polynomial> error_locations;
    vector<Element> auxiliary_discrepancies;
    vector<Polynomial> auxiliary_error_locations;
    vector<int> error_location_lengths;
    return this->bma_error_location(syndrome, discrepancies, error_locations, auxiliary_discrepancies, auxiliary_error_locations, error_location_lengths);
}

Polynomial RS::euclidean_error_location(const Polynomial& syndrome) const {
    #ifdef DEBUG
        if (syndrome.degree() >= 2 * this->correction_capability()) {
            throw "FEC\\RS\\euclidean_error_location(const Polynomial&)\\syndrome\\degree";
        } else if (*syndrome.field != *this->symbol_field) {
            throw "FEC\\RS\\euclidean_error_location(const Polynomial&)\\syndrome\\field";
        }
    #endif
    Polynomial temp_0_polynomial(this->symbol_field, vector<Element>(2 * this->correction_capability() + 1, this->symbol_field->zero_element()));
    temp_0_polynomial.coefficients[2 * this->correction_capability()] = this->symbol_field->unit_element();
    Polynomial temp_1_polynomial = syndrome;
    Polynomial temp_2_polynomial(this->symbol_field, {this->symbol_field->zero_element()});
    Polynomial temp_3_polynomial(this->symbol_field, {this->symbol_field->unit_element()});
    while (temp_1_polynomial.degree() >= this->correction_capability()) {
        temp_2_polynomial = temp_2_polynomial - (temp_0_polynomial / temp_1_polynomial) * temp_3_polynomial;
        swap(temp_2_polynomial, temp_3_polynomial);
        temp_0_polynomial = temp_0_polynomial % temp_1_polynomial;
        swap(temp_0_polynomial, temp_1_polynomial);
    }
    return temp_3_polynomial.redegree(this->correction_capability());
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::bma_error_location(const vector<Element>& syndrome, vector<Element>& discrepancies, vector<Polynomial>& error_locations, vector<Element>& auxiliary_discrepancies, vector<Polynomial>& auxiliary_error_locations, vector<int>& error_location_lengths) const {
    #ifdef DEBUG
        if (syndrome.size() != 2 * this->correction_capability()) {
            throw "FEC\\RS\\bma_error_location(const vector<Element>&, vector<Element>&, vector<Polynomial>&, vector<Element>&, vector<Polynomial>&, vector<int>&)\\syndrome\\size";
        } else if (true) {
            for (int i = 0; i < syndrome.size(); i ++) {
                if (*syndrome[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\bma_error_location(const vector<Element>&, vector<Element>&, vector<Polynomial>&, vector<Element>&, vector<Polynomial>&, vector<int>&)\\syndrome\\field";
                }
            }
        }
    #endif
    Element discrepancy = this->symbol_field->zero_element();
    Polynomial error_location(this->symbol_field, {this->symbol_field->unit_element()});
    Element auxiliary_discrepancy = this->symbol_field->unit_element();
    Polynomial auxiliary_error_location(this->symbol_field, {this->symbol_field->unit_element()});
    int error_location_length = 0;
    error_locations.push_back(error_location);
    auxiliary_discrepancies.push_back(auxiliary_discrepancy);
    auxiliary_error_locations.push_back(auxiliary_error_location);
    error_location_lengths.push_back(error_location_length);
    for (int i = 0; i < 2 * this->correction_capability(); i ++) {
        discrepancy = this->symbol_field->zero_element();
        for (int j = 0; j <= error_location_length; j ++) {
            discrepancy = discrepancy + (error_location.coefficients[j] * syndrome[i - j]);
        }
        discrepancies.push_back(discrepancy);
        Polynomial temp_0_polynomial(this->symbol_field, {this->symbol_field->zero_element(), this->symbol_field->unit_element()});
        Polynomial new_error_location = error_location + ((temp_0_polynomial * auxiliary_error_location) * (discrepancy * ~auxiliary_discrepancy));
        if (discrepancy != this->symbol_field->zero_element() && 2 * error_location_length <= i) {
            auxiliary_error_location = error_location;
            auxiliary_discrepancy = discrepancy;
            error_location_length = i + 1 - error_location_length;
        } else if (i - error_location_length == this->correction_capability() - 1) {
            auxiliary_error_location = auxiliary_error_location;
            auxiliary_discrepancy = auxiliary_discrepancy;
            error_location_length = error_location_length;
        } else {
            auxiliary_error_location = temp_0_polynomial * auxiliary_error_location;
            auxiliary_discrepancy = auxiliary_discrepancy;
            error_location_length = error_location_length;
        }
        error_location = new_error_location.redegree(error_location_length);
        error_locations.push_back(error_location);
        auxiliary_discrepancies.push_back(auxiliary_discrepancy);
        auxiliary_error_locations.push_back(auxiliary_error_location);
        error_location_lengths.push_back(error_location_length);
        if (i - error_location_length == this->correction_capability() - 1) {
            break;
        }
    }
    return error_location;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::error_evaluator(const Polynomial& syndrome, const Polynomial& error_location) const {
    #ifdef DEBUG
        if (syndrome.degree() >= 2 * this->correction_capability()) {
            throw "FEC\\RS\\error_evaluator(const Polynomial&, const Polynomial&)\\syndrome\\degree";
        } else if (error_location.degree() > this->correction_capability()) {
            throw "FEC\\RS\\error_evaluator(const Polynomial&, const Polynomial&)\\error_location\\degree";
        } else if (*syndrome.field != *this->symbol_field) {
            throw "FEC\\RS\\error_evaluator(const Polynomial&, const Polynomial&)\\syndrome\\field";
        } else if (*error_location.field != *this->symbol_field) {
            throw "FEC\\RS\\error_evaluator(const Polynomial&, const Polynomial&)\\error_location\\field";
        }
    #endif
    vector<Element> temp_0_coefficients(2 * this->correction_capability() + 1, this->symbol_field->zero_element());
    temp_0_coefficients[2 * this->correction_capability()] = this->symbol_field->unit_element();
    Polynomial temp_1_polynomial(this->symbol_field, temp_0_coefficients);
    return ((syndrome * error_location) % temp_1_polynomial).redegree(this->correction_capability() - 1);
}

vector<Element> RS::chien_search(const Polynomial& error_location) const {
    #ifdef DEBUG
        if (error_location.degree() > this->correction_capability()) {
            throw "FEC\\RS\\chien_search(const Polynomial&)\\error_location\\degree";
        } else if (*error_location.field != *this->symbol_field) {
            throw "FEC\\RS\\chien_search(const Polynomial&)\\error_location\\field";
        }
    #endif
    vector<Element> temp_0_elements;
    for (int i = 0; i < this->symbol_field->size() - 1; i ++) {
        if (error_location.evaluate(this->symbol_field->general_elements[i]) == this->symbol_field->zero_element()) {
            temp_0_elements.push_back(this->symbol_field->general_elements[i]);
        }
    }
    return temp_0_elements;
}

vector<Element> RS::forney_algorithm(const Polynomial& error_location, const Polynomial& error_evaluator, const vector<Element>& error_numbers) const {
    #ifdef DEBUG
        if (error_location.degree() > this->correction_capability()) {
            throw "FEC\\RS\\forney_algorithm(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_location\\degree";
        } else if (error_evaluator.degree() >= 2 * this->correction_capability()) {
            throw "FEC\\RS\\forney_algorithm(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\degree";
        } else if (error_numbers.size() > this->correction_capability()) {
            throw "FEC\\RS\\forney_algorithm(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_numbers\\size";
        } else if (*error_location.field != *this->symbol_field) {
            throw "FEC\\RS\\forney_algorithm(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_location\\field";
        } else if (*error_evaluator.field != *this->symbol_field) {
            throw "FEC\\RS\\forney_algorithm(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_evaluator\\field";
        } else if (true) {
            for (int i = 0; i < error_numbers.size(); i ++) {
                if (*error_numbers[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\forney_algorithm(const Polynomial&, const Polynomial&, const vector<Element>&)\\error_numbers\\field";
                }
            }
        }
    #endif
    vector<Element> temp_0_elements;
    for (int i = 0; i < error_numbers.size(); i ++) {
        temp_0_elements.push_back(error_evaluator.evaluate(error_numbers[i]) / error_location.derivative().evaluate(error_numbers[i]));
    }
    return temp_0_elements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

vector<Element> RS::chien_search(const Polynomial& error_location, vector<Element>& searched_odd_values, vector<Element>& searched_values) const {
    #ifdef DEBUG
        if (error_location.degree() > this->correction_capability()) {
            throw "FEC\\RS\\chien_search(const Polynomial&, vector<Element>&, vector<Element>&)\\error_location\\degree";
        } else if (*error_location.field != *this->symbol_field) {
            throw "FEC\\RS\\chien_search(const Polynomial&, vector<Element>&, vector<Element>&)\\error_location\\field";
        }
    #endif
    vector<Element> temp_0_elements;
    Polynomial temp_1_polynomial(this->symbol_field, {this->symbol_field->zero_element(), this->symbol_field->unit_element()});
    Polynomial temp_2_polynomial = temp_1_polynomial * error_location.derivative();
    for (int i = 0; i < this->symbol_field->size() - 1; i ++) {
        Element searched_odd_value = temp_2_polynomial.evaluate(this->symbol_field->general_elements[i]);
        searched_odd_values.push_back(searched_odd_value);
        Element searched_value = error_location.evaluate(this->symbol_field->general_elements[i]);
        searched_values.push_back(searched_value);
        if (searched_value == this->symbol_field->zero_element()) {
            temp_0_elements.push_back(this->symbol_field->general_elements[i]);
        }
    }
    return temp_0_elements;
}

vector<Element> RS::horiguchi_koetter_algorithm(const Element& auxiliary_discrepancy, const Polynomial& auxiliary_error_location, const int& error_location_length, const vector<Element>& searched_odd_values, const vector<Element>& searched_values, vector<Element>& evaluated_auxiliary_discrepancy_values, vector<Element>& evaluated_auxiliary_error_location_values, vector<Element>& evaluated_values) const {
    #ifdef DEBUG
        if (*auxiliary_discrepancy.field != *this->symbol_field) {
            throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\auxiliary_discrepancy\\field";
        } else if (auxiliary_error_location.degree() > this->correction_capability() - 1) {
            throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\auxiliary_error_location\\degree";
        } else if (*auxiliary_error_location.field != *this->symbol_field) {
            throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\auxiliary_error_location\\field";
        } else if (error_location_length < 0) {
            throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\error_location_length";
        } else if (searched_values.size() != this->symbol_field->size() - 1) {
            throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\searched_values\\size";
        } else if (true) {
            for (int i = 0; i < searched_values.size(); i ++) {
                if (*searched_values[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\searched_values\\field";
                }
            }
        } else if (searched_odd_values.size() != this->symbol_field->size() - 1) {
            throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\searched_odd_values\\size";
        } else if (true) {
            for (int i = 0; i < searched_odd_values.size(); i ++) {
                if (*searched_odd_values[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\horiguchi_koetter_algorithm(const Element&, const Polynomial&, const int&, const vector<Element>&, const vector<Element>&, vector<Element>&, vector<Element>&, vector<Element>&)\\searched_odd_values\\field";
                }
            }
        }
    #endif
    vector<Element> temp_0_elements;
    for (int i = 0; i < this->symbol_field->size() - 1; i ++) {
        Element evaluated_auxiliary_discrepancy_value = auxiliary_discrepancy * (this->symbol_field->general_elements[i] ^ (this->correction_capability() + error_location_length - 1));
        Element evaluated_auxiliary_error_location_value = auxiliary_error_location.evaluate(this->symbol_field->general_elements[i]);
        Element temp_1_element = evaluated_auxiliary_error_location_value * searched_odd_values[i];
        Element evaluated_value = this->symbol_field->zero_element();
        if (temp_1_element != this->symbol_field->zero_element()) {
            evaluated_value = evaluated_auxiliary_discrepancy_value / temp_1_element;
        }
        evaluated_auxiliary_discrepancy_values.push_back(evaluated_auxiliary_discrepancy_value);
        evaluated_auxiliary_error_location_values.push_back(evaluated_auxiliary_error_location_value);
        evaluated_values.push_back(evaluated_value);
        if (searched_values[i] == this->symbol_field->zero_element()) {
            temp_0_elements.push_back(evaluated_value);
        }
    }
    return temp_0_elements;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial RS::estimated_error(const vector<Element>& error_numbers, const vector<Element>& error_values) const {
    #ifdef DEBUG
        if (error_numbers.size() > this->correction_capability() || error_numbers.size() != error_values.size()) {
            throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\error_numbers\\size";
        } else if (error_values.size() > this->correction_capability() || error_values.size() != error_numbers.size()) {
            throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\error_values\\size";
        } else if (true) {
            for (int i = 0; i < error_numbers.size(); i ++) {
                if (*error_numbers[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\error_numbers\\field";
                }
            }
        } else if (true) {
            for (int i = 0; i < error_values.size(); i ++) {
                if (*error_values[i].field != *this->symbol_field) {
                    throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)\\error_values\\field";
                }
            }
        }
    #endif
    vector<Element> temp_0_coefficients(this->codeword_length(), this->symbol_field->zero_element());
    for (int i = 0; i < error_numbers.size(); i ++) {
        Element temp_1_element = ~error_numbers[i];
        for (int j = 0; j < this->symbol_field->size() - 1; j ++) {
            if (temp_1_element == this->symbol_field->general_elements[j]) {
                temp_0_coefficients[j] = error_values[i];
                break;
            }
            if (j == this->symbol_field->size() - 2) {
                throw "FEC\\RS\\estimated_error(const vector<Element>&, const vector<Element>&)";
            }
        }
    }
    return Polynomial(this->symbol_field, temp_0_coefficients);
}

Polynomial RS::estimated_codeword(const Polynomial& received, const Polynomial& estimated_error) const {
    #ifdef DEBUG
        if (received.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\received\\coefficients\\size";
        } else if (estimated_error.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\estimated_error\\coefficients\\size";
        } else if (*received.field != *this->symbol_field) {
            throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\received\\field";
        } else if (*estimated_error.field != *this->symbol_field) {
            throw "FEC\\RS\\estimated_codeword(const Polynomial&, const Polynomial&)\\estimated_error\\field";
        }
    #endif
    return (received - estimated_error).redegree(this->codeword_length() - 1);
}

Polynomial RS::estimated_message(const Polynomial& estimated_codeword) const {
    #ifdef DEBUG
        if (estimated_codeword.coefficients.size() != this->codeword_length()) {
            throw "FEC\\RS\\estimated_message(const Polynomial&)\\estimated_codeword\\coefficients\\size";
        } else if (*estimated_codeword.field != *this->symbol_field) {
            throw "FEC\\RS\\estimated_message(const Polynomial&)\\estimated_codeword\\field";
        }
    #endif
    Polynomial temp_0_polynomial(this->symbol_field, vector<Element>(this->message_length(), this->symbol_field->zero_element()));
    for (int i = 0; i < this->message_length(); i ++) {
        temp_0_polynomial.coefficients[i] = estimated_codeword.coefficients[i + this->parity_length()];
    }
    return temp_0_polynomial;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////