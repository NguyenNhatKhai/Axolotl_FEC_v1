////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: helper.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "simulator.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void separate(string file_name) {
    ifstream pre_separate_file("../Outputs/" + file_name + ".txt");
    string content;
    for (int i = 0; i < ITERATION; i ++) {
        string line;
        getline(pre_separate_file, line);
        content += line;
    }
    pre_separate_file.close();
    ofstream post_separate_file("../Outputs/" + file_name + "_" + to_string(INTERFACE) + ".txt");
    size_t start_index = 0;
    while (start_index < content.size()) {
        int end_index = min(start_index + INTERFACE * RS0.symbol_size(), content.size());
        string line = content.substr(start_index, end_index - start_index);
        line.append(INTERFACE * RS0.symbol_size() - line.size(), '0');
        start_index = end_index;
        post_separate_file << line << endl;
    }
    post_separate_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

string get_current_time() {
    time_t now_time_t = time(nullptr);
    tm now_tm = *localtime(&now_time_t);
    timespec now_timespec;
    clock_gettime(CLOCK_REALTIME, &now_timespec);
    stringstream now_stringstream;
    now_stringstream << put_time(&now_tm, "%Y_%m_%d_%H_%M_%S") << "_" << setw(9) << setfill('0') << now_timespec.tv_nsec;
    return now_stringstream.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

string get_ordinal_suffix(int number) {
    if (number % 10 == 1 && number % 100 != 11) {
        return "st";
    } else if (number % 10 == 2 && number % 100 != 12) {
        return "nd";
    } else if (number % 10 == 3 && number % 100 != 13) {
        return "rd";
    } else {
        return "th";
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

string convert_element_to_string(Element element) {
    #ifdef DEBUG
        if (element.field != &FIE0) {
            throw "Simulator\\Helper\\convert_element_to_string(Element)\\element\\field";
        }
    #endif
    string element_string = "";
    for (int i = 0; i < element.size(); i ++) {
        element_string += to_string(element.values[element.size() - i - 1].value);
    }
    return element_string;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

string convert_polynomial_to_string(Polynomial polynomial) {
    #ifdef DEBUG
        if (polynomial.field != &FIE0) {
            throw "Simulator\\Helper\\convert_polynomial_to_string(Polynomial)\\polynomial\\field";
        }
    #endif
    string polynomial_string = "";
    for (int i = 0; i < (int)polynomial.coefficients.size(); i ++) {
        polynomial_string += convert_element_to_string(polynomial.coefficients[(int)polynomial.coefficients.size() - i - 1]);
    }
    return polynomial_string;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Element convert_string_to_element(string element_string) {
    #ifdef DEBUG
        if (element_string.length() != FIE0.primitive_polynomial->degree()) {
            throw "Simulator\\Helper\\convert_string_to_element(string)\\element_string\\length";
        }
    #endif
    vector<Element> values;
    for (int i = 0; i < element_string.length(); i ++) {
        if (element_string[element_string.length() - i - 1] == '0') {
            values.push_back(BIN0);
        } else if (element_string[element_string.length() - i - 1] == '1') {
            values.push_back(BIN1);
        } else {
            throw "Simulator\\Helper\\convert_string_to_element(string)\\element_string";
        }
    }
    return Element(&FIE0, values);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial convert_string_to_polynomial(string polynomial_string) {
    #ifdef DEBUG
        if (polynomial_string.length() % FIE0.primitive_polynomial->degree() != 0) {
            throw "Simulator\\Helper\\convert_string_to_polynomial(string)\\polynomial_string\\length";
        }
    #endif
    vector<Element> coefficients;
    for (int i = 0; i < polynomial_string.length() / FIE0.primitive_polynomial->degree(); i ++) {
        coefficients.push_back(convert_string_to_element(polynomial_string.substr((polynomial_string.length() / FIE0.primitive_polynomial->degree() - i - 1) * FIE0.primitive_polynomial->degree(), FIE0.primitive_polynomial->degree())));
    }
    return Polynomial(&FIE0, coefficients);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void get_egf_pri_pol() {
    cout << "'{";
    for (int i = 0; i <= FIE0.primitive_polynomial->degree(); i ++) {
        cout << FIE0.primitive_polynomial->coefficients[FIE0.primitive_polynomial->degree() - i];
        if (i < FIE0.primitive_polynomial->degree()) {
            cout << ", ";
        }
    }
    cout << "}" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void get_egf_gen_ele() {
    cout << "'{";
    for (int i = 0; i < FIE0.size() - 1; i ++) {
        cout << "'b" << convert_element_to_string(FIE0.general_elements[FIE0.size() - 2 - i]);
        if (i != FIE0.size() - 2) {
            cout << ", ";
        }
    }
    cout << "}" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void get_egf_inv_map_ele() {
    vector<Element> temp_0_elements(FIE0.size(), FIE0.zero_element());
    for (int i = 0; i < FIE0.size() - 1; i ++) {
        temp_0_elements[stoi(convert_element_to_string(FIE0.general_elements[i]), nullptr, 2)] = ~FIE0.general_elements[i];
    }
    cout << "'{";
    for (int i = 0; i < temp_0_elements.size(); i ++) {
        cout << "'b" << convert_element_to_string(temp_0_elements[temp_0_elements.size() - 1 - i]);
        if (i != temp_0_elements.size() - 1) {
            cout << ", ";
        }
    }
    cout << "}" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void get_rsc_gen_pol() {
    cout << "'{";
    for (int i = 0; i <= RS0.generator_polynomial.degree(); i ++) {
        cout << "'b" << convert_element_to_string(RS0.generator_polynomial.coefficients[RS0.generator_polynomial.degree() - i]);
        if (i < RS0.generator_polynomial.degree()) {
            cout << ", ";
        }
    }
    cout << "}" << endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////