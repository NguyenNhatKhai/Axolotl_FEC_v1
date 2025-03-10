////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: generator.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "simulator.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

Polynomial get_random_message() {
    uniform_int_distribution<> random_element(0, RS0.symbol_field->size() - 1);
    vector<Element> message_coefficients;
    for (int i = 0; i < RS0.message_length(); i ++) {
        message_coefficients.push_back(ELE(random_element(random_engine)));
    }
    return Polynomial(&FIE0, message_coefficients);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void generate(ofstream& log_file) {
    string gen_file_name = "gen_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ofstream gen_file("../Outputs/" + gen_file_name + ".txt");
    for (int i = 0; i < ITERATION; i ++) {
        cout << get_current_time() << "\tGenerating at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        log_file << get_current_time() << "\tGenerating at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        Polynomial message = get_random_message();
        gen_file << convert_polynomial_to_string(message) << endl;
    }
    gen_file.close();
    separate(gen_file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////