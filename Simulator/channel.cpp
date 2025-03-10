////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: channel.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "simulator.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

vector<int> get_random_index(int size, int iteration) {
    uniform_int_distribution<> random_codeword_index(0, size);
    vector<int> indices(size);
    for (int i = 0; i < size; i ++) {
        indices[i] = i;
    }
    shuffle(indices.begin(), indices.end(), random_engine);
    return vector<int>(indices.begin(), indices.begin() + iteration);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void channel_random_1s1p(ofstream& log_file, vector<int> error_symbols) {
    uniform_int_distribution<> random_error(0, RS0.symbol_field->size() - 2);
    string enc_file_name = "enc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string cha_file_name = "cha_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream enc_file("../Outputs/" + enc_file_name + ".txt");
    ofstream cha_file("../Outputs/" + cha_file_name + ".txt");
    for (int i = 0; i < ITERATION; i ++) {
        cout << get_current_time() << "\tTransmitting at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        log_file << get_current_time() << "\tTransmitting at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        string line;
        getline(enc_file, line);
        Polynomial codeword = convert_string_to_polynomial(line);
        vector<Element> error_coefficients;
        vector<int> random_symbol_indexs = get_random_index(RS0.codeword_length(), error_symbols[i]);
        for (int j = 0; j < RS0.codeword_length(); j ++) {
            if (find(random_symbol_indexs.begin(), random_symbol_indexs.end(), j) != random_symbol_indexs.end()) {
                cout << "\t\t\t\tError occuring at symbol " << j << endl;
                log_file << "\t\t\t\tError occuring at symbol " << j << endl;
                error_coefficients.push_back(ELE(random_error(random_engine)));
            } else {
                error_coefficients.push_back(RS0.symbol_field->zero_element());
            }
        }
        codeword = RS0.add_error(codeword, Polynomial(&FIE0, error_coefficients));
        cha_file << convert_polynomial_to_string(codeword) << endl;
    }
    enc_file.close();
    cha_file.close();
    separate(cha_file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////