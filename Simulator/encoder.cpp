////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: encoder.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "simulator.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void encode(ofstream& log_file) {
    string gen_file_name = "gen_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string enc_file_name = "enc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream gen_file("../Outputs/" + gen_file_name + ".txt");
    ofstream enc_file("../Outputs/" + enc_file_name + ".txt");
    for (int i = 0; i < ITERATION; i ++) {
        cout << get_current_time() << "\tEncoding at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        log_file << get_current_time() << "\tEncoding at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        string gen_data_line;
        getline(gen_file, gen_data_line);
        Polynomial message = convert_string_to_polynomial(gen_data_line);
        Polynomial codeword = RS0.systematic_encode(message);
        enc_file << convert_polynomial_to_string(codeword) << endl;
    }
    gen_file.close();
    enc_file.close();
    separate(enc_file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////