////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: decoder.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "simulator.h"

////////////////////////////////////////////////////////////////////////////////////////////////////

void syndrome_calculate(ofstream& log_file) {
    string cha_file_name = "cha_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string syn_cal_file_name = "syn_cal_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream cha_file("../Outputs/" + cha_file_name + ".txt");
    ofstream syn_cal_file("../Outputs/" + syn_cal_file_name + ".txt");
    for (int i = 0; i < ITERATION; i ++) {
        cout << get_current_time() << "\tSyndrome calculating at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        log_file << get_current_time() << "\tSyndrome calculating at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        string cha_data_line;
        getline(cha_file, cha_data_line);
        Polynomial received = convert_string_to_polynomial(cha_data_line);
        Polynomial syndrome = RS0.syndrome(RS0.syndrome(received));
        syn_cal_file << convert_polynomial_to_string(syndrome) << endl;
    }
    cha_file.close();
    syn_cal_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void bma_solve(ofstream& log_file) {
    string syn_cal_file_name = "syn_cal_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_dis_file_name = "bma_sol_dis_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_err_loc_file_name = "bma_sol_err_loc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_aux_dis_file_name = "bma_sol_aux_dis_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_aux_err_loc_file_name = "bma_sol_aux_err_loc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_err_loc_len_file_name = "bma_sol_err_loc_len_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream syn_cal_file("../Outputs/" + syn_cal_file_name + ".txt");
    ofstream bma_sol_dis_file("../Outputs/" + bma_sol_dis_file_name + ".txt");
    ofstream bma_sol_err_loc_file("../Outputs/" + bma_sol_err_loc_file_name + ".txt");
    ofstream bma_sol_aux_dis_file("../Outputs/" + bma_sol_aux_dis_file_name + ".txt");
    ofstream bma_sol_aux_err_loc_file("../Outputs/" + bma_sol_aux_err_loc_file_name + ".txt");
    ofstream bma_sol_err_loc_len_file("../Outputs/" + bma_sol_err_loc_len_file_name + ".txt");
    vector<Element> discrepancies;
    vector<Polynomial> error_locations;
    vector<Element> auxiliary_discrepancies;
    vector<Polynomial> auxiliary_error_locations;
    vector<int> error_location_lengths;
    for (int i = 0; i < ITERATION; i ++) {
        discrepancies.clear();
        error_locations.clear();
        auxiliary_discrepancies.clear();
        auxiliary_error_locations.clear();
        error_location_lengths.clear();
        cout << get_current_time() << "\tBMA solving at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        log_file << get_current_time() << "\tBMA solving at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        string syn_cal_data_line;
        getline(syn_cal_file, syn_cal_data_line);
        Polynomial syndrome = convert_string_to_polynomial(syn_cal_data_line);
        Polynomial error_location = RS0.bma_error_location(syndrome.coefficients, discrepancies, error_locations, auxiliary_discrepancies, auxiliary_error_locations, error_location_lengths);
        #ifdef TRACK
            for (int j = 0; j < discrepancies.size(); j ++) {
                bma_sol_dis_file << convert_element_to_string(discrepancies[j]) << endl;
            }
            for (int j = 0; j < error_locations.size(); j ++) {
                bma_sol_err_loc_file << convert_polynomial_to_string(error_locations[j]) << endl;
            }
            for (int j = 0; j < auxiliary_discrepancies.size(); j ++) {
                bma_sol_aux_dis_file << convert_element_to_string(auxiliary_discrepancies[j]) << endl;
            }
            for (int j = 0; j < auxiliary_error_locations.size(); j ++) {
                bma_sol_aux_err_loc_file << convert_polynomial_to_string(auxiliary_error_locations[j]) << endl;
            }
            for (int j = 0; j < error_location_lengths.size(); j ++) {
                bma_sol_err_loc_len_file << error_location_lengths[j] << endl;
            }
        #else
            bma_sol_err_loc_file << convert_polynomial_to_string(error_locations[error_locations.size() - 1]) << endl;
            bma_sol_aux_dis_file << convert_element_to_string(auxiliary_discrepancies[auxiliary_discrepancies.size() - 1]) << endl;
            bma_sol_aux_err_loc_file << convert_polynomial_to_string(auxiliary_error_locations[auxiliary_error_locations.size() - 1]) << endl;
            bma_sol_err_loc_len_file << error_location_lengths[error_location_lengths.size() - 1] << endl;
        #endif
    }
    syn_cal_file.close();
    bma_sol_dis_file.close();
    bma_sol_err_loc_file.close();
    bma_sol_aux_dis_file.close();
    bma_sol_aux_err_loc_file.close();
    bma_sol_err_loc_len_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void error_evaluate(ofstream& log_file) {
    string bma_sol_err_loc_file_name = "bma_sol_err_loc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_aux_dis_file_name = "bma_sol_aux_dis_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_aux_err_loc_file_name = "bma_sol_aux_err_loc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string bma_sol_err_loc_len_file_name = "bma_sol_err_loc_len_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_sea_odd_file_name = "err_eva_sea_odd_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_sea_file_name = "err_eva_sea_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_num_file_name = "err_eva_num_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_aux_err_loc_file_name = "err_eva_aux_err_loc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_eva_file_name = "err_eva_eva_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_file_name = "err_eva_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream bma_sol_err_loc_file("../Outputs/" + bma_sol_err_loc_file_name + ".txt");
    ifstream bma_sol_aux_dis_file("../Outputs/" + bma_sol_aux_dis_file_name + ".txt");
    ifstream bma_sol_aux_err_loc_file("../Outputs/" + bma_sol_aux_err_loc_file_name + ".txt");
    ifstream bma_sol_err_loc_len_file("../Outputs/" + bma_sol_err_loc_len_file_name + ".txt");
    ofstream err_eva_sea_odd_file("../Outputs/" + err_eva_sea_odd_file_name + ".txt");
    ofstream err_eva_sea_file("../Outputs/" + err_eva_sea_file_name + ".txt");
    ofstream err_eva_num_file("../Outputs/" + err_eva_num_file_name + ".txt");
    ofstream err_eva_aux_err_loc_file("../Outputs/" + err_eva_aux_err_loc_file_name + ".txt");
    ofstream err_eva_eva_file("../Outputs/" + err_eva_eva_file_name + ".txt");
    ofstream err_eva_file("../Outputs/" + err_eva_file_name + ".txt");
    vector<Element> searched_odd_values;
    vector<Element> searched_values;
    vector<Element> evaluated_numerator_values;
    vector<Element> evaluated_auxiliary_error_location_values;
    vector<Element> evaluated_values;
    for (int i = 0; i < ITERATION; i ++) {
        searched_odd_values.clear();
        searched_values.clear();
        evaluated_numerator_values.clear();
        evaluated_auxiliary_error_location_values.clear();
        evaluated_values.clear();
        cout << get_current_time() << "\tError evaluating at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        log_file << get_current_time() << "\tError evaluating at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        string bma_sol_err_loc_data_line;
        getline(bma_sol_err_loc_file, bma_sol_err_loc_data_line);
        Polynomial error_location = convert_string_to_polynomial(bma_sol_err_loc_data_line);
        string bma_sol_aux_dis_data_line;
        getline(bma_sol_aux_dis_file, bma_sol_aux_dis_data_line);
        Element auxiliary_discrepancy = convert_string_to_element(bma_sol_aux_dis_data_line);
        string bma_sol_aux_err_loc_data_line;
        getline(bma_sol_aux_err_loc_file, bma_sol_aux_err_loc_data_line);
        Polynomial auxiliary_error_location = convert_string_to_polynomial(bma_sol_aux_err_loc_data_line);
        int error_location_length;
        bma_sol_err_loc_len_file >> error_location_length;
        vector<Element> error_numbers = RS0.chien_search(error_location, searched_odd_values, searched_values);
        vector<Element> error_values = RS0.horiguchi_koetter_algorithm(auxiliary_discrepancy, auxiliary_error_location, error_location_length, searched_odd_values, searched_values, evaluated_numerator_values, evaluated_auxiliary_error_location_values, evaluated_values);
        for (int j = 0; j < searched_odd_values.size(); j ++) {
            err_eva_sea_odd_file << convert_element_to_string(searched_odd_values[(j + 1) % searched_odd_values.size()]) << endl;
        }
        for (int j = 0; j < searched_values.size(); j ++) {
            err_eva_sea_file << convert_element_to_string(searched_values[(j + 1) % searched_values.size()]) << endl;
        }
        for (int j = 0; j < evaluated_numerator_values.size(); j ++) {
            err_eva_num_file << convert_element_to_string(evaluated_numerator_values[(j + 1) % evaluated_numerator_values.size()]) << endl;
        }
        for (int j = 0; j < evaluated_auxiliary_error_location_values.size(); j ++) {
            err_eva_aux_err_loc_file << convert_element_to_string(evaluated_auxiliary_error_location_values[(j + 1) % evaluated_auxiliary_error_location_values.size()]) << endl;
        }
        for (int j = 0; j < evaluated_values.size(); j ++) {
            err_eva_eva_file << convert_element_to_string(evaluated_values[(j + 1) % evaluated_values.size()]) << endl;
        }
        for (int j = 0; j < evaluated_values.size(); j ++) {
            if (searched_values[(j + 1) % evaluated_values.size()] == RS0.symbol_field->zero_element()) {
                err_eva_file << convert_element_to_string(evaluated_values[(j + 1) % evaluated_values.size()]) << endl;
            } else {
                err_eva_file << convert_element_to_string(RS0.symbol_field->zero_element()) << endl;
            }
        }
    }
    bma_sol_err_loc_file.close();
    bma_sol_aux_dis_file.close();
    bma_sol_aux_err_loc_file.close();
    bma_sol_err_loc_len_file.close();
    err_eva_sea_odd_file.close();
    err_eva_sea_file.close();
    err_eva_num_file.close();
    err_eva_aux_err_loc_file.close();
    err_eva_eva_file.close();
    err_eva_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void error_correct(ofstream& log_file) {
    string gen_file_name = "gen_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string cha_file_name = "cha_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_file_name = "err_eva_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_cor_file_name = "err_cor_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream gen_file("../Outputs/" + gen_file_name + ".txt");
    ifstream cha_file("../Outputs/" + cha_file_name + ".txt");
    ifstream err_eva_file("../Outputs/" + err_eva_file_name + ".txt");
    ofstream err_cor_file("../Outputs/" + err_cor_file_name + ".txt");
    vector<Element> evaluated_values;
    for (int i = 0; i < ITERATION; i ++) {
        evaluated_values.clear();
        cout << get_current_time() << "\tError correcting at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        log_file << get_current_time() << "\tError correcting at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
        string gen_data_line;
        getline(gen_file, gen_data_line);
        Polynomial message = convert_string_to_polynomial(gen_data_line);
        string cha_data_line;
        getline(cha_file, cha_data_line);
        Polynomial received = convert_string_to_polynomial(cha_data_line);
        for (int j = 0; j < RS0.symbol_field->size() - 1; j ++) {
            string err_eva_data_line;
            getline(err_eva_file, err_eva_data_line);
            evaluated_values.insert(evaluated_values.begin(), convert_string_to_element(err_eva_data_line));
        }
        Polynomial estimated_error(&FIE0, evaluated_values);
        Polynomial estimated_codeword = RS0.estimated_codeword(received, estimated_error);
        Polynomial estimated_message = RS0.estimated_message(estimated_codeword);
        err_cor_file << convert_polynomial_to_string(estimated_message) << endl;
        if (message == estimated_message) {
            cout << "\t\t\t\tError correcting passed" << endl;
        } else {
            cout << "\t\t\t\tError correcting failed" << endl;
        }
    }
    cha_file.close();
    err_eva_file.close();
    separate(err_cor_file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

void decode(ofstream& log_file) {
    string err_cor_file_name = "err_cor_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string dec_file_name = "dec_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream err_cor_file("../Outputs/" + err_cor_file_name + ".txt");
    ofstream dec_file("../Outputs/" + dec_file_name + ".txt");
    dec_file << err_cor_file.rdbuf();
    err_cor_file.close();
    dec_file.close();
    separate(dec_file_name);
}

////////////////////////////////////////////////////////////////////////////////////////////////////