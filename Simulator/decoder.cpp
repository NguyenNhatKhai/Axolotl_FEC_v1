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
    string err_eva_aux_dis_file_name = "err_eva_aux_dis_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_aux_err_loc_file_name = "err_eva_aux_err_loc_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    string err_eva_eva_file_name = "err_eva_eva_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
    ifstream bma_sol_err_loc_file("../Outputs/" + bma_sol_err_loc_file_name + ".txt");
    ifstream bma_sol_aux_dis_file("../Outputs/" + bma_sol_aux_dis_file_name + ".txt");
    ifstream bma_sol_aux_err_loc_file("../Outputs/" + bma_sol_aux_err_loc_file_name + ".txt");
    ifstream bma_sol_err_loc_len_file("../Outputs/" + bma_sol_err_loc_len_file_name + ".txt");
    ofstream err_eva_sea_odd_file("../Outputs/" + err_eva_sea_odd_file_name + ".txt");
    ofstream err_eva_sea_file("../Outputs/" + err_eva_sea_file_name + ".txt");
    ofstream err_eva_aux_dis_file("../Outputs/" + err_eva_aux_dis_file_name + ".txt");
    ofstream err_eva_aux_err_loc_file("../Outputs/" + err_eva_aux_err_loc_file_name + ".txt");
    ofstream err_eva_eva_file("../Outputs/" + err_eva_eva_file_name + ".txt");
    vector<Element> searched_odd_values;
    vector<Element> searched_values;
    vector<Element> evaluated_auxiliary_discrepancy_values;
    vector<Element> evaluated_auxiliary_error_location_values;
    vector<Element> evaluated_values;
    for (int i = 0; i < ITERATION; i ++) {
        searched_odd_values.clear();
        searched_values.clear();
        evaluated_auxiliary_discrepancy_values.clear();
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
        vector<Element> error_values = RS0.horiguchi_koetter_algorithm(auxiliary_discrepancy, auxiliary_error_location, error_location_length, searched_odd_values, searched_values, evaluated_auxiliary_discrepancy_values, evaluated_auxiliary_error_location_values, evaluated_values);
        for (int j = 0; j < searched_odd_values.size(); j ++) {
            err_eva_sea_odd_file << convert_element_to_string(searched_odd_values[(j + 1) % searched_odd_values.size()]) << endl;
        }
        for (int j = 0; j < searched_values.size(); j ++) {
            err_eva_sea_file << convert_element_to_string(searched_values[(j + 1) % searched_values.size()]) << endl;
        }
        for (int j = 0; j < evaluated_auxiliary_discrepancy_values.size(); j ++) {
            err_eva_aux_dis_file << convert_element_to_string(evaluated_auxiliary_discrepancy_values[(j + 1) % evaluated_auxiliary_discrepancy_values.size()]) << endl;
        }
        for (int j = 0; j < evaluated_auxiliary_error_location_values.size(); j ++) {
            err_eva_aux_err_loc_file << convert_element_to_string(evaluated_auxiliary_error_location_values[(j + 1) % evaluated_auxiliary_error_location_values.size()]) << endl;
        }
        for (int j = 0; j < evaluated_values.size(); j ++) {
            err_eva_eva_file << convert_element_to_string(evaluated_values[(j + 1) % evaluated_values.size()]) << endl;
        }
    }
    bma_sol_err_loc_file.close();
    bma_sol_aux_dis_file.close();
    bma_sol_aux_err_loc_file.close();
    bma_sol_err_loc_len_file.close();
    err_eva_sea_odd_file.close();
    err_eva_sea_file.close();
    err_eva_aux_dis_file.close();
    err_eva_aux_err_loc_file.close();
    err_eva_eva_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// void error_estimate() {
//     string bma_sol_err_loc_len_file_name = "bma_sol_err_loc_len_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
//     string rtc_sea_file_name = "rtc_sea_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
//     string err_eva_file_name = "err_eva_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
//     string err_est_file_name = "err_est_data_" + to_string(RS0.codeword_length()) + "_" + to_string(RS0.message_length());
//     ifstream bma_sol_err_loc_len_file("../Outputs/" + bma_sol_err_loc_len_file_name + ".txt");
//     ifstream rtc_sea_file("../Outputs/" + rtc_sea_file_name + ".txt");
//     ifstream err_eva_file("../Outputs/" + err_eva_file_name + ".txt");
//     ofstream err_est_file("../Outputs/" + err_est_file_name + ".txt");
//     for (int i = 0; i < ITERATION; i ++) {
//         cout << get_current_time() << "\tError estimating at the " << i + 1 << get_ordinal_suffix(i + 1) << " iteration" << endl;
//         vector<Element> error_numbers;
//         vector<Element> error_values;
//         string rtc_sea_data_line;
//         string err_eva_data_line;
//         for (int j = 0; j < RS0.symbol_field->size() - 1; j ++) {
//             getline(rtc_sea_file, rtc_sea_data_line);
//             getline(err_eva_file, err_eva_data_line);
//             if (convert_string_to_element(rtc_sea_data_line) == RS0.symbol_field->zero_element()) {
//                 error_numbers.push_back(~RS0.symbol_field->general_elements[RS0.symbol_field->size() - 2 - j]);
//                 error_values.push_back(convert_string_to_element(err_eva_data_line));
//             }
//         }
//         int error_location_length;
//         bma_sol_err_loc_len_file >> error_location_length;
//         Polynomial estimated_error(&FIE0, vector<Element>(RS0.codeword_length(), RS0.symbol_field->zero_element()));
//         if (error_numbers.size() == error_location_length) {
//             estimated_error = RS0.estimated_error(error_numbers, error_values);
//         }
//         err_est_file << convert_polynomial_to_string(estimated_error).substr(0, RS0.message_length() * RS0.symbol_size()) << endl;
//     }
//     bma_sol_err_loc_len_file.close();
//     rtc_sea_file.close();
//     err_eva_file.close();
//     err_est_file.close();
// }

////////////////////////////////////////////////////////////////////////////////////////////////////