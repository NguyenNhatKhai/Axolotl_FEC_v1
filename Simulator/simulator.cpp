////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: simulator.cpp
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "simulator.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    string log_file_name = "simulator_" + get_current_time();
    ofstream log_file("../Logs/" + log_file_name + ".txt");
    try {
        // generate(log_file);
        // encode(log_file);
        
        channel_random_1s1p(log_file, {
            0, 0, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 6, 4, 2, 0, 0,
            0, 0, 12, 12, 0, 0, 0, 0, 0, 0, 8, 6, 4, 2, 0, 0,
            0, 0, 18, 18, 0, 0, 0, 0, 0, 0, 0, 0
        });

        syndrome_calculate(log_file);
        bma_solve(log_file);
        error_evaluate(log_file);
        error_correct(log_file);
        decode(log_file);
    } catch (const char* error_message) {
        cout << error_message << endl;
    }
    log_file.close();
    cout << "Harw" << endl;
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// int main() {
//     try {
//         // get_egf_pri_pol();
//         // get_egf_gen_ele();
//         // get_egf_inv_map_ele();
//         // get_rsc_gen_pol();
//         report_general_element();
//         // report_generator_polynomial();
//     } catch (const char* error_message) {
//         cout << error_message << endl;
//     }
//     cout << "Harw" << endl;
//     return 0;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////