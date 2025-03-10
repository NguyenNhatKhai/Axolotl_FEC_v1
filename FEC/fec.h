////////////////////////////////////////////////////////////////////////////////////////////////////
//
// File: fec.h
// Author: Nhat Khai Nguyen
//
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _FEC_H_
#define _FEC_H_

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../FFA/ffa.h"
#include "../Maths/maths.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

class RS;

////////////////////////////////////////////////////////////////////////////////////////////////////

class RS {
    public:
    Field* symbol_field;
    Polynomial generator_polynomial;

    public:
    RS() = delete;
    RS(Field* symbol_field, int correction_capability);
    ~RS() = default;

    public:
    int codeword_length() const;
    int message_length() const;
    int parity_length() const;
    int symbol_size() const;
    int detection_capability() const;
    int correction_capability() const;

    public:
    Polynomial systematic_encode(const Polynomial& message) const;
    Polynomial nonsystematic_encode(const Polynomial& message) const;

    public:
    Polynomial add_error(const Polynomial& codeword, const Polynomial& error) const;

    public:
    Polynomial pgz_decode(const Polynomial& received) const;
    Polynomial bm_decode(const Polynomial& received) const;
    Polynomial euclidean_decode(const Polynomial& received) const;

    public:
    vector<Element> syndrome(const Polynomial& received) const;
    Polynomial syndrome(const vector<Element>& syndrome) const;

    public:
    Polynomial pgz_error_location(const vector<Element>& syndrome) const;
    Polynomial bm_error_location(const vector<Element>& syndrome) const;
    Polynomial bma_error_location(const vector<Element>& syndrome) const;
    Polynomial euclidean_error_location(const Polynomial& syndrome) const;

    public:
    Polynomial bma_error_location(const vector<Element>& syndrome, vector<Element>& discrepancies, vector<Polynomial>& error_locations, vector<Element>& auxiliary_discrepancies, vector<Polynomial>& auxiliary_error_locations, vector<int>& error_location_lengths) const;

    public:
    Polynomial error_evaluator(const Polynomial& syndrome, const Polynomial& error_location) const;
    vector<Element> chien_search(const Polynomial& error_location) const;
    vector<Element> forney_algorithm(const Polynomial& error_location, const Polynomial& error_evaluator, const vector<Element>& error_numbers) const;

    public:
    vector<Element> chien_search(const Polynomial& error_location, vector<Element>& searched_odd_values, vector<Element>& searched_values) const;
    vector<Element> horiguchi_koetter_algorithm(const Element& auxiliary_discrepancy, const Polynomial& auxiliary_error_location, const int& error_location_length, const vector<Element>& searched_odd_values, const vector<Element>& searched_values, vector<Element>& evaluated_auxiliary_discrepancy_values, vector<Element>& evaluated_auxiliary_error_location_values, vector<Element>& evaluated_values) const;

    public:
    Polynomial estimated_error(const vector<Element>& error_numbers, const vector<Element>& error_values) const;
    Polynomial estimated_codeword(const Polynomial& received, const Polynomial& estimated_error) const;
    Polynomial estimated_message(const Polynomial& estimated_codeword) const;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#include "rs.cpp"

#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////