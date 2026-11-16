#include "deamination_profile.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include <zlib.h>

namespace vg {

using namespace std;

namespace {

/// Returns the index of a base in the order A, C, G, T, counting anything else as A.
size_t base_index(char base) {
    switch (base) {
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 0;
    }
}

bool is_n(char base) {
    return base == 'N' || base == 'n';
}

/// Reads a whole file, which may be gzip-compressed. Throws std::runtime_error if it cannot.
string read_possibly_compressed_file(const string& filename) {
    gzFile file = gzopen(filename.c_str(), "rb");
    if (file == nullptr) {
        throw runtime_error("cannot open deamination profile " + filename);
    }
    string text;
    array<char, 1 << 16> buffer;
    int read_bytes;
    while ((read_bytes = gzread(file, buffer.data(), buffer.size())) > 0) {
        text.append(buffer.data(), read_bytes);
    }
    bool failed = read_bytes < 0;
    gzclose(file);
    if (failed) {
        throw runtime_error("cannot read deamination profile " + filename);
    }
    return text;
}

vector<string> split_tabs(string_view line) {
    vector<string> fields;
    size_t start = 0;
    while (true) {
        size_t tab = line.find('\t', start);
        fields.emplace_back(line.substr(start, tab == string_view::npos ? string_view::npos : tab - start));
        if (tab == string_view::npos) {
            break;
        }
        start = tab + 1;
    }
    return fields;
}

/// The number of columns in a profile, one for each change of a base into another.
constexpr size_t PROFILE_COLUMNS = 12;

/// Parses one rate of a profile. Throws std::runtime_error if the field is not a number.
double parse_rate(const string& field, const string& source) {
    double rate = 0.0;
    const char* end = field.data() + field.size();
    auto [parsed_end, error] = from_chars(field.data(), end, rate);
    if (field.empty() || error != errc() || parsed_end != end) {
        throw runtime_error("deamination profile " + source + " has a field that is not a number: " + field);
    }
    return rate;
}

/// Turns the fields of one row of a profile into the substitution matrix at its position. Throws
/// std::runtime_error if a field is not a number or the rates of changing one base sum to more
/// than 1.
SubstitutionMatrix parse_profile_row(const vector<string>& fields, const string& line, const string& source) {
    // Column c holds the rate of changing true base c / 3 into the c % 3-th of the other three
    // bases, in the order A, C, G, T. The rate of reading a base unchanged is what is left.
    SubstitutionMatrix matrix;
    size_t column = 0;
    for (size_t true_base = 0; true_base < 4; true_base++) {
        double unchanged = 1.0;
        for (size_t read_base = 0; read_base < 4; read_base++) {
            if (read_base != true_base) {
                matrix[true_base][read_base] = parse_rate(fields[column], source);
                unchanged -= matrix[true_base][read_base];
                column++;
            }
        }
        if (unchanged < 0) {
            throw runtime_error("deamination profile " + source
                                + " has rates for one base that sum to more than 1: " + line);
        }
        matrix[true_base][true_base] = unchanged;
    }
    return matrix;
}

}

vector<SubstitutionMatrix> DeaminationProfile::parse_profile(const string& text, const string& source) {
    vector<SubstitutionMatrix> result;
    istringstream lines(text);
    string line;
    bool header = true;
    while (getline(lines, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        vector<string> fields = split_tabs(line);
        if (fields.size() != PROFILE_COLUMNS) {
            throw runtime_error("deamination profile " + source + " has a line without "
                                + to_string(PROFILE_COLUMNS) + " tab-separated fields: " + line);
        }
        if (header) {
            header = false;
        } else {
            result.push_back(parse_profile_row(fields, line, source));
        }
    }
    if (result.empty()) {
        throw runtime_error("deamination profile " + source + " has no rows of rates");
    }
    return result;
}

DeaminationProfile::DeaminationProfile(const string& five_prime_filename, const string& three_prime_filename) :
    DeaminationProfile(parse_profile(read_possibly_compressed_file(five_prime_filename), five_prime_filename),
                       parse_profile(read_possibly_compressed_file(three_prime_filename), three_prime_filename)) {
    // Nothing else to do
}

DeaminationProfile::DeaminationProfile(vector<SubstitutionMatrix> five_prime, vector<SubstitutionMatrix> three_prime) :
    five_prime(std::move(five_prime)), three_prime(std::move(three_prime)) {
    if (this->five_prime.empty() || this->three_prime.empty()) {
        throw invalid_argument("DeaminationProfile needs rates for at least one position from each end");
    }
}

double DeaminationProfile::substitution_probability(size_t fragment_length, size_t position,
                                                    char true_base, char read_base) const {
    if (is_n(true_base) || is_n(read_base)) {
        return 0.25;
    }
    size_t from_three_prime = position < fragment_length ? fragment_length - 1 - position : 0;
    const SubstitutionMatrix& five = five_prime[min(position, five_prime.size() - 1)];
    const SubstitutionMatrix& three = three_prime[min(from_three_prime, three_prime.size() - 1)];
    size_t t = base_index(true_base);
    size_t r = base_index(read_base);
    const SubstitutionMatrix& chosen = five[t][t] <= three[t][t] ? five : three;
    return chosen[t][r];
}

double damage_posterior(const string& graph_kmer, const string& read_kmer, size_t read_length,
                        size_t kmer_start, const DeaminationProfile& profile,
                        const DamageModelParameters& parameters) {
    if (graph_kmer == read_kmer) {
        return 1.0;
    }

    // Log-likelihood of the read k-mer under the damage model.
    double damage_log_likelihood = 0.0;
    size_t compared = min(graph_kmer.size(), read_kmer.size());
    for (size_t i = 0; i < compared; i++) {
        double probability = profile.substitution_probability(read_length, kmer_start, graph_kmer[i], read_kmer[i]);
        if (probability == 0.0) {
            return 0.0;
        }
        damage_log_likelihood += log(probability);
    }

    // Log-likelihood of the read k-mer under the spurious model.
    size_t mismatches = 0;
    for (size_t i = 0; i < compared; i++) {
        if (graph_kmer[i] != read_kmer[i]) {
            mismatches++;
        }
    }
    double spurious_log_likelihood;
    if (read_kmer.size() < 22) {
        auto n = static_cast<double>(graph_kmer.size());
        double rate = min(1.0, max(0.0, 1.0014 * pow(n, -0.6628)));
        auto m = static_cast<double>(mismatches);
        spurious_log_likelihood = lgamma(n + 1) - lgamma(m + 1) - lgamma(n - m + 1)
                                  + m * log(rate) + (n - m) * log(1.0 - rate);
    } else {
        spurious_log_likelihood = log(0.01);
    }

    if (!isfinite(damage_log_likelihood) || !isfinite(spurious_log_likelihood)) {
        return 1e-8;
    }

    double damage_log_posterior = damage_log_likelihood + log(1.0 - parameters.spurious_prior);
    double spurious_log_posterior = spurious_log_likelihood + log(parameters.spurious_prior);
    double largest = max(damage_log_posterior, spurious_log_posterior);
    double log_evidence = largest + log(exp(damage_log_posterior - largest) + exp(spurious_log_posterior - largest));
    return exp(damage_log_posterior - log_evidence);
}

}
