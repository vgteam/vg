//
//  suffix_tree.cpp
//  
//
//  Created by Jordan Eizenga on 1/25/17.
//
//

#include <stdio.h>
#include <random>
#include <chrono>
#include <vector>
#include "suffix_tree.hpp"

#include "catch.hpp"

using namespace std;

namespace vg {
    namespace unittest {
        
        size_t longest_overlap(string& str1, string& str2) {
            size_t max_possible = min(str1.size(), str2.size());
            
            for (size_t i = str1.size() - max_possible; i < str1.size(); i++) {
                bool match = true;
                for (size_t j = i; j < str1.size(); j++) {
                    if (str1[j] != str2[j - i]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    return str1.size() - i;
                }
            }
            return 0;
        }
        
        vector<size_t> substring_locations(string& str, string& substr) {
            
            vector<size_t> locations;
            
            if (substr.empty()) {
                return locations;
            }
            
            for (size_t i = 0; i <= str.size() - substr.size(); i++) {
                bool match = true;
                for (size_t j = 0; j < substr.size(); j++) {
                    if (str[i + j] != substr[j]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    locations.push_back(i);
                }
            }
            return locations;
        }
        
        string random_string(string& alphabet, size_t length) {
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> distr(0, alphabet.size() - 1);
            
            string str;
            str.reserve(length);
            for (int i = 0; i < length; i++) {
                str.push_back(alphabet[distr(gen)]);
            }
            
            return str;
        }
        
        string random_repetive_string(vector<string>& chunks, size_t chunk_length,
                                      string& alphabet, double mismatch_rate) {
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> char_distr(0, alphabet.size() - 1);
            uniform_int_distribution<int> chunk_distr(0, chunks.size() - 1);
            uniform_real_distribution<double> err_distr(0.0, 1.0);
            
            string str;
            for (int i = 0; i < chunk_length; i++) {
                string chunk = chunks[chunk_distr(gen)];
                for (int j = 0; j < chunk.size(); j++) {
                    if (err_distr(gen) < mismatch_rate) {
                        chunk[j] = alphabet[char_distr(gen)];
                    }
                }
                str.append(chunk);
            }
            return str;
        }
        
        string random_substring(string& str, size_t substring_length, string& alphabet, double mismatch_rate) {
            
            if (str.size() < substring_length) {
                return "";
            }
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> start_distr(0, str.size() - substring_length);
            uniform_int_distribution<int> char_distr(0, alphabet.size() - 1);
            uniform_real_distribution<double> err_distr(0.0, 1.0);
            
            string substr;
            substr.reserve(substring_length);
            
            int substr_start = start_distr(gen);
            for (int i = substr_start; i < substr_start + substring_length; i++) {
                if (err_distr(gen) < mismatch_rate) {
                    substr.push_back(alphabet[char_distr(gen)]);
                }
                else {
                    substr.push_back(str[i]);
                }
            }
            return substr;
        }
        
        TEST_CASE("Suffix tree builds without crashing", "[suffix]") {
            
            string str = "GNANNNCGGATGCANNNCAGTGTCTTN";
            SuffixTree suffix_tree(str.begin(), str.end());
            
        }
        
        TEST_CASE("Suffix tree correctly finds longest overlap in hand-selected cases", "[suffix]") {
            
            SECTION("Suffix tree correctly finds longest overlap between two sequences") {
                
                string seq = "ACGTGACA";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                REQUIRE(suffix_tree.longest_overlap("ACAGCCT") == 3);
            }
            
            SECTION("Suffix tree correctly finds longest overlap in identical sequence") {
                
                string seq = "AATGGCATTNCGNAAGTACAGTG";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                REQUIRE(suffix_tree.longest_overlap(seq) == seq.size());
            }
            
            SECTION("Suffix tree correctly finds longest overlap onto empty sequence") {
                
                string seq = "AATGGCATTNCGNAAGTACAGTG";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                REQUIRE(suffix_tree.longest_overlap("") == 0);
            }
            
            SECTION("Suffix tree correctly finds longest overlap from empty sequence") {
                
                string seq = "";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                REQUIRE(suffix_tree.longest_overlap("AATGGCATTNCGNAAGTACAGTG") == 0);
            }
            
            SECTION("Suffix tree correctly finds longest overlap between two empty sequences") {
                
                string seq = "";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                REQUIRE(suffix_tree.longest_overlap("") == 0);
            }
        }
        
        TEST_CASE("Suffix tree correctly finds longest overlaps in 1000 random sequences", "[suffix]") {

            int max_str_len = 30;
            int num_seqs = 1000;
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> len_distr(0, max_str_len);

            string alphabet = "ACGTN";

            for (int i = 0; i < num_seqs; i++) {
                string str1 = random_string(alphabet, len_distr(gen));
                string str2 = random_string(alphabet, len_distr(gen));
                
                SuffixTree suffix_tree(str1.begin(), str1.end());

                size_t suffix_tree_longest_overlap = suffix_tree.longest_overlap(str2);
                size_t brute_longest_overlap = longest_overlap(str1, str2);

                if (suffix_tree_longest_overlap != brute_longest_overlap) {
                    // print out the failures since their random and we might have a hard time finding them again
                    cerr << "FAILURE: wrong overlap of " << suffix_tree_longest_overlap << " on " << str1 << " " << str2 << endl;
                }

                REQUIRE(suffix_tree_longest_overlap == brute_longest_overlap);
            }
        }
        
        TEST_CASE("Suffix tree correctly finds longest overlaps in 1000 low entropy random sequences", "[suffix]") {
            
            int max_num_chunks = 10;
            int max_chunk_len = 10;
            int num_chunks = 2;
            int num_seqs = 1000;
            double mismatch_rate = .01;
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> len_distr(0, max_chunk_len);
            uniform_int_distribution<int> num_chunks_distr(0, max_num_chunks);
            
            
            string alphabet = "ACGTN";
            
            for (int i = 0; i < num_seqs; i++) {
                
                vector<string> chunks;
                for (int j = 0; j < num_chunks; j++) {
                    chunks.push_back(random_string(alphabet, len_distr(gen)));
                }
                
                string str1 = random_repetive_string(chunks, num_chunks_distr(gen), alphabet, mismatch_rate);
                string str2 = random_repetive_string(chunks, num_chunks_distr(gen), alphabet, mismatch_rate);
                
                SuffixTree suffix_tree(str1.begin(), str1.end());
                
                size_t suffix_tree_longest_overlap = suffix_tree.longest_overlap(str2);
                size_t brute_longest_overlap = longest_overlap(str1, str2);
                
                if (suffix_tree_longest_overlap != brute_longest_overlap) {
                    // print out the failures since their random and we might have a hard time finding them again
                    cerr << "FAILURE: wrong overlap on " << str1 << " " << str2 << endl;
                }
                
                REQUIRE(suffix_tree_longest_overlap == brute_longest_overlap);
            }
        }
        
        TEST_CASE("Suffix tree correctly finds substring locations in hand-selected cases", "[suffix]") {
            
            SECTION("Suffix tree correctly finds substring locations in an example sequence") {
                
                string seq = "AGTGCGATAGATGATAGAAGATCGCTCGCTCCGCGATA";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                vector<size_t> locs = suffix_tree.substring_locations("GATA");
                sort(locs.begin(), locs.end());
                
                vector<size_t> correct_locs {5, 12, 34};
                
                REQUIRE(locs == correct_locs);
            }
            
            SECTION("Suffix tree correctly finds substring locations when substring does not occur") {
                
                string seq = "TACGGCAGATG";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                vector<size_t> locs = suffix_tree.substring_locations("GATA");
                sort(locs.begin(), locs.end());
                
                vector<size_t> correct_locs;
                
                REQUIRE(locs == correct_locs);
            }
            
            SECTION("Suffix tree correctly finds substring locations when substring is equal to suffix tree string") {
                
                string seq = "TACGGCAGATG";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                vector<size_t> locs = suffix_tree.substring_locations("TACGGCAGATG");
                sort(locs.begin(), locs.end());
                
                vector<size_t> correct_locs {0};
                
                REQUIRE(locs == correct_locs);
            }
            
            SECTION("Suffix tree correctly finds substring locations when substring is empty") {
                
                string seq = "TACGGCAGATG";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                vector<size_t> locs = suffix_tree.substring_locations("");
                sort(locs.begin(), locs.end());
                
                vector<size_t> correct_locs;
                
                REQUIRE(locs == correct_locs);
            }
            
            SECTION("Suffix tree correctly finds substring locations when substring is longer than suffix tree string") {
                
                string seq = "TACGGCAGATG";
                
                SuffixTree suffix_tree(seq.begin(), seq.end());
                
                vector<size_t> locs = suffix_tree.substring_locations("TACGGCAGATGA");
                sort(locs.begin(), locs.end());
                
                vector<size_t> correct_locs;
                
                REQUIRE(locs == correct_locs);
            }
        }
        
        TEST_CASE("Suffix tree correctly finds substring locations in 1000 randomly generated cases", "[suffix]") {
            
            int num_suffix_trees = 100;
            int num_substrings_per_tree = 10;
            int min_suffix_tree_length = 100;
            int max_suffix_tree_length = 300;
            int max_substring_length = 40;
            double mismatch_rate = .03;
            
            string alphabet = "ACGTN";
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> str_len_distr(min_suffix_tree_length, max_suffix_tree_length);
            uniform_int_distribution<int> substr_len_distr(0, max_substring_length);
            
            for (int i = 0; i < num_suffix_trees; i++) {
                
                string str = random_string(alphabet, substr_len_distr(gen));
                
                SuffixTree suffix_tree(str.begin(), str.end());
                
                for (int j = 0; j < num_substrings_per_tree; j++) {
                    
                    string substr = random_substring(str, substr_len_distr(gen), alphabet, mismatch_rate);
                    
                    vector<size_t> st_locations = suffix_tree.substring_locations(substr);
                    sort(st_locations.begin(), st_locations.end());
                    vector<size_t> direct_locations = substring_locations(str, substr);
                    
                    if (st_locations != direct_locations) {
                        // print out the failures since their random and we might have a hard time finding them again
                        cerr << "FAILURE: wrong substring locations on " << str << " " << substr << endl;
                    }
                    
                    REQUIRE(st_locations == direct_locations);
                }
            }
        }
        
        TEST_CASE("Suffix tree correctly finds substring locations in 1000 randomly generated low-entropy cases", "[suffix]") {
            
            int num_suffix_trees = 100;
            int num_substrings_per_tree = 10;
            int chunk_bank_size = 3;
            int min_num_chunks = 3;
            int max_num_chunks = 10;
            int min_chunk_length = 5;
            int max_chunk_length = 50;
            int max_substring_length = 40;
            double chunk_mismatch_rate = .05;
            double substring_mismatch_rate = .02;
            
            string alphabet = "ACGTN";
            
            random_device rd;
            default_random_engine gen(rd());
            uniform_int_distribution<int> num_chunks_distr(min_num_chunks, max_num_chunks);
            uniform_int_distribution<int> chunk_len_distr(min_chunk_length, max_chunk_length);
            uniform_int_distribution<int> substr_len_distr(0, max_substring_length);
            
            for (int i = 0; i < num_suffix_trees; i++) {
                
                vector<string> chunks;
                for (int j = 0; j < chunk_bank_size; j++) {
                    chunks.push_back(random_string(alphabet, chunk_len_distr(gen)));
                }
                
                string str = random_repetive_string(chunks, num_chunks_distr(gen), alphabet, chunk_mismatch_rate);
                
                SuffixTree suffix_tree(str.begin(), str.end());
                
                for (int j = 0; j < num_substrings_per_tree; j++) {
                    
                    string substr = random_substring(str, substr_len_distr(gen), alphabet, substring_mismatch_rate);
                    
                    vector<size_t> st_locations = suffix_tree.substring_locations(substr);
                    sort(st_locations.begin(), st_locations.end());
                    vector<size_t> direct_locations = substring_locations(str, substr);
                    
                    if (st_locations != direct_locations) {
                        // print out the failures since their random and we might have a hard time finding them again
                        cerr << "FAILURE: wrong substring locations on " << str << " " << substr << endl;
                    }
                    
                    REQUIRE(st_locations == direct_locations);
                }
            }
        }
    }
}
