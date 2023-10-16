// primer_filter.hpp
//
// Contains class Primer_finder for storing and filtering primers predicted 
// using Primer3. Also contains Primer struct and Primer_pair struct that stores
// information on primers and primer pairs.

#ifndef VG_PRIMER_FILTER_HPP_INCLUDED
#define VG_PRIMER_FILTER_HPP_INCLUDED

// Not sure what to include.. Will include everything from the unittest for now
#include <stdio.h>
#include <iostream>
#include <regex>
#include <fstream>
#include <vector>
#include <sstream>
#include <set>
#include <vg/vg.pb.h>
#include "snarl_distance_index.hpp"
#include "integrated_snarl_finder.hpp"
#include "genotypekit.hpp"
#include "traversal_finder.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>

using namespace std;

namespace vg {

struct Primer {
    string sequence;
    bool left = true;
    size_t position;
    size_t length;
    size_t offset;
    vector<size_t> mapped_nodes_ids;
};

struct Primer_pair {
    Primer left_primer;
    Primer right_primer;
    size_t linear_product_size;
    size_t min_product_size;
    size_t max_product_size;
};

class Primer_finder {

private:
    vector<Primer_pair> primer_pairs;
    vector<Primer_pair> selected_primer_pairs;
    PathPositionHandleGraph* graph;
    SnarlDistanceIndex* distance_index;
    path_handle_t reference_path_handle; 

public:
    Primer_finder() = default;
    
    /**
     * Construct Primer finder given PathPositionHandleGraph, reference graph name
     * and pointer to SnarlDistanceIndex
     */
    Primer_finder(unique_ptr<handlegraph::PathPositionHandleGraph>& graph_param,
                string reference_path_name, SnarlDistanceIndex* distance_index_param);

    /**
     * Destructor
     */
    ~Primer_finder();

    /**
     * Add a Primer_pair object given primers' starting node id, offset relative
     * to the starting node, and length, all in the POSTIVE orientation. The new
     * primer_pair object is automatically added to primer_pairs vector - and
     * selected_primer_pairs if conditions are met. Mainly used for unit testing.
     */
    void add_primer_pair(size_t left_primer_starting_node_id,
                    size_t left_primer_offset, size_t left_primer_length,
                    size_t right_primer_starting_node_id,
                    size_t right_primer_offset, size_t right_primer_length);

    /**
     * Read the path to the primer3 output. Primers information is parsed,
     * processed, and  stored in primer_pairs vector - and selected_primer_pairs
     * if conditions are met.
     */
    void load_primers(string path_to_primers);

    /**
     * return vector of Primer pairs
     */
    vector<Primer_pair> get_primer_pairs();

    /**
     * return vector selected primer pairs
     */
    vector<Primer_pair> get_selected_primer_pairs();

private:
    /**
     * Private functions used by public or private functions.
     */

    /**
     * Update minimum and maximum prodcut to a primer pair object.
     * Used in: add_primer_pair
     *          load_primers
     */
    void update_min_max_product_size(Primer_pair& primer_pair);

    /**
     * Update a Primer object given starting node id, offset relative to the starting node,
     * and the length of primer.
     * Used in: add_primer_pair
     */
    void make_primer(Primer& primer, size_t starting_node_id, size_t offset, size_t length, bool is_left);

    /**
     * Find and store corresponding node ids to Primer object.
     * Used in: make_primer
     *          load_primers
     */
    void map_to_nodes(Primer& primer);

    /**
     * Find the length of the longest match between two sequences. Also find and
     * store offset in Primer object.
     * Used in: map_to_nodes
     */
    size_t longest_match_len(Primer& primer, string const left_seq, string const right_seq,
        bool const first_node);

    /**
     * Strip empty spaces on the right side of a string.
     * Used in: load_primers
     */
    string rstrip(string const s);

    /**
     * Check if primers in a primer_pair object have variations on the pangenome.
     * Used in: add_primer_node
     *          load_primers
     */
    bool no_variation(const Primer_pair& primer_pair);

    /**
     * return the complement of a nucleotide.
     * Used in: revcomp
     */
    char complement(char nt);

    /**
     * return the reverse complement of a sequence.
     * Used in: make_primers
     *          map_to_nodes
     *      
     */
    string revcomp(string const seq);

    /**
     * Split a string into vectors.
     * Used in: load_priemrs
     */
    vector<string> split(string str, string const delim);

};

}

#endif /* primder_filter_hpp */