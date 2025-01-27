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
#include <vector>
#include <sstream>
#include <set>
#include <vg/vg.pb.h>
#include "utility.hpp"
#include "snarl_distance_index.hpp"
#include "minimizer_mapper.hpp"
#include "integrated_snarl_finder.hpp"
#include "genotypekit.hpp"
#include "traversal_finder.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>
#include "../primer_filter.hpp"
#include "../recombinator.hpp"

using namespace std;

namespace vg {

/**
 * Primer struct contains primer attributes, including sequence, left/right primer,
 * position on the reference genome, length, index offset in corresponding node on
 * sequence graph, and vector of corresponding nodes on the sequence graph. Everything
 * is in the positive/forward orientation.
 */
struct Primer {
    string sequence;
    bool left = true;
    size_t position_chromosome = numeric_limits<size_t>::max();
    size_t position_template   = numeric_limits<size_t>::max();
    size_t length              = numeric_limits<size_t>::max();
    size_t offset              = numeric_limits<size_t>::max();
    vector<size_t> mapped_nodes_ids;
};

/**
 * Primer_pair struct contains primer pair attributesm including left primer, right primer,
 * linear product size, minimum and maximum product size on the sequence graph, and boolean on
 * whether the primers locate in low variation region of the sequence graph.
 */
struct PrimerPair {
    Primer left_primer;
    Primer right_primer;
    string chromosome_name;
    string template_feature;
    size_t linear_product_size = numeric_limits<size_t>::max();
    size_t template_position   = numeric_limits<size_t>::max();
    size_t min_product_size    = numeric_limits<size_t>::max();
    size_t max_product_size    = numeric_limits<size_t>::max();
    double variation_level     = 0.0;
    vector<HaplotypePartitioner::sequence_type> sequence_visits;
};

class PrimerFinder {

private:
    unordered_map<string, vector<PrimerPair>> chroms; // map containing a vector of primer pairs for each chromosome
    const PathPositionHandleGraph* graph;
    const SnarlDistanceIndex* distance_index;
    MinimizerMapper* giraffe_mapper; 
    const gbwtgraph::GBWTGraph& gbwt_graph;
    const gbwt::GBWT& gbwt_index;
    const gbwt::FastLocate& r_index;


public:
    PrimerFinder() = default;
    
    /**
     * Construct Primer finder given PathPositionHandleGraph, reference graph name
     * and pointer to SnarlDistanceIndex
     */
    PrimerFinder(const handlegraph::PathPositionHandleGraph* graph_param,
        const SnarlDistanceIndex* distance_index_param,
        istream& primers_file_handle,
        const gbwtgraph::GBWTGraph& gbwt_graph, const gbwt::GBWT& gbwt_index,
        const gbwt::FastLocate& r_index, MinimizerMapper* giraffe_mapper_param=nullptr);

    /**
     * Destructor
     */
    ~PrimerFinder();

    /**
     * Add a Primer_pair object given primers' starting node id, offset relative
     * to the starting node, and length, all in the POSTIVE orientation. The new
     * primer_pair object is automatically added to primer_pairs vector - and
     * selected_primer_pairs if conditions are met. Mainly used for unit testing.
     */
    void add_primer_pair(const string& path_name, const size_t& left_primer_starting_node_id,
        const size_t& left_primer_offset, const size_t& left_primer_length,
        const size_t& right_primer_starting_node_id,
        const size_t& right_primer_offset, const size_t& right_primer_length);

    /**
     * Read the path to the primer3 output. Primers information is parsed,
     * processed, and  stored in primer_pairs vector - and selected_primer_pairs
     * if conditions are met.
     */
    void load_primers(istream& file_handle);

    /**
     * return vector of Primer pairs
     */
    const vector<PrimerPair>& get_primer_pairs_of_chrom(const string& chrom_name) const;

    /**
     * return the total number of reference paths
     */
    const size_t total_reference_paths() const;

    vector<string> get_reference_paths(); 

private:
    /**
     * Private functions used by public or private functions.
     */

    /**
     * Get the graph coordinates by mapping and surjecting the template
     * To be used if the chromosome_name isn't a valid path
     * Returns a pair of the path/chromosome name and the offset of the template in the path
     * Used in: load_primers
     */
     std::pair<string, size_t> get_graph_coordinates_from_sequence(const string& seq);


    /**
     * Update minimum and maximum prodcut to a primer pair object.
     * Used in: add_primer_pair
     *          load_primers
     */
    void update_min_max_product_size(PrimerPair& primer_pair);

    /**
     * Update a Primer object given starting node id, offset relative to the starting node,
     * and the length of primer.
     * Used in: add_primer_pair
     */
    void make_primer(Primer& primer, const string& path_name, const size_t& starting_node_id,
        const size_t& offset, const size_t& length, const bool& is_left);

    /**
     * Find and store corresponding node ids to Primer object.
     * Used in: make_primer
     *          load_primers
     */
    void map_to_nodes(Primer& primer, const string& path_name);

    /**
     * Find the length of the longest match between two sequences. Also find and
     * store offset in Primer object.
     * Used in: map_to_nodes
     */
    size_t longest_match_len(Primer& primer, const string& left_seq, const string& right_seq,
       const bool& first_node);

    /**
     * Strip empty spaces on the sides of a string.
     * Used in: load_primers
     */
    const string strip(const string& s) const;

    /**
     * Check if primers in a primer_pair object have variations on the pangenome.
     * Used in: add_primer_node
     *          load_primers
     */
    // const bool no_variation_at_primers(const PrimerPair& primer_pair) const;

    void update_variation(PrimerPair& primer_pair, const string& path_name);
    
    /**
     * Split a string into vectors.
     * Used in: load_primers
     */
    vector<string> split(const string& str);

    /**
     * Split a string into vectors given delimiter.
     */
    vector<string> split(const string& str, const char& delim);
    
    /**
     * Works like str.startswith(prefix) in python
     * Used in: load_primers
     */
    bool startswith(const string& str, const string& prefix);
};

}

#endif /* primer_filter_hpp */
