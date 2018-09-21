#ifndef VG_SIMULATOR_HPP_INCLUDED
#define VG_SIMULATOR_HPP_INCLUDED

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "xg.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "position.hpp"
#include "cached_position.hpp"
#include "xg_position.hpp"
#include "distributions.hpp"
#include "lru_cache.h"
#include "json2pb.h"

namespace vg {

using namespace std;

/// We have a utility function for turning positions along paths, with
/// orientations, into pos_ts. Remember that pos_t counts offset from the start
/// of the reoriented node, while here we count offset from the beginning of the
/// forward version of the path.
pos_t position_at(xg::XG* xgidx, const string& path_name, const size_t& path_offset, bool is_reverse);

/**
 * Generate Alignments (with or without mutations, and in pairs or alone) from
 * an XG index.
 */
class Sampler {

public:

    xg::XG* xgidx;
    // We need this so we don't re-load the node for every character we visit in
    // it.
    LRUCache<id_t, Node> node_cache;
    LRUCache<id_t, vector<Edge> > edge_cache;
    mt19937 rng;
    int64_t nonce;
    // If set, only sample positions/start reads on the forward strands of their
    // nodes.
    bool forward_only;
    // A flag that we set if we don't want to generate sequences with Ns (on by default)
    bool no_Ns;
    // A vector which, if nonempty, gives the names of the paths to restrict simulated reads to.
    vector<string> source_paths;
    vg::discrete_distribution<> path_sampler; // draw an index in source_paths
    inline Sampler(xg::XG* x,
            int seed = 0,
            bool forward_only = false,
            bool allow_Ns = false,
            const vector<string>& source_paths = {})
        : xgidx(x),
          node_cache(100),
          edge_cache(100),
          forward_only(forward_only),
          no_Ns(!allow_Ns),
          nonce(0),
          source_paths(source_paths) {
        if (!seed) {
            seed = time(NULL);
        }
        rng.seed(seed);
        set_source_paths(source_paths);
    }

    void set_source_paths(const vector<string>& source_paths);

    pos_t position(void);
    string sequence(size_t length);
    
    /// Get an alignment against the whole graph, or against the source path if
    /// one is selected.
    Alignment alignment(size_t length);
    
    /// Get an alignment against the whole graph.
    Alignment alignment_to_graph(size_t length);
    
    /// Get an alignment against the currently set source_path.
    Alignment alignment_to_path(const string& source_path, size_t length);
    
    Alignment alignment_with_error(size_t length,
                                   double base_error,
                                   double indel_error);
    vector<Alignment> alignment_pair(size_t read_length,
                                     size_t fragment_length,
                                     double fragment_std_dev,
                                     double base_error,
                                     double indel_error);
    size_t node_length(id_t id);
    char pos_char(pos_t pos);
    map<pos_t, char> next_pos_chars(pos_t pos);

    Alignment mutate(const Alignment& aln,
                     double base_error,
                     double indel_error);

    /**
     * Mutate the given edit, producing a vector of edits that should replace
     * it. Position is the position of the start of the edit, and is updated to
     * point to the next base after the mutated edit.
     */
    vector<Edit> mutate_edit(const Edit& edit,
                             const pos_t& position,
                             double base_error,
                             double indel_error,
                             const string& bases,
                             vg::uniform_real_distribution<double>& rprob,
                             vg::uniform_int_distribution<int>& rbase);

    string alignment_seq(const Alignment& aln);
    
    /// Return true if the alignment is semantically valid against the XG index
    /// we wrap, and false otherwise. Checks from_lengths on mappings to make
    /// sure all node bases are accounted for. Won't accept alignments with
    /// internal jumps between graph locations or regions; all skipped bases
    /// need to be accounted for by deletions.
    bool is_valid(const Alignment& aln);

};

/**
 * Class that simulates reads with alignments to a graph that mimic the error
 * profile of NGS sequencing data.
 */
class NGSSimulator {
public:
    /// Initialize simulator. FASTQ file will be used to train an error distribution.
    /// Most reads in the FASTQ should be the same length. Polymorphism rates apply
    /// uniformly along a read, whereas errors are distributed as indicated by the learned
    /// distribution. The simulation can also be restricted to named paths in the graph.
    NGSSimulator(xg::XG& xg_index,
                 const string& ngs_fastq_file,
                 bool interleaved_fastq = false,
                 const vector<string>& source_paths = {},
                 double substition_polymorphism_rate = 0.001,
                 double indel_polymorphism_rate = 0.0002,
                 double indel_error_proportion = 0.01,
                 double insert_length_mean = 1000.0,
                 double insert_length_stdev = 75.0,
                 double error_multiplier = 1.0,
                 bool retry_on_Ns = true,
                 size_t seed = 0);
    
    /// Sample an individual read and alignment
    Alignment sample_read();
    
    /// Sample a pair of reads an alignments
    pair<Alignment, Alignment> sample_read_pair();
    
private:
    template<class From, class To>
    class MarkovDistribution {
    public:
        MarkovDistribution(size_t seed);
        
        /// record a transition from the input data
        void record_transition(From from, To to);
        /// indicate that there is no more data and prepare for sampling
        void finalize();
        /// sample according to the training data
        To sample_transition(From from);
        
    private:
        
        default_random_engine prng;
        unordered_map<From, vg::uniform_int_distribution<size_t>> samplers;
        
        unordered_map<To, size_t> column_of;
        vector<To> value_at;
        unordered_map<From, vector<size_t>> cond_distrs;
        
    };
    
    NGSSimulator(void) = delete;
    
    /// DNA alphabet
    static const string alphabet;
    /// Remainder of the alphabet after removing a given character
    unordered_map<char, string> mutation_alphabets;
    
    /// Add a quality string to the training data
    void record_read_quality(const Alignment& aln, bool read_2 = false);
    /// Add a pair of quality strings to the training data
    void record_read_pair_quality(const Alignment& aln_1, const Alignment& aln_2);
    /// Indicate that there is no more training data
    void finalize();
    /// Get a quality string and a vector of 'N'-masks that mimics the training data
    pair<string, vector<bool>> sample_read_quality();
    /// Get a pair of quality strings and vectors of 'N'-masks that mimic the training data
    pair<pair<string, vector<bool>>, pair<string, vector<bool>>> sample_read_quality_pair();
    /// Wrapped internal function for quality sampling
    pair<string, vector<bool>> sample_read_quality_internal(pair<uint8_t, bool> first,
                                                            bool transitions_1);
    
    /// Internal method called by paired and unpaired samplers for both whole-
    /// graph and path sources. Offset and is_reverse are only used (and drive
    /// the iteration and update of curr_pos) in path node. Otherwise, in whole
    /// graph mode, they are ignored and curr_pos is used to traverse the graph
    /// directly.
    void sample_read_internal(Alignment& aln, size_t& offset, bool& is_reverse, pos_t& curr_pos,
                              const string& source_path);
    
    /// Sample an appropriate starting position according to the mode. Updates the arguments.
    void sample_start_pos(size_t& offset, bool& is_reverse, pos_t& pos, string& source_path);
    
    /// Get a random position in the graph
    pos_t sample_start_graph_pos();
    /// Get a random position along the source path
    tuple<size_t, bool, pos_t, string> sample_start_path_pos();
    
    /// Get an unclashing read name
    string get_read_name();
    
    /// Move forward one position in either the source path or the graph,
    /// depending on mode. Update the arguments. Return true if we can't because
    /// we hit a tip or false otherwise
    bool advance(size_t& offset, bool& is_reverse, pos_t& pos, char& graph_char, const string& source_path);
    /// Move forward a certain distance in either the source path or the graph,
    /// depending on mode. Update the arguments. Return true if we can't because
    /// we hit a tip or false otherwise
    bool advance_by_distance(size_t& offset, bool& is_reverse, pos_t& pos, size_t distance,
                             const string& source_path);
    
    /// Move forward one position in the source path, return true if we can't
    /// because we hit a tip or false otherwise
    bool advance_on_path(size_t& offset, bool& is_reverse, pos_t& pos, char& graph_char,
                         const string& source_path);
    /// Move forward a certain distance in the source path, return true if we
    /// can't because we hit a tip or false otherwise
    bool advance_on_path_by_distance(size_t& offset, bool& is_reverse, pos_t& pos, size_t distance,
        const string& source_path);
    
    /// Move forward one position in the graph along a random path, return true if we can't
    /// because we hit a tip or false otherwise
    bool advance_on_graph(pos_t& pos, char& graph_char);
    /// Move forward a certain distance in the graph along a random path, return true if we
    /// can't because we hit a tip or false otherwise
    bool advance_on_graph_by_distance(pos_t& pos, size_t distance);
    
    /// Mask out bases with 'N's if the mask is true
    void apply_N_mask(string& sequence, const vector<bool>& n_mask);
    
    /// Returns the position a given distance from the end of the path, walking backwards
    pos_t walk_backwards(const Path& path, size_t distance);
    /// Add a deletion to the alignment
    void apply_deletion(Alignment& aln, const pos_t& pos);
    /// Add an insertion to the alignment
    void apply_insertion(Alignment& aln, const pos_t& pos);
    /// Add a match/mismatch to the alignment
    void apply_aligned_base(Alignment& aln, const pos_t& pos, char graph_char, char read_char);
    
    /// Memo for Phred -> probability conversion
    vector<double> phred_prob;
    
    /// A Markov distribution for each read position indicating quality and whether the base is an 'N'
    vector<MarkovDistribution<pair<uint8_t, bool>, pair<uint8_t, bool>>> transition_distrs_1;
    /// A second set of Markov distributions for the second read in a pair
    vector<MarkovDistribution<pair<uint8_t, bool>, pair<uint8_t, bool>>> transition_distrs_2;
    /// A distribution for the joint initial qualities of a read pair
    MarkovDistribution<pair<uint8_t, bool>, pair<pair<uint8_t, bool>, pair<uint8_t, bool>>> joint_initial_distr;
    
    xg::XG& xg_index;
    
    LRUCache<id_t, Node> node_cache;
    LRUCache<id_t, vector<Edge> > edge_cache;
    
    default_random_engine prng;
    vg::discrete_distribution<> path_sampler;
    vector<vg::uniform_int_distribution<size_t>> start_pos_samplers;
    vg::uniform_int_distribution<uint8_t> strand_sampler;
    vg::uniform_int_distribution<size_t> background_sampler;
    vg::uniform_int_distribution<size_t> mut_sampler;
    vg::uniform_real_distribution<double> prob_sampler;
    vg::normal_distribution<double> insert_sampler;
    
    const double sub_poly_rate;
    const double indel_poly_rate;
    const double indel_error_prop;
    const double insert_mean;
    const double insert_sd;
    
    size_t sample_counter = 0;
    size_t seed;
    
    const bool retry_on_Ns;
    
    /// Restrict reads to just these paths (path-only mode) if nonempty.
    vector<string> source_paths;
};
    
/**
 * A finite state Markov distribution that supports sampling
 */


}

#endif
