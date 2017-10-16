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
#include "lru_cache.h"
#include "json2pb.h"

namespace vg {

using namespace std;

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
    // A string which, if nonempty, gives the name of the path to restrict simulated reads to.
    string source_path;
    inline Sampler(xg::XG* x,
            int seed = 0,
            bool forward_only = false,
            bool allow_Ns = false,
            const string& source_path = "")
        : xgidx(x),
          node_cache(100),
          edge_cache(100),
          forward_only(forward_only),
          no_Ns(!allow_Ns),
          nonce(0),
          source_path(source_path) {
        if (!seed) {
            seed = time(NULL);
        }
        rng.seed(seed);
    }

    pos_t position(void);
    string sequence(size_t length);
    
    /// Get an alignment against the whole graph, or against the source path if
    /// one is selected.
    Alignment alignment(size_t length);
    
    /// Get an alignment against the whole graph.
    Alignment alignment_to_graph(size_t length);
    
    /// Get an alignment against the currently set source_path.
    Alignment alignment_to_path(size_t length);
    
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
                             uniform_real_distribution<double>& rprob,
                             uniform_int_distribution<int>& rbase);

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
    /// Every read must be the same length. Polymorphism rates apply uniformly along a
    /// read, whereas errors are distributed as indicated by the learned distribution.
    NGSSimulator(xg::XG& xg_index,
                 const string& ngs_fastq_file,
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
    class MarkovDistribution;
    
    NGSSimulator(void) = delete;
    
    /// DNA alphabet
    static const string alphabet;
    /// Remainder of the alphabet after removing a given character
    unordered_map<char, string> mutation_alphabets;
    
    /// Add a quality string to the training data
    void record_read_quality(const Alignment& aln);
    /// Indicate that there is no more training data
    void finalize();
    /// Get a quality string that mimics the training data
    string sample_read_quality();
    
    /// Internal method called by paired and unpaired samplers
    void sample_read_internal(Alignment& aln, pos_t& curr_pos);
    /// Get a random position in the graph
    pos_t sample_start_pos();
    /// Get an unclashing read name
    string get_read_name();
    
    /// Move forward one position in the graph along a random path, return true if we can't
    /// because we hit a tip or false otherwise
    bool advance_on_graph(pos_t& pos, char& graph_char);
    /// Move forward a certain distance in the graph along a random path, return true if we
    /// can't because we hit a tip or false otherwise
    bool advance_on_graph_by_distance(pos_t& pos, size_t distance);
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
    
    /// A Markov distribution for each read position
    vector<MarkovDistribution> transition_distrs;
    
    xg::XG& xg_index;
    
    LRUCache<id_t, Node> node_cache;
    LRUCache<id_t, vector<Edge> > edge_cache;
    
    default_random_engine prng;
    uniform_int_distribution<size_t> start_pos_sampler;
    uniform_int_distribution<uint8_t> strand_sampler;
    uniform_int_distribution<size_t> background_sampler;
    uniform_int_distribution<size_t> mut_sampler;
    uniform_real_distribution<double> prob_sampler;
    normal_distribution<double> insert_sampler;
    
    const double sub_poly_rate;
    const double indel_poly_rate;
    const double indel_error_prop;
    const double insert_mean;
    const double insert_sd;
    
    size_t sample_counter = 0;
    size_t seed;
    
    const bool retry_on_Ns;
};
    
/**
 * A finite state Markov distribution that supports sampling
 */
class NGSSimulator::MarkovDistribution {
public:
    MarkovDistribution(size_t seed);
    
    /// record a transition from the input data
    void record_transition(uint8_t from, uint8_t to);
    /// indicate that there is no more data and prepare for sampling
    void finalize();
    /// sample according to the training data
    uint8_t sample_transition(uint8_t from);
    
private:
    
    default_random_engine prng;
    unordered_map<uint8_t, uniform_int_distribution<size_t>> samplers;
    
    unordered_map<uint8_t, size_t> column_of;
    vector<uint8_t> value_at;
    unordered_map<uint8_t, vector<size_t>> cond_distrs;
    
};

}

#endif
