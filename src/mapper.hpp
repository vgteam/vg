#ifndef MAPPER_H
#define MAPPER_H

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "xg.hpp"
#include "index.hpp"
#include "gcsa.h"
#include "lcp.h"
#include "alignment.hpp"
#include "path.hpp"
#include "json2pb.h"
#include "entropy.hpp"
#include "gssw_aligner.hpp"

namespace vg {

using namespace std;

class MaximalExactMatch {

public:

    string::const_iterator begin;
    string::const_iterator end;
    gcsa::range_type range;
    size_t match_count;
    std::vector<gcsa::node_type> nodes;
    MaximalExactMatch(string::const_iterator b,
                      string::const_iterator e,
                      gcsa::range_type r,
                      size_t m = 0)
        : begin(b), end(e), range(r), match_count(m) { }

    // construct the sequence of the MEM; useful in debugging
    string sequence(void) const {
        string seq; //seq.resize(end-begin);
        string::const_iterator c = begin;
        while (c != end) seq += *c++;
        return seq;
    }
    // uses GCSA to get the positions matching the range
    void fill_nodes(gcsa::GCSA* gcsa) {
        gcsa->locate(range, nodes);
    }
    // uses GCSA to get a count of the number of graph nodes in our range
    void fill_match_count(gcsa::GCSA* gcsa) {
        match_count = gcsa->count(range);
    }
};

class Mapper {


private:

    // Private constructor to delegate everything to. It might have all these
    // indexing structures null, for example if being called from the default
    // constructor.
    Mapper(Index* idex, xg::XG* xidex, gcsa::GCSA* g, gcsa::LCPArray* a);
    
    bool adjust_alignments_for_base_quality;

public:
    // Make a Mapper that pulls from a RocksDB index and optionally a GCSA2 kmer index.
    Mapper(Index* idex, gcsa::GCSA* g = nullptr, gcsa::LCPArray* a = nullptr);
    // Make a Mapper that pulls from an XG succinct graph and a GCSA2 kmer index + LCP array
    Mapper(xg::XG* xidex, gcsa::GCSA* g, gcsa::LCPArray* a);
    Mapper(void);
    ~Mapper(void);
    // rocksdb index
    Index* index;
    // xg index
    xg::XG* xindex;
    // GCSA index and its LCP array
    gcsa::GCSA* gcsa;
    gcsa::LCPArray* lcp;
    // GSSW aligner
    Aligner aligner;
    
    double estimate_gc_content();
    void init_aligner(int32_t match, int32_t mismatch, int32_t gap_open, int32_t gap_extend);
    void set_alignment_scores(int32_t match, int32_t mismatch, int32_t gap_open, int32_t gap_extend);
    void set_base_quality_adjusted_alignments(bool do_adjustments);

    // Align the given string and return an Alignment.
    Alignment align(const string& seq, int kmer_size = 0, int stride = 0, int band_width = 1000);
    // Align the given read and return an aligned copy. Does not modify the input Alignment.
    Alignment align(const Alignment& read, int kmer_size = 0, int stride = 0, int band_width = 1000);
    // Align the given read with multi-mapping. Returns the alignments in score
    // order, up to multimaps (or max_multimaps if multimaps is 0). Does not update the alignment passed in.
    // If the sequence is longer than the band_width, will only produce a single best banded alignment.
    // All alignments but the first are marked as secondary.
    vector<Alignment> align_multi(const Alignment& aln, int kmer_size = 0, int stride = 0, int band_width = 1000, int multimaps = 0);
    vector<Alignment> align_multi_kmers(const Alignment& aln, int kmer_size = 0, int stride = 0, int band_width = 1000, int multimaps = 0);

    // Return true of the two alignments are consistent for paired reads, and false otherwise
    bool alignments_consistent(const Alignment& aln1, const Alignment& aln2, int pair_window);

    // Align read2 to the subgraph near the alignment of read1.
    // TODO: support banded alignment and intelligently use orientation heuristics
    void align_mate_in_window(const Alignment& read1, Alignment& read2, int pair_window);

    // Return the one best banded alignment.
    Alignment align_banded(const Alignment& read,
                           int kmer_size = 0,
                           int stride = 0,
                           int band_width = 1000);
    
    vector<Alignment> resolve_banded_multi(vector<vector<Alignment>>& multi_alns);
    bool adjacent_positions(const Position& pos1, const Position& pos2);
    int64_t get_node_length(int64_t node_id);

    // paired-end based
    
    // Both vectors of alignments will be sorted in order of increasing score.
    // All alignments but the first in each vector are marked as secondary.
    // Alignments at corresponding positions in the two vectors may or may not
    // be corresponding paired alignments. If a read does not map, its vector
    // will be empty.
    pair<vector<Alignment>, vector<Alignment>> 
        align_paired_multi(const Alignment& read1,
                           const Alignment& read2,
                           int kmer_size = 0,
                           int stride = 0,
                           int band_width = 1000,
                           int pair_window = 64);
    
    // Paired-end alignment ignoring multi-mapping. Returns either the two
    // highest-scoring reads if no rescue was required, or the highest-scoring
    // read and its corresponding rescue result if rescue was used.
    pair<Alignment, Alignment> align_paired(const Alignment& read1,
                                            const Alignment& read2,
                                            int kmer_size = 0,
                                            int stride = 0,
                                            int band_width = 1000,
                                            int pair_window = 64);

    // base algorithm for above Update the passed-in Alignment with a highest-
    // score alignment, and return all good alignments sorted by score up to
    // max_multimaps. If the read does not map, the returned vector will be
    // empty. No alignments will be marked as secondary; the caller must do that
    // if they plan to produce GAM output.
    vector<Alignment> align_threaded(const Alignment& read,
                                     int& hit_count,
                                     int kmer_size = 0,
                                     int stride = 0,
                                     int attempt = 0);

    // MEM-based mapping
    // finds absolute super-maximal exact matches
    vector<MaximalExactMatch> find_smems(const string& seq);
    bool get_mem_hits_if_under_max(MaximalExactMatch& mem);
    // debugging, checking of mems using find interface to gcsa
    void check_mems(const vector<MaximalExactMatch>& mems);
    // finds "forward" maximal exact matches of the sequence using the GCSA index
    // stepping step between each one
    vector<MaximalExactMatch> find_forward_mems(const string& seq, size_t step = 1);
    // alignment based on the MEM approach
    vector<Alignment> align_mem(const Alignment& alignment, int multimaps = 0);
    Alignment align_mem_optimal(const Alignment& alignment, vector<MaximalExactMatch>& mems);
    vector<Alignment> align_mem_multi(const Alignment& alignment, vector<MaximalExactMatch>& mems, int max_multi = 0);
    // use BFS to expand the graph in an attempt to resolve soft clips
    void resolve_softclips(Alignment& aln, VG& graph);
    // use the xg index to get a character at a particular position (rc or foward)
    char pos_char(pos_t pos);
    // the next positions and their characters following the same strand of the graph
    map<pos_t, char> next_pos_chars(pos_t pos);
    // convert a single MEM hit into an alignment (by definition, a perfect one)
    Alignment walk_match(const string& seq, pos_t pos);
    // convert the set of hits of a MEM into a set of alignments
    vector<Alignment> mem_to_alignments(MaximalExactMatch& mem);
    // Use the GCSA index to look up the sequence
    set<pos_t> sequence_positions(const string& seq);
    
    bool debug;
    int alignment_threads; // how many threads will *this* mapper use when running banded alignmentsx

    // kmer/"threaded" mapper parameters
    //
    set<int> kmer_sizes; // taken from rocksdb index
    int best_clusters; // use up to this many clusters to build threads
    int cluster_min; // minimum number of hits nearby before we test local alignment
    int hit_size_threshold; // This is in bytes. TODO: Make it not in bytes, guessing at rocksdb records per byte.
    float min_kmer_entropy; // exclude kmers with less that this entropy/base
    int kmer_min; // don't decrease kmer size below this level when trying shorter kmers
    int max_thread_gap; // maximum number of nodes in id space to extend a thread (assumes semi partial order on graph ids)
    int kmer_sensitivity_step; // size to decrease the kmer length if we fail alignment
    bool prefer_forward; // attempt alignment of forward complement of the read against the graph (forward) first
    bool greedy_accept; // if we make an OK forward alignment, accept it
    float accept_identity; // for early bailout; target alignment score as a fraction of the score of a perfect match

    // mem mapper parameters (it is _much_ simpler)
    //
    int max_mem_length; // a mem must be <= this length
    int min_mem_length; // a mem must be >= this length

    // general parameters, applying to both types of mapping
    //
    int hit_max;       // ignore kmers or MEMs (TODO) with more than this many hits
    int context_depth; // how deeply the mapper will extend out the subgraph prior to alignment
    int max_attempts;  // maximum number of times to try to increase sensitivity or use a lower-hit subgraph
    int thread_extension; // add this many nodes in id space to the end of the thread when building thread into a subgraph
    int max_target_factor; // the maximum multiple of the read length we'll try to align to

    size_t max_query_graph_ratio;

    // multimapping
    int max_multimaps;
    // soft clip resolution
    int softclip_threshold; // if more than this many bp are clipped, try extension algorithm
    float min_identity; // require that alignment identity is at least this much to accept alignment
    // paired-end consistency enforcement
    bool promote_consistent_pairs; // Should consistent paired mappings be made primary over higher-scoring inconsistent ones?
    int extra_pairing_multimaps; // Extra mappings considered for finding consistent paired-end mappings

};

// utility
const vector<string> balanced_kmers(const string& seq, int kmer_size, int stride);
const string mems_to_json(const vector<MaximalExactMatch>& mems);
set<pos_t> gcsa_nodes_to_positions(const vector<gcsa::node_type>& nodes);

}

#endif
