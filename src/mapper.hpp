#ifndef MAPPER_H
#define MAPPER_H

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "vg.hpp"
#include "xg.hpp"
#include "index.hpp"
#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "alignment.hpp"
#include "path.hpp"
#include "position.hpp"
#include "lru_cache.h"
#include "json2pb.h"
#include "entropy.hpp"
#include "gssw_aligner.hpp"

namespace vg {

// uncomment to enable vg map --debug
//#define debug_mapper

using namespace std;
    
enum MappingQualityMethod { Approx, Exact, None };

class Mapper;

class MaximalExactMatch {

public:

    //const string* source;
    string::const_iterator begin;
    string::const_iterator end;
    gcsa::range_type range;
    size_t match_count;
    int fragment;
    std::vector<gcsa::node_type> nodes;
    map<string, vector<size_t> > positions;
    MaximalExactMatch(string::const_iterator b,
                      string::const_iterator e,
                      gcsa::range_type r,
                      size_t m = 0)
        : begin(b), end(e), range(r), match_count(m) { }

    // construct the sequence of the MEM; useful in debugging
    string sequence(void) const;
    // uses GCSA to get the positions matching the range
    void fill_nodes(gcsa::GCSA* gcsa);
    // uses GCSA to get a count of the number of graph nodes in our range
    void fill_match_count(gcsa::GCSA* gcsa);
    // get the length of the MEM
    int length(void) const;
    // uses an xgindex to fill out the MEM positions
    void fill_positions(Mapper* mapper);

    friend bool operator==(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend bool operator<(const MaximalExactMatch& m1, const MaximalExactMatch& m2);
    friend ostream& operator<<(ostream& out, const MaximalExactMatch& m);

    MaximalExactMatch(void) = default;                                      // Copy constructor
    MaximalExactMatch(const MaximalExactMatch&) = default;               // Copy constructor
    MaximalExactMatch(MaximalExactMatch&&) = default;                    // Move constructor
    MaximalExactMatch& operator=(const MaximalExactMatch&) & = default;  // MaximalExactMatchopy assignment operator
    MaximalExactMatch& operator=(MaximalExactMatch&&) & = default;       // Move assignment operator
    //virtual ~MaximalExactMatch() { }                     // Destructor
};


class MEMMarkovModelVertex {
public:
    MaximalExactMatch mem;
    vector<pair<MEMMarkovModelVertex*, double> > next_cost; // for forward
    vector<pair<MEMMarkovModelVertex*, double> > prev_cost; // for traceback
    double weight;
    double score;
    MEMMarkovModelVertex* prev;
    MEMMarkovModelVertex(void) = default;                                      // Copy constructor
    MEMMarkovModelVertex(const MEMMarkovModelVertex&) = default;               // Copy constructor
    MEMMarkovModelVertex(MEMMarkovModelVertex&&) = default;                    // Move constructor
    MEMMarkovModelVertex& operator=(const MEMMarkovModelVertex&) & = default;  // MEMMarkovModelVertexopy assignment operator
    MEMMarkovModelVertex& operator=(MEMMarkovModelVertex&&) & = default;       // Move assignment operator
    virtual ~MEMMarkovModelVertex() { }                     // Destructor
};

class MEMMarkovModel {
public:
    vector<MEMMarkovModelVertex> model;
    MEMMarkovModel(
        const vector<size_t>& aln_lengths,
        const vector<vector<MaximalExactMatch> >& matches,
        Mapper* mapper,
        const function<double(const MaximalExactMatch&, const MaximalExactMatch&)>& transition_weight,
        int band_width = 10);
    void score(const set<MEMMarkovModelVertex*>& exclude);
    MEMMarkovModelVertex* max_vertex(void);
    vector<vector<MaximalExactMatch> > traceback(int alt_alns, bool debug);
    void display(ostream& out);
    void clear_scores(void);
};


class Mapper {


private:

    // Private constructor to delegate everything to. It might have all these
    // indexing structures null, for example if being called from the default
    // constructor.
    Mapper(Index* idex, xg::XG* xidex, gcsa::GCSA* g, gcsa::LCPArray* a);
    
    Alignment align_to_graph(const Alignment& aln,
                             VG& vg,
                             size_t max_query_graph_ratio,
                             bool pinned_alignment = false,
                             bool pin_left = false,
                             int8_t full_length_bonus = 0,
                             bool global = false);
    vector<Alignment> align_multi_internal(bool compute_unpaired_qualities,
                                           const Alignment& aln,
                                           int kmer_size,
                                           int stride,
                                           int max_mem_length,
                                           int band_width,
                                           double& cluster_mq,
                                           int additional_multimaps = 0,
                                           vector<MaximalExactMatch>* restricted_mems = nullptr);
    void compute_mapping_qualities(vector<Alignment>& alns, double cluster_mq);
    void compute_mapping_qualities(pair<vector<Alignment>, vector<Alignment>>& pair_alns, double cluster_mq);
    vector<Alignment> score_sort_and_deduplicate_alignments(vector<Alignment>& all_alns, const Alignment& original_alignment);
    void filter_and_process_multimaps(vector<Alignment>& all_alns, int additional_multimaps);
    vector<Alignment> align_multi_kmers(const Alignment& aln, int kmer_size = 0, int stride = 0, int band_width = 1000);
    // Return the one best banded alignment.
    Alignment align_banded(const Alignment& read,
                           int kmer_size = 0,
                           int stride = 0,
                           int max_mem_length = 0,
                           int band_width = 1000);
    // alignment based on the MEM approach
    vector<Alignment> align_mem_multi(const Alignment& alignment, vector<MaximalExactMatch>& mems, double& cluster_mq, int additional_multimaps = 0);
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
    // GSSW aligner(s)
    vector<QualAdjAligner*> qual_adj_aligners;
    vector<Aligner*> regular_aligners;
    void clear_aligners(void);
    QualAdjAligner* get_qual_adj_aligner(void);
    Aligner* get_regular_aligner(void);

    // match walking support to prevent repeated calls to the xg index for the same node
    vector<LRUCache<id_t, Node>* > node_cache;
    LRUCache<id_t, Node>& get_node_cache(void);
    void init_node_cache(void);

    // match node traversals to path positions
    vector<LRUCache<gcsa::node_type, map<string, vector<size_t> > >* > node_pos_cache;
    LRUCache<gcsa::node_type, map<string, vector<size_t> > >& get_node_pos_cache(void);
    void init_node_pos_cache(void);
    map<string, vector<size_t> > node_positions_in_paths(gcsa::node_type node);

    // a collection of read pairs which we'd like to realign once we have estimated the fragment_size
    vector<pair<Alignment, Alignment> > imperfect_pairs_to_retry;

    // running estimation of fragment length distribution
    deque<double> fragment_lengths;
    void record_fragment_length(int length);
    double fragment_length_stdev(void);
    double fragment_length_mean(void);
    double fragment_length_pdf(double length);
    double cached_fragment_length_mean;
    double cached_fragment_length_stdev;
    int since_last_fragment_length_estimate;
    int fragment_length_estimate_interval;

    double estimate_gc_content();
    void init_aligner(int32_t match, int32_t mismatch, int32_t gap_open, int32_t gap_extend);
    void set_alignment_scores(int32_t match, int32_t mismatch, int32_t gap_open, int32_t gap_extend);

    // use the xg index to get the mean position of the nodes in the alignent for each reference that it corresponds to
    map<string, double> alignment_mean_path_positions(const Alignment& aln, bool first_hit_only = true);

    // Return true of the two alignments are consistent for paired reads, and false otherwise
    bool alignments_consistent(const map<string, double>& pos1,
                               const map<string, double>& pos2,
                               int fragment_size_bound);

    // Align read2 to the subgraph near the alignment of read1.
    // TODO: support banded alignment and intelligently use orientation heuristics
    void align_mate_in_window(const Alignment& read1, Alignment& read2, int pair_window);
    
    vector<Alignment> resolve_banded_multi(vector<vector<Alignment>>& multi_alns);
    set<MaximalExactMatch*> resolve_paired_mems(vector<MaximalExactMatch>& mems1,
                                                vector<MaximalExactMatch>& mems2);

    // uses heuristic clustering based on node id ranges to find alignment targets, and aligns
    vector<Alignment> mems_id_clusters_to_alignments(const Alignment& alignment, vector<MaximalExactMatch>& mems, int additional_multimaps);

    // uses approximate-positional clustering based on embedded paths in the xg index to find and align against alignment targets
    vector<Alignment> mems_pos_clusters_to_alignments(const Alignment& aln, vector<MaximalExactMatch>& mems, int additional_multimaps, double& cluster_mq);
    // helper for computing the number of bases in the query covered by a cluster
    int cluster_coverage(const vector<MaximalExactMatch>& cluster);

    // takes the input alignment (with seq, etc) so we have reference to the base sequence
    // for reconstruction the alignments from the SMEMs
    Alignment mems_to_alignment(const Alignment& aln, vector<MaximalExactMatch>& mems);
    Alignment mem_to_alignment(MaximalExactMatch& mem);
    // fix up a SMEM-threaded exact match alignment by locally aligning small pieces against gaps in alignment
    Alignment patch_alignment(const Alignment& aln);
    // try to locally align portions of the alignment with "knots" where we slip backwards in the graph
    Alignment smooth_alignment(const Alignment& aln);
    // use the scoring provided by the internal aligner to re-score the alignment, scoring gaps using graph distance
    int32_t score_alignment(const Alignment& aln);

    bool adjacent_positions(const Position& pos1, const Position& pos2);
    int64_t get_node_length(int64_t node_id);
    bool check_alignment(const Alignment& aln);
    VG alignment_subgraph(const Alignment& aln, int context_size = 1);
    
    // Align the given string and return an Alignment.
    Alignment align(const string& seq,
                    int kmer_size = 0,
                    int stride = 0,
                    int max_mem_length = 0,
                    int band_width = 1000);

    // Align the given read and return an aligned copy. Does not modify the input Alignment.
    Alignment align(const Alignment& read,
                    int kmer_size = 0,
                    int stride = 0,
                    int max_mem_length = 0,
                    int band_width = 1000);

    // Align the given read with multi-mapping. Returns the alignments in score
    // order, up to multimaps (or max_multimaps if multimaps is 0). Does not update the alignment passed in.
    // If the sequence is longer than the band_width, will only produce a single best banded alignment.
    // All alignments but the first are marked as secondary.
    vector<Alignment> align_multi(const Alignment& aln,
                                  int kmer_size = 0,
                                  int stride = 0,
                                  int max_mem_length = 0,
                                  int band_width = 1000);
    
    // paired-end based
    
    // Both vectors of alignments will be sorted in order of increasing score.
    // All alignments but the first in each vector are marked as secondary.
    // Alignments at corresponding positions in the two vectors may or may not
    // be corresponding paired alignments. If a read does not map, its vector
    // will be empty.
    // If only_top_scoring_pair is set, then the vectors will be empty unless
    // the primary pair of alignments each have top scores individually as well. 
    pair<vector<Alignment>, vector<Alignment>> 
        align_paired_multi(const Alignment& read1,
                           const Alignment& read2,
                           bool& queued_resolve_later,
                           int kmer_size = 0,
                           int stride = 0,
                           int max_mem_length = 0,
                           int band_width = 1000,
                           int pair_window = 64,
                           bool only_top_scoring_pair = false);

    // align each fragment separately, then find consistent results using various heuristics
    pair<vector<Alignment>, vector<Alignment>> 
        align_paired_multi_sep(const Alignment& read1,
                               const Alignment& read2,
                               bool& queued_resolve_later,
                               int kmer_size = 0,
                               int stride = 0,
                               int max_mem_length = 0,
                               int band_width = 1000,
                               int pair_window = 64,
                               bool only_top_scoring_pair = false);

    // align the pair as a single component using MEM threading and patching on the pair simultaneously
    pair<vector<Alignment>, vector<Alignment>> 
        align_paired_multi_simul(const Alignment& read1,
                                 const Alignment& read2,
                                 bool& queued_resolve_later,
                                 int max_mem_length = 0,
                                 bool only_top_scoring_pair = false);

    // lossily project an alignment into a particular path space of a graph
    // the resulting alignment is equivalent to a SAM record against the chosen path
    Alignment surject_alignment(const Alignment& source,
                                set<string>& path_names,
                                string& path_name,
                                int64_t& path_pos,
                                bool& path_reverse,
                                int window);

    // MEM-based mapping
    // find maximal exact matches
    // These are SMEMs by definition when shorter than the max_mem_length or GCSA2 order.
    // Enables reseeding of long matches when reseed_length > 0.
    // Reseeding requires an occurrence count for the MEM that we are reseeding.
    // Where k is the number of hits for the MEM we are reseeding,
    // we return the MEMs shorter than the MEM to reseed that have count > k.
    vector<MaximalExactMatch> find_mems(string::const_iterator seq_begin,
                                        string::const_iterator seq_end,
                                        int max_mem_length,
                                        int reseed_length = 0);
    bool get_mem_hits_if_under_max(MaximalExactMatch& mem);
    // debugging, checking of mems using find interface to gcsa
    void check_mems(const vector<MaximalExactMatch>& mems);
    // compute a mapping quality component based only on the MEMs we've obtained
    double compute_cluster_mapping_quality(const vector<vector<MaximalExactMatch> >& clusters, int read_length);
    // use BFS to expand the graph in an attempt to resolve soft clips
    void resolve_softclips(Alignment& aln, VG& graph);
    // walks the graph one base at a time from pos1 until we find pos2
    int graph_distance(pos_t pos1, pos_t pos2, int maximum = 1e3);
    // use the xg index to get a character at a particular position (rc or foward)
    char pos_char(pos_t pos);
    // the next positions and their characters following the same strand of the graph
    map<pos_t, char> next_pos_chars(pos_t pos);
    // get the positions some specific distance from the given position (in the forward direction)
    set<pos_t> positions_bp_from(pos_t pos, int distance, bool rev);
    // convert a single MEM hit into an alignment (by definition, a perfect one)
    Alignment walk_match(const string& seq, pos_t pos);
    vector<Alignment> walk_match(const Alignment& base, const string& seq, pos_t pos);
    // convert the set of hits of a MEM into a set of alignments
    vector<Alignment> mem_to_alignments(MaximalExactMatch& mem);
    // Use the GCSA index to look up the sequence
    set<pos_t> sequence_positions(const string& seq);

    // fargment length estimation
    map<string, int> approx_pair_fragment_length(const Alignment& aln1, const Alignment& aln2);
    // uses the cached information about the graph in the xg index to get an approximate node length
    double average_node_length(void);
    
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
    //int max_mem_length; // a mem must be <= this length
    int min_mem_length; // a mem must be >= this length
    bool mem_threading; // whether to use the mem threading mapper or not
    int mem_reseed_length; // the length above which we reseed MEMs to get potentially missed hits

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
    int max_softclip_iterations; // Extend no more than this many times (while softclips are getting shorter)
    float min_identity; // require that alignment identity is at least this much to accept alignment
    // paired-end consistency enforcement
    int extra_multimaps; // Extra mappings considered
    
    bool adjust_alignments_for_base_quality; // use base quality adjusted alignments
    MappingQualityMethod mapping_quality_method; // how to compute mapping qualities
    int max_mapping_quality; // the cap for mapping quality

    bool always_rescue; // Should rescue be attempted for all imperfect alignments?
    int fragment_max; // the maximum length fragment which we will consider when estimating fragment lengths
    int fragment_size; // Used to bound clustering of MEMs during paired end mapping, also acts as sentinel to determine
                       // if consistent pairs should be reported; dynamically estimated at runtime
    double fragment_sigma; // the number of times the standard deviation above the mean to set the fragment_size
    int fragment_length_cache_size;
    float perfect_pair_identity_threshold;
    int8_t full_length_alignment_bonus;

};

// utility
const vector<string> balanced_kmers(const string& seq, int kmer_size, int stride);
const string mems_to_json(const vector<MaximalExactMatch>& mems);
set<pos_t> gcsa_nodes_to_positions(const vector<gcsa::node_type>& nodes);

}

#endif
