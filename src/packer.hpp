#ifndef VG_PACKER_HPP_INCLUDED
#define VG_PACKER_HPP_INCLUDED

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include <mutex>
#include "omp.h"
#include "lru_cache.h"
#include "alignment.hpp"
#include "path.hpp"
#include "position.hpp"
#include "vg/io/json2pb.h"
#include "graph.hpp"
#include "gcsa/internal.h"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "utility.hpp"
#include <bdsg/overlays/overlay_helper.hpp>

namespace vg {

using namespace sdsl;

/// Packer collects coverage of a GAM using compressed indexes
/// Any combination of these 3 types of information can be stored
/// - base coverage : number of reads aligning to a given base (offset in node) in the graph
/// - edge coverage : number of reads aligning to a given edge in the graph
/// - edits : a list of edits at a given base in the graph
/// In memory, the coverages are stored in SDSL int vectors (dynamic) and on disk they are compressed int vectors
class Packer {
public:
    
    /// Some helper functions to heuristically estimate input parameters for constructor
    static size_t estimate_data_width(size_t expected_coverage);
    static size_t estimate_batch_size(size_t num_threads);
    static size_t estimate_bin_count(size_t num_threads);

    /// Create a Packer (to read from a file)
    Packer(const HandleGraph* base_graph = nullptr);
    
    /// Create a Packer (to write to)
    /// graph : Must implement the VectorizableHandleGraph interface
    /// record_bases : Store the base coverage
    /// record_edges : Store the edge coverage
    /// record_edits : Store the edits
    /// record_qualities : Store the average MAPQ for each node rank    
    /// bin_size : Bin coverage into bins
    /// coverage_bins : Use this many coverage objects.  Using one / thread allows faster merge
    /// coverage_locks : Number of mutexes to use for each of node and edge coverage.
    /// data_width : Number of bits per entry in the dynamic coverage vector.  Higher values get stored in a map
    Packer(const HandleGraph* base_graph, bool record_bases, bool record_edges, bool record_edits, bool record_qualities,
           size_t bin_size = 0, size_t coverage_bins = 1, size_t data_width = 8);
    ~Packer();
    void clear();

    /// Add coverage from given alignment to the indexes
    /// aln : given alignemnt
    /// min_mapq : ignore alignments with mapping_quality below this value
    /// min_baseq : ignore bases in the alignment if their read quality is below this value
    /// trim_ends : ignore first and last <trim_ends> bases
    void add(const Alignment& aln, int min_mapq = 0, int min_baseq = 0, int trim_ends = 0);

    void merge_from_files(const vector<string>& file_names);
    void merge_from_dynamic(vector<Packer*>& packers);
    void load_from_file(const string& file_name);
    void save_to_file(const string& file_name);
    void load(istream& in);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* s = NULL,
                     std::string name = "");
    void make_compact(void);
    void make_dynamic(void);
    size_t position_in_basis(const Position& pos) const;
    string pos_key(size_t i) const;
    string edit_value(const Edit& edit, bool revcomp) const;
    vector<Edit> edits_at_position(size_t i) const;
    size_t coverage_at_position(size_t i) const;
    void collect_coverage(const vector<Packer*>& packers);
    ostream& as_table(ostream& out, bool show_edits, vector<vg::id_t> node_ids);
    ostream& as_edge_table(ostream& out, vector<vg::id_t> node_ids);
    ostream& as_quality_table(ostream& out, vector<vg::id_t> node_ids);
    ostream& show_structure(ostream& out); // debugging
    void write_edits(vector<ofstream*>& out) const; // for merge
    void write_edits(ostream& out, size_t bin) const; // for merge
    size_t get_bin_size(void) const;
    size_t get_n_bins(void) const;
    bool is_dynamic(void) const;
    const HandleGraph* get_graph() const;
    const VectorizableHandleGraph* get_vec_graph() const;
    size_t coverage_size(void) const ;
    void increment_coverage(size_t i);
    void increment_coverage(size_t i, size_t v);

    size_t edge_coverage(Edge& e) const;
    size_t edge_coverage(size_t i) const;
    size_t edge_vector_size(void) const;
    /// Return the 1-based index of the given edge in our vectorization order,
    /// or 0 if the edge does not exist in the graph.
    size_t edge_index(const Edge& e) const;
    void increment_edge_coverage(size_t i);
    void increment_edge_coverage(size_t i, size_t v);

    /// total node quality (faster from dynamimc)
    size_t total_node_quality(size_t i) const;
    /// average node quality (faster from static)
    size_t average_node_quality(size_t i) const;
    size_t node_quality_vector_size(void) const;
    /// Return the 1-based node rank or 0 if node not in graph
    size_t node_index(nid_t node_id) const;
    /// and back
    nid_t index_to_node(size_t i) const;
    void increment_node_quality(size_t i, size_t v);
    /// return true if there's at least one nonzero quality in the structure
    bool has_qualities() const;
    
private:
    /// map from absolute postion to positions in the binned arrays
    pair<size_t, size_t> coverage_bin_offset(size_t i) const;
    pair<size_t, size_t> edge_coverage_bin_offset(size_t i) const;
    pair<size_t, size_t> node_quality_bin_offset(size_t i) const;
    /// get the size of a bin
    size_t coverage_bin_size(size_t i) const;
    size_t edge_coverage_bin_size(size_t i) const;
    size_t node_quality_bin_size(size_t i) const;
    /// initialize coverage bins on demand
    void init_coverage_bin(size_t i);
    void init_edge_coverage_bin(size_t i);
    void init_node_quality_bin(size_t i);
    
    void ensure_edit_tmpfiles_open(void);
    void close_edit_tmpfiles(void);
    void remove_edit_tmpfiles(void);
    bool is_compacted = false;
    
    // base graph
    const HandleGraph* base_graph;
    // vectorizable graph
    bdsg::VectorizableOverlayHelper overlay_helper;
    const VectorizableHandleGraph* graph;

    // data with for counter arrays
    size_t data_width;
    // bin sizes (last bins may be a bit bigger)
    size_t cov_bin_size;
    size_t edge_cov_bin_size;
    size_t node_qual_bin_size;

    // dynamic model
    // base coverage.  we bin to make merging faster
    vector<gcsa::CounterArray*> coverage_dynamic;
    // total length of above vectors
    size_t num_bases_dynamic;
    // one mutex per element of coverage_dynamic
    std::mutex* base_locks;
    // edge coverage.  we bin to make merging faster
    vector<gcsa::CounterArray*> edge_coverage_dynamic;
    // total length of above
    size_t num_edges_dynamic;
    // one mutex per element of edge_coverage_dynamic
    std::mutex* edge_locks;
    // node qualities
    // in the dynamic structure, we keep count of total quality for a node
    // and it gets averaged out in the static, so slightly different
    // semantics than the node and edge coverage which store the same information
    vector<gcsa::CounterArray*> node_quality_dynamic;
    // total length of above
    size_t num_nodes_dynamic;
    // one mutex per element of node_quality_dynamic
    std::mutex* node_quality_locks;
    
    vector<string> edit_tmpfile_names;
    vector<ofstream*> tmpfstreams;
    std::mutex* tmpfstream_locks;
    // which bin should we use
    size_t bin_for_position(size_t i) const;
    size_t n_bins = 1;
    size_t bin_size = 0;
    size_t edit_length = 0;
    size_t edit_count = 0;
    dac_vector<> coverage_civ; // graph coverage (compacted coverage_dynamic)
    vlc_vector<> edge_coverage_civ; // edge coverage (compacted edge_coverage_dynamic)
    vlc_vector<> node_quality_civ; // averge mapq for each node rank (compacted node_quality_dynamic)
    // edits
    vector<csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<> > > edit_csas;
    // make separators that are somewhat unusual, as we escape these
    char delim1 = '\xff';
    char delim2 = '\xfe';
    // double the delimiter in the string
    string escape_delim(const string& s, char d) const;
    string escape_delims(const string& s) const;
    // take each double delimiter back to a single
    string unescape_delim(const string& s, char d) const;
    string unescape_delims(const string& s) const;

    // toggles:
    bool record_bases;
    bool record_edges;
    bool record_edits;
    bool record_qualities;
    
    // Combine the MAPQ and base quality (if available) for a given position in the read
    int compute_quality(const Alignment& aln, size_t position_in_read) const;
    int combine_qualities(int map_quality, int base_quality) const;
    
    // Avoid recomputing qualities in above (one per thread)
    mutable vector<LRUCache<pair<int, int>, int>*> quality_cache;
    static const int maximum_quality;
    static const int lru_cache_size;
    
};

// for making a combined matrix output and maybe doing other fun operations
class Packers : public vector<Packer> {
    void load(const vector<string>& file_names);
    ostream& as_table(ostream& out);
};

}

#endif
