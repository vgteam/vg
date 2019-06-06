#ifndef VG_PACKER_HPP_INCLUDED
#define VG_PACKER_HPP_INCLUDED

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "omp.h"
#include "lru_cache.h"
#include "xg.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "position.hpp"
#include "json2pb.h"
#include "graph.hpp"
#include "gcsa/internal.h"
#include "xg_position.hpp"
#include "utility.hpp"

namespace vg {

using namespace sdsl;

class Packer {
public:
    Packer(void);
    Packer(xg::XG* xidx, size_t bin_size = 0, bool qual_adjust = false);
    ~Packer(void);
    xg::XG* xgidx;
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
    void add(const Alignment& aln, bool record_edits = true);
    size_t graph_length(void) const;
    size_t position_in_basis(const Position& pos) const;
    string pos_key(size_t i) const;
    string edit_value(const Edit& edit, bool revcomp) const;
    vector<Edit> edits_at_position(size_t i) const;
    size_t coverage_at_position(size_t i) const;
    void collect_coverage(const Packer& c);
    ostream& as_table(ostream& out, bool show_edits, vector<vg::id_t> node_ids);
    ostream& as_edge_table(ostream& out, vector<vg::id_t> node_ids);
    ostream& show_structure(ostream& out); // debugging
    void write_edits(vector<ofstream*>& out) const; // for merge
    void write_edits(ostream& out, size_t bin) const; // for merge
    size_t get_bin_size(void) const;
    size_t get_n_bins(void) const;
    bool is_dynamic(void);
    size_t coverage_size(void);

    size_t edge_coverage(Edge& e) const;
    size_t edge_coverage(size_t i) const;
    size_t edge_count(void) const;
    size_t edge_vector_size(void) const;
private:
    void ensure_edit_tmpfiles_open(void);
    void close_edit_tmpfiles(void);
    void remove_edit_tmpfiles(void);
    bool is_compacted = false;
    // dynamic model
    gcsa::CounterArray coverage_dynamic;
    gcsa::CounterArray edge_coverage_dynamic;
    vector<string> edit_tmpfile_names;
    vector<ofstream*> tmpfstreams;
    // which bin should we use
    size_t bin_for_position(size_t i) const;
    size_t n_bins = 1;
    size_t bin_size = 0;
    size_t edit_length = 0;
    size_t edit_count = 0;
    dac_vector<> coverage_civ; // graph coverage (compacted coverage_dynamic)
    vlc_vector<> edge_coverage_civ; // edge coverage (compacted edge_coverage_dynamic)
    //
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

    // toggle quality adjusted mode
    bool qual_adjust;
    
    // Combine the MAPQ and base quality (if available) for a given position in the read
    int compute_quality(const Alignment& aln, size_t position_in_read) const;
    int combine_qualities(int map_quality, int base_quality) const;
    
    // Avoid recomputing qualities in above
    mutable LRUCache<pair<int, int>, int>* quality_cache;
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
