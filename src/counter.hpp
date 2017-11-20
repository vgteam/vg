#ifndef VG_MAPPER_HPP_INCLUDED
#define VG_MAPPER_HPP_INCLUDED

#include <iostream>
#include <map>
#include <chrono>
#include <ctime>
#include "omp.h"
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

class Counter {
public:
    Counter(void);
    Counter(xg::XG* xidx);
    ~Counter(void);
    xg::XG* xgidx;
    //void load(const vector<string>& file_names);
    void load_from_file(const string& file_name);
    void save_to_file(const string& file_name);
    void load(istream& in);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* s = NULL,
                     std::string name = "");
    void ensure_edit_tmpfile_open(void);
    void close_edit_tmpfile(void);
    void remove_edit_tmpfile(void);
    void make_compact(void);
    void make_dynamic(void);
    void add(const Alignment& aln, bool record_edits = true);
    size_t position_in_basis(const Position& pos);
    string pos_key(size_t i);
    vector<Edit> edits_at_position(size_t i);
    ostream& as_table(ostream& out);
    ostream& show_structure(ostream& out); // debugging
private:
    bool is_compacted;
    // dynamic model
    gcsa::CounterArray coverage_dynamic;
    string edit_tmpfile_name;
    fstream tmpfstream;
    //map<size_t, map<string, int32_t> > edits;
    // compact model
    //size_t total_count; // sum of counts in coverage and edit coverage
    size_t edit_length;
    size_t edit_count;
    vlc_vector<> coverage_civ; // graph coverage (compacted coverage_dynamic)
    csa_sada<enc_vector<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<> > edit_csa;
    // make separators that are like improbable varints
    // so that they are unlikely to occur in structured data we write
    string sep_start = string(2, '\xff');
    string sep_end = string(2, '\xff');
};

// for making a combined matrix output and maybe doing other fun operations
class Counters : public vector<Counter> {
    void load(const vector<string>& file_names);
    ostream& as_table(ostream& out);
};

}

#endif
