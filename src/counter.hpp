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

namespace vg {

class Counter {
public:
    Counter(void);
    Counter(xg::XG* xidx);
    ~Counter(void);
    xg::XG* xgidx;
    void load(const vector<string>& file_names);
    void write(const string& file_name);
    void make_compact(void);
    void make_dynamic(void);
    void add(const Alignment& aln);
    ostream& as_table(ostream& out);
private:
    bool is_compacted;
    // dynamic model
    gcsa::CounterArray coverage_dynamic;
    map<int64_t, map<string, int32_t> > edit_coverage_dynamic;
    // compact model
    vlc_vector<> coverage_civ; // graph coverage (compacted coverage_dynamic)
    wt_gmr<> edit_coverage_wt; // map from graph positions to edits
    vlc_vector<> edits_coverage_civ; // counts of the edits
    sd_vector<> edit_cbv; // mark beginning of edit record in edits_civ
    vlc_vector<> edits_civ; // serialized protobuf of edits
};

}

#endif
