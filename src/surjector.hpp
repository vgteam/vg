#ifndef VG_SURJECTOR_HPP_INCLUDED
#define VG_SURJECTOR_HPP_INCLUDED

/** \file
 *
 *  A class to hold surjection algorithms that do lossy realignment restricted to paths in the graph
 */

#include <set>

#include "alignment.hpp"
#include "mapper.hpp"
#include "xg.hpp"
#include "vg.hpp"
#include "translator.hpp"
#include "vg.pb.h"
#include "multipath_alignment_graph.hpp"

namespace vg {

using namespace std;

    class Surjector : Mapper {
    public:
        
        Surjector(xg::XG* xg_index);
        ~Surjector();
        
        // lossily project an alignment into a particular path space of a graph
        // the resulting alignment is equivalent to a SAM record against the chosen path
        Alignment surject_classic(const Alignment& source,
                                  const set<string>& path_names,
                                  string& path_name,
                                  int64_t& path_pos,
                                  bool& path_reverse);
        
        Alignment path_anchored_surject(const Alignment& source, const string& path_name);
        
        
        
    private:
        
        vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>
        extract_overlapping_paths(const Alignment& source, size_t path_rank,
                                  unordered_map<int64_t, vector<size_t>>* paths_of_node_memo = nullptr);
        
        pair<size_t, size_t>
        compute_path_interval(const Alignment& source, size_t path_rank, const xg::XGPath& xpath,
                              const vector<pair<pair<string::const_iterator, string::const_iterator>, Path>>& path_chunks,
                              unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo = nullptr);
        
        VG extract_linearized_path_graph(size_t begin, size_t end, const xg::XGPath& xpath,
                                         unordered_map<id_t, pair<id_t, bool>>& node_trans);
        
    };
}
 
#endif
