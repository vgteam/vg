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

#include "algorithms/topological_sort.hpp"
#include "algorithms/split_strands.hpp"

namespace vg {

using namespace std;

    class Surjector : BaseMapper {
    public:
        
        Surjector(xg::XG* xg_index);
        ~Surjector();
        
        /// Extract the portions of an alignment that are on a chosen set of paths and try to
        /// align realign the portions that are off of the chosen paths to the intervening
        /// path segments to obtain an alignment that is fully restricted to the paths.
        ///
        /// Also returns the path name, position, and strand of the new alignment.
        ///
        /// Optionally either allow softclips so that the alignment has a nonnegative score on
        /// the path or require the full-length alignment, possibly creating a negative score.
        Alignment surject(const Alignment& source,
                          const set<string>& path_names,
                          string& path_name_out,
                          int64_t& path_pos_out,
                          bool& path_rev_out,
                          bool allow_negative_scores = false);
        
        /// a local type that represents a read interval matched to a portion of the alignment path
        using path_chunk_t = pair<pair<string::const_iterator, string::const_iterator>, Path>;
        
    private:
        
        /// get the chunks of the alignment path that follow the given reference paths
        unordered_map<size_t, vector<path_chunk_t>>
        extract_overlapping_paths(const Alignment& source, const unordered_map<size_t, string>& path_rank_to_name,
                                  unordered_map<int64_t, vector<size_t>>* paths_of_node_memo = nullptr,
                                  unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo = nullptr);
        
        /// compute the widest interval of path positions that the realigned sequence could align to
        pair<size_t, size_t>
        compute_path_interval(const Alignment& source, size_t path_rank, const xg::XGPath& xpath, const vector<path_chunk_t>& path_chunks,
                              unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo = nullptr);
        
        /// make a linear graph that corresponds to a path interval, possibly duplicating nodes in case of cycles
        VG extract_linearized_path_graph(size_t first, size_t last, const xg::XGPath& xpath,
                                         unordered_map<id_t, pair<id_t, bool>>& node_trans);
        
        
        /// associate a path position and strand to a surjected alignment against this path
        void set_path_position(const Alignment& surjected, size_t best_path_rank, const xg::XGPath& xpath,
                               string& path_name_out, int64_t& path_pos_out, bool& path_rev_out,
                               unordered_map<pair<int64_t, size_t>, vector<pair<size_t, bool>>>* oriented_occurrences_memo = nullptr);
        
        // make a sentinel meant to indicate an unmapped read
        Alignment make_null_alignment(const Alignment& source);
    };
}

#endif
