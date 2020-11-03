#ifndef VG_PRUNE_HPP_INCLUDED
#define VG_PRUNE_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "handle.hpp"
#include "hash_map.hpp"
#include "position.hpp"

/** \file 
 * Functions for working with `kmers_t`'s in HandleGraphs.
 */

namespace vg {

using namespace std;

/// Record a <=k-length walk in the context of a graph.
struct walk_t {
    walk_t(uint16_t l,
           const pos_t& b,
           const pos_t& e,
           const handle_t& c,
           uint16_t f)
        : length(l), begin(b), end(e), curr(c), forks(f) { };
    /// our start position
    pos_t begin;
    pos_t end; /// one past the (current) end of the kmer
    handle_t curr; /// the next handle we extend into
    uint16_t forks; /// how many branching edge crossings we took to get here
    uint16_t length; /// how far we've been
};

/// Iterate over all the walks up to length k, adding edges which 
pair_hash_set<edge_t> find_edges_to_prune(const HandleGraph& graph, size_t k, size_t edge_max);

}

#endif
