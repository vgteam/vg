#pragma once

#include <iostream>
#include <string>
#include <list>
#include <handlegraph/util.hpp>
#include <handlegraph/handle_graph.hpp>
#include <gbwt/gbwt.h>
#include "position.hpp"
#include "gbwt_helper.hpp"

/** \file 
 * Functions for working with `kmers_t`'s in HandleGraphs.
 */

namespace vg {

namespace algorithms {

using namespace handlegraph;

/// Stores a walk in the context of a graph.
struct walk_t {
    walk_t(const std::string& s,
           const pos_t& b,
           const pos_t& e,
           const handle_t& c)
        : seq(s), begin(b), end(e), curr(c), path({c}) { };
    /// the walk
    std::vector<handle_t> path;
    /// the sequence
    std::string seq;
    /// our start position
    pos_t begin;
    /// Used in construction
    pos_t end; /// one past the (current) end of the walk
    handle_t curr; /// the next handle we extend into
    uint16_t forks; /// how many branching edge crossings we took to get here
};

/// Iterate over all the walks in the graph, running lambda on each
void for_each_walk(const HandleGraph& graph, size_t k, size_t edge_max,
                   const std::function<void(const walk_t&)>& lambda);

/// Print a walk_t to a stream.
std::ostream& operator<<(std::ostream& out, const walk_t& walk);

uint64_t walk_haplotype_frequency(const HandleGraph& graph,
                                  const gbwt::GBWT& haplotypes,
                                  const walk_t& walk);

std::vector<std::string> walk_haplotype_names(const HandleGraph& graph,
                                              const gbwt::GBWT& haplotypes,
                                              const walk_t& walk);

}

}
