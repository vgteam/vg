#ifndef VG_CONSTRUCTOR_HPP
#define VG_CONSTRUCTOR_HPP

/**
 * constructor.hpp: defines a tool class used for constructing VG graphs from
 * VCF files.
 */

#include <vector>
#include <set>

#include "types.hpp"

// We need vcflib
#include "Variant.h"
#include "Fasta.h"

namespace vg {

using namespace std;

/**
 * Represents a constructed region of the graph alogn a single linear sequence.
 * Contains the protobuf Graph holding all the created components (which may be
 * too large to serialize), a set of node IDs whose left sides need to be
 * connected to when you connect to the start of the chunk, and a set of node
 * IDs whose right sides need to be connected to when you connect to the end of
 * the chunk.
 */
struct ConstructedChunk {
    // What nodes, edges, and mappings exist?
    Graph graph;
    
    // What nodes have left sides that match up with the left edge of the chunk?
    set<id_t> left_ends;
    // And similarly for right sides on the right edge of the chunk?
    set<id_t> right_ends;
};

class Constructor {

public:

    // Should alts be interpreted as flat (false) or aligned back to the
    // reference by vcflib (true)?
    bool flat = false;
    
    // What's the maximum node size we should allow?
    size_t max_node_size = 1024;

    /**
     * Construct a ConstructedChunk of graph from the given piece of sequence,
     * with the given name, applying the given variants. The variants need to be
     * sorted by start position, have their start positions set relative to the
     * first base (0) of the given sequence, and not overlap with any variants
     * not in the vector we have (i.e. we need access to all overlapping
     * variants for this region). The variants must not extend beyond the given
     * sequence, though they can abut its edges.
     */
    ConstructedChunk construct_chunk(string reference_sequence, string reference_path_name,
        vector<vcflib::Variant> variants) const;

};

}

#endif
