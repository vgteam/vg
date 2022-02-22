/**
 * \file flat_file_back_translation.hpp
 * Defines a back-translation from graph node ID space to named node space,
 * backed by a flat text file.
 */

#ifndef VG_FLAT_FILE_BACK_TRANSLATION_HPP_INCLUDED
#define VG_FLAT_FILE_BACK_TRANSLATION_HPP_INCLUDED

#include "handle.hpp"

#include <unordered_map>
#include <istream>
#include <vector>

namespace vg {

/**
 * A NamedNodeBackTranslation loadable from a text file.
 *
 * The file is a GFA-like tab-separated file with types of lines identified by
 * a letter in the first field. It consists of 0 or more T lines, each giving a
 * segment name and an assigned number for it. This is followed by 0 or more K
 * lines, each giving a segment number, offsets along that segment in forward
 * and then reverse orientations, and then the graph node ID that begins at
 * that offset.
 *
 * Note that an empty file is allowed, and that this format can only represent
 * translations where nodes are broken up (and not merged) and where
 * orientation does not change.
 *
 * Many applications (such as loading the translation into a GBWTGraph) will
 * expect the graph node IDs alogn a segment to be dense, contiguous, and
 * increasing.
 */
class FlatFileBackTranslation : public NamedNodeBackTranslation {

public:
    /**
     * Create a FlatFileBackTranslation by reading it from an open file.
     */
    FlatFileBackTranslation(std::istream& stream);
    
    virtual ~FlatFileBackTranslation() = default;
    
    /**
     * Translate the given range of bases on the given orientation of the given
     * node in the current graph, to zero or more ranges on orientations of
     * nodes in some prior graph.
     */
    virtual std::vector<oriented_node_range_t> translate_back(const oriented_node_range_t& range) const;
    
    /**
     * Get the name of a node in the graph that translate_back() translates
     * into, given its number.
     */
    virtual std::string get_back_graph_node_name(const nid_t& back_node_id) const;
    
protected:
    /**
     * This holds, for each node ID, the segment number and starting offset, if
     * it is not offset 0, on each orientation of the segment with the same
     * number as the node ID.
     */
    std::unordered_map<nid_t, std::tuple<nid_t, size_t, size_t>> node_to_segment_and_offsets;
    
    /**
     * This holds, for each segment ID, the segment name, if it is not the
     * string version of the segment number.
     */
    std::unordered_map<nid_t, std::string> segment_to_name;

};

}

#endif
