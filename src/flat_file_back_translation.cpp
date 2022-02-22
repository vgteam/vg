/**
 * \file flat_file_back_translation.cpp
 * Implementation for flat-file-backed named-node back-translation.
 */

#include "flat_file_back_translation.hpp"
#include "utility.hpp"

namespace vg {

FlatFileBackTranslation::FlatFileBackTranslation(std::istream& stream) {
    if (!stream) {
        // File didn't open properly or something.
        throw std::runtime_error("Could not read translation from stream");
    }
    while (stream) {
        // Get each line.
        std::string line;
        std::getline(stream, line);
        if (line.empty()) {
            // Skip blank lines.
            continue;
        }
        // Split on tabs.
        auto parts = split_delims(line, "\t");
        if (parts[0] == "T") {
            // This is T <tab> segment name <tab> segment number
            if (parts.size() != 3) {
                throw std::runtime_error("Encountered unparseable T line: " + line);
            }
            segment_to_name[parse<nid_t>(parts[2])] = parts[1];
        } else if (parts[0] == "K") {
            // This is K <tab> old ID <tab> forward offset <tab> reverse offset <tab> new ID
            if (parts.size() != 5) {
                throw std::runtime_error("Encountered unparseable K line: " + line);
            }
            // Save, under the new node ID, the old node ID and the offsets.
            node_to_segment_and_offsets[parse<nid_t>(parts[4])] =  {parse<nid_t>(parts[1]), parse<size_t>(parts[2]), parse<size_t>(parts[3])};
        } else {
            // This shouldn't be here.
            throw std::runtime_error("Encountered unrecognized line: " + line);
        }
    }
}

std::vector<oriented_node_range_t> FlatFileBackTranslation::translate_back(const oriented_node_range_t& range) const {
    // Look up the node ID
    auto it = node_to_segment_and_offsets.find(std::get<0>(range));
    
    if (it == node_to_segment_and_offsets.end()) {
        // This doesn't have to go anywhere else.
        return {range};
    }
    
    // Otherwise this goes somewhere else.
    return {{
             // The destination segment
             std::get<0>(it->second),
             // In the requested orientation
             std::get<1>(range),
             // Starting at the correct offset along that orientation of the segment
             std::get<2>(range) + (std::get<1>(range) ? std::get<2>(it->second) : std::get<1>(it->second)),
             // And running for the specified length
             std::get<3>(range)
           }};
}

std::string FlatFileBackTranslation::get_back_graph_node_name(const nid_t& back_node_id) const {
    auto it = segment_to_name.find(back_node_id);
    if (it == segment_to_name.end()) {
        // There's no non-default name for this segment.
        return std::to_string(back_node_id);
    }
    
    // Otherwise, we found a name.
    return it->second;
}

}
