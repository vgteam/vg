/**
 * \file md5_sum_path.cpp
 * Implementations of algorithms for MD5-summing paths
 */

#include "md5_sum_path.hpp"
#include <htslib/hts.h>

namespace vg {
namespace algorithms {

std::pair<std::string, size_t> md5_sum_path_with_length(const PathHandleGraph& graph, const path_handle_t& path) {
    size_t length = 0;
    hts_md5_context* context = hts_md5_init();
    std::string digest_string;
    try {
        for (handle_t handle : graph.scan_path(path)) {
            // Get the sequence of each visit along the path
            std::string visit_sequence = graph.get_sequence(handle);
            // Count its length
            length += visit_sequence.size();
            // And add it to the hash
            hts_md5_update(context, (const void*)visit_sequence.c_str(), visit_sequence.size() * sizeof(char));
        }
        // Get a binary MD5 digest
        unsigned char binary_digest[16];
        hts_md5_final(binary_digest, context);
        // Convert to a null-terminated hex string
        char hex_digest[33];
        hts_md5_hex(hex_digest, binary_digest);
        // Copy to a C++ string
        digest_string = std::string(hex_digest);
    } catch (...) {
        // Free alloated memory on exception
        hts_md5_destroy(context);
        throw;
    }
    // Free allocated memory when successful
    hts_md5_destroy(context);
    return std::make_pair(digest_string, length);
}

}
}
