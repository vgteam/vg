#pragma once

#include "../handle.hpp"
#include "vg/io/edit.hpp"
#include <vg/vg.pb.h>
#include <string>
#include <functional>

namespace vg {
namespace algorithms {

using namespace std;

/// use the given oriented node sequence and the mapping to reconstruct the sequence represented by the mapping
void append_mapping_sequence(const Mapping& m, const string& node_seq, string& seq);
    
/// use the given graph and the path to determine our path string
std::string path_string(const HandleGraph& graph, const Path& path);

}
}
