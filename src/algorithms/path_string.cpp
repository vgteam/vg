#include "path_string.hpp"

namespace vg {
namespace algorithms {

using namespace std;
using namespace vg::io;

void append_mapping_sequence(const Mapping& m, const string& node_seq, string& seq) {
    size_t t = 0;
    size_t f = m.position().offset();
    for (size_t i = 0; i < m.edit_size(); ++i) {
        auto& e = m.edit(i);
        if (edit_is_match(e)) {
            seq.append(node_seq.substr(f, e.from_length()));
        } else if (edit_is_sub(e)) {
            seq.append(e.sequence());
        } else if (edit_is_insertion(e)) {
            seq.append(e.sequence());
        } else if (edit_is_deletion(e)) {
            // no-op
        }
        t += e.to_length();
        f += e.from_length();
    }
}

string path_string(const HandleGraph& graph, const Path& path) {
    string seq;
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        append_mapping_sequence(m,
                                graph.get_sequence(
                                    graph.get_handle(m.position().node_id(),
                                                     m.position().is_reverse())),
                                seq);
    }
    return seq;
}
    
}
}
