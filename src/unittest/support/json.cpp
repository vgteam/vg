#include "json.hpp"

#include "vg/io/json2pb.h"
#include "vg.hpp"

namespace vg {
namespace unittest {

std::unique_ptr<MutablePathMutableHandleGraph> json_to_graph(const std::string& json) {
    // Load into a Protobuf object
    Graph source;
    json2pb(source, json.c_str(), json.size());
   
    // Make a HandleGraph that knows how to load from Protobuf
    std::unique_ptr<MutablePathMutableHandleGraph> to_return = new vg::VG();

    // Load it from Protobuf
    to_return->extend(source);
    return to_return;
}
    

}
}
