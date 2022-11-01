/// \file protobuf.cpp
///  
/// unit tests for Protobuf linking and usage
///

#include <vg/vg.pb.h>
#include <google/protobuf/descriptor.h>

#include "../version.hpp"

#include "catch.hpp"

#include <iostream>

namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Protobuf reflection works", "[protobuf][io]") {

    std::cerr << "Testing libvgio..." << std::endl;
    
    // We need to actually instantiate a vg Protobuf (or otherwise depend on
    // its symbols at link time) or the vg Protobuf symbols might not come
    // along in the binary.
#ifdef debug
    std::cerr << "Creating Graph..." << std::endl;
#endif
    vg::Graph g;
#ifdef debug
    std::cerr << "Graph exists at " << &g << std::endl;
#endif

    // We also need to make sure to actually pull in something from libvg so it gets linked in if it is dynamic.
    auto version = vg::Version::get_version();
#ifdef debug
    std::cerr << "libvg version: " << version << std::endl;
#endif

#ifdef debug
    std::cerr << "Checking for Protobuf descriptor pool..." << std::endl;
#endif
    const google::protobuf::DescriptorPool* pool =  google::protobuf::DescriptorPool::generated_pool();
    if (pool == nullptr) {
        throw std::runtime_error("Cound not find Protobuf descriptor pool: is libvgio working?");
    }
#ifdef debug
    std::cerr << "Found descriptor pool at " << pool << std::endl;
#endif
    
    for (auto message_name : {"vg.Graph", "vg.Alignment", "vg.Position"}) {
#ifdef debug
        std::cerr << "Checking for Protobuf message type " << message_name << "..." << std::endl;
#endif
        const google::protobuf::Descriptor* descriptor = pool->FindMessageTypeByName(message_name);
        if (descriptor == nullptr) {
            throw std::runtime_error(std::string("Cound not find Protobuf descriptor for message type ") + message_name + ": is libvgio working?");
        }
#ifdef debug
        std::cerr << "Found " << message_name << " as " << descriptor->full_name() << " at " << descriptor << std::endl;
#endif
    }


}


}
}
