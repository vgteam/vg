/**
 * \file protobuf_tags.cpp
 * Defines tag registrations for common Protobuf types.
 * Gets applied by being linked into the binary.
 */
 
#include "registry.hpp"

#include "../vg.pb.h"

namespace vg {

namespace stream {

using namespace std;

static ProtobufTag<Alignment> gam_tag("GAM");
static ProtobufTag<MultipathAlignment> mgam_tag("MGAM");
static ProtobufTag<Graph> vg_tag("VG");
static ProtobufTag<Snarl> snarl_tag("SNARL");


}

}
