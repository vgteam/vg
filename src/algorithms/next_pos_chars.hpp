#pragma once

#include "../handle.hpp"
#include <vg/vg.pb.h>
#include <functional>
#include "../position.hpp"

namespace vg {
namespace algorithms {

using namespace std;

map<pos_t, char> next_pos_chars(const PathPositionHandleGraph& graph, pos_t pos);

}

}
