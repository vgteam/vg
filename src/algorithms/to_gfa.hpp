#pragma once

#include "../handle.hpp"
#include <ostream>
#include <string>
#include <functional>

namespace vg {
namespace algorithms {

using namespace std;

void to_gfa(const PathHandleGraph& graph, ostream& out);

}
}
