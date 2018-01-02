#ifndef VG_POSITION_HPP_INCLUDED
#define VG_POSITION_HPP_INCLUDED

#include "vg.pb.h"
#include <iostream>
#include "json2pb.h"
#include "handle.hpp"
#include "dfs.hpp"

namespace vg {

using namespace std;

void for_each_kmer(const HandleGraph& graph, size_t k, const function<void(const string&, const vector<handle_t>&, const size_t)>& lambda);

}

#endif
