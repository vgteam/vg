/**
 * \file index_manager.cpp: implementations of common indexing functionality
 */

#include <iostream>
#include <vector>
#include <string>

#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>

#include <gbwtgraph/minimizer.h>

#include "index_manager.hpp"

using namespace std;

namespace vg {

unique_ptr<gbwt::GBWT> IndexManager::make_gbwt(PathHandleGraph* xg_index, bool index_paths, const vector<string>& gam_file_names,  bool show_progress) {
    unique_ptr<gbwt::GBWT> to_return;
    return to_return;
}

}
