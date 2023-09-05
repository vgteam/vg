/**
 * \file register_libvg_io.hpp
 * Includes calls to register all libvg types with libvgio.
 */

// Keep these includes in alphabetical order.

#include "register_loader_saver_distance_index.hpp"
#include "register_loader_saver_gbwt.hpp"
#include "register_loader_saver_r_index.hpp"
#include "register_loader_saver_gbwtgraph.hpp"
#include "register_loader_saver_gbz.hpp"
#include "register_loader_saver_gbzgraph.hpp"
#include "register_loader_saver_gcsa.hpp"
#include "register_loader_saver_lcp.hpp"
#include "register_loader_saver_minimizer.hpp"
#include "register_loader_saver_snarl_manager.hpp"
#include "register_loader_saver_vg.hpp"
#include "register_loader_saver_xg.hpp"
#include "register_loader_saver_packed_graph.hpp"
#include "register_loader_saver_hash_graph.hpp"
#include "register_loader_saver_gfa.hpp"

#include "register_libvg_io.hpp"


namespace vg {

namespace io {

using namespace std;

bool register_libvg_io() {
    register_loader_saver_distance_index();
    register_loader_saver_gbwt();
    register_loader_saver_r_index();
    register_loader_saver_gbwtgraph();
    register_loader_saver_gbz();
    register_loader_saver_gbzgraph();
    register_loader_saver_gcsa();
    register_loader_saver_lcp();
    register_loader_saver_minimizer();
    register_loader_saver_snarl_manager();
    register_loader_saver_vg();
    register_loader_saver_gfa();
    register_loader_saver_xg();
    register_loader_saver_packed_graph();
    register_loader_saver_hash_graph();
    return true;
}
    
}

}
