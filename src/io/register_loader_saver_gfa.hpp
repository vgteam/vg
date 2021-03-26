#ifndef VG_IO_REGISTER_LOADER_SAVER_GFA_HPP_INCLUDED
#define VG_IO_REGISTER_LOADER_SAVER_GFA_HPP_INCLUDED

/**
 * \file register_loader_saver_gfa.hpp
 * Defines IO for a graph in GFA format.  It's best if the GFA is "canonical"
 * with S lines before L lines before P lines.
 */

namespace vg {

namespace io {

using namespace std;

void register_loader_saver_gfa();

}

}

#endif
