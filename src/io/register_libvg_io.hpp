#ifndef VG_IO_REGISTER_LIBVG_IO_HPP_INCLUDED
#define VG_IO_REGISTER_LIBVG_IO_HPP_INCLUDED

/**
 * \file register_libvg_io.hpp
 * Includes a function to call to register IO handlers for libvg types.
 */

namespace vg {

namespace io {

using namespace std;

/**
 * Register libvg types with libvgio.
 * Must be called by library users before doing IO.
 * Does not magically statically call itself.
 * TODO: work out a way it can.
 * Returns true on success.
 */
bool register_libvg_io();

}

}

#endif
