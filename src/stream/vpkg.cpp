/**
 * \file vpkg.cpp: Implementations for VPKG loader functions.
 */


#include "vpkg.hpp"
#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>

namespace vg {

namespace stream {

using namespace std;

void do_stuff() {
    VPKG::load_all<gcsa::GCSA, gcsa::LCPArray>(std::cin);
}

}

}
