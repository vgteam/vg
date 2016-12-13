#ifndef VG_UNITTEST_DRIVER_H
#define VG_UNITTEST_DRIVER_H

namespace vg {
namespace unittest {

using namespace std;
using namespace vg;
    
/**
 * Take the original argc and argv from a `vg unittest` command-line call and
 * run the unit tests. We keep this in its own CPP/HPP to keep our unit test
 * library from being a dependency of main.o and other real application code.
 *
 * Passes the args along to the unit test system.
 * 
 * Returns exit code 0 on success, other codes on failure.
 */
int run_unit_tests(int argc, char** argv);

}
}

#endif
