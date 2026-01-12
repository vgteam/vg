#define CATCH_CONFIG_RUNNER
// Make sure to grab catch relative to us
#include "../catch.hpp"
#include "../../io/register_libvg_io.hpp"
#include "../../utility.hpp"
#include "../../log.hpp"

using namespace vg;

/**
 * Main function for unit test suite binaries.
 */
int main(int argc, char** argv) {
    Logger logger("driver.cpp");

    // Determine temp directory from environment variables.
    temp_file::set_system_dir();

    // Tell the IO library about libvg types.
    // TODO: Make a more generic libvg startup function?
    if (!vg::io::register_libvg_io()) {
        logger.error() << "Could not register libvg types with libvgio" << endl;
    }

    // Make a Catch session
    Catch::Session session;

    int return_code = session.applyCommandLine(argc, argv);
    if(return_code != 0) {
        // Complain the user didn't specify good arguments
        return return_code;
    }
    
    // Actually run the tests
    return session.run();
}
