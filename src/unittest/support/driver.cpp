#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

/**
 * Main function for unit test suite binaries.
 */
int main(int argc, char** argv) {
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
