#include "driver.hpp"

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

namespace vg {
namespace unittest {

using namespace std;
using namespace vg;

int run_unit_tests(int argc, char** argv) {
    // argc and argv are going to have command and subcommand, which I don't
    // think argv speaks.
    
    assert(argc >= 2);
    
    // writing to session.configData() or session.Config() here 
    // overrides command line args
    // only do this if you know you need to

    // We're going to trick it by making a fake program name with a space in it
    auto new_program_name = string(argv[0]) + " " + string(argv[1]);
    
    // Delete an argument
    int fixed_argc = argc - 1;
    // Advance the pointer to the next char*
    char** fixed_argv = argv + 1;
    fixed_argv[0] = &new_program_name[0];
    
    // Make a Catch session
    Catch::Session session;

    int return_code = session.applyCommandLine(fixed_argc, fixed_argv);
    if(return_code != 0) {
        // Complain the user didn't specify good arguments
        return return_code;
    }
    
    // Actually run the tests
    return session.run();
    
}

}
}
