/** \file test_main.cpp
 *
 * Defines the "vg test" subcommand, which runs unit tests.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#define CATCH_CONFIG_RUNNER
#include "unittest/catch.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

// No help_test is necessary because the unit testing library takes care of
// complaining about missing options.

/**
 * Take the original argc and argv from a `vg unittest` command-line call and
 * run the unit tests. We keep this in its own CPP/HPP to keep our unit test
 * library from being a dependency of main.o and other real application code.
 *
 * Passes the args along to the unit test system.
 * 
 * Returns exit code 0 on success, other codes on failure.
 */
static int run_unit_tests(int argc, char** argv) {
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

int main_test(int argc, char** argv){
    // Forward arguments along to the main unit test driver
    return run_unit_tests(argc, argv);
}

// Register subcommand
static Subcommand vg_test("test", "run unit tests", DEVELOPMENT, main_test);

