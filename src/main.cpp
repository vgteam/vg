#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <csignal>
#include <getopt.h>
#include <sys/stat.h>

#include "google/protobuf/stubs/common.h"
#include "version.hpp"
#include "utility.hpp"
#include "crash.hpp"
#include "preflight.hpp"
#include "config/allocator_config.hpp"
#include "io/register_libvg_io.hpp"

// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"

using namespace std;
using namespace google::protobuf;
using namespace vg;


void vg_help(char** argv) {
    cerr << "vg: variation graph tool, version " << Version::get_short() << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << vg::subcommand::PIPELINE << ":" << endl;
         
     vg::subcommand::Subcommand::for_each(vg::subcommand::PIPELINE, [](const vg::subcommand::Subcommand& command) {
        // Announce every subcommand we have
        
        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });
     
     cerr << endl << "For more commands, type `vg help`." << endl;
     cerr << "For technical support, please visit: https://www.biostars.org/tag/vg/" << endl << endl;
 }

/// Main entry point once we know we're on a supported CPU.
int vg_main(int argc, char *argv[]) {

    // Make sure we configure the memory allocator appropriately for our environment
    AllocatorConfig::configure();
    
    // Set up stack trace support from crash.hpp
    enable_crash_handling();
    set_crash_context("Starting up");

    // Determine a sensible default number of threads and apply it.
    choose_good_thread_count();

    // Determine temp directory from environment variables.
    temp_file::set_system_dir();

    // Tell the IO library about libvg types.
    // TODO: Make a more generic libvg startup function?
    if (!vg::io::register_libvg_io()) {
        cerr << "error[vg]: Could not register libvg types with libvgio" << endl;
        return 1;
    }
    
    if (argc == 1) {
        vg_help(argv);
        return 1;
    }
    
    auto* subcommand = vg::subcommand::Subcommand::get(argc, argv);
    if (subcommand != nullptr) {
        // We found a matching subcommand, so run it
        if (subcommand->get_category() == vg::subcommand::CommandCategory::DEPRECATED) {
            cerr << endl << "WARNING:[vg] Subcommand '" << argv[1] << "' is deprecated and is no longer being actively maintained. Future releases may eliminate it entirely." << endl << endl;
        }
        set_crash_context("Starting '" +  std::string(argv[1]) + "' subcommand");
        return (*subcommand)(argc, argv);
    } else {
        // No subcommand found
        string command = argv[1];
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

}

// We make sure to compile main for the lowest common denominator architecture.
// This macro is defined in the preflight header on supported compiler setups.
// But to use it we have to declare and then define main.
// Note that on GCC 13.1 the always-inline allocator functions can't be inlined
// into code for architectures this old, causing an error if we try and
// allocate or use std::string. So the real main() function can't use C++
// allocators.
int main(int argc, char *argv[]) VG_PREFLIGHT_EVERYWHERE;

// TODO: What about static initialization code? It might use instructions not
// supported on the current CPU!

/// Make sure the system meets system requirements (i.e. has all the
/// instructions we need), then call vg_main
int main(int argc, char** argv) {
    preflight_check();
    return vg_main(argc, argv);
}
