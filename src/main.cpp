// Needed for crash.hpp to work because it uses newer types
#define _POSIX_C_SOURCE 200809L

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
     cerr << "For technical support, please visit: https://www.biostars.org/t/vg/" << endl << endl;
 }

// We make sure to compile main for the lowest common denominator architecture.
// This macro is defined in the preflight header on supported compiler setups.
// But to use it we have to declare and then define main.
int main(int argc, char *argv[]) VG_PREFLIGHT_EVERYWHERE;

int main(int argc, char *argv[]) {

    // Make sure the system meets system requirements (i.e. has all the instructions we need)
    preflight_check();
    
    // Make sure we configure the memory allocator appropriately for our environment
    //configure_memory_allocator();
    
    // Set up stack trace support from crash.hpp
    enable_crash_handling();

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
        return (*subcommand)(argc, argv);
    } else {
        // No subcommand found
        string command = argv[1];
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

}
