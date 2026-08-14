/** \file help_main.cpp
 *
 * Defines the "vg help" subcommand, which describes subcommands.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include <iostream>

#include "subcommand.hpp"
#include "../utility.hpp"
#include "../version.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_help(char** argv) {
    cerr << "usage: " << argv[0] << " help <cmd>" << endl
         << "Print out information about all vg subcommands" << endl
         << "  -m, --man          output full Markdown-formatted manpage" << endl
         << "  -h, --help         print this help message to stderr and exit" << endl
         << endl;
}

int main_help(int argc, char** argv) {
    Logger logger("vg help");

    bool manpage = false;
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"man", no_argument, 0, 'm'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?m", long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'm':
                manpage = true;
                break;
            case 'h':
            case '?':
            default:
                help_help(argv);
                exit(1);
                break;
        }
    }

    if (optind < argc) {
        // We have a positional argument
        if (manpage) {
            logger.error() << "Cannot request both manpage and individual help at once" << endl;
        }
        string subcommand_name(argv[optind++]);
        if (optind < argc) {
            logger.error() << "Cannot request multiple helptexts at once" << endl;
        }

        // We have been asked for a specific subcommand's helptext
        argv[1] = &subcommand_name[0];
        auto* subcommand = vg::subcommand::Subcommand::get(argc, argv);
        if (subcommand != nullptr) {
            subcommand->run_help(argv);
            return 0;
        } else if (vg::subcommand::REMOVED_CMD_MESSAGES.count(subcommand_name)) {
            // We have a special message prepared!
            logger.error() << subcommand_name << " has been removed:\n"
                           << vg::subcommand::REMOVED_CMD_MESSAGES.at(subcommand_name) << endl;
            return 1;
        } else {
            // No subcommand found
            logger.error() << "command " << argv[1] << " not found" << endl;
            return 1;
        }
    }

    if (manpage) {
        cerr << "## NAME" << endl
             << endl;
    }

    cerr << "vg: variation graph tool, version " << Version::get_short() << endl
         << endl;

    if (manpage) {
        cerr << "## DESCRIPTION" << endl
             << endl
             << "[vg](https://github.com/vgteam/vg) is a toolkit for variation graph data structures, "
             << "interchange formats, alignment, genotyping, and variant calling methods." << endl
             << endl
             << "For more in-depth explanations of tools and workflows, "
             << "see the [general wiki page](https://github.com/vgteam/vg/wiki)." << endl
             << endl
             << "## SYNOPSIS" << endl
             << endl;

        cerr << "## COMMANDS" << endl
             << endl;
        vg::subcommand::Subcommand::for_each([&argv](const vg::subcommand::Subcommand& command) {
                // Print out all helptext
                cerr << "### " << command.get_name() << ": " << command.get_description() << endl
                     << endl
                     << "```" << endl;
                command.run_help(argv);
                cerr << "```" << endl
                     << endl;
            });
    } else {
        cerr << "usage: " << argv[0] << " <command> [options]" << endl
             << endl;
        // We're just going through a list of names
        for (auto category : {PIPELINE, TOOLKIT, WIDGET, DEVELOPMENT}) {

            cerr << category << ":" << endl;
            
            vg::subcommand::Subcommand::for_each(category, [](const vg::subcommand::Subcommand& command) {
                
                // Announce every subcommand we have
                
                // Pad all the names so the descriptions line up
                string name = command.get_name();
                name.resize(18, ' ');
                cerr << "  -- " << name << command.get_description() << endl;
            });
            
            cerr << endl;
        }
    }

    if (manpage) {
        cerr << "## BUGS" << endl
             << endl
             << "Bugs can be reported at: https://github.com/vgteam/vg/issues" << endl
             << endl;
    }
     
    cerr << "For technical support, please visit: https://www.biostars.org/tag/vg/" << endl << endl;
     
    
    return 0;
}

// Register subcommand
static Subcommand vg_help("help", "show all subcommands", PIPELINE, help_help, main_help);

