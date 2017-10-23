/** \file bugs_main.cpp
 *
 * Defines the "vg bugs" subcommand, which shows and opens Github issues.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>

#include <iostream>

#include "subcommand.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_bugs(char** argv){
    cerr << "usage: " << argv[0] << " bugs" << endl
         << "options: " << endl
         << "--new, -n      make a new bug report" << endl
         << endl;
}

int main_bugs(int argc, char** argv){

    bool new_bug = false;
    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"new", no_argument, 0, 'n'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "n",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'n':
            new_bug = true;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_bugs(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    // What should we open?
    string url = new_bug ? string("https://github.com/vgteam/vg/issues/new") : string("https://github.com/vgteam/vg/issues");
    
    // What should we run?
#ifdef __APPLE__
    string open_command = "open";
#else
    string open_command = "xdg-open";
#endif

    // Open the URL in the appropriate browser (which may be lynx or similar)
    return system((open_command + " " + url).c_str());
    
}

// Register subcommand
static Subcommand vg_bugs("bugs", "show or create bugs", main_bugs);

