#include "gamsorter.hpp"
#include "stream.hpp"
#include <getopt.h>
#include "subcommand.hpp"

/**
* GAM sort main
*/

using namespace std;
using namespace vg;
using namespace vg::subcommand;
void help_gamsort(char **argv)
{
    cerr << "gamsort: sort a GAM file (or index it) without Rocksdb" << endl
         << "Usage: " << argv[1] << " [Options] gamfile" << endl
         << "Options:" << endl
         << "  -i / --index <INDEX>    produce a node-to-alignment index"
         << "  -d / --dumb-sort      use naive sorting algorithm (no tmp files, faster for small GAMs)"
         << endl;
}

int main_gamsort(int argc, char **argv)
{
    string gamfile;
    string index_file;
    bool dumb_sort = false;
    bool is_paired = false;
    bool just_use_rocks = false;
    int c;
    optind = 2; // force optind past command positional argument
    while (true)
    {
        static struct option long_options[] =
            {
                {"index", required_argument, 0, 'i'},
                {"dumb-sort", no_argument, 0, 'd'},
                {"paired", no_argument, 0, 'p'},
                {"rocks", no_argument, 0, 'R'},
                {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "i:dhR",
                        long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            index_file = optarg;
            break;
        case 'd':
            dumb_sort = true;
            break;
        case 'R':
            just_use_rocks = true;
            break;
        case 'h':
        case '?':
        default:
            help_gamsort(argv);
            abort();
        }
    }

    gamfile = argv[optind];

    GAMSorter gs;

    if (just_use_rocks)
    {
        // Do the sort the old way - write a big ol'
        // RocksDB index of alignments, then dump them
        // from that DB. Loses unmapped reads.
    }

    if (dumb_sort)
    {
        gs.dumb_sort(gamfile);
    }
    else if (index_file.empty() && !dumb_sort)
    {
        cerr << "Smart sort not implemented :' ( " << endl;
        exit(11);
    }

    if (!index_file.empty())
    {
        // Write a super simple index file of format
        // node_d \tab firstAlnLineNumber \tab secondAlnLineNumber ...
        gs.write_index(gamfile, gamfile + ".gai", false);
    }

    return 1;
}

static Subcommand vg_gamsort("gamsort", "Perform naive sorts and indexing on a GAM file.", main_gamsort);
