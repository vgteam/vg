#include "gamsorter.hpp"
#include "stream.hpp"
#include <getopt.h>
#include "subcommand.hpp"
#include "index.hpp"
#include "stream.hpp"

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
         << "  -p / --paired           Index a paired-end GAM." << endl
         << "  -s / --sorted           Input GAM is already sorted." << endl
         << "  -i / --index <INDEX>    produce a node-to-alignment index" << endl
         << "  -d / --dumb-sort        use naive sorting algorithm (no tmp files, faster for small GAMs)" << endl
         << "  -r / --rocks            Just use the old RocksDB-style indexing scheme for sorting." << endl
         << "  -a / --aln-index        Create the old RocksDB-style node-to-alignment index." << endl
         << endl;
}

int main_gamsort(int argc, char **argv)
{
    string gamfile;
    bool do_index = false;
    bool dumb_sort = false;
    bool is_paired = false;
    bool is_sorted = false;
    bool just_use_rocks = false;
    bool do_aln_index = false;
    int c;
    optind = 2; // force optind past command positional argument
    while (true)
    {
        static struct option long_options[] =
            {
                {"index", required_argument, 0, 'i'},
                {"dumb-sort", no_argument, 0, 'd'},
                {"paired", no_argument, 0, 'p'},
                {"rocks", no_argument, 0, 'r'},
                {"aln-index", no_argument, 0, 'a'},
                {"is-sorted", no_argument, 0, 's'},
                {0, 0, 0, 0}};
        int option_index = 0;
        c = getopt_long(argc, argv, "idhraps",
                        long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            do_index = true;
            break;
        case 'd':
            dumb_sort = true;
            break;
        case 's':
            is_sorted = true;
            break;
        case 'r':
            just_use_rocks = true;
            break;
        case 'a':
            do_aln_index = true;
            break;
        case 'p':
            is_paired = true;
            break;
        case 'h':
        case '?':
        default:
            help_gamsort(argv);
            abort();
        }
    }


    if (argc <= 2){
        help_gamsort(argv);
        exit(11);
    }

    gamfile = argv[optind];

    GAMSorter gs;

    if (just_use_rocks && !do_index)
    {
        // Do the sort the old way - write a big ol'
        // RocksDB index of alignments, then dump them
        // from that DB. Loses unmapped reads.
        string dbname = gamfile + ".grai";
        Index index;

        index.open_for_bulk_load(dbname);
        int64_t aln_idx = 0;
        function<void(Alignment&)> lambda_reader = [&index](Alignment& aln) {
                index.put_alignment(aln);
        };
        get_input_file(gamfile, [&](istream& in) {
            stream::for_each_parallel(in, lambda_reader);
        });

        vector<Alignment> output_buf;
        auto lambda_writer = [&output_buf](const Alignment& aln) {
                output_buf.push_back(aln);
                stream::write_buffered(cout, output_buf, 100);
            };
        index.for_each_alignment(lambda_writer);
        ofstream outstream;
        outstream.open(dbname);
        //stream::write_buffered(outstream, output_buf, 0);
        index.flush();
        index.close();
    }

    if (dumb_sort)
    {
        gs.dumb_sort(gamfile);
    }
    else if (!dumb_sort)
    {
        cerr << "Smart sort not implemented :' ( " << endl;
        exit(11);
    }

    if (do_index && just_use_rocks)
    {
        string dbname = gamfile + ".grai";
        Index index;

        index.open_for_bulk_load(dbname);
        int64_t aln_idx = 0;
        function<void(Alignment&)> lambda_reader = [&index](Alignment& aln) {
                index.put_alignment(aln);
        };
        get_input_file(gamfile, [&](istream& in) {
            stream::for_each_parallel(in, lambda_reader);
        });

        vector<Alignment> output_buf;
        auto lambda_writer = [&output_buf](const Alignment& aln) {
                output_buf.push_back(aln);
                stream::write_buffered(cout, output_buf, 100);
            };
        index.for_each_alignment(lambda_writer);
        ofstream outstream;
        outstream.open(dbname);
        //stream::write_buffered(outstream, output_buf, 0);
        index.flush();
       
    }
    
    else if (do_index){
        // Write a super simple index file of format
        // node_d \tab firstAlnLineNumber \tab secondAlnLineNumber ...
        gs.write_index(gamfile, gamfile + ".gai", false);
    }

    return 1;
}

static Subcommand vg_gamsort("gamsort", "Perform naive sorts and indexing on a GAM file.", main_gamsort);
