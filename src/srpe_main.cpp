#include <iostream>
#include <vector>
#include <getopt.h>
#include "stream.hpp"
#include "index.hpp"
#include "vg.pb.h"
#include "srpe.hpp"

using namespace std;
using namespace vg;

void help_srpe(char** argv){
    cerr << "Usage: " << argv[0] << " srpe [options] <data.gam>" << endl
        << "Options: " << endl
        << ""
        << endl;

}

int main_srpe(int argc, char** argv){
    string gam_name = "";
    string index_name = "";

    /*
     * int sc_cutoff
     * int max_dist between reads
     *
     */

    if (argc <= 2) {
        help_srpe(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"index", required_argument, 0, 'i'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "i:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'i':
                index_name = optarg;
                break;
            case 'h':
            case '?':
            default:
                help_srpe(argv);
                abort();
        }

    }

    vector<Alignment> buffer;
    static const int buffer_size = 1000; // we let this be off by 1
    function<Alignment&(uint64_t)> write_buffer = [&buffer](uint64_t i) -> Alignment& {
            return buffer[i];
    };


    SRPE srpe;
    vector<Alignment> alns;
    vector<string> positions;
    std::function<void(Alignment&)> lambda = [&srpe, &alns, &positions, &buffer](Alignment& aln){
        
        //srpe.apply(aln, alns);
        srpe.apply(aln, positions);
        /* For each alignment in the gam:
         * Find the node the alignment maps to
         * find its mate read, if it has one
         * Check if alignment/mate fail any filters
         * if they do:
         *  find all the alignments that map to that node
         *  Look for matching signatures on other reads at the node
         *  Check the depth at the Locus, and update it as well
         *
         */

    };

    
    /*
     * now that we've been through every read, check for erroneous depth signals
     */


    gam_name = argv[optind];

    if (gam_name == "-"){
        stream::for_each_parallel(cin, lambda);
    }
    else{
        ifstream in;
        in.open(gam_name);
        if (in.good()){
            stream::for_each_parallel(in, lambda);
        }
        else{
            cerr << "Could not open " << gam_name << endl;
            help_srpe(argv);
        }
    }

    if (buffer.size() > 0) {
        stream::write(cout, buffer.size(), write_buffer);
        buffer.clear();
    }


}
