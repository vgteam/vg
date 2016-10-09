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
    cerr << "Usage: " << argv[0] << " srpe [options] <data.gam> <data.gam.index> <graph.vg>" << endl
        << "Options: " << endl
        << "-m / --max-iter <MaximumIterations> max number of iterations for homogenization." << endl
        
        << endl;
        //<< "-S / --SV-TYPE comma separated list of SV types to detect (default: all)." << endl
        


}



int main_srpe(int argc, char** argv){
    string gam_name = "";
    string gam_index_name = "";
    string graph_name = "";
    int max_iter = 0;
    int max_frag_len = -1;
    int max_softclip = -1;

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
            {"max-iter", required_argument, 0, 'm'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "m:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'm':
                max_iter = atoi(optarg);
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
    Filter ff;
    vector<Alignment> alns;
    vector<string> sigtypes;
    vector<string> searches = {"DEL"};
    std::function<void(Alignment&, Alignment&)> flag_signature_read_support = [&ff, &alns, &sigtypes, &buffer](Alignment& aln_first, Alignment& aln_second){
        
        //srpe.apply(aln, alns);
        //srpe.apply(aln, positions);

        
    };


    gam_name = argv[optind];
    gam_index_name = argv[++optind];
    graph_name = argv[++optind];

    // Read in GAM and filter for Sigs
    
    // Convert sigs to edges and nodes

    // Find newly incorporated nodes / edges within X nodes or Y bp of each other
    // and remap reads to refine them.
    // Refinement: for N variant edges/nodes within X bp/nodes of one another, calculate the
    // sig quality ~ (reads mapped, mapping qual of reads, soft clipping, fragment length consistency,
    // )
    // Take the X highest-scoring edges / nodes and remap reads to just these.
    // Repeat this process until convergence, or after max-iter iterations

    // Can now emit refined variant calls.



}
