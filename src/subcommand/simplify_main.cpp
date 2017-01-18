// simplify_main.cpp: define the "vg simplify" subcommand, which removes small variation

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../simplifier.hpp"




using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_simplify(char** argv) {
    cerr << "usage: " << argv[0] << " simplify [options] old.vg >new.vg" << endl
         << "options:" << endl
         << "    -m, --min-size N       remove leaf sites with fewer than N bases involved (default: 10)" << endl
         << "    -i, --max-iterations N perform up to N iterations of simplification (default: 10)" << endl
         << "    -p, --progress         show progress" << endl
         << "    -b, --bed-in           read in the given BED file in the cordinates of the original paths" << endl
         << "    -B, --bed-out          output transformed features in the coordinates of the new paths" << endl
         << "    -t, --threads N        use N threads to construct graph (defaults to numCPUs)" << endl;
}

int main_simplify(int argc, char** argv) {

    if (argc == 2) {
        help_simplify(argv);
        return 1;
    }


    // TODO: The simplifier needs the graph when we make it, so we can't store
    // our settings in it directly.
    size_t min_size = 10;
    size_t max_iterations = 10;
    string bed_in_filename;
    string bed_out_filename;
    bool show_progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"min-size", required_argument, 0, 'm'},
                {"max-iterations", required_argument, 0, 'i'},
                {"progress",  no_argument, 0, 'p'},
                {"bed-in", required_argument, 0, 'b'},
                {"bed-out", required_argument, 0, 'B'},
                {"threads", required_argument, 0, 't'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "m:i:pb:B:t:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'm':
            min_size = atoi(optarg);
            break;
            
        case 'i':
            max_iterations = atoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;
            
        case 'b':
            bed_in_filename = optarg;
            break;
        
        case 'B':
            bed_out_filename = optarg;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_simplify(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }
    
    // TODO: move all this to a simplifier object
    
    // Load the graph
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in, show_progress);
    });
    
    if (graph == nullptr) {
        cerr << "error:[vg simplify]: Could not load graph" << endl;
        exit(1);
    }
    
    {
        // Make a Simplifier for the graph and copy over settings. Need sto be
        // in a block so that the graph doesn't get deleted before the
        // simplifier goes out of scope.
        Simplifier simplifier(*graph);
        simplifier.show_progress = show_progress;
        simplifier.max_iterations = max_iterations;
        simplifier.min_size = min_size;
        
        if (!bed_in_filename.empty()) {
            // Load BED features
            ifstream bed_stream(bed_in_filename.c_str());
            simplifier.features.load_bed(bed_stream);
        }
        
        // Do the simplification
        simplifier.simplify();
        
        // Serialize the graph
        graph->serialize_to_ostream(std::cout);
        
        if (!bed_out_filename.empty()) {
            // Save BED features
            ofstream bed_stream(bed_out_filename.c_str());
            simplifier.features.save_bed(bed_stream);
        }
        
    }
    
    delete graph;

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_construct("simplify", "graph simplification", main_simplify);

