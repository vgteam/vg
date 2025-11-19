#include <getopt.h>

#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include <vg/io/vpkg.hpp>
#include "../haplotype_extracter.hpp"
#include "../algorithms/find_gbwt.hpp"
#include <bdsg/overlays/overlay_helper.hpp>

using namespace vg;
using namespace std;
using namespace vg::subcommand;

void help_trace(char** argv) {
    cerr << "usage: " << argv[0] << " trace [options]" << endl
         << "Trace and extract haplotypes from an index" << endl
         << endl
         << "options:" << endl
         << "  -x, --index FILE            use this XG index or graph" << endl
         << "  -G, --gbwt-name FILE        use GBWT haplotype index instead of any in graph" << endl
         << "  -n, --start-node INT        start at this node ID" << endl
        //TODO: implement backwards iteration over graph
        // << "  -b, --backwards             iterate backwards over graph" << endl
         << "  -d, --extend-distance INT   extend search this many nodes [50]" << endl
         << "  -a, --annotation-path FILE  output file for haplotype frequency annotations" << endl
         << "  -j, --json                  output subgraph in json instead of protobuf" << endl
         << "  -h, --help                  print this help message to stderr and exit" << endl;
}

int main_trace(int argc, char** argv) {
    Logger logger("vg trace");
    if (argc == 2) {
        help_trace(argv);
        return 1;
    }

    string xg_name;
    string gbwt_name;
    string annotation_path;
    int64_t start_node = 0;
    int extend_distance = 50;
    bool backwards = false;
    bool json = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"index", required_argument, 0, 'x'},
                {"gbwt-name", required_argument, 0, 'G'},
                {"annotation-path", required_argument, 0, 'a'},
                {"start-node", required_argument, 0, 'n'},
                {"extend-distance", required_argument, 0, 'd'},
                {"json", no_argument, 0, 'j'},
                {"help", no_argument, 0, 'h'},
                //{"backwards", no_argument, 0, 'b'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "x:G:a:n:d:jh?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'x':
            xg_name = require_exists(logger, optarg);
            break;

        case 'G':
            gbwt_name = require_exists(logger, optarg);
            break;

        case 'a':
            annotation_path = ensure_writable(logger, optarg);
            break;

        case 'n':
            start_node = parse<int>(optarg);
            break;

        case 'd':
            extend_distance = parse<int>(optarg);
            break;

        case 'j':
            json = true;
            break;

        //case 'b':
            //backwards = true;
            //break;

        case '?':
        case 'h':
            help_trace(argv);
            return 1;

        default:
            help_trace(argv);
            abort();
        }
    }

    if (xg_name.empty()) {
        logger.error() << "XG index must be specified with -x" << endl;
    }
    if (start_node < 1) {
        logger.error() << "start node must be specified with -n" << endl;
    }
    unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* xindex = overlay_helper.apply(path_handle_graph.get());    

    // Now load the haplotype data
    unique_ptr<gbwt::GBWT> gbwt_index_holder;
    const gbwt::GBWT* gbwt_index = vg::algorithms::find_gbwt(path_handle_graph.get(), gbwt_index_holder, gbwt_name);

    if (gbwt_index == nullptr) {
        // Complain if we couldn't.
        logger.error() << "unable to find GBWT index in graph or separate file" << endl;
    }
    
    // trace out our graph and paths from the start node
    Graph trace_graph;
    map<string, int> haplotype_frequences;
    trace_haplotypes_and_paths(*xindex, *gbwt_index, start_node, extend_distance,
                               trace_graph, haplotype_frequences);

    // dump our graph to stdout
    if (json) {
        cout << pb2json(trace_graph);
    } else {
        VG vg_graph;
        vg_graph.extend(trace_graph);
        vg_graph.serialize_to_ostream(cout);
    }

    // if requested, write thread frequencies to a file
    if (!annotation_path.empty()) {
        ofstream annotation_file(annotation_path);
        for (auto tf : haplotype_frequences) {
            annotation_file << tf.first << "\t" << tf.second << endl;
        }
    }

    return 0;
}

static Subcommand vg_trace("trace", "trace haplotypes", main_trace);
