#include <getopt.h>

#include <string>
#include <vector>

#include <subcommand.hpp>

#include "../primer_filter.hpp"
#include "../snarl_distance_index.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_primers(char** argv) {
    cerr << "usage: " << argv[0] << " primers [options] input.primer3 > filtered_primers.out" << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-path FILE         use this xg graph" << endl
         << "    -s, --snarl-index FILE     use this snarl index" << endl
         << "    -z, --zero-variance        allow no variance in the product" << endl
         << "    -l, --tolerance INT        allow this much difference between minimum and maximum sizes compared to the linear product size (default: 10)" << endl
         << "    -n, --minimum-size INT     minimum product size allowed (has precedence over --tolerance)" << endl
         << "    -m, --maximum-size INT     maximum product size allowed (has precedence over --tolerance)" << endl
         << "    -a, --all-primers          output all primers" << endl;
}

int main_primers(int argc, char** argv) {
    
    if (argc == 2) {
        help_primers(argv);
        return 1;
    }

    string xg_path;
    string snarl_index_path;
    bool zero_variance = false;
    int tolerance = 20;
    int minimum_product_size = numeric_limits<int>::max();
    int maximum_product_size = numeric_limits<int>::max();

    int c;
    optind = 2;

    while (true) {
        static struct option long_options[] =
        {
          {"help",          no_argument,       0, 'h'},
          {"xg-path",       required_argument, 0, 'x'},
          {"snarl-index",   required_argument, 0, 's'},
          {"zero-variance", required_argument, 0, 'z'},
          {"tolerance",     required_argument, 0, 'l'},
          {"minimum-size",  required_argument, 0, 'n'},
          {"maximum-size",  required_argument, 0, 'm'},
          {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:s:zl:n:m:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {
        case 'x':
            xg_path = optarg;
            break;
        
        case 's':
            snarl_index_path = optarg;
            break;

        case 'z':
            zero_variance = true;
            break;
        
        case 'l':
            tolerance = parse<int>(optarg);
            break;

        case 'n':
            minimum_product_size = parse<int>(optarg);
            break;

        case 'm':
            maximum_product_size = parse<int>(optarg);
            break;

        case 'h':
        case '?':
            help_primers(argv);
            exit(1);
            break;

        default:
          abort ();
        }
    }

    if (xg_path.empty()) {
        cerr << "error:[vg primers] xg file (-x) is required" << endl;
        exit(1);
    }

    if (snarl_index_path.empty()) {
        cerr << "error:[vg priemrs] snarl index file (-s) is required" << endl;
        exit(1);
    }

    string primers_path = get_input_file_name(optind, argc, argv);

    cout << "primer file name: "      << primers_path     << endl
         << "xg file name: "          << xg_path          << endl
         << "snarl index file name: " << snarl_index_path << endl;
    
    SnarlDistanceIndex distance_index;
    unique_ptr<handlegraph::PathPositionHandleGraph> graph;
    distance_index.deserialize(snarl_index_path);
    graph = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_path);
    PrimerFinder primer_finder(graph, &distance_index);
    primer_finder.load_primers(primers_path);

    

    return 0;
}

static Subcommand vg_primers("primers", "filter primers for low variation", main_primers);