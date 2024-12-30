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
         << "    -x, --xg-path FILE                 use this xg graph" << endl
         << "    -s, --snarl-index FILE             use this snarl index" << endl
         << "    -r, --r-index FILE                 use this r index" << endl
         << "    -g, --gbz FILE                     use this gbz file" << endl
         << "    -v, --variation-threshold DOUBLE   output primers that work for at least this percentage of haplotypes (default: 0.8)"
         << "    -l, --tolerance INT                allow this much difference between minimum and maximum sizes compared to the linear product size (default: 10)" << endl
         << "    -n, --minimum-size INT             minimum product size allowed (has precedence over --tolerance)" << endl
         << "    -m, --maximum-size INT             maximum product size allowed (has precedence over --tolerance)" << endl
         << "    -a, --all-primers                  output all primers" << endl;
}

size_t difference(const size_t& a, const size_t& b) {
    size_t diff;
    if (a == b) {
        diff = 0;
    } else if (a > b) {
        diff = a - b;
    } else {
        diff = b - a;
    }
    return diff;
}

void print_tabular(const string& genome_name, const PrimerPair& primer_pair) {
    const Primer& left_primer  = primer_pair.left_primer;
    const Primer& right_primer = primer_pair.right_primer;
    size_t rln = right_primer.mapped_nodes_ids.size() - 1; //right primer last node index
    cout << genome_name                          << "\t"
         << primer_pair.template_feature         << "\t"
         << primer_pair.template_position        << "\t"
         << left_primer.sequence                 << "\t" 
         << right_primer.sequence                << "\t"
         << left_primer.position_template        << "\t"
         << right_primer.position_template       << "\t"
         << left_primer.position_chromosome      << "\t"
         << right_primer.position_chromosome     << "\t"
         << left_primer.mapped_nodes_ids[0]      << "\t"
         << right_primer.mapped_nodes_ids[rln]   << "\t"
         << left_primer.length                   << "\t"
         << right_primer.length                  << "\t"
         << primer_pair.linear_product_size      << "\t"
         << primer_pair.min_product_size         << "\t"
         << primer_pair.max_product_size         << "\t"
         << primer_pair.variation_level          << endl;
}

int main_primers(int argc, char** argv) {
    
    if (argc == 2) {
        help_primers(argv);
        return 1;
    }

    string xg_path;
    string snarl_index_path;
    string ri_path;
    string gbz_path;
    bool zero_variation = false;
    bool all_primers   = false;
    int  tolerance     = 10;
    double variation_threshold = 0.8;
    int  minimum_product_size = numeric_limits<int>::max();
    int  maximum_product_size = numeric_limits<int>::max();

    int c;
    optind = 2;

    while (true) {
        static struct option long_options[] =
        {
          {"help",                no_argument,       0, 'h'},
          {"xg-path",             required_argument, 0, 'x'},
          {"snarl-index",         required_argument, 0, 's'},
          {"ri-path",             required_argument, 0, 'r'},
          {"gbz-path",            required_argument, 0, 'g'},
          {"variation-threshold", required_argument, 0, 'v'},
          {"tolerance",           required_argument, 0, 'l'},
          {"minimum-size",        required_argument, 0, 'n'},
          {"maximum-size",        required_argument, 0, 'm'},
          {"all-primers",         required_argument, 0, 'a'},
          {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:s:r:g:l:n:m:a", long_options, &option_index);

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
        
        case 'r':
            ri_path = optarg;
            break;
        
        case 'g':
            gbz_path = optarg;
            break;
        
        case 'v':
            variation_threshold = parse<double>(optarg);
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
        
        case 'a':
            all_primers = true;
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
        cerr << "error:[vg primers] snarl index file (-s) is required" << endl;
        exit(1);
    }

    if (ri_path.empty()) {
        cerr << "error:[vg primers] r index file (-r) is required" << endl;
        exit(1);
    }

    if (gbz_path.empty()) {
        cerr << "error:[vg primers] gbz file (-g) is required" << endl;
        exit(1);
    }

    string primers_path = get_input_file_name(optind, argc, argv);
    
    SnarlDistanceIndex distance_index;
    unique_ptr<handlegraph::PathPositionHandleGraph> graph;
    gbwtgraph::GBWTGraph gbwt_graph;
    gbwt::GBWT gbwt_index;
    gbwt::FastLocate r_index;
    load_r_index(r_index, ri_path);
    load_gbz(gbwt_index, gbwt_graph, gbz_path);
    gbwt_graph.set_gbwt(gbwt_index);
    r_index.setGBWT(gbwt_index);
    distance_index.deserialize(snarl_index_path);
    graph = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_path);
    ifstream file_handle(primers_path);
    PrimerFinder primer_finder(graph, &distance_index, file_handle, gbwt_graph, gbwt_index, r_index);

    cout << "chrom\ttplfeat\ttplpos\tlpseq\trpseq\tlppostpl\trppostmp\tlpposchrom\trpposchrom\t"
         << "lpnid\trpnid\tlplen\trplen\tlinsize\tminsize\tmaxsize\tvarlevel" << endl;
         
    vector<string> reference_paths = primer_finder.get_reference_paths();
    for (size_t i = 0; i < reference_paths.size(); ++i) {
        string path_name = reference_paths[i];
        const vector<PrimerPair>& primer_pairs = primer_finder.get_primer_pairs_of_chrom(path_name);
        for (size_t j = 0; j < primer_pairs.size(); ++j) {
            const PrimerPair& primer_pair = primer_pairs[j];
            if (all_primers) {
                print_tabular(path_name, primer_pair);
            } else {
                if (minimum_product_size != numeric_limits<int>::max() &&
                    primer_pair.min_product_size < minimum_product_size) {
                    continue;
                }
                if (maximum_product_size != numeric_limits<int>::max() &&
                    primer_pair.max_product_size > maximum_product_size) {
                    continue;
                }

                if (difference(primer_pair.linear_product_size, primer_pair.min_product_size) > tolerance
                || difference(primer_pair.linear_product_size, primer_pair.max_product_size) > tolerance) {
                    continue; 
                }

                if (primer_pair.variation_level < variation_threshold) {
                    continue;
                }
                print_tabular(path_name, primer_pair);
            }
        }
    }

    return 0;
}

static Subcommand vg_primers("primers", "filter primers for low variation", main_primers);
