#include "subcommand.hpp"
#include "../io/save_handle_graph.hpp"
#include "../handle.hpp"
#include "../xg.hpp"
#include "../utility.hpp"
#include "../snarls.hpp"
#include "../masker.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <gbwtgraph/gbz.h>

#include <unistd.h>
#include <getopt.h>
#include <sstream>

using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;
using namespace std;

vector<tuple<string, size_t, size_t>> parse_bed(istream& in){
    
    vector<tuple<string, size_t, size_t>> regions;
    
    string line_buffer;
    string token_buffer;
    while (getline(in, line_buffer)) {
        
        if (line_buffer.empty() || line_buffer.front() == '#' || line_buffer.substr(0, 5) == "track" || line_buffer.substr(0, 7) == "browser") {
            // header line
            continue;
        }
        istringstream line_strm(line_buffer);
        
        regions.emplace_back();
        
        getline(line_strm, token_buffer, '\t');
        get<0>(regions.back()) = token_buffer;
        getline(line_strm, token_buffer, '\t');
        get<1>(regions.back()) = parse<int64_t>(token_buffer);
        getline(line_strm, token_buffer, '\t');
        get<2>(regions.back()) = parse<int64_t>(token_buffer);
    }
    
    return regions;
}

void help_mask(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph>" << endl
       << "Mask out specified regions of a graph with N's" << endl
       << endl
       << "input options: " << endl
       << "    -b, --bed FILE       BED regions corresponding to path intervals of the graph to target (required)" << endl
       << "    -g, --gbz-input      Input graph is in GBZ format" << endl
       << "    -s, --snarls FILE    Snarls from vg snarls (computed directly if not provided)" << endl
    
       << endl;
}    

int main_mask(int argc, char** argv) {
    
    string bed_filepath;
    string snarls_filepath;
    bool gbz_input = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"bed", required_argument, 0, 'b'},
            {"gbz-input", no_argument, 0, 'g'},
            {"snarls", required_argument, 0, 's'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hb:gs:", long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case '?':
            case 'h':
                help_mask(argv);
                return 0;
            case 'b':
                bed_filepath = optarg;
                break;
            case 'g':
                gbz_input = true;
                break;
            case 's':
                snarls_filepath = optarg;
                break;
            default:
                help_mask(argv);
                return 1;
        }
    }
    
    
    if (argc - optind != 1) {
        cerr << "error: vg mask requires exactly 1 positional argument" << endl;
        help_mask(argv);
        return 1;
    }

    if (bed_filepath.empty()) {
        cerr << "error: vg mask requires an input BED file from -b / --bed" << endl;
        return 1;
    }
    if (bed_filepath != "-") {
        if (!ifstream(bed_filepath)) {
            cerr << "error: could not open BED file from " << bed_filepath << endl;
            return 1;
        }
    }
    if (!snarls_filepath.empty() && snarls_filepath != "-") {
        if (!ifstream(snarls_filepath)) {
            cerr << "error: could not open snarls file from " << snarls_filepath << endl;
            return 1;
        }
    }
    
    // load the graph
    string graph_path = get_input_file_name(optind, argc, argv);
    unique_ptr<MutablePathDeletableHandleGraph> graph;
    unique_ptr<gbwtgraph::GBZ> gbz;
    if (gbz_input) {
        gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(graph_path);;
    }
    else {
        graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(graph_path);
    }
    
    // load the BED
    vector<tuple<string, size_t, size_t>> regions;
    get_input_file(bed_filepath, [&](istream& in) {
        regions = std::move(parse_bed(in));
    });

    // load the snarls
    unique_ptr<SnarlManager> snarl_manager;
    if (!snarls_filepath.empty()) {
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarls_filepath);
    }
    
    Masker masker;
    if (gbz.get()) {
        masker = std::move(Masker(gbz->graph, snarl_manager.get()));
    }
    else {
        masker = std::move(Masker(*graph, snarl_manager.get()));
    }
    
    masker.mask_sequences(regions);

    // write the graph
    if (gbz.get()) {
        gbz->simple_sds_serialize(cout);
    }
    else {
        vg::io::save_handle_graph(graph.get(), cout);
    }
    
    return 0;
}


// Register subcommand
static Subcommand vg_mask("mask", "Mask out sequences in a graph with N's", main_mask);
