// inject_main.cpp: define the "vg inject" subcommand, which lifts over alignments from the linear space

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include <subcommand.hpp>

#include "../utility.hpp"
#include "../alignment.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../hts_alignment_emitter.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_inject(char** argv) {
    cerr << "usage: " << argv[0] << " inject -x graph.xg [options] input.[bam|sam|cram] >output.gam" << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE       use this graph or xg index (required, non-XG formats also accepted)" << endl
         << "    -o, --output-format NAME output the alignments in NAME format (gam / gaf / json) [gam]" << endl
         << "    -t, --threads N          number of threads to use" << endl;
}

int main_inject(int argc, char** argv) {

    if (argc == 2) {
        help_inject(argv);
    }

    string xg_name;
    string output_format = "GAM";
    std::set<std::string> output_formats = { "GAM", "GAF", "JSON" };
    int threads = get_thread_count();

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
          {"help", no_argument, 0, 'h'},
          {"xg-name", required_argument, 0, 'x'},
          {"output-format", required_argument, 0, 'o'},
          {"threads", required_argument, 0, 't'},
          {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:o:t:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;
        
        case 'o':
            {
                output_format = optarg;
                for (char& c : output_format) {
                    c = std::toupper(c);
                }
                if (output_formats.find(output_format) == output_formats.end()) {
                    std::cerr << "error: [vg inject] Invalid output format: " << optarg << std::endl;
                    std::exit(1);
                }
            }
            break;

        case 't':
          threads = parse<int>(optarg);
          break;

        case 'h':
        case '?':
          help_inject(argv);
          exit(1);
          break;

        default:
          abort ();
        }
    }
    
    omp_set_num_threads(threads);

    string file_name = get_input_file_name(optind, argc, argv);

    // We require an XG index
    if (xg_name.empty()) {
        cerr << "error[vg inject]: XG index (-x) is required" << endl;
        exit(1);
    }
    unique_ptr<PathHandleGraph> path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* xgidx = overlay_helper.apply(path_handle_graph.get());

    // We don't do HTS output formats but we do need an empty paths collection to make an alignment emitter
    vector<tuple<path_handle_t, size_t, size_t>> paths;
    unique_ptr<vg::io::AlignmentEmitter> alignment_emitter = get_alignment_emitter("-", output_format, paths, threads, xgidx);

    function<void(Alignment&)> lambda = [&](Alignment& aln) {
        alignment_emitter->emit_mapped_single({std::move(aln)});
    };
    if (threads > 1) {
        hts_for_each_parallel(file_name, lambda, xgidx);
    } else {
        hts_for_each(file_name, lambda, xgidx);
    }
    return 0;
}

// Register subcommand
static Subcommand vg_inject("inject", "lift over alignments for the graph", main_inject);
