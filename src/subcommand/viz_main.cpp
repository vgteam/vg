#include "subcommand.hpp"
#include "../utility.hpp"
#include "../viz.hpp"
#include "../xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_viz(char** argv) {
    cerr << "usage: " << argv[0] << " viz [options]" << endl
         << "options:" << endl
         << "    -x, --xg FILE         use this basis graph" << endl
         << "    -i, --pack-in FILE    use this compressed coverage format (multiple allowed)" << endl
         << "    -n, --name NAME       apply name to the previous .pack (multiple allowed)" << endl
         << "    -o, --out FILE        write to file (could be .png or .svg)" << endl
         << "    -X, --width N         write an image N pixels wide (default 1024)" << endl
         << "    -Y, --height N        write an image N pixels high (default 1024)" << endl
         << "    -C, --show-cnv        visualize CNVs in paths on new rows (default uses text)" << endl
         << "    -P, --hide-paths      hide reference paths in the graph" << endl
         << "    -D, --hide-dna        suppress the visualization of DNA sequences" << endl;
}

int main_viz(int argc, char** argv) {

    string xg_name;
    vector<string> packs_in;
    vector<string> pack_names;
    string image_out;
    int image_height = 0;
    int image_width = 0;
    double scale_factor = 1.0;
    bool show_cnv = false;
    bool show_dna = true;
    bool show_paths = true;
    
    if (argc == 2) {
        help_viz(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg", required_argument,0, 'x'},
            {"pack-in", required_argument,0, 'i'},
            {"name", required_argument, 0, 'n'},
            {"width", required_argument, 0, 'X'},
            {"height", required_argument, 0, 'Y'},
            {"out", required_argument, 0, 'o'},
            {"scale", required_argument, 0, 's'},
            {"hide-cnv", no_argument, 0, 'C'},
            {"hide-dna", no_argument, 0, 'D'},
            {"hide-paths", no_argument, 0, 'P'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:i:n:o:X:Y:s:CDP",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_viz(argv);
            return 1;
        case 'x':
            xg_name = optarg;
            break;
        case 'n':
            pack_names.push_back(optarg);
            break;
        case 'i':
            packs_in.push_back(optarg);
            break;
        case 'o':
            image_out = optarg;
            break;
        case 'X':
            image_width = parse<int>(optarg);
            break;
        case 'Y':
            image_height = parse<int>(optarg);
            break;
        case 's':
            scale_factor = parse<double>(optarg);
            break;
        case 'C':
            show_cnv = true;
            break;
        case 'D':
            show_dna = false;
            break;
        case 'P':
            show_paths = false;
            break;
        default:
            abort();
        }
    }

    PathPositionHandleGraph* xgidx = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionVectorizableOverlayHelper overlay_helper;
    if (xg_name.empty()) {
        cerr << "No input graph given. An input graph (-x) must be provided." << endl;
        exit(1);
    } else {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        // We know the PathPositionVectorizableOverlayHelper produces a PathPositionVectorizableOverlay which implements PathPositionHandleGraph.
        // TODO: Make the types actually work out here.
        xgidx = dynamic_cast<PathPositionHandleGraph*>(overlay_helper.apply(path_handle_graph.get()));
        assert(xgidx != nullptr);
    }

    // todo one packer per thread and merge
    vector<vg::Packer> packs;
    for (auto& f : packs_in) {
        packs.emplace_back();
        auto& p = packs.back();
        p.load_from_file(f);
    }

    // default to using the file names as the row names for the packs
    if (pack_names.empty()) {
        pack_names = packs_in;
    }

    {
        Viz viz(xgidx, &packs, pack_names, image_out, image_width, image_height, show_cnv, show_dna, show_paths);
        viz.draw();
    }

    return 0;
}

// Register subcommand
static Subcommand vg_viz("viz", "render visualizations of indexed graphs and read sets", main_viz);
