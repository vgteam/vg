#include "subcommand.hpp"
#include "../utility.hpp"
#include "../viz.hpp"
#include "../stream.hpp"

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
         << "    -C, --hide-cnv        suppress the visualization of CNVs in paths" << endl;
}

int main_viz(int argc, char** argv) {

    string xg_name;
    vector<string> packs_in;
    vector<string> names;
    string image_out;
    int image_height = 1024;
    int image_width = 0;
    double scale_factor = 1.0;
    bool show_cnv = true;
    
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
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:i:n:o:X:Y:s:C",
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
            names.push_back(optarg);
            break;
        case 'i':
            packs_in.push_back(optarg);
            break;
        case 'o':
            image_out = optarg;
            break;
        case 'X':
            image_width = atoi(optarg);
            break;
        case 'Y':
            image_height = atoi(optarg);
            break;
        case 's':
            scale_factor = atof(optarg);
            break;
        case 'C':
            show_cnv = false;
            break;
        default:
            abort();
        }
    }

    xg::XG xgidx;
    if (xg_name.empty()) {
        cerr << "No XG index given. An XG index must be provided." << endl;
        exit(1);
    } else {
        ifstream in(xg_name.c_str());
        xgidx.load(in);
    }

    // todo one packer per thread and merge
    vector<vg::Packer> packs;
    for (auto& f : packs_in) {
        packs.emplace_back();
        auto& p = packs.back();
        p.load_from_file(f);
    }

    Viz viz(&xgidx, &packs, image_out, image_width, image_height, show_cnv);
    viz.draw();

    return 0;
}

// Register subcommand
static Subcommand vg_viz("viz", "render visualizations of indexed graphs and read sets", main_viz);
