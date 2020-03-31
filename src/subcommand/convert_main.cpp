#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "xg.hpp"
#include "../convert_handle.hpp"
#include "../io/save_handle_graph.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_convert(char** argv) {
    cerr << "usage: " << argv[0] << " convert [options] <input-graph>" << endl
         << "output options:" << endl
         << "    -v, --vg-out           output in VG format [default]" << endl
         << "    -a, --hash-out         output in HashGraph format" << endl
         << "    -p, --packed-out       output in PackedGraph format" << endl
         << "    -x, --xg-out           output in XG format" << endl
         << "    -o, --odgi-out         output in ODGI format" << endl;
}

int main_convert(int argc, char** argv) {

    string output_format = "vg";

    if (argc == 2) {
        help_convert(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"vg-out", no_argument, 0, 'v'},
            {"hash-out", no_argument, 0, 'a'},
            {"packed-out", no_argument, 0, 'p'},
            {"xg-out", no_argument, 0, 'x'},
            {"odgi-out", no_argument, 0, 'o'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hvxapxo",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_convert(argv);
            return 1;
        case 'v':
            output_format = "vg";
            break;
        case 'a':
            output_format = "hash";
            break;
        case 'p':
            output_format = "packed";
            break;
        case 'x':
            output_format = "xg";
            break;
        case 'o':
            output_format = "odgi";
            break;

        default:
            abort();
        }
    }

    unique_ptr<HandleGraph> input_graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        input_graph = vg::io::VPKG::load_one<HandleGraph>(in);
      });

    // allocate a graph using the graph_type string to decide a class
    unique_ptr<HandleGraph> output_graph;
    if (output_format == "vg") {
      output_graph = unique_ptr<HandleGraph>(new VG());
    } else if (output_format == "hash") {
      output_graph = unique_ptr<HandleGraph>(new bdsg::HashGraph());
    } else if (output_format == "packed") {
      output_graph = unique_ptr<HandleGraph>(new bdsg::PackedGraph());
    } else if (output_format == "xg") {
      output_graph = unique_ptr<HandleGraph>(new xg::XG());
    } else if (output_format == "odgi") {
      output_graph = unique_ptr<HandleGraph>(new bdsg::ODGI());
    }

    PathHandleGraph* input_path_graph = dynamic_cast<PathHandleGraph*>(input_graph.get());
    PathHandleGraph* output_path_graph = dynamic_cast<PathHandleGraph*>(output_graph.get());

    if (output_format == "xg") {
        if (input_path_graph != nullptr) {
            dynamic_cast<xg::XG*>(output_graph.get())->from_path_handle_graph(*input_path_graph);
        } else {
            // not that this will happen any time soon, but we can relax it once xg0 is in,
            // as it has a from_handle_graph
            cerr << "error [vg convert]: cannot output to xg as input graph does not support paths" << endl;
            return 1;
        }
    }
    else if (input_path_graph != nullptr && output_path_graph != nullptr) {
        MutablePathMutableHandleGraph* mutable_output_graph = dynamic_cast<MutablePathMutableHandleGraph*>(output_path_graph);
        assert(mutable_output_graph != nullptr);
        convert_path_handle_graph(input_path_graph, mutable_output_graph);
    } else if (input_path_graph != nullptr) {
        MutableHandleGraph* mutable_output_graph = dynamic_cast<MutableHandleGraph*>(output_graph.get());
        assert(mutable_output_graph != nullptr);
        cerr << "warning [vg convert]: output format does not support paths" << endl;
        convert_handle_graph(input_graph.get(), mutable_output_graph);
    }

    // Serialize the graph using VPKG.
    vg::io::save_handle_graph(output_graph.get(), cout);

    return 0;
}

// Register subcommand
static Subcommand vg_convert("convert", "convert graphs between handle-graph compiant formats", main_convert);
