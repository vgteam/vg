#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../xg.hpp"
#include "../convert_handle.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

#include "sglib/packed_graph.hpp"
#include "sglib/hash_graph.hpp"

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_convert(char** argv) {
    cerr << "usage: " << argv[0] << " convert [options] <input-graph>" << endl
         << "input options:" << endl
         << "    -v, --vg-in            input in VG format [default]" << endl
         << "    -a, --hash-in          input in HashGraph format" << endl
         << "    -p, --packed-in        input in PackedGraph format" << endl
         << "    -x, --xg-in            input in XG format" << endl
         << "output options:" << endl
         << "    -V, --vg-out           output in VG format [default]" << endl
         << "    -A, --hash-out         output in HashGraph format" << endl
         << "    -P, --packed-out       output in PackedGraph format" << endl;
}

int main_convert(int argc, char** argv) {

    string input_format = "vg";
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
            {"vg-in", no_argument, 0, 'v'},
            {"hash-in", no_argument, 0, 'a'},
            {"packed-in", no_argument, 0, 'p'},
            {"xg-in", no_argument, 0, 'x'},
            {"vg-out", no_argument, 0, 'V'},
            {"hash-out", no_argument, 0, 'A'},
            {"packed-out", no_argument, 0, 'P'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hvxapxVAP",
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
            input_format = "vg";
            break;
        case 'a':
            input_format = "hash";
            break;
        case 'p':
            input_format = "packed";
            break;
        case 'x':
            input_format = "xg";
            break;            
        case 'V':
            output_format = "vg";
            break;
        case 'A':
            output_format = "hash";
            break;
        case 'P':
            output_format = "packed";
            break;

        default:
            abort();
        }
    }

    string input_path = get_input_file_name(optind, argc, argv);

    // allocate a graph using the graph_type string to decide a class
    function<PathHandleGraph*(const string&)> graph_factory = [&](const string& graph_type) -> PathHandleGraph*{
        if (graph_type == "vg") {
            return new VG();
        } else if (graph_type == "hash") {
            return new sglib::HashGraph();
        } else if (graph_type == "packed") {
            return new sglib::PackedGraph();
        } else if (graph_type == "xg") {
            return new XG();
        }
        return nullptr;
    };

    ifstream input_stream;
    if (input_path != "-") {
        input_stream.open(input_path);
    }

    // Make an empty input graph
    PathHandleGraph* input_graph = graph_factory(input_format);
    assert(input_graph != nullptr);
    // Read the graph from the stream
    if (input_format != "xg") {
        dynamic_cast<SerializableHandleGraph*>(input_graph)->deserialize(input_path == "-" ? cin : input_stream);
    } else {
        //todo: XG::deserialize() doesn't work.  Need to go through vpkg
        unique_ptr<XG> xindex = vg::io::VPKG::load_one<XG>(input_path);
        delete input_graph;
        input_graph = xindex.release();
    }
    

    // Make an empty output graph
    MutablePathMutableHandleGraph* output_graph = dynamic_cast<MutablePathMutableHandleGraph*>(graph_factory(output_format));
    assert(output_graph != nullptr);
    // Copy over the input graph
    convert_path_handle_graph(input_graph, output_graph);
    // Write the output graph to the stream
    dynamic_cast<SerializableHandleGraph*>(output_graph)->serialize(cout);

    delete input_graph;
    delete output_graph;

    return 0;
}

// Register subcommand
static Subcommand vg_convert("convert", "convert graphs between handle-graph compiant formats", main_convert);
