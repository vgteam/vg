/** \file version_main.cpp
 *
 * Defines the "vg version" subcommand, which evaluates graphs and alignments.
 */


#include <omp.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <getopt.h>


#include "subcommand.hpp"

#include "sdsl/bit_vectors.hpp"
#include <vg/vg.pb.h>
#include "../version.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../xg.hpp"
#include "../region.hpp"
#include "../convert_handle.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_xg(char** argv) {
    cerr << "usage: " << argv[0] << " xg [options]" << endl
         << "Manipluate succinct representations of queryable sequence graphs" << endl
         << endl
         << "options:" << endl
         << "    -v, --vg FILE              compress graph in vg FILE" << endl
         << "    -V, --validate             validate compression" << endl
         << "    -o, --out FILE             serialize graph to FILE in xg format" << endl
         << "    -i, --in FILE              use index in FILE" << endl
         << "    -X, --extract-vg FILE      serialize graph to FILE in vg format" << endl
         << "    -n, --node ID              graph neighborhood around node with ID" << endl
         << "    -c, --context N            steps of context to extract when building neighborhood" << endl
         << "    -s, --node-seq ID          provide node sequence for ID" << endl
         << "    -P, --char POS             give the character at a given position in the graph" << endl
         << "    -F, --substr POS:LEN       extract the substr of LEN on the node at the position" << endl
         << "    -f, --edges-from ID        list edges from node with ID" << endl
         << "    -t, --edges-to ID          list edges to node with ID" << endl
         << "    -O, --edges-of ID          list all edges related to node with ID" << endl
         << "    -S, --edges-on-start ID    list all edges on start of node with ID" << endl
         << "    -E, --edges-on-end ID      list all edges on start of node with ID" << endl
         << "    -p, --path TARGET          gets the region of the graph @ TARGET (chr:start-end)" << endl
         << "    -x, --extract-threads      extract succinct threads as paths" << endl
         << "    -r, --store-threads        store perfect match paths as succinct threads" << endl
         << "    -d, --is-sorted-dag        graph is a sorted dag; use fast thread insert" << endl
         << "    -R, --report FILE          save an HTML space usage report to FILE when serializing" << endl
         << "    -D, --debug                show debugging output" << endl
         << "    -T, --text-output          write text instead of vg protobuf" << endl
         << "    -b, --dump-bs FILE         dump the gPBWT to the given file" << endl
         << "    -h, --help                 this text" << endl;
}

int main_xg(int argc, char** argv) {

    if (argc == 2) {
        help_xg(argv);
        return 1;
    }

    string vg_in;
    string vg_out;
    string out_name;
    string in_name;
    int64_t node_id;
    bool edges_from = false;
    bool edges_to = false;
    bool edges_of = false;
    bool edges_on_start = false;
    bool edges_on_end = false;
    bool node_sequence = false;
    string pos_for_char;
    string pos_for_substr;
    int context_steps = 0;
    bool node_context = false;
    string target;
    bool print_graph = false;
    bool text_output = false;
    bool validate_graph = false;
    bool extract_threads = false;
    bool store_threads = false;
    bool is_sorted_dag = false;
    string report_name;
    string b_array_name;
    
    int c;
    optind = 2; // force optind past "xg" positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"vg", required_argument, 0, 'v'},
                {"out", required_argument, 0, 'o'},
                {"in", required_argument, 0, 'i'},
                {"extract-vg", required_argument, 0, 'X'},
                {"node", required_argument, 0, 'n'},
                {"char", required_argument, 0, 'P'},
                {"substr", required_argument, 0, 'F'},
                //{"range", required_argument, 0, 'r'},
                {"context", required_argument, 0, 'c'},
                {"edges-from", required_argument, 0, 'f'},
                {"edges-to", required_argument, 0, 't'},
                {"edges-of", required_argument, 0, 'O'},
                {"edges-on-start", required_argument, 0, 'S'},
                {"edges-on-end", required_argument, 0, 'E'},
                {"node-seq", required_argument, 0, 's'},
                {"path", required_argument, 0, 'p'},
                {"extract-threads", no_argument, 0, 'x'},
                {"store-threads", no_argument, 0, 'r'},
                {"is-sorted-dag", no_argument, 0, 'd'},
                {"report", required_argument, 0, 'R'},
                {"debug", no_argument, 0, 'D'},
                {"text-output", no_argument, 0, 'T'},
                {"validate", no_argument, 0, 'V'},
                {"dump-bs", required_argument, 0, 'b'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hv:o:i:X:f:t:s:c:n:p:DxrdTO:S:E:VR:P:F:b:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'v':
            vg_in = optarg;
            break;

        case 'V':
            validate_graph = true;
            break;

        case 'o':
            out_name = optarg;
            break;

        case 'D':
            print_graph = true;
            break;

        case 'T':
            text_output = true;
            break;
            
        case 'x':
            extract_threads = true;
            break;
            
        case 'r':
            store_threads = true;
            break;
            
        case 'd':
            is_sorted_dag = true;
            break;

        case 'i':
            in_name = optarg;
            break;

        case 'X':
            vg_out = optarg;
            break;

        case 'n':
            node_id = parse<int64_t>(optarg);
            node_context = true;
            break;

        case 'c':
            context_steps = parse<int>(optarg);
            break;

        case 'f':
            node_id = parse<int64_t>(optarg);
            edges_from = true;
            break;
            
        case 't':
            node_id = parse<int64_t>(optarg);
            edges_to = true;
            break;

        case 'O':
            node_id = parse<int64_t>(optarg);
            edges_of = true;
            break;

        case 'S':
            node_id = parse<int64_t>(optarg);
            edges_on_start = true;
            break;

        case 'E':
            node_id = parse<int64_t>(optarg);
            edges_on_end = true;
            break;

        case 's':
            node_id = parse<int64_t>(optarg);
            node_sequence = true;
            break;

        case 'p':
            target = optarg;
            break;

        case 'P':
            pos_for_char = optarg;
            break;
            
        case 'F':
            pos_for_substr = optarg;
            break;
            
        case 'R':
            report_name = optarg;
            break;
            
        case 'b':
            b_array_name = optarg;
            break;
            
        case 'h':
        case '?':
            help_xg(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    unique_ptr<XG> graph;
    //string file_name = argv[optind];
    if (in_name.empty()) assert(!vg_in.empty());
    if (vg_in == "-") {
        // Read VG from stdin
        graph = unique_ptr<XG>(new XG());
        graph->from_stream(std::cin, validate_graph, print_graph, store_threads, is_sorted_dag);
    } else if (vg_in.size()) {
        // Read VG from a file
        ifstream in;
        in.open(vg_in.c_str());
        graph = unique_ptr<XG>(new XG());
        graph->from_stream(in, validate_graph, print_graph, store_threads, is_sorted_dag);
    }

    if (in_name.size()) {
        get_input_file(in_name, [&](istream& in) {
            // Load from an XG file or - (stdin)
            graph = vg::io::VPKG::load_one<XG>(in);
        });
    }

    // Prepare structure tree for serialization
    unique_ptr<sdsl::structure_tree_node> structure;
    
    if (!report_name.empty()) {
        // We need to make a report, so we need the structure. Make a real tree
        // node. The unique_ptr handles deleting.
        structure = unique_ptr<sdsl::structure_tree_node>(new sdsl::structure_tree_node("name", "type"));
    }

    if(!vg_out.empty()) {
        if (graph.get() == nullptr) {
             cerr << "error [vg xg] no xg graph exists to convert; Try: vg xg -i graph.xg -X graph.vg" << endl;
             return 1;
        }
        
        VG converted;
        // Convert the xg graph to vg format
        convert_path_handle_graph(graph.get(), &converted);
        
        if (vg_out == "-") {
            converted.serialize_to_ostream(std::cout);
        } else {
            converted.serialize_to_file(vg_out);
        }
    }

    if (!out_name.empty()) {
        // Open a destination file if it is a file we want to write to
        ofstream out_file;
        if (out_name != "-") {
            out_file.open(out_name);
        }
        // Work out where to save to
        ostream& out = (out_name == "-") ? std::cout : out_file;
        
        // Encapsulate output in VPKG
        vg::io::VPKG::with_save_stream(out, "XG", [&](ostream& tagged) {
            // Serialize to the file while recording space usage to the structure.
            graph->serialize_and_measure(tagged, structure.get(), "xg");
        });
        
        out.flush();
    }

    if (!report_name.empty()) {
        // Save the report
        ofstream out;
        out.open(report_name.c_str());
        sdsl::write_structure_tree<HTML_FORMAT>(structure.get(), out, 0);
    }

    // queries
    if (node_sequence) {
        cout << node_id << ": " << graph->node_sequence(node_id) << endl;
    }
    if (!pos_for_char.empty()) {
        // extract the position from the string
        int64_t id;
        bool is_rev;
        size_t off;
        extract_pos(pos_for_char, id, is_rev, off);
        // then pick it up from the graph
        cout << graph->pos_char(id, is_rev, off) << endl;
    }
    if (!pos_for_substr.empty()) {
        int64_t id;
        bool is_rev;
        size_t off;
        size_t len;
        extract_pos_substr(pos_for_substr, id, is_rev, off, len);
        cout << graph->pos_substr(id, is_rev, off, len) << endl;
    }
    
    if (edges_from) {
        vector<Edge> edges = graph->edges_from(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_to) {
        vector<Edge> edges = graph->edges_to(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_of) {
        vector<Edge> edges = graph->edges_of(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_on_start) {
        vector<Edge> edges = graph->edges_on_start(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_on_end) {
        vector<Edge> edges = graph->edges_on_end(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }

    if (node_context) {
        Graph g;
        graph->neighborhood(node_id, context_steps, g);
        if (text_output) {
            to_text(cout, g);
        } else {
            vector<Graph> gb = { g };
            vg::io::write_buffered(cout, gb, 0);
        }
    }

    if (!target.empty()) {
        string name;
        int64_t start, end;
        Graph g;
        parse_region(target, name, start, end);
        graph->get_path_range(name, start, end, g);
        graph->expand_context(g, context_steps);
        if (text_output) {
            to_text(cout, g);
        } else {
            vector<Graph> gb = { g };
            vg::io::write_buffered(cout, gb, 0);
        }
    }
    
    if (extract_threads) {
        list<XG::thread_t> threads;
        for (auto& p : graph->extract_threads(false)) {
            for (auto& t : p.second) {
                threads.push_back(t);
            }
        }
        for (auto& p : graph->extract_threads(true)) {
            for (auto& t : p.second) {
                threads.push_back(t);
            }
        }

        size_t thread_number = 0;
        for(XG::thread_t& thread : threads) {
            // Convert to a Path
            Path path;
            for(XG::ThreadMapping& m : thread) {
                // Convert all the mappings
                Mapping mapping;
                mapping.mutable_position()->set_node_id(m.node_id);
                mapping.mutable_position()->set_is_reverse(m.is_reverse);
                
                *(path.add_mapping()) = mapping;
            }
        
        
            // Give each thread a name
            path.set_name("_thread_" + to_string(thread_number++));
            
            // We need a Graph for serialization purposes. We do one chunk per
            // thread in case the threads are long.
            Graph g;
            
            *(g.add_path()) = path;
            
            // Dump the graph with its mappings. TODO: can we restrict these to
            // mappings to nodes we have already pulled out? Or pull out the
            // whole compressed graph?
            if (text_output) {
                to_text(cout, g);
            } else {
                vector<Graph> gb = { g };
                vg::io::write_buffered(cout, gb, 0);
            }
            
        }
    }

    if (!b_array_name.empty()) {
        // Dump B array
        ofstream out;
        out.open(b_array_name.c_str());
        graph->bs_dump(out);
    }

    return 0;
}

// Register subcommand
static Subcommand vg_xg("xg", "manipulate xg files", main_xg);
