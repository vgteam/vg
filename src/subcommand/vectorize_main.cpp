/** \file vectorize_main.cpp
 *
 * Defines the "vg vectorize" subcommand, which turns graphs into ML vectors.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../vectorizer.hpp"
#include "../mapper.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

const string context = "[vg vectorize]";

void help_vectorize(char** argv) {
    cerr << "usage: " << argv[0] << " vectorize [options] -x <index.xg> <alignments.gam>" << endl
         << "Vectorize a set of alignments to a variety of vector formats." << endl
         << endl
         << "options: " << endl
         << "  -x, --xg FILE              an xg index or graph of interest" << endl
         << "  -g, --gcsa FILE            a GCSA2 index to use if generating MEM sketches" << endl
         << "  -l, --aln-label LABEL      output all alignments with name LABEL" << endl
         << "  -f, --format               tab-delimit output so it can be used in R." << endl
         << "  -A, --annotate             create a header with each node/edge's name" << endl
         << "                             and a column with alignment names." << endl
         << "  -a, --a-hot                instead of a 1-hot, output a vector of {0|1|2}" << endl
         << "                             for covered, reference, or alt." << endl
         << "  -w, --wabbit               output a format that's friendly to vowpal wabbit" << endl
         << "  -M, --wabbit-mapping FILE  output the mappings used for vowpal wabbit classes" << endl
         << "                             (default: print to stderr)" << endl
         << "  -m, --mem-sketch           generate a MEM sketch of a given read based on GCSA" << endl
         << "  -p, --mem-positions        add the positions to the MEM sketch of a given read" << endl
         << "                             based on the GCSA" << endl
         << "  -H, --mem-hit-max N        ignore MEMs with more than this many hits when" << endl
         << "                             extracting poisitions" << endl
         << "  -i, --identity-hot         output a score vector for each alignment based on" << endl
         << "                             percent identity and coverage" << endl
         << "  -h, --help                 print this help message to stderr and exit" << endl;
}

int main_vectorize(int argc, char** argv){

    string xg_name;
    string aln_label = "";
    string gcsa_name;
    string wabbit_mapping_file = "";
    bool format = false;
    bool show_header = false;
    bool map_alns = false;
    bool annotate = false;
    bool a_hot = false;
    bool output_wabbit = false;
    bool use_identity_hot = false;
    bool mem_sketch = false;
    bool mem_positions = false;
    bool mem_hit_max = 0;
    int max_mem_length = 0;

    if (argc <= 2) {
        help_vectorize(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"annotate", no_argument, 0, 'A'},
            {"xg", required_argument,0, 'x'},
            {"gcsa", required_argument,0, 'g'},
            {"format", no_argument, 0, 'f'},
            {"a-hot", no_argument, 0, 'a'},
            {"wabbit", no_argument, 0, 'w'},
            {"wabbit-mapping", required_argument, 0, 'M'},
            {"mem-sketch", no_argument, 0, 'm'},
            {"mem-positions", no_argument, 0, 'p'},
            {"mem-hit-max", required_argument, 0, 'H'},
            {"identity-hot", no_argument, 0, 'i'},
            {"aln-label", required_argument, 0, 'l'},
            {"reads", required_argument, 0, 'r'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "Aaih?wM:fmpx:g:l:H:r:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'l':
                aln_label = optarg;
                break;
            case 'r':
                emit_warning(context, "The --reads option is deprecated.");
                break;
            case '?':
            case 'h':
                help_vectorize(argv);
                return 1;
            case 'x':
                xg_name = require_exists(context, optarg);
                break;
            case 'g':
                gcsa_name = require_exists(context, optarg);
                // We also need the LCP index
                require_exists(context, gcsa_name + ".lcp");
                break;
            case 'm':
                mem_sketch = true;
                break;
            case 'p':
                mem_positions = true;
                break;
            case 'H':
                mem_hit_max = parse<int>(optarg);
                break;
            case 'a':
                a_hot = true;
                break;
            case 'w':
                output_wabbit = true;
                break;
            case 'i':
                use_identity_hot = true;
                break;
            case 'f':
                format = true;
                break;
            case 'A':
                annotate = true;
                format = true;
                break;
            case 'M':
                wabbit_mapping_file = ensure_writable(context, optarg);
                break;
            default:
                abort();
        }
    }

    PathPositionHandleGraph* xg_index = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionOverlayHelper overlay_helper;
    if (!xg_name.empty()) {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        xg_index = overlay_helper.apply(path_handle_graph.get());
    }
    else{
        error_and_exit(context, "No XG index given. An XG index must be provided.");
    }

    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);

    unique_ptr<gcsa::GCSA> gcsa_index;
    unique_ptr<gcsa::LCPArray> lcp_index;
    if (!gcsa_name.empty()) {
        gcsa_index = vg::io::VPKG::load_one<gcsa::GCSA>(gcsa_name);
        
        // default LCP is the gcsa base name +.lcp
        string lcp_name = gcsa_name + ".lcp";
        lcp_index = vg::io::VPKG::load_one<gcsa::LCPArray>(lcp_name);
    }

    Mapper* mapper = nullptr;
    if (mem_sketch) {
        if (gcsa_name.empty()) {
            error_and_exit(context, "an XG index and GCSA index are required when making MEM sketches");
        } else {
            mapper = new Mapper(xg_index, gcsa_index.get(), lcp_index.get());
        }
        if (mem_hit_max) {
            mapper->hit_max = mem_hit_max;
        }
    }

    Vectorizer vz(xg_index);

    // write the header if needed
    if (format) {
        cout << "aln.name";
        xg_index->for_each_handle([&](handle_t handle) {
                cout << "\tnode." << xg_index->get_id(handle);
            });
        cout << endl;
    }

    //Generate a 1-hot coverage vector for graph entities.
    function<void(Alignment&)> lambda = [&vz, &mapper, use_identity_hot, output_wabbit, aln_label, mem_sketch, 
                                         mem_positions, format, a_hot, max_mem_length](Alignment& a){
        //vz.add_bv(vz.alignment_to_onehot(a));
        //vz.add_name(a.name());
        if (a_hot) {
            vector<int> v = vz.alignment_to_a_hot(a);
            if (output_wabbit){
                cout << vz.wabbitize(aln_label == "" ? a.name() : aln_label, v) << endl;
            }
            else if (format){
                cout << a.name() << "\t" << vz.format(v) << endl;
            } else{
                cout << v << endl;
            }
        }
        else if (use_identity_hot){
            vector<double> v = vz.alignment_to_identity_hot(a);
            if (output_wabbit){
                cout << vz.wabbitize(aln_label == "" ? a.name() : aln_label, v) << endl;
            }
            else if (format){
                cout << a.name() << "\t" << vz.format(v) << endl;
            }
            else {
                cout << vz.format(v) << endl;
            }

        } else if (mem_sketch) {
            // get the mems
            map<string, int> mem_to_count;
            auto mems = mapper->find_mems_simple(a.sequence().begin(), a.sequence().end(),
                                                 max_mem_length, mapper->min_mem_length);
            for (auto& mem : mems) {
                mem_to_count[mem.sequence()]++;
            }
            cout << " |info count:" << mems.size() << " unique:" << mem_to_count.size();
            cout << " |mems";
            for (auto m : mem_to_count) {
                cout << " " << m.first << ":" << m.second;
            }
            if (mem_positions) {
                cout << " |positions";
                for (auto& mem : mems) {
                    for (auto& node : mem.nodes) {
                        cout << " " << gcsa::Node::id(node);
                        if (gcsa::Node::rc(node)) {
                            cout << "-";
                        } else {
                            cout << "+";
                        }
                        cout << ":" << mem.end - mem.begin;
                    }
                }
            }
            cout << endl;
        } else {
            bit_vector v = vz.alignment_to_onehot(a);
            if (output_wabbit){
                cout << vz.wabbitize(aln_label == "" ? a.name() : aln_label, v) << endl;
            } else if (format) {
                cout << a.name() << "\t" << vz.format(v) << endl;
            } else{
                cout << v << endl;
            }
        }
    };
    
    get_input_file(optind, argc, argv, [&](istream& in) {
        vg::io::for_each(in, lambda);
    });

    string mapping_str = vz.output_wabbit_map();
    if (output_wabbit){
        if (!wabbit_mapping_file.empty()){
            ofstream ofi;
            ofi.open(wabbit_mapping_file);
            ofi << mapping_str;
            ofi.close();
        }
        else{

            cerr << mapping_str;
        }
    }


    delete mapper;

    return 0;
}

// Register subcommand
static Subcommand vg_vectorize("vectorize", "transform alignments to simple ML-compatible vectors", main_vectorize);

