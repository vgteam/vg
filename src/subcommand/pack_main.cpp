#include "subcommand.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../utility.hpp"
#include "../packer.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <handlegraph/handle_graph.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_pack(char** argv) {
    cerr << "usage: " << argv[0] << " pack [options]" << endl
         << "options:" << endl
         << "    -x, --xg FILE          use this basis graph (any format accepted, does not have to be xg)" << endl
         << "    -o, --packs-out FILE   write compressed coverage packs to this output file" << endl
         << "    -i, --packs-in FILE    begin by summing coverage packs from each provided FILE" << endl
         << "    -g, --gam FILE         read alignments from this GAM file (could be '-' for stdin)" << endl
         << "    -a, --gaf FILE         read alignments from this GAF file (could be '-' for stdin)" << endl
         << "    -d, --as-table         write table on stdout representing packs" << endl
         << "    -D, --as-edge-table    write table on stdout representing edge coverage" << endl
         << "    -u, --as-qual-table    write table on stdout representing average node mapqs" << endl
         << "    -e, --with-edits       record and write edits rather than only recording graph-matching coverage" << endl
         << "    -b, --bin-size N       number of sequence bases per CSA bin [default: inf]" << endl
         << "    -n, --node ID          write table for only specified node(s)" << endl
         << "    -N, --node-list FILE   a white space or line delimited list of nodes to collect" << endl
         << "    -Q, --min-mapq N       ignore reads with MAPQ < N and positions with base quality < N [default: 0]" << endl
         << "    -c, --expected-cov N   expected coverage.  used only for memory tuning [default : 128]" << endl
         << "    -s, --trim-ends N      ignore the first and last N aligned bases of each read" << endl 
         << "    -t, --threads N        use N threads (defaults to numCPUs)" << endl;
}


int main_pack(int argc, char** argv) {

    string xg_name;
    vector<string> packs_in;
    string packs_out;
    string gam_in;
    string gaf_in;
    bool write_table = false;
    bool write_edge_table = false;
    bool write_qual_table = false;
    bool record_edits = false;
    size_t bin_size = 0;
    vector<vg::id_t> node_ids;
    string node_list_file;
    int min_mapq = 0;
    int min_baseq = 0;
    size_t expected_coverage = 128;
    int trim_ends = 0;

    if (argc == 2) {
        help_pack(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg", required_argument,0, 'x'},
            {"packs-out", required_argument,0, 'o'},
            {"count-in", required_argument, 0, 'i'},
            {"gam", required_argument, 0, 'g'},
            {"gaf", required_argument, 0, 'a'},
            {"as-table", no_argument, 0, 'd'},
            {"as-edge-table", no_argument, 0, 'D'},
            {"as-qual-table", no_argument, 0, 'u'},
            {"threads", required_argument, 0, 't'},
            {"with-edits", no_argument, 0, 'e'},
            {"node", required_argument, 0, 'n'},
            {"node-list", required_argument, 0, 'N'},
            {"bin-size", required_argument, 0, 'b'},
            {"min-mapq", required_argument, 0, 'Q'},
            {"expected-cov", required_argument, 0, 'c'},
            {"trim-ends", required_argument, 0, 's'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:o:i:g:a:dDut:eb:n:N:Q:c:s:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_pack(argv);
            return 1;
        case 'x':
            xg_name = optarg;
            break;
        case 'o':
            packs_out = optarg;
            break;
        case 'i':
            packs_in.push_back(optarg);
            break;
        case 'g':
            gam_in = optarg;
            break;
        case 'a':
            gaf_in = optarg;
            break;
        case 'd':
            write_table = true;
            break;
        case 'D':
            write_edge_table = true;
            break;
        case 'u':
            write_qual_table = true;
            break;
        case 'e':
            record_edits = true;
            break;
        case 'b':
            bin_size = atoll(optarg);
            break;            
        case 't':
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg pack] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }
        case 'n':
            node_ids.push_back(parse<int>(optarg));
            break;
        case 'N':
            node_list_file = optarg;
            break;
        case 'Q':
            min_mapq = parse<int>(optarg);
            min_baseq = min_mapq;
            break;
        case 'c':
            expected_coverage = parse<size_t>(optarg);
            break;
        case 's':
            trim_ends = parse<int>(optarg);
            break;
        default:
            abort();
        }
    }

    unique_ptr<HandleGraph> handle_graph;
    HandleGraph* graph = nullptr;
    if (xg_name.empty()) {
        cerr << "error [vg pack]: No basis graph given. One must be provided with -x." << endl;
        exit(1);
    } else {
        handle_graph = vg::io::VPKG::load_one<HandleGraph>(xg_name);
    }
    bdsg::VectorizableOverlayHelper overlay_helper;
    graph = dynamic_cast<HandleGraph*>(overlay_helper.apply(handle_graph.get()));

    if (gam_in.empty() && packs_in.empty() && gaf_in.empty()) {
        cerr << "error [vg pack]: Input must be provided with -g, -a or -i" << endl;
        exit(1);
    }

    if (!gam_in.empty() && !gaf_in.empty()) {
        cerr << "error [vg pack]: -g cannot be used with -a" << endl;
        exit(1);
    }

    if (packs_out.empty() && write_table == false && write_edge_table == false && write_qual_table == false) {
        cerr << "error [vg pack]: Output must be selected with -o, -d or -D" << endl;
        exit(1);
    }

    // process input node list
    if (!node_list_file.empty()) {
        ifstream nli;
        nli.open(node_list_file);
        if (!nli.good()){
            cerr << "[vg pack] error, unable to open the node list input file." << endl;
            exit(1);
        }
        string line;
        while (getline(nli, line)){
            for (auto& idstr : split_delims(line, " \t")) {
                node_ids.push_back(parse<int64_t>(idstr.c_str()));
            }
        }
        nli.close();
    }

    // get a data width from our expected coverage, using simple heuristic of counting
    // bits needed to store double the coverage
    size_t data_width = Packer::estimate_data_width(expected_coverage);

    // use some naive heuristics to come up with bin count and batch size based on thread count
    // more bins: finer grained parallelism at cost of more mutexes and allocations
    // bigger batch size: more robustness to sorted input at cost of less parallelism
    size_t num_threads = get_thread_count();
    size_t batch_size = Packer::estimate_batch_size(num_threads);
    size_t bin_count = Packer::estimate_bin_count(num_threads);

    // create our packer
    Packer packer(graph, true, true, record_edits, true, bin_size, bin_count, data_width);
    
    // todo one packer per thread and merge
    if (packs_in.size() == 1) {
        packer.load_from_file(packs_in.front());
    } else if (packs_in.size() > 1) {
        packer.merge_from_files(packs_in);
    }

    std::function<void(Alignment&)> lambda = [&packer,&min_mapq,&min_baseq,&trim_ends](Alignment& aln) {
        packer.add(aln, min_mapq, min_baseq, trim_ends);
    };

    if (!gam_in.empty()) {
        get_input_file(gam_in, [&](istream& in) {
                vg::io::for_each_parallel(in, lambda, batch_size);
            });
    } else if (!gaf_in.empty()) {
        // we use this interface so we can ignore sequence, which takes a lot of time to parse
        // and is unused by pack
        function<size_t(nid_t)> node_to_length = [&graph](nid_t node_id) {
            return graph->get_length(graph->get_handle(node_id));
        };
        function<string(nid_t, bool)> node_to_sequence = [&graph](nid_t node_id, bool is_reversed) {
            return graph->get_sequence(graph->get_handle(node_id, is_reversed));
        };

        // computed batch size was tuned for GAM performance.  some small tests show that
        // gaf benefits from a slightly larger one. 
        vg::io::gaf_unpaired_for_each_parallel(node_to_length, record_edits ? node_to_sequence : nullptr,
                                               gaf_in, lambda, batch_size * 4);
    }

    if (!packs_out.empty()) {
        packer.save_to_file(packs_out);
    }
    if (write_table || write_edge_table || write_qual_table) {
        packer.make_compact();
        if (write_table) {
            packer.as_table(cout, record_edits, node_ids);
        }
        if (write_edge_table) {
            packer.as_edge_table(cout, node_ids);
        }
        if (write_qual_table) {
            packer.as_quality_table(cout, node_ids);
        }
    }

    return 0;
}

// Register subcommand
static Subcommand vg_pack("pack", "convert alignments to a compact coverage index", PIPELINE, 9, main_pack);
