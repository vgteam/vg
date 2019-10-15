#include "subcommand.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../utility.hpp"
#include "../packer.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <handlegraph/handle_graph.hpp>
#include <bdsg/overlay_helper.hpp>

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
         << "    -g, --gam FILE         read alignments from this file (could be '-' for stdin)" << endl
         << "    -d, --as-table         write table on stdout representing packs" << endl
         << "    -D, --as-edge-table    write table on stdout representing edge coverage" << endl
         << "    -e, --with-edits       record and write edits rather than only recording graph-matching coverage" << endl
         << "    -b, --bin-size N       number of sequence bases per CSA bin [default: inf]" << endl
         << "    -n, --node ID          write table for only specified node(s)" << endl
         << "    -N, --node-list FILE   a white space or line delimited list of nodes to collect" << endl
         << "    -q, --qual-adjust      scale coverage by phred quality (combined from mapq and base quality)" << endl
         << "    -Q, --min-mapq N       ignore reads with MAPQ < N and positions with base quality < N [default: 0]" << endl
         << "    -c, --expected-cov N   expected coverage.  used only for memory tuning [default : 128]" << endl
         << "    -t, --threads N        use N threads (defaults to numCPUs)" << endl;
}


int main_pack(int argc, char** argv) {

    string xg_name;
    vector<string> packs_in;
    string packs_out;
    string gam_in;
    bool write_table = false;
    bool write_edge_table = false;
    int thread_count = 1;
    bool record_edits = false;
    size_t bin_size = 0;
    vector<vg::id_t> node_ids;
    string node_list_file;
    bool qual_adjust = false;
    int min_mapq = 0;
    int min_baseq = 0;
    size_t expected_coverage = 128;

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
            {"as-table", no_argument, 0, 'd'},
            {"as-edge-table", no_argument, 0, 'D'},
            {"threads", required_argument, 0, 't'},
            {"with-edits", no_argument, 0, 'e'},
            {"node", required_argument, 0, 'n'},
            {"node-list", required_argument, 0, 'N'},
            {"bin-size", required_argument, 0, 'b'},
            {"qual-adjust", no_argument, 0, 'q'},
            {"min-mapq", required_argument, 0, 'Q'},
            {"expected-cov", required_argument, 0, 'c'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:o:i:g:dDt:eb:n:N:qQ:c:",
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
        case 'd':
            write_table = true;
            break;
        case 'D':
            write_edge_table = true;
            break;
        case 'e':
            record_edits = true;
            break;
        case 'b':
            bin_size = atoll(optarg);
            break;
        case 't':
            thread_count = parse<int>(optarg);
            break;
        case 'n':
            node_ids.push_back(parse<int>(optarg));
            break;
        case 'N':
            node_list_file = optarg;
            break;
        case 'q':
            qual_adjust = true;
            break;
        case 'Q':
            min_mapq = parse<int>(optarg);
            min_baseq = min_mapq;
            break;
        case 'c':
            expected_coverage = parse<size_t>(optarg);
            break;
        default:
            abort();
        }
    }

    omp_set_num_threads(thread_count);

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

    if (gam_in.empty() && packs_in.empty()) {
        cerr << "error [vg pack]: Input must be provided with -g or -i" << endl;
        exit(1);
    }

    if (packs_out.empty() && write_table == false && write_edge_table == false) {
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
    size_t data_width = std::ceil(std::log2(2 * expected_coverage));

    // create our packers
    size_t num_packers = !gam_in.empty() ? thread_count : 1;
    vector<Packer*> packers(num_packers, nullptr);
#pragma omp parallel for
    for (int i = 0; i < packers.size(); ++i) {
        packers[i] = new Packer(graph, bin_size, thread_count, data_width, true, true, record_edits);
    }

    Packer& packer = *packers[0];
    
    // todo one packer per thread and merge
    if (packs_in.size() == 1) {
        packer.load_from_file(packs_in.front());
    } else if (packs_in.size() > 1) {
        packer.merge_from_files(packs_in);
    }

    if (!gam_in.empty()) {
        std::function<void(Alignment&)> lambda = [&packer,&min_mapq,&min_baseq,&qual_adjust,&packers](Alignment& aln) {
            packers[omp_get_thread_num()]->add(aln, min_mapq, min_baseq, qual_adjust);
        };
        if (gam_in == "-") {
            vg::io::for_each_parallel(std::cin, lambda);
        } else {
            ifstream gam_stream(gam_in);
            if (!gam_stream) {
                cerr << "[vg pack] error reading gam file: " << gam_in << endl;
                return 1;
            }
            vg::io::for_each_parallel(gam_stream, lambda);
            gam_stream.close();
        }
        if (packers.size() > 1) {
            vector<Packer*> others(packers.begin() + 1, packers.end());
            packer.merge_from_dynamic(others);
            for (auto other : others) {
                delete other;
            }
            packers.resize(1);
        }
    }

    if (!packs_out.empty()) {
        packer.save_to_file(packs_out);
    }
    if (write_table || write_edge_table) {
        packer.make_compact();
        if (write_table) {
            packer.as_table(cout, record_edits, node_ids);
        }
        if (write_edge_table) {
            packer.as_edge_table(cout, node_ids);
        }
    }

    for (int i = 0; i < packers.size(); ++i) {
        delete packers[i];
    }        
    packers.clear();
    
    return 0;
}

// Register subcommand
static Subcommand vg_pack("pack", "convert alignments to a compact coverage index", PIPELINE, 6, main_pack);
