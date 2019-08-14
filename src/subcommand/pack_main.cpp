#include "subcommand.hpp"
#include "../vg.hpp"
#include "../xg.hpp"
#include "../utility.hpp"
#include "../packer.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_pack(char** argv) {
    cerr << "usage: " << argv[0] << " pack [options]" << endl
         << "options:" << endl
         << "    -x, --xg FILE          use this basis graph" << endl
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
         << "    -Q, --min-mapq N       ignore read mappings with Mapping Quality < N [default: 0]" << endl
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
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:o:i:g:dDt:eb:n:N:qQ:",
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
            break;
        default:
            abort();
        }
    }

    omp_set_num_threads(thread_count);

    unique_ptr<PathPositionHandleGraph> xgidx;
    if (xg_name.empty()) {
        cerr << "No XG index given. An XG index must be provided." << endl;
        exit(1);
    } else {
        xgidx = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_name);
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

    // todo one packer per thread and merge

    vg::Packer packer(xgidx.get(), bin_size, qual_adjust);
    if (packs_in.size() == 1) {
        packer.load_from_file(packs_in.front());
    } else if (packs_in.size() > 1) {
        packer.merge_from_files(packs_in);
    }

    if (!gam_in.empty()) {
        vector<vg::Packer*> packers;
        if (thread_count == 1) {
            packers.push_back(&packer);
        } else {
            for (size_t i = 0; i < thread_count; ++i) {
                packers.push_back(new Packer(xgidx.get(), bin_size, qual_adjust));
            }
        }
        std::function<void(Alignment&)> lambda = [&packer,&record_edits,&packers,&min_mapq](Alignment& aln) {
            if (aln.mapping_quality() >= min_mapq) {
                packers[omp_get_thread_num()]->add(aln, record_edits);
            }
        };
        if (gam_in == "-") {
            vg::io::for_each_parallel(std::cin, lambda);
        } else {
            ifstream gam_stream(gam_in);
            vg::io::for_each_parallel(gam_stream, lambda);
            gam_stream.close();
        }
        if (thread_count == 1) {
            packers.clear();
        } else {
            packer.merge_from_dynamic(packers);
            for (auto& p : packers) {
                delete p;
            }
            packers.clear();
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

    return 0;
}

// Register subcommand
static Subcommand vg_pack("pack", "convert alignments to a compact coverage, edit, and path index", main_pack);
