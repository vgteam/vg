#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../packer.hpp"
#include "../stream.hpp"

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
         << "    -n, --no-edits         don't record or write edits, just graph-matching coverage" << endl
         << "    -b, --bin-size N       number of sequence bases per CSA bin" << endl
         << "    -t, --threads N        use N threads (defaults to numCPUs)" << endl;
}

int main_pack(int argc, char** argv) {

    string xg_name;
    vector<string> packs_in;
    string packs_out;
    string gam_in;
    bool write_table = false;
    int thread_count = 1;
    bool record_edits = true;
    size_t bin_size = 1 << 20;

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
            {"threads", required_argument, 0, 't'},
            {"no-edits", no_argument, 0, 'n'},
            {"bin-size", required_argument, 0, 'b'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:o:i:g:dt:nb:",
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
        case 'n':
            record_edits = false;
            break;
        case 'b':
            bin_size = atoll(optarg);
            break;
        case 't':
            thread_count = atoi(optarg);
            break;

        default:
            abort();
        }
    }

    omp_set_num_threads(thread_count);

    xg::XG xgidx;
    if (xg_name.empty()) {
        cerr << "No XG index given. An XG index must be provided." << endl;
        exit(1);
    } else {
        ifstream in(xg_name.c_str());
        xgidx.load(in);
    }

    // todo one packer per thread and merge
    vg::Packer packer(&xgidx, bin_size);
    if (packs_in.size() == 1) {
        packer.load_from_file(packs_in.front());
    } else if (packs_in.size() > 1) {
        packer.merge_from_files(packs_in);
    }

    if (!gam_in.empty()) {
        std::function<void(Alignment&)> lambda = [&packer,&record_edits](Alignment& aln) { packer.add(aln, record_edits); };
        if (gam_in == "-") {
            stream::for_each(std::cin, lambda);
        } else {
            ifstream gam_stream(gam_in);
            stream::for_each(gam_stream, lambda);
            gam_stream.close();
        }
    }
    if (!packs_out.empty()) {
        packer.save_to_file(packs_out);
    }
    if (write_table) {
        packer.make_compact();
        packer.as_table(cout, record_edits);
    }

    return 0;
}

// Register subcommand
static Subcommand vg_pack("pack", "convert alignments to a compact coverage, edit, and path index", main_pack);
