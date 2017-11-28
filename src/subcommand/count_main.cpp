#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "../counter.hpp"
#include "../stream.hpp"

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;

void help_count(char** argv) {
    cerr << "usage: " << argv[0] << " count [options]" << endl
         << "options:" << endl
         << "    -x, --xg FILE          use this basis graph" << endl
         << "    -o, --counts-out FILE  write compressed coverage counts to this output file" << endl
         << "    -i, --counts-in FILE   begin by summing coverage counts from each provided FILE" << endl
         << "    -g, --gam FILE         read alignments from this file (could be '-' for stdin)" << endl
         << "    -d, --as-table         write table on stdout representing counts" << endl
         << "    -n, --no-edits         don't record or write edits, just graph-matching coverage" << endl
         << "    -t, --threads N        use N threads (defaults to numCPUs)" << endl;
}

int main_count(int argc, char** argv) {

    string xg_name;
    vector<string> counts_in;
    string counts_out;
    string gam_in;
    bool write_table = false;
    int thread_count = 1;
    bool record_edits = true;

    if (argc == 2) {
        help_count(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg", required_argument,0, 'x'},
            {"counts-out", required_argument,0, 'o'},
            {"count-in", required_argument, 0, 'i'},
            {"gam", required_argument, 0, 'g'},
            {"as-table", no_argument, 0, 'd'},
            {"threads", required_argument, 0, 't'},
            {"no-edits", no_argument, 0, 'n'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hx:o:i:g:dt:n",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_count(argv);
            return 1;
        case 'x':
            xg_name = optarg;
            break;
        case 'o':
            counts_out = optarg;
            break;
        case 'i':
            counts_in.push_back(optarg);
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

    // todo one counter per thread and merge
    vg::Counter counter(&xgidx);
    if (counts_in.size() == 1) {
        counter.load_from_file(counts_in.front());
    } else if (counts_in.size() > 1) {
        counter.merge_from_files(counts_in);
    }

    if (!gam_in.empty()) {
        std::function<void(Alignment&)> lambda = [&counter,&record_edits](Alignment& aln) { counter.add(aln, record_edits); };
        if (gam_in == "-") {
            stream::for_each(std::cin, lambda);
        } else {
            ifstream gam_stream(gam_in);
            stream::for_each(gam_stream, lambda);
            gam_stream.close();
        }
    }
    if (!counts_out.empty()) {
        counter.save_to_file(counts_out);
    }
    if (write_table) {
        counter.make_compact();
        counter.as_table(cout, record_edits);
    }

    return 0;
}

// Register subcommand
static Subcommand vg_count("count", "count features on the graph", main_count);
