// inject_main.cpp: define the "vg inject" subcommand, which lifts over alignments from the linear space

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include <subcommand.hpp>

#include "../alignment.hpp"
#include "../vg.hpp"
#include "../stream.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_inject(char** argv) {
    cerr << "usage: " << argv[0] << " inject [options] input.[bam|sam|cram] >output.gam" << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE       use the graph in this xg index" << endl
         << "    -t, --threads N          number of threads to use" << endl;
}

int main_inject(int argc, char** argv) {

    if (argc == 2) {
        help_inject(argv);
    }

    string xg_name;
    int threads = 1;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
          {"help", no_argument, 0, 'h'},
          {"xg-name", required_argument, 0, 'x'},
          {"threads", required_argument, 0, 't'},
          {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:t:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;

        case 't':
          threads = parse<int>(optarg);
          omp_set_num_threads(threads);
          break;

        case 'h':
        case '?':
          help_inject(argv);
          exit(1);
          break;

        default:
          abort ();
        }
    }

    string file_name = get_input_file_name(optind, argc, argv);

    xg::XG* xgidx = nullptr;
    ifstream xg_stream(xg_name);
    if(xg_stream) {
      xgidx = new xg::XG(xg_stream);
    }
    if (!xg_stream || xgidx == nullptr) {
      cerr << "[vg inject] error: could not open xg index" << endl;
      return 1;
    }

    stream::ProtobufEmitter<Alignment> buf(cout);
    function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
#pragma omp critical (buf)
        {
            buf.write(std::move(aln));
        }
    };
    if (threads > 1) {
        hts_for_each_parallel(file_name, lambda, xgidx);
    } else {
        hts_for_each(file_name, lambda, xgidx);
    }
    return 0;
}

// Register subcommand
static Subcommand vg_inject("inject", "lift over alignments for the graph", main_inject);
