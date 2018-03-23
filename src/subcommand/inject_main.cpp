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

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_inject(char** argv) {
    cerr << "usage: " << argv[0] << " inject [options] input.bam >output.gam" << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE       use the graph in this xg index" << endl
         //<< "    -c, --cram-name FILE     input CRAM file" << endl
         //<< "    -b, --bam-name FILE      input BAM file" << endl
         //<< "    -s, --sam-name FILE      input SAM file" << endl
         << "    -t, --threads N          number of threads to use" << endl
         << "    -p, --into-path NAME     inject from this path (many allowed, default: all in bam)" << endl
         << "    -F, --into-paths FILE    inject from nonoverlapping path names listed in FILE (one per line)" << endl;
}

int main_inject(int argc, char** argv) {

    if (argc == 2) {
        help_inject(argv);
    }

    string xg_name;
    set<string> path_names;
    string path_prefix;
    string path_file;
    string output_type = "gam";
    string input_type = "bam";
    int threads = 1;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
          {"help", no_argument, 0, 'h'},
          {"xg-name", required_argument, 0, 'x'},
          {"threads", required_argument, 0, 't'},
          {"into-path", required_argument, 0, 'p'},
          {"into-paths", required_argument, 0, 'F'},
          {"into-prefix", required_argument, 0, 'P'},
          {"cram-name", no_argument, 0, 'c'},
          {"bam-name", no_argument, 0, 'b'},
          {"sam-name", no_argument, 0, 's'},
          {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:p:F:P:cbst:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;

        case 'p':
            path_names.insert(optarg);
            break;

        case 'F':
            path_file = optarg;
            break;

        case 'c':
          output_type = "cram";
          break;

        case 'b':
          output_type = "bam";
          break;

        case 's':
          output_type = "sam";
          break;

        case 't':
          threads = atoi(optarg);
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

    if (!path_file.empty()){
      // open the file
      ifstream in(path_file);
      string line;
      while (std::getline(in,line)) {
        path_names.insert(line);
      }
    }

    xg::XG* xgidx = nullptr;
    ifstream xg_stream(xg_name);
    if(xg_stream) {
      xgidx = new xg::XG(xg_stream);
    }
    if (!xg_stream || xgidx == nullptr) {
      cerr << "[vg surject] error: could not open xg index" << endl;
      return 1;
    }

    // if no paths were given take all of those in the index
    if (path_names.empty()) {
      for (size_t i = 1; i <= xgidx->path_count; ++i) {
        path_names.insert(xgidx->path_name(i));
      }
    }

    //function<void(const Alignment&)>& lambda) {
    // todo write buffering procedure in alignment.cpp
    vector<Alignment> buf;
    function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
      buf.push_back(aln);
      if (buf.size() > 1000) {
        write_alignments(std::cout, buf);
        buf.clear();
      }
    };
    if (threads > 1) {
        hts_for_each_parallel(file_name, lambda);
    } else {
        hts_for_each(file_name, lambda);
    }
    write_alignments(std::cout, buf);
    buf.clear();
    cout.flush();
    return 0;
}

// Register subcommand
static Subcommand vg_inject("inject", "lift over alignments for the graph", main_inject);
