#include <getopt.h>

#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../haplotype_extracter.hpp"

using namespace vg;
using namespace std;
using namespace vg::subcommand;

using thread_t = vector<xg::XG::ThreadMapping>;

void help_trace(char** argv) {
    cerr << "usage: " << argv[0] << " trace [options]" << endl
         << "options:" << endl
         << "    -x, --index FILE           use this xg index" << endl
         << "    -n, --start-node INT       start at this node" << endl
        //TODO: implement backwards iteration over graph
        // << "    -b, --backwards            iterate backwards over graph" << endl
         << "    -d, --extend-distance INT  extend search this many nodes" << endl
         << "    -j, --output-json          path to which to output json" << endl;
}

int main_trace(int argc, char** argv) {
  if (argc == 2) {
      help_trace(argv);
      return 1;
  }

  string xg_name;
  string output_path;
  int64_t start_node;
  int extend_distance = 50;
  bool backwards = false;

  int c;
  optind = 2; // force optind past command positional argument
  while (true) {
    static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"index", required_argument, 0, 'x'},
            {"output-json", required_argument, 0, 'j'},
            {"start-node", required_argument, 0, 'n'},
            {"extend-distance", required_argument, 0, 'd'},
            //{"backwards", no_argument, 0, 'b'},
            {0, 0, 0, 0}
        };

    int option_index = 0;
    c = getopt_long (argc, argv, "x:j:n:d:h",
                     long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
        break;

    switch (c)
    {
    case 'x':
        xg_name = optarg;
        break;

    case 'j':
        output_path = optarg;
        break;

    case 'n':
        start_node = atoi(optarg);
        break;

    case 'd':
        extend_distance = atoi(optarg);
        break;

    //case 'b':
        //backwards = true;
        //break;

    case '?':
    case 'h':
        help_trace(argv);
        return 1;

    default:
        help_trace(argv);
        abort();
    }
  }
  
  xg::XG xindex;
  if (!xg_name.empty()) {
      ifstream in(xg_name.c_str());
      xindex.load(in);
  }
  
  xg::XG::ThreadMapping n;
  n.node_id = start_node;
  n.is_reverse = backwards;
  int64_t offset = 0;
  
  if(!output_path.empty()) {
    ofstream annotation_ofstream(output_path+".annotation");
    ofstream json_ofstream(output_path);
    // what subhaplotypes of length extend_distance starting at n are embedded
    // in xindex? How many identical copies of each subhaplotype?
    vector<pair<thread_t,int> > haplotype_list =
              list_haplotypes(xindex, n, extend_distance);
    // haplotype counts go into annotation_ofstream
    output_haplotype_counts(annotation_ofstream, haplotype_list, xindex);
    // json containing graph and paths goes into json_ofstream
    output_graph_with_embedded_paths(json_ofstream, haplotype_list, xindex);
  }

  return 0;
}

static Subcommand vg_trace("trace", "trace haplotypes", main_trace);
