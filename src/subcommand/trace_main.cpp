#include <getopt.h>

#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../haplotypes.hpp"
#include "../haplotype_extracter.hpp"

using namespace vg;
using namespace std;


void help_trace(char** argv) {
    cerr << "usage: " << argv[0] << " trace [options]" << endl
         << "options:" << endl
         << "    -x, --index FILE           use this xg index" << endl
         << "    -n, --start-node INT       start at this node" << endl
         << "    -b, --backwards            iterate backwards over graph" << endl
         << "    -d, --extend-distance INT  extend search this many nodes" << endl
         << "    -j, --output-json          path to which to output json" << endl
         << "    -l, --likelihoods          print haplotype likelihoods to annotation file" << endl;
}

int main_trace(int argc, char** argv) {
  if (argc == 2) {
      help_trace(argv);
      return 1;
  }

  string xg_name;
  string json_out;
  int64_t start_node;
  int extend_distance = 50;
  bool backwards = false;
  bool likelihoods = false;
  bool test_decomp = false;
  bool test_decomp_log = false;

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
            {"backwards", no_argument, 0, 'b'},
            {"likelihoods", no_argument, 0, 'l'},
            {"test-haplo-decomps", no_argument, 0, 'T'},
            {"test-haplo-decomps-log", no_argument, 0, 'L'},
            {0, 0, 0, 0}
        };

    int option_index = 0;
    c = getopt_long (argc, argv, "x:j:n:d:blLTh",
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
        json_out = optarg;
        break;

    case 'n':
        start_node = atoi(optarg);
        break;

    case 'd':
        extend_distance = atoi(optarg);
        break;

    case 'b':
        backwards = true;
        break;

    case 'l':
        likelihoods = true;
        break;

    case 'T':
        test_decomp = true;
        break;

    case 'L':
        test_decomp_log = true;
        break;

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
  if(!json_out.empty() && !test_decomp && !test_decomp_log) {
    vector<pair<thread_t,int> > haplotype_list = list_haplotypes(xindex, n, extend_distance);
    output_weighted_haplotype_list(json_out, haplotype_list, xindex, likelihoods);
    output_graph_with_embedded_paths(json_out, haplotype_list, xindex);
  }

  if(test_decomp) {
    thread_t t = extract_thread(xindex, n, offset, extend_distance);
    haplo_d h = haplo_d(t, xindex);
    h.calculate_Is(xindex);
    h.print_detailed_searchstates(cout);
  }

  if(test_decomp_log) {
    thread_t t = extract_thread(xindex, n, offset, extend_distance);
    cerr << "made thread" << endl;
    haplo_d h = haplo_d(t, xindex);
    h.log_calculate_Is(xindex);
    cerr << "Generated I's by binary search method" << endl;
    h.print_detailed_searchstates(cout);
  }

  return 0;
}
