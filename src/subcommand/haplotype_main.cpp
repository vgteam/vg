#include <getopt.h>

#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../haplotypes.hpp"
#include "../haplotype_extracter.hpp"

using namespace vg;
using namespace std;
using namespace vg::subcommand;

void help_haplo(char** argv) {
    cerr << "usage: " << argv[0] << " <command> [options]" << endl
         << "options:" << endl
         << "    -x, --xg-name FILE         input xg index" << endl
         << "    -n, --display-path-names   lists path in xg by name" << endl
         << "    -l, --path_labels PATH     lists node labels for path specified by PATH" << endl
         << "    -o, --output-path PATH     folder to which to output statistics" << endl
         << "    -q, --query-haplo PATH     outputs statistics for PATH" << endl
         << "    -A, --all-haplos           outputs statistics for all paths in xg index" << endl
         << "    -s, --start-node node      node (in xg index) at which to start assessment of whole xg index" << endl
         << "    -e, --end-node node        node (in xg index) at which to end assessment of whole xg index, inclusive" << endl
         << "    -i, --inner-index integer        where (within each node) to start enumeration" << endl
         << "    -g, --graphical            output a (non-compact) array showing all strips" << endl;
}

int main_haplo(int argc, char** argv) {
  if (argc == 2) {
      help_haplo(argv);
      return 1;
  }

  string xg_name;
  string named_path;

  bool use_embedded_thread = false;
  int64_t start_node;
  int extend_distance = -1;
  bool backwards = false;
  int64_t offset = 0;

  double recombination_penalty = 9;

  bool display_path_names = false;
  bool decomposition_only = false;
  bool print_haplod = false;
  bool print_haplod_detailed = false;

  bool tests = false;
  bool alt_calc_I = false;

  bool recombine = false;
  int64_t target = -2;
  int cutpoint = -2;
  int joinpoint = -2;
  bool findnode = false;

  int c;
  optind = 2;
  while (true) {
    static struct option long_options[] =
    {
      /* These options set a flag. */
      //{"verbose", no_argument,       &verbose_flag, 1},
      {"xg-name", required_argument, 0, 'x'},
      {"display-path-names", no_argument, 0, 'N'},
      {"named_path", required_argument, 0, 'P'},
      {"start-node", required_argument, 0, 'n'},
      {"offset", required_argument, 0, 'i'},
      {"backwards", no_argument, 0, 'b'},
      {"extend-distance", required_argument, 0, 'd'},
      {"recombination-penalty", required_argument, 0, 'r'},
      {"decomposition-only", no_argument, 0, 'o'},
      {"print-haplo-decomposition", no_argument, 0, 'y'},
      {"print-detailed", no_argument, 0, 'z'},
      {"alt_calc_I", no_argument, 0, 's'},
      {"recombine", no_argument, 0, 'R'},
      {"target", required_argument, 0, 'T'},
      {"cutpoint", required_argument, 0, 'C'},
      {"joinpoint", required_argument, 0, 'J'},
      {"find-node", required_argument, 0, 'F'},
      {0, 0, 0, 0}
    };

    int option_index = 0;
    c = getopt_long (argc, argv, "x:NPn:bi:d:r:oyzsRT:C:J:F:h",
    long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1)
    break;

    switch (c)
    {
      case 'x':
      xg_name = optarg;
      break;

      case 'N':
      display_path_names = true;
      break;

      case 'P':
      named_path = optarg;
      break;

      case 'n':
      use_embedded_thread = true;
      start_node = atoi(optarg);
      break;

      case 'b':
      backwards = true;
      break;

      case 'i':
      offset = atoi(optarg);
      break;

      case 'd':
      extend_distance = atoi(optarg);
      break;

      case 'r':
      recombination_penalty = atof(optarg);
      break;

      case 'o':
      decomposition_only = true;
      break;

      case 'y':
      print_haplod = true;
      break;

      case 'z':
      print_haplod_detailed = true;
      break;

      case 's':
      alt_calc_I = true;
      break;

      case 'R':
      recombine = true;
      break;

      case 'T':
      target = atoi(optarg);
      break;

      case 'C':
      cutpoint = atoi(optarg);
      break;

      case 'J':
      joinpoint = atoi(optarg);
      break;

      case 'F':
      findnode = true;
      break;

      case '?':
      case 'h':
          help_haplo(argv);
          return 1;

      default:
          help_haplo(argv);
          abort();
    }
  }

  xg::XG xindex;
  if (!xg_name.empty()) {
      ifstream in(xg_name.c_str());
      xindex.load(in);
  } else {
    cerr << "[xg haplo] error: need to specify an xg index to work on" << endl;
  }

  if(tests) {
    cerr << "running log-probability calculation tests" << endl;
    logRR_tests(recombination_penalty);
  }

  if(display_path_names){
    cerr << "The xg index contains the following named paths, by rank:" << endl;
    for(size_t path_rank = 1; path_rank <= xindex.max_path_rank(); path_rank++) {
      cout << path_rank << "\t" << xindex.path_name(path_rank) << endl;
    }
  }

  if(!named_path.empty() && use_embedded_thread) {
    cerr << "[xg haplo] error: can't both extract named thread and unnamed embedded thread" << endl;
    return 0;
  }

  if(!recombine) {
    thread_t t;
    if(!named_path.empty()) {
      Path path = xindex.path(named_path);
      t = path_to_thread_t(path);
    } else if(use_embedded_thread) {
      xg::XG::ThreadMapping n;
      n.node_id = start_node;
      n.is_reverse = backwards;
      t = extract_thread(xindex, n, offset, extend_distance);
    }
    if(t.size() == 0) {
      cerr << "[xg haplo] error: thread of length zero queried" << endl;
    } else {
      haplo_d h = haplo_d(t, xindex);
      if(!alt_calc_I) {
        h.log_calculate_Is(xindex);
      } else {
        h.seeded_log_calculate_Is(xindex);
      }
      if(!decomposition_only) {
        RRMemo memo(recombination_penalty);
        cout << h.log_probability(memo) << endl;
      }
      if(print_haplod) {
        h.print(cout);
      }
      if(print_haplod_detailed) {
        h.print_detailed(cout);
      }
    }
  } else {
    thread_t l;
    thread_t r;
    if(!use_embedded_thread) {
      cerr << "[xg haplo] error: need to specify LHS path (use -n <node> -i <offset> -d <length> (-b))" << endl;
    } else {
      xg::XG::ThreadMapping n;
      n.node_id = start_node;
      n.is_reverse = backwards;
      l = extract_thread(xindex, n, offset, extend_distance);
    }
    if(target == -2) {
      cerr << "[xg haplo] error: need to specify RHS haplotype offset (use -T <offset>)" << endl;
    } else {
      xg::XG::ThreadMapping n;
      n.node_id = start_node;
      n.is_reverse = backwards;
      r = extract_thread(xindex, n, target, extend_distance);
    }
    // if(l.size() == 0 || r.size() == 0) {
    //   cerr << "[xg haplo] error: need two nonempty haplotypes to recombine" << endl;
    // } else {
    //   if(findnode) {
    //     if(cutpoint == -2) {
    //       cerr << "[xg haplo] error: need specify cutpoint (use -C <cutpoint>)" << endl;
    //     } else {
    //       haplo_d hl = haplo_d(l, xindex);
    //       haplo_d hr = haplo_d(r, xindex);
    //       xg::XG::ThreadMapping to_find = hl.cs[cutpoint].get_node();
    //       int result = find_node(hr, to_find, cutpoint);
    //       cout << result << endl;
    //       cout << start_node << "\t" << r[result].node_id << endl;
    //     }
    //   } else {
    //     if(cutpoint == -2 || joinpoint == -2) {
    //       cerr << "[xg haplo] error: need specify cutpoint and joinpoint (use -C <cutpoint> -J <joinpoint>)" << endl;
    //     } else {
    //       haplo_d hl = haplo_d(l, xindex);
    //       hl.log_calculate_Is(xindex);
    //       haplo_d hr = haplo_d(r, xindex);
    //       hr.log_calculate_Is(xindex);
    //       if(joinpoint > 0) {
    //         haplo_d h = recombine_arms(hl, hr, cutpoint, joinpoint, xindex);
    //         h.print_detailed(cout);
    //       }
    //     }
    //   }
    // }
  }

  return 0;
}

static Subcommand vg_haplo("haplo", "work with haplotype decompositions", main_haplo);
