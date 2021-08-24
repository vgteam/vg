// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce
// more efficient representations of snarls.

#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include "../gbwt_helper.hpp"

#include "subcommand.hpp"

// todo: should be able to remove '../../include/...' and replace with e.g.
// <bdsg/hash...>
#include "../../include/bdsg/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"

#include "../io/save_handle_graph.hpp"

#include "../algorithms/0_oo_normalize_snarls.hpp"
#include "../algorithms/0_snarl_analyzer.hpp"

#include <chrono> // for high_resolution_clock

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_normalize(char **argv) {
  cerr
      << "usage: " << argv[0] << " normalize [options] <graph.vg> >[normalized.vg]"
      << endl
      << "Modifies snarls, outputs modified on stdout." << endl
      << endl
      << "options:" << endl
      << "    -g, --gbwt       gbwt corresponding to hashgraph." << endl
      << "    -s, --snarls       snarls file corresponding to hashgraph."
      << endl
      << "    -o, --output_gbwt   name for the gbwt corresponding to the normalized graph. Default is normalized.gbwt."
      << "    -n, --normalize_type        can be either 'all' or 'one'. "
         "Default is 'all'. If 'all', all top-level snarls in the graph are "
         "normalized. If 'one', only the one specified snarl is normalized."
      << endl
      << "    -a, --source        only used when --normalize_type is 'one'. "
         "The source of the single node to be normalized."
      << endl
      << "    -b, --sink      only used when --normalize_type is 'one'. The "
         "sink of the single node to be normalized."
      << endl
      << "    -p, --paths_right_to_left       only used when --normalize_type "
         "is 'one'. This must be passed to the normalizer when the gbwt "
         "embedded paths are moving through the snarl from right to left, "
         "rather than left to right. This is frequently the case when the "
         "node_ids are larger on the left than on the right in the graph."
      << endl
      << "    -m, --max_alignment_size       limits the size of the snarl that "
         "will be normalized. If the number of gbwt paths through the snarl "
         "exceeds the max size, that snarl will be skipped for normalization. "
         "(Any additional embedded paths stretching from source to sink that "
         "aren't represented in the gbwt will also contribute to the max size.)"
      << endl
      << "    -i, --snarl_sizes       identifies the size of every top-level "
         "snarl of the inputted graph, outputs in a document specified with "
         "format 'source\tsink\tsize\n'. When passed this argument, there is "
         "no normalization."
      << endl
      << "    -k, --snarl_sizes_skip_source_sink      when passed in "
         "conjunction with --snarl_sizes/-i, the snarl size counting ignores "
         "the sequence in the source and sink. This avoids double-counting the "
         "sequence in a source/sink node between two adjacent snarls."
      << endl
      << "    -h, --max_handle_size       currently, default is INT_MAX. This "
         "is for compatibility with changes in graph size measures (file "
         "0_snarl_analyzer, arg '-i'). Eventually should change to handle size "
         "standard."
      << endl
      << "    -x, --handles_in_snarl      used in conjunction with arguments "
         "source and sink. Will print all the node ids in between source and "
         "sink, inclusive."
      << endl
      << "    -h, --help      print this help info." << endl;
}

int main_normalize(int argc, char **argv) {

  if (argc == 2) {
    help_normalize(argv);
    return 1;
  }

  int max_alignment_size =
      INT_MAX; // default cutoff used to be 200 threads in a snarl.
  string gbwt;
  string snarls;
  string output_gbwt = "normalized.gbwt";
  string normalize_type = "all";
  int source = NULL; // todo: do something other than NULL to avoid the compiler
                     // warnings.
  int sink = NULL;
  bool paths_right_to_left = false;
  bool evaluate = false;
  string snarl_sizes;
  bool snarl_sizes_skip_source_sink = false;
  int max_handle_size = INT_MAX;
  bool handles_in_snarl = false;

  int c;
  optind = 2; // force optind past command positional argument
  while (true) {
    static struct option long_options[] =

        {{"help", no_argument, 0, 'h'},
         {"gbwt", required_argument, 0, 'g'},
         {"snarls", required_argument, 0, 's'},
         {"output_gbwt", required_argument, 0, 'o'},
         {"normalize_type", required_argument, 0, 'n'},
         {"source", required_argument, 0, 'a'},
         {"sink", required_argument, 0, 'b'},
         {"paths_right_to_left", no_argument, 0, 'p'},
         {"max_alignment_size", required_argument, 0, 'm'},
         {"snarl_sizes", required_argument, 0, 'i'},
         {"snarl_sizes_skip_source_sink", no_argument, 0, 'k'},
         {"max_handle_size", required_argument, 0, 'h'},
         {"handles_in_snarl", no_argument, 0, 'x'},
         {0, 0, 0, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "hg:s:o:n:a:b:pm:i:kh:x", long_options,
                    &option_index);

    // Detect the end of the options.
    if (c == -1)
      break;

    switch (c) {

    case 'g':
      gbwt = optarg;
      break;

    case 's':
      snarls = optarg;
      break;

    case 'o':
      output_gbwt = optarg;
      break;

    case 'n':
      normalize_type = optarg;
      break;

    case 'a':
      source = parse<int>(optarg);
      break;

    case 'b':
      sink = parse<int>(optarg);
      break;

    case 'p':
      paths_right_to_left = true;
      break;

    case 'm':
      max_alignment_size = parse<int>(optarg);
      // if max_alignment_size is 0, then that signifies that it should actually
      // be infinite, i.e. that we should not exclude any snarls.
      if (max_alignment_size == 0) {
        max_alignment_size = INT_MAX;
      }
      break;

    case 'i':
      snarl_sizes = optarg;
      // normalize_type = "none";
      break;

    case 'k':
      snarl_sizes_skip_source_sink = true;
      break;

    case 'h':
      max_handle_size = parse<int>(optarg);
      break;

    case 'x':
      handles_in_snarl = true;
      normalize_type = "none";
      break;

    default:
      cerr << "error:[vg normalize] abort" << endl;
      abort();
    }
  }

  // getting graph of any type, except non-mutable graphs (e.g., xg)
  unique_ptr<MutablePathDeletableHandleGraph> graph;
  get_input_file(optind, argc, argv, [&](istream &in) {
    graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
  });

  if (normalize_type != "all" && normalize_type != "one" &&  normalize_type != "none") {
    cerr << "please enter a valid normalize_type: all or one." << endl;
  }

  if (normalize_type == "all" || normalize_type == "one") {
    cerr << "running normalize!" << endl;

    /// Build the gbwt:
    ifstream gbwt_stream;
    gbwt_stream.open(gbwt);

    // Load the GBWT from its container
    unique_ptr<gbwt::GBWT> gbwt;
    gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
    gbwtgraph::GBWTGraph gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);

    std::ifstream snarl_stream;
    string snarl_file = snarls;
    snarl_stream.open(snarl_file);

    if (!snarl_stream) {
      cerr << "error:[vg normalize] Cannot open Snarls file " << snarl_file
           << endl;
      exit(1);
    }
    // Record start time
    auto start = chrono::high_resolution_clock::now();

    algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
        *graph, *gbwt, gbwt_graph, max_alignment_size, max_handle_size);

    if (normalize_type == "all") 
    {

      // gbwt_graph.get_handle()
      
      gbwt::GBWT normalized_gbwt = normalizer.normalize_top_level_snarls(snarl_stream);
      save_gbwt(normalized_gbwt, output_gbwt, true);

  //     //todo: delete this secondary normalize:
  //     algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
  //       *graph, *gbwt, gbwt_graph, max_alignment_size, max_handle_size);
  // // getting graph of any type, except non-mutable graphs (e.g., xg)
  // unique_ptr<MutablePathDeletableHandleGraph> graph;
  // get_input_file(optind, argc, argv, [&](istream &in) {
  //   graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
  // });


    }
    else if (normalize_type == "one")
    {
      if (source == NULL && sink == NULL) 
      {
        cerr << "ERROR: please provide a source and sink for the snarl you "
                "want to normalize."
             << endl;
        return 0;
      }
      vector<int> error_record =
          normalizer.normalize_snarl(source, sink, paths_right_to_left);
      if (!(error_record[0] || error_record[1] || error_record[2] ||
            error_record[3] || error_record[6])) 
      {
        cerr << "snarl starting at " << source << " and ending at " << sink
             << " normalized." << endl;
        cerr << "amount of sequence in normalized snarl before normalization: "
             << error_record[4] << endl;
        cerr << "amount of sequence in normalized snarl after normalization: "
             << error_record[5] << endl;
      } 
      else 
      {
        // todo: make it so it only prints the relevant message:
        cerr << "snarl skipped because...\nthey exceeded the size limit ("
             << error_record[0] << " snarls),\n"
             << "had haplotypes starting/ending in the middle of the snarl ("
             << error_record[1] << "),\n"
             << "the snarl was cyclic (" << error_record[3] << " snarls),\n"
             << " there were handles not connected by the gbwt info ("
             << error_record[2] << " snarls),\n"
             << "the snarl was cyclic (" << error_record[3] << " snarls),\n"
             << "or the snarl was trivial - composed of only one or two nodes ("
             << error_record[6] << " snarls)." << endl;
      }
    }
    // // run test code on all snarls in graph. (non obj-oriented code)
    // disambiguate_top_level_snarls(*graph, gbwt_graph, snarl_stream,
    // max_alignment_size);

    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cerr << "Elapsed time: " << elapsed.count() << " s\n";

    // Save the modified graph
    vg::io::save_handle_graph(graph.get(), std::cout);
  }

  // snarl_analyzer identifies the size of every top-level snarl, outputs in a
  // document specified with format "source\tsink\tsize\n"
  if (snarl_sizes.size() != 0) {
    std::ifstream snarl_stream;
    snarl_stream.open(snarls);
    if (!snarl_stream) {
      cerr << "error:[vg normalize] Cannot open Snarls file " << snarls << endl;
      exit(1);
    }

    algorithms::SnarlAnalyzer sizes = algorithms::SnarlAnalyzer(
        *graph, snarl_stream, snarl_sizes_skip_source_sink);

    sizes.output_snarl_sizes(snarl_sizes);
  }

  if (handles_in_snarl) {
    if (source == NULL && sink == NULL) {
      cerr << "error:[vg normalize] please enter a values for source and sink "
              "to define the snarl."
           << endl;
    } else {
      algorithms::print_handles_in_snarl(*graph, source, sink);
    }
  }

  return 0;
}

// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication",
                               TOOLKIT, main_normalize);