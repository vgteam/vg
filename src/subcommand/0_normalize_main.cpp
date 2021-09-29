// mod.cpp: define the "normalize" subcommand, which realigns snarls to produce
// more efficient representations of snarls.

#include <getopt.h>
#include <omp.h>
#include <unistd.h>

#include <gbwtgraph/gbwtgraph.h>
#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"

#include "subcommand.hpp"

// todo: should be able to remove '../../include/...' and replace with e.g.
// <bdsg/hash...>
#include "../../include/bdsg/hash_graph.hpp"
#include "../../include/vg/io/vpkg.hpp"

#include "../io/save_handle_graph.hpp"

#include "../algorithms/0_oo_normalize_snarls.hpp"
#include "../algorithms/0_snarl_analyzer.hpp"

#include "../snarls.hpp"

#include <chrono> // for high_resolution_clock

using namespace std;
using namespace vg;
using namespace vg::subcommand;

pair<shared_ptr<MutablePathDeletableHandleGraph>, gbwt::GBWT> run_norm(vector<const Snarl *> snarl_roots, int optind, int argc, char** argv, string gbwt_file, string gbwt_graph_file, int max_handle_size, int max_alignment_size){
  // getting graph of any type, except non-mutable graphs (e.g., xg)
  shared_ptr<MutablePathDeletableHandleGraph> graph;
  get_input_file(optind, argc, argv, [&](istream &in) {
    graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
  });


    /// Build the gbwt:
  ifstream gbwt_stream;
  gbwt_stream.open(gbwt_file);

  // Load the GBWT from its container
  unique_ptr<gbwt::GBWT> gbwt;
  gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_stream);
  
  // gbwtgraph::GBWTGraph gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);
  // algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
  // *graph, *gbwt, gbwt_graph, max_handle_size, max_alignment_size);
  
  // gbwtgraph::GBWTGraph gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);
  // save_gbwtgraph(gbwt_graph, gbwt_file+".gg");
  // unique_ptr<gbwtgraph::GBWTGraph> gbwt_graph_p;
  // gbwt_graph_p = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_file+".gg");
  // algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
  // *graph, *gbwt, *gbwt_graph_p, max_handle_size, max_alignment_size);

  unique_ptr<gbwtgraph::GBWTGraph> gbwt_graph;
  if (gbwt_graph_file.size() == 0)
  {
    string gbwt_graph_output_file = gbwt_file + ".gg";
    cerr << "gbwt_graph option is empty. Making new GBWTGraph from gbwt and graph. Saving as " << gbwt_graph_output_file << endl;
    gbwtgraph::GBWTGraph new_gbwt_graph = gbwtgraph::GBWTGraph(*gbwt, *graph);
    save_gbwtgraph(new_gbwt_graph, gbwt_graph_output_file);
    //todo: find way to load gbwtgraph's unique pointer without saving and then reloading file.
    gbwt_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_graph_output_file);
    gbwt_graph->set_gbwt(*gbwt);
  }
  else 
  {
    gbwt_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(gbwt_graph_file);
    gbwt_graph->set_gbwt(*gbwt);
  }

  algorithms::SnarlNormalizer normalizer = algorithms::SnarlNormalizer(
    *graph, *gbwt, *gbwt_graph, max_handle_size, max_alignment_size);

  gbwt::GBWT normalized_gbwt = normalizer.normalize_snarls(snarl_roots);
  return make_pair(graph, normalized_gbwt);
}

//binary search:
void binary_search_norm(vector<const Snarl *> chosen_snarls, int snarl_start, int snarl_end, int optind, int argc, char** argv, string gbwt_file, string gbwt_graph, int max_handle_size, int max_alignment_size) // pass 0 and chosen_snarls.size() for first snarl_start/end. 
{
  cerr << "\nnormalizing snarls between: " << snarl_start << " and " << snarl_end << endl;
  if (snarl_end - snarl_start  == 0)
  {
    return;
  }
  else
  {
    try
    {
      run_norm(chosen_snarls, optind, argc, argv, gbwt_file, gbwt_graph, max_handle_size, max_alignment_size);
    } catch (const std::out_of_range& e) {
      vector<const Snarl *>::const_iterator first = chosen_snarls.begin();
      vector<const Snarl *>::const_iterator mid = chosen_snarls.begin() + chosen_snarls.size()/2;
      vector<const Snarl *>::const_iterator last = chosen_snarls.end();
      vector<const Snarl *> first_half(first, mid);
      vector<const Snarl *> second_half(mid, last);
      cerr << "normalizing snarls between: " << snarl_start << " and " << snarl_end << " failed. Splitting." << endl;
      binary_search_norm(first_half, snarl_start, snarl_start + chosen_snarls.size()/2, optind, argc, argv, gbwt_file, gbwt_graph, max_handle_size, max_alignment_size);
      binary_search_norm(second_half, snarl_start + chosen_snarls.size()/2, snarl_end, optind, argc, argv, gbwt_file, gbwt_graph, max_handle_size, max_alignment_size);
    } catch (...) {
      cerr << "snarls between: " << snarl_start << " and " << snarl_end << " ended successfully." << endl;
    }
  }
  return;
} 

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
      << "    -a, --source        only used when --normalize_type is 'one', or when using -x option. "
         "The source of the single node to be normalized."
      << endl
      << "    -b, --sink      only used when --normalize_type is 'one' or when using -x option. The "
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
      << "    -c, --start_snarl_num      counting starting from 1, determines which snarls from the .snarls.pb will be normalized if normalize_type='all'."
      << endl
      << "    -v, --end_snarl_num      counting starting from 1, determines which snarls from the .snarls.pb will be normalized if normalize_type='all'."
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
  string gbwt_graph;
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
  int max_handle_size = 32;
  bool handles_in_snarl = false;
  //todo: note: options mostly for debugging:
  int start_snarl_num = 0;
  int end_snarl_num = 0;


  int c;
  optind = 2; // force optind past command positional argument
  while (true) {
    static struct option long_options[] =

        {{"help", no_argument, 0, 'h'},
         {"gbwt", required_argument, 0, 'g'},
         {"gbwt_graph", required_argument, 0, 'r'},
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
         {"start_snarl_num", required_argument, 0, 'c'},
         {"end_snarl_num", required_argument, 0, 'v'},
         {0, 0, 0, 0}};

    int option_index = 0;
    c = getopt_long(argc, argv, "hg:r:s:o:n:a:b:pm:i:kh:xc:v:", long_options,
                    &option_index);

    // Detect the end of the options.
    if (c == -1)
      break;

    switch (c) {

    case 'g':
      gbwt = optarg;
      break;

    case 'r':
      gbwt_graph = optarg;
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

    case 'c':
      start_snarl_num = parse<int>(optarg);
      break;

    case 'v':
      end_snarl_num = parse<int>(optarg);
      break;

    default:
      cerr << "error:[vg normalize] abort" << endl;
      abort();
    }
  }


  
  if (normalize_type != "all" && normalize_type != "none") 
  {
    cerr << "please enter a valid normalize_type: all or one." << endl;
  }

  if (normalize_type == "all") 
  {
    // prep snarl_roots:
    std::ifstream snarl_stream;
    string snarl_file = snarls;
    snarl_stream.open(snarl_file);
    if (!snarl_stream) {
      cerr << "error:[vg normalize] Cannot open Snarls file " << snarl_file
            << endl;
      exit(1);
    }

    SnarlManager *snarl_manager = new SnarlManager(snarl_stream);
    vector<const Snarl *> snarl_roots = snarl_manager->top_level_snarls();
      
    cerr << "snarl_roots.size()" << snarl_roots.size() << endl;

    //todo: use for debugging:
    /// for select snarl:
    // vector<const Snarl *> chosen_snarls(snarl_roots.begin() + 1370135, snarl_roots.begin() + 1370136); 
    // run_norm(chosen_snarls, optind, argc, argv, gbwt, max_handle_size, max_alignment_size);
    /// for binary search:
    // binary_search_norm(snarl_roots, 0, snarl_roots.size(), optind, argc, argv, gbwt, max_handle_size, max_alignment_size);

    
    cerr << "running normalize!" << endl;

    // Record start time
    auto start = chrono::high_resolution_clock::now();
    pair<shared_ptr<MutablePathDeletableHandleGraph>, gbwt::GBWT> output;

    // unique_ptr<MutablePathDeletableHandleGraph> graph;
    // gbwt::GBWT normalized_gbwt;
    // pair<MutablePathDeletableHandleGraph, gbwt::GBWT> output = make_pair(*graph, normalized_gbwt);
    
    //standard, normalize all snarls:
    if (start_snarl_num == 0 && end_snarl_num == 0)
    {
      output = run_norm(snarl_roots, optind, argc, argv, gbwt, gbwt_graph, max_handle_size, max_alignment_size);
    }
    //normalize select snarls:
    else
    {
      //check that start and end are valid:
      if (start_snarl_num >= end_snarl_num)
      {
        cerr << "error:[vg normalize] start_snarl_num >= end_snarl_num."  << endl;
        exit(1);
      }
      else if (end_snarl_num > snarl_roots.size())
      {
        cerr << "WARNING:[vg normalize] end_snarl_num greater than snarl_roots.size(). Will normalize starting at start_snarl_num and ending at last snarl in snarl_roots." << endl;
      }
      
      vector<const Snarl *> chosen_snarls(snarl_roots.begin() + start_snarl_num, snarl_roots.begin() + end_snarl_num); 

      cerr << "of snarl selection that is " << chosen_snarls.size() << " long, first snarl selected has source: " << (*chosen_snarls.front()).start().node_id() << " and sink: " << (*chosen_snarls.front()).end().node_id() << endl; 
      output = run_norm(chosen_snarls, optind, argc, argv, gbwt, gbwt_graph, max_handle_size, max_alignment_size);
    }
    save_gbwt(output.second, output_gbwt, true);

    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;
    cerr << "Elapsed time: " << elapsed.count() << " s" << endl;

    // Save the modified graph
    vg::io::save_handle_graph(output.first.get(), std::cout);
  }

  // snarl_analyzer identifies the size of every top-level snarl, outputs in a
  // document specified with format "source\tsink\tsize\n"
  if (snarl_sizes.size() != 0) {
    shared_ptr<MutablePathDeletableHandleGraph> graph;
    get_input_file(optind, argc, argv, [&](istream &in) {
      graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
    });
    
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

  if (handles_in_snarl) 
  {
    if (source == NULL && sink == NULL) 
    {
      cerr << "error:[vg normalize] please enter a values for source and sink "
              "to define the snarl."
           << endl;
    } 
    else 
    {
      shared_ptr<MutablePathDeletableHandleGraph> graph;
      get_input_file(optind, argc, argv, [&](istream &in) {
        graph = vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(in);
      });
      algorithms::print_handles_in_snarl(*graph, source, sink);
    }
  }

  return 0;
}

// Register subcommand
static Subcommand vg_normalize("normalize",
                               "edit snarls to reduce information duplication",
                               TOOLKIT, main_normalize) ;