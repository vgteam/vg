#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include "gcsa.h"
// From gcsa2
#include "files.h"
#include "json2pb.h"
#include "vg.hpp"
#include "vg.pb.h"
#include "vg_set.hpp"
#include "index.hpp"
#include "mapper.hpp"
#include "Variant.h"
#include "Fasta.h"
#include "stream.hpp"
#include "alignment.hpp"
#include "convert.hpp"
#include "pileup.hpp"
#include "caller.hpp"
#include "deconstructor.hpp"
#include "vectorizer.hpp"
#include "google/protobuf/stubs/common.h"
#include "progress_bar.hpp"

using namespace std;
using namespace google::protobuf;
using namespace vg;

void help_validate(char** argv) {
    cerr << "usage: " << argv[0] << " validate [options] graph" << endl
        << "Validate the graph." << endl
        << endl
        << "options:" << endl
        << "    default: check all aspects of the graph, if options are specified do only those" << endl
        << "    -n, --nodes    verify that we have the expected number of nodes" << endl
        << "    -e, --edges    verify that the graph contains all nodes that are referred to by edges" << endl
        << "    -p, --paths    verify that contiguous path segments are connected by edges" << endl
        << "    -o, --orphans  verify that all nodes have edges" << endl;
}

int main_validate(int argc, char** argv) {

    if (argc <= 2) {
        help_validate(argv);
        return 1;
    }

    bool check_nodes = false;
    bool check_edges = false;
    bool check_orphans = false;
    bool check_paths = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"nodes", no_argument, 0, 'n'},
            {"edges", no_argument, 0, 'e'},
            {"paths", no_argument, 0, 'o'},
            {"orphans", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hneop",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'n':
                check_nodes = true;
                break;

            case 'e':
                check_edges = true;
                break;

            case 'o':
                check_orphans = true;
                break;

            case 'p':
                check_paths = true;
                break;

            case 'h':
            case '?':
                help_validate(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    // if we chose a specific subset, do just them
    if (check_nodes || check_edges || check_orphans || check_paths) {
        if (graph->is_valid(check_nodes, check_edges, check_orphans, check_paths)) {
            return 0;
        } else {
            return 1;
        }
        // otherwise do everything
    } else if (graph->is_valid()) {
        return 0;
    } else {
        return 1;
    }
}

void help_vectorize(char** argv){
    cerr << "usage: " << argv[0] << " vectorize [options] -x <index.xg> <alignments.gam>" << endl
        << "Vectorize a set of alignments to a variety of vector formats." << endl
        << endl
        << "options: " << endl
        << "  -x --xg-name FILE  An xg index for the graph of interest" << endl
        << "  -f --format        Tab-delimit output so it can be used in R." << endl
        << "  -A --annotate      Create a header with each node/edge's name and a column with alignment names." << endl
        << endl;
}

int main_vectorize(int argc, char** argv){

    string xg_name;
    bool format = false;
    bool show_header = false;
    bool map_alns = false;
    bool annotate = false;

    if (argc <= 2) {
        help_vectorize(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"annotate", no_argument, 0, 'A'},
            {"xg-name", required_argument,0, 'x'},
            {"threads", required_argument, 0, 't'},
            {"format", no_argument, 0, 'f'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "Ahfx:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case '?':
            case 'h':
                help_vectorize(argv);
                return 1;
            case 'x':
                xg_name = optarg;
                break;
            case 'f':
                format = true;
                break;
            case 'A':
                annotate = true;
                format = true;
                break;
            default:
                abort();
        }
    }

    xg::XG* xindex;
    if (!xg_name.empty()) {
        ifstream in(xg_name);
        xindex = new xg::XG(in);
    }
    else{
        cerr << "No XG index given. An XG index must be provided." << endl;
        exit(1);
    }

    Vectorizer vz(xindex);
    string alignment_file = argv[optind];

    //Generate a 1-hot coverage vector for graph entities.
    function<void(Alignment&)> lambda = [&vz, format](Alignment& a){
      //vz.add_bv(vz.alignment_to_onehot(a));
      //vz.add_name(a.name());
      bit_vector v = vz.alignment_to_onehot(a);
      if (format){
        cout << a.name() << "\t";
        cout << vz.format(v) << endl;
      }
      else{
        cout << v << endl;
      }

    };
    if (alignment_file == "-"){
        stream::for_each_parallel(cin, lambda);
    }
    else{
        ifstream in;
        in.open(alignment_file);
        if (in.good()){
            stream::for_each(in, lambda);
        }
    }



    //TODO handle custom scores settings.

    vz.emit(cout, format, annotate);



    return 0;
}

void help_compare(char** argv) {
    cerr << "usage: " << argv[0] << " compare [options] graph1 graph2" << endl
        << "Compare kmer sets of two graphs" << endl
        << endl
        << "options:" << endl
        << "    -d, --db-name1 FILE  use this db for graph1 (defaults to <graph1>.index/)" << endl
        << "    -e, --db-name2 FILE  use this db for graph2 (defaults to <graph1>.index/)" << endl
        << "    -t, --threads N      number of threads to use" << endl;
}

int main_compare(int argc, char** argv) {

    if (argc <= 3) {
        help_compare(argv);
        return 1;
    }

    string db_name1;
    string db_name2;
    int num_threads = 1;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"db-name1", required_argument, 0, 'd'},
            {"db-name2", required_argument, 0, 'e'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hd:e:t:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'd':
                db_name1 = optarg;
                break;

            case 'e':
                db_name2 = optarg;
                break;

            case 't':
                num_threads = atoi(optarg);
                break;

            case 'h':
            case '?':
                help_compare(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    omp_set_num_threads(num_threads);

    string file_name1 = argv[optind++];
    string file_name2 = argv[optind];

    if (db_name1.empty()) {
        db_name1 = file_name1;
    }
    if (db_name2.empty()) {
        db_name2 = file_name2;
    }

    // Note: only supporting rocksdb index for now.

    Index index1;
    index1.open_read_only(db_name1);

    Index index2;
    index2.open_read_only(db_name2);

    pair<int64_t, int64_t> index1_vs_index2;
    pair<int64_t, int64_t> index2_vs_index1;

    // Index::compare is not parallel, but at least we can do the
    // two directions at the same time...
#pragma omp parallel sections
    {
#pragma omp section
        {
            index1_vs_index2 = index1.compare_kmers(index2);
        }
#pragma omp section
        {
            index2_vs_index1 = index2.compare_kmers(index1);
        }
    }
    {// <-- for emacs
        assert(index1_vs_index2.first == index2_vs_index1.first);

        int64_t db1_count = index1_vs_index2.first + index1_vs_index2.second;
        int64_t db2_count = index2_vs_index1.first + index2_vs_index1.second;
        int64_t db1_only = index1_vs_index2.second;
        int64_t db2_only = index2_vs_index1.second;
        int64_t db1_and_db2 = index1_vs_index2.first;
        int64_t db1_or_db2 = db1_only + db2_only + db1_and_db2;

        cout << "{\n"
            << "\"db1_path\": " << "\"" << db_name1 << "\"" << ",\n"
            << "\"db2_path\": " << "\"" << db_name2 << "\"" << ",\n"
            << "\"db1_total\": " << db1_count << ",\n"
            << "\"db2_total\": " << db2_count << ",\n"
            << "\"db1_only\": " << db1_only << ",\n"
            << "\"db2_only\": " << db2_only << ",\n"
            << "\"intersection\": " << db1_and_db2 << ",\n"
            << "\"union\": " << db1_or_db2 << "\n"
            << "}" << endl;
    }
    return 0;
}

void help_call(char** argv) {
    cerr << "usage: " << argv[0] << " call [options] <graph.vg> <pileup.vgpu> > sample_graph.vg" << endl
        << "Compute SNPs from pilup data (prototype! for evaluation only). " << endl
        << endl
        << "options:" << endl
        << "    -d, --min_depth         minimum depth of pileup (default=" << Caller::Default_min_depth <<")" << endl
        << "    -e, --max_depth         maximum depth of pileup (default=" << Caller::Default_max_depth <<")" << endl
        << "    -s, --min_support       minimum number of reads required to support snp (default=" << Caller::Default_min_support <<")" << endl
        << "    -r, --het_prior         prior for heterozygous genotype (default=" << Caller::Default_het_prior <<")" << endl
        << "    -q, --default_read_qual phred quality score to use if none found in the pileup (default="
        << (int)Caller::Default_default_quality << ")" << endl
        << "    -l, --leave_uncalled    leave un-called graph regions in output" << endl
        << "    -j, --json              output in JSON" << endl
        << "    -p, --progress          show progress" << endl
        << "    -t, --threads N         number of threads to use" << endl;
}

int main_call(int argc, char** argv) {

    if (argc <= 3) {
        help_call(argv);
        return 1;
    }

    double het_prior = Caller::Default_het_prior;
    int min_depth = Caller::Default_min_depth;
    int max_depth = Caller::Default_max_depth;
    int min_support = Caller::Default_min_support;
    int default_read_qual = Caller::Default_default_quality;
    bool leave_uncalled = false;
    bool output_json = false;
    bool show_progress = false;
    int thread_count = 1;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
        {
            {"min_depth", required_argument, 0, 'd'},
            {"max_depth", required_argument, 0, 'e'},
            {"min_support", required_argument, 0, 's'},
            {"default_read_qual", required_argument, 0, 'q'},
            {"leave_uncalled", no_argument, 0, 'l'},
            {"json", no_argument, 0, 'j'},
            {"progress", no_argument, 0, 'p'},
            {"het_prior", required_argument, 0, 'r'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:e:s:q:ljpr:t:",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 'd':
                min_depth = atoi(optarg);
                break;
            case 'e':
                max_depth = atoi(optarg);
                break;
            case 's':
                min_support = atoi(optarg);
                break;
            case 'q':
                default_read_qual = atoi(optarg);
                break;
            case 'l':
                leave_uncalled = true;
                break;
            case 'j':
                output_json = true;
                break;
            case 'p':
                show_progress = true;
                break;
            case 'r':
                het_prior = atof(optarg);
                break;
            case 't':
                thread_count = atoi(optarg);
                break;
            case 'h':
            case '?':
                /* getopt_long already printed an error message. */
                help_call(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }
    omp_set_num_threads(thread_count);
    thread_count = get_thread_count();

    // read the graph
    if (optind >= argc) {
        help_call(argv);
        return 1;
    }
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }
    VG* graph;
    string graph_file_name = argv[optind++];
    if (graph_file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(graph_file_name.c_str());
        if (!in) {
            cerr << "error: input file " << graph_file_name << " not found." << endl;
            exit(1);
        }
        graph = new VG(in);
    }

    // setup pileup stream
    if (optind >= argc) {
        help_call(argv);
        return 1;
    }
    string pileup_file_name = argv[optind];
    istream* pileup_stream = NULL;
    ifstream in;
    if (pileup_file_name == "-") {
        if (graph_file_name == "-") {
            cerr << "error: graph and pileup can't both be from stdin." << endl;
            exit(1);
        }
        pileup_stream = &std::cin;
    } else {
        in.open(pileup_file_name);
        if (!in) {
            cerr << "error: input file " << pileup_file_name << " not found." << endl;
            exit(1);
        }
        pileup_stream = &in;
    }

    // compute the variants.
    if (show_progress) {
        cerr << "Computing variants" << endl;
    }
    Caller caller(graph,
            het_prior, min_depth, max_depth, min_support,
            Caller::Default_min_frac, Caller::Default_min_likelihood,
            leave_uncalled, default_read_qual);

    function<void(NodePileup&)> lambda = [&caller](NodePileup& pileup) {
        caller.call_node_pileup(pileup);
    };
    stream::for_each(*pileup_stream, lambda);

    // map the edges from original graph
    if (show_progress) {
        cerr << "Mapping edges into call graph" << endl;
    }
    caller.update_call_graph();

    // write the call graph
    if (show_progress) {
        cerr << "Writing call graph" << endl;
    }
    caller.write_call_graph(cout, output_json);

    return 0;
}

void help_pileup(char** argv) {
    cerr << "usage: " << argv[0] << " pileup [options] <graph.vg> <alignment.gam> > out.vgpu" << endl
        << "Calculate pileup for each position in graph and output in VG Pileup format (list of protobuf NodePileups)." << endl
        << endl
        << "options:" << endl
        << "    -j, --json           output in JSON" << endl
        << "    -p, --progress       show progress" << endl
        << "    -t, --threads N      number of threads to use" << endl;
}

int main_pileup(int argc, char** argv) {

    if (argc <= 3) {
        help_pileup(argv);
        return 1;
    }

    bool output_json = false;
    bool show_progress = false;
    int thread_count = 1;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
        {
            {"json", required_argument, 0, 'j'},
            {"progress", required_argument, 0, 'p'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "jpt:",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 'j':
                output_json = true;
                break;
            case 'p':
                show_progress = true;
                break;
            case 't':
                thread_count = atoi(optarg);
                break;
            case 'h':
            case '?':
                /* getopt_long already printed an error message. */
                help_pileup(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }
    omp_set_num_threads(thread_count);
    thread_count = get_thread_count();

    // read the graph
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }
    VG* graph;
    string graph_file_name = argv[optind++];
    if (graph_file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(graph_file_name.c_str());
        if (!in) {
            cerr << "error: input file " << graph_file_name << " not found." << endl;
            exit(1);
        }
        graph = new VG(in);
    }

    // setup alignment stream
    string alignments_file_name = argv[optind++];
    istream* alignment_stream = NULL;
    ifstream in;
    if (alignments_file_name == "-") {
        if (graph_file_name == "-") {
            cerr << "error: graph and alignments can't both be from stdin." << endl;
            exit(1);
        }
        alignment_stream = &std::cin;
    } else {
        in.open(alignments_file_name);
        if (!in) {
            cerr << "error: input file " << graph_file_name << " not found." << endl;
            exit(1);
        }
        alignment_stream = &in;
    }

    // compute the pileups.
    if (show_progress) {
        cerr << "Computing pileups" << endl;
    }
    vector<Pileups> pileups(thread_count);
    function<void(Alignment&)> lambda = [&pileups, &graph](Alignment& aln) {
        int tid = omp_get_thread_num();
        pileups[tid].compute_from_alignment(*graph, aln);
    };
    stream::for_each_parallel(*alignment_stream, lambda);

    // single-threaded (!) merge
    if (show_progress && pileups.size() > 1) {
        cerr << "Merging pileups" << endl;
    }
    for (int i = 1; i < pileups.size(); ++i) {
        pileups[0].merge(pileups[i]);
    }
    // spit out the pileup
    if (show_progress) {
        cerr << "Writing pileups" << endl;
    }
    if (output_json == false) {
        pileups[0].write(std::cout);
    } else {
        pileups[0].to_json(std::cout);
    }

    delete graph;
    return 0;
}

void help_msga(char** argv) {
    cerr << "usage: " << argv[0] << " msga [options] >graph.vg" << endl
         << "Multiple sequence / graph aligner." << endl
         << endl
         << "options:" << endl
         << "inputs:" << endl
         << "    -f, --from FILE         use sequences in (fasta) FILE" << endl
         << "    -n, --name NAME         include this sequence" << endl
         << "                             (If any --name is specified, use only" << endl
         << "                              specified sequences from FASTA files.)" << endl
         << "    -b, --base NAME         use this sequence as the graph basis if graph is empty" << endl
         << "    -s, --seq SEQUENCE      literally include this sequence" << endl
         << "    -g, --graph FILE        include this graph" << endl
         << "mem mapping:" << endl
         << "    -L, --min-mem-length N  ignore SMEMs shorter than this length (default: 0/unset)" << endl
         << "    -Y, --max-mem-length N  ignore SMEMs longer than this length by stopping backward search (default: 0/unset)" << endl
         << "    -H, --hit-max N         SMEMs which have >N hits in our index (default: 100)" << endl
         << "    -c, --context-depth N   follow this many steps out from each subgraph for alignment (default: 7)" << endl
         << "    -T, --thread-ex N       cluster nodes when successive ids are within this distance (default: 10)" << endl
         << "    -P, --min-score N       accept alignment only if the normalized alignment score is >N (default: 0.75)" << endl
         << "    -B, --band-width N      use this bandwidth when mapping (default: 256)" << endl
         << "    -G, --greedy-accept     if a tested alignment achieves -S score/bp don't try clusters with fewer hits" << endl
         << "    -S, --accept-score N    accept early alignment if the normalized alignment score is > N and -G is set" << endl
         << "    -M, --max-attempts N    only attempt the N best subgraphs ranked by SMEM support (default: 10)" << endl
         << "    -q, --max-target-x N    skip cluster subgraphs with length > N*read_length (default: 100; 0=unset)" << endl
         << "    -I, --max-multimaps N   if N>1, keep N best mappings of each band, resolve alignment by DP (default: 1)" << endl
         << "index generation:" << endl
         << "    -K, --idx-kmer-size N   use kmers of this size for building the GCSA indexes (default: 16)" << endl
         << "    -O, --idx-no-recomb     index only embedded paths, not recombinations of them" << endl
         << "    -E, --idx-edge-max N    reduce complexity of graph indexed by GCSA using this edge max (default: off)" << endl
         << "    -Q, --idx-prune-subs N  prune subgraphs shorter than this length from input graph to GCSA (default: off)" << endl
         << "    -m, --node-max N        chop nodes to be shorter than this length (default: 2* --idx-kmer-size)" << endl
         << "    -X, --idx-doublings N   use this many doublings when building the GCSA indexes (default: 2)" << endl
         << "graph normalization:" << endl
         << "    -N, --normalize         normalize the graph after assembly" << endl
         << "    -z, --allow-nonpath     don't remove parts of the graph that aren't in the paths of the inputs" << endl
         << "    -D, --debug             print debugging information about construction to stderr" << endl
         << "    -A, --debug-align       print debugging information about alignment to stderr" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << endl
         << "Construct a multiple sequence alignment from all sequences in the" << endl
         << "input fasta-format files, graphs, and sequences. Uses the MEM mapping algorithm." << endl
         << endl
         << "Emits the resulting MSA as a (vg-format) graph." << endl;
}

int main_msga(int argc, char** argv) {

    if (argc == 2) {
        help_msga(argv);
        return 1;
    }

    vector<string> fasta_files;
    set<string> seq_names;
    vector<string> sequences;
    vector<string> graph_files;
    string base_seq_name;
    int idx_kmer_size = 16;
    int idx_doublings = 2;
    int best_clusters = 0;
    int hit_max = 100;
    int max_attempts = 10;
    // if we set this above 1, we use a dynamic programming process to determine the
    // optimal alignment through a series of bands based on a proximity metric
    int max_multimaps = 1;
    // if this is set too low, we may miss optimal alignments
    int context_depth = 7;
    // same here; initial clustering
    int thread_extension = 10;
    float min_norm_score = 0.75;
    int band_width = 256;
    size_t doubling_steps = 2;
    bool debug = false;
    bool debug_align = false;
    size_t node_max = 0;
    int alignment_threads = 1;
    int edge_max = 0;
    int subgraph_prune = 0;
    bool normalize = false;
    bool allow_nonpath = false;
    int iter_max = 1;
    int max_mem_length = 0;
    int min_mem_length = 0;
    bool greedy_accept = false;
    float accept_score = 0;
    int max_target_factor = 100;
    bool idx_path_only = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                {"help", no_argument, 0, 'h'},
                {"from", required_argument, 0, 'f'},
                {"name", required_argument, 0, 'n'},
                {"seq", required_argument, 0, 's'},
                {"graph", required_argument, 0, 'g'},
                {"base", required_argument, 0, 'b'},
                {"idx-kmer-size", required_argument, 0, 'K'},
                {"idx-no-recomb", no_argument, 0, 'O'},
                {"idx-doublings", required_argument, 0, 'X'},
                {"band-width", required_argument, 0, 'B'},
                {"debug", no_argument, 0, 'D'},
                {"debug-align", no_argument, 0, 'A'},
                {"context-depth", required_argument, 0, 'c'},
                {"min-score", required_argument, 0, 'P'},
                {"idx-edge-max", required_argument, 0, 'E'},
                {"idx-prune-subs", required_argument, 0, 'Q'},
                {"normalize", no_argument, 0, 'N'},
                {"allow-nonpath", no_argument, 0, 'z'},
                {"min-mem-length", required_argument, 0, 'L'},
                {"max-mem-length", required_argument, 0, 'Y'},
                {"hit-max", required_argument, 0, 'H'},
                {"threads", required_argument, 0, 't'},
                {"node-max", required_argument, 0, 'm'},
                {"greedy-accept", no_argument, 0, 'G'},
                {"accept-score", required_argument, 0, 'S'},
                {"max-attempts", required_argument, 0, 'M'},
                {"thread-ex", required_argument, 0, 'T'},
                {"max-target-x", required_argument, 0, 'q'},
                {"max-multimaps", required_argument, 0, 'I'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hf:n:s:g:b:K:X:B:DAc:P:E:Q:NzI:L:Y:H:t:m:GS:M:T:q:OI:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {


        case 'L':
            min_mem_length = atoi(optarg);
            break;

        case 'Y':
            max_mem_length = atoi(optarg);
            break;

        case 'H':
            hit_max = atoi(optarg);
            break;

        case 'I':
            max_multimaps = atoi(optarg);
            break;

        case 'q':
            max_target_factor = atoi(optarg);
            break;

        case 'M':
            max_attempts = atoi(optarg);
            break;

        case 'T':
            thread_extension = atoi(optarg);
            break;

        case 'G':
            greedy_accept = true;
            break;

        case 'S':
            accept_score = atof(optarg);
            break;

        case 'c':
            context_depth = atoi(optarg);
            break;


        case 'f':
            fasta_files.push_back(optarg);
            break;

            case 'n':
                seq_names.insert(optarg);
                break;

            case 's':
                sequences.push_back(optarg);
                break;

            case 'b':
                base_seq_name = optarg;
                break;

            case 'g':
                if (graph_files.size() != 0) {
                    cerr << "[vg msga] Error: graph-graph alignment is not yet implemented." << endl
                        << "We can only use one input graph." << endl;
                    return 1;
                }
                graph_files.push_back(optarg);
                break;

        case 'B':
            band_width = atoi(optarg);
            break;

            case 'D':
                debug = true;
                break;

            case 'A':
                debug_align = true;
                break;

            case 'X':
                doubling_steps = atoi(optarg);
                break;

            case 'K':
                idx_kmer_size = atoi(optarg);
                break;


        case 'O':
            idx_path_only = true;
            break;

        case 'm':
            node_max = atoi(optarg);
            break;

        case 'N':
            normalize = true;
            break;


        case 'P':
            min_norm_score = atof(optarg);
            break;

            case 't':
                omp_set_num_threads(atoi(optarg));
                alignment_threads = atoi(optarg);
                break;

            case 'Q':
                subgraph_prune = atoi(optarg);
                break;

            case 'E':
                edge_max = atoi(optarg);
                break;

            case 'z':
                allow_nonpath = true;
                break;

            case 'h':
            case '?':
                help_msga(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    // build the graph or read it in from input
    VG* graph;
    if (graph_files.size() == 1) {
        string file_name = graph_files.front();
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
    } else {
        graph = new VG;
    }

    // we should chop up the inputs into bits
    // then run a kind of alignment/overlap assembly on them
    // to generate the new graph/msa
    // TODO refactor into class

    // map from name to sequence, just a transformation of FASTA records
    map<string, string> strings;

    // open the fasta files, read in the sequences
    vector<string> names_in_order;

    for (auto& fasta_file_name : fasta_files) {
        FastaReference ref;
        ref.open(fasta_file_name);
        if (debug) cerr << "loading " << fasta_file_name << endl;
        for (auto& name : ref.index->sequenceNames) {
            // only use the sequence if we have whitelisted it
            string seq = ref.getSequence(name);
            strings[name] = seq;
        }
    }

    // give a null label to sequences passed on the command line
    // thus far there is no way to specify these (use fasta instead)
    for (auto& s : sequences) {
        strings[sha1head(s, 8)] = s;
    }

    // align, include, repeat

    if (debug) cerr << "preparing initial graph" << endl;

    // if our graph is empty, we need to take the first sequence and build a graph from it
    if (graph->empty()) {
        // what's the first sequence?
        if (base_seq_name.empty()) {
            graph->create_node(strings.begin()->second);
        } else {
            // we specified one we wanted to use as the first
            graph->create_node(strings[base_seq_name]);
        }
    }

    size_t max_query_size = pow(2, doubling_steps) * idx_kmer_size;
    // limit max node size
    if (!node_max) node_max = 2*idx_kmer_size;
    graph->dice_nodes(node_max);
    graph->sort();
    graph->compact_ids();

    // questions:
    // should we preferentially use sequences from fasta files in the order they were given?
    // (considering this a todo)
    // reverse complement?
    Mapper* mapper = nullptr;
    gcsa::GCSA* gcsaidx = nullptr;
    gcsa::LCPArray* lcpidx = nullptr;
    xg::XG* xgidx = nullptr;
    size_t iter = 0;

    auto rebuild = [&](VG* graph) {
        //stringstream s; s << iter++ << ".vg";

        if (mapper) delete mapper;
        if (xgidx) delete xgidx;
        if (gcsaidx) delete gcsaidx;
        if (lcpidx) delete lcpidx;

        if (debug) cerr << "building xg index" << endl;
        xgidx = new xg::XG(graph->graph);
        if (debug) cerr << "building GCSA2 index" << endl;
        if (edge_max) {
            VG gcsa_graph = *graph; // copy the graph
            // remove complex components
            gcsa_graph.prune_complex_with_head_tail(idx_kmer_size, edge_max);
            if (subgraph_prune) gcsa_graph.prune_short_subgraphs(subgraph_prune);
            // then index
            gcsa_graph.build_gcsa_lcp(gcsaidx, lcpidx, idx_kmer_size, idx_path_only, false, doubling_steps);
        } else {
            // if no complexity reduction is requested, just build the index
            graph->build_gcsa_lcp(gcsaidx, lcpidx, idx_kmer_size, idx_path_only, false, doubling_steps);
        }
        mapper = new Mapper(xgidx, gcsaidx, lcpidx);
        { // set mapper variables
            mapper->debug = debug_align;
            mapper->context_depth = context_depth;
            mapper->thread_extension = thread_extension;
            mapper->max_attempts = max_attempts;
            mapper->min_norm_score = min_norm_score;
            mapper->alignment_threads = alignment_threads;
            mapper->max_mem_length = max_mem_length;
            mapper->min_mem_length = min_mem_length;
            mapper->hit_max = hit_max;
            mapper->greedy_accept = greedy_accept;
            mapper->max_target_factor = max_target_factor;
            mapper->max_multimaps = max_multimaps;
            if (accept_score) mapper->accept_norm_score = accept_score;
        }
    };

    // set up the graph for mapping
    rebuild(graph);

    // todo restructure so that we are trying to map everything
    // add alignment score/bp bounds to catch when we get a good alignment
    for (auto& sp : strings) {
        bool incomplete = true; // complete when we've fully included the sequence set
        int iter = 0;
        auto& name = sp.first;
        //graph->serialize_to_file("msga-pre-" + name + ".vg");
        while (incomplete && iter++ < iter_max) {
            stringstream s; s << iter; string iterstr = s.str();
            if (debug) cerr << name << ": adding to graph" << iter << endl;
            vector<Path> paths;
            vector<Alignment> alns;
            int j = 0;
            // TODO we should use just one seq
            auto& seq = sp.second;
            // align to the graph
            if (debug) cerr << name << ": aligning sequence of " << seq.size() << "bp against " <<
                           graph->node_count() << " nodes" << endl;
            Alignment aln = simplify(mapper->align(seq, 0, 0, band_width));
            auto aln_seq = graph->path_string(aln.path());
            if (aln_seq != seq) {
                cerr << "[vg msga] alignment corrupted, failed to obtain correct banded alignment (alignment seq != input seq)" << endl;
                cerr << "expected " << seq << endl;
                cerr << "got      " << aln_seq << endl;
                ofstream f(name + "-failed-alignment-" + convert(j) + ".gam");
                stream::write(f, 1, (std::function<Alignment(uint64_t)>)([&aln](uint64_t n) { return aln; }));
                f.close();
                graph->serialize_to_file(name + "-corrupted-alignment.vg");
                exit(1);
            }
            alns.push_back(aln);
            //if (debug) cerr << pb2json(aln) << endl; // huge in some cases
            paths.push_back(aln.path());
            paths.back().set_name(name); // cache name to trigger inclusion of path elements in graph by edit

            /*
               ofstream f(name + "-pre-edit-" + convert(j) + ".gam");
               stream::write(f, 1, (std::function<Alignment(uint64_t)>)([&aln](uint64_t n) { return aln; }));
               f.close();
               */

            ++j;

            // now take the alignment and modify the graph with it
            if (debug) cerr << name << ": editing graph" << endl;
            //graph->serialize_to_file(name + "-pre-edit.vg");
            graph->edit(paths);
            //if (!graph->is_valid()) cerr << "invalid after edit" << endl;
            //graph->serialize_to_file(name + "-immed-post-edit.vg");
            graph->dice_nodes(node_max);
            //if (!graph->is_valid()) cerr << "invalid after dice" << endl;
            //graph->serialize_to_file(name + "-post-dice.vg");
            if (debug) cerr << name << ": sorting and compacting ids" << endl;
            graph->sort();
            //if (!graph->is_valid()) cerr << "invalid after sort" << endl;
            graph->compact_ids(); // xg can't work unless IDs are compacted.
            //if (!graph->is_valid()) cerr << "invalid after compact" << endl;

            // the edit needs to cut nodes at mapping starts and ends
            // thus allowing paths to be included that map directly to entire nodes
            // XXX

            // check that all is well
            //graph->serialize_to_file(name + "-pre-index.vg");
            rebuild(graph);

            // verfy validity of path
            bool is_valid = graph->is_valid();
            auto path_seq = graph->path_string(graph->paths.path(name));
            incomplete = !(path_seq == seq) || !is_valid;
            if (incomplete) {
                cerr << "[vg msga] failed to include alignment, retrying " << endl
                    << "expected " << seq << endl
                    << "got " << path_seq << endl
                    << pb2json(aln.path()) << endl
                    << pb2json(graph->paths.path(name)) << endl;
                graph->serialize_to_file(name + "-post-edit.vg");
            }
        }
        // if (debug && !graph->is_valid()) cerr << "graph is invalid" << endl;
        if (incomplete && iter >= iter_max) {
            cerr << "[vg msga] Error: failed to include path " << name << endl;
            exit(1);
        }
    }

    // auto include_paths = [&mapper,
    //      kmer_size,
    //      kmer_stride,
    //      band_width,
    //      debug,
    //      &strings](VG* graph) {
    //          // include the paths in the graph
    //          if (debug) cerr << "including paths" << endl;
    //          for (auto& group : strings) {
    //              auto& name = group.first;
    //              if (debug) cerr << name << ": tracing path through graph" << endl;
    //              auto& seq = group.second;
    //              if (debug) cerr << name << ": aligning sequence of " << seq.size() << "bp" << endl;
    //              Alignment aln = mapper->align(seq, kmer_size, kmer_stride, band_width);
    //              //if (debug) cerr << "alignment score: " << aln.score() << endl;
    //              aln.mutable_path()->set_name(name);
    //              //if (debug) cerr << "alignment: " << pb2json(aln) << endl;
    //              // todo simplify in the mapper itself when merging the banded bits
    //              if (debug) cerr << name << ": labeling" << endl;
    //              graph->include(aln.path());
    //              // now repeat back the path
    //          }
    //      };

    if (normalize) {
        if (debug) cerr << "normalizing graph" << endl;
        graph->remove_non_path();
        graph->normalize();
        graph->dice_nodes(node_max);
        graph->sort();
        graph->compact_ids();
        if (!graph->is_valid()) {
            cerr << "[vg msga] warning! graph is not valid after normalization" << endl;
        }
    }

    // remove nodes in the graph that have no assigned paths
    // this should be pretty minimal now that we've made one iteration
    if (!allow_nonpath) {
        graph->remove_non_path();
    }

    // finally, validate the included paths
    set<string> failures;
    for (auto& sp : strings) {
        auto& name = sp.first;
        auto& seq = sp.second;
        if (seq != graph->path_string(graph->paths.path(name))) {
            /*
               cerr << "failed inclusion" << endl
               << "expected " << graph->path_string(graph->paths.path(name)) << endl
               << "got      " << seq << endl;
               */
            failures.insert(name);
        }
    }

    if (!failures.empty()) {
        stringstream ss;
        ss << "vg-msga-failed-include_";
        for (auto& s : failures) {
            cerr << "[vg msga] Error: failed to include path " << s << endl;
            ss << s << "_";
        }
        ss << ".vg";
        graph->serialize_to_file(ss.str());
        exit(1);
    }

    // return the graph
    graph->serialize_to_ostream(std::cout);
    delete graph;

    // todo....
    //
    // strategy for graph/graph alignment
    // ---------------
    // multiple graphs can be aligned by converting them into collections of named sequences
    // e.g. using a strided sampling of a long kmer space
    // of sufficient length to align to a second graph
    //
    // the multi-graph alignment is a graph which contains both of the
    // with homologous sequences merged and the paths from the input graphs retained
    //
    // a progressive approach can be used, where we first attempt to construct a graph using a particular
    // sequence size
    // then, by labeling the first graph where it is shown to map to the new graph, we can retain only
    // the portion which was not included, then attempt to include the remaining fragments using
    // more compute-intensive parameters
    //
    // to limit path complexity, random walks should be used to sample the path space of the first graph
    // we can erode the graph we are aligning as regions of it become completely aligned,
    // so as to avoid over-sampling the graph
    // we already have functionality for this in `vg sim`
    //
    // this is an elaborate but easily-written and flexible approach to aligning large graphs
    // efficiently

    return 0;
}

void help_surject(char** argv) {
    cerr << "usage: " << argv[0] << " surject [options] <aln.gam> >[proj.cram]" << endl
        << "Transforms alignments to be relative to particular paths." << endl
        << endl
        << "options:" << endl
        << "    -d, --db-name DIR       use the graph in this database" << endl
        << "    -t, --threads N         number of threads to use" << endl
        << "    -p, --into-path NAME    surject into just this path" << endl
        << "    -i, --into-paths FILE   surject into nonoverlapping path names listed in FILE (one per line)" << endl
        << "    -P, --into-prefix NAME  surject into all paths with NAME as their prefix" << endl
        //<< "    -H, --header-from FILE  use the header in the SAM/CRAM/BAM file for the output" << endl
        // todo, reenable
        // << "    -c, --cram-output       write CRAM to stdout (default is vg::Aligment/GAM format)" << endl
        // << "    -f, --reference FILE    use this file when writing CRAM to rebuild sequence header" << endl
        << "    -b, --bam-output        write BAM to stdout" << endl
        << "    -s, --sam-output        write SAM to stdout" << endl
        << "    -C, --compression N     level for compression [0-9]" << endl
        << "    -w, --window N          use N nodes on either side of the alignment to surject (default 5)" << endl;
}

int main_surject(int argc, char** argv) {

    if (argc == 2) {
        help_surject(argv);
        return 1;
    }

    string db_name;
    string path_name;
    string path_prefix;
    string path_file;
    string output_type = "gam";
    string input_type = "gam";
    string header_file;
    int compress_level = 9;
    int window = 5;
    string fasta_filename;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"db-name", required_argument, 0, 'd'},
            {"threads", required_argument, 0, 't'},
            {"into-path", required_argument, 0, 'p'},
            {"into-paths", required_argument, 0, 'i'},
            {"into-prefix", required_argument, 0, 'P'},
            {"cram-output", no_argument, 0, 'c'},
            {"reference", required_argument, 0, 'f'},
            {"bam-output", no_argument, 0, 'b'},
            {"sam-output", no_argument, 0, 's'},
            {"header-from", required_argument, 0, 'H'},
            {"compress", required_argument, 0, 'C'},
            {"window", required_argument, 0, 'w'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hd:p:i:P:cbsH:C:t:w:f:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'd':
                db_name = optarg;
                break;

            case 'p':
                path_name = optarg;
                break;

            case 'i':
                path_file = optarg;
                break;

            case 'P':
                path_prefix = optarg;
                break;

            case 'H':
                header_file = optarg;
                break;

            case 'c':
                output_type = "cram";
                break;

            case 'f':
                fasta_filename = optarg;
                break;

            case 'b':
                output_type = "bam";
                break;

            case 's':
                compress_level = -1;
                output_type = "sam";
                break;

            case 't':
                omp_set_num_threads(atoi(optarg));
                break;

            case 'C':
                compress_level = atoi(optarg);
                break;

            case 'w':
                window = atoi(optarg);
                break;

            case 'h':
            case '?':
                help_surject(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    string file_name = argv[optind];

    Index index;
    // open index
    index.open_read_only(db_name);

    set<string> path_names;
    if (!path_file.empty()){
        // open the file
        ifstream in(path_file);
        string line;
        while (std::getline(in,line)) {
            path_names.insert(line);
        }
    } else {
        path_names.insert(path_name);
    }

    if (input_type == "gam") {
        if (output_type == "gam") {
            int thread_count = get_thread_count();
            vector<vector<Alignment> > buffer;
            buffer.resize(thread_count);
            function<void(Alignment&)> lambda = [&index, &path_names, &buffer, &window](Alignment& src) {
                int tid = omp_get_thread_num();
                Alignment surj;
                string path_name;
                int64_t path_pos;
                index.surject_alignment(src, path_names, surj, path_name, path_pos, window);
                buffer[tid].push_back(surj);
                stream::write_buffered(cout, buffer[tid], 1000);
            };
            if (file_name == "-") {
                stream::for_each_parallel(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each_parallel(in, lambda);
            }
            for (int i = 0; i < thread_count; ++i) {
                stream::write_buffered(cout, buffer[i], 0); // flush
            }
        } else {
            char out_mode[5];
            string out_format = "";
            strcpy(out_mode, "w");
            if (output_type == "bam") { out_format = "b"; }
            else if (output_type == "cram") { out_format = "c"; }
            else { out_format = ""; }
            strcat(out_mode, out_format.c_str());
            if (compress_level >= 0) {
                char tmp[2];
                tmp[0] = compress_level + '0'; tmp[1] = '\0';
                strcat(out_mode, tmp);
            }
            // get the header
            /*
               if (header_file.empty()) {
               cerr << "[vg surject] error: --header-from must be specified for SAM/BAM/CRAM output" << endl;
               return 1;
               }
               */
            string header;
            map<string, int64_t> path_by_id = index.paths_by_id();
            map<string, pair<pair<int64_t, bool>, pair<int64_t, bool>>> path_layout;
            map<string, int64_t> path_length;
            index.path_layout(path_layout, path_length);
            int thread_count = get_thread_count();
            vector<vector<tuple<string, int64_t, Alignment> > > buffer;
            buffer.resize(thread_count);
            map<string, string> rg_sample;

            // bam/sam/cram output
            samFile* out = 0;
            int buffer_limit = 100;

            bam_hdr_t* hdr = NULL;
            int64_t count = 0;
            omp_lock_t output_lock;
            omp_init_lock(&output_lock);

            // handles buffers, possibly opening the output file if we're on the first record
            auto handle_buffer =
                [&hdr, &header, &path_length, &rg_sample, &buffer_limit,
                &out_mode, &out, &output_lock, &fasta_filename](vector<tuple<string, int64_t, Alignment> >& buf) {
                    if (buf.size() >= buffer_limit) {
                        // do we have enough data to open the file?
#pragma omp critical (hts_header)
                        {
                            if (!hdr) {
                                hdr = hts_string_header(header, path_length, rg_sample);
                                if ((out = sam_open("-", out_mode)) == 0) {
                                    /*
                                       if (!fasta_filename.empty()) {
                                       string fai_filename = fasta_filename + ".fai";
                                       hts_set_fai_filename(out, fai_filename.c_str());
                                       }
                                       */
                                    cerr << "[vg surject] failed to open stdout for writing HTS output" << endl;
                                    exit(1);
                                } else {
                                    // write the header
                                    if (sam_hdr_write(out, hdr) != 0) {
                                        cerr << "[vg surject] error: failed to write the SAM header" << endl;
                                    }
                                }
                            }
                        }
                        // try to get a lock, and force things if we've built up a huge buffer waiting
                        if (omp_test_lock(&output_lock) || buf.size() > 10*buffer_limit) {
                            for (auto& s : buf) {
                                auto& path_nom = get<0>(s);
                                auto& path_pos = get<1>(s);
                                auto& surj = get<2>(s);
                                string cigar = cigar_against_path(surj);
                                bam1_t* b = alignment_to_bam(header,
                                        surj,
                                        path_nom,
                                        path_pos,
                                        cigar,
                                        "=",
                                        path_pos,
                                        0);
                                int r = 0;
#pragma omp critical (cout)
                                r = sam_write1(out, hdr, b);
                                if (r == 0) { cerr << "[vg surject] error: writing to stdout failed" << endl; exit(1); }
                                bam_destroy1(b);
                            }
                            omp_unset_lock(&output_lock);
                            buf.clear();
                        }
                    }
                };

            function<void(Alignment&)> lambda = [&index,
                &path_names,
                &path_length,
                &window,
                &rg_sample,
                &header,
                &out,
                &buffer,
                &count,
                &hdr,
                &out_mode,
                &handle_buffer](Alignment& src) {
                    int tid = omp_get_thread_num();
                    Alignment surj;
                    string path_name;
                    int64_t path_pos;
                    index.surject_alignment(src, path_names, surj, path_name, path_pos, window);
                    if (!surj.path().mapping_size()) {
                        surj = src;
                    }
                    // record
                    if (!hdr && !surj.read_group().empty() && !surj.sample_name().empty()) {
#pragma omp critical (hts_header)
                        rg_sample[surj.read_group()] = surj.sample_name();
                    }

                    buffer[tid].push_back(make_tuple(path_name, path_pos, surj));
                    handle_buffer(buffer[tid]);

                };

            // now apply the alignment processor to the stream
            if (file_name == "-") {
                stream::for_each_parallel(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each_parallel(in, lambda);
            }
            buffer_limit = 0;
            for (auto& buf : buffer) {
                handle_buffer(buf);
            }
            bam_hdr_destroy(hdr);
            sam_close(out);
            omp_destroy_lock(&output_lock);
        }
    }
    cout.flush();

    return 0;
}

void help_mod(char** argv) {
    cerr << "usage: " << argv[0] << " mod [options] <graph.vg> >[mod.vg]" << endl
         << "Modifies graph, outputs modified on stdout." << endl
         << endl
         << "options:" << endl
         << "    -i, --include-aln FILE  merge the paths implied by alignments into the graph" << endl
         << "    -P, --label-paths       don't edit with -i alignments, just use them for labeling the graph" << endl
         << "    -c, --compact-ids       should we sort and compact the id space? (default false)" << endl
         << "    -C, --compact-ranks     compact mapping ranks in paths" << endl
         << "    -z, --sort              sort the graph using an approximate topological sort" << endl
         << "    -b, --break-cycles      use an approximate topological sort to break cycles in the graph" << endl
         << "    -n, --normalize         normalize the graph so that edges are always non-redundant" << endl
         << "                            (nodes have unique starting and ending bases relative to neighbors," << endl
         << "                            and edges that do not introduce new paths are removed and neighboring" << endl
         << "                            nodes are merged)" << endl
         << "    -s, --simplify          remove redundancy from the graph that will not change its path space" << endl
         << "    -T, --strong-connect    outputs the strongly-connected components of the graph" << endl
         << "    -d, --dagify-step N     copy strongly connected components of the graph N times, forwarding" << endl
         << "                            edges from old to new copies to convert the graph into a DAG" << endl
         << "    -w, --dagify-to N       copy strongly connected components of the graph forwarding" << endl
         << "                            edges from old to new copies to convert the graph into a DAG" << endl
         << "                            until the shortest path through each SCC is N bases long" << endl
         << "    -L, --dagify-len-max N  stop a dagification step if the unrolling component has this much sequence" << endl
         << "    -U, --unroll N          using backtracking to unroll cycles in the graph, preserving paths of length N" << endl
         << "    -B, --max-branch N      maximum number of branchings to consider when unrolling" << endl
         << "    -f, --unfold N          represent inversions accesible up to N from the forward" << endl
         << "                            component of the graph" << endl
         << "    -O, --orient-forward    orient the nodes in the graph forward" << endl
         << "    -D, --drop-paths        remove the paths of the graph" << endl
         << "    -r, --retain-path NAME  remove any path not specified for retention" << endl
         << "    -k, --keep-path NAME    keep only nodes and edges in the path" << endl
         << "    -N, --remove-non-path   keep only nodes and edges which are part of paths" << endl
         << "    -o, --remove-orphans    remove orphan edges from graph (edge specified but node missing)" << endl
         << "    -R, --remove-null       removes nodes that have no sequence, forwarding their edges" << endl
         << "    -g, --subgraph ID       gets the subgraph rooted at node ID, multiple allowed" << endl
         << "    -x, --context N         steps the subgraph out by N steps (default: 1)" << endl
         << "    -p, --prune-complex     remove nodes that are reached by paths of --path-length which" << endl
         << "                            cross more than --edge-max edges" << endl
         << "    -S, --prune-subgraphs   remove subgraphs which are shorter than --length" << endl
         << "    -l, --length N          for pruning complex regions and short subgraphs" << endl
         << "    -X, --chop N            chop nodes in the graph so they are not more than N bp long" << endl
         << "    -u, --unchop            where two nodes are only connected to each other and by one edge" << endl
         << "                            replace the pair with a single node that is the concatenation of their labels" << endl
         << "    -K, --kill-labels       delete the labels from the graph, resulting in empty nodes" << endl
         << "    -e, --edge-max N        only consider paths which make edge choices at <= this many points" << endl
         << "    -m, --markers           join all head and tails nodes to marker nodes" << endl
         << "                            ('###' starts and '$$$' ends) of --path-length, for debugging" << endl
         << "    -F, --force-path-match  sets path edits explicitly equal to the nodes they traverse" << endl
         << "    -t, --threads N         for tasks that can be done in parallel, use this many threads" << endl;
}

int main_mod(int argc, char** argv) {

    if (argc == 2) {
        help_mod(argv);
        return 1;
    }

    string path_name;
    bool remove_orphans = false;
    string aln_file;
    bool label_paths = false;
    bool compact_ids = false;
    bool prune_complex = false;
    int path_length = 0;
    int edge_max = 0;
    int chop_to = 0;
    bool add_start_and_end_markers = false;
    bool prune_subgraphs = false;
    bool kill_labels = false;
    bool simplify_graph = false;
    bool unchop = false;
    bool normalize_graph = false;
    bool sort_graph = false;
    bool remove_non_path = false;
    bool compact_ranks = false;
    bool drop_paths = false;
    bool force_path_match = false;
    set<string> paths_to_retain;
    vector<int64_t> root_nodes;
    int32_t context_steps;
    bool remove_null;
    bool strong_connect = false;
    uint32_t unroll_to = 0;
    uint32_t unfold_to = 0;
    uint32_t unroll_max_branch = 0;
    bool break_cycles = false;
    uint32_t dagify_steps = 0;
    uint32_t dagify_to = 0;
    uint32_t dagify_component_length_max = 0;
    bool orient_forward = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                {"help", no_argument, 0, 'h'},
                {"include-aln", required_argument, 0, 'i'},
                {"compact-ids", no_argument, 0, 'c'},
                {"compact-ranks", no_argument, 0, 'C'},
                {"drop-paths", no_argument, 0, 'D'},
                {"keep-path", required_argument, 0, 'k'},
                {"remove-orphans", no_argument, 0, 'o'},
                {"prune-complex", no_argument, 0, 'p'},
                {"prune-subgraphs", no_argument, 0, 'S'},
                {"length", required_argument, 0, 'l'},
                {"edge-max", required_argument, 0, 'e'},
                {"chop", required_argument, 0, 'X'},
                {"kill-labels", no_argument, 0, 'K'},
                {"markers", no_argument, 0, 'm'},
                {"threads", no_argument, 0, 't'},
                {"label-paths", no_argument, 0, 'P'},
                {"simplify", no_argument, 0, 's'},
                {"unchop", no_argument, 0, 'u'},
                {"normalize", no_argument, 0, 'n'},
                {"sort", no_argument, 0, 'z'},
                {"remove-non-path", no_argument, 0, 'N'},
                {"orient-forward", no_argument, 0, 'O'},
                {"unfold", required_argument, 0, 'f'},
                {"force-path-match", no_argument, 0, 'F'},
                {"retain-path", required_argument, 0, 'r'},
                {"subgraph", required_argument, 0, 'g'},
                {"context", required_argument, 0, 'x'},
                {"remove-null", no_argument, 0, 'R'},
                {"strong-connect", no_argument, 0, 'T'},
                {"dagify-steps", required_argument, 0, 'd'},
                {"dagify-to", required_argument, 0, 'w'},
                {"dagify-len-max", required_argument, 0, 'L'},
                {"unroll", required_argument, 0, 'U'},
                {"max-branch", required_argument, 0, 'B'},
                {"break-cycles", no_argument, 0, 'b'},
                {"orient-forward", no_argument, 0, 'O'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:oi:cpl:e:mt:SX:KPsunzNf:CDFr:g:x:RTU:B:bd:Ow:L:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'i':
                aln_file = optarg;
                break;

            case 'c':
                compact_ids = true;
                break;

            case 'C':
                compact_ranks = true;
                break;

            case 'k':
                path_name = optarg;
                break;

            case 'r':
                paths_to_retain.insert(optarg);
                break;

            case 'o':
                remove_orphans = true;
                break;

            case 'p':
                prune_complex = true;
                break;

            case 'S':
                prune_subgraphs = true;
                break;

            case 'l':
                path_length = atoi(optarg);
                break;

            case 'X':
                chop_to = atoi(optarg);
                break;

            case 'u':
                unchop = true;
                break;

            case 'K':
                kill_labels = true;
                break;

            case 'e':
                edge_max = atoi(optarg);
                break;

            case 'm':
                add_start_and_end_markers = true;
                break;

            case 't':
                omp_set_num_threads(atoi(optarg));
                break;

            case 'f':
                unfold_to = atoi(optarg);
                break;

            case 'O':
                orient_forward = true;
                break;

            case 'F':
                force_path_match = true;
                break;

            case 'P':
                label_paths = true;
                break;

            case 'D':
                drop_paths = true;
                break;

            case 's':
                simplify_graph = true;
                break;

            case 'n':
                normalize_graph = true;
                break;

            case 'N':
                remove_non_path = true;
                break;

            case 'T':
                strong_connect = true;
                break;

            case 'U':
                unroll_to = atoi(optarg);
                break;

            case 'd':
                dagify_steps = atoi(optarg);
                break;

            case 'w':
                dagify_to = atoi(optarg);
                break;


        case 'L':
            dagify_component_length_max = atoi(optarg);
            break;

        case 'B':
            unroll_max_branch = atoi(optarg);
            break;


            case 'z':
                sort_graph = true;
                break;

            case 'b':
                break_cycles = true;
                break;

            case 'g':
                root_nodes.push_back(atoi(optarg));
                break;

            case 'x':
                context_steps = atoi(optarg);
                break;

            case 'R':
                remove_null = true;
                break;

            case 'h':
            case '?':
                help_mod(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    if (!path_name.empty()) {
        graph->keep_path(path_name);
    }

    if (!paths_to_retain.empty()) {
        graph->paths.keep_paths(paths_to_retain);
    }

    if (drop_paths) {
        graph->paths.clear();
    }

    if (remove_orphans) {
        graph->remove_orphan_edges();
    }

    if (unchop) {
        graph->unchop();
    }

    if (simplify_graph) {
        graph->simplify_siblings();
    }

    if (normalize_graph) {
        graph->normalize();
    }

    if (strong_connect) {
        graph->keep_multinode_strongly_connected_components();
    }

    if (remove_non_path) {
        graph->remove_non_path();
    }

    if (force_path_match) {
        graph->force_path_match();
    }

    if (orient_forward) {
        set<int64_t> flipped;
        graph->orient_nodes_forward(flipped);
    }

    if (dagify_steps) {
        map<int64_t, pair<int64_t, bool> > node_translation;
        *graph = graph->dagify(dagify_steps, node_translation, 0, dagify_component_length_max);
    }

    if (dagify_to) {
        map<int64_t, pair<int64_t, bool> > node_translation;
        // use the walk as our maximum number of steps; it's the worst case
        *graph = graph->dagify(dagify_to, node_translation, dagify_to, dagify_component_length_max);
    }

    if (unroll_to) {
        map<int64_t, pair<int64_t, bool> > node_translation;
        *graph = graph->backtracking_unroll(unroll_to, unroll_max_branch, node_translation);
    }

    if (unfold_to) {
        map<int64_t, pair<int64_t, bool> > node_translation;
        *graph = graph->unfold(unfold_to, node_translation);
    }

    if (remove_null) {
        graph->remove_null_nodes_forwarding_edges();
    }

    if (sort_graph) {
        graph->sort();
    }

    if (break_cycles) {
        graph->break_cycles();
    }

    // to subset the graph
    if (!root_nodes.empty()) {
        VG g;
        for (auto root : root_nodes) {
            graph->nonoverlapping_node_context_without_paths(graph->get_node(root), g);
            graph->expand_context(g, max(context_steps, 1));
            g.remove_orphan_edges();
        }
        *graph = g;
    }

    if (!aln_file.empty()) {
        // read in the alignments and save their paths
        vector<Path> paths;
        function<void(Alignment&)> lambda = [&graph, &paths](Alignment& aln) {
            Path path = simplify(aln.path());
            path.set_name(aln.name());
            paths.push_back(path);
        };
        if (aln_file == "-") {
            stream::for_each(std::cin, lambda);
        } else {
            ifstream in;
            in.open(aln_file.c_str());
            stream::for_each(in, lambda);
        }
        if (!label_paths) {
            // execute the edits
            graph->edit(paths);
        } else {
            // just add the path labels to the graph
            for (auto& path : paths) {
                graph->paths.extend(path);
            }
        }
    }

    // and optionally compact ids
    if (compact_ids) {
        graph->sort();
        graph->compact_ids();
    }

    if (compact_ranks) {
        graph->paths.compact_ranks();
    }

    if (prune_complex) {
        if (!(path_length > 0 && edge_max > 0)) {
            cerr << "[vg mod]: when pruning complex regions you must specify a --path-length and --edge-max" << endl;
            return 1;
        }
        graph->prune_complex_with_head_tail(path_length, edge_max);
    }

    if (prune_subgraphs) {
        graph->prune_short_subgraphs(path_length);
    }

    if (chop_to) {
        graph->dice_nodes(chop_to);
        graph->paths.compact_ranks();
    }

    if (kill_labels) {
        graph->for_each_node([](Node* n) { n->clear_sequence(); });
    }

    if (add_start_and_end_markers) {
        if (!(path_length > 0)) {
            cerr << "[vg mod]: when adding start and end markers you must provide a --path-length" << endl;
            return 1;
        }
        Node* head_node = NULL;
        Node* tail_node = NULL;
        graph->add_start_end_markers(path_length, '#', '$', head_node, tail_node);
    }

    graph->serialize_to_ostream(std::cout);

    delete graph;

    return 0;
}

void help_sim(char** argv) {
    cerr << "usage: " << argv[0] << " sim [options] <graph.vg>" << endl
        << "Simulates reads from the graph(s). Output is a list of reads." << endl
        << endl
        << "options:" << endl
        << "    -l, --read-length N   write reads of length N" << endl
        << "    -n, --num-reads N     simulate N reads" << endl
        << "    -s, --random-seed N   use this specific seed for the PRNG" << endl
        << "    -e, --base-error N    base substitution error rate (default 0.0)" << endl
        << "    -i, --indel-error N   indel error rate (default 0.0)" << endl
        << "    -f, --forward-only    don't simulate from the reverse strand" << endl
        << "    -a, --align-out       generate true alignments on stdout rather than reads" << endl;
}

int main_sim(int argc, char** argv) {

    if (argc == 2) {
        help_sim(argv);
        return 1;
    }

    int read_length = 100;
    int num_reads = 1;
    int seed_val = time(NULL);
    double base_error = 0;
    double indel_error = 0;
    bool forward_only = false;
    bool align_out = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"read-length", required_argument, 0, 'l'},
            {"num-reads", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"forward-only", no_argument, 0, 'f'},
            {"align-out", no_argument, 0, 'a'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:e:i:fa",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'l':
                read_length = atoi(optarg);
                break;

            case 'n':
                num_reads = atoi(optarg);
                break;

            case 's':
                seed_val = atoi(optarg);
                break;

            case 'e':
                base_error = atof(optarg);
                break;

            case 'i':
                indel_error = atof(optarg);
                break;

            case 'f':
                forward_only = true;
                break;

            case 'a':
                align_out = true;
                break;

            case 'h':
            case '?':
                help_sim(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    int64_t max_id = graph->max_node_id();
    int64_t min_id = graph->min_node_id();

    mt19937 rng;
    rng.seed(seed_val);

    string bases = "ATGC";
    uniform_real_distribution<double> rprob(0, 1);
    uniform_int_distribution<int> rbase(0, 3);

    function<string(const string&)> introduce_read_errors
        = [&rng, &rprob, &rbase, &bases, base_error, indel_error](const string& perfect_read) {

            if (base_error == 0 && indel_error == 0) return perfect_read;
            string read;
            for (auto c : perfect_read) {
                if (rprob(rng) <= base_error) {
                    // pick another base than what c is
                    char e;
                    do {
                        e = bases[rbase(rng)];
                    } while (e == c);
                    c = e;
                }
                if (rprob(rng) <= indel_error) {
                    if (rprob(rng) < 0.5) {
                        read.push_back(bases[rbase(rng)]);
                    } // else do nothing, == deletion of base
                } else {
                    read.push_back(c);
                }
            }
            return read;
        };

    size_t max_iter = 1000;
    for (int i = 0; i < num_reads; ++i) {
        auto perfect_read = graph->random_read(read_length, rng, min_id, max_id, !forward_only);
        // avoid short reads at the end of the graph by retrying
        int iter = 0;
        while (perfect_read.sequence().size() < read_length && ++iter < max_iter) {
            perfect_read = graph->random_read(read_length, rng, min_id, max_id, !forward_only);
            // if we can't make a suitable read in 1000 tries, then maybe the graph is too small?
        }
        if (iter == max_iter) {
            cerr << "couldn't simulate read, perhaps the chosen length is too long for this graph?" << endl;
        } else {
            // apply errors
            if (!align_out) {
                string readseq = introduce_read_errors(perfect_read.sequence());
                cout << readseq << endl;
            } else {
                function<Alignment(uint64_t)> lambda =
                    [&perfect_read] (uint64_t n) {
                        return perfect_read;
                    };
                stream::write(cout, 1, lambda);
            }
        }
    }
    delete graph;

    return 0;
}

void help_kmers(char** argv) {
    cerr << "usage: " << argv[0] << " kmers [options] <graph1.vg> [graph2.vg ...] >kmers.tsv" << endl

         << "Generates kmers of the graph(s). Output is: kmer id pos" << endl
         << endl
         << "options:" << endl
         << "    -k, --kmer-size N     print kmers of size N in the graph" << endl
         << "    -e, --edge-max N      only consider paths which make edge choices at <= this many points" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers in paths (default 1)" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -d, --ignore-dups     filter out duplicated kmers in normal output" << endl
         << "    -n, --allow-negs      don't filter out relative negative positions of kmers in normal output" << endl
         << "    -g, --gcsa-out        output a table suitable for input to GCSA2:" << endl
         << "                          kmer, starting position, previous characters," << endl
         << "                          successive characters, successive positions." << endl
         << "                          Forward and reverse strand kmers are reported." << endl
         << "    -B, --gcsa-binary     Write the GCSA graph in binary format." << endl
         << "    -F, --forward-only    When producing GCSA2 output, don't describe the reverse strand" << endl
         << "    -P, --path-only       Only consider kmers if they occur in a path embedded in the graph" << endl
         << "    -H, --head-id N       use the specified ID for the GCSA2 head sentinel node" << endl
         << "    -T, --tail-id N       use the specified ID for the GCSA2 tail sentinel node" << endl
         << "    -p, --progress        show progress" << endl;
}

int main_kmers(int argc, char** argv) {

    if (argc == 2) {
        help_kmers(argv);
        return 1;
    }

    int kmer_size = 0;
    bool path_only = false;
    int edge_max = 0;
    int kmer_stride = 1;
    bool show_progress = false;
    bool gcsa_out = false;
    bool allow_dups = true;
    bool allow_negs = false;
    // for distributed GCSA2 kmer generation
    int64_t head_id = 0;
    int64_t tail_id = 0;
    bool forward_only = false;
    bool gcsa_binary = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                {"help", no_argument, 0, 'h'},
                {"kmer-size", required_argument, 0, 'k'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"edge-max", required_argument, 0, 'e'},
                {"threads", required_argument, 0, 't'},
                {"gcsa-out", no_argument, 0, 'g'},
                {"ignore-dups", no_argument, 0, 'd'},
                {"allow-negs", no_argument, 0, 'n'},
                {"progress",  no_argument, 0, 'p'},
                {"head-id", required_argument, 0, 'H'},
                {"tail-id", required_argument, 0, 'T'},
                {"forward-only", no_argument, 0, 'F'},
                {"gcsa-binary", no_argument, 0, 'B'},
                {"path-only", no_argument, 0, 'P'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:j:pt:e:gdnH:T:FBP",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'k':
                kmer_size = atoi(optarg);
                break;

            case 'j':
                kmer_stride = atoi(optarg);
                break;

            case 'e':
                edge_max = atoi(optarg);
                break;

            case 't':
                omp_set_num_threads(atoi(optarg));
                break;

            case 'g':
                gcsa_out = true;
                break;

            case 'F':
                forward_only = true;
                break;


        case 'P':
            path_only = true;
            break;

        case 'd':
            allow_dups = false;
            break;

            case 'n':
                allow_negs = true;
                break;

            case 'p':
                show_progress = true;
                break;

            case 'H':
                head_id = atoi(optarg);
                break;

            case 'T':
                tail_id = atoi(optarg);
                break;

            case 'B':
                gcsa_binary = true;
                break;

            case 'h':
            case '?':
                help_kmers(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    vector<string> graph_file_names;
    while (optind < argc) {
        string file_name = argv[optind++];
        graph_file_names.push_back(file_name);
    }

    VGset graphs(graph_file_names);

    graphs.show_progress = show_progress;

    if (gcsa_out) {
        if (!gcsa_binary) {
            graphs.write_gcsa_out(cout, kmer_size, path_only, forward_only, head_id, tail_id);
        } else {
            graphs.write_gcsa_kmers_binary(cout, kmer_size, path_only, forward_only, head_id, tail_id);
        }
    } else {
        function<void(string&, list<NodeTraversal>::iterator, int, list<NodeTraversal>&, VG& graph)>
            lambda = [](string& kmer, list<NodeTraversal>::iterator n, int p, list<NodeTraversal>& path, VG& graph) {
                // We encode orientation by negating the IDs for backward nodes.
                // Their offsets are from the end of the node in its local forward
                // orientation, and are negated in the output.
                int sign = (*n).backward ? -1 : 1;
#pragma omp critical (cout)

            cout << kmer << '\t' << (*n).node->id() * sign << '\t' << p * sign << '\n';
        };
        graphs.for_each_kmer_parallel(lambda, kmer_size, path_only, edge_max, kmer_stride, allow_dups, allow_negs);
    }
    cout.flush();

    return 0;
}

void help_concat(char** argv) {
    cerr << "usage: " << argv[0] << " concat [options] <graph1.vg> [graph2.vg ...] >merged.vg" << endl
        << "Concatenates graphs in order by adding edges from the tail nodes of the" << endl
        << "predecessor to the head nodes of the following graph. Node IDs are" << endl
        << "compacted, so care should be taken if consistent IDs are required." << endl;
}

int main_concat(int argc, char** argv) {

    if (argc == 2) {
        help_concat(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
            case '?':
                help_concat(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    list<VG*> graphs;

    while (optind < argc) {
        VG* graph;
        string file_name = argv[optind++];
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
        graphs.push_back(graph);
    }

    VG merged;
    for (list<VG*>::iterator g = graphs.begin(); g != graphs.end(); ++g) {
        merged.append(**g);
    }

    // output
    merged.serialize_to_ostream(std::cout);

    return 0;
}

void help_ids(char** argv) {
    cerr << "usage: " << argv[0] << " ids [options] <graph1.vg> [graph2.vg ...] >new.vg" << endl
        << "options:" << endl
        << "    -c, --compact        minimize the space of integers used by the ids" << endl
        << "    -i, --increment N    increase ids by N" << endl
        << "    -d, --decrement N    decrease ids by N" << endl
        << "    -j, --join           make a joint id space for all the graphs that are supplied" << endl
        << "                         by iterating through the supplied graphs and incrementing" << endl
        << "                         their ids to be non-conflicting" << endl
        << "    -s, --sort           assign new node IDs in (generalized) topological sort order" << endl;
}

int main_ids(int argc, char** argv) {

    if (argc == 2) {
        help_ids(argv);
        return 1;
    }

    bool join = false;
    bool compact = false;
    bool sort = false;
    int64_t increment = 0;
    int64_t decrement = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"compact", no_argument, 0, 'c'},
            {"increment", required_argument, 0, 'i'},
            {"decrement", required_argument, 0, 'd'},
            {"join", no_argument, 0, 'j'},
            {"sort", no_argument, 0, 's'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hci:d:js",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'c':
                compact = true;
                break;

            case 'i':
                increment = atoi(optarg);
                break;

            case 'd':
                decrement = atoi(optarg);
                break;

            case 'j':
                join = true;
                break;

            case 's':
                sort = true;
                break;

            case 'h':
            case '?':
                help_ids(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    if (!join) {
        VG* graph;
        string file_name = argv[optind];
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }

        if (sort) {
            // Set up the nodes so we go through them in topological order
            graph->sort();
        }

        if (compact || sort) {
            // Compact only, or compact to re-assign IDs after sort
            graph->compact_ids();
        }

        if (increment != 0) {
            graph->increment_node_ids(increment);
        }

        if (decrement != 0) {
            graph->decrement_node_ids(decrement);
        }

        graph->serialize_to_ostream(std::cout);
        delete graph;
    } else {

        vector<string> graph_file_names;
        while (optind < argc) {
            VG* graph;
            string file_name = argv[optind++];
            graph_file_names.push_back(file_name);
        }

        VGset graphs(graph_file_names);
        graphs.merge_id_space();

    }

    return 0;

}

void help_join(char** argv) {
    cerr << "usage: " << argv[0] << " join [options] <graph1.vg> [graph2.vg ...] >joined.vg" << endl
        << "Joins graphs and sub-graphs into a single variant graph by connecting their" << endl
        << "heads to a single root node with sequence 'N'." << endl
        << "Assumes a single id namespace for all graphs to join." << endl;
}

int main_join(int argc, char** argv) {

    if (argc == 2) {
        help_join(argv);
        return 1;
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'h':
            case '?':
                help_join(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    list<VG*> graphs;

    while (optind < argc) {
        VG* graph;
        string file_name = argv[optind++];
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
        graphs.push_back(graph);
    }

    VG joined;
    for (list<VG*>::iterator g = graphs.begin(); g != graphs.end(); ++g) {
        // Stick all the graphs together, complaining if they use the same node IDs (since they probably shouldn't).
        joined.extend(**g, true);
    }

    // combine all subgraphs
    joined.join_heads();

    // output
    joined.serialize_to_ostream(std::cout);

    return 0;
}

void help_stats(char** argv) {
    cerr << "usage: " << argv[0] << " stats [options] <graph.vg>" << endl
        << "options:" << endl
        << "    -z, --size            size of graph" << endl
        << "    -N, --node-count      number of nodes in graph" << endl
        << "    -E, --edge-count      number of edges in graph" << endl
        << "    -l, --length          length of sequences in graph" << endl
        << "    -s, --subgraphs       describe subgraphs of graph" << endl
        << "    -H, --heads           list the head nodes of the graph" << endl
        << "    -T, --tails           list the tail nodes of the graph" << endl
        << "    -S, --siblings        describe the siblings of each node" << endl
        << "    -c, --components      print the strongly connected components of the graph" << endl
        << "    -n, --node ID         consider node with the given id" << endl
        << "    -d, --to-head         show distance to head for each provided node" << endl
        << "    -t, --to-tail         show distance to head for each provided node" << endl;
}

int main_stats(int argc, char** argv) {

    if (argc == 2) {
        help_stats(argv);
        return 1;
    }

    bool stats_size = false;
    bool stats_length = false;
    bool stats_subgraphs = false;
    bool stats_heads = false;
    bool stats_tails = false;
    bool show_sibs = false;
    bool show_components = false;
    bool distance_to_head = false;
    bool distance_to_tail = false;
    bool node_count = false;
    bool edge_count = false;
    set<vg::id_t> ids;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"size", no_argument, 0, 'z'},
            {"node-count", no_argument, 0, 'N'},
            {"edge-count", no_argument, 0, 'E'},
            {"length", no_argument, 0, 'l'},
            {"subgraphs", no_argument, 0, 's'},
            {"heads", no_argument, 0, 'H'},
            {"tails", no_argument, 0, 'T'},
            {"help", no_argument, 0, 'h'},
            {"siblings", no_argument, 0, 'S'},
            {"components", no_argument, 0, 'c'},
            {"to-head", no_argument, 0, 'd'},
            {"to-tail", no_argument, 0, 't'},
            {"node", required_argument, 0, 'n'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hzlsHTScdtn:NE",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'z':
                stats_size = true;
                break;

            case 'N':
                node_count = true;
                break;

            case 'E':
                edge_count = true;
                break;

            case 'l':
                stats_length = true;
                break;

            case 's':
                stats_subgraphs = true;
                break;

            case 'H':
                stats_heads = true;
                break;

            case 'T':
                stats_tails = true;
                break;

            case 'S':
                show_sibs = true;
                break;

            case 'c':
                show_components = true;
                break;

            case 'd':
                distance_to_head = true;
                break;

            case 't':
                distance_to_tail = true;
                break;

            case 'n':
                ids.insert(atoi(optarg));
                break;

            case 'h':
            case '?':
                help_stats(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    if (stats_size) {
        cout << "nodes" << "\t" << graph->node_count() << endl
            << "edges" << "\t" << graph->edge_count() << endl;
    }

    if (node_count) {
        cout << graph->node_count() << endl;
    }

    if (edge_count) {
        cout << graph->edge_count() << endl;
    }

    if (stats_length) {
        cout << "length" << "\t" << graph->total_length_of_nodes() << endl;
    }

    if (stats_heads) {
        vector<Node*> heads;
        graph->head_nodes(heads);
        cout << "heads" << "\t";
        for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
            cout << (*h)->id() << " ";
        }
        cout << endl;
    }

    if (stats_tails) {
        vector<Node*> tails;
        graph->tail_nodes(tails);
        cout << "tails" << "\t";
        for (vector<Node*>::iterator t = tails.begin(); t != tails.end(); ++t) {
            cout << (*t)->id() << " ";
        }
        cout << endl;
    }

    if (stats_subgraphs) {
        list<VG> subgraphs;
        graph->disjoint_subgraphs(subgraphs);
        // these are topologically-sorted
        for (list<VG>::iterator s = subgraphs.begin(); s != subgraphs.end(); ++s) {
            VG& subgraph = *s;
            vector<Node*> heads;
            subgraph.head_nodes(heads);
            int64_t length = subgraph.total_length_of_nodes();
            for (vector<Node*>::iterator h = heads.begin(); h != heads.end(); ++h) {
                cout << (h==heads.begin()?"":",") << (*h)->id();
            }
            cout << "\t" << length << endl;
        }
    }

    if (show_sibs) {
        graph->for_each_node([graph](Node* n) {
                for (auto trav : graph->full_siblings_to(NodeTraversal(n, false))) {
                cout << n->id() << "\t" << "to-sib" << "\t" << trav.node->id() << endl;
                }
                for (auto trav : graph->full_siblings_from(NodeTraversal(n, false))) {
                cout << n->id() << "\t" << "from-sib" << "\t" << trav.node->id() << endl;
                }
                });
    }

    if (show_components) {
        for (auto& c : graph->strongly_connected_components()) {
            for (auto& id : c) {
                cout << id << ", ";
            }
            cout << endl;
        }
    }

    if (distance_to_head) {
        for (auto id : ids) {
            cout << id << " to head:\t"
                << graph->distance_to_head(NodeTraversal(graph->get_node(id), false)) << endl;
        }
    }

    if (distance_to_tail) {
        for (auto id : ids) {
            cout << id << " to tail:\t"
                << graph->distance_to_tail(NodeTraversal(graph->get_node(id), false)) << endl;
        }
    }

    delete graph;

    return 0;

}

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options] <graph.vg>" << endl
         << "options:" << endl
         << "  obtain paths in GAM:" << endl
         << "    -x, --extract         return (as alignments) the stored paths in the graph" << endl
         << "  generation:" << endl
         << "    -n, --node ID         starting at node with ID" << endl
         << "    -l, --max-length N    generate paths of at most length N" << endl
         << "    -e, --edge-max N      only consider paths which make edge choices at this many points" << endl
         << "    -s, --as-seqs         write each path as a sequence" << endl
         << "    -p, --path-only       only write kpaths from the graph if they traverse embedded paths" << endl;
}

int main_paths(int argc, char** argv) {

    if (argc == 2) {
        help_paths(argv);
        return 1;
    }

    int max_length = 0;
    int edge_max = 0;
    int64_t node_id = 0;
    bool as_seqs = false;
    bool extract = false;
    bool path_only = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                {"extract", no_argument, 0, 'x'},
                {"node", required_argument, 0, 'n'},
                {"max-length", required_argument, 0, 'l'},
                {"edge-max", required_argument, 0, 'e'},
                {"as-seqs", no_argument, 0, 's'},
                {"path-only", no_argument, 0, 'p'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "n:l:hse:xp",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'x':
                extract = true;
                break;

            case 'n':
                node_id = atoll(optarg);
                break;

            case 'l':
                max_length = atoi(optarg);
                break;

            case 'e':
                edge_max = atoi(optarg);
                break;

            case 's':
                as_seqs = true;
                break;

        case 'p':
            path_only = true;
            break;

        case 'h':
        case '?':
            help_paths(argv);
            exit(1);
            break;

            default:
                abort ();
        }
    }

    if (edge_max == 0) edge_max = max_length + 1;

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    if (extract) {
        vector<Alignment> alns = graph->paths_as_alignments();
        write_alignments(cout, alns);
        delete graph;
        return 0;
    }

    if (max_length == 0) {
        cerr << "error:[vg paths] a --max-length is required when generating paths" << endl;
    }

    function<void(size_t,Path&)> paths_to_seqs = [graph](size_t mapping_index, Path& p) {
        string seq = graph->path_sequence(p);
#pragma omp critical(cout)
        cout << seq << endl;
    };

    function<void(size_t,Path&)> paths_to_json = [](size_t mapping_index, Path& p) {
        string json2 = pb2json(p);
#pragma omp critical(cout)
        cout<<json2<<endl;
    };

    function<void(size_t, Path&)>* callback = &paths_to_seqs;
    if (!as_seqs) {
        callback = &paths_to_json;
    }

    auto noop = [](NodeTraversal) { }; // don't handle the failed regions of the graph yet

    if (node_id) {

        graph->for_each_kpath_of_node(graph->get_node(node_id),
                                      path_only,
                                      max_length,
                                      edge_max,
                                      noop, noop,
                                      *callback);
    } else {
        graph->for_each_kpath_parallel(max_length,
                                       path_only,
                                       edge_max,
                                       noop, noop,
                                       *callback);
    }

    delete graph;

    return 0;

}

void help_find(char** argv) {
    cerr << "usage: " << argv[0] << " find [options] <graph.vg> >sub.vg" << endl
        << "options:" << endl
        << "    -d, --db-name DIR      use this db (defaults to <graph>.index/)" << endl
        // TODO, dump these from the index
        //<< "    -a, --alignments       write all stored alignments in sorted order (in GAM)" << endl
        //<< "    -m, --mappings         write stored mappings in sorted order (in json)" << endl
         << "    -x, --xg-name FILE     use this xg index (instead of rocksdb db)" << endl
         << "graph features:" << endl
         << "    -n, --node ID          find node, return 1-hop context as graph" << endl
         << "    -e, --edges-end ID     return edges on end of node with ID" << endl
         << "    -s, --edges-start ID   return edges on start of node with ID" << endl
         << "    -c, --context STEPS    expand the context of the kmer hit subgraphs" << endl
         << "    -p, --path TARGET      find the node(s) in the specified path range TARGET=path[:pos1[-pos2]]" << endl
         << "    -P, --position-in PATH find the position of the node (specified by -n) in the given path" << endl
         << "    -r, --node-range N:M   get nodes from N to M" << endl
         << "sequences:" << endl
         << "    -g, --gcsa FILE        use this GCSA2 index of the sequence space of the graph" << endl
         << "    -z, --kmer-size N      split up --sequence into kmers of size N" << endl
         << "    -j, --kmer-stride N    step distance between succesive kmers in sequence (default 1)" << endl
         << "    -S, --sequence STR     search for sequence STR using --kmer-size kmers" << endl
         << "    -M, --mems STR         describe the super-maximal exact matches of the STR (gcsa2) in JSON" << endl
         << "    -k, --kmer STR         return a graph of edges and nodes matching this kmer" << endl
         << "    -T, --table            instead of a graph, return a table of kmers" << endl
         << "                           (works only with kmers in the index)" << endl
         << "    -C, --kmer-count       report approximate count of kmer (-k) in db" << endl;


}

int main_find(int argc, char** argv) {

    if (argc == 2) {
        help_find(argv);
        return 1;
    }

    string db_name;
    string sequence;
    int kmer_size=0;
    int kmer_stride = 1;
    vector<string> kmers;
    string output_format;
    int64_t end_id=0, start_id=0;
    vector<int64_t> node_ids;
    int context_size=0;
    bool count_kmers = false;
    bool kmer_table = false;
    string target;
    string path_name;
    string range;
    bool get_alignments = false;
    bool get_mappings = false;
    string gcsa_in;
    string xg_name;
    bool get_mems = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"db-name", required_argument, 0, 'd'},
                {"xg-name", required_argument, 0, 'x'},
                {"gcsa", required_argument, 0, 'g'},
                {"node", required_argument, 0, 'n'},
                {"edges-end", required_argument, 0, 'e'},
                {"edges-start", required_argument, 0, 's'},
                {"kmer", required_argument, 0, 'k'},
                {"table", no_argument, 0, 'T'},
                {"sequence", required_argument, 0, 'S'},
                {"mems", required_argument, 0, 'M'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"kmer-size", required_argument, 0, 'z'},
                {"output", required_argument, 0, 'o'},
                {"context", required_argument, 0, 'c'},
                {"kmer-count", no_argument, 0, 'C'},
                {"path", required_argument, 0, 'p'},
                {"position-in", required_argument, 0, 'P'},
                {"node-range", required_argument, 0, 'r'},
                {"alignments", no_argument, 0, 'a'},
                {"mappings", no_argument, 0, 'm'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:x:n:e:s:o:k:hc:S:z:j:CTp:P:r:amg:M:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'd':
                db_name = optarg;
                break;

            case 'x':
                xg_name = optarg;
                break;

            case 'g':
                gcsa_in = optarg;
                break;

            case 'k':
                kmers.push_back(optarg);
                break;

            case 'S':
                sequence = optarg;
                break;

        case 'M':
            sequence = optarg;
            get_mems = true;
            break;

        case 'j':
            kmer_stride = atoi(optarg);
            break;

            case 'z':
                kmer_size = atoi(optarg);
                break;

            case 'C':
                count_kmers = true;
                break;

            case 'p':
                target = optarg;
                break;

            case 'P':
                path_name = optarg;
                break;

            case 'c':
                context_size = atoi(optarg);
                break;

            case 'n':
                node_ids.push_back(atoi(optarg));
                break;

            case 'e':
                end_id = atoi(optarg);
                break;

            case 's':
                start_id = atoi(optarg);
                break;

            case 'T':
                kmer_table = true;
                break;

            case 'r':
                range = optarg;
                break;

            case 'a':
                get_alignments = true;
                break;

            case 'm':
                get_mappings = true;
                break;

            case 'o':
                output_format = optarg;
                break;

            case 'h':
            case '?':
                help_find(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }
    if (optind < argc) {
        //string file_name = argv[optind];
        cerr << "[vg find] find requires -d, -g, or -x to know where to find its database" << endl;
        return 1;
    }

    // open index
    Index* vindex = nullptr;
    if (db_name.empty()) {
        assert(!gcsa_in.empty() || !xg_name.empty());
    } else {
        vindex = new Index;
        vindex->open_read_only(db_name);
    }

    xg::XG xindex;
    if (!xg_name.empty()) {
        ifstream in(xg_name.c_str());
        xindex.load(in);
    }

    if (get_alignments) {
        // todo
    }

    if (!xg_name.empty()) {
        if (!node_ids.empty() && path_name.empty()) {
            // get the context of the node
            vector<Graph> graphs;
            for (auto node_id : node_ids) {
                Graph g;
                xindex.neighborhood(node_id, context_size, g);
                graphs.push_back(g);
            }
            VG result_graph;
            for (auto& graph : graphs) {
                // Allow duplicate nodes and edges (from e.g. multiple -n options); silently collapse them.
                result_graph.extend(graph);
            }
            result_graph.remove_orphan_edges();
            // return it
            result_graph.serialize_to_ostream(cout);
        } else if (end_id != 0) {
            for (auto& e : xindex.edges_on_end(end_id)) {
                cout << (e.from_start() ? -1 : 1) * e.from() << "\t" <<  (e.to_end() ? -1 : 1) * e.to() << endl;
            }
        } else if (start_id != 0) {
            for (auto& e : xindex.edges_on_start(start_id)) {
                cout << (e.from_start() ? -1 : 1) * e.from() << "\t" <<  (e.to_end() ? -1 : 1) * e.to() << endl;
            }
        }
        if (!target.empty()) {
            string name;
            int64_t start, end;
            Graph graph;
            parse_region(target, name, start, end);
            xindex.get_path_range(name, start, end, graph);
            if (context_size > 0) {
                xindex.expand_context(graph, context_size);
            }
            VG vgg; vgg.extend(graph); // removes dupes
            vgg.serialize_to_ostream(cout);
        }
        if (!range.empty()) {
            Graph graph;
            int64_t id_start=0, id_end=0;
            vector<string> parts = split_delims(range, ":");
            if (parts.size() == 1) {
                cerr << "[vg find] error, format of range must be \"N:M\" where start id is N and end id is M, got " << range << endl;
                exit(1);
            }
            convert(parts.front(), id_start);
            convert(parts.back(), id_end);
            xindex.get_id_range(id_start, id_end, graph);
            if (context_size > 0) {
                xindex.expand_context(graph, context_size);
            }
            VG vgg; vgg.extend(graph); // removes dupes
            vgg.remove_orphan_edges();
            vgg.serialize_to_ostream(cout);
        }
    } else if (!db_name.empty()) {
        if (!node_ids.empty() && path_name.empty()) {
            // get the context of the node
            vector<VG> graphs;
            for (auto node_id : node_ids) {
                VG g;
                vindex->get_context(node_id, g);
                if (context_size > 0) {
                    vindex->expand_context(g, context_size);
                }
                graphs.push_back(g);
            }
            VG result_graph;
            for (auto& graph : graphs) {
                // Allow duplicate nodes and edges (from e.g. multiple -n options); silently collapse them.
                result_graph.extend(graph);
            }
            result_graph.remove_orphan_edges();
            // return it
            result_graph.serialize_to_ostream(cout);
        } else if (end_id != 0) {
            vector<Edge> edges;
            vindex->get_edges_on_end(end_id, edges);
            for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
                cout << (e->from_start() ? -1 : 1) * e->from() << "\t" <<  (e->to_end() ? -1 : 1) * e->to() << endl;
            }
        } else if (start_id != 0) {
            vector<Edge> edges;
            vindex->get_edges_on_start(start_id, edges);
            for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
                cout << (e->from_start() ? -1 : 1) * e->from() << "\t" <<  (e->to_end() ? -1 : 1) * e->to() << endl;
            }
        }
        if (!node_ids.empty() && !path_name.empty()) {
            int64_t path_id = vindex->get_path_id(path_name);
            for (auto node_id : node_ids) {
                list<pair<int64_t, bool>> path_prev, path_next;
                int64_t prev_pos=0, next_pos=0;
                bool prev_backward, next_backward;
                if (vindex->get_node_path_relative_position(node_id, false, path_id,
                                                          path_prev, prev_pos, prev_backward,
                                                          path_next, next_pos, next_backward)) {

                    // Negate IDs for backward nodes
                    cout << node_id << "\t" << path_prev.front().first * (path_prev.front().second ? -1 : 1) << "\t" << prev_pos
                        << "\t" << path_next.back().first * (path_next.back().second ? -1 : 1) << "\t" << next_pos << "\t";

                    Mapping m = vindex->path_relative_mapping(node_id, false, path_id,
                                                            path_prev, prev_pos, prev_backward,
                                                            path_next, next_pos, next_backward);
                    cout << pb2json(m) << endl;
                }
            }
        }
        if (!target.empty()) {
            string name;
            int64_t start, end;
            VG graph;
            parse_region(target, name, start, end);
            vindex->get_path(graph, name, start, end);
            if (context_size > 0) {
                vindex->expand_context(graph, context_size);
            }
            graph.remove_orphan_edges();
            graph.serialize_to_ostream(cout);
        }
        if (!range.empty()) {
            VG graph;
            int64_t id_start=0, id_end=0;
            vector<string> parts = split_delims(range, ":");
            if (parts.size() == 1) {
                cerr << "[vg find] error, format of range must be \"N:M\" where start id is N and end id is M, got " << range << endl;
                exit(1);
            }
            convert(parts.front(), id_start);
            convert(parts.back(), id_end);
            vindex->get_range(id_start, id_end, graph);
            if (context_size > 0) {
                vindex->expand_context(graph, context_size);
            }
            graph.remove_orphan_edges();
            graph.serialize_to_ostream(cout);
        }
    }

    // todo cleanup if/else logic to allow only one function

    if (!sequence.empty()) {
        if (gcsa_in.empty()) {
            if (get_mems) {
                cerr << "error:[vg find] a GCSA index must be passed to get MEMs" << endl;
                return 1;
            }
            set<int> kmer_sizes = vindex->stored_kmer_sizes();
            if (kmer_sizes.empty()) {
                cerr << "error:[vg find] index does not include kmers, add with vg index -k" << endl;
                return 1;
            }
            if (kmer_size == 0) {
                kmer_size = *kmer_sizes.begin();
            }
            for (int i = 0; i <= sequence.size()-kmer_size; i+=kmer_stride) {
                kmers.push_back(sequence.substr(i,kmer_size));
            }
        } else {
            // let's use the GCSA index
            // first open it
            ifstream in_gcsa(gcsa_in.c_str());
            gcsa::GCSA gcsa_index;
            gcsa_index.load(in_gcsa);
            gcsa::LCPArray lcp_index;
            // default LCP is the gcsa base name +.lcp
            string lcp_in = gcsa_in + ".lcp";
            ifstream in_lcp(lcp_in.c_str());
            lcp_index.load(in_lcp);
            //range_type find(const char* pattern, size_type length) const;
            //void locate(size_type path, std::vector<node_type>& results, bool append = false, bool sort = true) const;
            //locate(i, results);
            if (!get_mems) {
                auto paths = gcsa_index.find(sequence.c_str(), sequence.length());
                //cerr << paths.first << " - " << paths.second << endl;
                for (gcsa::size_type i = paths.first; i <= paths.second; ++i) {
                    std::vector<gcsa::node_type> ids;
                    gcsa_index.locate(i, ids);
                    for (auto id : ids) {
                        cout << gcsa::Node::decode(id) << endl;
                    }
                }
            } else {
                // for mems we need to load up the gcsa and lcp structures into the mapper
                Mapper mapper;
                mapper.gcsa = &gcsa_index;
                mapper.lcp = &lcp_index;
                // get the mems
                auto mems = mapper.find_smems(sequence);
                // then fill the nodes that they match
                for (auto& mem : mems) mem.fill_nodes(&gcsa_index);
                // dump them to stdout
                cout << mems_to_json(mems) << endl;
            }
        }
    }

    if (!kmers.empty()) {
        if (count_kmers) {
            for (auto& kmer : kmers) {
                cout << kmer << "\t" << vindex->approx_size_of_kmer_matches(kmer) << endl;
            }
        } else if (kmer_table) {
            for (auto& kmer : kmers) {
                map<string, vector<pair<int64_t, int32_t> > > positions;
                vindex->get_kmer_positions(kmer, positions);
                for (auto& k : positions) {
                    for (auto& p : k.second) {
                        cout << k.first << "\t" << p.first << "\t" << p.second << endl;
                    }
                }
            }
        } else {
            vector<VG> graphs;
            for (auto& kmer : kmers) {
                VG g;
                vindex->get_kmer_subgraph(kmer, g);
                if (context_size > 0) {
                    vindex->expand_context(g, context_size);
                }
                graphs.push_back(g);
            }

            VG result_graph;
            for (auto& graph : graphs) {
                // Allow duplicate nodes and edges (from multiple kmers); silently collapse them.
                result_graph.extend(graph);
            }
            result_graph.remove_orphan_edges();
            result_graph.serialize_to_ostream(cout);
        }
    }

    if (vindex) delete vindex;

    return 0;

}

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space, and the graph containing the highest-ID node " << endl
         << "must come first." << endl
         << "xg options:" << endl
         << "    -x, --xg-name FILE     use this file to store a succinct, queryable version of" << endl
         << "                           the graph(s) (effectively replaces rocksdb)" << endl
         << "    -v, --vcf-phasing FILE import phasing blocks from the given VCF file as threads" << endl
         << "    -T, --store-threads    use gPBWT to store the embedded paths as threads" << endl
         << "gcsa options:" << endl
         << "    -g, --gcsa-out FILE    output a GCSA2 index instead of a rocksdb index" << endl
         << "    -i, --dbg-in FILE      optionally use deBruijn graph encoded in FILE rather than an input VG (multiple allowed" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph" << endl
         << "    -X, --doubling-steps N use this number of doubling steps for GCSA2 construction" << endl
         << "    -Z, --size-limit N     limit of memory to use for GCSA2 construction in gigabytes" << endl
         << "    -O, --path-only        only index the kmers in paths embedded in the graph" << endl
         << "    -F, --forward-only     omit the reverse complement of the graph from indexing" << endl
         << "    -e, --edge-max N       only consider paths which make edge choices at <= this many points" << endl
         << "    -j, --kmer-stride N    step distance between succesive kmers in paths (default 1)" << endl
         << "    -d, --db-name PATH     create rocksdb in PATH directory (default: <graph>.index/)" << endl
         << "                           or GCSA2 index in PATH file (default: <graph>" << gcsa::GCSA::EXTENSION << ")" << endl
         << "                           (this is required if you are using multiple graphs files)" << endl
         << "    -t, --threads N        number of threads to use" << endl
         << "    -p, --progress         show progress" << endl
         << "    -V, --verify-index     validate the GCSA2 index using the input kmers (important for testing)" << endl
         << "rocksdb options (ignored with -g):" << endl
         << "    -s, --store-graph      store graph as xg" << endl
         << "    -m, --store-mappings   input is .gam format, store the mappings in alignments by node" << endl
         << "    -a, --store-alignments input is .gam format, store the alignments by node" << endl
         << "    -A, --dump-alignments  graph contains alignments, output them in sorted order" << endl
         << "    -P, --prune KB         remove kmer entries which use more than KB kilobytes" << endl
         << "    -n, --allow-negs       don't filter out relative negative positions of kmers" << endl
         << "    -D, --dump             print the contents of the db to stdout" << endl
         << "    -M, --metadata         describe aspects of the db stored in metadata" << endl
         << "    -L, --path-layout      describes the path layout of the graph" << endl
         << "    -S, --set-kmer         assert that the kmer size (-k) is in the db" << endl
        //<< "    -b, --tmp-db-base S    use this base name for temporary indexes" << endl
        << "    -C, --compact          compact the index into a single level (improves performance)" << endl
        << "    -Q, --use-snappy       use snappy compression (faster, larger) rather than zlib" << endl;

}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    string rocksdb_name;
    string gcsa_name;
    string xg_name;
    // Where should we import haplotype phasing paths from, if anywhere?
    string vcf_name;
    vector<string> dbg_names;
    int kmer_size = 0;
    bool path_only = false;
    int edge_max = 0;
    int kmer_stride = 1;
    int prune_kb = -1;
    bool store_graph = false;
    bool dump_index = false;
    bool describe_index = false;
    bool show_progress = false;
    bool set_kmer_size = false;
    bool path_layout = false;
    bool store_alignments = false;
    bool store_mappings = false;
    bool allow_negs = false;
    bool compact = false;
    bool dump_alignments = false;
    bool use_snappy = false;
    int doubling_steps = 2;
    bool verify_index = false;
    bool forward_only = false;
    size_t size_limit = 200; // in gigabytes
    bool store_threads = false; // use gPBWT to store paths

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"db-name", required_argument, 0, 'd'},
                {"kmer-size", required_argument, 0, 'k'},
                {"edge-max", required_argument, 0, 'e'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"store-graph", no_argument, 0, 's'},
                {"store-alignments", no_argument, 0, 'a'},
                {"dump-alignments", no_argument, 0, 'A'},
                {"store-mappings", no_argument, 0, 'm'},
                {"dump", no_argument, 0, 'D'},
                {"metadata", no_argument, 0, 'M'},
                {"set-kmer", no_argument, 0, 'S'},
                {"threads", required_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"prune",  required_argument, 0, 'P'},
                {"path-layout", no_argument, 0, 'L'},
                {"compact", no_argument, 0, 'C'},
                {"allow-negs", no_argument, 0, 'n'},
                {"use-snappy", no_argument, 0, 'Q'},
                {"gcsa-name", required_argument, 0, 'g'},
                {"xg-name", required_argument, 0, 'x'},
                {"vcf-phasing", required_argument, 0, 'v'},
                {"verify-index", no_argument, 0, 'V'},
                {"forward-only", no_argument, 0, 'F'},
                {"size-limit", no_argument, 0, 'Z'},
                {"path-only", no_argument, 0, 'O'},
                {"store-threads", no_argument, 0, 'T'},
                {"dbg-in", required_argument, 0, 'i'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:k:j:pDshMt:b:e:SP:LmaCnAQg:X:x:v:VFZ:Oi:T",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'd':
            rocksdb_name = optarg;
            break;

        case 'x':
            xg_name = optarg;
            break;

        case 'v':
            vcf_name = optarg;
            break;

        case 'P':
            prune_kb = atoi(optarg);
            break;

        case 'k':
            kmer_size = atoi(optarg);
            break;


        case 'O':
            path_only = true;
            break;

        case 'e':
            edge_max = atoi(optarg);
            break;


        case 'j':
            kmer_stride = atoi(optarg);
            break;

        case 'p':
            show_progress = true;
            break;

        case 'D':
            dump_index = true;
            break;

        case 'M':
            describe_index = true;
            break;

        case 'L':
            path_layout = true;
            break;

        case 'S':
            set_kmer_size = true;
            break;

        case 's':
            store_graph = true;
            break;

        case 'a':
            store_alignments = true;
            break;

        case 'A':
            dump_alignments = true;
            break;

        case 'm':
            store_mappings = true;
            break;

        case 'n':
            allow_negs = true;
            break;

        case 'C':
            compact = true;
            break;

        case 'Q':
            use_snappy = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'g':
            gcsa_name = optarg;
            break;

        case 'V':
            verify_index = true;
            break;
        case 'i':
            dbg_names.push_back(optarg);
            break;
        case 'F':
            forward_only = true;
            break;

        case 'X':
            doubling_steps = atoi(optarg);
            break;

        case 'Z':
            size_limit = atoi(optarg);
            break;

        case 'T':
            store_threads = true;
            break;

        case 'h':
        case '?':
            help_index(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (edge_max == 0) edge_max = kmer_size + 1;

    vector<string> file_names;
    while (optind < argc) {
        string file_name = argv[optind++];
        file_names.push_back(file_name);
    }

    if(kmer_size == 0 && !gcsa_name.empty() && dbg_names.empty()) {
        // gcsa doesn't do anything if we tell it a kmer size of 0.
        cerr << "error:[vg index] kmer size for GCSA2 index must be >0" << endl;
        return 1;
    }

    if(kmer_size < 0) {
        cerr << "error:[vg index] kmer size cannot be negative" << endl;
        return 1;
    }

    if(kmer_stride <= 0) {
        // kmer strides of 0 (or negative) are silly.
        cerr << "error:[vg index] kmer stride must be positive and nonzero" << endl;
        return 1;
    }

    if (!xg_name.empty()) {
        // We need to build an xg index

        // We'll fill this with the opened VCF file if we need one.
        vcflib::VariantCallFile variant_file;

        if(!vcf_name.empty()) {
            // There's a VCF we should load haplotype info from

            if (!vcf_name.empty()) {
                variant_file.open(vcf_name);
                if (!variant_file.is_open()) {
                    cerr << "error:[vg index] could not open" << vcf_name << endl;
                    return 1;
                }
            }

        }

        // We want to siphon off the "_alt_<variant>_<number>" paths from "vg
        // construct -a" and not index them, and use them for creating haplotype
        // threads.
        // TODO: a better way to store path metadata
        map<string, Path> alt_paths;
        // This is matched against the entire string.
        //regex is_alt("_alt_.+_[0-9]+");

        // store the graphs
        VGset graphs(file_names);
        // Turn into an XG index, except for the alt paths which we pull out and load into RAM instead.
        xg::XG index = graphs.to_xg(store_threads); //TODO graphs.to_xg(is_alt, alt_paths);

        if(variant_file.is_open()) {
            // Now go through and add the varaints.

            // How many phases are there?
            size_t num_samples = variant_file.sampleNames.size();
            // And how many phases?
            size_t num_phases = num_samples * 2;

            for(size_t path_rank = 1; path_rank <= index.max_path_rank(); path_rank++) {
                // Find all the reference paths and loop over them. We'll just
                // assume paths that don't start with "_" might appear in the
                // VCF. We need to use the xg path functions, since we didn't
                // load up the whole vg graph.

                // What path is this?
                string path_name = index.path_name(path_rank);

                // We already know it's not a variant's alt, since those were
                // removed, so it might be a primary contig.

                // How many bases is it?
                size_t path_length = index.path_length(path_name);

                // Allocate some Paths to store phase threads
                vector<Path> active_phase_paths{num_phases};
                // We need to remember how many paths of a particular phase have
                // already been generated.
                vector<int> saved_phase_paths(num_phases, 0);

                // What's the first reference position after the last variant?
                size_t nonvariant_start = 0;

                // Completed ones just get dumped into the index
                auto finish_phase = [&](size_t phase_number) {
                    // We have finished a phase (because an unphased variant
                    // came up or we ran out of variants); dump it into the
                    // index under a name and make a new Path for that phase.

                    // Find where this path is in our vector
                    Path& to_save = active_phase_paths[phase_number];

                    if(to_save.mapping_size() > 0) {
                        // Only actually do anything if we put in some mappings.

                        // Make sure the path has a name, and increment the
                        // number of saved paths for this phase so the next path
                        // will get a different name.
                        to_save.set_name("_phase_" + to_string(phase_number) +
                            "_" + to_string(saved_phase_paths[phase_number]++));

                        // Actually send the path off to XG
                        index.insert_thread(to_save);

                        // Clear it out for re-use
                        to_save = Path();
                    }
                };

                // We need a way to dump mappings into pahse threads. The
                // mapping edits and rank and offset info will be ignored; the
                // Mapping just represents an oriented node traversal.
                auto append_mapping = [&](size_t phase_number, const Mapping& mapping) {
                    // Find the path to add to
                    Path& to_extend = active_phase_paths[phase_number];

                    // See if the edge we need to follow exists
                    if(to_extend.mapping_size() > 0) {
                        // If there's a previous mapping, go find it
                        const Mapping& previous = to_extend.mapping(to_extend.mapping_size() - 1);

                        // Break out the IDs and flags we need to check for the edge
                        int64_t last_node = previous.position().node_id();
                        bool last_from_start = previous.position().is_reverse();

                        int64_t new_node = mapping.position().node_id();
                        bool new_to_end = mapping.position().is_reverse();

                        if(!index.has_edge(last_node, last_from_start, new_node, new_to_end)) {
                            // We can't have a thread take this edge. Split ane
                            // emit the current mappings and start a new path.
                            cerr << "warning:[vg index] phase " << phase_number << " wants edge "
                                << last_node << (last_from_start ? "L" : "R") << " - "
                                << new_node << (new_to_end ? "R" : "L")
                                << " which does not exist. Splitting!" << endl;

                            finish_phase(phase_number);
                        }
                    }

                    // Make a new Mapping in the Path
                    Mapping& new_place = *(active_phase_paths[phase_number].add_mapping());

                    // Set it
                    new_place = mapping;

                    // Make sure to clear out the rank and edits
                    new_place.set_rank(0);
                    new_place.clear_edit();
                };

                // We need an easy way to append any reference mappings from the
                // last variant up until a certain position (which may be past
                // the end of the entire reference path).
                auto append_reference_mappings_until = [&](size_t phase_number, size_t end) {
                    // We need to look and add in the mappings to the
                    // intervening reference nodes from the last variant, if
                    // any. For which we need access to the last variant's past-
                    // the-end reference position.
                    size_t ref_pos = nonvariant_start;
                    while(ref_pos < end) {
                        // While there is intervening reference
                        // sequence, add it to our phase.

                        // What mapping is here?
                        Mapping ref_mapping = index.mapping_at_path_position(path_name, ref_pos);

                        // Stick it in the phase path
                        append_mapping(phase_number, ref_mapping);

                        // Advance to what's after that mapping
                        ref_pos += index.node_length(ref_mapping.position().node_id());
                    }
                };

                // We also have another function to handle each variant as it comes in.
                auto handle_variant = [&](vcflib::Variant& variant) {
                    // So we have a variant

                    // Grab its id, or make one by hashing stuff if it doesn't
                    // have an ID.
                    string var_name = get_or_make_variant_id(variant);

                    if(alt_paths.count("_alt_" + var_name + "_0") == 0) {
                        // There isn't a reference alt path for this variant.
#ifdef debug
                        cerr << "Reference alt for " << var_name << " not in VG set! Skipping!" << endl;
#endif
                        // Don't bother with this variant
                        return;
                    }

                    for(int sample_number = 0; sample_number < num_samples; sample_number++) {
                        // For each sample

                        // What sample is it?
                        string& sample_name = variant_file.sampleNames[sample_number];

                        // Parse it out and see if it's phased.
                        string genotype = variant.getGenotype(sample_name);

                        // Find the phasing bar
                        auto bar_pos = genotype.find('|');

                        for(int phase_offset = 0; phase_offset < 2; phase_offset++) {
                            // For both the phases for the sample, add mappings
                            // through all the fixed reference nodes between the
                            // last variant and here.
                            append_reference_mappings_until(sample_number * 2 + phase_offset, variant.position);

                            // If this variant isn't phased, this will just be a
                            // reference-matching piece of thread after the last
                            // variant. If that wasn't phased either, it's just
                            // a floating perfect reference match.
                        }

                        if(bar_pos == string::npos || bar_pos == 0 || bar_pos + 1 >= genotype.size()) {
                            // If it isn't phased, or we otherwise don't like
                            // it, we need to break phasing paths.
                            for(int phase_offset = 0; phase_offset < 2; phase_offset++) {
                                // Finish both the phases for this sample.
                                finish_phase(sample_number * 2 + phase_offset);
                            }
                        }

                        // If it is phased, parse out the two alleles and handle
                        // each separately.
                        vector<int> alt_indices({stoi(genotype.substr(0, bar_pos)),
                            stoi(genotype.substr(bar_pos + 1))});

                        for(int phase_offset = 0; phase_offset < 2; phase_offset++) {
                            // Handle each phase and its alt
                            int& alt_index = alt_indices[phase_offset];

                            // We need to find the path for this alt of this
                            // variant. We can pull out the whole thing since it
                            // should be short.
                            Path alt_path = alt_paths.at("_alt_" + var_name + "_" + to_string(alt_index));
                            // TODO: if we can't find this path, it probaby
                            // means we mismatched the vg file and the vcf file.
                            // Maybe we should complain to the user instead of
                            // just failing an assert in at()?


                            for(size_t i = 0; i < alt_path.mapping_size(); i++) {
                                // Then blit mappings from the alt over to the phase thread
                                append_mapping(sample_number * 2 + phase_offset, alt_path.mapping(i));
                            }

                            // TODO: We can't really land anywhere on the other
                            // side of a deletion if the phasing breaks right at
                            // it, because we don't know that the first
                            // reference base after the deletion hasn't been
                            // replaced. TODO: can we inspect the next reference
                            // node and see if any alt paths touch it?
                        }

                        // Now we have processed both phasinbgs for this sample.
                    }

                    // Update the past-the-last-variant position, globally,
                    // after we do all the samples.
                    nonvariant_start = variant.position + variant.ref.size();
                };

                // Look for variants only on this path
                variant_file.setRegion(path_name);

                // Set up progress bar
                ProgressBar* progress = nullptr;
                // Message needs to last as long as the bar itself.
                string progress_message = "loading variants for " + path_name;
                if(show_progress) {
                    progress = new ProgressBar(path_length, progress_message.c_str());
                    progress->Progressed(0);
                }

                // TODO: For a first attempt, let's assume we can actually store
                // all the Path objects for all the phases.

                // Allocate a place to store actual variants
                vcflib::Variant var(variant_file);

                // How many variants have we done?
                size_t variants_processed = 0;
                while (variant_file.is_open() && variant_file.getNextVariant(var)) {
                    // this ... maybe we should remove it as for when we have calls against N
                    bool isDNA = allATGC(var.ref);
                    for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                        if (!allATGC(*a)) isDNA = false;
                    }
                    // only work with DNA sequences
                    if (!isDNA) {
                        continue;
                    }

                    var.position -= 1; // convert to 0-based

                    // Handle the variant
                    handle_variant(var);


                    if (variants_processed++ % 1000 == 0 && progress != nullptr) {
                        // Say we made progress
                        progress->Progressed(var.position);
                    }
                }

                // Now finish up all the threads
                for(size_t i = 0; i < num_phases; i++) {
                    // Each thread runs out until the end of the reference path
                    append_reference_mappings_until(i, path_length);

                    // And then we save all the threads
                    finish_phase(i);
                }

                if(progress != nullptr) {
                    // Throw out our progress bar
                    delete progress;
                }

            }


        }


        // save the xg version to the file name we've been given
        ofstream db_out(xg_name);
        index.serialize(db_out);
        db_out.close();
    }

    if(!gcsa_name.empty()) {
        // We need to make a gcsa index.

        // Load up the graphs
        vector<string> tmpfiles;
        if (dbg_names.empty()) {
            VGset graphs(file_names);
            graphs.show_progress = show_progress;
            // Go get the kmers of the correct size
            tmpfiles = graphs.write_gcsa_kmers_binary(kmer_size, path_only, forward_only);
        } else {
            tmpfiles = dbg_names;
        }
        // Make the index with the kmers
        gcsa::InputGraph input_graph(tmpfiles, true);
        gcsa::ConstructionParameters params;
        params.setSteps(doubling_steps);
        params.setLimit(size_limit);

        // build the GCSA index
        gcsa::GCSA* gcsa_index = new gcsa::GCSA(input_graph, params);

        if (verify_index) {
            //cerr << "verifying index" << endl;
            if (!gcsa_index->verifyIndex(input_graph)) {
                cerr << "[vg::main]: GCSA2 index verification failed" << endl;
            }
        }

        // build the LCP array
        string lcp_name = gcsa_name + ".lcp";
        gcsa::LCPArray* lcp_array = new gcsa::LCPArray(input_graph, params);

        // clean up input graph temp files
        if (dbg_names.empty()) {
            for (auto& tfn : tmpfiles) {
                remove(tfn.c_str());
            }
        }

        // Save the GCSA2 index
        sdsl::store_to_file(*gcsa_index, gcsa_name);
        delete gcsa_index;

        // Save the LCP array
        sdsl::store_to_file(*lcp_array, lcp_name);
        delete lcp_array;

    }

    if (!rocksdb_name.empty()) {

        Index index;
        index.use_snappy = use_snappy;

        if (compact) {
            index.open_for_write(rocksdb_name);
            index.compact();
            index.flush();
            index.close();
        }

        // todo, switch to xg for graph storage
        // index should write and load index/xg or such
        // then a handful of functions used in main.cpp and mapper.cpp need to be rewritten to use the xg index
        if (store_graph && file_names.size() > 0) {
            index.open_for_write(rocksdb_name);
            VGset graphs(file_names);
            graphs.show_progress = show_progress;
            graphs.store_in_index(index);
            //index.flush();
            //index.close();
            // reopen to index paths
            // this requires the index to be queryable
            //index.open_for_write(db_name);
            graphs.store_paths_in_index(index);
            index.compact();
            index.flush();
            index.close();
        }

        if (store_alignments && file_names.size() > 0) {
            index.open_for_write(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                index.put_alignment(aln);
            };
            for (auto& file_name : file_names) {
                if (file_name == "-") {
                    stream::for_each(std::cin, lambda);
                } else {
                    ifstream in;
                    in.open(file_name.c_str());
                    stream::for_each(in, lambda);
                }
            }
            index.flush();
            index.close();
        }

        if (dump_alignments) {
            vector<Alignment> output_buf;
            index.open_read_only(rocksdb_name);
            auto lambda = [&output_buf](const Alignment& aln) {
                output_buf.push_back(aln);
                stream::write_buffered(cout, output_buf, 1000);
            };
            index.for_each_alignment(lambda);
            stream::write_buffered(cout, output_buf, 0);
            index.close();
        }

        if (store_mappings && file_names.size() > 0) {
            index.open_for_write(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                const Path& path = aln.path();
                for (int i = 0; i < path.mapping_size(); ++i) {
                    index.put_mapping(path.mapping(i));
                }
            };
            for (auto& file_name : file_names) {
                if (file_name == "-") {
                    stream::for_each(std::cin, lambda);
                } else {
                    ifstream in;
                    in.open(file_name.c_str());
                    stream::for_each(in, lambda);
                }
            }
            index.flush();
            index.close();
        }

        if (kmer_size != 0 && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            VGset graphs(file_names);
            graphs.show_progress = show_progress;
            graphs.index_kmers(index, kmer_size, path_only, edge_max, kmer_stride, allow_negs);
            index.flush();
            index.close();
            // forces compaction
            index.open_for_write(rocksdb_name);
            index.flush();
            index.compact();
            index.close();
        }

        if (prune_kb >= 0) {
            if (show_progress) {
                cerr << "pruning kmers > " << prune_kb << " on disk from " << rocksdb_name << endl;
            }
            index.open_for_write(rocksdb_name);
            index.prune_kmers(prune_kb);
            index.compact();
            index.close();
        }

        if (set_kmer_size) {
            assert(kmer_size != 0);
            index.open_for_write(rocksdb_name);
            index.remember_kmer_size(kmer_size);
            index.close();
        }

        if (dump_index) {
            index.open_read_only(rocksdb_name);
            index.dump(cout);
            index.close();
        }

        if (describe_index) {
            index.open_read_only(rocksdb_name);
            set<int> kmer_sizes = index.stored_kmer_sizes();
            cout << "kmer sizes: ";
            for (auto kmer_size : kmer_sizes) {
                cout << kmer_size << " ";
            }
            cout << endl;
            index.close();
        }

        if (path_layout) {
            index.open_read_only(rocksdb_name);
            //index.path_layout();
            map<string, int64_t> path_by_id = index.paths_by_id();
            map<string, pair<pair<int64_t, bool>, pair<int64_t, bool>>> layout;
            map<string, int64_t> length;
            index.path_layout(layout, length);
            for (auto& p : layout) {
                // Negate IDs for backward nodes
                cout << p.first << " " << p.second.first.first * (p.second.first.second ? -1 : 1) << " "
                     << p.second.second.first * (p.second.second.second ? -1 : 1) << " " << length[p.first] << endl;
            }
            index.close();
        }
    }

    return 0;

}

void help_align(char** argv) {
    cerr << "usage: " << argv[0] << " align [options] <graph.vg> >alignments.gam" << endl
        << "options:" << endl
        << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
        << "    -Q, --seq-name STR    name the sequence using this value" << endl
        << "    -j, --json            output alignments in JSON format (default GAM)" << endl;
}

int main_align(int argc, char** argv) {

    string seq;
    string seq_name;

    if (argc == 2) {
        help_align(argv);
        return 1;
    }

    bool print_cigar = false;
    bool output_json = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"sequence", required_argument, 0, 's'},
            {"seq-name", no_argument, 0, 'Q'},
            {"json", no_argument, 0, 'j'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:jhQ:",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 's':
                seq = optarg;
                break;

            case 'Q':
                seq_name = optarg;
                break;

            case 'j':
                output_json = true;
                break;

            case 'h':
            case '?':
                /* getopt_long already printed an error message. */
                help_align(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    Alignment alignment = graph->align(seq);

    if (output_json) {
        cout << pb2json(alignment) << endl;
    } else {
        if (!seq_name.empty()) {
            alignment.set_name(seq_name);
        }
        function<Alignment(uint64_t)> lambda =
            [&alignment] (uint64_t n) {
                return alignment;
            };
        stream::write(cout, 1, lambda);
    }

    delete graph;

    return 0;

}

void help_map(char** argv) {
    cerr << "usage: " << argv[0] << " map [options] <graph.vg> >alignments.vga" << endl
         << "options:" << endl
         << "    -d, --db-name DIR     use this db (defaults to <graph>.index/)" << endl
         << "                          A graph is not required. But GCSA/xg take precedence if available." << endl
         << "    -x, --xg-name FILE    use this xg index (defaults to <graph>.xg)" << endl
         << "    -g, --gcsa-name FILE  use this GCSA2 index (defaults to <graph>" << gcsa::GCSA::EXTENSION << ")" << endl
         << "    -V, --in-memory       build the XG and GCSA2 indexes in-memory from the provided .vg file" << endl
         << "    -O, --in-mem-path-only  when making the in-memory temporary index, only look at embedded paths" << endl
         << "input:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -Q, --seq-name STR    name the sequence using this value (for graph modification with new named paths)" << endl
         << "    -r, --reads FILE      take reads (one per line) from FILE, write alignments to stdout" << endl
         << "    -b, --hts-input FILE  align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout" << endl
         << "    -K, --keep-secondary  produce alignments for secondary input alignments in addition to primary ones" << endl
         << "    -f, --fastq FILE      input fastq (possibly compressed), two are allowed, one for each mate" << endl
         << "    -i, --interleaved     fastq is interleaved paired-ended" << endl
         << "    -p, --pair-window N   align to a graph up to N ids away from the mapping location of one mate for the other" << endl
         << "    -N, --sample NAME     for --reads input, add this sample" << endl
         << "    -R, --read-group NAME for --reads input, add this read group" << endl
         << "output:" << endl
         << "    -J, --output-json     output JSON rather than an alignment stream (helpful for debugging)" << endl
         << "    -B, --band-width N    for very long sequences, align in chunks then merge paths (default 1000bp)" << endl
         << "    -D, --debug           print debugging information about alignment to stderr" << endl
         << "generic mapping parameters:" << endl
         << "    -P, --min-score N     accept alignment only if the normalized alignment score is > N (default: 0)" << endl
         << "    -n, --context-depth N follow this many edges out from each thread for alignment (default: 5)" << endl
         << "    -M, --max-multimaps N produce up to N alignments for each read (default: 1)" << endl
         << "    -T, --softclip-trig N trigger graph extension and realignment when either end has softclips (default: 0)" << endl
         << "    -m, --hit-max N       ignore kmers or MEMs who have >N hits in our index (default: 100)" << endl
         << "maximal exact match (MEM) mapper:" << endl
         << "  This algorithm is used when --kmer-size is not specified and a GCSA index is given" << endl
         << "    -L, --min-mem-length N   ignore MEMs shorter than this length (default: 0/unset)" << endl
         << "    -Y, --max-mem-length N   ignore MEMs longer than this length by stopping backward search (default: 0/unset)" << endl
         << "kmer-based mapper:" << endl
         << "  This algorithm is used when --kmer-size is specified or a rocksdb index is given" << endl
         << "    -k, --kmer-size N     use this kmer size, it must be < kmer size in db (default: from index)" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers to use for seeding (default: kmer size)" << endl
         << "    -E, --min-kmer-entropy N  require shannon entropy of this in order to use kmer (default: no limit)" << endl
         << "    -S, --sens-step N     decrease kmer size by N bp until alignment succeeds (default: 5)" << endl
         << "    -A, --max-attempts N  try to improve sensitivity and align this many times (default: 7)" << endl
         << "    -l, --kmer-min N      give up aligning if kmer size gets below this threshold (default: 8)" << endl
         << "    -e, --thread-ex N     grab this many nodes in id space around each thread for alignment (default: 10)" << endl
         << "    -c, --clusters N      use at most the largest N ordered clusters of the kmer graph for alignment (default: all)" << endl
         << "    -C, --cluster-min N   require at least this many kmer hits in a cluster to attempt alignment (default: 2)" << endl
         << "    -H, --max-target-x N  skip cluster subgraphs with length > N*read_length (default: 100; unset: 0)" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -F, --prefer-forward  if the forward alignment of the read works, accept it" << endl
         << "    -G, --greedy-accept   if a tested alignment achieves -X score/bp don't try worse seeds" << endl
         << "    -X, --accept-score N  accept early alignment if the normalized alignment score is > N and -F or -G is set" << endl;
}

int main_map(int argc, char** argv) {

    if (argc == 2) {
        help_map(argv);
        return 1;
    }

    string seq;
    string seq_name;
    string db_name;
    string xg_name;
    string gcsa_name;
    int kmer_size = 0;
    int kmer_stride = 0;
    int sens_step = 0;
    int best_clusters = 0;
    int cluster_min = 2;
    int max_attempts = 7;
    string read_file;
    string hts_file;
    bool keep_secondary = false;
    int hit_max = 100;
    int max_multimaps = 1;
    int thread_count = 1;
    int thread_ex = 10;
    int context_depth = 5;
    bool output_json = false;
    bool debug = false;
    bool prefer_forward = false;
    bool greedy_accept = false;
    float min_score = 0;
    string sample_name;
    string read_group;
    string fastq1, fastq2;
    bool interleaved_fastq = false;
    int pair_window = 64; // ~11bp/node
    int band_width = 1000; // anything > 1000bp sequences is difficult to align efficiently
    bool try_both_mates_first = false;
    float min_kmer_entropy = 0;
    float accept_score = 0;
    size_t kmer_min = 8;
    int softclip_threshold = 0;
    bool build_in_memory = false;
    int max_mem_length = 0;
    int min_mem_length = 0;
    int max_target_factor = 100;
    bool in_mem_path_only = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {"seq-name", required_argument, 0, 'Q'},
                {"db-name", required_argument, 0, 'd'},
                {"xg-name", required_argument, 0, 'x'},
                {"gcsa-name", required_argument, 0, 'g'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"kmer-size", required_argument, 0, 'k'},
                {"min-kmer-entropy", required_argument, 0, 'E'},
                {"clusters", required_argument, 0, 'c'},
                {"cluster-min", required_argument, 0, 'C'},
                {"max-attempts", required_argument, 0, 'A'},
                {"reads", required_argument, 0, 'r'},
                {"sample", required_argument, 0, 'N'},
                {"read-group", required_argument, 0, 'R'},
                {"hit-max", required_argument, 0, 'm'},
                {"max-multimaps", required_argument, 0, 'N'},
                {"threads", required_argument, 0, 't'},
                {"prefer-forward", no_argument, 0, 'F'},
                {"greedy-accept", no_argument, 0, 'G'},
                {"accept-score", required_argument, 0, 'X'},
                {"sens-step", required_argument, 0, 'S'},
                {"thread-ex", required_argument, 0, 'e'},
                {"context-depth", required_argument, 0, 'n'},
                {"output-json", no_argument, 0, 'J'},
                {"hts-input", required_argument, 0, 'b'},
                {"keep-secondary", no_argument, 0, 'K'},
                {"fastq", no_argument, 0, 'f'},
                {"interleaved", no_argument, 0, 'i'},
                {"pair-window", required_argument, 0, 'p'},
                {"band-width", required_argument, 0, 'B'},
                {"min-score", required_argument, 0, 'P'},
                {"kmer-min", required_argument, 0, 'l'},
                {"softclip-trig", required_argument, 0, 'T'},
                {"in-memory", no_argument, 0, 'V'},
                {"in-mem-path-only", no_argument, 0, 'O'},
                {"debug", no_argument, 0, 'D'},
                {"min-mem-length", required_argument, 0, 'L'},
                {"max-mem-length", required_argument, 0, 'Y'},
                {"max-target-x", required_argument, 0, 'H'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:j:hd:x:g:c:r:m:k:M:t:DX:FS:Jb:KR:N:if:p:B:h:GC:A:E:Q:n:P:l:e:T:VL:Y:H:O",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 's':
            seq = optarg;
            break;

        case 'V':
            build_in_memory = true;
            break;

        case 'd':
            db_name = optarg;
            break;

        case 'x':
            xg_name = optarg;
            break;

        case 'g':
            gcsa_name = optarg;
            break;

        case 'j':
            kmer_stride = atoi(optarg);
            break;


        case 'O':
            in_mem_path_only = true;
            break;

        case 'Q':
            seq_name = optarg;
            break;

        case 'S':
            sens_step = atoi(optarg);
            break;

        case 'c':
            best_clusters = atoi(optarg);
            break;

        case 'C':
            cluster_min = atoi(optarg);
            break;

        case 'E':
            min_kmer_entropy = atof(optarg);
            break;

        case 'A':
            max_attempts = atoi(optarg);
            break;

        case 'm':
            hit_max = atoi(optarg);
            break;

        case 'M':
            max_multimaps = atoi(optarg);
            break;

        case 'k':
            kmer_size = atoi(optarg);
            break;

        case 'e':
            thread_ex = atoi(optarg);
            break;

        case 'n':
            context_depth = atoi(optarg);
            break;

        case 'T':
            softclip_threshold = atoi(optarg);
            break;

        case 'r':
            read_file = optarg;
            break;

        case 'R':
            read_group = optarg;
            break;

        case 'N':
            sample_name = optarg;
            break;

        case 'b':
            hts_file = optarg;
            break;

        case 'K':
            keep_secondary = true;
            break;

        case 'f':
            if (fastq1.empty()) fastq1 = optarg;
            else if (fastq2.empty()) fastq2 = optarg;
            else { cerr << "[vg map] error: more than two fastqs specified" << endl; exit(1); }
            break;

        case 'i':
            interleaved_fastq = true;
            break;

        case 'p':
            pair_window = atoi(optarg);
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'D':
            debug = true;
            break;

        case 'F':
            prefer_forward = true;
            break;

        case 'G':
            greedy_accept = true;
            break;

        case 'X':
            accept_score = atof(optarg);
            break;

        case 'J':
            output_json = true;
            break;

        case 'B':
            band_width = atoi(optarg);
            break;

        case 'P':
            min_score = atof(optarg);
            break;

        case 'l':
            kmer_min = atoi(optarg);
            break;

        case 'L':
            min_mem_length = atoi(optarg);
            break;

        case 'Y':
            max_mem_length = atoi(optarg);
            break;

        case 'H':
            max_target_factor = atoi(optarg);
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_map(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (seq.empty() && read_file.empty() && hts_file.empty() && fastq1.empty()) {
        cerr << "error:[vg map] a sequence or read file is required when mapping" << endl;
        return 1;
    }

    // should probably disable this
    string file_name;
    if (optind < argc) {
        file_name = argv[optind];
    }

    if (gcsa_name.empty() && !file_name.empty()) {
        gcsa_name = file_name + gcsa::GCSA::EXTENSION;
    }

    // Load up our indexes.
    xg::XG* xindex = nullptr;
    gcsa::GCSA* gcsa = nullptr;
    gcsa::LCPArray* lcp = nullptr;

    // for testing, we sometimes want to run the mapper on indexes we build in memory
    if (build_in_memory) {
        VG* graph;
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
        xindex = new xg::XG(graph->graph);
        assert(kmer_size);
        int doubling_steps = 2;
        graph->build_gcsa_lcp(gcsa, lcp, kmer_size, in_mem_path_only, false, 2);
        delete graph;
    } else {
        // We try opening the file, and then see if it worked
        ifstream xg_stream(xg_name);

        if(xg_stream) {
            // We have an xg index!
            if(debug) {
                cerr << "Loading xg index " << xg_name << "..." << endl;
            }
            xindex = new xg::XG(xg_stream);
        }

        ifstream gcsa_stream(gcsa_name);
        if(gcsa_stream) {
            // We have a GCSA index too!
            if(debug) {
                cerr << "Loading GCSA2 index " << gcsa_name << "..." << endl;
            }
            gcsa = new gcsa::GCSA();
            gcsa->load(gcsa_stream);
        }

        string lcp_name = gcsa_name + ".lcp";
        ifstream lcp_stream(lcp_name);
        if (lcp_stream) {
            if(debug) {
                cerr << "Loading LCP index " << gcsa_name << "..." << endl;
            }
            lcp = new gcsa::LCPArray();
            lcp->load(lcp_stream);
        }
    }

    Index* idx = nullptr;

    if(!xindex || !gcsa) {
        // We only need a Rocksdb index if we don't have the others.
        if(debug) {
            cerr << "Loading RocksDB index " << db_name << "..." << endl;
        }
        idx = new Index();
        idx->open_read_only(db_name);
    }

    thread_count = get_thread_count();

    vector<Mapper*> mapper;
    mapper.resize(thread_count);
    vector<vector<Alignment> > output_buffer;
    output_buffer.resize(thread_count);

    // We have one function to dump alignments into
    // Make sure to flush the buffer at the end of the program!
    auto output_alignments = [&output_buffer, &output_json](vector<Alignment>& alignments) {
        // for(auto& alignment : alignments){
        // 		cerr << "This is in output_alignments" << alignment.DebugString() << endl;
        // }

        if (output_json) {
            // If we want to convert to JSON, convert them all to JSON and dump them to cout.
            for(auto& alignment : alignments) {
                string json = pb2json(alignment);
#pragma omp critical (cout)
                cout << json << "\n";
            }
        } else {
            // Otherwise write them through the buffer for our thread
            int tid = omp_get_thread_num();
            auto& output_buf = output_buffer[tid];

            // Copy all the alignments over to the output buffer
            copy(alignments.begin(), alignments.end(), back_inserter(output_buf));

            stream::write_buffered(cout, output_buf, 1000);
        }
    };

    for (int i = 0; i < thread_count; ++i) {
        Mapper* m;
        if(xindex && gcsa && lcp) {
            // We have the xg and GCSA indexes, so use them
            m = new Mapper(xindex, gcsa, lcp);
        } else {
            // Use the Rocksdb index and maybe the GCSA one
            m = new Mapper(idx, gcsa);
        }
        m->best_clusters = best_clusters;
        m->hit_max = hit_max;
        m->max_multimaps = max_multimaps;
        m->debug = debug;
        if (accept_score) m->accept_norm_score = accept_score;
        if (sens_step) m->kmer_sensitivity_step = sens_step;
        m->prefer_forward = prefer_forward;
        m->greedy_accept = greedy_accept;
        m->thread_extension = thread_ex;
        m->cluster_min = cluster_min;
        m->context_depth = context_depth;
        m->max_attempts = max_attempts;
        m->min_kmer_entropy = min_kmer_entropy;
        m->kmer_min = kmer_min;
        m->min_norm_score = min_score;
        m->softclip_threshold = softclip_threshold;
        m->min_mem_length = min_mem_length;
        m->max_mem_length = max_mem_length;
        m->max_target_factor = max_target_factor;
        mapper[i] = m;
    }

    if (!seq.empty()) {
        int tid = omp_get_thread_num();

        Alignment unaligned;
        unaligned.set_sequence(seq);
        vector<Alignment> alignments = mapper[tid]->align_multi(unaligned, kmer_size, kmer_stride, band_width);
        if(alignments.size() == 0) {
            // If we didn't have any alignments, report the unaligned alignment
            alignments.push_back(unaligned);
        }


        for(auto& alignment : alignments) {
            if (!sample_name.empty()) alignment.set_sample_name(sample_name);
            if (!read_group.empty()) alignment.set_read_group(read_group);
            if (!seq_name.empty()) alignment.set_name(seq_name);
        }

        // Output the alignments in JSON or protobuf as appropriate.
        output_alignments(alignments);
    }

    if (!read_file.empty()) {
        ifstream in(read_file);
        bool more_data = in.good();
#pragma omp parallel shared(in)
        {
            string line;
            int tid = omp_get_thread_num();
            while (in.good()) {
                line.clear();
#pragma omp critical (readq)
                {
                    std::getline(in,line);
                }
                if (!line.empty()) {
                    // Make an alignment
                    Alignment unaligned;
                    unaligned.set_sequence(line);

                    vector<Alignment> alignments = mapper[tid]->align_multi(unaligned, kmer_size, kmer_stride, band_width);
                    if(alignments.empty()) {
                        alignments.push_back(unaligned);
                    }

                    for(auto& alignment : alignments) {
                        // Set the alignment metadata
                        if (!sample_name.empty()) alignment.set_sample_name(sample_name);
                        if (!read_group.empty()) alignment.set_read_group(read_group);
                    }


                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alignments);
                }
            }
        }
    }

    if (!hts_file.empty()) {
        function<void(Alignment&)> lambda =
            [&mapper,
            &output_alignments,
            &keep_secondary,
            &kmer_size,
            &kmer_stride,
            &band_width]
                (Alignment& alignment) {

                    if(alignment.is_secondary() && !keep_secondary) {
                        // Skip over secondary alignments in the input; we don't want several output mappings for each input *mapping*.
                        return;
                    }

                    int tid = omp_get_thread_num();
                    vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, band_width);
                    if(alignments.empty()) {
                        alignments.push_back(alignment);
                    }

                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alignments);
                };
        // run
        hts_for_each_parallel(hts_file, lambda);
    }

    if (!fastq1.empty()) {
        if (interleaved_fastq) {
            // paired interleaved
            function<void(Alignment&, Alignment&)> lambda =
                [&mapper,
                &output_alignments,
                &kmer_size,
                &kmer_stride,
                &band_width,
                &pair_window]
                    (Alignment& aln1, Alignment& aln2) {

                        int tid = omp_get_thread_num();
                        auto alnp = mapper[tid]->align_paired_multi(aln1, aln2, kmer_size, kmer_stride, band_width, pair_window);

                        // Make sure we have unaligned "alignments" for things that don't align.
                        if(alnp.first.empty()) {
                            alnp.first.push_back(aln1);
                        }
                        if(alnp.second.empty()) {
                            alnp.second.push_back(aln2);
                        }

                        // Output the alignments in JSON or protobuf as appropriate.
                        output_alignments(alnp.first);
                        output_alignments(alnp.second);
                    };
            fastq_paired_interleaved_for_each_parallel(fastq1, lambda);
        } else if (fastq2.empty()) {
            // single
            function<void(Alignment&)> lambda =
                [&mapper,
                &output_alignments,
                &kmer_size,
                &kmer_stride,
                &band_width]
                    (Alignment& alignment) {

                        int tid = omp_get_thread_num();
                        vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, band_width);

                        if(alignments.empty()) {
                            // Make sure we have a "no alignment" alignment
                            alignments.push_back(alignment);
                        }

                        //cerr << "This is just before output_alignments" << alignment.DebugString() << endl;
                        output_alignments(alignments);
                    };
            fastq_unpaired_for_each_parallel(fastq1, lambda);
        } else {
            // paired two-file
            function<void(Alignment&, Alignment&)> lambda =
                [&mapper,
                &output_alignments,
                &kmer_size,
                &kmer_stride,
                &band_width,
                &pair_window]
                    (Alignment& aln1, Alignment& aln2) {

                        int tid = omp_get_thread_num();
                        auto alnp = mapper[tid]->align_paired_multi(aln1, aln2, kmer_size, kmer_stride, band_width, pair_window);

                        // Make sure we have unaligned "alignments" for things that don't align.
                        if(alnp.first.empty()) {
                            alnp.first.push_back(aln1);
                        }
                        if(alnp.second.empty()) {
                            alnp.second.push_back(aln2);
                        }

                        output_alignments(alnp.first);
                        output_alignments(alnp.second);
                    };
            fastq_paired_two_files_for_each_parallel(fastq1, fastq2, lambda);
        }
    }

    // clean up
    for (int i = 0; i < thread_count; ++i) {
        delete mapper[i];
        auto& output_buf = output_buffer[i];
        if (!output_json) {
            stream::write_buffered(cout, output_buf, 0);
        }
    }

    if(idx)  {
        delete idx;
        idx = nullptr;
    }
    if(gcsa) {
        delete gcsa;
        gcsa = nullptr;
    }
    if(xindex) {
        delete xindex;
        xindex = nullptr;
    }

    cout.flush();

    return 0;

}

void help_view(char** argv) {
    cerr << "usage: " << argv[0] << " view [options] [ <graph.vg> | <graph.json> | <aln.gam> | <read1.fq> [<read2.fq>] ]" << endl
        << "options:" << endl
        << "    -g, --gfa            output GFA format (default)" << endl
        << "    -F, --gfa-in         input GFA format" << endl

        << "    -v, --vg             output VG format" << endl
        << "    -V, --vg-in          input VG format (default)" << endl

        << "    -j, --json           output JSON format" << endl
        << "    -J, --json-in        input JSON format" << endl
        << "    -c, --json-stream    streaming conversion of a VG format graph in line delimited JSON format" << endl
        << "                         (this cannot be loaded directly via -J)" << endl
        << "    -G, --gam            output GAM format (vg alignment format: Graph " << endl
        << "                         Alignment/Map)" << endl
        << "    -t, --turtle         output RDF/turtle format (can not be loaded by VG)" << endl
        << "    -T  --turtle-in      input turtle format." << endl
        << "    -r, --rdf_base_uri   set base uri for the RDF output" << endl

        << "    -a, --align-in       input GAM format" << endl
        << "    -A, --aln-graph GAM  add alignments from GAM to the graph" << endl

        << "    -d, --dot            output dot format" << endl
        << "    -S, --simple-dot     simplify the dot output; remove node labels, simplify alignments" << endl
        << "    -C, --color          color nodes that are not in the reference path (DOT OUTPUT ONLY)" << endl
        << "    -p, --show-paths     show paths in dot output" << endl
        << "    -w, --walk-paths     add labeled edges to represent paths in dot output" << endl
        << "    -n, --annotate-paths add labels to normal edges to represent paths in dot output" << endl
        << "    -M, --show-mappings  with -p print the mappings in each path in JSON" << endl
        << "    -I, --invert-ports   invert the edge ports in dot so that ne->nw is reversed" << endl
        << "    -s, --random-seed N  use this seed when assigning path symbols in dot output" << endl

        << "    -b, --bam            input BAM or other htslib-parseable alignments" << endl

        << "    -f, --fastq          input fastq (output defaults to GAM). Takes two " << endl
        << "                         positional file arguments if paired" << endl
        << "    -i, --interleaved    fastq is interleaved paired-ended" << endl

        << "    -L, --pileup         ouput VG Pileup format" << endl
        << "    -l, --pileup-in      input VG Pileup format" << endl;
// TODO: Can we regularize the option names for input and output types?

}

int main_view(int argc, char** argv) {

    if (argc == 2) {
        help_view(argv);
        return 1;
    }

    // Supported conversions:
    //      TO  vg  json    gfa gam bam fastq   dot
    // FROM
    // vg       Y   Y       Y   N   N   N       Y
    // json     Y   Y       Y   N   N   N       Y
    // gfa      Y   Y       Y   N   N   N       Y
    // gam      N   Y       N   N   N   N       N
    // bam      N   N       N   Y   N   N       N
    // fastq    N   N       N   Y   N   N       N
    // dot      N   N       N   N   N   N       N
    //
    // and json-gam -> gam
    //     json-pileup -> pileup

    string output_type;
    string input_type;
    string rdf_base_uri;
    bool input_json = false;
    string alignments;
    string fastq1, fastq2;
    bool interleaved_fastq = false;
    bool show_paths_in_dot = false;
    bool walk_paths_in_dot = false;
    bool annotate_paths_in_dot = false;
    bool invert_edge_ports_in_dot = false;
    bool show_mappings_in_dot = false;
    bool simple_dot = false;
    int seed_val = time(NULL);
    bool color_variants = false;

    int c;
    optind = 2; // force optind past "view" argument
    while (true) {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            //{"verbose", no_argument,       &verbose_flag, 1},
            {"dot", no_argument, 0, 'd'},
            {"gfa", no_argument, 0, 'g'},
            {"turtle", no_argument, 0, 't'},
            {"rdf-base-uri", no_argument, 0, 'r'},
            {"gfa-in", no_argument, 0, 'F'},
            {"json",  no_argument, 0, 'j'},
            {"json-in",  no_argument, 0, 'J'},
            {"json-stream", no_argument, 0, 'c'},
            {"vg", no_argument, 0, 'v'},
            {"vg-in", no_argument, 0, 'V'},
            {"align-in", no_argument, 0, 'a'},
            {"gam", no_argument, 0, 'G'},
            {"bam", no_argument, 0, 'b'},
            {"fastq", no_argument, 0, 'f'},
            {"interleaved", no_argument, 0, 'i'},
            {"aln-graph", required_argument, 0, 'A'},
            {"show-paths", no_argument, 0, 'p'},
            {"turtle-in", no_argument, 0, 'T'},
            {"walk-paths", no_argument, 0, 'w'},
            {"annotate-paths", no_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"pileup", no_argument, 0, 'L'},
            {"pileup-in", no_argument, 0, 'l'},
            {"invert-ports", no_argument, 0, 'I'},
            {"show-mappings", no_argument, 0, 'M'},
            {"simple-dot", no_argument, 0, 'S'},
            {"color", no_argument, 0, 'C'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgFjJhvVpaGbifA:s:wnlLIMcTtr:SC",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 'C':
                color_variants = true;
                break;

            case 'd':
                output_type = "dot";
                break;

            case 'S':
                simple_dot = true;
                break;

            case 'p':
                show_paths_in_dot = true;
                break;

            case 'M':
                show_mappings_in_dot = true;
                break;

            case 'w':
                walk_paths_in_dot = true;
                break;


            case 'n':
                annotate_paths_in_dot = true;
                break;

            case 's':
                seed_val = atoi(optarg);
                break;

            case 'g':
                output_type = "gfa";
                break;

            case 'F':
                input_type = "gfa";
                break;

            case 'j':
                output_type = "json";
                break;

            case 'J':
                // -J can complement input GAM/Pileup, hence the extra logic here.
                if (input_type.empty()) {
                    input_type = "json";
                }
                input_json = true;
                break;

            case 'c':
                input_type = "vg";
                output_type = "stream";
                break;

            case 'v':
                output_type = "vg";
                break;

            case 'V':
                input_type = "vg";
                break;

            case 'G':
                output_type = "gam";
                break;

            case 't':
                output_type = "turtle";
                break;

        case 'r':
            rdf_base_uri = optarg;
            break;
        case 'T':
            input_type= "turtle-in";
            break;
        case 'a':
            input_type = "gam";
            if(output_type.empty()) {
                // Default to GAM -> JSON
                output_type = "json";
            }
            break;

            case 'b':
                input_type = "bam";
                if(output_type.empty()) {
                    // Default to BAM -> GAM, since BAM isn't convertable to our normal default.
                    output_type = "gam";
                }
                break;

            case 'f':
                input_type = "fastq";
                if(output_type.empty()) {
                    // Default to FASTQ -> GAM
                    output_type = "gam";
                }
                break;

            case 'i':
                interleaved_fastq = true;
                break;

            case 'A':
                alignments = optarg;
                break;

            case 'I':
                invert_edge_ports_in_dot = true;
                break;

            case 'L':
                output_type = "pileup";
                break;

            case 'l':
                input_type = "pileup";
                if (output_type.empty()) {
                    // Default to Pileup -> JSON
                    output_type = "json";
                }
                break;

            case 'h':
            case '?':
                /* getopt_long already printed an error message. */
                help_view(argv);
                exit(1);
                break;

            default:
                abort ();
        }
    }

    // If the user specified nothing else, we default to VG in and GFA out.
    if (input_type.empty()) {
        input_type = "vg";
    }
    if (output_type.empty()) {
        output_type = "gfa";
    }
    if (rdf_base_uri.empty()) {
        rdf_base_uri = "http://example.org/vg/";
    }
    vector<Alignment> alns;
    if (!alignments.empty()) {
        function<void(Alignment&)> lambda = [&alns](Alignment& aln) { alns.push_back(aln); };
        ifstream in;
        in.open(alignments.c_str());
        stream::for_each(in, lambda);
    }

    VG* graph = nullptr;
    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }
    string file_name = argv[optind];
    if (input_type == "vg") {
        if (output_type == "stream") {
            function<void(Graph&)> lambda = [&](Graph& g) { cout << pb2json(g) << endl; };
            if (file_name == "-") {
                stream::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each(in, lambda);
            }
            return 0;
        } else {
            if (file_name == "-") {
                graph = new VG(std::cin);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                graph = new VG(in);
            }
        }
        // VG can convert to any of the graph formats, so keep going
    } else if (input_type == "gfa") {
        if (file_name == "-") {
            graph = new VG;
            graph->from_gfa(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG;
            graph->from_gfa(in);
        }
        // GFA can convert to any of the graph formats, so keep going
    } else if(input_type == "json") {
        assert(input_json == true);
        JSONStreamHelper<Graph> json_helper(file_name);
        function<bool(Graph&)> get_next_graph = json_helper.get_read_fn();
        graph = new VG(get_next_graph, false);
    } else if(input_type == "turtle-in") {
        graph = new VG;
        bool pre_compress=color_variants;
        if (file_name == "-") {
             graph->from_turtle("/dev/stdin", rdf_base_uri);
        } else {
             graph->from_turtle(file_name, rdf_base_uri);
	}
    } else if (input_type == "gam") {
        if (input_json == false) {
            if (output_type == "json") {
                // convert values to printable ones
                function<void(Alignment&)> lambda = [](Alignment& a) {
                    alignment_quality_short_to_char(a);
                    cout << pb2json(a) << "\n";
                };
                if (file_name == "-") {
                    stream::for_each(std::cin, lambda);
                } else {
                    ifstream in;
                    in.open(file_name.c_str());
                    stream::for_each(in, lambda);
                }
            } else {
                // todo
                cerr << "[vg view] error: (binary) GAM can only be converted to JSON" << endl;
                return 1;
            }
        } else {
            if (output_type == "json" || output_type == "gam") {
                JSONStreamHelper<Alignment> json_helper(file_name);
                json_helper.write(cout, output_type == "json");
            } else {
                cerr << "[vg view] error: JSON GAM can only be converted to GAM or JSON" << endl;
                return 1;
            }
        }
        cout.flush();
        return 0;
    } else if (input_type == "bam") {
        if (output_type == "gam") {
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
            hts_for_each(file_name, lambda);
            write_alignments(std::cout, buf);
            buf.clear();
            cout.flush();
            return 0;
        } else if (output_type == "json") {
            // todo
            cerr << "[vg view] error: BAM to JSON conversion not yet implemented" << endl;
            return 0;
        } else {
            cerr << "[vg view] error: BAM can only be converted to GAM" << endl;
            return 1;
        }
        } else if (input_type == "fastq") {
            fastq1 = argv[optind++];
            if (optind < argc) {
                fastq2 = argv[optind];
            }
            if (output_type == "gam") {
                vector<Alignment> buf;
                if (!interleaved_fastq && fastq2.empty()) {
                    function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
                        buf.push_back(aln);
                        if (buf.size() > 1000) {
                            write_alignments(std::cout, buf);
                            buf.clear();
                        }
                    };
                    fastq_unpaired_for_each(fastq1, lambda);
                } else if (interleaved_fastq && fastq2.empty()) {
                    function<void(Alignment&, Alignment&)> lambda = [&buf](Alignment& aln1, Alignment& aln2) {
                        buf.push_back(aln1);
                        buf.push_back(aln2);
                        if (buf.size() > 1000) {
                            write_alignments(std::cout, buf);
                            buf.clear();
                        }
                    };
                    fastq_paired_interleaved_for_each(fastq1, lambda);
                } else if (!fastq2.empty()) {
                    function<void(Alignment&, Alignment&)> lambda = [&buf](Alignment& aln1, Alignment& aln2) {
                        buf.push_back(aln1);
                        buf.push_back(aln2);
                        if (buf.size() > 1000) {
                            write_alignments(std::cout, buf);
                            buf.clear();
                        }
                    };
                    fastq_paired_two_files_for_each(fastq1, fastq2, lambda);
                }
                write_alignments(std::cout, buf);
                buf.clear();
            } else {
                // We can't convert fastq to the other graph formats
                cerr << "[vg view] error: FASTQ can only be converted to GAM" << endl;
                return 1;
            }
            cout.flush();
            return 0;
        } else if (input_type == "pileup") {
            if (input_json == false) {
                if (output_type == "json") {
                    // convert values to printable ones
                    function<void(NodePileup&)> lambda = [](NodePileup& p) {
                        cout << pb2json(p) << "\n";
                    };
                    if (file_name == "-") {
                        stream::for_each(std::cin, lambda);
                    } else {
                        ifstream in;
                        in.open(file_name.c_str());
                        stream::for_each(in, lambda);
                    }
                } else {
                    // todo
                    cerr << "[vg view] error: (binary) Pileup can only be converted to JSON" << endl;
                    return 1;
                }
            } else {
                if (output_type == "json" || output_type == "pileup") {
                    JSONStreamHelper<NodePileup> json_helper(file_name);
                    json_helper.write(cout, output_type == "json");
                } else {
                    cerr << "[vg view] error: JSON Pileup can only be converted to Pileup or JSON" << endl;
                    return 1;
                }
            }
            cout.flush();
            return 0;
        }

        if(graph == nullptr) {
            // Make sure we didn't forget to implement an input format.
            cerr << "[vg view] error: cannot load graph in " << input_type << " format" << endl;
            return 1;
        }

        if(!graph->is_valid()) {
            // If we're converting the graph, we might as well make sure it's valid.
            // This is especially useful for JSON import.
            cerr << "[vg view] warning: graph is invalid!" << endl;
        }

        // Now we know graph was filled in from the input format. Spit it out in the
        // requested output format.

        if (output_type == "dot") {
            graph->to_dot(std::cout,
                    alns,
                    show_paths_in_dot,
                    walk_paths_in_dot,
                    annotate_paths_in_dot,
                    show_mappings_in_dot,
                    simple_dot,
                    invert_edge_ports_in_dot,
                    color_variants,
                    seed_val);
        } else if (output_type == "json") {
            cout << pb2json(graph->graph) << endl;
        } else if (output_type == "gfa") {
            graph->to_gfa(std::cout);
        } else if (output_type == "turtle") {
            graph->to_turtle(std::cout, rdf_base_uri, color_variants);
        } else if (output_type == "vg") {
            graph->serialize_to_ostream(cout);
        } else {
            // We somehow got here with a bad output format.
            cerr << "[vg view] error: cannot save a graph in " << output_type << " format" << endl;
            return 1;
        }

        cout.flush();
        delete graph;

        return 0;
    }

    void help_deconstruct(char** argv){
        cerr << "usage: " << argv[0] << " deconstruct [options] <my_graph>.vg" << endl
            << "options: " << endl
            << endl;
    }
int main_deconstruct(int argc, char** argv){
        cerr << "WARNING: EXPERIMENTAL" << endl;
            if (argc <= 2) {
            help_deconstruct(argv);
            return 1;
        }

        int c;
        optind = 2; // force optind past command positional argument
        while (true) {
            static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"xg-name", required_argument,0, 'g'},
                {0, 0, 0, 0}

            };
            int option_index = 0;
            c = getopt_long (argc, argv, "hg:",
                    long_options, &option_index);

            // Detect the end of the options.
            if (c == -1)
                break;

            switch (c)
            {
                case '?':
                case 'h':
                    help_deconstruct(argv);
                    return 1;
                default:
                    abort();
            }
        }

        VG* graph;
        string file_name = argv[optind];
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
        //This struct contains the number of vertices in the graph
        //and the edges in a format that the superbubble package can understand.
        //Alright, superbubble time.
        for (auto i : graph->get_superbubbles()){
            cerr << i.first << " " << i.second << endl;
        }

          /* Find superbubbles */

        return 1;
    }
    void help_construct(char** argv) {
        cerr << "usage: " << argv[0] << " construct [options] >new.vg" << endl
             << "options:" << endl
             << "    -v, --vcf FILE        input VCF" << endl
             << "    -r, --reference FILE  input FASTA reference" << endl
             << "    -P, --ref-paths FILE  write reference paths in protobuf/gzip format to FILE" << endl
             << "    -B, --phase-blocks    save paths for phased blocks with the ref paths" << endl
             << "    -a, --alt-paths       save paths for alts of variants by variant ID" << endl
             << "    -R, --region REGION   specify a particular chromosome" << endl
             << "    -C, --region-is-chrom don't attempt to parse the region (use when the reference" << endl
             << "                          sequence name could be inadvertently parsed as a region)" << endl
             << "    -z, --region-size N   variants per region to parallelize" << endl
             << "    -m, --node-max N      limit the maximum allowable node sequence size" << endl
             << "                          nodes greater than this threshold will be divided" << endl
             << "    -p, --progress        show progress" << endl
             << "    -t, --threads N       use N threads to construct graph (defaults to numCPUs)" << endl
             << "    -f, --flat-alts N     don't chop up alternate alleles from input vcf" << endl;
    }

    int main_construct(int argc, char** argv) {

        if (argc == 2) {
            help_construct(argv);
            return 1;
        }

        string fasta_file_name, vcf_file_name, json_filename;
        string region;
        bool region_is_chrom = false;
        string output_type = "VG";
        bool progress = false;
        int vars_per_region = 25000;
        int max_node_size = 0;
        string ref_paths_file;
        bool flat_alts = false;
        // Should we make paths out of phasing blocks in the called samples?
        bool load_phasing_paths = false;
        // Should we make alt paths for variants?
        bool load_alt_paths = false;

        int c;
        while (true) {
            static struct option long_options[] =
                {
                    /* These options set a flag. */
                    //{"verbose", no_argument,       &verbose_flag, 1},
                    {"vcf", required_argument, 0, 'v'},
                    {"reference", required_argument, 0, 'r'},
                    // TODO: change the long option here?
                    {"ref-paths", required_argument, 0, 'P'},
                    {"phase-blocks", no_argument, 0, 'B'},
                    {"alt-paths", no_argument, 0, 'a'},
                    {"progress",  no_argument, 0, 'p'},
                    {"region-size", required_argument, 0, 'z'},
                    {"threads", required_argument, 0, 't'},
                    {"region", required_argument, 0, 'R'},
                    {"region-is-chrom", no_argument, 0, 'C'},
                    {"node-max", required_argument, 0, 'm'},\
                    {"flat-alts", no_argument, 0, 'f'},
                    {0, 0, 0, 0}
                };

            int option_index = 0;
            c = getopt_long (argc, argv, "v:r:phz:t:R:m:P:Bas:Cf",
                             long_options, &option_index);

            /* Detect the end of the options. */
            if (c == -1)
                break;

            switch (c)
            {
            case 'v':
                vcf_file_name = optarg;
                break;

            case 'r':
                fasta_file_name = optarg;
                break;

            case 'P':
                ref_paths_file = optarg;
                break;

            case 'B':
                load_phasing_paths = true;
                break;

            case 'a':
                load_alt_paths = true;
                break;

            case 'p':
                progress = true;
                break;

            case 'z':
                vars_per_region = atoi(optarg);
                break;

            case 'R':
                region = optarg;
                break;

            case 'C':
                region_is_chrom = true;
                break;

            case 't':
                omp_set_num_threads(atoi(optarg));
                break;

            case 'm':
                max_node_size = atoi(optarg);
                break;

            case 'f':
                flat_alts = true;
                break;

            case 'h':
            case '?':
                /* getopt_long already printed an error message. */
                help_construct(argv);
                exit(1);
                break;

            default:
                abort ();
            }
        }

        if(load_phasing_paths && ref_paths_file.empty()) {
            cerr << "error:[vg construct] cannot save phasing paths without a paths file name" << endl;
            return 1;
        }

        vcflib::VariantCallFile variant_file;
        if (!vcf_file_name.empty()) {
            variant_file.open(vcf_file_name);
            if (!variant_file.is_open()) {
                cerr << "error:[vg construct] could not open" << vcf_file_name << endl;
                return 1;
            }
        }

        FastaReference reference;
        if (fasta_file_name.empty()) {
            cerr << "error:[vg construct] a reference is required for graph construction" << endl;
            return 1;
        }
        reference.open(fasta_file_name);

        // store our reference sequence paths
        // TODO: use this. Maybe dump paths here instead of in the graph?
        Paths ref_paths;

        VG graph(variant_file, reference, region, region_is_chrom, vars_per_region,
            max_node_size, flat_alts, load_phasing_paths, load_alt_paths, progress);

        if (!ref_paths_file.empty()) {
            ofstream paths_out(ref_paths_file);
            graph.paths.write(paths_out);
            if(load_phasing_paths) {
                // Keep only the non-phasing paths in the graph. If you keep too
                // many paths in a graph, you'll make chunks that are too large.
                // TODO: dynamically deliniate the chunks in the serializer so you
                // won't write vg files you can't read.

                set<string> non_phase_paths;
                string phase_prefix = "_phase";
                graph.paths.for_each_name([&](string path_name) {
                    if(!equal(phase_prefix.begin(), phase_prefix.end(), path_name.begin())) {
                        // Path is not a phase path
                        non_phase_paths.insert(path_name);
                    }
                });

                // Keep only the non-phase paths
                graph.paths.keep_paths(non_phase_paths);
            }
        }

        graph.serialize_to_ostream(std::cout);

        // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
        // this would free all the memory used by protobuf:
        //ShutdownProtobufLibrary();

        return 0;
    }

    void vg_help(char** argv) {
        cerr << "usage: " << argv[0] << " <command> [options]" << endl
             << endl
             << "commands:" << endl
             << "  -- construct     graph construction" << endl
             << "  -- deconstruct   convert a graph into VCF relative to a reference." << endl
             << "  -- view          format conversions for graphs and alignments" << endl
             << "  -- vectorize     Transform alignments to one-hot vectors." << endl
             << "  -- index         index features of the graph in a disk-backed key/value store" << endl
             << "  -- find          use an index to find nodes, edges, kmers, or positions" << endl
             << "  -- paths         traverse paths in the graph" << endl
             << "  -- align         local alignment" << endl
             << "  -- map           global alignment" << endl
             << "  -- stats         metrics describing graph properties" << endl
             << "  -- join          combine graphs via a new head" << endl
             << "  -- ids           manipulate node ids" << endl
             << "  -- concat        concatenate graphs tail-to-head" << endl
             << "  -- kmers         enumerate kmers of the graph" << endl
             << "  -- sim           simulate reads from the graph" << endl
             << "  -- mod           filter, transform, and edit the graph" << endl
             << "  -- surject       map alignments onto specific paths" << endl
             << "  -- msga          multiple sequence graph alignment" << endl
             << "  -- pileup        build a pileup from a set of alignments" << endl
             << "  -- call          prune the graph by genotyping a pileup" << endl
             << "  -- compare       compare the kmer space of two graphs" << endl
             << "  -- validate      validate the semantics of a graph" << endl;
    }

    int main(int argc, char *argv[])
    {

        if (argc == 1) {
            vg_help(argv);
            return 1;
        }

        //omp_set_dynamic(1); // use dynamic scheduling

        string command = argv[1];
        if (command == "construct") {
            return main_construct(argc, argv);
        } else if (command == "deconstruct"){
            return main_deconstruct(argc, argv);
        } else if (command == "view") {
            return main_view(argc, argv);
        } else if (command == "align") {
            return main_align(argc, argv);
        } else if (command == "map") {
            return main_map(argc, argv);
        } else if (command == "index") {
            return main_index(argc, argv);
        } else if (command == "find") {
            return main_find(argc, argv);
        } else if (command == "paths") {
            return main_paths(argc, argv);
        } else if (command == "stats") {
            return main_stats(argc, argv);
        } else if (command == "join") {
            return main_join(argc, argv);
        } else if (command == "ids") {
            return main_ids(argc, argv);
        } else if (command == "concat") {
            return main_concat(argc, argv);
        } else if (command == "kmers") {
            return main_kmers(argc, argv);
        } else if (command == "sim") {
            return main_sim(argc, argv);
        } else if (command == "mod") {
            return main_mod(argc, argv);
        } else if (command == "surject") {
            return main_surject(argc, argv);
        } else if (command == "msga") {
            return main_msga(argc, argv);
        } else if (command == "pileup") {
            return main_pileup(argc, argv);
        } else if (command == "call") {
            return main_call(argc, argv);
        } else if (command == "compare") {
            return main_compare(argc, argv);
        } else if (command == "validate") {
            return main_validate(argc, argv);
        } else if (command == "vectorize") {
            return main_vectorize(argc, argv);
        }else {
            cerr << "error:[vg] command " << command << " not found" << endl;
            vg_help(argv);
            return 1;
        }

        return 0;

    }
