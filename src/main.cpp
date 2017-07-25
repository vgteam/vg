#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdio>
#include <getopt.h>
#include <sys/stat.h>
#include "gcsa/gcsa.h"
#include "gcsa/algorithms.h"
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
#include "filter.hpp"
#include "google/protobuf/stubs/common.h"
#include "progress_bar.hpp"
#include "version.hpp"
#include "genotyper.hpp"
#include "bubbles.hpp"
#include "translator.hpp"
#include "readfilter.hpp"
#include "distributions.hpp"
#include "unittest/driver.hpp"
// New subcommand system provides all the subcommands that used to live here
#include "subcommand/subcommand.hpp"
#include "flow_sort.hpp"


using namespace std;
using namespace google::protobuf;
using namespace vg;

void help_translate(char** argv) {
    cerr << "usage: " << argv[0] << " translate [options] translation" << endl
         << "Translate alignments or paths using the translation map." << endl
         << endl
         << "options:" << endl
         << "    -p, --paths FILE      project the input paths into the from-graph" << endl
         << "    -a, --alns FILE       project the input alignments into the from-graph" << endl
         << "    -l, --loci FILE       project the input locus descriptions into the from-graph" << endl
         << "    -m, --mapping JSON    print the from-mapping corresponding to the given JSON mapping" << endl
         << "    -P, --position JSON   print the from-position corresponding to the given JSON position" << endl
         << "    -o, --overlay FILE    overlay this translation on top of the one we are given" << endl;
}

int main_translate(int argc, char** argv) {

    if (argc <= 2) {
        help_translate(argv);
        return 1;
    }

    string position_string;
    string mapping_string;
    string path_file;
    string aln_file;
    string loci_file;
    string overlay_file;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"position", required_argument, 0, 'P'},
            {"mapping", required_argument, 0, 'm'},
            {"paths", required_argument, 0, 'p'},
            {"alns", required_argument, 0, 'a'},
            {"loci", required_argument, 0, 'l'},
            {"overlay", required_argument, 0, 'o'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp:m:P:a:o:l:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'P':
            position_string = optarg;
            break;

        case 'm':
            mapping_string = optarg;
            break;

        case 'p':
            path_file = optarg;
            break;

        case 'a':
            aln_file = optarg;
            break;

        case 'l':
            loci_file = optarg;
            break;

        case 'o':
            overlay_file = optarg;
            break;

        case 'h':
        case '?':
            help_translate(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    Translator* translator;
    get_input_file(optind, argc, argv, [&](istream& in) {
        translator = new Translator(in);
    });

    // test the position translation
    if (!position_string.empty()) {
        Position position;
        json2pb(position, position_string.c_str(), position_string.size());
        cout << pb2json(translator->translate(position)) << endl;
    }

    // test the mapping translation
    if (!mapping_string.empty()) {
        Mapping mapping;
        json2pb(mapping, mapping_string.c_str(), mapping_string.size());
        cout << pb2json(translator->translate(mapping)) << endl;
    }

    if (!path_file.empty()) {
        vector<Path> buffer;
        function<void(Path&)> lambda = [&](Path& path) {
            buffer.push_back(translator->translate(path));
            stream::write_buffered(cout, buffer, 100);
        };
        ifstream path_in(path_file);
        stream::for_each(path_in, lambda);
        stream::write_buffered(cout, buffer, 0);
    } else if (!aln_file.empty()) {
        vector<Alignment> buffer;
        function<void(Alignment&)> lambda = [&](Alignment& aln) {
            buffer.push_back(translator->translate(aln));
            stream::write_buffered(cout, buffer, 100);
        };
        ifstream aln_in(aln_file);
        stream::for_each(aln_in, lambda);
        stream::write_buffered(cout, buffer, 0);
    } else if (!loci_file.empty()) {
        vector<Locus> buffer;
        function<void(Locus&)> lambda = [&](Locus& locus) {
            buffer.push_back(translator->translate(locus));
            stream::write_buffered(cout, buffer, 100);
        };
        ifstream loci_in(loci_file);
        stream::for_each(loci_in, lambda);
        stream::write_buffered(cout, buffer, 0);
    }

    if (!overlay_file.empty()) {
        vector<Translation> buffer;
        function<void(Translation&)> lambda = [&](Translation& trans) {
            buffer.push_back(translator->overlay(trans));
            stream::write_buffered(cout, buffer, 100);
        };
        ifstream overlay_in(overlay_file);
        stream::for_each(overlay_in, lambda);
        stream::write_buffered(cout, buffer, 0);
    }

    return 0;
}

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
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

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

    if (db_name1.empty()) {
        db_name1 = get_input_file_name(optind, argc, argv);
    }
    if (db_name2.empty()) {
        db_name2 = get_input_file_name(optind, argc, argv);
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


void help_pileup(char** argv) {
    cerr << "usage: " << argv[0] << " pileup [options] <graph.vg> <alignment.gam> > out.vgpu" << endl
         << "Calculate pileup for each position in graph and output in VG Pileup format (list of protobuf NodePileups)." << endl
         << endl
         << "options:" << endl
         << "    -j, --json              output in JSON" << endl
         << "    -q, --min-quality N     ignore bases with PHRED quality < N (default=10)" << endl
         << "    -m, --max-mismatches N  ignore bases with > N mismatches within window centered on read (default=1)" << endl
         << "    -w, --window-size N     size of window to apply -m option (default=0)" << endl
         << "    -d, --max-depth N       maximum depth pileup to create (further maps ignored) (default=1000)" << endl
         << "    -M, --ignore-mapq       do not combine mapping qualities with base qualities" << endl
         << "    -p, --progress          show progress" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << "    -v, --verbose           print stats on bases filtered" << endl;
}

int main_pileup(int argc, char** argv) {

    if (argc <= 3) {
        help_pileup(argv);
        return 1;
    }

    bool output_json = false;
    bool show_progress = false;
    int thread_count = 1;
    int min_quality = 10;
    int max_mismatches = 1;
    int window_size = 0;
    int max_depth = 1000; // used to prevent protobuf messages getting to big
    bool verbose = false;
    bool use_mapq = true;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"json", required_argument, 0, 'j'},
                {"min-quality", required_argument, 0, 'q'},
                {"max-mismatches", required_argument, 0, 'm'},
                {"window-size", required_argument, 0, 'w'},
                {"progress", required_argument, 0, 'p'},
                {"max-depth", required_argument, 0, 'd'},
                {"ignore-mapq", no_argument, 0, 'M'},
                {"threads", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "jq:m:w:pd:at:v",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'j':
            output_json = true;
            break;
        case 'q':
            min_quality = atoi(optarg);
            break;
        case 'm':
            max_mismatches = atoi(optarg);
            break;
        case 'w':
            window_size = atoi(optarg);
            break;
        case 'd':
            max_depth = atoi(optarg);
            break;
        case 'M':
            use_mapq = false;
            break;
        case 'p':
            show_progress = true;
            break;
        case 't':
            thread_count = atoi(optarg);
            break;
        case 'v':
            verbose = true;
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

    // Parse the arguments
    if (optind >= argc) {
        help_pileup(argv);
        return 1;
    }
    string graph_file_name = get_input_file_name(optind, argc, argv);
    if (optind >= argc) {
        help_pileup(argv);
        return 1;
    }
    string alignments_file_name = get_input_file_name(optind, argc, argv);
    
    if (alignments_file_name == "-" && graph_file_name == "-") {
        cerr << "error: graph and alignments can't both be from stdin." << endl;
        exit(1);
    }

    // read the graph
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }
    VG* graph;
    get_input_file(graph_file_name, [&](istream& in) {
        graph = new VG(in);
    });

    // Make Pileups makers for each thread.
    vector<Pileups> pileups(thread_count, Pileups(graph, min_quality, max_mismatches, window_size, max_depth, use_mapq));
    
    // setup alignment stream
    get_input_file(alignments_file_name, [&](istream& alignment_stream) {
        // compute the pileups.
        if (show_progress) {
            cerr << "Computing pileups" << endl;
        }
        
        function<void(Alignment&)> lambda = [&pileups, &graph](Alignment& aln) {
            int tid = omp_get_thread_num();
            pileups[tid].compute_from_alignment(aln);
        };
        stream::for_each_parallel(alignment_stream, lambda);
    });

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

    // number of bases filtered
    if (verbose) {
        cerr << "Bases filtered by min. quality: " << pileups[0]._min_quality_count << endl
             << "Bases filtered by max mismatch: " << pileups[0]._max_mismatch_count << endl
             << "Total bases:                    " << pileups[0]._bases_count << endl << endl;
    }

    return 0;
}

void help_circularize(char** argv){
    cerr << "usage: " << argv[0] << " circularize [options] <graph.vg> > [circularized.vg]" << endl
        << "Makes specific paths or nodes in a graph circular." << endl
        << endl
        << "options:" << endl
        << "    -p  --path  <PATHNAME>  circularize the path by connecting its head/tail node." << endl
        << "    -P, --pathfile <PATHSFILE> circularize all paths in the provided file." << endl
        << "    -a, --head  <node_id>   circularize a head and tail node (must provide a tail)." << endl
        << "    -z, --tail  <tail_id>   circularize a head and tail node (must provide a head)." << endl
        << "    -d  --describe          list all the paths in the graph."   << endl
        << endl;
    exit(1);
}

int main_circularize(int argc, char** argv){
    if (argc == 2){
        help_circularize(argv);
        exit(1);
    }

    string path = "";
    string pathfile = "";
    bool describe = false;
    vg::id_t head = -1;
    vg::id_t tail = -1;


    int c;
    optind = 2;
    while (true){
        static struct option long_options[] =
        {
            {"path", required_argument, 0, 'p'},
            {"pathfile", required_argument, 0, 'P'},
            {"head", required_argument, 0, 'a'},
            {"tail", required_argument, 0, 'z'},
            {"describe", required_argument, 0, 'd'},
            {0,0,0,0}
        };


    int option_index = 0;
    c = getopt_long (argc, argv, "hdp:P:a:z:",
            long_options, &option_index);
    if (c == -1){
        break;
    }

        switch(c){
            case 'a':
                head = atoi(optarg);
                break;
            case 'z':
                tail = atoi(optarg);
                break;
            case 'p':
                path = optarg;
                break;
            case 'P':
                pathfile = optarg;
                break;
            case 'd':
                describe = true;
                break;
            case 'h':
            case '?':
                help_circularize(argv);
                exit(1);
                break;

            default:
                abort();
        }
    }

    vector<string> paths_to_circularize;
    if (!((head * tail) > 0)){
        cerr << "Both a head and tail node must be provided" << endl;
        help_circularize(argv);
        exit(1);
    }
    if  (pathfile != ""){
        string line;
        ifstream pfi;
        pfi.open(pathfile);
        if (!pfi.good()){
            cerr << "There is an error with the input file." << endl;
            help_circularize(argv);
        }
        while (getline(pfi, line)){
            paths_to_circularize.push_back(line);
        }
        pfi.close();

    }
    else if (path != ""){
        paths_to_circularize.push_back(path);
    }

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    // Check if paths are in graph:
    for (string p : paths_to_circularize){
        bool paths_in_graph = true;
        if (!graph->paths.has_path(p)){
            cerr << "ERROR: PATH NOT IN GRAPH - " << p << endl;
            paths_in_graph = false;
        }

        if (!paths_in_graph){
            exit(1);
        }

    }

    if (describe){
       for (pair<string, list<Mapping> > p : graph->paths._paths){
            cout << p.first << endl;
       }
       exit(0);
    }

    if (head > 0 && tail > head){
        graph->circularize(head, tail);
    }
    else{
        graph->circularize(paths_to_circularize);
    }

    graph->serialize_to_ostream(std::cout);
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
        string file_name = get_input_file_name(optind, argc, argv);
        graph_file_names.push_back(file_name);
    }

    VGset graphs(graph_file_names);

    graphs.show_progress = show_progress;

    if (gcsa_out) {
        if (edge_max != 0) {
            // I have been passing this option to vg index -g for months
            // thinking it worked. But it can't work. So we should tell the user
            // they're wrong.
            cerr << "error:[vg kmers] Cannot limit edge crossing (-e) when generating GCSA kmers (-g)."
                << " Use vg mod -p to prune the graph instead." << endl;
            exit(1);
        }
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
        get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });
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
        get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });

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
            string file_name = get_input_file_name(optind, argc, argv);
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
        get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });
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

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options] <graph.vg>" << endl
        << "options:" << endl
        << "  inspection:" << endl
        << "    -x, --extract         return (as GAM alignments) the stored paths in the graph" << endl
        << "    -L, --list            return (as a list of names, one per line) the path names" << endl
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
    bool list_paths = false;
    bool path_only = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            {"extract", no_argument, 0, 'x'},
            {"list", no_argument, 0, 'L'},
            {"node", required_argument, 0, 'n'},
            {"max-length", required_argument, 0, 'l'},
            {"edge-max", required_argument, 0, 'e'},
            {"as-seqs", no_argument, 0, 's'},
            {"path-only", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "n:l:hse:xLp",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

            case 'x':
                extract = true;
                break;
                
            case 'L':
                list_paths = true;
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
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    if (extract) {
        vector<Alignment> alns = graph->paths_as_alignments();
        write_alignments(cout, alns);
        delete graph;
        return 0;
    }
    
    if (list_paths) {
        graph->paths.for_each_name([&](const string& name) {
            cout << name << endl;
        });
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

    void help_sv(char** argv){
        cerr << "usage: " << argv[0] << " sv [options] <aln.gam>" << endl
            << "options: " << endl
            << " -g --graph <graph>.vg " << endl
            << " -m --mask <vcf>.vcf" << endl
            << endl;
    }

void help_locify(char** argv){
    cerr << "usage: " << argv[0] << " locify [options] " << endl
         << "    -l, --loci FILE      input loci over which to locify the alignments" << endl
         << "    -a, --aln-idx DIR    use this rocksdb alignment index (from vg index -N)" << endl
         << "    -x, --xg-idx FILE    use this xg index" << endl
         << "    -n, --name-alleles   generate names for each allele rather than using full Paths" << endl
         << "    -f, --forwardize     flip alignments on the reverse strand to the forward" << endl
         << "    -s, --sorted-loci FILE  write the non-nested loci out in their sorted order" << endl
         << "    -b, --n-best N       keep only the N-best alleles by alignment support" << endl
         << "    -o, --out-loci FILE  rewrite the loci with only N-best alleles kept" << endl;
        // TODO -- add some basic filters that are useful downstream in whatshap
}

int main_locify(int argc, char** argv){
    string gam_idx_name;
    string loci_file;
    Index gam_idx;
    string xg_idx_name;
    bool name_alleles = false;
    bool forwardize = false;
    string loci_out, sorted_loci;
    int n_best = 0;

    if (argc <= 2){
        help_locify(argv);
        exit(1);
    }

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"gam-idx", required_argument, 0, 'g'},
            {"loci", required_argument, 0, 'l'},
            {"xg-idx", required_argument, 0, 'x'},
            {"name-alleles", no_argument, 0, 'n'},
            {"forwardize", no_argument, 0, 'f'},
            {"sorted-loci", required_argument, 0, 's'},
            {"loci-out", required_argument, 0, 'o'},
            {"n-best", required_argument, 0, 'b'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:x:g:nfo:b:s:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            gam_idx_name = optarg;
            break;

        case 'l':
            loci_file = optarg;
            break;

        case 'x':
            xg_idx_name = optarg;
            break;

        case 'n':
            name_alleles = true;
            break;

        case 'f':
            forwardize = true;
            break;

        case 'o':
            loci_out = optarg;
            break;

        case 's':
            sorted_loci = optarg;
            break;

        case 'b':
            n_best = atoi(optarg);
            name_alleles = true;
            break;

        case 'h':
        case '?':
            help_locify(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    if (!gam_idx_name.empty()) {
        gam_idx.open_read_only(gam_idx_name);
    }

    if (xg_idx_name.empty()) {
        cerr << "[vg locify] Error: no xg index provided" << endl;
        return 1;
    }
    ifstream xgstream(xg_idx_name);
    xg::XG xgidx(xgstream);

    std::function<vector<string>(string, char)> strsplit = [&](string x, char delim){

        vector<string> ret;
        stringstream ss;
        std::string tok;
        while (getline(ss, tok, delim)){
            ret.push_back(tok);
        }
        return ret;

    };

    vector<string> locus_names;
    map<string, map<string, int > > locus_allele_names;
    map<string, Alignment> alignments_with_loci;
    map<pos_t, set<string> > pos_to_loci;
    map<string, set<pos_t> > locus_to_pos;
    map<string, map<int, int> > locus_allele_support;
    map<string, vector<int> > locus_to_best_n_alleles;
    map<string, set<int> > locus_to_keep;
    int count = 0;

    std::function<void(Locus&)> lambda = [&](Locus& l){
        locus_names.push_back(l.name());
        set<vg::id_t> nodes_in_locus;
        for (int i = 0; i < l.allele_size(); ++i) {
            auto& allele = l.allele(i);
            for (int j = 0; j < allele.mapping_size(); ++j) {
                auto& position = allele.mapping(j).position();
                nodes_in_locus.insert(position.node_id());
            }
            // for position in mapping
            map<pos_t, int> ref_positions;
            map<pos_t, Edit> edits;
            decompose(allele, ref_positions, edits);
            // warning: uses only reference positions!!!
            for (auto& pos : ref_positions) {
                pos_to_loci[pos.first].insert(l.name());
                locus_to_pos[l.name()].insert(pos.first);
            }
        }
        // void for_alignment_in_range(int64_t id1, int64_t id2, std::function<void(const Alignment&)> lambda);
        std::function<void(const Alignment&)> fill_alns = [&](const Alignment& a){
            // TODO reverse complementing alleles ?
            // overlap is stranded
            //matching
            // find the most-matching allele
            map<double, vector<int> > matches;
            for (int i = 0; i < l.allele_size(); ++i) {
                auto& allele = l.allele(i);
                matches[overlap(a.path(), allele)].push_back(i);
            }
            assert(l.allele_size());
            int best = matches.rbegin()->second.front();
            Locus matching;
            matching.set_name(l.name());
            if (name_alleles) {
                //map<string, map<string, int > > locus_allele_names;
                auto& allele = l.allele(best);
                string s;
                allele.SerializeToString(&s);
                auto& l_names = locus_allele_names[l.name()];
                auto f = l_names.find(s);
                int name_int = 0;
                if (f == l_names.end()) {
                    int next_id = l_names.size() + 1;
                    l_names[s] = next_id;
                    name_int = next_id;
                } else {
                    name_int = f->second;
                }
                string allele_name = vg::convert(name_int);
                Path p;
                p.set_name(allele_name);
                *matching.add_allele() = p;
                if (n_best) {
                    // record support for this allele
                    // we'll use to filter the locus records later
                    locus_allele_support[l.name()][name_int]++;
                }
            } else {
                *matching.add_allele() = l.allele(best);
                // TODO get quality score relative to this specific allele / alignment
                // record in the alignment we'll save
            }
            if (alignments_with_loci.find(a.name()) == alignments_with_loci.end()) {
                alignments_with_loci[a.name()] = a;
            }
            Alignment& aln = alignments_with_loci[a.name()];
            *aln.add_locus() = matching;
        };
        vector<vg::id_t> nodes_vec;
        for (auto& id : nodes_in_locus) nodes_vec.push_back(id);
        gam_idx.for_alignment_to_nodes(nodes_vec, fill_alns);
    };

    if (!loci_file.empty()){
        ifstream ifi(loci_file);
        stream::for_each(ifi, lambda);
    } else {
        cerr << "[vg locify] Warning: empty locus file given, could not annotate alignments with loci." << endl;
    }

    // find the non-nested loci
    vector<string> non_nested_loci;
    for (auto& name : locus_names) {
        // is it nested?
        auto& positions = locus_to_pos[name];
        int min_loci = 0;
        for (auto& pos : positions) {
            auto& loci = pos_to_loci[pos];
            min_loci = (min_loci == 0 ? (int)loci.size() : min(min_loci, (int)loci.size()));
        }
        if (min_loci == 1) {
            // not fully contained in any other locus
            non_nested_loci.push_back(name);
        }
    }

    // filter out the non-best alleles
    if (n_best) {
        // find the n-best
        for (auto& supp : locus_allele_support) {
            auto& name = supp.first;
            auto& alleles = supp.second;
            map<int, int> ranked;
            for (auto& allele : alleles) {
                ranked[allele.second] = allele.first;
            }
            auto& to_keep = locus_to_keep[name];
            for (auto r = ranked.rbegin(); r != ranked.rend(); ++r) {
                to_keep.insert(r->second);
                if (to_keep.size() == n_best) {
                    break;
                }
            }
        }
        // filter out non-n-best from the alignments
        for (auto& a : alignments_with_loci) {
            auto& aln = a.second;
            vector<Locus> kept;
            for (int i = 0; i < aln.locus_size(); ++i) {
                auto& allele = aln.locus(i).allele(0);
                if (locus_to_keep[aln.locus(i).name()].count(atoi(allele.name().c_str()))) {
                    kept.push_back(aln.locus(i));
                }
            }
            aln.clear_locus();
            for (auto& l : kept) {
                *aln.add_locus() = l;
            }
        }
    }

    if (n_best && !loci_out.empty()) {
        // filter out non-n-best from the loci
        if (!loci_file.empty()){
            ofstream outloci(loci_out);
            vector<Locus> buffer;
            std::function<void(Locus&)> lambda = [&](Locus& l){
                // remove the alleles which are to filter
                //map<string, map<string, int > > locus_allele_names;
                auto& allele_names = locus_allele_names[l.name()];
                auto& to_keep = locus_to_keep[l.name()];
                vector<Path> alleles_to_keep;
                for (int i = 0; i < l.allele_size(); ++i) {
                    auto allele = l.allele(i);
                    string s; allele.SerializeToString(&s);
                    auto& name = allele_names[s];
                    if (to_keep.count(name)) {
                        allele.set_name(vg::convert(name));
                        alleles_to_keep.push_back(allele);
                    }
                }
                l.clear_allele();
                for (auto& allele : alleles_to_keep) {
                    *l.add_allele() = allele;
                }
                buffer.push_back(l);
                stream::write_buffered(outloci, buffer, 100);
            };
            ifstream ifi(loci_file);
            stream::for_each(ifi, lambda);
            stream::write_buffered(outloci, buffer, 0);
            outloci.close();
        } else {
            cerr << "[vg locify] Warning: empty locus file given, could not update loci." << endl;
        }
    }

    // sort them using... ? ids?
    sort(non_nested_loci.begin(), non_nested_loci.end(),
         [&locus_to_pos](const string& s1, const string& s2) {
             return *locus_to_pos[s1].begin() < *locus_to_pos[s2].begin();
         });

    if (!sorted_loci.empty()) {
        ofstream outsorted(sorted_loci);
        for (auto& name : non_nested_loci) {
            outsorted << name << endl;
        }
        outsorted.close();
    }

    vector<Alignment> output_buf;
    for (auto& aln : alignments_with_loci) {
        // TODO order the loci by their order in the alignments
        if (forwardize) {
            if (aln.second.path().mapping_size() && aln.second.path().mapping(0).position().is_reverse()) {
                output_buf.push_back(reverse_complement_alignment(aln.second,
                                                                  [&xgidx](int64_t id) { return xgidx.node_length(id); }));
            } else {
                output_buf.push_back(aln.second);
            }
        } else {
            output_buf.push_back(aln.second);
        }
        stream::write_buffered(cout, output_buf, 100);
    }
    stream::write_buffered(cout, output_buf, 0);        
    
    return 0;
}

void help_deconstruct(char** argv){
    cerr << "usage: " << argv[0] << " deconstruct [options] -p <PATH> <my_graph>.vg" << endl
         << "Outputs VCF records for Snarls present in a graph (relative to a chosen reference path)." << endl
         << "options: " << endl
         << "--path / -p     REQUIRED: A reference path to deconstruct against." << endl
         << endl;
}

int main_deconstruct(int argc, char** argv){
    //cerr << "WARNING: EXPERIMENTAL" << endl;
    if (argc <= 2) {
        help_deconstruct(argv);
        return 1;
    }

    vector<string> refpaths;
    string graphname;
    string outfile = "";
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"path", required_argument, 0, 'p'},
                {0, 0, 0, 0}

            };

            int option_index = 0;
            c = getopt_long (argc, argv, "hp:",
                    long_options, &option_index);

            // Detect the end of the options.
            if (c == -1)
                break;

            switch (c)
            {
                case 'p':
                    refpaths = split(optarg, ",");
                    break;
                case '?':
                case 'h':
                    help_deconstruct(argv);
                    return 1;
                default:
                    help_deconstruct(argv);
                    abort();
            }

        }
        graphname = argv[optind];
        vg::VG* graph;
        if (!graphname.empty()){
            ifstream gstream(graphname);
            graph = new vg::VG(gstream);
        }

        // load graph

        // Deconstruct
        Deconstructor dd;
        dd.deconstruct(refpaths, graph);
    return 0;
}


void help_sort(char** argv){
    cerr << "usage: " << argv[0] << " sort [options] -i <input_file> -r <reference_name> > sorted.vg " << endl
         << "options: " << endl
         << "           -g, --gfa              input in GFA format" << endl
         << "           -i, --in               input file" << endl
         << "           -r, --ref              reference name" << endl
         << "           -w, --without-grooming no grooming mode" << endl
         << "           -f, --fast             sort using Eades algorithm, otherwise max-flow sorting is used" << endl   
         << endl;
}

int main_sort(int argc, char *argv[]) {

    //default input format is vg
    bool gfa_input = false;
    string file_name = "";
    string reference_name = "";
    bool without_grooming = false;
    bool use_fast_algorithm = false;
    int c;
    while (true) {
        static struct option long_options[] =
            {
                {"gfa", no_argument, 0, 'g'},
                {"in", required_argument, 0, 'i'},
                {"ref", required_argument, 0, 'r'},
                {"without-grooming", no_argument, 0, 'w'},
                {"fast", no_argument, 0, 'f'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "i:r:gwf",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'g':
            gfa_input = true;
            break;
        case 'r':
            reference_name = optarg;
            break;
        case 'i':
            file_name = optarg;
            break;
        case 'w':
            without_grooming = true;
            break;
        case 'f':
            use_fast_algorithm = true;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_sort(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }
  
    if (reference_name.empty() || file_name.empty()) {
        help_sort(argv);
        exit(1);
    }
    
    ifstream in;
    std::unique_ptr<VG> graph;
    {
        in.open(file_name.c_str());        
        if (gfa_input) {
            graph.reset(new VG());
            graph->from_gfa(in);
        } else {
            graph.reset(new VG(in));
        }
    }
    FlowSort flow_sort(*graph.get());
    if (use_fast_algorithm) {
        flow_sort.fast_linear_sort(reference_name, !without_grooming);
    } else {
        flow_sort.max_flow_sort(reference_name);
    }
    
    graph->serialize_to_ostream(std::cout);
    in.close();
    return 0;
}

void help_version(char** argv){
    cerr << "usage: " << argv[0] << " version" << endl
         << "options: " << endl
         << endl;
}

int main_version(int argc, char** argv){

    if (argc != 2) {
        help_version(argv);
        return 1;
    }

    cout << VG_VERSION_STRING << endl;
    return 0;
}

// No help_test is necessary because the unit testing library takes care of
// complaining about missing options.

int main_test(int argc, char** argv){
    // Forward arguments along to the main unit test driver
    return vg::unittest::run_unit_tests(argc, argv);
}

void vg_help(char** argv) {
    cerr << "vg: variation graph tool, version " << VG_VERSION_STRING << endl
         << endl
         << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << "commands:" << endl;
         
     vg::subcommand::Subcommand::for_each([](const vg::subcommand::Subcommand& command) {
        // Announce every subcommand we have
        
        // Pad all the names so the descriptions line up
        string name = command.get_name();
        name.resize(14, ' ');
        cerr << "  -- " << name << command.get_description() << endl;
     });
         
     // Also announce all the old-style hardcoded commands
     cerr << "  -- deconstruct   convert a graph into VCF relative to a reference." << endl
         << "  -- paths         traverse paths in the graph" << endl
         << "  -- join          combine graphs via a new head" << endl
         << "  -- ids           manipulate node ids" << endl
         << "  -- concat        concatenate graphs tail-to-head" << endl
         << "  -- kmers         enumerate kmers of the graph" << endl
         << "  -- pileup        build a pileup from a set of alignments" << endl
         << "  -- genotype      compute genotypes from aligned reads" << endl
         << "  -- compare       compare the kmer space of two graphs" << endl
         << "  -- circularize   circularize a path within a graph." << endl
         << "  -- translate     project alignments and paths through a graph translation" << endl
         << "  -- validate      validate the semantics of a graph" << endl
         << "  -- sort          sort variant graph using max flow algorithm or Eades fast heuristic algorithm" << endl
         << "  -- test          run unit tests" << endl
         << "  -- version       version information" << endl;
}

int main(int argc, char *argv[])
{

    // set a higher value for tcmalloc warnings
    setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);

    if (argc == 1) {
        vg_help(argv);
        return 1;
    }
    
    auto* subcommand = vg::subcommand::Subcommand::get(argc, argv);
    if (subcommand != nullptr) {
        // We found a matching subcommand, so run it
        return (*subcommand)(argc, argv);
    }
    
    // Otherwise, fall abck on the old chain of if statements.

    //omp_set_dynamic(1); // use dynamic scheduling

    string command = argv[1];
    if (command == "deconstruct"){
        return main_deconstruct(argc, argv);
    } else if (command == "paths") {
        return main_paths(argc, argv);
    } else if (command == "join") {
        return main_join(argc, argv);
    } else if (command == "ids") {
        return main_ids(argc, argv);
    } else if (command == "concat") {
        return main_concat(argc, argv);
    } else if (command == "kmers") {
        return main_kmers(argc, argv);
    } else if (command == "pileup") {
        return main_pileup(argc, argv);
    } else if (command == "compare") {
        return main_compare(argc, argv);
    } else if (command == "validate") {
        return main_validate(argc, argv);
    } else if (command == "circularize"){
        return main_circularize(argc, argv);
    }  else if (command == "translate") {
        return main_translate(argc, argv);
    }  else if (command == "version") {
        return main_version(argc, argv);
    } else if (command == "test") {
        return main_test(argc, argv);
    } else if (command == "locify"){
        return main_locify(argc, argv);
    } else if (command == "sort") {
        return main_sort(argc, argv);
    } else {
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

    return 0;

}
