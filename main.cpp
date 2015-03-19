#include <iostream>
#include <fstream>
#include <ctime>
#include <getopt.h>
#include "pb2json.h"
#include "vg.hpp"
#include "vg.pb.h"
#include "vg_set.hpp"
#include "index.hpp"
#include "mapper.hpp"
#include "Variant.h"
#include "Fasta.h"
#include "stream.hpp"
#include "alignment.hpp"
#include <google/protobuf/stubs/common.h>

using namespace std;
using namespace google::protobuf;
using namespace vg;

void help_sim(char** argv) {
    cerr << "usage: " << argv[0] << " sim [options] <graph.vg>" << endl
         << "Simulates reads from the graph(s). Output is a list of reads." << endl
         << endl
         << "options:" << endl
         << "    -l, --read-length N   write reads of length N" << endl
         << "    -n, --num-reads N     simulate N reads" << endl
         << "    -s, --random-seed N   use this specific seed for the PRNG" << endl;
}

int main_sim(int argc, char** argv) {

    if (argc == 2) {
        help_sim(argv);
        return 1;
    }

    int read_length = 100;
    int num_reads = 1;
    int seed_val = time(NULL);

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"read-length", required_argument, 0, 'l'},
                {"num-reads", required_argument, 0, 'n'},
                {"random-seed", required_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:",
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

    for (int i = 0; i < num_reads; ++i) {
        string readseq = graph->random_read(read_length, rng, min_id, max_id, true);
        // avoid short reads at the end of the graph by retrying
        while (readseq.size() < read_length) {
            readseq = graph->random_read(read_length, rng, min_id, max_id, true);
        }
        cout << readseq << endl;
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
         << "    -e, --edge-max N     cross no more than N edges when determining k-paths" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers in paths (default 1)" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -p, --progress        show progress" << endl;
}

int main_kmers(int argc, char** argv) {

    if (argc == 2) {
        help_kmers(argv);
        return 1;
    }

    int kmer_size = 0;
    int edge_max = 0;
    int kmer_stride = 1;
    bool show_progress = false;

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
                {"progress",  no_argument, 0, 'p'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hk:j:pt:e:",
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

        case 'p':
            show_progress = true;
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

    if (edge_max == 0) edge_max = kmer_size + 1;

    vector<string> graph_file_names;
    while (optind < argc) {
        string file_name = argv[optind++];
        graph_file_names.push_back(file_name);
    }

    VGset graphs(graph_file_names);
    function<void(string&, Node*, int)>
        lambda = [](string& kmer, Node* n, int p) {
#pragma omp critical (cout)
        cout << kmer << '\t' << n->id() << '\t' << p << '\n';
    };
    graphs.show_progress = show_progress;
    graphs.for_each_kmer_parallel(lambda, kmer_size, edge_max, kmer_stride);
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
         << "                         their ids to be non-conflicting" << endl;
}

int main_ids(int argc, char** argv) {

    if (argc == 2) {
        help_ids(argv);
        return 1;
    }

    bool join = false;
    bool compact = false;
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
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hci:d:j",
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

        if (compact) {
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
        joined.extend(**g);
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
         << "    -l, --length          length of sequences in graph" << endl
         << "    -s, --subgraphs       describe subgraphs of graph" << endl
         << "    -H, --heads           list the head nodes of the graph" << endl
         << "    -T, --tails           list the tail nodes of the graph" << endl;
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

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"size", no_argument, 0, 'z'},
                {"length", no_argument, 0, 'l'},
                {"subgraphs", no_argument, 0, 's'},
                {"heads", no_argument, 0, 'H'},
                {"tails", no_argument, 0, 'T'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hzlsHT",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'z':
            stats_size = true;
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

    delete graph;

    return 0;

}

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -n, --node ID         starting at node with ID" << endl
         << "    -l, --max-length N    generate paths of at most length N" << endl
         << "    -e, --edge-max N      cross no more than N edges when determining k-paths" << endl
         << "    -s, --as-seqs         write each path as a sequence" << endl;
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

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"node", required_argument, 0, 'n'},
                {"max-length", required_argument, 0, 'l'},
                {"edge-max", required_argument, 0, 'e'},
                {"as-seqs", no_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "n:l:hse:",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
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

    if (max_length == 0) {
        cerr << "error:[vg paths] a --max-length is required when generating paths" << endl;
    }

    function<void(Node*,Path&)> paths_to_seqs = [graph](Node* n, Path& p) {
        string seq = graph->path_sequence(p);
#pragma omp critical(cout)
        cout << seq << endl;
    };

    function<void(Node*,Path&)> paths_to_json = [](Node* n, Path& p) {
        char *json2 = pb2json(p);
#pragma omp critical(cout)
        cout<<json2<<endl;
        free(json2);
    };

    function<void(Node*, Path&)>* callback = &paths_to_seqs;
    if (!as_seqs) {
        callback = &paths_to_json;
    }

    if (node_id) {
        graph->for_each_kpath_of_node(graph->get_node(node_id), max_length, edge_max, *callback);
    } else {
        graph->for_each_kpath_parallel(max_length, edge_max, *callback);
    }

    delete graph;

    return 0;

}

void help_find(char** argv) {
    cerr << "usage: " << argv[0] << " find [options] <graph.vg> >sub.vg" << endl
         << "options:" << endl
         << "    -n, --node ID          find node, return 1-hop context as graph" << endl
         << "    -f, --edges-from ID    return edges from node with ID" << endl
         << "    -t, --edges-to ID      return edges from node with ID" << endl
         << "    -k, --kmer STR         return a graph of edges and nodes matching this kmer" << endl
         << "    -T, --table            instead of a graph, return a table of kmers" << endl
         << "    -c, --context STEPS    expand the context of the kmer hit subgraphs" << endl
         << "    -s, --sequence STR     search for sequence STR using --kmer-size kmers" << endl
         << "    -j, --kmer-stride N    step distance between succesive kmers in sequence (default 1)" << endl
         << "    -z, --kmer-size N      split up --sequence into kmers of size N" << endl
         << "    -C, --kmer-count       report approximate count of kmer (-k) in db" << endl
         << "    -p, --path TARGET      find the node(s) in the specified path range TARGET=path[:pos1[-pos2]]" << endl
         << "    -P, --position-in PATH find the position of the node (specified by -n) in the given path" << endl
         << "    -d, --db-name DIR      use this db (defaults to <graph>.index/)" << endl;
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
    int64_t from_id=0, to_id=0;
    vector<int64_t> node_ids;
    int context_size=0;
    bool count_kmers = false;
    bool kmer_table = false;
    string target;
    string path_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"db-name", required_argument, 0, 'd'},
                {"node", required_argument, 0, 'n'},
                {"edges-from", required_argument, 0, 'f'},
                {"edges-to", required_argument, 0, 't'},
                {"kmer", required_argument, 0, 'k'},
                {"table", no_argument, 0, 'T'},
                {"sequence", required_argument, 0, 's'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"kmer-size", required_argument, 0, 'z'},
                {"output", required_argument, 0, 'o'},
                {"context", required_argument, 0, 'c'},
                {"kmer-count", no_argument, 0, 'C'},
                {"path", required_argument, 0, 'p'},
                {"position-in", required_argument, 0, 'P'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:n:f:t:o:k:hc:s:z:j:CTp:P:",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'd':
            db_name = optarg;
            break;

        case 'k':
            kmers.push_back(optarg);
            break;

        case 's':
            sequence = optarg;
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

        case 'f':
            from_id = atoi(optarg);
            break;

        case 't':
            to_id = atoi(optarg);
            break;

        case 'T':
            kmer_table = true;
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
        string file_name = argv[optind];
        if (file_name == "-") {
            if (db_name.empty()) {
                cerr << "error:[vg find] reading variant graph from stdin and no db name (-d) given, exiting" << endl;
                return 1;
            }
        }
        ifstream in;
        if (db_name.empty()) {
            db_name = file_name + ".index";
        }
        in.open(file_name.c_str());
    }

    Index index;
    index.open_read_only(db_name);

    if (!node_ids.empty() && path_name.empty()) {
        // open index
        // our result
        // get the context of the node
        vector<VG> graphs;
        for (auto node_id : node_ids) {
            VG g;
            index.get_context(node_id, g);
            if (context_size > 0) {
                index.expand_context(g, context_size);
            }
            graphs.push_back(g);
        }
        VG result_graph;
        for (auto& graph : graphs) {
            result_graph.extend(graph);
        }
        result_graph.remove_orphan_edges();
        // return it
        result_graph.serialize_to_ostream(cout);
    } else if (from_id != 0) {
        vector<Edge> edges;
        index.get_edges_from(from_id, edges);
        for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
            cout << e->from() << "\t" << e->to() << endl;
        }
    } else if (to_id != 0) {
        vector<Edge> edges;
        index.get_edges_to(to_id, edges);
        for (vector<Edge>::iterator e = edges.begin(); e != edges.end(); ++e) {
            cout << e->from() << "\t" << e->to() << endl;
        }
    }

    if (!node_ids.empty() && !path_name.empty()) {
        int64_t path_id = index.get_path_id(path_name);
        for (auto node_id : node_ids) {
            list<int64_t> path_prev, path_next;
            int64_t prev_pos=0, next_pos=0;
            if (index.get_node_path_relative_position(node_id, path_id,
                                                      path_prev, prev_pos, path_next, next_pos)) {
                cout << node_id << " " << path_prev.front() << " " << prev_pos
                     << " " << path_next.back() << " " << next_pos << endl;
            }
        }
    }

    if (!target.empty()) {
        string name;
        int64_t start, end;
        VG graph;
        parse_region(target, name, start, end);
        index.get_path(graph, name, start, end);
        if (context_size > 0) {
            index.expand_context(graph, context_size);
        }
        graph.serialize_to_ostream(cout);
    }

    // todo cleanup if/else logic to allow only one function
    
    if (!sequence.empty()) {
        set<int> kmer_sizes = index.stored_kmer_sizes();
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
    }

    if (!kmers.empty()) {
        if (count_kmers) {
            for (auto& kmer : kmers) {
                cout << kmer << "\t" << index.approx_size_of_kmer_matches(kmer) << endl;
            }
        } else if (kmer_table) {
            for (auto& kmer : kmers) {
                map<string, vector<pair<int64_t, int32_t> > > positions;
                index.get_kmer_positions(kmer, positions);
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
                index.get_kmer_subgraph(kmer, g);
                if (context_size > 0) {
                    index.expand_context(g, context_size);
                }
                graphs.push_back(g);
            }

            VG result_graph;
            for (auto& graph : graphs) {
                result_graph.extend(graph);
            }
            result_graph.remove_orphan_edges();
            result_graph.serialize_to_ostream(cout);
        }
    }

    return 0;

}

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "options:" << endl
         << "    -s, --store           store graph (do this first to build db!)" << endl
         << "    -k, --kmer-size N     index kmers of size N in the graph" << endl
         << "    -e, --edge-max N     cross no more than N edges when determining k-paths" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers in paths (default 1)" << endl
         << "    -P, --prune KB        remove kmer entries which use more than KB kilobytes" << endl
         << "    -D, --dump            print the contents of the db to stdout" << endl
         << "    -M, --metadata        describe aspects of the db stored in metadata" << endl
         << "    -S, --set-kmer        assert that the kmer size (-k) is in the db" << endl
         << "    -d, --db-name DIR     create rocksdb in DIR (defaults to <graph>.index/)" << endl
         << "                          (this is required if you are using multiple graphs files" << endl
        //<< "    -b, --tmp-db-base S   use this base name for temporary indexes" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -p, --progress        show progress" << endl;
}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    string db_name;
    int kmer_size = 0;
    int edge_max = 0;
    int kmer_stride = 1;
    int prune_kb = -1;
    bool store_graph = false;
    bool dump_index = false;
    bool describe_index = false;
    bool show_progress = false;
    bool set_kmer_size = false;

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
                {"store", no_argument, 0, 's'},
                {"dump", no_argument, 0, 'D'},
                {"metadata", no_argument, 0, 'M'},
                {"set-kmer", no_argument, 0, 'S'},
                {"threads", required_argument, 0, 't'},
                {"progress",  no_argument, 0, 'p'},
                {"prune",  required_argument, 0, 'P'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:k:j:pDshMt:b:e:SP:",
                         long_options, &option_index);
        
        // Detect the end of the options.
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'd':
            db_name = optarg;
            break;

        case 'P':
            prune_kb = atoi(optarg);
            break;

        case 'k':
            kmer_size = atoi(optarg);
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

        case 'S':
            set_kmer_size = true;
            break;

        case 's':
            store_graph = true;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
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

    vector<string> graph_file_names;
    while (optind < argc) {
        string file_name = argv[optind++];
        graph_file_names.push_back(file_name);
    }

    if (db_name.empty()) {
        if (graph_file_names.size() > 1) {
            cerr << "error:[vg index] working on multiple graphs and no db name (-d) given, exiting" << endl;
            return 1;
        } else if (graph_file_names.size() == 1) {
            db_name = *graph_file_names.begin() + ".index";
        } else {
            cerr << "error:[vg index] no graph or db given, exiting" << endl;
            return 1;
        }
    }

    Index index;

    if (store_graph && graph_file_names.size() > 0) {
        index.open_for_bulk_load(db_name);
        VGset graphs(graph_file_names);
        graphs.show_progress = show_progress;
        graphs.store_in_index(index);
        index.flush();
        index.close();
        // reopen to index paths
        // this requires the index to be queryable
        index.open_for_write(db_name);
        index.compact();
        graphs.store_paths_in_index(index);
        index.compact();
        index.flush();
        index.close();
    }

    if (kmer_size != 0 && graph_file_names.size() > 0) {
        index.open_for_bulk_load(db_name);
        VGset graphs(graph_file_names);
        graphs.show_progress = show_progress;
        graphs.index_kmers(index, kmer_size, edge_max, kmer_stride);
        index.flush();
        index.close();
        index.open_for_write(db_name);
        index.compact();
        index.flush();
        index.close();
    }

    if (prune_kb >= 0) {
        if (show_progress) {
            cerr << "pruning kmers > " << prune_kb << " on disk from " << db_name << endl;
        }
        index.open_for_write(db_name);
        index.prune_kmers(prune_kb);
        index.compact();
        index.close();
    }

    if (set_kmer_size) {
        assert(kmer_size != 0);
        index.open_for_write(db_name);
        index.remember_kmer_size(kmer_size);
        index.close();
    }

    if (dump_index) {
        index.open_read_only(db_name);
        index.dump(cout);
        index.close();
    }

    if (describe_index) {
        index.open_read_only(db_name);
        set<int> kmer_sizes = index.stored_kmer_sizes();
        cout << "kmer sizes: ";
        for (auto kmer_size : kmer_sizes) {
            cout << kmer_size << " ";
        }
        cout << endl;
        index.close();
    }

    return 0;

}

void help_align(char** argv) {
    cerr << "usage: " << argv[0] << " align [options] <graph.vg> >alignments.vga" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
        //<< "    -p, --print-cigar     output graph cigar for alignments" << endl
         << "    -j, --json            output alignments in JSON format (default)" << endl;
}

int main_align(int argc, char** argv) {

    string seq;

    if (argc == 2) {
        help_align(argv);
        return 1;
    }

    bool print_cigar = false;
    bool output_json = true;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:jhd:",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 's':
            seq = optarg;
            break;

            /*
        case 'p':
            print_cigar = true;
            break;
            */

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

    char *json2 = pb2json(alignment);
    cout<<json2<<endl;
    free(json2);

    delete graph;

    return 0;

}

void help_map(char** argv) {
    cerr << "usage: " << argv[0] << " map [options] <graph.vg> >alignments.vga" << endl
         << "options:" << endl
         << "    -d, --db-name DIR     use this db (defaults to <graph>.index/)" << endl
         << "                          a graph is not required" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -r, --reads FILE      take reads (one per line) from FILE, write alignments to stdout" << endl
         << "    -b, --hts-input FILE  align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout" << endl
         << "    -k, --kmer-size N     use this kmer size, it must be < kmer size in db (default: from index)" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers to use for seeding (default: kmer size)" << endl
         << "    -S, --sens-step N     decrease kmer size by N bp until alignment succeeds (default 5)" << endl
         << "    -c, --clusters N      use at most the largest N ordered clusters of the kmer graph for alignment" << endl
         << "    -m, --hit-max N       ignore kmers who have >N hits in our index (default 100)" << endl
         << "    -t, --threads N       number of threads to use" << endl
         << "    -F, --prefer-forward  if the forward alignment of the read works, accept it" << endl
         << "    -X, --score-per-bp N  accept forward if the alignment score per base is > N and -F is set" << endl
         << "    -J, --output-json     output JSON rather than an alignment stream (helpful for debugging)" << endl
         << "    -D, --debug           print debugging information about alignment to stderr" << endl;
}

int main_map(int argc, char** argv) {

    if (argc == 2) {
        help_map(argv);
        return 1;
    }

    string seq;
    string db_name;
    int kmer_size = 0;
    int kmer_stride = 0;
    int sens_step = 0;
    int best_clusters = 0;
    string read_file;
    string hts_file;
    int hit_max = 100;
    int thread_count = 1;
    bool output_json = false;
    bool debug = false;
    bool prefer_forward = false;
    float score_per_bp = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {"db-name", required_argument, 0, 'd'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"kmer-size", required_argument, 0, 'k'},
                {"clusters", required_argument, 0, 'c'},
                {"reads", required_argument, 0, 'r'},
                {"hit-max", required_argument, 0, 'm'},
                {"threads", required_argument, 0, 't'},
                {"prefer-forward", no_argument, 0, 'F'},
                {"score-per-bp", required_argument, 0, 'X'},
                {"sens-step", required_argument, 0, 'S'},
                {"output-json", no_argument, 0, 'J'},
                {"hts-input", no_argument, 0, 'b'},
                {"debug", no_argument, 0, 'D'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:j:hd:c:r:m:k:t:DX:FS:Jb:",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 's':
            seq = optarg;
            break;

        case 'd':
            db_name = optarg;
            break;

        case 'j':
            kmer_stride = atoi(optarg);
            break;

        case 'k':
            kmer_size = atoi(optarg);
            break;

        case 'S':
            sens_step = atoi(optarg);
            break;

        case 'c':
            best_clusters = atoi(optarg);
            break;

        case 'm':
            hit_max = atoi(optarg);
            break;

        case 'r':
            read_file = optarg;
            break;

        case 'b':
            hts_file = optarg;
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

        case 'X':
            score_per_bp = atof(optarg);
            break;

        case 'J':
            output_json = true;
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

    if (seq.empty() && read_file.empty() && hts_file.empty()) {
        cerr << "error:[vg map] a sequence or read file is required when mapping" << endl;
        return 1;
    }

    string file_name;
    if (optind < argc) {
        file_name = argv[optind];
    }

    if (db_name.empty()) {
        if (file_name.empty()) {
            cerr << "error:[vg map] no graph or db given, exiting" << endl;
            return 1;
        } else {
            db_name = file_name + ".index";
        }
    }

    thread_count = get_thread_count();

    vector<Mapper*> mapper;
    mapper.resize(thread_count);
    vector<vector<Alignment> > output_buffer;
    output_buffer.resize(thread_count);

    Index idx;
    idx.open_read_only(db_name);

    for (int i = 0; i < thread_count; ++i) {
        Mapper* m = new Mapper(&idx);
        m->best_clusters = best_clusters;
        m->hit_max = hit_max;
        m->debug = debug;
        if (score_per_bp) m->target_score_per_bp = score_per_bp;
        if (sens_step) m->kmer_sensitivity_step = sens_step;
        m->prefer_forward = prefer_forward;
        mapper[i] = m;
    }

    if (!seq.empty()) {
        int tid = omp_get_thread_num();
        Alignment alignment = mapper[tid]->align(seq, kmer_size, kmer_stride);
        if (output_json) {
            char *json2 = pb2json(alignment);
            cout<<json2<<endl;
            free(json2);
        } else {
            function<Alignment(uint64_t)> lambda =
                [&alignment] (uint64_t n) {
                return alignment;
            };
            stream::write(cout, 1, lambda);
        }
    }

    if (!read_file.empty()) {
        ifstream in(read_file);
        bool more_data = true;
#pragma omp parallel shared(in, more_data)
        {
            string line;
            int tid = omp_get_thread_num();
            while (more_data) {
                line.clear();
#pragma omp critical (readq)
                {
                    more_data = std::getline(in,line);
                }
                if (!line.empty()) {
                    Alignment alignment = mapper[tid]->align(line, kmer_size, kmer_stride);
                    if (output_json) {
                        char *json2 = pb2json(alignment);
#pragma omp critical (cout)
                        cout << json2 << "\n";
                        free(json2);
                    } else {
                        auto& output_buf = output_buffer[tid];
                        output_buf.push_back(alignment);
                        if (output_buf.size() >= 1000) {
                            function<Alignment(uint64_t)> lambda =
                                [&output_buf] (uint64_t n) {
                                return output_buf[n];
                            };
#pragma omp critical (cout)
                            stream::write(cout, output_buf.size(), lambda);
                            output_buf.clear();
                        }
                    }
                }
            }
        }
    }


    if (!hts_file.empty()) {
        function<void(Alignment&)> lambda =
            [&mapper,
             &output_buffer,
             &output_json,
             &kmer_size,
             &kmer_stride]
            (Alignment& alignment) {
            int tid = omp_get_thread_num();
            alignment = mapper[tid]->align(alignment, kmer_size, kmer_stride);
            if (output_json) {
                char *json2 = pb2json(alignment);
#pragma omp critical (cout)
                cout << json2 << "\n";
                free(json2);
            } else {
                auto& output_buf = output_buffer[tid];
                output_buf.push_back(alignment);
                if (output_buf.size() >= 1000) {
                    function<Alignment(uint64_t)> lambda =
                        [&output_buf] (uint64_t n) {
                        return output_buf[n];
                    };
#pragma omp critical (cout)
                    stream::write(cout, output_buf.size(), lambda);
                    output_buf.clear();
                }
            }
        };
        // run
        hts_for_each_parallel(hts_file, lambda);
    }

    // clean up
    for (int i = 0; i < thread_count; ++i) {
        delete mapper[i];
        auto& output_buf = output_buffer[i];
        if (!output_buf.empty()) {
            function<Alignment(uint64_t)> lambda =
                [&output_buf] (uint64_t n) {
                return output_buf[n];
            };
            stream::write(cout, output_buf.size(), lambda);
            output_buf.clear();
        }
    }

    cout.flush();

    return 0;

}

void help_view(char** argv) {
    cerr << "usage: " << argv[0] << " view [options] [<graph.vg>|<paths.vgp>]" << endl
         << "options:" << endl
         << "    -g, --gfa          output GFA format (default)" << endl
         << "    -d, --dot          output dot format" << endl
         << "    -j, --json         output JSON format" << endl
         << "    -v, --vg           write VG format (input is GFA)" << endl
         << "    -a, --align-in     input is GAM (vg alignment format: Graph Alignment/Map)" << endl
         << "    -G, --gam          output GAM format alignments" << endl
         << "    -b, --bam          input is htslib-parseable alignments" << endl;
    //<< "    -p, --paths           extract paths from graph in VG format" << endl;
}

int main_view(int argc, char** argv) {

    if (argc == 2) {
        help_view(argv);
        return 1;
    }

    string output_type;
    string input_type;

    int c;
    optind = 2; // force optind past "view" argument
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"dot", no_argument, 0, 'd'},
                {"gfa", no_argument, 0, 'g'},
                {"json",  no_argument, 0, 'j'},
                {"vg", no_argument, 0, 'v'},
                {"paths", no_argument, 0, 'p'},
                {"align-in", no_argument, 0, 'a'},
                {"gam", no_argument, 0, 'G'},
                {"bam", no_argument, 0, 'b'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgjhvpaGb",
                         long_options, &option_index);
        
        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c)
        {
        case 'd':
            output_type = "dot";
            break;
 
        case 'g':
            output_type = "gfa";
            break;
 
        case 'j':
            output_type = "json";
            break;

        case 'v':
            input_type = "gfa";
            output_type = "vg";
            break;

        case 'a':
            input_type = "gam";
            break;

        case 'G':
            output_type = "gam";
            break;

        case 'b':
            input_type = "bam";
            break;

        case 'p':
            input_type = "paths";
            output_type = "paths";
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

    if (output_type.empty()) {
        if (input_type == "gam") {
            output_type = "json";
        } else {
            output_type = "gfa";
        }
    }

    if (input_type.empty()) {
        input_type = "vg";
    }

    VG* graph;
    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }
    string file_name = argv[optind];
    if (input_type == "vg") {
        if (file_name == "-") {
            graph = new VG(std::cin);
        } else {
            ifstream in;
            in.open(file_name.c_str());
            graph = new VG(in);
        }
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
    } else if (input_type == "paths") {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG;
        graph->paths.load(in);
    } else if (input_type == "gam") {
        if (output_type == "json") {
            // convert values to printable ones
            function<void(Alignment&)> lambda = [](Alignment& a) {
                alignment_quality_short_to_char(a);
                char *json2 = pb2json(a);
                cout << json2 << "\n";
                free(json2);
            };
            if (file_name == "-") {
                stream::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(file_name.c_str());
                stream::for_each(in, lambda);
            }
            return 0;
        } else {
            // todo
            return 0;
        }
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
            return 0;
        } else if (output_type == "json") {
            // todo
            return 0;
        }
    }

    if (output_type == "dot") {
        graph->to_dot(std::cout);
    } else if (output_type == "json") {
        char *json2 = pb2json(graph->graph);
        cout<<json2<<endl;
        free(json2);
    } else if (output_type == "gfa") {
        graph->to_gfa(std::cout);
    } else if (output_type == "vg") {
        graph->serialize_to_ostream(cout);
    } else if (output_type == "paths") {
        function<void(Path&)> dump_path = [](Path& p) {
            for (int i = 0; i < p.mapping_size(); ++i) {
                cout << p.name() << "\t" << p.mapping(i).node_id() << endl;
            }
        };
        graph->paths.for_each(dump_path);
    }

    delete graph;

    return 0;
}

void help_construct(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options] >new.vg" << endl
         << "options:" << endl
         << "    -v, --vcf FILE        input VCF" << endl
         << "    -r, --reference FILE  input FASTA reference" << endl
         << "    -P, --ref-paths FILE  write reference paths in protobuf/gzip format to FILE" << endl
         << "    -R, --region REGION   specify a particular chromosome" << endl
         << "    -z, --region-size N   variants per region to parallelize" << endl
         << "    -m, --node-max N      limit the maximum allowable node sequence size" << endl
         << "                          nodes greater than this threshold will be divided" << endl
         << "    -p, --progress        show progress" << endl
         << "    -t, --threads N       use N threads to construct graph (defaults to numCPUs)" << endl;
}

int main_construct(int argc, char** argv) {

    if (argc == 2) {
        help_construct(argv);
        return 1;
    }

    string fasta_file_name, vcf_file_name;
    string region;
    string output_type = "VG";
    bool progress = false;
    int vars_per_region = 25000;
    int max_node_size = 0;
    string ref_paths_file;

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"ref-paths", required_argument, 0, 'P'},
                {"progress",  no_argument, 0, 'p'},
                {"region-size", required_argument, 0, 'z'},
                {"threads", required_argument, 0, 't'},
                {"region", required_argument, 0, 'R'},
                {"node-max", required_argument, 0, 'm'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:phz:t:R:m:P:",
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

        case 'p':
            progress = true;
            break;
 
        case 'z':
            vars_per_region = atoi(optarg);
            break;

        case 'R':
            region = optarg;
            break;

        case 't':
            omp_set_num_threads(atoi(optarg));
            break;

        case 'm':
            max_node_size = atoi(optarg);
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

    // set up our inputs

    vcflib::VariantCallFile variant_file;
    if (vcf_file_name.empty()) {
        cerr << "error:[vg construct] a VCF file is required for graph construction" << endl;
        return 1;
    }
    variant_file.open(vcf_file_name);
    if (!variant_file.is_open()) {
        cerr << "error:[vg construct] could not open" << vcf_file_name << endl;
        return 1;
    }

    FastaReference reference;
    if (fasta_file_name.empty()) {
        cerr << "error:[vg construct] a reference is required for graph construction" << endl;
        return 1;
    }
    reference.open(fasta_file_name);

    // store our reference sequence paths
    Paths ref_paths;

    VG graph(variant_file, reference, region, vars_per_region, max_node_size, progress);

    if (!ref_paths_file.empty()) {
        ofstream paths_out(ref_paths_file);
        graph.paths.write(paths_out);
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
         << "  -- view          format conversions for graphs and alignments" << endl
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
         << "  -- sim           simulate reads from the graph" << endl;
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
    } else {
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

    return 0;

}
