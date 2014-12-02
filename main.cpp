#include <iostream>
#include <fstream>
#include <getopt.h>
#include "pb2json.h"
#include "vg.pb.h"
#include "Variant.h"
#include "Fasta.h"
#include "index.h"
#include "vg.h"
#include "leveldb/db.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;

void vg_help(char** argv) {
    cerr << "usage: " << argv[0] << " <command> [options]" << endl
         << endl
         << "commands:" << endl 
         << "  -- construct     graph construction" << endl
         << "  -- view          conversion (protobuf/json/GFA)" << endl
         << "  -- index         index features of the graph in a disk-backed key/value store" << endl
         << "  -- find          use an index to find nodes, edges, kmers, or positions" << endl
         << "  -- paths         traverse paths in the graph" << endl
         << "  -- align         local alignment" << endl
         << "  -- stats         metrics describing graph properties" << endl
         << "  -- join          combine graphs" << endl;
}

void help_align(char** argv) {
    cerr << "usage: " << argv[0] << " align [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
        //<< "    -p, --print-cigar     output graph cigar for alignments" << endl
         << "    -j, --json            output alignments in JSON format (default)" << endl;
}

void help_view(char** argv) {
    cerr << "usage: " << argv[0] << " view [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -d, --dot             output dot format (default)" << endl
         << "    -g, --gfa             output GFA format" << endl
         << "    -j, --json            output VG JSON format" << endl;
}

void help_construct(char** argv) {
    cerr << "usage: " << argv[0] << " construct [options]" << endl
         << "options:" << endl
         << "    -v, --vcf FILE        input VCF" << endl
         << "    -r, --reference FILE  input FASTA reference" << endl
         << "    -p, --protobuf        output VG protobuf format (default)" << endl
         << "    -g, --gfa             output GFA format" << endl
         << "    -j, --json            output VG JSON format" << endl;
}

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -s, --store           store graph (do this first to build db!)" << endl
         << "    -k, --kmer-size N     index kmers of size N in the graph" << endl
         << "    -p, --positions       index nodes and edges by position" << endl
         << "    -D, --dump            print the contents of the db to stdout" << endl
         << "    -d, --db-name DIR     create leveldb in DIR (defaults to <graph>.index/)" << endl;
}

void help_find(char** argv) {
    cerr << "usage: " << argv[0] << " find [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -n, --node ID         find node, return 1-hop context as graph" << endl
        // << "    -f, --edges-from ID   return edges from node with ID" << endl
        // << "    -t, --edges-to ID     return edges from node with ID" << endl
         << "    -k, --kmer STR        return a list of edges and nodes matching this kmer" << endl
         << "    -c, --context STEPS   expand the context of the kmer hit subgraphs" << endl
         << "    -s, --sequence STR    search for sequence STR using --kmer-size kmers" << endl
         << "    -z, --kmer-size N     split up --sequence into kmers of size N" << endl
        // << "    -o, --output FORMAT   use this output format for found elements (default: JSON)" << endl
         << "    -d, --db-name DIR     use this db (defaults to <graph>.index/)" << endl;
}

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -n, --node ID         starting at node with ID" << endl
         << "    -l, --max-length N    generate paths of at most length N" << endl;
}

void help_stats(char** argv) {
    cerr << "usage: " << argv[0] << " stats [options] <graph.vg>" << endl
         << "options:" << endl
         << "    -z, --size            size of graph" << endl
         << "    -l, --length          length of sequences in graph" << endl
         << "    -s, --subgraphs       describe subgraphs of graph" << endl;
}

void help_join(char** argv) {
    cerr << "usage: " << argv[0] << " join [options] <graph1.vg> [graph2.vg ...]" << endl
         << "joins graphs and sub-graphs into a single variant graph" << endl
         << "a single root node with sequence 'N' is added" << endl
         << "assumes a single id namespace for all graphs to join" << endl;
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
    joined.graph.SerializeToOstream(&std::cout);

    return 0;
}

int main_stats(int argc, char** argv) {

    if (argc == 2) {
        help_stats(argv);
        return 1;
    }

    bool stats_size = false;
    bool stats_length = false;
    bool stats_subgraphs = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"size", no_argument, 0, 'z'},
                {"length", no_argument, 0, 'l'},
                {"subgraphs", no_argument, 0, 's'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hzls",
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

int main_paths(int argc, char** argv) {

    if (argc == 2) {
        help_paths(argv);
        return 1;
    }

    int max_length = 0;
    int64_t node_id = 0;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"node", required_argument, 0, 'n'},
                {"max-length", required_argument, 0, 'l'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "n:l:h",
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

        case 'h':
        case '?':
            help_paths(argv);
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

    vector<Path> paths;
    if (node_id) {
        graph->bounded_paths(node_id, paths, max_length);
    } else {
        graph->bounded_paths(paths, max_length);
    }

    for (vector<Path>::iterator p = paths.begin(); p != paths.end(); ++p) {
        char *json2 = pb2json(*p);
        cout<<json2<<endl;
        free(json2);
    }

    delete graph;

    return 0;

}

int main_find(int argc, char** argv) {

    if (argc == 2) {
        help_find(argv);
        return 1;
    }

    string db_name;
    string sequence;
    int kmer_size=0;
    vector<string> kmers;
    string output_format;
    int64_t node_id=0, from_id=0, to_id=0;
    int context_size=0;

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
                {"sequence", required_argument, 0, 's'},
                {"kmer-size", required_argument, 0, 'z'},
                {"output", required_argument, 0, 'o'},
                {"context", required_argument, 0, 'c'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:n:f:t:o:k:hc:s:z:",
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

        case 'z':
            kmer_size = atoi(optarg);
            break;

        case 'c':
            context_size = atoi(optarg);
            break;

        case 'n':
            node_id = atoi(optarg);
            break;

        case 'f':
            from_id = atoi(optarg);
            break;

        case 't':
            to_id = atoi(optarg);
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

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        if (db_name.empty()) {
            cerr << "error:[vg find] reading variant graph from stdin and no db name (-d) given, exiting" << endl;
            return 1;
        }
        graph = new VG(std::cin);
    } else {
        ifstream in;
        if (db_name.empty()) {
            db_name = file_name + ".index";
        }
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    Index index(db_name);

    if (node_id != 0) {
        // open index
        // our result
        VG result_graph;
        // get the context of the node
        index.get_context(node_id, result_graph);
        if (context_size > 0) {
            index.expand_context(result_graph, context_size);
        }
        // return it
        result_graph.graph.SerializeToOstream(&cout);
    }

    if (!sequence.empty()) {
        if (kmer_size == 0) {
            cerr << "--kmer-size is required when supplying --sequence" << endl;
            return 1;
        }
        for (int i = 0; i < sequence.size()-kmer_size; ++i) {
            kmers.push_back(sequence.substr(i,kmer_size));
        }
    }

    if (!kmers.empty()) {
        VG result_graph;
        for (vector<string>::iterator k = kmers.begin(); k != kmers.end(); ++k) {
            VG g;
            index.get_kmer_subgraph(*k, g);
            if (context_size > 0) {
                index.expand_context(g, context_size);
            }
            result_graph.extend(g);
        }
        result_graph.graph.SerializeToOstream(&cout);
    }
    
    delete graph;

    return 0;

}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    string db_name;
    bool index_by_position = false;
    int kmer_size = 0;
    bool store_graph = false;
    bool dump_index = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"db-name", required_argument, 0, 'd'},
                {"kmer-size", required_argument, 0, 'k'},
                {"positions", no_argument, 0, 'p'},
                {"store", no_argument, 0, 's'},
                {"dump", no_argument, 0, 'D'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:k:pDsh",
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
            kmer_size = atoi(optarg);
            break;

        case 'p':
            index_by_position = true;
            break;

        case 'D':
            dump_index = true;
            break;

        case 's':
            store_graph = true;
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

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        if (db_name.empty()) {
            cerr << "error:[vg index] reading variant graph from stdin and no db name (-d) given, exiting" << endl;
            return 1;
        }
        graph = new VG(std::cin);
    } else {
        ifstream in;
        if (db_name.empty()) {
            db_name = file_name + ".index";
        }
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    Index index(db_name);

    if (store_graph) {
        index.load_graph(*graph);
    }

    if (kmer_size != 0) {
        map<string, map<Node*, int> > kmer_map;
        graph->kmers_of(kmer_map, kmer_size);
        index.store_kmers(kmer_map);
    }

    if (dump_index) {
        index.dump(cout);
    }

    delete graph;

    return 0;

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
        c = getopt_long (argc, argv, "s:jh",
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

int main_view(int argc, char** argv) {

    if (argc == 2) {
        help_view(argv);
        return 1;
    }

    string output_type = "dot";

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
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgjh",
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
            output_type = "GFA";
            break;
 
        case 'j':
            output_type = "JSON";
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

    VG* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new VG(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new VG(in);
    }

    if (output_type == "dot") {
        graph->to_dot(std::cout);
    } else if (output_type == "JSON") {
        char *json2 = pb2json(graph->graph);
        cout<<json2<<endl;
        free(json2);
    }

    delete graph;

    return 0;
}

int main_construct(int argc, char** argv) {

    if (argc == 2) {
        help_construct(argv);
        return 1;
    }

    string fasta_file_name, vcf_file_name;
    string output_type = "VG";

    int c;
    while (true) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"vcf", required_argument, 0, 'v'},
                {"reference", required_argument, 0, 'r'},
                {"protobuf",  no_argument, 0, 'p'},
                {"gfa",  no_argument, 0, 'g'},
                {"json",  no_argument, 0, 'j'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:r:pgjh",
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

        case 'p':
            output_type = "VG";
            break;
 
        case 'g':
            output_type = "GFA";
            break;
 
        case 'j':
            output_type = "JSON";
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

    vcf::VariantCallFile variant_file;
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

    VG graph(variant_file, reference);

    if (output_type == "VG") {
        //ofstream of("test.vg");
        graph.graph.SerializeToOstream(&std::cout);
    } else if (output_type == "JSON") {
        char *json2 = pb2json(graph.graph);
        cout<<json2<<endl;
        free(json2);
    }

    return 0;
}

int main(int argc, char *argv[])
{

    if (argc == 1) {
        vg_help(argv);
        return 1;
    }

    string command = argv[1];
    if (command == "construct") {
        return main_construct(argc, argv);
    } else if (command == "view") {
        return main_view(argc, argv);
    } else if (command == "align") {
        return main_align(argc, argv);
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
    } else {
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

    return 0;

}
