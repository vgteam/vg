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
#include "vectorizer.hpp"
#include "sampler.hpp"
#include "filter.hpp"
#include "google/protobuf/stubs/common.h"
#include "progress_bar.hpp"
#include "version.hpp"
#include "genotyper.hpp"
#include "bubbles.hpp"
#include "translator.hpp"
#include "homogenize_main.cpp"
#include "sift_main.cpp"
#include "srpe_main.cpp"
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

void help_filter(char** argv) {
    cerr << "usage: " << argv[0] << " filter [options] <alignment.gam> > out.gam" << endl
         << "Filter low-scoring alignments using different heuristics." << endl
         << endl
         << "options:" << endl
         << "    -s, --min-secondary N   minimum score to keep secondary alignment [default=0]" << endl
         << "    -r, --min-primary N     minimum score to keep primary alignment [default=0]" << endl
         << "    -f, --frac-score        normalize score based on length" << endl
         << "    -u, --substitutions     use substitution count instead of score" << endl
         << "    -o, --max-overhang N    filter reads whose alignments begin or end with an insert > N [default=99999]" << endl
         << "    -S, --drop-split        remove split reads taking nonexistent edges" << endl
         << "    -x, --xg-name FILE      use this xg index (required for -R, -S, and -D)" << endl
         << "    -R, --regions-file      only output alignments that intersect regions (BED file with 0-based coordinates expected)" << endl
         << "    -B, --output-basename   output to file(s) (required for -R).  The ith file will correspond to the ith BED region" << endl
         << "    -A, --append-regions    append to alignments created with -RB" << endl
         << "    -c, --context STEPS     expand the context of the subgraph this many steps when looking up chunks" << endl
         << "    -v, --verbose           print out statistics on numbers of reads filtered by what." << endl
         << "    -q, --min-mapq N        filter alignments with mapping quality < N" << endl
         << "    -E, --repeat-ends N     filter reads with tandem repeat (motif size <= 2N, spanning >= N bases) at either end" << endl
         << "    -D, --defray-ends N     clip back the ends of reads that are ambiguously aligned, up to N bases" << endl
         << "    -C, --defray-count N    stop defraying after N nodes visited (used to keep runtime in check) [default=99999]" << endl
         << "    -t, --threads N         number of threads [1]" << endl;
}

int main_filter(int argc, char** argv) {

    if (argc <= 2) {
        help_filter(argv);
        return 1;
    }

    // This is the better design for a subcommand: we have a class that
    // implements it and encapsulates all the default parameters, and then we
    // just feed in overrides in the option parsing code. Thsi way we don't have
    // multiple defaults all over the place.
    ReadFilter filter;

    // What XG index, if any, should we load to support the other options?
    string xg_name;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"min-secondary", required_argument, 0, 's'},
                {"min-primary", required_argument, 0, 'r'},
                {"frac-score", required_argument, 0, 'f'},
                {"substitutions", required_argument, 0, 'u'},
                {"max-overhang", required_argument, 0, 'o'},
                {"drop-split",  no_argument, 0, 'S'},
                {"xg-name", required_argument, 0, 'x'},
                {"regions-file",  required_argument, 0, 'R'},
                {"output-basename",  required_argument, 0, 'B'},
                {"append-regions", no_argument, 0, 'A'},
                {"context",  required_argument, 0, 'c'},
                {"verbose",  no_argument, 0, 'v'},
                {"min-mapq", required_argument, 0, 'q'},
                {"repeat-ends", required_argument, 0, 'E'},
                {"defray-ends", required_argument, 0, 'D'},
                {"defray-count", required_argument, 0, 'C'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:r:d:e:fauo:Sx:R:B:Ac:vq:E:D:C:t:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 's':
            filter.min_secondary = atof(optarg);
            break;
        case 'r':
            filter.min_primary = atof(optarg);
            break;
        case 'f':
            filter.frac_score = true;
            break;
        case 'u':
            filter.sub_score = true;
            break;
        case 'o':
            filter.max_overhang = atoi(optarg);
            break;
        case 'S':
            filter.drop_split = true;
        case 'x':
            xg_name = optarg;
            break;
        case 'R':
            filter.regions_file = optarg;
            break;
        case 'B':
            filter.outbase = optarg;
            break;
        case 'A':
            filter.append_regions = true;
            break;
        case 'c':
            filter.context_size = atoi(optarg);
            break;
        case 'q':
            filter.min_mapq = atof(optarg);
            break;
        case 'v':
            filter.verbose = true;
            break;
        case 'E':
            filter.repeat_size = atoi(optarg);
            break;
        case 'D':
            filter.defray_length = atoi(optarg);
            break;
        case 'C':
            filter.defray_count = atoi(optarg);
            break;          
        case 't':
            filter.threads = atoi(optarg);
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_filter(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    omp_set_num_threads(filter.threads);
    
    // setup alignment stream
    if (optind >= argc) {
        help_filter(argv);
        return 1;
    }

    // What should our return code be?
    int error_code = 0;

    get_input_file(optind, argc, argv, [&](istream& in) {
        // Open up the alignment stream
        
        // If the user gave us an XG index, we probably ought to load it up.
        // TODO: make sure if we add any other error exits from this function we
        // remember to delete this!
        xg::XG* xindex = nullptr;
        if (!xg_name.empty()) {
            // read the xg index
            ifstream xg_stream(xg_name);
            if(!xg_stream) {
                cerr << "Unable to open xg index: " << xg_name << endl;
                error_code = 1;
                return;
            }
            xindex = new xg::XG(xg_stream);
        }
    
        // Read in the alignments and filter them.
        error_code = filter.filter(&in, xindex);
        
        if(xindex != nullptr) {
            delete xindex;
        }
    });

    return error_code;
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

void help_vectorize(char** argv){
    cerr << "usage: " << argv[0] << " vectorize [options] -x <index.xg> <alignments.gam>" << endl

         << "Vectorize a set of alignments to a variety of vector formats." << endl
         << endl
         << "options: " << endl
         << "  -x --xg FILE       An xg index for the graph of interest" << endl
         << "  -g --gcsa FILE     A gcsa2 index to use if generating MEM sketches" << endl
         << "  -l --aln-label LABEL   Rename every alignment to LABEL when outputting alignment name." << endl
         << "  -f --format        Tab-delimit output so it can be used in R." << endl
         << "  -A --annotate      Create a header with each node/edge's name and a column with alignment names." << endl
         << "  -a --a-hot         Instead of a 1-hot, output a vector of {0|1|2} for covered, reference, or alt." << endl
         << "  -w --wabbit        Output a format that's friendly to vowpal wabbit" << endl
         << "  -M --wabbit-mapping <FILE> output the mappings used for vowpal wabbit classes (default: print to stderr)" << endl
         << "  -m --mem-sketch    Generate a MEM sketch of a given read based on the GCSA" << endl
         << "  -p --mem-positions Add the positions to the MEM sketch of a given read based on the GCSA" << endl
         << "  -H --mem-hit-max N Ignore MEMs with this many hits when extracting poisitions" << endl
         << "  -i --identity-hot  Output a score vector based on percent identity and coverage" << endl
         << endl;
}

int main_vectorize(int argc, char** argv){

    string xg_name;
    string read_file = "";
    string aln_label = "";
    string gcsa_name;
    string wabbit_mapping_file = "";
    bool format = false;
    bool show_header = false;
    bool map_alns = false;
    bool annotate = false;
    bool a_hot = false;
    bool output_wabbit = false;
    bool use_identity_hot = false;
    bool mem_sketch = false;
    bool mem_positions = false;
    bool mem_hit_max = 0;
    int max_mem_length = 0;

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
            {"xg", required_argument,0, 'x'},
            {"gcsa", required_argument,0, 'g'},
            {"threads", required_argument, 0, 't'},
            {"format", no_argument, 0, 'f'},
            {"a-hot", no_argument, 0, 'a'},
            {"wabbit", no_argument, 0, 'w'},
            {"wabbit-mapping", required_argument, 0, 'M'},
            {"mem-sketch", no_argument, 0, 'm'},
            {"mem-positions", no_argument, 0, 'p'},
            {"mem-hit-max", required_argument, 0, 'H'},
            {"identity-hot", no_argument, 0, 'i'},
            {"aln-label", required_argument, 0, 'l'},
            {"reads", required_argument, 0, 'r'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "AaihwM:fmpx:g:l:H:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

 /*           case '?':
            case 'h':
                help_vectorize(argv);
                return 1;
            case 'x':
                xg_name = optarg;
                break;
            case 'a':
                a_hot = true;
                break;
            case 'w':
                output_wabbit = true;
                break; */
            case 'l':
                aln_label = optarg;
                break;
            case 'r':
                read_file = optarg;
                break;
                /*

            case 'i':
                use_identity_hot = true;
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
*/

        case '?':
        case 'h':
            help_vectorize(argv);
            return 1;
        case 'x':
            xg_name = optarg;
            break;
        case 'g':
            gcsa_name = optarg;
            break;
        case 'm':
            mem_sketch = true;
            break;
        case 'p':
            mem_positions = true;
            break;
        case 'H':
            mem_hit_max = atoi(optarg);
            break;
        case 'a':
            a_hot = true;
            break;
        case 'w':
            output_wabbit = true;
            break;
        case 'i':
            use_identity_hot = true;
            break;
        case 'f':
            format = true;
            break;
        case 'A':
            annotate = true;
            format = true;
            break;
        case 'M':
            wabbit_mapping_file = optarg;
            break;
        default:
            abort();
        }
    }

    xg::XG* xg_index;
    if (!xg_name.empty()) {
        ifstream in(xg_name);
        xg_index = new xg::XG(in);
    }
    else{
        cerr << "No XG index given. An XG index must be provided." << endl;
        exit(1);
    }

    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);

    gcsa::GCSA gcsa_index;
    gcsa::LCPArray lcp_index;
    if (!gcsa_name.empty()) {
        ifstream in_gcsa(gcsa_name.c_str());
        gcsa_index.load(in_gcsa);
        // default LCP is the gcsa base name +.lcp
        string lcp_in = gcsa_name + ".lcp";
        ifstream in_lcp(lcp_in.c_str());
        lcp_index.load(in_lcp);
    }

    Mapper mapper;
    if (mem_sketch) {
        if (gcsa_name.empty()) {
            cerr << "[vg vectorize] error : an xg index and gcsa index are required when making MEM sketches" << endl;
            return 1;
        } else {
            mapper.gcsa = &gcsa_index;
            mapper.lcp = &lcp_index;
        }
        if (mem_hit_max) {
            mapper.hit_max = mem_hit_max;
        }
    }

    Vectorizer vz(xg_index);

    //Generate a 1-hot coverage vector for graph entities.
    function<void(Alignment&)> lambda = [&vz, &mapper, use_identity_hot, output_wabbit, aln_label, mem_sketch, mem_positions, format, a_hot, max_mem_length](Alignment& a){
        //vz.add_bv(vz.alignment_to_onehot(a));
        //vz.add_name(a.name());
        if (a_hot) {
            vector<int> v = vz.alignment_to_a_hot(a);
            if (output_wabbit){
                cout << vz.wabbitize(aln_label == "" ? a.name() : aln_label, v) << endl;
            }
            else if (format){
                cout << a.name() << "\t" << vz.format(v) << endl;
            } else{
                cout << v << endl;
            }
        }
        else if (use_identity_hot){
            vector<double> v = vz.alignment_to_identity_hot(a);
            if (output_wabbit){
                cout << vz.wabbitize(aln_label == "" ? a.name() : aln_label, v) << endl;
            }
            else if (format){
                cout << a.name() << "\t" << vz.format(v) << endl;
            }
            else {
                cout << vz.format(v) << endl;
            }

        } else if (mem_sketch) {
            // get the mems
            map<string, int> mem_to_count;
            auto mems = mapper.find_mems(a.sequence().begin(), a.sequence().end(), max_mem_length);
            for (auto& mem : mems) {
                mem_to_count[mem.sequence()]++;
            }
            cout << " |info count:" << mems.size() << " unique:" << mem_to_count.size();
            cout << " |mems";
            for (auto m : mem_to_count) {
                cout << " " << m.first << ":" << m.second;
            }
            if (mem_positions) {
                cout << " |positions";
                for (auto& mem : mems) {
                    mapper.get_mem_hits_if_under_max(mem);
                    for (auto& node : mem.nodes) {
                        cout << " " << gcsa::Node::id(node);
                        if (gcsa::Node::rc(node)) {
                            cout << "-";
                        } else {
                            cout << "+";
                        }
                        cout << ":" << mem.end - mem.begin;
                    }
                }
            }
            cout << endl;
        } else {
            bit_vector v = vz.alignment_to_onehot(a);
            if (output_wabbit){
                cout << vz.wabbitize(aln_label == "" ? a.name() : aln_label, v) << endl;
            } else if (format) {
                cout << a.name() << "\t" << vz.format(v) << endl;
            } else{
                cout << v << endl;
            }
        }
    };
    
    get_input_file(optind, argc, argv, [&](istream& in) {
        stream::for_each(in, lambda);
    });

    string mapping_str = vz.output_wabbit_map();
    if (output_wabbit){
        if (!wabbit_mapping_file.empty()){
            ofstream ofi;
            ofi.open(wabbit_mapping_file);
            if (!ofi.good()){
                cerr << "Error with outputting wabbit mapping file. Make sure the filename is a valid string" << endl;
                return 1;
            }
            ofi << mapping_str;
            ofi.close();
        }
        else{

            cerr << mapping_str;
        }
    }



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

void help_call(char** argv) {
    cerr << "usage: " << argv[0] << " call [options] <graph.vg> <pileup.vgpu> > output.vcf" << endl
         << "Output variant calls in VCF format given a graph and pileup" << endl
         << endl
         << "options:" << endl
         << "    -d, --min_depth INT        minimum depth of pileup [" << Caller::Default_min_depth <<"]" << endl
         << "    -e, --max_depth INT        maximum depth of pileup [" << Caller::Default_max_depth <<"]" << endl
         << "    -s, --min_support INT      minimum number of reads required to support snp [" << Caller::Default_min_support <<"]" << endl
         << "    -f, --min_frac FLOAT       minimum percentage of reads required to support snp[" << Caller::Default_min_frac <<"]" << endl
         << "    -q, --default_read_qual N  phred quality score to use if none found in the pileup ["
         << (int)Caller::Default_default_quality << "]" << endl
         << "    -b, --max_strand_bias FLOAT limit to absolute difference between 0.5 and proportion of supporting reads on reverse strand. [" << Caller::Default_max_strand_bias << "]" << endl
         << "    -a, --link-alts            add all possible edges between adjacent alts" << endl
         << "    -A, --aug-graph FILE       write out the agumented graph in vg format" << endl
         << "    -r, --ref PATH             use the given path name as the reference path" << endl
         << "    -c, --contig NAME          use the given name as the VCF contig name" << endl
         << "    -S, --sample NAME          name the sample in the VCF with the given name [SAMPLE]" << endl
         << "    -o, --offset INT           offset variant positions by this amount in VCF [0]" << endl
         << "    -l, --length INT           override total sequence length in VCF" << endl
         << "    -P, --pileup               write pileup under VCF lines (for debugging, output not valid VCF)" << endl
         << "    -D, --depth INT            maximum depth for path search [default 10 nodes]" << endl
         << "    -F, --min_cov_frac FLOAT   min fraction of average coverage at which to call [0.0]" << endl
         << "    -H, --max_het_bias FLOAT   max imbalance factor between alts to call heterozygous [3]" << endl
         << "    -R, --max_ref_bias FLOAT   max imbalance factor between ref and alts to call heterozygous ref [4]" << endl
         << "    -M, --bias_mult FLOAT      multiplier for bias limits for indels as opposed to substitutions [1]" << endl
         << "    -n, --min_count INT        min total supporting read count to call a variant [1]" << endl
         << "    -B, --bin_size  INT        bin size used for counting coverage [250]" << endl
         << "    -C, --exp_coverage INT     specify expected coverage (instead of computing on reference)" << endl
         << "    -O, --no_overlap           don't emit new variants that overlap old ones" << endl
         << "    -u, --use_avg_support      use average instead of minimum support" << endl
         << "    -I, --singleallelic        disable support for multiallelic sites" << endl
         << "    -E, --min_mad              min. minimum allele depth required to PASS filter [5]" << endl
         << "    -h, --help                 print this help message" << endl
         << "    -p, --progress             show progress" << endl
         << "    -v, --verbose              print information and warnings about vcf generation" << endl
         << "    -t, --threads N            number of threads to use" << endl;
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
    double min_frac = Caller::Default_min_frac;
    int default_read_qual = Caller::Default_default_quality;
    double max_strand_bias = Caller::Default_max_strand_bias;
    string aug_file;
    bool bridge_alts = false;
    // Option variables (formerly from glenn2vcf)
    // What's the name of the reference path in the graph?
    string refPathName = "";
    // What name should we give the contig in the VCF file?
    string contigName = "";
    // What name should we use for the sample in the VCF file?
    string sampleName = "SAMPLE";
    // How far should we offset positions of variants?
    int64_t variantOffset = 0;
    // How many nodes should we be willing to look at on our path back to the
    // primary path? Keep in mind we need to look at all valid paths (and all
    // combinations thereof) until we find a valid pair.
    int64_t maxDepth = 10;
    // Should we write pileup information to the VCF for debugging?
    bool pileupAnnotate = false;
    // What should the total sequence length reported in the VCF header be?
    int64_t lengthOverride = -1;
    // What fraction of average coverage should be the minimum to call a variant (or a single copy)?
    double minFractionForCall = 0;
    // What fraction of the reads supporting an alt are we willing to discount?
    // At 2, if twice the reads support one allele as the other, we'll call
    // homozygous instead of heterozygous. At infinity, every call will be
    // heterozygous if even one read supports each allele.
    double maxHetBias = 3;
    // Like above, but applied to ref / alt ratio (instead of alt / ref)
    double maxRefHetBias = 4;
    // How many times more bias do we allow for indels?
    double indelBiasMultiple = 1;
    // What's the minimum integer number of reads that must support a call? We
    // don't necessarily want to call a SNP as het because we have a single
    // supporting read, even if there are only 10 reads on the site.
    size_t minTotalSupportForCall = 1;
    // Bin size used for counting coverage along the reference path.  The
    // bin coverage is used for computing the probability of an allele
    // of a certain depth
    size_t refBinSize = 250;
    // On some graphs, we can't get the coverage because it's split over
    // parallel paths.  Allow overriding here
    size_t expCoverage = 0;
    // Should we drop variants that would overlap old ones? TODO: we really need
    // a proper system for accounting for usage of graph material.
    bool suppress_overlaps = false;
    // Should we use average support instead minimum support for our calculations?
    bool useAverageSupport = false;
    // Should we go by sites and thus support multiallelic sites (true), or use
    // the old single-branch-bubble method (false)?
    bool multiallelic_support = true;
    // How big a site should we try to type all at once instead of replacing
    // with its children if it has any?
    size_t max_ref_length = 100;
    // What's the maximum number of bubble path combinations we can explore
    // while finding one with maximum support?
    size_t max_bubble_paths = 100;
    // what's the minimum minimum allele depth to give a PASS in the filter column
    // (anything below gets FAIL)
    size_t min_mad_for_filter = 5;

    bool show_progress = false;
    bool verbose = false;
    int thread_count = 1;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"min_depth", required_argument, 0, 'd'},
                {"max_depth", required_argument, 0, 'e'},
                {"min_support", required_argument, 0, 's'},
                {"min_frac", required_argument, 0, 'f'},
                {"default_read_qual", required_argument, 0, 'q'},
                {"max_strand_bias", required_argument, 0, 'b'},
                {"aug_graph", required_argument, 0, 'A'},
                {"link-alts", no_argument, 0, 'a'},
                {"progress", no_argument, 0, 'p'},
                {"verbose", no_argument, 0, 'v'},
                {"threads", required_argument, 0, 't'},
                {"ref", required_argument, 0, 'r'},
                {"contig", required_argument, 0, 'c'},
                {"sample", required_argument, 0, 'S'},
                {"offset", required_argument, 0, 'o'},
                {"depth", required_argument, 0, 'D'},
                {"length", required_argument, 0, 'l'},
                {"pileup", no_argument, 0, 'P'},
                {"min_cov_frac", required_argument, 0, 'F'},
                {"max_het_bias", required_argument, 0, 'H'},
                {"max_ref_bias", required_argument, 0, 'R'},
                {"bias_mult", required_argument, 0, 'M'},
                {"min_count", required_argument, 0, 'n'},
                {"bin_size", required_argument, 0, 'B'},
                {"avg_coverage", required_argument, 0, 'C'},
                {"no_overlap", no_argument, 0, 'O'},
                {"use_avg_support", no_argument, 0, 'u'},
                {"singleallelic", no_argument, 0, 'I'},
                {"min_mad", required_argument, 0, 'E'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:e:s:f:q:b:A:apvt:r:c:S:o:D:l:PF:H:R:M:n:B:C:OuIE:h",
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
        case 'f':
            min_frac = atof(optarg);
            break;
        case 'q':
            default_read_qual = atoi(optarg);
            break;
        case 'b':
            max_strand_bias = atof(optarg);
            break;
        case 'A':
            aug_file = optarg;
            break;
        case 'a':
            bridge_alts = true;
            break;
        // old glenn2vcf opts start here
        case 'r':
            // Set the reference path name
            refPathName = optarg;
            break;
        case 'c':
            // Set the contig name
            contigName = optarg;
            break;
        case 'S':
            // Set the sample name
            sampleName = optarg;
            break;
        case 'o':
            // Offset variants
            variantOffset = std::stoll(optarg);
            break;
        case 'D':
            // Limit max depth for pathing to primary path
            maxDepth = std::stoll(optarg);
            break;
        case 'l':
            // Set a length override
            lengthOverride = std::stoll(optarg);
            break;
        case 'P':
            pileupAnnotate = true;
            break;
        case 'F':
            // Set min fraction of average coverage for a call
            minFractionForCall = std::stod(optarg);
            break;
        case 'H':
            // Set max factor between reads on one alt and reads on the other
            // alt for calling a het.
            maxHetBias = std::stod(optarg);
            break;
        case 'R':
            // Set max factor between reads on ref and reads on the other
            // alt for calling a homo ref.
            maxRefHetBias = std::stod(optarg);
            break;
        case 'M':
            // Set multiplier for bias limits for indels
            indelBiasMultiple = std::stod(optarg);
            break;
        case 'n':
            // How many reads need to touch an allele before we are willing to
            // call it?
            minTotalSupportForCall = std::stoll(optarg);
            break;
        case 'B':
            // Set the reference bin size
            refBinSize = std::stoll(optarg);
            break;
        case 'C':
            // Override expected coverage
            expCoverage = std::stoll(optarg);
            break;
        case 'O':
            // Suppress variants that overlap others
            suppress_overlaps = true;
            break;
        case 'u':
            // Average (isntead of min) support
            useAverageSupport = true;
            break;
        case 'I':
            // Disallow for multiallelic sites by using a different algorithm
            multiallelic_support = false;
            break;
        case 'E':
            // Minimum min-allele-depth required to give Filter column a PASS
            min_mad_for_filter = std::stoi(optarg);
            break;
        case 'p':
            show_progress = true;
            break;
        case 'v':
            verbose = true;
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

    // Parse the arguments
    if (optind >= argc) {
        help_call(argv);
        return 1;
    }
    string graph_file_name = get_input_file_name(optind, argc, argv);
    if (optind >= argc) {
        help_call(argv);
        return 1;
    }
    string pileup_file_name = get_input_file_name(optind, argc, argv);
    
    if (pileup_file_name == "-" && graph_file_name == "-") {
        cerr << "error: graph and pileup can't both be from stdin." << endl;
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

    // this is the call tsv file that was used to communicate with glenn2vcf
    // it's still here for the time being but never actually written
    // (just passed as a string to the caller)
    stringstream text_file_stream;

    if (show_progress) {
        cerr << "Computing augmented graph" << endl;
    }
    Caller caller(graph,
                  het_prior, min_depth, max_depth, min_support,
                  min_frac, Caller::Default_min_log_likelihood,
                  true, default_read_qual, max_strand_bias,
                  &text_file_stream, bridge_alts);

    // setup pileup stream
    get_input_file(pileup_file_name, [&](istream& pileup_stream) {
        // compute the augmented graph
        function<void(Pileup&)> lambda = [&caller](Pileup& pileup) {
            for (int i = 0; i < pileup.node_pileups_size(); ++i) {
                caller.call_node_pileup(pileup.node_pileups(i));
            }
            for (int i = 0; i < pileup.edge_pileups_size(); ++i) {
                caller.call_edge_pileup(pileup.edge_pileups(i));
            }
        };
        stream::for_each(pileup_stream, lambda);
    });
    
    // map the edges from original graph
    if (show_progress) {
        cerr << "Mapping edges into augmented graph" << endl;
    }
    caller.update_call_graph();

    // map the paths from the original graph
    if (show_progress) {
        cerr << "Mapping paths into augmented graph" << endl;
    }
    caller.map_paths();

    if (!aug_file.empty()) {
        // write the augmented graph
        if (show_progress) {
            cerr << "Writing augmented graph" << endl;
        }
        ofstream aug_stream(aug_file.c_str());
        caller.write_call_graph(aug_stream, false);
    }

    if (show_progress) {
        cerr << "Calling variants" << endl;
    }

    // project the augmented graph to a reference path
    // in order to create a VCF of calls.  this
    // was once a separate tool called glenn2vcf
    glenn2vcf::call2vcf(caller._call_graph,
                        text_file_stream.str(),
                        refPathName,
                        contigName,
                        sampleName,
                        variantOffset,
                        maxDepth,
                        lengthOverride,
                        pileupAnnotate ? pileup_file_name : string(),
                        minFractionForCall,
                        maxHetBias,
                        maxRefHetBias,
                        indelBiasMultiple,
                        minTotalSupportForCall,
                        refBinSize,
                        expCoverage,
                        suppress_overlaps,
                        useAverageSupport,
                        multiallelic_support,
                        max_ref_length,
                        max_bubble_paths,
                        min_mad_for_filter,
                        verbose);

    return 0;
}

void help_genotype(char** argv) {
    cerr << "usage: " << argv[0] << " genotype [options] <graph.vg> <reads.index/> > <calls.vcf>" << endl
         << "Compute genotypes from a graph and an indexed collection of reads" << endl
         << endl
         << "options:" << endl
         << "    -j, --json              output in JSON" << endl
         << "    -v, --vcf               output in VCF" << endl
         << "    -r, --ref PATH          use the given path name as the reference path" << std::endl
         << "    -c, --contig NAME       use the given name as the VCF contig name" << std::endl
         << "    -s, --sample NAME       name the sample in the VCF with the given name" << std::endl
         << "    -o, --offset INT        offset variant positions by this amount" << std::endl
         << "    -l, --length INT        override total sequence length" << std::endl
         << "    -a, --augmented FILE    dump augmented graph to FILE" << std::endl
         << "    -q, --use_mapq          use mapping qualities" << std::endl
         << "    -C, --cactus            use cactus ultrabubbles for site finding" << std::endl
         << "    -S, --subset-graph      only use the reference and areas of the graph with read support" << std::endl
         << "    -i, --realign_indels    realign at indels" << std::endl
         << "    -d, --het_prior_denom   denominator for prior probability of heterozygousness" << std::endl
         << "    -P, --min_per_strand    min consistent reads per strand for an allele" << std::endl
         << "    -p, --progress          show progress" << endl
         << "    -t, --threads N         number of threads to use" << endl;
}

int main_genotype(int argc, char** argv) {

    if (argc <= 3) {
        help_genotype(argv);
        return 1;
    }
    // Should we output genotypes in JSON (true) or Protobuf (false)?
    bool output_json = false;
    // Should we output VCF instead of protobuf?
    bool output_vcf = false;
    // Should we show progress with a progress bar?
    bool show_progress = false;
    // How many threads should we use?
    int thread_count = 0;

    // What reference path should we use
    string ref_path_name;
    // What sample name should we use for output
    string sample_name;
    // What contig name override do we want?
    string contig_name;
    // What offset should we add to coordinates
    int64_t variant_offset = 0;
    // What length override should we use
    int64_t length_override = 0;

    // Should we use mapping qualities?
    bool use_mapq = false;
    // Should we do indel realignment?
    bool realign_indels = false;

    // Should we dump the augmented graph to a file?
    string augmented_file_name;

    // Should we do superbubbles/sites with Cactus (true) or supbub (false)
    bool use_cactus = false;
    // Should we find superbubbles on the supported subset (true) or the whole graph (false)?
    bool subset_graph = false;
    // What should the heterozygous genotype prior be? (1/this)
    double het_prior_denominator = 10.0;
    // At least how many reads must be consistent per strand for a call?
    size_t min_consistent_per_strand = 2;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"json", no_argument, 0, 'j'},
                {"vcf", no_argument, 0, 'v'},
                {"ref", required_argument, 0, 'r'},
                {"contig", required_argument, 0, 'c'},
                {"sample", required_argument, 0, 's'},
                {"offset", required_argument, 0, 'o'},
                {"length", required_argument, 0, 'l'},
                {"augmented", required_argument, 0, 'a'},
                {"use_mapq", no_argument, 0, 'q'},
                {"cactus", no_argument, 0, 'C'},
                {"subset-graph", no_argument, 0, 'S'},
                {"realign_indels", no_argument, 0, 'i'},
                {"het_prior_denom", required_argument, 0, 'd'},
                {"min_per_strand", required_argument, 0, 'P'},
                {"progress", no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hjvr:c:s:o:l:a:qCSid:P:pt:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'j':
            output_json = true;
            break;
        case 'v':
            output_vcf = true;
            break;
        case 'r':
            // Set the reference path name
            ref_path_name = optarg;
            break;
        case 'c':
            // Set the contig name for output
            contig_name = optarg;
            break;
        case 's':
            // Set the sample name
            sample_name = optarg;
            break;
        case 'o':
            // Offset variants
            variant_offset = std::stoll(optarg);
            break;
        case 'l':
            // Set a length override
            length_override = std::stoll(optarg);
            break;
        case 'a':
            // Dump augmented graph
            augmented_file_name = optarg;
            break;
        case 'q':
            // Use mapping qualities
            use_mapq = true;
            break;
        case 'C':
            // Use Cactus to find sites
            use_cactus = true;
            break;
        case 'S':
            // Find sites on the graph subset with any read support
            subset_graph = true;
            break;
        case 'i':
            // Do indel realignment
            realign_indels = true;
            break;
        case 'd':
            // Set heterozygous genotype prior denominator
            het_prior_denominator = std::stod(optarg);
            break;
        case 'P':
            // Set min consistent reads per strand required to keep an allele
            min_consistent_per_strand = std::stoll(optarg);
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
            help_genotype(argv);
            exit(1);
            break;
        default:
          abort ();
        }
    }

    if(thread_count > 0) {
        omp_set_num_threads(thread_count);
    }

    // read the graph
    if (optind >= argc) {
        help_genotype(argv);
        return 1;
    }
    if (show_progress) {
        cerr << "Reading input graph..." << endl;
    }
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    // setup reads index
    if (optind >= argc) {
        help_genotype(argv);
        return 1;
    }

    string reads_index_name = get_input_file_name(optind, argc, argv);
    // This holds the RocksDB index that has all our reads, indexed by the nodes they visit.
    Index index;
    index.open_read_only(reads_index_name);

    // Build the set of all the node IDs to operate on
    vector<vg::id_t> graph_ids;
    graph->for_each_node([&](Node* node) {
        // Put all the ids in the set
        graph_ids.push_back(node->id());
    });

    // Load all the reads matching the graph into memory
    vector<Alignment> alignments;

    if(show_progress) {
        cerr << "Loading reads..." << endl;
    }

    index.for_alignment_to_nodes(graph_ids, [&](const Alignment& alignment) {
        // Extract all the alignments

        // Only take alignments that don't visit nodes not in the graph
        bool contained = true;
        for(size_t i = 0; i < alignment.path().mapping_size(); i++) {
            if(!graph->has_node(alignment.path().mapping(i).position().node_id())) {
                // Throw out the read
                contained = false;
            }
        }

        if(contained) {
            // This alignment completely falls within the graph
            alignments.push_back(alignment);
        }
    });

    if(show_progress) {
        cerr << "Loaded " << alignments.size() << " alignments" << endl;
    }

    // Make a Genotyper to do the genotyping
    Genotyper genotyper;
    // Configure it
    genotyper.use_mapq = use_mapq;
    genotyper.realign_indels = realign_indels;
    assert(het_prior_denominator > 0);
    genotyper.het_prior_logprob = prob_to_logprob(1.0/het_prior_denominator);
    genotyper.min_consistent_per_strand = min_consistent_per_strand;
    // TODO: move arguments below up into configuration
    genotyper.run(*graph,
                  alignments,
                  cout,
                  ref_path_name,
                  contig_name,
                  sample_name,
                  augmented_file_name,
                  use_cactus,
                  subset_graph,
                  show_progress,
                  output_vcf,
                  output_json,
                  length_override,
                  variant_offset);

    delete graph;

    return 0;
}

void help_pileup(char** argv) {
    cerr << "usage: " << argv[0] << " pileup [options] <graph.vg> <alignment.gam> > out.vgpu" << endl
         << "Calculate pileup for each position in graph and output in VG Pileup format (list of protobuf NodePileups)." << endl
         << endl
         << "options:" << endl
         << "    -j, --json              output in JSON" << endl
         << "    -q, --min-quality N     ignore bases with PHRED quality < N (default=0)" << endl
         << "    -m, --max-mismatches N  ignore bases with > N mismatches within window centered on read (default=1)" << endl
         << "    -w, --window-size N     size of window to apply -m option (default=0)" << endl
         << "    -d, --max-depth N       maximum depth pileup to create (further maps ignored) (default=1000)" << endl
         << "    -a, --use-mapq          combine mapping qualities with base qualities" << endl
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
    int min_quality = 0;
    int max_mismatches = 1;
    int window_size = 0;
    int max_depth = 1000; // used to prevent protobuf messages getting to big
    bool verbose = false;
    bool use_mapq = false;

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
                {"use-mapq", no_argument, 0, 'a'},
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
        case 'a':
            use_mapq = true;
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
        help_call(argv);
        return 1;
    }
    string graph_file_name = get_input_file_name(optind, argc, argv);
    if (optind >= argc) {
        help_call(argv);
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
         << "local alignment parameters:" << endl
         << "    -a, --match N         use this match score (default: 1)" << endl
         << "    -i, --mismatch N      use this mismatch penalty (default: 4)" << endl
         << "    -o, --gap-open N      use this gap open penalty (default: 6)" << endl
         << "    -e, --gap-extend N    use this gap extension penalty (default: 1)" << endl
         << "mem mapping:" << endl
         << "    -L, --min-mem-length N  ignore SMEMs shorter than this length (default: 0/unset)" << endl
         << "    -Y, --max-mem-length N  ignore SMEMs longer than this length by stopping backward search (default: 0/unset)" << endl
         << "    -H, --hit-max N         SMEMs which have >N hits in our index (default: 100)" << endl
         << "    -c, --context-depth N   follow this many steps out from each subgraph for alignment (default: 7)" << endl
         << "    -T, --thread-ex N       cluster nodes when successive ids are within this distance (default: 10)" << endl
         << "    -P, --min-identity N    accept alignment only if the alignment-based identity is >= N (default: 0.75)" << endl
         << "    -B, --band-width N      use this bandwidth when mapping (default: 256)" << endl
         << "    -G, --greedy-accept     if a tested alignment achieves -S identity don't try clusters with fewer hits" << endl
         << "    -S, --accept-identity N accept early alignment if the alignment identity is >= N and -G is set" << endl
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
         << "    -C, --circularize       the input sequences are from circular genomes, circularize them after inclusion" << endl
         << "generic parameters:" << endl
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
    int hit_max = 100;
    int max_attempts = 10;
    // if we set this above 1, we use a dynamic programming process to determine the
    // optimal alignment through a series of bands based on a proximity metric
    int max_multimaps = 1;
    // if this is set too low, we may miss optimal alignments
    int context_depth = 7;
    // same here; initial clustering
    int thread_extension = 10;
    float min_identity = 0.0;
    int band_width = 256;
    size_t doubling_steps = 3;
    bool debug = false;
    bool debug_align = false;
    size_t node_max = 0;
    int alignment_threads = get_thread_count();
    int edge_max = 0;
    int subgraph_prune = 0;
    bool normalize = false;
    bool allow_nonpath = false;
    int iter_max = 1;
    int max_mem_length = 0;
    int min_mem_length = 8;
    bool greedy_accept = false;
    float accept_identity = 0;
    int max_target_factor = 100;
    bool idx_path_only = false;
    int match = 1;
    int mismatch = 4;
    int gap_open = 6;
    int gap_extend = 1;
    bool circularize = false;
    int sens_step = 5;

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
                {"min-identity", required_argument, 0, 'P'},
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
                {"accept-identity", required_argument, 0, 'S'},
                {"max-attempts", required_argument, 0, 'M'},
                {"thread-ex", required_argument, 0, 'T'},
                {"max-target-x", required_argument, 0, 'q'},
                {"max-multimaps", required_argument, 0, 'I'},
                {"match", required_argument, 0, 'a'},
                {"mismatch", required_argument, 0, 'i'},
                {"gap-open", required_argument, 0, 'o'},
                {"gap-extend", required_argument, 0, 'e'},
                {"circularize", no_argument, 0, 'C'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hf:n:s:g:b:K:X:B:DAc:P:E:Q:NzI:L:Y:H:t:m:GS:M:T:q:OI:a:i:o:e:C",
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
            accept_identity = atof(optarg);
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

        case 'C':
            circularize = true;
            break;

        case 'P':
            min_identity = atof(optarg);
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

        case 'a':
            match = atoi(optarg);
            break;

        case 'i':
            mismatch = atoi(optarg);
            break;

        case 'o':
            gap_open = atoi(optarg);
            break;

        case 'e':
            gap_extend = atoi(optarg);
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
        get_input_file(file_name, [&](istream& in) {
            graph = new VG(in);
        });
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
            if (!seq_names.empty() && seq_names.count(name) == 0) continue;
            // only use the sequence if we have whitelisted it
            string seq = ref.getSequence(name);
            strings[name] = seq;
            names_in_order.push_back(name);
        }
    }

    // give a label to sequences passed on the command line
    // use the sha1sum, take the head
    // collision avoidance with nonce ensures we get the same names for the same sequences across multiple runs
    for (auto& s : sequences) {
        auto name = sha1head(s, 8);
        int nonce = 0;
        while (strings.find(name) != strings.end()) {
            stringstream ss;
            ss << s << ++nonce;
            name = sha1head(ss.str(), 8);
        }
        strings[name] = s;
        names_in_order.push_back(name);
    }

    // always go from biggest to smallest in progression
    sort(names_in_order.begin(), names_in_order.end(),
         [&strings](const string& s1, const string& s2) {
             return strings[s1].size() > strings[s2].size(); });

    // align, include, repeat

    if (debug) cerr << "preparing initial graph" << endl;

    size_t max_query_size = pow(2, doubling_steps) * idx_kmer_size;
    // limit max node size
    if (!node_max) node_max = 2*idx_kmer_size;

    // if our graph is empty, we need to take the first sequence and build a graph from it
    if (graph->empty()) {
        auto build_graph = [&graph,&node_max](const string& seq, const string& name) {
            graph->create_node(seq);
            graph->dice_nodes(node_max);
            graph->sort();
            graph->compact_ids();
            // the graph will have a single embedded path in it
            Path& path = *graph->graph.add_path();
            path.set_name(name);
            graph->for_each_node([&path](Node* node) {
                    auto mapping = path.add_mapping();
                    mapping->mutable_position()->set_node_id(node->id());
                    auto edit = mapping->add_edit();
                    edit->set_from_length(node->sequence().size());
                    edit->set_to_length(node->sequence().size());
                });
            graph->paths.extend(path);
            graph->sync_paths(); // sync the paths
        };
        // what's the first sequence?
        if (base_seq_name.empty()) {
            base_seq_name = names_in_order.front();
        }
        // we specified one we wanted to use as the first
        build_graph(strings[base_seq_name], base_seq_name);
    }

#ifdef debug
    cerr << "path going in " << pb2json(graph->graph.path(0)) << endl;
#endif

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
        if (mapper) delete mapper;
        if (xgidx) delete xgidx;
        if (gcsaidx) delete gcsaidx;
        if (lcpidx) delete lcpidx;
    
        //stringstream s; s << iter++ << ".vg";
        graph->sort();
        graph->sync_paths();
        graph->graph.clear_path();
        graph->paths.to_graph(graph->graph);
        graph->rebuild_indexes();

        if (debug) cerr << "building xg index" << endl;
        xgidx = new xg::XG(graph->graph);

        if (debug) cerr << "building GCSA2 index" << endl;
        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if(!debug) gcsa::Verbosity::set(gcsa::Verbosity::SILENT);

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
            mapper->min_identity = min_identity;
            mapper->min_mem_length = min_mem_length;
            mapper->hit_max = hit_max;
            mapper->greedy_accept = greedy_accept;
            mapper->max_target_factor = max_target_factor;
            mapper->max_multimaps = max_multimaps;
            mapper->accept_identity = accept_identity;
            mapper->mem_threading = true;

            // set up the multi-threaded alignment interface
            // TODO abstract this into a single call!!
            mapper->alignment_threads = alignment_threads;
            mapper->clear_aligners(); // number of aligners per mapper depends on thread count
                                      // we have to reset this here to re-init scores to the right number
            mapper->set_alignment_scores(match, mismatch, gap_open, gap_extend);
            mapper->init_node_cache();
            mapper->init_node_pos_cache();
            mapper->mem_threading = true;
        }
    };

    // set up the graph for mapping
    rebuild(graph);

    // todo restructure so that we are trying to map everything
    // add alignment score/bp bounds to catch when we get a good alignment
    for (auto& name : names_in_order) {
        //cerr << "do... " << name << " ?" << endl;
        if (!base_seq_name.empty() && name == base_seq_name) continue; // already embedded
        bool incomplete = true; // complete when we've fully included the sequence set
        int iter = 0;
        auto& seq = strings[name];
        //cerr << "doing... " << name << endl;
#ifdef debug
        {
            graph->serialize_to_file("msga-pre-" + name + ".vg");
            ofstream db_out("msga-pre-" + name + ".xg");
            xgidx->serialize(db_out);
            db_out.close();
        }
#endif
        while (incomplete && iter++ < iter_max) {
            stringstream s; s << iter; string iterstr = s.str();
            if (debug) cerr << name << ": adding to graph" << iter << endl;
            vector<Path> paths;
            vector<Alignment> alns;
            int j = 0;
            // align to the graph
            if (debug) cerr << name << ": aligning sequence of " << seq.size() << "bp against " <<
                graph->node_count() << " nodes" << endl;
            Alignment aln = simplify(mapper->align(seq, 0, sens_step, max_mem_length, band_width));
            if (aln.path().mapping_size()) {
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
            } else {
                Edit* edit = aln.mutable_path()->add_mapping()->add_edit();
                edit->set_sequence(aln.sequence());
                edit->set_to_length(aln.sequence().size());
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
            if (circularize) {
                if (debug) cerr << name << ": circularizing" << endl;
                graph->circularize({name});
                //graph->serialize_to_file(name + "-post-circularize.vg");
            }

            // the edit needs to cut nodes at mapping starts and ends
            // thus allowing paths to be included that map directly to entire nodes
            // XXX

            //graph->serialize_to_file(name + "-pre-index.vg");
            // update the paths
            graph->graph.clear_path();
            graph->paths.to_graph(graph->graph);
            // and rebuild the indexes
            rebuild(graph);
            //graph->serialize_to_file(name + "-post-index.vg");

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
        << "    -x, --xg-name FILE      use the graph in this xg index" << endl
        << "    -t, --threads N         number of threads to use" << endl
        << "    -p, --into-path NAME    surject into just this path" << endl
        << "    -i, --into-paths FILE   surject into nonoverlapping path names listed in FILE (one per line)" << endl
        << "    -P, --into-prefix NAME  surject into all paths with NAME as their prefix" << endl
        //<< "    -H, --header-from FILE  use the header in the SAM/CRAM/BAM file for the output" << endl
        // todo, reenable
        // << "    -c, --cram-output       write CRAM to stdout (default is vg::Aligment/GAM format)" << endl
        // << "    -f, --reference FILE    use this file when writing CRAM to rebuild sequence header" << endl
         << "    -n, --context-depth N     expand this many steps when preparing graph for surjection (default: 3)" << endl
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

    string xg_name;
    string path_name;
    string path_prefix;
    string path_file;
    string output_type = "gam";
    string input_type = "gam";
    string header_file;
    int compress_level = 9;
    int window = 5;
    string fasta_filename;
    int context_depth = 3;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xb-name", required_argument, 0, 'x'},
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
            {"context-depth", required_argument, 0, 'n'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:p:i:P:cbsH:C:t:w:f:n:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
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

        case 'n':
            context_depth = atoi(optarg);
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

    string file_name = get_input_file_name(optind, argc, argv);

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

    xg::XG* xgidx = nullptr;
    ifstream xg_stream(xg_name);
    if(xg_stream) {
        xgidx = new xg::XG(xg_stream);
    }
    if (!xg_stream || xgidx == nullptr) {
        cerr << "[vg sim] error: could not open xg index" << endl;
        return 1;
    }

    map<string, int64_t> path_by_id;// = index.paths_by_id();
    map<string, int64_t> path_length;
    int num_paths = xgidx->max_path_rank();
    for (int i = 1; i <= num_paths; ++i) {
        auto name = xgidx->path_name(i);
        path_by_id[name] = i;
        path_length[name] = xgidx->path_length(name);
    }

    int thread_count = get_thread_count();
    vector<Mapper*> mapper;
    mapper.resize(thread_count);
    for (int i = 0; i < thread_count; ++i) {
        Mapper* m = new Mapper;
        m->xindex = xgidx;
        m->context_depth = context_depth;
        mapper[i] = m;
    }

    if (input_type == "gam") {
        if (output_type == "gam") {
            int thread_count = get_thread_count();
            vector<vector<Alignment> > buffer;
            buffer.resize(thread_count);
            function<void(Alignment&)> lambda = [&xgidx, &path_names, &buffer, &window, &mapper](Alignment& src) {
                int tid = omp_get_thread_num();
                Alignment surj;
                // Since we're outputting full GAM, we ignore all this info
                // about where on the path the alignment falls. But we need to
                // provide the space to the surject call anyway.
                string path_name;
                int64_t path_pos;
                bool path_reverse;
                buffer[tid].push_back(mapper[tid]->surject_alignment(src, path_names,path_name, path_pos, path_reverse, window));
                stream::write_buffered(cout, buffer[tid], 100);
            };
            get_input_file(file_name, [&](istream& in) {
                stream::for_each_parallel(in, lambda);
            });
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
            int thread_count = get_thread_count();
            vector<vector<tuple<string, int64_t, bool, Alignment> > > buffer;
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
                &out_mode, &out, &output_lock, &fasta_filename](vector<tuple<string, int64_t, bool, Alignment> >& buf) {
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
                                auto& path_reverse = get<2>(s);
                                auto& surj = get<3>(s);
                                string cigar = cigar_against_path(surj, path_reverse);
                                bam1_t* b = alignment_to_bam(header,
                                        surj,
                                        path_nom,
                                        path_pos,
                                        path_reverse,
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

            function<void(Alignment&)> lambda = [&xgidx,
                                                 &mapper,
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
                    string path_name;
                    int64_t path_pos;
                    bool path_reverse;
                    int tid = omp_get_thread_num();
                    auto surj = mapper[tid]->surject_alignment(src, path_names, path_name, path_pos, path_reverse, window);
                    if (!surj.path().mapping_size()) {
                        surj = src;
                    }
                    // record
                    if (!hdr && !surj.read_group().empty() && !surj.sample_name().empty()) {
#pragma omp critical (hts_header)
                        rg_sample[surj.read_group()] = surj.sample_name();
                    }

                    buffer[tid].push_back(make_tuple(path_name, path_pos, path_reverse, surj));
                    handle_buffer(buffer[tid]);

                };


            // now apply the alignment processor to the stream
            get_input_file(file_name, [&](istream& in) {
                stream::for_each_parallel(in, lambda);
            });
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

void help_sim(char** argv) {
    cerr << "usage: " << argv[0] << " sim [options]" << endl
         << "Samples sequences from the xg-indexed graph." << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE    use the xg index in FILE" << endl
         << "    -l, --read-length N   write reads of length N" << endl
         << "    -n, --num-reads N     simulate N reads" << endl
         << "    -s, --random-seed N   use this specific seed for the PRNG" << endl
         << "    -e, --base-error N    base substitution error rate (default 0.0)" << endl
         << "    -i, --indel-error N   indel error rate (default 0.0)" << endl
         << "    -f, --forward-only    don't simulate from the reverse strand" << endl
         << "    -p, --frag-len N      make paired end reads with given fragment length N" << endl
         << "    -v, --frag-std-dev N  use this standard deviation for fragment length estimation" << endl
         << "    -a, --align-out       generate true alignments on stdout rather than reads" << endl
         << "    -J, --json-out        write alignments in json" << endl;
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
    bool json_out = false;
    int fragment_length = 0;
    double fragment_std_dev = 0;
    string xg_name;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"read-length", required_argument, 0, 'l'},
            {"num-reads", required_argument, 0, 'n'},
            {"random-seed", required_argument, 0, 's'},
            {"forward-only", no_argument, 0, 'f'},
            {"align-out", no_argument, 0, 'a'},
            {"json-out", no_argument, 0, 'J'},
            {"base-error", required_argument, 0, 'e'},
            {"indel-error", required_argument, 0, 'i'},
            {"frag-len", required_argument, 0, 'p'},
            {"frag-std-dev", required_argument, 0, 'v'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:e:i:fax:Jp:v:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;

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

        case 'J':
            json_out = true;
            align_out = true;
            break;

        case 'p':
            fragment_length = atoi(optarg);
            break;

        case 'v':
            fragment_std_dev = atof(optarg);
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

    if (xg_name.empty()) {
        cerr << "[vg sim] error: we need an xg index to sample reads from" << endl;
        return 1;
    }

    mt19937 rng;
    rng.seed(seed_val);

    xg::XG* xgidx = nullptr;
    ifstream xg_stream(xg_name);
    if(xg_stream) {
        xgidx = new xg::XG(xg_stream);
    }
    if (!xg_stream || xgidx == nullptr) {
        cerr << "[vg sim] error: could not open xg index" << endl;
        return 1;
    }

    Sampler sampler(xgidx, seed_val, forward_only);
    size_t max_iter = 1000;
    int nonce = 1;
    for (int i = 0; i < num_reads; ++i) {
        if (fragment_length) {
            auto alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
            size_t iter = 0;
            while (iter++ < max_iter) {
                if (alns.front().sequence().size() < read_length
                    || alns.back().sequence().size() < read_length) {
                    alns = sampler.alignment_pair(read_length, fragment_length, fragment_std_dev, base_error, indel_error);
                }
            }
            // write the alignment or its string
            if (align_out) {
                // write it out as requested
                if (json_out) {
                    cout << pb2json(alns.front()) << endl;
                    cout << pb2json(alns.back()) << endl;
                } else {
                    function<Alignment(uint64_t)> lambda = [&alns](uint64_t n) { return alns[n]; };
                    stream::write(cout, 2, lambda);
                }
            } else {
                cout << alns.front().sequence() << "\t" << alns.back().sequence() << endl;
            }
        } else {
            auto aln = sampler.alignment_with_error(read_length, base_error, indel_error);
            size_t iter = 0;
            while (iter++ < max_iter) {
                if (aln.sequence().size() < read_length) {
                    auto aln_prime = sampler.alignment_with_error(read_length, base_error, indel_error);
                    if (aln_prime.sequence().size() > aln.sequence().size()) {
                        aln = aln_prime;
                    }
                }
            }
            // write the alignment or its string
            if (align_out) {
                // write it out as requested
                if (json_out) {
                    cout << pb2json(aln) << endl;
                } else {
                    function<Alignment(uint64_t)> lambda = [&aln](uint64_t n) { return aln; };
                    stream::write(cout, 1, lambda);
                }
            } else {
                cout << aln.sequence() << endl;
            }
        }
    }

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
         << "    -b, --superbubbles    describe the superbubbles of the graph" << endl
         << "    -u, --ultrabubbles    describe the ultrabubbles of the graph" << endl
         << "    -c, --components      print the strongly connected components of the graph" << endl
         << "    -A, --is-acyclic      print if the graph is acyclic or not" << endl
         << "    -n, --node ID         consider node with the given id" << endl
         << "    -d, --to-head         show distance to head for each provided node" << endl
         << "    -t, --to-tail         show distance to head for each provided node" << endl
         << "    -a, --alignments FILE compute stats for reads aligned to the graph" << endl
         << "    -v, --verbose         output longer reports" << endl;
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
    bool superbubbles = false;
    bool ultrabubbles = false;
    bool verbose = false;
    bool is_acyclic = false;
    set<vg::id_t> ids;
    // What alignments GAM file should we read and compute stats on with the
    // graph?
    string alignments_filename;

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
            {"superbubbles", no_argument, 0, 'b'},
            {"ultrabubbles", no_argument, 0, 'u'},
            {"alignments", required_argument, 0, 'a'},
            {"is-acyclic", no_argument, 0, 'A'},
            {"verbose", no_argument, 0, 'v'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hzlsHTScdtn:NEbua:vA",
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

        case 'b':
            superbubbles = true;
            break;

        case 'u':
            ultrabubbles = true;
            break;

        case 'A':
            is_acyclic = true;
            break;

        case 'a':
            alignments_filename = optarg;
            break;

        case 'v':
            verbose = true;
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
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

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

    if (superbubbles || ultrabubbles) {
        auto bubbles = superbubbles ? vg::superbubbles(*graph) : vg::ultrabubbles(*graph);
        for (auto& i : bubbles) {
            auto b = i.first;
            auto v = i.second;
            // sort output for now, to help do diffs in testing
            sort(v.begin(), v.end());
            cout << b.first << "\t" << b.second << "\t";
            for (auto& n : v) {
                cout << n << ",";
            }
            cout << endl;
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

    if (is_acyclic) {
        if (graph->is_acyclic()) {
            cout << "acyclic" << endl;
        } else {
            cout << "cyclic" << endl;
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

    if (!alignments_filename.empty()) {
        // Read in the given GAM
        ifstream alignment_stream(alignments_filename);

        // We need some allele parsing functions

        // This one decided if a path is really an allele path
        auto path_name_is_allele = [](const string path_name) -> bool {
            string prefix = "_alt_";
            // It needs to start with "_alt_" and have another separating
            // underscore between site name and allele number
            return(prefix.size() < path_name.size() &&
                count(path_name.begin(), path_name.end(), '_') >= 3 &&
                equal(prefix.begin(), prefix.end(), path_name.begin()));
        };

        // This one gets the site name from an allele path name
        auto path_name_to_site = [](const string& path_name) -> string {
            auto last_underscore = path_name.rfind('_');
            assert(last_underscore != string::npos);
            return path_name.substr(0, last_underscore);
        };

        // This one gets the allele name from an allele path name
        auto path_name_to_allele = [](const string& path_name) -> string {
            auto last_underscore = path_name.rfind('_');
            assert(last_underscore != string::npos);
            return path_name.substr(last_underscore + 1);
        };

        // Before we go over the reads, we need to make a map that tells us what
        // nodes are unique to what allele paths. Stores site and allele parts
        // separately.
        map<vg::id_t, pair<string, string>> allele_path_for_node;

        // This is what we really care about: for each pair of allele paths in
        // the graph, we need to find out whether the coverage imbalance between
        // them among primary alignments is statistically significant. For this,
        // we need to track how many reads overlap the distinct parts of allele
        // paths.

        // This is going to be indexed by site
        // ("_alt_f6d951572f9c664d5d388375aa8b018492224533") and then by allele
        // ("0"). A read only counts if it visits a node that's on one allele
        // and not any others in that site.

        // We need to pre-populate it with 0s so we know which sites actually
        // have 2 alleles and which only have 1 in the graph.
        map<string, map<string, size_t>> reads_on_allele;

        graph->for_each_node_parallel([&](Node* node) {
            // For every node

            if(!graph->paths.has_node_mapping(node)) {
                // No paths to go over. If we try and get them we'll be
                // modifying the paths in parallel, which will explode.
                return;
            }

            // We want an allele path on it
            string allele_path;
            for(auto& name_and_mappings : graph->paths.get_node_mapping(node)) {
                // For each path on it
                if(path_name_is_allele(name_and_mappings.first)) {
                    // If it's an allele path
                    if(allele_path.empty()) {
                        // It's the first. Take it.
                        allele_path = name_and_mappings.first;
                    } else {
                        // It's a subsequent one. This node is not uniquely part
                        // of any allele path.
                        return;
                    }
                }
            }

            if(!allele_path.empty()) {
                // We found an allele path for this node

                // Get its site and allele so we can count it as a biallelic
                // site. Note that sites where an allele has no unique nodes
                // (pure indels, for example) can't be handled and will be
                // ignored.
                auto site = path_name_to_site(allele_path);
                auto allele = path_name_to_allele(allele_path);


                #pragma omp critical (allele_path_for_node)
                allele_path_for_node[node->id()] = make_pair(site, allele);

                #pragma omp critical (reads_on_allele)
                reads_on_allele[site][allele] = 0;
            }
        });


        // These are the general stats we will compute.
        size_t total_alignments = 0;
        size_t total_aligned = 0;
        size_t total_primary = 0;
        size_t total_secondary = 0;

        // These are for counting significantly allele-biased hets
        size_t total_hets = 0;
        size_t significantly_biased_hets = 0;

        // These are for tracking which nodes are covered and which are not
        map<vg::id_t, size_t> node_visit_counts;

        // And for counting indels
        // Inserted bases also counts softclips
        size_t total_insertions = 0;
        size_t total_inserted_bases = 0;
        size_t total_deletions = 0;
        size_t total_deleted_bases = 0;
        // And substitutions
        size_t total_substitutions = 0;
        size_t total_substituted_bases = 0;
        // And softclips
        size_t total_softclips = 0;
        size_t total_softclipped_bases = 0;

        // In verbose mode we want to report details of insertions, deletions,
        // and substitutions, and soft clips.
        vector<pair<vg::id_t, Edit>> insertions;
        vector<pair<vg::id_t, Edit>> deletions;
        vector<pair<vg::id_t, Edit>> substitutions;
        vector<pair<vg::id_t, Edit>> softclips;

        function<void(Alignment&)> lambda = [&](Alignment& aln) {
            int tid = omp_get_thread_num();

            // We ought to be able to do many stats on the alignments.

            // Now do all the non-mapping stats
            #pragma omp critical (total_alignments)
            total_alignments++;
            if(aln.is_secondary()) {
                #pragma omp critical (total_secondary)
                total_secondary++;
            } else {
                #pragma omp critical (total_primary)
                total_primary++;
                if(aln.score() > 0) {
                    // We only count aligned primary reads in "total aligned";
                    // the primary can't be unaligned if the secondary is
                    // aligned.
                    #pragma omp critical (total_aligned)
                    total_aligned++;
                }

                // Which sites and alleles does this read support. TODO: if we hit
                // unique nodes from multiple alleles of the same site, we should...
                // do something. Discard the read? Not just count it on both sides
                // like we do now.
                set<pair<string, string>> alleles_supported;

                for(size_t i = 0; i < aln.path().mapping_size(); i++) {
                    // For every mapping...
                    auto& mapping = aln.path().mapping(i);
                    vg::id_t node_id = mapping.position().node_id();

                    if(allele_path_for_node.count(node_id)) {
                        // We hit a unique node for this allele. Add it to the set,
                        // in case we hit another unique node for it later in the
                        // read.
                        alleles_supported.insert(allele_path_for_node.at(node_id));
                    }

                    // Record that there was a visit to this node.
                    #pragma omp critical (node_visit_counts)
                    node_visit_counts[node_id]++;

                    for(size_t j = 0; j < mapping.edit_size(); j++) {
                        // Go through edits and look for each type.
                        auto& edit = mapping.edit(j);

                        if(edit.to_length() > edit.from_length()) {
                            if((j == 0 && i == 0) || (j == mapping.edit_size() - 1 && i == aln.path().mapping_size() - 1)) {
                                // We're at the very end of the path, so this is a soft clip.
                                #pragma omp critical (total_softclipped_bases)
                                total_softclipped_bases += edit.to_length() - edit.from_length();
                                #pragma omp critical (total_softclips)
                                total_softclips++;
                                if(verbose) {
                                    // Record the actual insertion
                                    #pragma omp critical (softclips)
                                    softclips.push_back(make_pair(node_id, edit));
                                }
                            } else {
                                // Record this insertion
                                #pragma omp critical (total_inserted_bases)
                                total_inserted_bases += edit.to_length() - edit.from_length();
                                #pragma omp critical (total_insertions)
                                total_insertions++;
                                if(verbose) {
                                    // Record the actual insertion
                                    #pragma omp critical (insertions)
                                    insertions.push_back(make_pair(node_id, edit));
                                }
                            }

                        } else if(edit.from_length() > edit.to_length()) {
                            // Record this deletion
                            #pragma omp critical (total_deleted_bases)
                            total_deleted_bases += edit.from_length() - edit.to_length();
                            #pragma omp critical (total_deletions)
                            total_deletions++;
                            if(verbose) {
                                // Record the actual deletion
                                #pragma omp critical (deletions)
                                deletions.push_back(make_pair(node_id, edit));
                            }
                        } else if(!edit.sequence().empty()) {
                            // Record this substitution
                            // TODO: a substitution might also occur as part of a deletion/insertion above!
                            #pragma omp critical (total_substituted_bases)
                            total_substituted_bases += edit.from_length();
                            #pragma omp critical (total_substitutions)
                            total_substitutions++;
                            if(verbose) {
                                // Record the actual substitution
                                #pragma omp critical (substitutions)
                                substitutions.push_back(make_pair(node_id, edit));
                            }
                        }

                    }
                }

                for(auto& site_and_allele : alleles_supported) {
                    // This read is informative for an allele of a site.
                    // Up the reads on that allele of that site.
                    #pragma omp critical (reads_on_allele)
                    reads_on_allele[site_and_allele.first][site_and_allele.second]++;
                }
            }

        };

        // Actually go through all the reads and count stuff up.
        stream::for_each_parallel(alignment_stream, lambda);

        // Calculate stats about the reads per allele data
        for(auto& site_and_alleles : reads_on_allele) {
            // For every site
            if(site_and_alleles.second.size() == 2) {
                // If it actually has 2 alleles with unique nodes in the
                // graph (so we can use the binomial)

                // We'll fill this with the counts for the two present alleles.
                vector<size_t> counts;

                for(auto& allele_and_count : site_and_alleles.second) {
                    // Collect all the counts
                    counts.push_back(allele_and_count.second);
                }

                if(counts[0] > counts[1]) {
                    // We have a 50% underlying probability so we can just put
                    // the rarer allele first.
                    swap(counts[0], counts[1]);
                }

                // What's the log prob for the smaller tail?
                auto tail_logprob = binomial_cmf_ln(prob_to_logprob(0.5),  counts[1] + counts[0], counts[0]);

                // Double it to get the two-tailed test
                tail_logprob += prob_to_logprob(2);

#ifdef debug
                cerr << "Site " << site_and_alleles.first << " has " << counts[0]
                    << " and " << counts[1] << " p=" << logprob_to_prob(tail_logprob) << endl;
#endif

                if(tail_logprob < prob_to_logprob(0.05)) {
                    significantly_biased_hets++;
                }
                total_hets++;

            }
        }

        // Go through all the nodes again and sum up unvisited nodes
        size_t unvisited_nodes = 0;
        // And unvisited base count
        size_t unvisited_node_bases = 0;
        // And nodes that are visited by only one thing (which is useful if
        // we're checking diploid assembly pairs).
        size_t single_visited_nodes = 0;
        size_t single_visited_node_bases = 0;
        // If we're in verbose mode, collect IDs too.
        set<vg::id_t> unvisited_ids;
        set<vg::id_t> single_visited_ids;
        // Note that you need to subtract out substituted-away and deleted bases
        // from the sum of 2 * double- and single-visited bases to get the bases
        // actually present in reads, because deleted bases are still "visited"
        // as many times as their nodes are touched. Also note that we ignore
        // edge effects and a read that stops before the end of a node will
        // visit the whole node.
        graph->for_each_node_parallel([&](Node* node) {
            // For every node
            if(!node_visit_counts.count(node->id()) || node_visit_counts.at(node->id()) == 0) {
                // If we never visited it with a read, count it.
                #pragma omp critical (unvisited_nodes)
                unvisited_nodes++;
                #pragma omp critical (unvisited_node_bases)
                unvisited_node_bases += node->sequence().size();
                if(verbose) {
                    #pragma omp critical (unvisited_ids)
                    unvisited_ids.insert(node->id());
                }
            } else if(node_visit_counts.at(node->id()) == 1) {
                // If we visited it with only one read, count it.
                #pragma omp critical (single_visited_nodes)
                single_visited_nodes++;
                #pragma omp critical (single_visited_node_bases)
                single_visited_node_bases += node->sequence().size();
                if(verbose) {
                    #pragma omp critical (single_visited_ids)
                    single_visited_ids.insert(node->id());
                }
            }
        });

        cout << "Total alignments: " << total_alignments << endl;
        cout << "Total primary: " << total_primary << endl;
        cout << "Total secondary: " << total_secondary << endl;
        cout << "Total aligned: " << total_aligned << endl;

        cout << "Insertions: " << total_inserted_bases << " bp in " << total_insertions << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : insertions) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.sequence()
                    << " on " << id_and_edit.first << endl;
            }
        }
        cout << "Deletions: " << total_deleted_bases << " bp in " << total_deletions << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : deletions) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.to_length()
                    << " on " << id_and_edit.first << endl;
            }
        }
        cout << "Substitutions: " << total_substituted_bases << " bp in " << total_substitutions << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : substitutions) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.sequence()
                    << " on " << id_and_edit.first << endl;
            }
        }
        cout << "Softclips: " << total_softclipped_bases << " bp in " << total_softclips << " read events" << endl;
        if(verbose) {
            for(auto& id_and_edit : softclips) {
                cout << "\t" << id_and_edit.second.from_length() << " -> " << id_and_edit.second.sequence()
                    << " on " << id_and_edit.first << endl;
            }
        }

        cout << "Unvisited nodes: " << unvisited_nodes << "/" << graph->node_count()
            << " (" << unvisited_node_bases << " bp)" << endl;
        if(verbose) {
            for(auto& id : unvisited_ids) {
                cout << "\t" << id << endl;
            }
        }

        cout << "Single-visited nodes: " << single_visited_nodes << "/" << graph->node_count()
            << " (" << single_visited_node_bases << " bp)" << endl;
        if(verbose) {
            for(auto& id : single_visited_ids) {
                cout << "\t" << id << endl;
            }
        }

        cout << "Significantly biased heterozygous sites: " << significantly_biased_hets << "/" << total_hets;
        if(total_hets > 0) {
            cout << " (" << (double)significantly_biased_hets / total_hets * 100 << "%)";
        }
        cout << endl;


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
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

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
    cerr << "usage: " << argv[0] << " find [options] >sub.vg" << endl
         << "options:" << endl
         << "    -d, --db-name DIR      use this db (defaults to <graph>.index/)" << endl
         << "    -x, --xg-name FILE     use this xg index (instead of rocksdb db)" << endl
         << "graph features:" << endl
         << "    -n, --node ID          find node, return 1-hop context as graph" << endl
         << "    -N, --node-list FILE   a white space or line delimited list of nodes to collect" << endl
         << "    -e, --edges-end ID     return edges on end of node with ID" << endl
         << "    -s, --edges-start ID   return edges on start of node with ID" << endl
         << "    -c, --context STEPS    expand the context of the subgraph this many steps" << endl
         << "    -L, --use-length       treat STEPS in -c or M in -r as a length in bases" << endl
         << "    -p, --path TARGET      find the node(s) in the specified path range TARGET=path[:pos1[-pos2]]" << endl
         << "    -P, --position-in PATH find the position of the node (specified by -n) in the given path" << endl
         << "    -r, --node-range N:M   get nodes from N to M" << endl
         << "    -G, --gam GAM          accumulate the graph touched by the alignments in the GAM" << endl
         << "alignments: (rocksdb only)" << endl
         << "    -a, --alignments       writes alignments from index, sorted by node id" << endl
         << "    -i, --alns-in N:M      writes alignments whose start nodes is between N and M (inclusive)" << endl
         << "    -o, --alns-on N:M      writes alignments which align to any of the nodes between N and M (inclusive)" << endl
         << "    -A, --to-graph VG      get alignments to the provided subgraph" << endl
         << "sequences:" << endl
         << "    -g, --gcsa FILE        use this GCSA2 index of the sequence space of the graph" << endl
         << "    -z, --kmer-size N      split up --sequence into kmers of size N" << endl
         << "    -j, --kmer-stride N    step distance between succesive kmers in sequence (default 1)" << endl
         << "    -S, --sequence STR     search for sequence STR using --kmer-size kmers" << endl
         << "    -M, --mems STR         describe the super-maximal exact matches of the STR (gcsa2) in JSON" << endl
         << "    -Y, --max-mem N        the maximum length of the MEM (default: GCSA2 order)" << endl
         << "    -k, --kmer STR         return a graph of edges and nodes matching this kmer" << endl
         << "    -T, --table            instead of a graph, return a table of kmers" << endl
         << "                           (works only with kmers in the index)" << endl
         << "    -C, --kmer-count       report approximate count of kmer (-k) in db" << endl
         << "    -D, --distance         return distance on path between pair of nodes (-n). if -P not used, best path chosen heurstically" << endl
         << "haplotypes:" << endl
         << "    -H, --haplotypes FILE  count xg threads in agreement with alignments in the GAM" << endl;

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
    vector<vg::id_t> node_ids;
    string node_list_file;
    int context_size=0;
    bool use_length = false;
    bool count_kmers = false;
    bool kmer_table = false;
    string target;
    string path_name;
    string range;
    string gcsa_in;
    string xg_name;
    bool get_mems = false;
    bool get_alignments = false;
    bool get_mappings = false;
    string node_id_range;
    string aln_on_id_range;
    vg::id_t start_id = 0;
    vg::id_t end_id = 0;
    bool pairwise_distance = false;
    string haplotype_alignments;
    string gam_file;
    int max_mem_length = 0;
    string to_graph_file;

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
                {"node-list", required_argument, 0, 'N'},
                {"edges-end", required_argument, 0, 'e'},
                {"edges-start", required_argument, 0, 's'},
                {"kmer", required_argument, 0, 'k'},
                {"table", no_argument, 0, 'T'},
                {"sequence", required_argument, 0, 'S'},
                {"mems", required_argument, 0, 'M'},
                {"kmer-stride", required_argument, 0, 'j'},
                {"kmer-size", required_argument, 0, 'z'},
                {"context", required_argument, 0, 'c'},
                {"use-length", no_argument, 0, 'L'},
                {"kmer-count", no_argument, 0, 'C'},
                {"path", required_argument, 0, 'p'},
                {"position-in", required_argument, 0, 'P'},
                {"node-range", required_argument, 0, 'r'},
                {"alignments", no_argument, 0, 'a'},
                {"mappings", no_argument, 0, 'm'},
                {"alns-in", required_argument, 0, 'i'},
                {"alns-on", required_argument, 0, 'o'},
                {"distance", no_argument, 0, 'D'},
                {"haplotypes", required_argument, 0, 'H'},
                {"gam", required_argument, 0, 'G'},
                {"to-graph", required_argument, 0, 'A'},
                {"max-mem", required_argument, 0, 'Y'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "d:x:n:e:s:o:k:hc:LS:z:j:CTp:P:r:amg:M:i:DH:G:N:A:Y:",
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

        case 'Y':
            max_mem_length = atoi(optarg);
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

        case 'L':
            use_length = true;
            break;

        case 'n':
            node_ids.push_back(atoi(optarg));
            break;

        case 'N':
            node_list_file = optarg;
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

        case 'i':
            node_id_range = optarg;
            break;

        case 'm':
            get_mappings = true;
            break;

        case 'o':
            aln_on_id_range = optarg;
            break;

        case 'D':
            pairwise_distance = true;
            break;

        case 'H':
            haplotype_alignments = optarg;
            break;

        case 'G':
            gam_file = optarg;
            break;

        case 'A':
            to_graph_file = optarg;
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
        cerr << "[vg find] find does not accept positional arguments" << endl;
        return 1;
    }

    if (db_name.empty() && gcsa_in.empty() && xg_name.empty()) {
        cerr << "[vg find] find requires -d, -g, or -x to know where to find its database" << endl;
        return 1;
    }

    if (context_size > 0 && use_length == true && xg_name.empty()) {
        cerr << "[vg find] error, -L not supported without -x" << endl;
        exit(1);
    }

    // process input node list
    if (!node_list_file.empty()) {
        ifstream nli;
        nli.open(node_list_file);
        if (!nli.good()){
            cerr << "[vg find] error, unable to open the node list input file." << endl;
            exit(1);
        }
        string line;
        while (getline(nli, line)){
            for (auto& idstr : split_delims(line, " \t")) {
                node_ids.push_back(atol(idstr.c_str()));
            }
        }
        nli.close();
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
        assert(!db_name.empty());
        vector<Alignment> output_buf;
        auto lambda = [&output_buf](const Alignment& aln) {
            output_buf.push_back(aln);
            stream::write_buffered(cout, output_buf, 100);
        };
        vindex->for_each_alignment(lambda);
        stream::write_buffered(cout, output_buf, 0);
    }

    if (!node_id_range.empty()) {
        assert(!db_name.empty());
        vector<string> parts = split_delims(node_id_range, ":");
        if (parts.size() == 1) {
            convert(parts.front(), start_id);
            end_id = start_id;
        } else {
            convert(parts.front(), start_id);
            convert(parts.back(), end_id);
        }
        vector<Alignment> output_buf;
        auto lambda = [&output_buf](const Alignment& aln) {
            output_buf.push_back(aln);
            stream::write_buffered(cout, output_buf, 100);
        };
        vindex->for_alignment_in_range(start_id, end_id, lambda);
        stream::write_buffered(cout, output_buf, 0);
    }

    if (!aln_on_id_range.empty()) {
        assert(!db_name.empty());
        vector<string> parts = split_delims(aln_on_id_range, ":");
        if (parts.size() == 1) {
            convert(parts.front(), start_id);
            end_id = start_id;
        } else {
            convert(parts.front(), start_id);
            convert(parts.back(), end_id);
        }
        vector<vg::id_t> ids;
        for (auto i = start_id; i <= end_id; ++i) {
            ids.push_back(i);
        }
        vector<Alignment> output_buf;
        auto lambda = [&output_buf](const Alignment& aln) {
            output_buf.push_back(aln);
            stream::write_buffered(cout, output_buf, 100);
        };
        vindex->for_alignment_to_nodes(ids, lambda);
        stream::write_buffered(cout, output_buf, 0);
    }

    if (!to_graph_file.empty()) {
        assert(vindex != nullptr);
        ifstream tgi(to_graph_file);
        VG graph(tgi);
        vector<vg::id_t> ids;
        graph.for_each_node([&](Node* n) { ids.push_back(n->id()); });
        vector<Alignment> output_buf;
        auto lambda = [&output_buf](const Alignment& aln) {
            output_buf.push_back(aln);
            stream::write_buffered(cout, output_buf, 100);
        };
        vindex->for_alignment_to_nodes(ids, lambda);
        stream::write_buffered(cout, output_buf, 0);
    }

    if (!xg_name.empty()) {
        if (!node_ids.empty() && path_name.empty() && !pairwise_distance) {
            // get the context of the node
            vector<Graph> graphs;
            set<vg::id_t> ids;
            for (auto node_id : node_ids) ids.insert(node_id);
            for (auto node_id : node_ids) {
                Graph g;
                xindex.neighborhood(node_id, context_size, g, !use_length);
                if (context_size == 0) {
                    for (auto& edge : xindex.edges_of(node_id)) {
                        // if both ends of the edge are in our targets, keep them
                        if (ids.count(edge.to()) && ids.count(edge.from())) {
                            *g.add_edge() = edge;
                        }
                    }
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
            for (auto& e : xindex.edges_on_end(end_id)) {
                cout << (e.from_start() ? -1 : 1) * e.from() << "\t" <<  (e.to_end() ? -1 : 1) * e.to() << endl;
            }
        } else if (start_id != 0) {
            for (auto& e : xindex.edges_on_start(start_id)) {
                cout << (e.from_start() ? -1 : 1) * e.from() << "\t" <<  (e.to_end() ? -1 : 1) * e.to() << endl;
            }
        }
        if (!node_ids.empty() && !path_name.empty() && !pairwise_distance) {
            // Note: this isn't at all consistent with -P option with rocksdb, which couts a range
            // and then mapping, but need this info right now for scripts/chunked_call
            for (auto node_id : node_ids) {
                cout << node_id;
                vector<size_t> positions = xindex.position_in_path(node_id, path_name);
                for (auto pos : positions) {
                    cout << "\t" << pos;
                }
                cout << endl;
            }
        }
        if (pairwise_distance) {
            if (node_ids.size() != 2) {
                cerr << "[vg find] error, exactly 2 nodes (-n) required with -D" << endl;
                exit(1);
            }
            if (!path_name.empty()) {
                cout << xindex.approx_path_distance(path_name, node_ids[0], node_ids[1]) << endl;
            } else {
                cout << xindex.min_approx_path_distance(vector<string>(), node_ids[0], node_ids[1]) << endl;
            }
            return 0;
        }
        if (!target.empty()) {
            string name;
            int64_t start, end;
            Graph graph;
            parse_region(target, name, start, end);
            if(xindex.path_rank(name) == 0) {
                // Passing a nonexistent path to get_path_range produces Undefined Behavior
                cerr << "[vg find] error, path " << name << " not found in index" << endl;
                exit(1);
            }
            xindex.get_path_range(name, start, end, graph);
            if (context_size > 0) {
                xindex.expand_context(graph, context_size, true, !use_length);
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
            if (!use_length) {
                xindex.get_id_range(id_start, id_end, graph);
            } else {
                // treat id_end as length instead.
                xindex.get_id_range_by_length(id_start, id_end, graph, true);
            }
            if (context_size > 0) {
                xindex.expand_context(graph, context_size, true, !use_length);
            }
            VG vgg; vgg.extend(graph); // removes dupes
            vgg.remove_orphan_edges();
            vgg.serialize_to_ostream(cout);
        }
        if(!haplotype_alignments.empty()) {
            // What should we do with each alignment?
            function<void(Alignment&)> lambda = [&xindex](Alignment& aln) {
                // Count the amtches to the path. The path might be empty, in
                // which case it will yield the biggest size_t you can have.
                size_t matches = xindex.count_matches(aln.path());

                // We do this single-threaded, at least for now, so we don't
                // need to worry about coordinating output, and we can just
                // spit out the counts as bare numbers.
                cout << matches << endl;
            };
            if (haplotype_alignments == "-") {
                stream::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(haplotype_alignments.c_str());
                if(!in.is_open()) {
                    cerr << "[vg find] error: could not open alignments file " << haplotype_alignments << endl;
                    exit(1);
                }
                stream::for_each(in, lambda);
            }

        }
        if (!gam_file.empty()) {
            set<vg::id_t> nodes;
            function<void(Alignment&)> lambda = [&nodes](Alignment& aln) {
                // accumulate nodes matched by the path
                auto& path = aln.path();
                for (int i = 0; i < path.mapping_size(); ++i) {
                    nodes.insert(path.mapping(i).position().node_id());
                }
            };
            if (gam_file == "-") {
                stream::for_each(std::cin, lambda);
            } else {
                ifstream in;
                in.open(gam_file.c_str());
                if(!in.is_open()) {
                    cerr << "[vg find] error: could not open alignments file " << gam_file << endl;
                    exit(1);
                }
                stream::for_each(in, lambda);
            }
            // now we have the nodes to get
            VG graph;
            for (auto& node : nodes) {
                *graph.graph.add_node() = xindex.node(node);
            }
            xindex.expand_context(graph.graph, max(1, context_size)); // get connected edges
            graph.rebuild_indexes();
            graph.serialize_to_ostream(cout);
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

            // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);

            // Open it
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
                auto mems = mapper.find_mems(sequence.begin(), sequence.end(), max_mem_length);
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

void help_align(char** argv) {
    cerr << "usage: " << argv[0] << " align [options] <graph.vg> >alignments.gam" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -Q, --seq-name STR    name the sequence using this value" << endl
         << "    -j, --json            output alignments in JSON format (default GAM)" << endl
         << "    -m, --match N         use this match score (default: 1)" << endl
         << "    -M, --mismatch N      use this mismatch penalty (default: 4)" << endl
         << "    -g, --gap-open N      use this gap open penalty (default: 6)" << endl
         << "    -e, --gap-extend N    use this gap extension penalty (default: 1)" << endl
         << "    -T, --full-l-bonus N  provide this bonus for alignments that are full length (default: 5)" << endl
         << "    -b, --banded-global   use the banded global alignment algorithm" << endl
         << "    -p, --pinned          pin the (local) alignment traceback to the optimal edge of the graph" << endl
         << "    -L, --pin-left        pin the first rather than last bases of the graph and sequence" << endl
         << "    -D, --debug           print out score matrices and other debugging info" << endl
         << "options:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -Q, --seq-name STR    name the sequence using this value" << endl
         << "    -r, --reference STR   don't use an input graph--- run SSW alignment between -s and -r" << endl
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
    int match = 1;
    int mismatch = 4;
    int gap_open = 6;
    int gap_extend = 1;
    int full_length_bonus = 5;
    string ref_seq;
    bool debug = false;
    bool banded_global = false;
    bool pinned_alignment = false;
    bool pin_left = false;

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
            {"match", required_argument, 0, 'm'},
            {"mismatch", required_argument, 0, 'M'},
            {"gap-open", required_argument, 0, 'g'},
            {"gap-extend", required_argument, 0, 'e'},
            {"reference", required_argument, 0, 'r'},
            {"debug", no_argument, 0, 'D'},
            {"banded-global", no_argument, 0, 'b'},
            {"full-l-bonus", required_argument, 0, 'T'},
            {"pinned", no_argument, 0, 'p'},
            {"pinned-left", no_argument, 0, 'L'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:jhQ:m:M:g:e:Dr:F:O:bT:pL",
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

        case 'm':
            match = atoi(optarg);
            break;

        case 'M':
            mismatch = atoi(optarg);
            break;

        case 'g':
            gap_open = atoi(optarg);
            break;

        case 'e':
            gap_extend = atoi(optarg);
            break;

        case 'T':
            full_length_bonus = atoi(optarg);
            break;

        case 'r':
            ref_seq = optarg;
            break;

        case 'D':
            debug = true;
            break;

        case 'b':
            banded_global = true;
            break;

        case 'p':
            pinned_alignment = true;
            break;

        case 'L':
            pin_left = true;
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

    VG* graph = nullptr;
    if (ref_seq.empty()) {
        // Only look at a filename if we don't have an explicit reference
        // sequence.
        get_input_file(optind, argc, argv, [&](istream& in) {
            graph = new VG(in);
        });
    }
    
    Alignment alignment;
    if (!ref_seq.empty()) {
        SSWAligner ssw = SSWAligner(match, mismatch, gap_open, gap_extend);
        alignment = ssw.align(seq, ref_seq);
    } else {
        Aligner aligner = Aligner(match, mismatch, gap_open, gap_extend);
        alignment = graph->align(seq, &aligner, 0, pinned_alignment, pin_left, full_length_bonus, banded_global, debug);
    }

    if (!seq_name.empty()) {
        alignment.set_name(seq_name);
    }

    if (output_json) {
        cout << pb2json(alignment) << endl;
    } else {
        function<Alignment(uint64_t)> lambda =
            [&alignment] (uint64_t n) {
                return alignment;
            };
        stream::write(cout, 1, lambda);
    }

    if (graph != nullptr) {
        delete graph;
    }

    return 0;

}

void help_map(char** argv) {
    cerr << "usage: " << argv[0] << " map [options] <graph.vg> >alignments.vga" << endl
         << "options:" << endl
         << "    -d, --db-name DIR     use this db (defaults to <graph>.vg.index/)" << endl
         << "                          A graph is not required. But GCSA/xg take precedence if available." << endl
         << "    -x, --xg-name FILE    use this xg index (defaults to <graph>.vg.xg)" << endl
         << "    -g, --gcsa-name FILE  use this GCSA2 index (defaults to <graph>" << gcsa::GCSA::EXTENSION << ")" << endl
         << "input:" << endl
         << "    -s, --sequence STR    align a string to the graph in graph.vg using partial order alignment" << endl
         << "    -I, --quality STR     Phred+33 base quality of sequence (for base quality adjusted alignment)" << endl
         << "    -Q, --seq-name STR    name the sequence using this value (for graph modification with new named paths)" << endl
         << "    -r, --reads FILE      take reads (one per line) from FILE, write alignments to stdout" << endl
         << "    -b, --hts-input FILE  align reads from htslib-compatible FILE (BAM/CRAM/SAM) stdin (-), alignments to stdout" << endl
         << "    -G, --gam-input FILE  realign GAM input" << endl
         << "    -K, --keep-secondary  produce alignments for secondary input alignments in addition to primary ones" << endl
         << "    -f, --fastq FILE      input fastq (possibly compressed), two are allowed, one for each mate" << endl
         << "    -i, --interleaved     fastq is interleaved paired-ended" << endl
         << "    -N, --sample NAME     for --reads input, add this sample" << endl
         << "    -R, --read-group NAME for --reads input, add this read group" << endl
         << "output:" << endl
         << "    -J, --output-json     output JSON rather than an alignment stream (helpful for debugging)" << endl
         << "    -Z, --buffer-size N   buffer this many alignments together before outputting in GAM (default: 100)" << endl
         << "    -w, --compare         consider GAM input (-G) as thruth, table of name, overlap with truth, identity, score, mapqual" << endl
         << "    -D, --debug           print debugging information about alignment to stderr" << endl
         << "local alignment parameters:" << endl
         << "    -q, --match N         use this match score (default: 1)" << endl
         << "    -z, --mismatch N      use this mismatch penalty (default: 4)" << endl
         << "    -o, --gap-open N      use this gap open penalty (default: 6)" << endl
         << "    -y, --gap-extend N    use this gap extension penalty (default: 1)" << endl
         << "    -T, --full-l-bonus N  the full-length alignment bonus (default: 5)" << endl
         << "    -1, --qual-adjust     perform base quality adjusted alignments (requires base quality input)" << endl
         << "paired end alignment parameters:" << endl
         << "    -W, --fragment-max N       maximum fragment size to be used for estimating the fragment length distribution (default: 1e4)" << endl
         << "    -2, --fragment-sigma N     calculate fragment size as mean(buf)+sd(buf)*N where buf is the buffer of perfect pairs we use (default: 10)" << endl
         << "    -p, --pair-window N        maximum distance between properly paired reads in node ID space" << endl
         << "    -u, --extra-multimaps N    examine N extra mappings looking for a consistent read pairing (default: 2)" << endl
         << "    -U, --always-rescue        rescue each imperfectly-mapped read in a pair off the other" << endl
         << "    -O, --top-pairs-only       only produce paired alignments if both sides of the pair are top-scoring individually" << endl
         << "generic mapping parameters:" << endl
         << "    -B, --band-width N        for very long sequences, align in chunks then merge paths, no mapping quality (default 1000bp)" << endl
         << "    -P, --min-identity N      accept alignment only if the alignment identity to ref is >= N (default: 0)" << endl
         << "    -n, --context-depth N     follow this many edges out from each thread for alignment (default: 7)" << endl
         << "    -M, --max-multimaps N     produce up to N alignments for each read (default: 1)" << endl
         << "    -3, --softclip-trig N     trigger graph extension and realignment when either end has softclips (default: 0)" << endl
         << "    -m, --hit-max N           ignore kmers or MEMs who have >N hits in our index (default: 100)" << endl
         << "    -c, --clusters N          use at most the largest N ordered clusters of the kmer graph for alignment (default: all)" << endl
         << "    -C, --cluster-min N       require at least this many kmer hits in a cluster to attempt alignment (default: 1)" << endl
         << "    -H, --max-target-x N      skip cluster subgraphs with length > N*read_length (default: 100; unset: 0)" << endl
         << "    -e, --thread-ex N         grab this many nodes in id space around each thread for alignment (default: 10)" << endl
         << "    -t, --threads N           number of threads to use" << endl
         << "    -X, --accept-identity N   accept early alignment if the normalized alignment score is >= N and -F or -G is set" << endl
         << "    -A, --max-attempts N      try to improve sensitivity and align this many times (default: 7)" << endl
         << "    -v  --map-qual-method OPT mapping quality method: 0 - none, 1 - fast approximation, 2 - exact (default 1)" << endl
         << "    -S, --sens-step N     decrease maximum MEM size or kmer size by N bp until alignment succeeds (default: 0/off)" << endl
         << "maximal exact match (MEM) mapper:" << endl
         << "  This algorithm is used when --kmer-size is not specified and a GCSA index is given" << endl
         << "    -L, --min-mem-length N   ignore MEMs shorter than this length (default: 8)" << endl
         << "    -Y, --max-mem-length N   ignore MEMs longer than this length by stopping backward search (default: 0/unset)" << endl
         << "    -V, --mem-reseed N       reseed MEMs longer than this length (default: 64)" << endl
         << "    -a, --id-clustering      use id clustering to drive the mapper, rather than MEM-threading" << endl
         << "kmer-based mapper:" << endl
         << "  This algorithm is used when --kmer-size is specified or a rocksdb index is given" << endl
         << "    -k, --kmer-size N     use this kmer size, it must be < kmer size in db (default: from index)" << endl
         << "    -j, --kmer-stride N   step distance between succesive kmers to use for seeding (default: kmer size)" << endl
         << "    -E, --min-kmer-entropy N  require shannon entropy of this in order to use kmer (default: no limit)" << endl
         << "    -l, --kmer-min N      give up aligning if kmer size gets below this threshold (default: 8)" << endl
         << "    -F, --prefer-forward  if the forward alignment of the read works, accept it" << endl;
}

int main_map(int argc, char** argv) {

    if (argc == 2) {
        help_map(argv);
        return 1;
    }

    string seq;
    string qual;
    string seq_name;
    string db_name;
    string xg_name;
    string gcsa_name;
    int kmer_size = 0;
    int kmer_stride = 0;
    int sens_step = 0;
    int best_clusters = 0;
    int cluster_min = 1;
    int max_attempts = 7;
    string read_file;
    string hts_file;
    bool keep_secondary = false;
    int hit_max = 100;
    int max_multimaps = 1;
    int thread_count = 1;
    int thread_ex = 10;
    int context_depth = 7;
    bool output_json = false;
    bool debug = false;
    bool prefer_forward = false;
    bool greedy_accept = false;
    float min_score = 0;
    string sample_name;
    string read_group;
    string fastq1, fastq2;
    bool interleaved_input = false;
    int pair_window = 64; // ~11bp/node
    int band_width = 1000; // anything > 1000bp sequences is difficult to align efficiently
    bool try_both_mates_first = false;
    bool always_rescue = false;
    bool top_pairs_only = false;
    float min_kmer_entropy = 0;
    float accept_identity = 0;
    size_t kmer_min = 8;
    int softclip_threshold = 0;
    int max_mem_length = 0;
    int min_mem_length = 8;
    int mem_reseed_length = 64;
    bool mem_threading = true;
    int max_target_factor = 100;
    int buffer_size = 100;
    int match = 1;
    int mismatch = 4;
    int gap_open = 6;
    int gap_extend = 1;
    int full_length_bonus = 5;
    bool qual_adjust_alignments = false;
    int extra_multimaps = 1;
    int max_mapping_quality = 60;
    int method_code = 1;
    string gam_input;
    bool compare_gam = false;
    int fragment_max = 1e4;
    double fragment_sigma = 10;
    

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

            {
                /* These options set a flag. */
                //{"verbose", no_argument,       &verbose_flag, 1},
                {"sequence", required_argument, 0, 's'},
                {"quality", required_argument, 0, 'I'},
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
                {"max-multimaps", required_argument, 0, 'M'},
                {"threads", required_argument, 0, 't'},
                {"prefer-forward", no_argument, 0, 'F'},
                {"gam-input", required_argument, 0, 'G'},
                {"accept-identity", required_argument, 0, 'X'},
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
                {"min-identity", required_argument, 0, 'P'},
                {"always-rescue", no_argument, 0, 'U'},
                {"top-pairs-only", no_argument, 0, 'O'},
                {"kmer-min", required_argument, 0, 'l'},
                {"softclip-trig", required_argument, 0, '3'},
                {"debug", no_argument, 0, 'D'},
                {"min-mem-length", required_argument, 0, 'L'},
                {"max-mem-length", required_argument, 0, 'Y'},
                {"mem-reseed-length", required_argument, 0, 'V'},
                {"id-clustering", no_argument, 0, 'a'},
                {"max-target-x", required_argument, 0, 'H'},
                {"buffer-size", required_argument, 0, 'Z'},
                {"match", required_argument, 0, 'q'},
                {"mismatch", required_argument, 0, 'z'},
                {"gap-open", required_argument, 0, 'o'},
                {"gap-extend", required_argument, 0, 'y'},
                {"qual-adjust", no_argument, 0, '1'},
                {"extra-multimaps", required_argument, 0, 'u'},
                {"map-qual-method", required_argument, 0, 'v'},
                {"compare", no_argument, 0, 'w'},
                {"fragment-max", required_argument, 0, 'W'},
                {"fragment-sigma", required_argument, 0, '2'},
                {"full-l-bonus", required_argument, 0, 'T'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "s:I:j:hd:x:g:c:r:m:k:M:t:DX:FS:Jb:KR:N:if:p:B:h:G:C:A:E:Q:n:P:UOl:e:T:L:Y:H:Z:q:z:o:y:1u:v:wW:a2:3:V:",
                         long_options, &option_index);


        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 's':
            seq = optarg;
            break;

        case 'I':
            qual = string_quality_char_to_short(string(optarg));
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
            full_length_bonus = atoi(optarg);
            break;

        case '3':
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
            interleaved_input = true;
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
            gam_input = optarg;
            break;

        case 'X':
            accept_identity = atof(optarg);
            greedy_accept = true;
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

        case 'U':
            always_rescue = true;
            break;

        case 'O':
            top_pairs_only = true;
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

        case 'V':
            mem_reseed_length = atoi(optarg);
            break;

        case 'a':
            mem_threading = false;
            break;

        case 'H':
            max_target_factor = atoi(optarg);
            break;

        case 'Z':
            buffer_size = atoi(optarg);
            break;

        case 'q':
            match = atoi(optarg);
            break;

        case 'z':
            mismatch = atoi(optarg);
            break;

        case 'o':
            gap_open = atoi(optarg);
            break;

        case 'y':
            gap_extend = atoi(optarg);
            break;

        case '1':
            qual_adjust_alignments = true;
            break;

        case 'u':
            extra_multimaps = atoi(optarg);
            break;

        case 'v':
            method_code = atoi(optarg);
            break;

        case 'w':
            compare_gam = true;
            output_json = true;
            break;

        case 'W':
            fragment_max = atoi(optarg);
            break;

        case '2':
            fragment_sigma = atof(optarg);
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

    if (seq.empty() && read_file.empty() && hts_file.empty() && fastq1.empty() && gam_input.empty()) {
        cerr << "error:[vg map] a sequence or read file is required when mapping" << endl;
        return 1;
    }

    if (!qual.empty() && (seq.length() != qual.length())) {
        cerr << "error:[vg map] sequence and base quality string must be the same length" << endl;
        return 1;
    }

    if (qual_adjust_alignments && ((fastq1.empty() && hts_file.empty() && qual.empty()) // must have some quality input
                                   || (!seq.empty() && qual.empty())                    // can't provide sequence without quality
                                   || !read_file.empty()))                              // can't provide sequence list without qualities
    {
        cerr << "error:[vg map] quality adjusted alignments require base quality scores for all sequences" << endl;
        return 1;
    }
    // note: still possible that hts file types don't have quality, but have to check the file to know

    MappingQualityMethod mapping_quality_method;
    if (method_code == 0) {
        mapping_quality_method = None;
    }
    else if (method_code == 1) {
        mapping_quality_method = Approx;
    }
    else if (method_code == 2) {
        mapping_quality_method = Exact;
    }
    else {
        cerr << "error:[vg map] unrecognized mapping quality method command line arg '" << method_code << "'" << endl;
        return 1;
    }


    string file_name;
    if (optind < argc) {
        file_name = get_input_file_name(optind, argc, argv);
    }

    if (gcsa_name.empty() && !file_name.empty()) {
        gcsa_name = file_name + gcsa::GCSA::EXTENSION;
    }

    if (xg_name.empty() && !file_name.empty()) {
        xg_name = file_name + ".xg";
    }

    if (db_name.empty() && !file_name.empty()) {
        db_name = file_name + ".index";
    }

    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);

    // Load up our indexes.
    xg::XG* xindex = nullptr;
    gcsa::GCSA* gcsa = nullptr;
    gcsa::LCPArray* lcp = nullptr;

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
    auto output_alignments = [&output_buffer, &output_json, &buffer_size](vector<Alignment>& alignments) {
        // for(auto& alignment : alignments){
        //     cerr << "This is in output_alignments" << alignment.DebugString() << endl;
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

            stream::write_buffered(cout, output_buf, buffer_size);
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
        m->accept_identity = accept_identity;
        m->kmer_sensitivity_step = sens_step;
        m->prefer_forward = prefer_forward;
        m->greedy_accept = greedy_accept;
        m->thread_extension = thread_ex;
        m->cluster_min = cluster_min;
        m->context_depth = context_depth;
        m->max_attempts = max_attempts;
        m->min_kmer_entropy = min_kmer_entropy;
        m->kmer_min = kmer_min;
        m->min_identity = min_score;
        m->softclip_threshold = softclip_threshold;
        m->min_mem_length = min_mem_length;
        m->mem_reseed_length = mem_reseed_length;
        m->mem_threading = mem_threading;
        m->max_target_factor = max_target_factor;
        m->set_alignment_scores(match, mismatch, gap_open, gap_extend);
        m->adjust_alignments_for_base_quality = qual_adjust_alignments;
        m->extra_multimaps = extra_multimaps;
        m->mapping_quality_method = mapping_quality_method;
        m->always_rescue = always_rescue;
        m->fragment_max = fragment_max;
        m->fragment_sigma = fragment_sigma;
        m->full_length_alignment_bonus = full_length_bonus;
        m->max_mapping_quality = max_mapping_quality;
        mapper[i] = m;
    }

    if (!seq.empty()) {
        int tid = omp_get_thread_num();

        Alignment unaligned;
        unaligned.set_sequence(seq);

        if (!qual.empty()) {
            unaligned.set_quality(qual);
        }

        vector<Alignment> alignments = mapper[tid]->align_multi(unaligned, kmer_size, kmer_stride, max_mem_length, band_width);
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

                    vector<Alignment> alignments = mapper[tid]->align_multi(unaligned, kmer_size, kmer_stride, max_mem_length, band_width);
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
             &max_mem_length,
             &band_width]
                (Alignment& alignment) {

                    if(alignment.is_secondary() && !keep_secondary) {
                        // Skip over secondary alignments in the input; we don't want several output mappings for each input *mapping*.
                        return;
                    }

                    int tid = omp_get_thread_num();
                    vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);
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
        if (interleaved_input) {
            // paired interleaved
            auto output_func = [&output_alignments,
                                &compare_gam]
                (Alignment& aln1,
                 Alignment& aln2,
                 pair<vector<Alignment>, vector<Alignment>>& alnp) {
                // Output the alignments in JSON or protobuf as appropriate.
                output_alignments(alnp.first);
                output_alignments(alnp.second);
            };
            function<void(Alignment&,Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &pair_window,
                 &top_pairs_only,
                 &output_func](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1, aln2, queued_resolve_later, kmer_size, kmer_stride, max_mem_length, band_width, pair_window, top_pairs_only);
                if (!queued_resolve_later) {
                    output_func(aln1, aln2, alnp);
                    // check if we should try to align the queued alignments
                    if (our_mapper->fragment_size != 0
                        && !our_mapper->imperfect_pairs_to_retry.empty()) {
                        int i = 0;
                        for (auto p : our_mapper->imperfect_pairs_to_retry) {
                            auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                                       queued_resolve_later, kmer_size,
                                                                       kmer_stride, max_mem_length,
                                                                       band_width, pair_window,
                                                                       top_pairs_only);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }
            };
            fastq_paired_interleaved_for_each_parallel(fastq1, lambda);
#pragma omp parallel
            { // clean up buffered alignments that weren't perfect
                auto our_mapper = mapper[omp_get_thread_num()];
                // if we haven't yet computed these, assume we couldn't get an estimate for fragment size
                our_mapper->fragment_size = fragment_max;
                for (auto p : our_mapper->imperfect_pairs_to_retry) {
                    bool queued_resolve_later = false;
                    auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                               queued_resolve_later, kmer_size,
                                                               kmer_stride, max_mem_length,
                                                               band_width, pair_window,
                                                               top_pairs_only);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        } else if (fastq2.empty()) {
            // single
            function<void(Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width]
                    (Alignment& alignment) {

                        int tid = omp_get_thread_num();
                        vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);

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
            auto output_func = [&output_alignments]
                (Alignment& aln1,
                 Alignment& aln2,
                 pair<vector<Alignment>, vector<Alignment>>& alnp) {
                // Make sure we have unaligned "alignments" for things that don't align.
                // Output the alignments in JSON or protobuf as appropriate.
                output_alignments(alnp.first);
                output_alignments(alnp.second);
            };
            function<void(Alignment&,Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &pair_window,
                 &top_pairs_only,
                 &output_func](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1, aln2, queued_resolve_later, kmer_size, kmer_stride, max_mem_length, band_width, pair_window, top_pairs_only);
                if (!queued_resolve_later) {
                    output_func(aln1, aln2, alnp);
                    // check if we should try to align the queued alignments
                    if (our_mapper->fragment_size != 0
                        && !our_mapper->imperfect_pairs_to_retry.empty()) {
                        int i = 0;
                        for (auto p : our_mapper->imperfect_pairs_to_retry) {
                            auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                                       queued_resolve_later, kmer_size,
                                                                       kmer_stride, max_mem_length,
                                                                       band_width, pair_window,
                                                                       top_pairs_only);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }
            };
            fastq_paired_two_files_for_each_parallel(fastq1, fastq2, lambda);
#pragma omp parallel
            {
                auto our_mapper = mapper[omp_get_thread_num()];
                our_mapper->fragment_size = fragment_max;
                for (auto p : our_mapper->imperfect_pairs_to_retry) {
                    bool queued_resolve_later = false;
                    auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                               queued_resolve_later, kmer_size,
                                                               kmer_stride, max_mem_length,
                                                               band_width, pair_window,
                                                               top_pairs_only);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        }
    }

    if (!gam_input.empty()) {
        ifstream gam_in(gam_input);
        if (interleaved_input) {
            auto output_func = [&output_alignments,
                                &compare_gam]
                (Alignment& aln1,
                 Alignment& aln2,
                 pair<vector<Alignment>, vector<Alignment>>& alnp) {
                if (compare_gam) {
#pragma omp critical (cout)
                    {
                        cout << aln1.name() << "\t" << overlap(aln1.path(), alnp.first.front().path())
                                            << "\t" << alnp.first.front().identity()
                                            << "\t" << alnp.first.front().score()
                                            << "\t" << alnp.first.front().mapping_quality() << endl
                             << aln2.name() << "\t" << overlap(aln2.path(), alnp.second.front().path())
                                            << "\t" << alnp.second.front().identity()
                                            << "\t" << alnp.second.front().score()
                                            << "\t" << alnp.second.front().mapping_quality() << endl;
                    }
                } else {
                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alnp.first);
                    output_alignments(alnp.second);
                }
            };
            function<void(Alignment&,Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &compare_gam,
                 &pair_window,
                 &top_pairs_only,
                 &output_func](Alignment& aln1, Alignment& aln2) {
                auto our_mapper = mapper[omp_get_thread_num()];
                bool queued_resolve_later = false;
                auto alnp = our_mapper->align_paired_multi(aln1, aln2, queued_resolve_later, kmer_size, kmer_stride, max_mem_length, band_width, pair_window, top_pairs_only);
                if (!queued_resolve_later) {
                    output_func(aln1, aln2, alnp);
                    // check if we should try to align the queued alignments
                    if (our_mapper->fragment_size != 0
                        && !our_mapper->imperfect_pairs_to_retry.empty()) {
                        int i = 0;
                        for (auto p : our_mapper->imperfect_pairs_to_retry) {
                            auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                                       queued_resolve_later, kmer_size,
                                                                       kmer_stride, max_mem_length,
                                                                       band_width, pair_window,
                                                                       top_pairs_only);
                            output_func(p.first, p.second, alnp);
                        }
                        our_mapper->imperfect_pairs_to_retry.clear();
                    }
                }
            };
            stream::for_each_interleaved_pair_parallel(gam_in, lambda);
#pragma omp parallel
            {
                auto our_mapper = mapper[omp_get_thread_num()];
                our_mapper->fragment_size = fragment_max;
                for (auto p : our_mapper->imperfect_pairs_to_retry) {
                    bool queued_resolve_later = false;
                    auto alnp = our_mapper->align_paired_multi(p.first, p.second,
                                                               queued_resolve_later, kmer_size,
                                                               kmer_stride, max_mem_length,
                                                               band_width, pair_window,
                                                               top_pairs_only);
                    output_func(p.first, p.second, alnp);
                }
                our_mapper->imperfect_pairs_to_retry.clear();
            }
        } else {
            function<void(Alignment&)> lambda =
                [&mapper,
                 &output_alignments,
                 &keep_secondary,
                 &kmer_size,
                 &kmer_stride,
                 &max_mem_length,
                 &band_width,
                 &compare_gam]
                (Alignment& alignment) {
                int tid = omp_get_thread_num();
                vector<Alignment> alignments = mapper[tid]->align_multi(alignment, kmer_size, kmer_stride, max_mem_length, band_width);
                if(alignments.empty()) {
                    alignments.push_back(alignment);
                }
                if (compare_gam) {
#pragma omp critical (cout)
                    cout << alignment.name() << "\t" << overlap(alignment.path(), alignments.front().path())
                                             << "\t" << alignments.front().identity()
                                             << "\t" << alignments.front().score()
                                             << "\t" << alignments.front().mapping_quality() << endl;
                } else {
                    // Output the alignments in JSON or protobuf as appropriate.
                    output_alignments(alignments);
                }
            };
            stream::for_each_parallel(gam_in, lambda);
        }
        gam_in.close();
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
         << "    -Z, --translation-in input is a graph translation description" << endl

         << "    -t, --turtle         output RDF/turtle format (can not be loaded by VG)" << endl
         << "    -T, --turtle-in      input turtle format." << endl
         << "    -r, --rdf_base_uri   set base uri for the RDF output" << endl

         << "    -a, --align-in       input GAM format" << endl
         << "    -A, --aln-graph GAM  add alignments from GAM to the graph" << endl

         << "    -q, --locus-in       input stream is Locus format" << endl
         << "    -z, --locus-out      output stream Locus format" << endl
         << "    -Q, --loci FILE      input is Locus format for use by dot output" << endl

         << "    -d, --dot            output dot format" << endl
         << "    -S, --simple-dot     simplify the dot output; remove node labels, simplify alignments" << endl
         << "    -B, --bubble-label   label nodes with emoji/colors that correspond to superbubbles" << endl
         << "    -Y, --ultra-label    same as -Y but using ultrabubbles" << endl
         << "    -m, --skip-missing   skip mappings to nodes not in the graph when drawing alignments" << endl
         << "    -C, --color          color nodes that are not in the reference path (DOT OUTPUT ONLY)" << endl
         << "    -p, --show-paths     show paths in dot output" << endl
         << "    -w, --walk-paths     add labeled edges to represent paths in dot output" << endl
         << "    -n, --annotate-paths add labels to normal edges to represent paths in dot output" << endl
         << "    -M, --show-mappings  with -p print the mappings in each path in JSON" << endl
         << "    -I, --invert-ports   invert the edge ports in dot so that ne->nw is reversed" << endl
         << "    -s, --random-seed N  use this seed when assigning path symbols in dot output" << endl

         << "    -b, --bam            input BAM or other htslib-parseable alignments" << endl

         << "    -f, --fastq-in       input fastq (output defaults to GAM). Takes two " << endl
         << "                         positional file arguments if paired" << endl
         << "    -X, --fastq-out      output fastq (input defaults to GAM)" << endl
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
    string loci_file;
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
    bool superbubble_ranking = false;
    bool superbubble_labeling = false;
    bool ultrabubble_labeling = false;
    bool skip_missing_nodes = false;

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
                {"fastq-in", no_argument, 0, 'f'},
                {"fastq-out", no_argument, 0, 'X'},
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
                {"translation-in", no_argument, 0, 'Z'},
                {"ultra-label", no_argument, 0, 'Y'},
                {"bubble-label", no_argument, 0, 'B'},
                {"skip-missing", no_argument, 0, 'm'},
                {"locus-in", no_argument, 0, 'q'},
                {"loci", no_argument, 0, 'Q'},
                {"locus-out", no_argument, 0, 'z'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgFjJhvVpaGbifA:s:wnlLIMcTtr:SCZBYmqQ:zX",
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

        case 'Y':
            ultrabubble_labeling = true;
            break;

        case 'B':
            superbubble_labeling = true;
            break;

        case 'm':
            skip_missing_nodes = true;
            break;

        case 'Z':
            input_type = "translation";
            output_type = "json";
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

        case 'X':
            output_type = "fastq";
            if(input_type.empty()) {
                // Default to FASTQ -> GAM
                input_type = "gam";
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

        case 'q':
            input_type = "locus";
            if (output_type.empty()) {
                // Default to Locus -> JSON
                output_type = "json";
            }
            break;

        case 'z':
            output_type = "locus";
            break;

        case 'Q':
            loci_file = optarg;
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
    vector<Locus> loci;
    if (!loci_file.empty()) {
        function<void(Locus&)> lambda = [&loci](Locus& locus) { loci.push_back(locus); };
        ifstream in;
        in.open(loci_file.c_str());
        stream::for_each(in, lambda);
    }

    VG* graph = nullptr;
    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }
    string file_name = get_input_file_name(optind, argc, argv);
    if (input_type == "vg") {
        if (output_type == "stream") {
            function<void(Graph&)> lambda = [&](Graph& g) { cout << pb2json(g) << endl; };
            get_input_file(file_name, [&](istream& in) {
                stream::for_each(in, lambda);
            });
            return 0;
        } else {
            get_input_file(file_name, [&](istream& in) {
                graph = new VG(in);
            });
        }
        // VG can convert to any of the graph formats, so keep going
    } else if (input_type == "gfa") {
        get_input_file(file_name, [&](istream& in) {
            graph = new VG;
            graph->from_gfa(in);
        });
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
                    if(std::isnan(a.identity())) {
                        // Fix up NAN identities that can't be serialized in
                        // JSON. We shouldn't generate these any more, and they
                        // are out of spec, but they can be in files.
                        a.set_identity(0);
                    }
                    cout << pb2json(a) << "\n";
                };
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each(in, lambda);
                });
            } else if (output_type == "fastq") {
                function<void(Alignment&)> lambda = [](Alignment& a) {
                    cout << "@" << a.name() << endl
                         << a.sequence() << endl
                    << "+" << endl;
                    if (a.quality().empty()) {
                        cout << string(a.sequence().size(), quality_short_to_char(30)) << endl;
                    } else {
                        cout << string_quality_short_to_char(a.quality()) << endl;
                    }
                };
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each(in, lambda);
                });
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
        // The first FASTQ is the filename we already grabbed
        fastq1 = file_name;
        if (optind < argc) {
            // There may be a second one
            fastq2 = get_input_file_name(optind, argc, argv);
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
                function<void(Pileup&)> lambda = [](Pileup& p) {
                    cout << pb2json(p) << "\n";
                };
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each(in, lambda);
                });
            } else {
                // todo
                cerr << "[vg view] error: (binary) Pileup can only be converted to JSON" << endl;
                return 1;
            }
        } else {
            if (output_type == "json" || output_type == "pileup") {
                JSONStreamHelper<Pileup> json_helper(file_name);
                json_helper.write(cout, output_type == "json");
            } else {
                cerr << "[vg view] error: JSON Pileup can only be converted to Pileup or JSON" << endl;
                return 1;
            }
        }
        cout.flush();
        return 0;
    } else if (input_type == "translation") {
        if (output_type == "json") {
            function<void(Translation&)> lambda = [](Translation& t) {
                cout << pb2json(t) << "\n";
            };
            get_input_file(file_name, [&](istream& in) {
                stream::for_each(in, lambda);
            });
        } else {
            cerr << "[vg view] error: (binary) Translation can only be converted to JSON" << endl;
            return 1;
        }
        return 0;
    } else if (input_type == "locus") {
        if (input_json == false) {
            if (output_type == "json") {
                // convert values to printable ones
                function<void(Locus&)> lambda = [](Locus& l) {
                    cout << pb2json(l) << "\n";
                };
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each(in, lambda);
                });
            } else {
                // todo
                cerr << "[vg view] error: (binary) Locus can only be converted to JSON" << endl;
                return 1;
            }
        } else {
            if (output_type == "json" || output_type == "locus") {
                JSONStreamHelper<Locus> json_helper(file_name);
                json_helper.write(cout, output_type == "json");
            } else {
                cerr << "[vg view] error: JSON Locus can only be converted to Locus or JSON" << endl;
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
                      loci,
                      show_paths_in_dot,
                      walk_paths_in_dot,
                      annotate_paths_in_dot,
                      show_mappings_in_dot,
                      simple_dot,
                      invert_edge_ports_in_dot,
                      color_variants,
                      superbubble_ranking,
                      superbubble_labeling,
                      ultrabubble_labeling,
                      skip_missing_nodes,
                      seed_val);
    } else if (output_type == "json") {
        cout << pb2json(graph->graph) << endl;
    } else if (output_type == "gfa") {
        graph->to_gfa(std::cout);
    } else if (output_type == "turtle") {
        graph->to_turtle(std::cout, rdf_base_uri, color_variants);
    } else if (output_type == "vg") {
        graph->serialize_to_ostream(cout);
    } else if (output_type == "locus") {

    } else {
        // We somehow got here with a bad output format.
        cerr << "[vg view] error: cannot save a graph in " << output_type << " format" << endl;
        return 1;
    }

    // We made it to the end and nothing broke.
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
            map<int, Edit> edits;
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
                string allele_name = convert(name_int);
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
                        allele.set_name(convert(name));
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
    cerr << "usage: " << argv[0] << " deconstruct [options] <my_graph>.vg" << endl
         << "options: " << endl
         << " -x --xg-name  <XG>.xg an XG index from which to extract distance information." << endl
         << " -s --superbubbles  Print the superbubbles of the graph and exit." << endl
         << " -o --output <FILE>      Save output to <FILE> rather than STDOUT." << endl
         << " -d --dagify             DAGify the graph before enumeratign superbubbles" << endl
         << " -u --unroll <STEPS>    Unroll the graph <STEPS> steps before calling variation." << endl
         << " -c --compact <ROUNDS>   Perform <ROUNDS> rounds of superbubble compaction on the graph." << endl
         << " -m --mask <vcf>.vcf    Look for variants not in <vcf> in the graph" << endl
         << " -i --invert           Invert the mask (i.e. find only variants present in <vcf>.vcf. Requires -m. " << endl
         << endl;
}

int main_deconstruct(int argc, char** argv){
    //cerr << "WARNING: EXPERIMENTAL" << endl;
    if (argc <= 2) {
        help_deconstruct(argv);
        return 1;
    }

    bool print_sbs = false;
    string outfile = "";
    bool dagify = false;
    int unroll_steps = 0;
    int compact_steps = 0;
    bool invert = false;
    string mask_file = "";
    string xg_name;
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"xg-name", required_argument,0, 'x'},
                {"output", required_argument, 0, 'o'},
                {"unroll", required_argument, 0, 'u'},
                {"compact", required_argument, 0, 'c'},
                {"mask", required_argument, 0, 'm'},
                {"dagify", no_argument, 0, 'd'},
                {"superbubbles", no_argument, 0, 's'},
                {"invert", no_argument, 0, 'v'},
                {0, 0, 0, 0}

            };

            int option_index = 0;
            c = getopt_long (argc, argv, "dho:u:c:vm:sx:",
                    long_options, &option_index);

            // Detect the end of the options.
            if (c == -1)
                break;

            switch (c)
            {
                case 's':
                    print_sbs = true;
                    break;
                case 'x':
                    xg_name = optarg;
                    break;
                case 'o':
                    outfile = optarg;
                    break;
                case 'u':
                    unroll_steps = atoi(optarg);
                    break;
                case 'c':
                    compact_steps = atoi(optarg);
                    break;
                case 'm':
                    mask_file = optarg;
                    break;
                case 'd':
                    dagify = true;
                    break;
                case 'v':
                    invert = true;
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

    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in);
    });

    Deconstructor decon = Deconstructor(graph);
    if (!xg_name.empty()){
        ifstream xg_stream(xg_name);
        if(!xg_stream) {
            cerr << "Unable to open xg index: " << xg_name << endl;
            exit(1);
        }

        xg::XG* xindex = new  xg::XG(xg_stream);
        decon.set_xg(xindex);
    }

		if (unroll_steps > 0){
			cerr << "Unrolling " << unroll_steps << " steps..." << endl;
            decon.unroll_my_vg(unroll_steps);
			cerr << "Done." << endl;
		}

        if (dagify){
            int dagify_steps = 3;
            cerr << "DAGifying..." << endl;
            decon.dagify_my_vg(dagify_steps);
            cerr << "Done." << endl;
        }



    // At this point, we can detect the superbubbles

    map<pair<vg::id_t, vg::id_t>, vector<vg::id_t> > sbs = decon.get_all_superbubbles();


    if (compact_steps > 0){
        cerr << "Compacting superbubbles of graph " << compact_steps << " steps..." << endl;
        decon.compact(compact_steps);
        cerr << "Done." << endl;
    }
    if (print_sbs){
        for (auto s: sbs){
            cout << s.first.first << "\t";
            cout << "\t" << s.first.second << endl;
        }
    }
    else{
        if (xg_name.empty()){
            cerr << "An xg index must be provided for full deconstruction." << endl;
            exit(1);
        }
        decon.sb2vcf( outfile);
    }
    /* Find superbubbles */

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
         << "  -- view          format conversions for graphs and alignments" << endl
         << "  -- vectorize     transform alignments to simple ML-compatible vectors" << endl
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
         << "  -- homogenize    homogenize long variants in the graph to improve genotyping" << endl
         << "  -- surject       map alignments onto specific paths" << endl
         << "  -- msga          multiple sequence graph alignment" << endl
         << "  -- pileup        build a pileup from a set of alignments" << endl
         << "  -- call          prune the graph by genotyping a pileup" << endl
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
    } else if (command == "view") {
        return main_view(argc, argv);
    } else if (command == "align") {
        return main_align(argc, argv);
    } else if (command == "map") {
        return main_map(argc, argv);
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
    } else if (command == "surject") {
        return main_surject(argc, argv);
    } else if (command == "msga") {
        return main_msga(argc, argv);
    } else if (command == "pileup") {
        return main_pileup(argc, argv);
    } else if (command == "call") {
        return main_call(argc, argv);
    } else if (command == "genotype") {
        return main_genotype(argc, argv);
    } else if (command == "compare") {
        return main_compare(argc, argv);
    } else if (command == "validate") {
        return main_validate(argc, argv);
    } else if (command == "filter") {
        return main_filter(argc, argv);
    } else if (command == "vectorize") {
        return main_vectorize(argc, argv);
    } else if (command == "circularize"){
        return main_circularize(argc, argv);
    }  else if (command == "translate") {
        return main_translate(argc, argv);
    }  else if (command == "version") {
        return main_version(argc, argv);
    }  else if (command == "homogenize"){
        return main_homogenize(argc, argv);
    } else if (command == "sift"){
        return main_sift(argc, argv);  
    } else if (command == "test") {
        return main_test(argc, argv);
    } else if (command == "srpe"){
        return main_srpe(argc, argv);
    } else if (command == "locify"){
        return main_locify(argc, argv);
    } else if (command == "sort") {
        return main_sort(argc, argv);
    }else {
        cerr << "error:[vg] command " << command << " not found" << endl;
        vg_help(argv);
        return 1;
    }

    return 0;

}
