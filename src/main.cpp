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
         << "    -m, --min-end-matches N filter reads that don't begin with at least N matches on each end" << endl
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
                {"min-end-matches", required_argument, 0, 'm'},
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
        c = getopt_long (argc, argv, "s:r:d:e:fauo:m:Sx:R:B:Ac:vq:E:D:C:t:",
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
        case 'm':
            filter.min_end_matches = atoi(optarg);
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
    
    // Configure its temp directory to the system temp directory
    gcsa::TempFile::setDirectory(find_temp_dir());

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

    Mapper* mapper = nullptr;
    if (mem_sketch) {
        if (gcsa_name.empty()) {
            cerr << "[vg vectorize] error : an xg index and gcsa index are required when making MEM sketches" << endl;
            return 1;
        } else {
            mapper = new Mapper(xg_index, &gcsa_index, &lcp_index);
        }
        if (mem_hit_max) {
            mapper->hit_max = mem_hit_max;
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
            auto mems = mapper->find_mems_simple(a.sequence().begin(), a.sequence().end(),
                                                 max_mem_length, mapper->min_mem_length);
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


    delete mapper;

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
         << "    -N, --allow-Ns        allow reads to be sampled from the graph with Ns in them" << endl
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
    bool reads_may_contain_Ns = false;
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
            {"allow-Ns", no_argument, 0, 'N'},
            {"base-error", required_argument, 0, 'e'},
            {"indel-error", required_argument, 0, 'i'},
            {"frag-len", required_argument, 0, 'p'},
            {"frag-std-dev", required_argument, 0, 'v'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hl:n:s:e:i:fax:Jp:v:N",
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

        case 'N':
            reads_may_contain_Ns = true;
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

    Sampler sampler(xgidx, seed_val, forward_only, reads_may_contain_Ns);
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
         << "    -A, --is-acyclic      print if the graph is acyclic or not" << endl
         << "    -n, --node ID         consider node with the given id" << endl
         << "    -d, --to-head         show distance to head for each provided node" << endl
         << "    -t, --to-tail         show distance to head for each provided node" << endl
         << "    -a, --alignments FILE compute stats for reads aligned to the graph" << endl
         << "    -r, --node-id-range   X:Y where X and Y are the smallest and largest "
        "node id in the graph, respectively" << endl
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
    bool verbose = false;
    bool is_acyclic = false;
    bool stats_range = false;
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
            {"alignments", required_argument, 0, 'a'},
            {"is-acyclic", no_argument, 0, 'A'},
            {"node-id-range", no_argument, 0, 'r'},
            {"verbose", no_argument, 0, 'v'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hzlsHTScdtn:NEa:vAr",
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

        case 'A':
            is_acyclic = true;
            break;

        case 'a':
            alignments_filename = optarg;
            break;

        case 'r':
            stats_range = true;
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

    if (stats_range) {
        cout << "node-id-range\t" << graph->min_node_id() << ":" << graph->max_node_id() << endl;
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
        Aligner aligner = Aligner(match, mismatch, gap_open, gap_extend, full_length_bonus);
        alignment = graph->align(seq, &aligner, 0, pinned_alignment, pin_left,
            banded_global, 0, max(seq.size(), graph->length()), debug);
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
         << "    -l, --pileup-in      input VG Pileup format" << endl

         << "    -R, --snarl-in       input VG Snarl format" << endl
         << "    -E, --snarl-traversal-in input VG SnarlTraversal format" << endl;
    
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
                {"snarls", no_argument, 0, 'R'},
                {"snarltraversals", no_argument, 0, 'E'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgFjJhvVpaGbifA:s:wnlLIMcTtr:SCZBYmqQ:zXRE",
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

        case 'R':
            input_type = "snarls";
            if (output_type.empty()) {
                // Default to Snarl -> JSON
                output_type = "json";
            }
            break;

        case 'E':
            input_type = "snarltraversals";
            if (output_type.empty()) {
                // Default to Locus -> JSON
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
    } else if (input_type == "snarls") {
        if (output_type == "json") {
            function<void(Snarl&)> lambda = [](Snarl& s) {
                cout << pb2json(s) << "\n";
            };
            get_input_file(file_name, [&](istream& in) {
                stream::for_each(in, lambda);
            });
        } else {
            cerr << "[vg view] error: (binary) Snarls can only be converted to JSON" << endl;
            return 1;
        }
        return 0;
    } else if (input_type == "snarltraversals") {
        if (output_type == "json") {
            function<void(SnarlTraversal&)> lambda = [](SnarlTraversal& s) {
                cout << pb2json(s) << "\n";
            };
            get_input_file(file_name, [&](istream& in) {
                stream::for_each(in, lambda);
            });
        } else {
            cerr << "[vg view] error: (binary) SnarlTraversals can only be converted to JSON" << endl;
            return 1;
        }
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
         << "Outputs VCF records for Snarls present in a graph." << endl
         << "options: " << endl
         << "--path / -p     A reference path to deconstruct against." << endl
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
         << "  -- view          format conversions for graphs and alignments" << endl
         << "  -- vectorize     transform alignments to simple ML-compatible vectors" << endl
         << "  -- paths         traverse paths in the graph" << endl
         << "  -- align         local alignment" << endl
         << "  -- stats         metrics describing graph properties" << endl
         << "  -- join          combine graphs via a new head" << endl
         << "  -- ids           manipulate node ids" << endl
         << "  -- concat        concatenate graphs tail-to-head" << endl
         << "  -- kmers         enumerate kmers of the graph" << endl
         << "  -- sim           simulate reads from the graph" << endl
         << "  -- mod           filter, transform, and edit the graph" << endl
         << "  -- homogenize    homogenize long variants in the graph to improve genotyping" << endl
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
    } else if (command == "view") {
        return main_view(argc, argv);
    } else if (command == "align") {
        return main_align(argc, argv);
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
    } else if (command == "pileup") {
        return main_pileup(argc, argv);
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
