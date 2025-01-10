/** \file view_main.cpp
 *
 * Defines the "vg view" subcommand, which converts formats.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../multipath_alignment.hpp"
#include "../vg.hpp"
#include "../snarl_distance_index.hpp"
#include "../gfa.hpp"
#include "../io/json_stream_helper.hpp"
#include "../handle.hpp"
#include "../algorithms/gfa_to_handle.hpp"

#include <vg/io/message_iterator.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/hash_graph.hpp>


using namespace std;
using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

void help_view(char** argv) {
    cerr << "usage: " << argv[0] << " view [options] [ <graph.vg> | <graph.json> | <aln.gam> | <read1.fq> [<read2.fq>] ]" << endl
         << "options:" << endl
         << "    -g, --gfa                  output GFA format (default)" << endl
         << "    -F, --gfa-in               input GFA format, reducing overlaps if they occur" << endl

         << "    -v, --vg                   output VG format [DEPRECATED, use vg convert instead]" << endl
         << "    -V, --vg-in                input VG format only" << endl

         << "    -j, --json                 output JSON format" << endl
         << "    -J, --json-in              input JSON format (use with e.g. -a as necessary) << endl
         << "    -c, --json-stream          streaming conversion of a VG format graph in line delimited JSON format" << endl
         << "                               (this cannot be loaded directly via -J)" << endl

         << "    -G, --gam                  output GAM format (vg alignment format: Graph Alignment/Map)" << endl
         << "    -Z, --translation-in       input is a graph translation description" << endl

         << "    -t, --turtle               output RDF/turtle format (can not be loaded by VG)" << endl
         << "    -T, --turtle-in            input turtle format." << endl
         << "    -r, --rdf_base_uri         set base uri for the RDF output" << endl

         << "    -a, --align-in             input GAM format, or JSON version of GAM format" << endl
         << "    -A, --aln-graph GAM        add alignments from GAM to the graph" << endl

         << "    -q, --locus-in             input is Locus format, or JSON version of Locus format" << endl
         << "    -z, --locus-out            output is Locus format" << endl
         << "    -Q, --loci FILE            input is Locus format for use by dot output" << endl

         << "    -d, --dot                  output dot format" << endl
         << "    -S, --simple-dot           simplify the dot output; remove node labels, simplify alignments" << endl
         << "    -u, --noseq-dot            shows size information instead of sequence in the dot output" << endl
         << "    -e, --ascii-labels         use labels for paths or superbubbles with char/colors rather than emoji" << endl
         << "    -Y, --ultra-label          label nodes with emoji/colors that correspond to ultrabubbles" << endl
         << "    -m, --skip-missing         skip mappings to nodes not in the graph when drawing alignments" << endl
         << "    -C, --color                color nodes that are not in the reference path (DOT OUTPUT ONLY)" << endl
         << "    -p, --show-paths           show paths in dot output" << endl
         << "    -w, --walk-paths           add labeled edges to represent paths in dot output" << endl
         << "    -n, --annotate-paths       add labels to normal edges to represent paths in dot output" << endl
         << "    -M, --show-mappings        with -p print the mappings in each path in JSON" << endl
         << "    -I, --invert-ports         invert the edge ports in dot so that ne->nw is reversed" << endl
         << "    -s, --random-seed N        use this seed when assigning path symbols in dot output" << endl

         << "    -b, --bam                  input BAM or other htslib-parseable alignments" << endl

         << "    -f, --fastq-in             input fastq (output defaults to GAM). Takes two " << endl
         << "                               positional file arguments if paired" << endl
         << "    -X, --fastq-out            output fastq (input defaults to GAM)" << endl
         << "    -i, --interleaved          fastq is interleaved paired-ended" << endl

         << "    -L, --pileup               output VG Pileup format" << endl
         << "    -l, --pileup-in            input VG Pileup format, or JSON version of VG Pileup format" << endl

         << "    -B, --distance-in          input distance index" << endl
         << "    -R, --snarl-in             input VG Snarl format" << endl
         << "    -E, --snarl-traversal-in   input VG SnarlTraversal format" << endl
         << "    -K, --multipath-in         input VG MultipathAlignment format (GAMP), or JSON version of GAMP format" << endl
         << "    -k, --multipath            output VG MultipathAlignment format (GAMP)" << endl
         << "    -D, --expect-duplicates    don't warn if encountering the same node or edge multiple times" << endl
         << "    -x, --extract-tag TAG      extract and concatenate messages with the given tag" << endl
         << "    --verbose                  explain the file being read with --extract-tag" << endl
         << "    --threads N                for parallel operations use this many threads [1]" << endl;
    
    // TODO: Can we regularize the option names for input and output types?

}

int main_view(int argc, char** argv) {

    if (argc == 2) {
        help_view(argv);
        return 1;
    }

    // Supported conversions:
    //      TO  hg  vg  json    gfa gam bam fastq   dot
    // FROM
    // hg       N   Y   Y       Y   N   N   N       Y
    // vg       N   Y   Y       Y   N   N   N       Y
    // json     N   Y   Y       Y   N   N   N       Y
    // gfa      N   Y   Y       Y   N   N   N       Y
    // gam      N   N   Y       N   N   N   N       N
    // bam      N   N   N       N   Y   N   N       N
    // fastq    N   N   N       N   Y   N   N       N
    // dot      N   N   N       N   N   N   N       N
    //
    // and json-gam -> gam
    //     json-pileup -> pileup

    bool which_read_in_pair = true;
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
    bool noseq_dot = false;
    int seed_val = time(NULL);
    bool color_variants = false;
    bool ultrabubble_labeling = false;
    bool skip_missing_nodes = false;
    bool expect_duplicates = false;
    string extract_tag;
    bool verbose;
    bool ascii_labels = false;
    omp_set_num_threads(1); // default to 1 thread
    
    #define OPT_VERBOSE 1000

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
                {"rdf-base-uri", required_argument, 0, 'r'},
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
                {"noseq-dot", no_argument, 0, 'u'},
                {"color", no_argument, 0, 'C'},
                {"translation-in", no_argument, 0, 'Z'},
                {"ultra-label", no_argument, 0, 'Y'},
                {"skip-missing", no_argument, 0, 'm'},
                {"locus-in", no_argument, 0, 'q'},
                {"loci", no_argument, 0, 'Q'},
                {"locus-out", no_argument, 0, 'z'},
                {"distance-in", no_argument, 0, 'B'},
                {"snarl-in", no_argument, 0, 'R'},
                {"snarl-traversal-in", no_argument, 0, 'E'},
                {"expect-duplicates", no_argument, 0, 'D'},
                {"extract-tag", required_argument, 0, 'x'},
                {"verbose", no_argument, 0, OPT_VERBOSE},
                {"multipath", no_argument, 0, 'k'},
                {"multipath-in", no_argument, 0, 'K'},
                {"ascii-labels", no_argument, 0, 'e'},
                {"threads", required_argument, 0, '7'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "dgFjJhvVpaGbifA:s:wnlLIMcTtr:SuCZYmqQ:zXBREDx:kKe7:",
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

        case 'u':
            noseq_dot = true;
            break;

        case 'Y':
            ultrabubble_labeling = true;
            break;

        case 'e':
            ascii_labels = true;
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
            seed_val = parse<int>(optarg);
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

        case 'B':
            input_type = "distance";
            if (output_type.empty()) {
                // Default to DistanceIndex -> JSON
                output_type = "json";
            }
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
            
        case 'K':
            input_type = "multipath";
            break;
            
        case 'k':
            output_type = "multipath";
            break;
            
        case 'D':
            expect_duplicates = true;
            break;
            
        case 'x':
            extract_tag = optarg;
            break;

        case OPT_VERBOSE:
            verbose = true;
            break;

        case '7':
            omp_set_num_threads(parse<int>(optarg));
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

    // If the user specified nothing else, we default to handle graph in and GFA out.
    if (input_type.empty()) {
        input_type = "handlegraph";
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
        vg::io::for_each(in, lambda);
    }
    vector<Locus> loci;
    if (!loci_file.empty()) {
        function<void(Locus&)> lambda = [&loci](Locus& locus) { loci.push_back(locus); };
        ifstream in;
        in.open(loci_file.c_str());
        vg::io::for_each(in, lambda);
    }
    
    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }
    if (output_type == "vg") {
        cerr << "[vg view] warning: vg-protobuf output (-v / --v) is deprecated. please use vg convert instead." << endl;
    }
    
    string file_name = get_input_file_name(optind, argc, argv);
    
    // Tag extraction has to be handled specially.
    // TODO: We can't dump just untagged messages.
    if (!extract_tag.empty()) {
        get_input_file(file_name, [&](istream& in) {
            // Iterate over the input as tagged messages.
            vg::io::MessageIterator it(in, verbose);
            while(it.has_current()) {
                if ((*it).first == extract_tag && (*it).second.get() != nullptr) {
                    // We match the tag, so dump this message.
                    if (verbose) {
                        cerr << "Message of " << (*it).second->size() << " bytes in matches tag to extract" << endl;
                    }
                    cout << *((*it).second.get());
                } else {
                    if (verbose) {
                        cerr << "Message of " << (*it).second->size() << " bytes does not match tag; skip" << endl;
                    }
                }
                ++it;
            }
            if (verbose) {
                cerr << "Iterator no longer has messages" << endl;
            }
        });
        return 0;
    }
    

    unique_ptr<PathHandleGraph> graph;
     
    if (input_type == "vg") {
        if (output_type == "stream") {
            function<void(Graph&)> lambda = [&](Graph& g) { cout << pb2json(g) << endl; };
            get_input_file(file_name, [&](istream& in) {
                vg::io::for_each(in, lambda);
            });
            return 0;
        } else {
            get_input_file(file_name, [&](istream& in) {
                graph = make_unique<VG>(in, false, !expect_duplicates);
            });
        }
        // VG can convert to any of the graph formats, so keep going
    } else if (input_type == "handlegraph") {
        if (output_type == "stream") {
            cerr << "[vg view] error: Cannot stream a generic HandleGraph to JSON" << endl;
            exit(1);
        } else {
            graph = vg::io::VPKG::load_one<PathHandleGraph>(file_name);
        }
    } else if (input_type == "gfa") {
        graph = make_unique<bdsg::HashGraph>();
       
        try {
            // Use the disk-backed GFA loader that `vg convert` also uses.
            vg::algorithms::gfa_to_path_handle_graph(file_name,
                                                 dynamic_cast<MutablePathMutableHandleGraph*>(graph.get()),
                                                 nullptr,
                                                 0); // set rgfa path rank to 0 to be consistent with vg convert's default logic
        } catch (vg::algorithms::GFAFormatError& e) {
            cerr << "error:[vg view] Input GFA is not acceptable." << endl;
            cerr << e.what() << endl;
            exit(1);
        } catch (std::ios_base::failure& e) {
            cerr << "error:[vg view] IO error processing input GFA." << endl;
            cerr << e.what() << endl;
            exit(1);
        }
        
        // GFA can convert to any of the graph formats, so keep going
    } else if(input_type == "json") {
        assert(input_json);
        vg::io::JSONStreamHelper<Graph> json_helper(file_name);
        function<bool(Graph&)> get_next_graph = json_helper.get_read_fn();
        // TODO: This is less inversion of control and more putting control in the middle.
        graph = make_unique<VG>([&](const function<void(Graph&)> use_graph) {
            Graph g;
            while (get_next_graph(g)) {
                use_graph(g);
            }
        }, false, !expect_duplicates);
    } else if(input_type == "turtle-in") {
        // TODO: Only vg::VG can load Turtle right now
        graph = make_unique<VG>();
        bool pre_compress=color_variants;
        if (file_name == "-") {
            dynamic_cast<vg::VG*>(graph.get())->from_turtle("/dev/stdin", rdf_base_uri);
        } else {
            dynamic_cast<vg::VG*>(graph.get())->from_turtle(file_name, rdf_base_uri);
        }
    } else if (input_type == "gam") {
        if (!input_json) {
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
                    vg::io::for_each(in, lambda);
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
                    vg::io::for_each(in, lambda);
                });
            }
            else if (output_type == "multipath") {
                vector<MultipathAlignment> buf;
                function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
                    multipath_alignment_t mp_aln;
                    to_multipath_alignment(aln, mp_aln);
                    buf.emplace_back();
                    to_proto_multipath_alignment(mp_aln, buf.back());
                    transfer_proto_metadata(aln, buf.back());
                    vg::io::write_buffered(cout, buf, 1000);
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
                vg::io::write_buffered(cout, buf, 0);
            }
            else {
                // todo
                cerr << "[vg view] error: (binary) GAM can only be converted to JSON, GAMP or FASTQ" << endl;
                return 1;
            }
        } else {
            vg::io::JSONStreamHelper<Alignment> json_helper(file_name);
            if (output_type == "json" || output_type == "gam") {
                json_helper.write(cout, output_type == "json");
            }
            else if (output_type == "multipath") {
                vector<MultipathAlignment> buf;
                Alignment aln;
                while (json_helper.get_read_fn()(aln)) {
                    multipath_alignment_t mp_aln;
                    to_multipath_alignment(aln, mp_aln);
                    buf.emplace_back();
                    to_proto_multipath_alignment(mp_aln, buf.back());
                    vg::io::write_buffered(std::cout, buf, 1000);
                }
                vg::io::write_buffered(cout, buf, 0);
            }
            else {
                cerr << "[vg view] error: JSON GAM can only be converted to GAM, GAMP, or JSON" << endl;
                return 1;
            }
        }
        cout.flush();
        return 0;
    } else if (input_type == "bam") {
        if (output_type == "gam") {
            vg::io::ProtobufEmitter<Alignment> buf(std::cout);
            function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
                buf.write(std::move(aln));
            };
            hts_for_each(file_name, lambda);
            return 0;
        } else if (output_type == "json") {
            // todo
            cerr << "[vg view] error: BAM to JSON conversion not yet implemented" << endl;
            return 0;
        } else {
            cerr << "[vg view] error: BAM can only be converted to GAM" << endl;
            return 1;
        }
    } else if (input_type == "multipath") {
        if (input_json) {
            vg::io::JSONStreamHelper<MultipathAlignment> json_helper(file_name);
            if (output_type == "multipath") {
                json_helper.write(cout, false);
            }
            else if (output_type == "fastq") {
                MultipathAlignment mp_aln;
                while (json_helper.get_read_fn()(mp_aln)) {
                    cout << "@" << mp_aln.name() << endl
                         << mp_aln.sequence() << endl
                         << "+" << endl;
                    if (mp_aln.quality().empty()) {
                        cout << string(mp_aln.sequence().size(), quality_short_to_char(30)) << endl;
                    } else {
                        cout << string_quality_short_to_char(mp_aln.quality()) << endl;
                    }
                }
                
            }
            else if (output_type == "gam") {
                vector<Alignment> buf;
                MultipathAlignment proto_mp_aln;
                while (json_helper.get_read_fn()(proto_mp_aln)) {
                    multipath_alignment_t mp_aln;
                    from_proto_multipath_alignment(proto_mp_aln, mp_aln);
                    buf.emplace_back();
                    optimal_alignment(mp_aln, buf.back());
                    transfer_proto_metadata(proto_mp_aln, buf.back());
                    if (!proto_mp_aln.paired_read_name().empty()) {
                        // alternate using next/prev
                        if (which_read_in_pair) {
                            buf.back().mutable_fragment_next()->set_name(proto_mp_aln.paired_read_name());
                        }
                        else {
                            buf.back().mutable_fragment_prev()->set_name(proto_mp_aln.paired_read_name());
                        }
                    }
                    which_read_in_pair = !which_read_in_pair;
                    vg::io::write_buffered(std::cout, buf, 1000);
                }
                vg::io::write_buffered(cout, buf, 0);
            }
            else if (output_type == "json") {
                json_helper.write(cout, true);
            }
            else if (output_type == "dot") {
                MultipathAlignment proto_mp_aln;
                while (json_helper.get_read_fn()(proto_mp_aln)) {
                    multipath_alignment_t mp_aln;
                    from_proto_multipath_alignment(proto_mp_aln, mp_aln);
                    view_multipath_alignment_as_dot(std::cout, mp_aln, true);
                }
            }
            else {
                cerr << "[vg view] error: Unrecognized output format for MultipathAlignment (GAMP)" << endl;
                return 1;
            }
            return 0;
        }
        else {
            if (output_type == "multipath") {
                vector<MultipathAlignment> buf;
                function<void(MultipathAlignment&)> lambda = [&buf](MultipathAlignment& mp_aln) {
                    buf.push_back(mp_aln);
                    vg::io::write_buffered(cout, buf, 1000);
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
                vg::io::write_buffered(std::cout, buf, 0);
            }
            else if (output_type == "fastq") {
                function<void(MultipathAlignment&)> lambda = [](MultipathAlignment& mp_aln) {
                    cout << "@" << mp_aln.name() << endl
                         << mp_aln.sequence() << endl
                         << "+" << endl;
                    if (mp_aln.quality().empty()) {
                        cout << string(mp_aln.sequence().size(), quality_short_to_char(30)) << endl;
                    } else {
                        cout << string_quality_short_to_char(mp_aln.quality()) << endl;
                    }
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
            else if (output_type == "gam") {
                vector<Alignment> buf;
                function<void(MultipathAlignment&)> lambda = [&buf](MultipathAlignment& proto_mp_aln) {
                    multipath_alignment_t mp_aln;
                    from_proto_multipath_alignment(proto_mp_aln, mp_aln);
                    buf.emplace_back();
                    optimal_alignment(mp_aln, buf.back());
                    transfer_proto_metadata(proto_mp_aln, buf.back());
                    vg::io::write_buffered(cout, buf, 1000);
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
                vg::io::write_buffered(std::cout, buf, 0);
            }
            else if (output_type == "json") {
                function<void(MultipathAlignment&)> lambda = [&](MultipathAlignment& mp_aln) {
                    cout << pb2json(mp_aln) << '\n';
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
            else if (output_type == "dot") {
                function<void(MultipathAlignment&)> lambda = [&](MultipathAlignment& proto_mp_aln) {
                    multipath_alignment_t mp_aln;
                    from_proto_multipath_alignment(proto_mp_aln, mp_aln);
                    view_multipath_alignment_as_dot(std::cout, mp_aln, true);
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
            else {
                cerr << "[vg view] error: Unrecognized output format for MultipathAlignment (GAMP)" << endl;
                return 1;
            }
            return 0;
        }
    } else if (input_type == "fastq") {
        // The first FASTQ is the filename we already grabbed
        fastq1 = file_name;
        if (optind < argc) {
            // There may be a second one
            fastq2 = get_input_file_name(optind, argc, argv);
        }
        if (output_type == "gam") {
            vg::io::ProtobufEmitter<Alignment> buf(std::cout);
            if (!interleaved_fastq && fastq2.empty()) {
                function<void(Alignment&)> lambda = [&buf](Alignment& aln) {
                    buf.write(std::move(aln));
                };
                fastq_unpaired_for_each(fastq1, lambda);
            } else if (interleaved_fastq && fastq2.empty()) {
                function<void(Alignment&, Alignment&)> lambda = [&buf](Alignment& aln1, Alignment& aln2) {
                    buf.write(std::move(aln1));
                    buf.write(std::move(aln2));
                };
                fastq_paired_interleaved_for_each(fastq1, lambda);
            } else if (!fastq2.empty()) {
                function<void(Alignment&, Alignment&)> lambda = [&buf](Alignment& aln1, Alignment& aln2) {
                    buf.write(std::move(aln1));
                    buf.write(std::move(aln2));
                };
                fastq_paired_two_files_for_each(fastq1, fastq2, lambda);
            }
        } else {
            // We can't convert fastq to the other graph formats
            cerr << "[vg view] error: FASTQ can only be converted to GAM" << endl;
            return 1;
        }
        cout.flush();
        return 0;
    } else if (input_type == "pileup") {
        if (!input_json) {
            if (output_type == "json") {
                // convert values to printable ones
                function<void(Pileup&)> lambda = [](Pileup& p) {
                    cout << pb2json(p) << "\n";
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            } else {
                // todo
                cerr << "[vg view] error: (binary) Pileup can only be converted to JSON" << endl;
                return 1;
            }
        } else {
            if (output_type == "json" || output_type == "pileup") {
                vg::io::JSONStreamHelper<Pileup> json_helper(file_name);
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
                vg::io::for_each(in, lambda);
            });
        } else {
            cerr << "[vg view] error: (binary) Translation can only be converted to JSON" << endl;
            return 1;
        }
        return 0;
    } else if (input_type == "locus") {
        if (!input_json) {
            if (output_type == "json") {
                // convert values to printable ones
                function<void(Locus&)> lambda = [](Locus& l) {
                    cout << pb2json(l) << "\n";
                };
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            } else {
                // todo
                cerr << "[vg view] error: (binary) Locus can only be converted to JSON" << endl;
                return 1;
            }
        } else {
            if (output_type == "json" || output_type == "locus") {
                vg::io::JSONStreamHelper<Locus> json_helper(file_name);
                json_helper.write(cout, output_type == "json");
            } else {
                cerr << "[vg view] error: JSON Locus can only be converted to Locus or JSON" << endl;
                return 1;
            }
        }
        cout.flush();
        return 0;
    } else if (input_type == "distance") {
        if (output_type == "json") {
            get_input_file(file_name, [&](istream& in) {
                auto distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(in);
                distance_index->write_snarls_to_json();
            });
        } else {
            cerr << "[vg view] error: (binary) Distance index can only be converted to JSON" << endl;
            return 1;
        }
        return 0;
    } else if (input_type == "snarls") {
        if (output_type == "json") {
            function<void(Snarl&)> lambda = [](Snarl& s) {
                cout << pb2json(s) << "\n";
            };
            get_input_file(file_name, [&](istream& in) {
                vg::io::for_each(in, lambda);
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
                vg::io::for_each(in, lambda);
            });
        } else {
            cerr << "[vg view] error: (binary) SnarlTraversals can only be converted to JSON" << endl;
            return 1;
        }
        return 0;
    }

    if(!graph) {
        // Make sure we didn't forget to implement an input format.
        cerr << "[vg view] error: cannot load graph in " << input_type << " format" << endl;
        return 1;
    }

    if (output_type == "gfa") {
        graph_to_gfa(graph.get(), std::cout);
        return 0;
    } 

    // Now we know graph was filled in from the input format. Spit it out in the
    // requested output format.
    
    // TODO: for now, all our output formats require a copy through vg:VG.
    // Look at the graph as a vg if possible
    VG* vg_graph = dynamic_cast<vg::VG*>(graph.get());
    
    if (vg_graph == nullptr) {
        // Copy instead. Should be fine because we on;y ever want to run this on small graphs anyway.
        vg_graph = new vg::VG();
        handlealgs::copy_path_handle_graph(graph.get(), vg_graph);
        
        // Make sure the new VG has its Proto right
        // TODO: if we didn't reach into vg.graph we wouldn't need to do this.
        vg_graph->paths.to_graph(vg_graph->graph);
        
#ifdef debug
        cerr << "Paths before conversion: " << endl;
        graph->for_each_path_handle([&](const path_handle_t& p) {
            cerr << graph->get_path_name(p) << endl;
        });
        
        cerr << "Paths after conversion: " << endl;
        vg_graph->for_each_path_handle([&](const path_handle_t& p) {
            cerr << vg_graph->get_path_name(p) << endl;
        });
        
        cerr << "VG Protobuf paths:" << endl;
        for (auto& p : vg_graph->graph.path()) {
            cerr << p.name() << endl;
        }
#endif
        
        // Give the unique_ptr ownership and delete the graph we loaded.
        graph.reset(vg_graph);
    }
    
    if(!vg_graph->is_valid()) {
        // If we're converting the graph via VG, we might as well make sure it's valid.
        // This is especially useful for JSON import.
        cerr << "[vg view] warning: graph is invalid!" << endl;
    }
    if (output_type == "dot") {
        vg_graph->to_dot(std::cout,
                        alns,
                        loci,
                        show_paths_in_dot,
                        walk_paths_in_dot,
                        annotate_paths_in_dot,
                        show_mappings_in_dot,
                        simple_dot,
                        noseq_dot,
                        invert_edge_ports_in_dot,
                        color_variants,
                        ultrabubble_labeling,
                        skip_missing_nodes,
                        ascii_labels,
                        seed_val);
    } else if (output_type == "json") {
        cout << pb2json(vg_graph->graph) << endl;
    } else if (output_type == "turtle") {
        vg_graph->to_turtle(std::cout, rdf_base_uri, color_variants);
    } else if (output_type == "vg") {
        vg_graph->serialize_to_ostream(cout);
    } else if (output_type != "gfa") {
        // We somehow got here with a bad output format.
        cerr << "[vg view] error: cannot save a graph in " << output_type << " format" << endl;
        return 1;
    }
    
    // We made it to the end and nothing broke.
    return 0;
}

// Register subcommand
static Subcommand vg_view("view", "format conversions for graphs and alignments", TOOLKIT, main_view);

