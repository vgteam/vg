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

#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;


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
                cerr << "[vg view] error: (binary) GAM can only be converted to JSON or FASTQ" << endl;
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

// Register subcommand
static Subcommand vg_view("view", "format conversions for graphs and alignments", main_view);

