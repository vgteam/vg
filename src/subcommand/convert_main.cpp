#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "xg.hpp"
#include "../algorithms/gfa_to_handle.hpp"
#include "../io/save_handle_graph.hpp"
#include "../gfa.hpp"
#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/io/alignment_emitter.hpp>

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

//------------------------------------------------------------------------------

enum input_type { input_default, input_gam, input_gaf, input_gfa, input_gbwtgraph, input_gbz };

void help_convert(char** argv);
void no_multiple_inputs(input_type input);
void gbwtgraph_to_xg(const gbwtgraph::GBWTGraph* input, xg::XG* output, const std::string& ref_sample);
void add_paths(const gbwtgraph::GBWTGraph* input, MutablePathHandleGraph* output, const std::string& ref_sample);


//------------------------------------------------------------------------------

int main_convert(int argc, char** argv) {

    string output_format;
    input_type input = input_default;
    int64_t input_rgfa_rank = 0;
    string gfa_trans_path;
    string input_aln;
    string gbwt_name, ref_sample = gbwtgraph::REFERENCE_PATH_SAMPLE_NAME;
    set<string> rgfa_paths;
    vector<string> rgfa_prefixes;
    bool rgfa_pline = false;
    string wline_sep;

    if (argc == 2) {
        help_convert(argv);
        return 1;
    }

    constexpr int OPT_REF_SAMPLE = 1000;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"gfa-in", no_argument, 0, 'g'},
            {"in-rgfa-rank", required_argument, 0, 'r'},
            {"gbwt-in", required_argument, 0, 'b'},
            {"gbz-in", no_argument, 0, 'Z'},
            {"ref-sample", required_argument, 0, OPT_REF_SAMPLE},
            {"vg-out", no_argument, 0, 'v'},
            {"hash-out", no_argument, 0, 'a'},
            {"packed-out", no_argument, 0, 'p'},
            {"xg-out", no_argument, 0, 'x'},
            {"odgi-out", no_argument, 0, 'o'},
            {"gfa-out", no_argument, 0, 'f'},
            {"rgfa-path", required_argument, 0, 'P'},
            {"rgfa-prefix", required_argument, 0, 'Q'},
            {"rgfa-pline", no_argument, 0, 'B'},
            {"gfa-trans", required_argument, 0, 'T'},
            {"wline-sep", required_argument, 0, 'w'},
            {"gam-to-gaf", required_argument, 0, 'G'},
            {"gaf-to-gam", required_argument, 0, 'F'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "hgr:b:ZvxapxofP:Q:BT:w:G:F:t:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_convert(argv);
            return 0;
        case 'g':
            no_multiple_inputs(input);
            input = input_gfa;
            break;
        case 'r':
            input_rgfa_rank = stol(optarg);
            break;
        case 'b':
            no_multiple_inputs(input);
            input = input_gbwtgraph;
            gbwt_name = optarg;
            break;
        case 'Z':
            no_multiple_inputs(input);
            input = input_gbz;
            break;
        case OPT_REF_SAMPLE:
            ref_sample = optarg;
            break;
        case 'v':
            output_format = "vg";
            break;
        case 'a':
            output_format = "hash";
            break;
        case 'p':
            output_format = "packed";
            break;
        case 'x':
            output_format = "xg";
            break;
        case 'o':
            output_format = "odgi";
            break;
        case 'f':
            output_format = "gfa";
            break;
        case 'P':
            rgfa_paths.insert(optarg);
            break;
        case 'Q':
            rgfa_prefixes.push_back(optarg);
            break;
        case 'B':
            rgfa_pline = true;
            break;
        case 'T':
            gfa_trans_path = optarg;
            break;
        case 'w':
            wline_sep = optarg;
            break;
        case 'G':
            no_multiple_inputs(input);
            input = input_gam;
            input_aln = optarg;
            break;
        case 'F':
            no_multiple_inputs(input);
            input = input_gaf;
            input_aln = optarg;
            break;
        case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg mpmap] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
            break;
        default:
            abort();
        }
    }

    if (!gfa_trans_path.empty() && input != input_gfa) {
        cerr << "error [vg convert]: -T can only be used with -g" << endl;
        return 1;
    }
    if (output_format != "gfa" && (!rgfa_paths.empty() || !rgfa_prefixes.empty() || !wline_sep.empty())) {
        cerr << "error [vg convert]: -P, -Q, and -w can only be used with -f" << endl;
        return 1;
    }
    if (input == input_gbwtgraph || input == input_gbz) {
        if (output_format == "gfa" && !(rgfa_paths.empty() && rgfa_prefixes.empty() && wline_sep.empty())) {
            cerr << "error [vg convert]: GFA output options (-P, -Q, -w) cannot be used with a GBWTGraph" << endl;
            return 1;
        }
    }

    // with -F or -G we convert an alignment and not a graph
    if (input == input_gam || input == input_gaf) {
        if (!output_format.empty()) {
            cerr << "error [vg convert]: Alignment conversion options (-F and -G) cannot be used "
                 << "with any graph conversion options" << endl;
            return 1;
        }

        unique_ptr<HandleGraph> input_graph;
        string input_graph_filename = get_input_file_name(optind, argc, argv);
        input_graph = vg::io::VPKG::load_one<HandleGraph>(input_graph_filename);

        unique_ptr<AlignmentEmitter> emitter = get_non_hts_alignment_emitter("-", (input == input_gam) ? "GAF" : "GAM", {}, get_thread_count(),
                                                                             input_graph.get());
        std::function<void(Alignment&)> lambda = [&] (Alignment& aln) {
            emitter->emit_singles({aln});
        };                
        if (input == input_gam) {
            get_input_file(input_aln, [&](istream& in) {
                    vg::io::for_each_parallel(in, lambda);
            });
        } else {
            gaf_unpaired_for_each_parallel(*input_graph, input_aln, lambda);
        }
        return 0;
    }

    if (output_format.empty()) {
        // default to HashGraph
        output_format = "hash";
    }
        
    // allocate a graph using the graph_type string to decide a class
    unique_ptr<HandleGraph> output_graph;
    if (output_format == "vg") {
        output_graph = unique_ptr<HandleGraph>(new VG());
    } else if (output_format == "hash") {
        output_graph = unique_ptr<HandleGraph>(new bdsg::HashGraph());
    } else if (output_format == "packed") {
        output_graph = unique_ptr<HandleGraph>(new bdsg::PackedGraph());
    } else if (output_format == "xg") {
        output_graph = unique_ptr<HandleGraph>(new xg::XG());
    } else if (output_format == "odgi") {
        output_graph = unique_ptr<HandleGraph>(new bdsg::ODGI());
    } else if (output_format == "gfa") {
        // we need an intermediary for going gfa to gfa, use packed graph
        output_graph = unique_ptr<HandleGraph>(new bdsg::PackedGraph());
    }

    unique_ptr<HandleGraph> input_graph;
    unique_ptr<gbwt::GBWT> input_gbwt;
    PathHandleGraph* output_path_graph = dynamic_cast<PathHandleGraph*>(output_graph.get());

    if (input == input_gfa) {
        // we have to check this manually since we're not using the istream-based loading
        // functions in order to be able to use the disk-backed loading algorithm
        if (optind >= argc) {
            cerr << "error [vg convert]: no input graph supplied" << endl;
            return 1;
        }
        string input_stream_name = argv[optind];
        if (output_format == "xg") {
            if (input_stream_name == "-") {
                cerr << "error [vg convert]: currently cannot convert GFA from stdin to XG, try loading from disk instead" << endl;
                return 1;
            }
            dynamic_cast<xg::XG*>(output_graph.get())->from_gfa(input_stream_name);
        }
        else {
            if (input_stream_name == "-") {
                cerr << "warning [vg convert]: Converting a GFA to an indexed sequence graph from piped input can be memory intensive because the entire GFA must be loaded into memory. If memory usage is too high, consider serializing the GFA to disk to allow multiple passes over the file." << endl;
            }
            try {
                if (output_path_graph != nullptr) {
                    MutablePathMutableHandleGraph* mutable_output_graph = dynamic_cast<MutablePathMutableHandleGraph*>(output_path_graph);
                    assert(mutable_output_graph != nullptr);
                    algorithms::gfa_to_path_handle_graph(input_stream_name, mutable_output_graph,
                                                         input_stream_name != "-",
                                                         input_rgfa_rank, gfa_trans_path);
                }
                else {
                    MutableHandleGraph* mutable_output_graph = dynamic_cast<MutableHandleGraph*>(output_graph.get());
                    assert(mutable_output_graph != nullptr);
                    algorithms::gfa_to_handle_graph(input_stream_name, mutable_output_graph,
                                                    input_stream_name != "-",
                                                    gfa_trans_path);
                }
            } catch (algorithms::GFAFormatError& e) {
                cerr << "error [vg convert]: Input GFA is not acceptable." << endl;
                cerr << e.what() << endl;
                exit(1);
            } catch (std::ios_base::failure& e) {
                cerr << "error [vg convert]: IO error processing input GFA." << endl;
                cerr << e.what() << endl;
                exit(1);
            }
        }
    }
    else {
        if (input == input_gbwtgraph) {
            get_input_file(optind, argc, argv, [&](istream& in) {
                input_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(in);
            });
            gbwtgraph::GBWTGraph* gbwt_graph = dynamic_cast<gbwtgraph::GBWTGraph*>(input_graph.get());
            if (gbwt_graph == nullptr) {
                cerr << "error [vg convert]: input graph is not a GBWTGraph" << endl;
                exit(1);
            }
            input_gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name);
            gbwt_graph->set_gbwt(*input_gbwt);
        } else if (input == input_gbz) {
            std::unique_ptr<gbwtgraph::GBZ> gbz;
            get_input_file(optind, argc, argv, [&](istream& in) {
                gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(in);
            });
            input_gbwt = std::make_unique<gbwt::GBWT>(std::move(gbz->index));
            gbz->graph.set_gbwt(*input_gbwt); // The GBWT pointer will get moved to the input graph.
            input_graph = std::make_unique<gbwtgraph::GBWTGraph>(std::move(gbz->graph));
        } else {
            string input_graph_filename = get_input_file_name(optind, argc, argv);
            input_graph = vg::io::VPKG::load_one<HandleGraph>(input_graph_filename);
        }
        
        PathHandleGraph* input_path_graph = dynamic_cast<PathHandleGraph*>(input_graph.get());

        // Convert HandleGraph to HandleGraph.
        if (output_format != "gfa") {
            // XG output.
            if (output_format == "xg") {
                xg::XG* xg_graph = dynamic_cast<xg::XG*>(output_graph.get());
                if (input == input_gbwtgraph || input == input_gbz) {
                    gbwtgraph_to_xg(dynamic_cast<gbwtgraph::GBWTGraph*>(input_graph.get()), xg_graph, ref_sample);
                } else if (input_path_graph != nullptr) {
                    xg_graph->from_path_handle_graph(*input_path_graph);
                } else {
                    xg_graph->from_handle_graph(*input_graph);
                }
            }
            // GBWTGraph to PathHandleGraph.
            else if ((input == input_gbwtgraph || input == input_gbz) && output_path_graph != nullptr) {
                MutablePathMutableHandleGraph* mutable_output_graph = dynamic_cast<MutablePathMutableHandleGraph*>(output_path_graph);
                assert(mutable_output_graph != nullptr);
                // ID hint (currently only for odgi)
                mutable_output_graph->set_id_increment(input_graph->min_node_id());
                handlealgs::copy_handle_graph(input_graph.get(), mutable_output_graph);
                add_paths(dynamic_cast<gbwtgraph::GBWTGraph*>(input_graph.get()), mutable_output_graph, ref_sample);
            }
            // Normal graph to PathHandleGraph.
            else if (input_path_graph != nullptr && output_path_graph != nullptr) {
                MutablePathMutableHandleGraph* mutable_output_graph = dynamic_cast<MutablePathMutableHandleGraph*>(output_path_graph);
                assert(mutable_output_graph != nullptr);
                // ID hint (currently only for odgi)
                mutable_output_graph->set_id_increment(input_graph->min_node_id());
                handlealgs::copy_path_handle_graph(input_path_graph, mutable_output_graph);
            }
            // HandleGraph output.
            else {
                if (input_path_graph != nullptr) {
                    cerr << "warning [vg convert]: output format does not support paths, they are being dropped from the input" << endl;
                }
                MutableHandleGraph* mutable_output_graph = dynamic_cast<MutableHandleGraph*>(output_graph.get());
                assert(mutable_output_graph != nullptr);
                // ID hint (currently only for odgi)
                mutable_output_graph->set_id_increment(input_graph->min_node_id());
                handlealgs::copy_handle_graph(input_graph.get(), mutable_output_graph);
            }
        }
    }

    // GFA output.
    if (output_format == "gfa") {
        if (input == input_gbwtgraph || input == input_gbz) {
            gbwtgraph::gbwt_to_gfa(*dynamic_cast<const gbwtgraph::GBWTGraph*>(input_graph.get()), std::cout, false);
        }
        else {
            const PathHandleGraph* graph_to_write;
            if (input == input_gfa) {
                graph_to_write = dynamic_cast<const PathHandleGraph*>(output_graph.get());
            } else {
                graph_to_write = dynamic_cast<const PathHandleGraph*>(input_graph.get());
            }
            for (const string& path_name : rgfa_paths) {
                if (!graph_to_write->has_path(path_name)) {
                    cerr << "error [vg convert]: specified path, " << " not found in graph" << path_name << endl;
                    return 1;
                }
            }
            graph_to_write->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = graph_to_write->get_path_name(path_handle);
                    for (const string& prefix : rgfa_prefixes) {
                        if (path_name.substr(0, prefix.length()) == prefix) {
                            rgfa_paths.insert(path_name);
                            continue;
                        }
                    }
                });
            graph_to_gfa(graph_to_write, std::cout, rgfa_paths, rgfa_pline, wline_sep);
        }
    }
    // Serialize the output graph.
    else {
        vg::io::save_handle_graph(output_graph.get(), cout);
    }

    return 0;
}

//------------------------------------------------------------------------------

void help_convert(char** argv) {
    cerr << "usage: " << argv[0] << " convert [options] <input-graph>" << endl
         << "input options:" << endl
         << "    -g, --gfa-in           input in GFA format" << endl
         << "    -r, --in-rgfa-rank N   import rgfa tags with rank <= N as paths [default=0]" << endl
         << "    -b, --gbwt-in FILE     input graph is a GBWTGraph using the GBWT in FILE" << endl
         << "    -Z, --gbz-in           input in GBZ format (GBWT + GBWTGraph)" << endl
         << "        --ref-sample STR   convert the threads for GBWT sample STR into paths (default " << gbwtgraph::REFERENCE_PATH_SAMPLE_NAME << ")" << endl
         << "output options:" << endl
         << "    -v, --vg-out           output in VG format" << endl
         << "    -a, --hash-out         output in HashGraph format [default]" << endl
         << "    -p, --packed-out       output in PackedGraph format" << endl
         << "    -x, --xg-out           output in XG format" << endl
         << "    -o, --odgi-out         output in ODGI format" << endl
         << "    -f, --gfa-out          output in GFA format" << endl
         << "gfa options:" << endl
         << "    -P, --rgfa-path STR    write given path as rGFA tags instead of P-line (use with -f, multiple allowed, only rank-0 supported)" << endl
         << "    -Q, --rgfa-prefix STR  write paths with given prefix as rGFA tags instead of P-lines (use with -f, multiple allowed, only rank-0 supported)" << endl
         << "    -B, --rgfa-pline       paths written as rGFA tags also written as P-lines (or W-lines if selected by -w)" << endl
         << "    -T, --gfa-trans FILE   write gfa id conversions to FILE (use with -g)" << endl
         << "    -w, --wline-sep SEP    write paths with names that can be parsed as <sample><SEP><hap><SEP><contig> as GFA W-lines. (use with -f)" << endl
         << "                           suffixes of the form [start] or [start-end] will be converted into start and end coordinates if found." << endl
         << "                           ex: using \"-w .\" will convert path HG00735.1.chr1[10000-20000] to \"W HG00735 1 chr1 10000 20000\"" << endl
         << "                           multiple characters allowed (they will all be treated as separators)" << endl
         << "alignment options:" << endl
         << "    -G, --gam-to-gaf FILE  convert GAM FILE to GAF" << endl
         << "    -F, --gaf-to-gam FILE  convert GAF FILE to GAM" << endl
         << "general options:" << endl
         << "    -t, --threads N        use N threads (defaults to numCPUs)" << endl;
}

void no_multiple_inputs(input_type input) {
    if (input != input_default) {
        std::cerr << "error [vg convert]: cannot combine input types (GFA, GBWTGraph, GBZ, GAM, GAF)" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//------------------------------------------------------------------------------

void check_duplicate_contigs(const gbwt::GBWT& index, const std::vector<gbwt::size_type>& thread_ids, const std::string& ref_sample) {
    std::unordered_set<gbwt::size_type> found_contigs;
    for (auto thread_id : thread_ids) {
        gbwt::size_type contig = index.metadata.path(thread_id).contig;
        if (found_contigs.find(contig) != found_contigs.end()) {
            std::cerr << "error [vg convert]: duplicate reference path name " << index.metadata.contig(contig) << " in the GBWT index" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        found_contigs.insert(contig);
    }
}

void gbwtgraph_to_xg(const gbwtgraph::GBWTGraph* input, xg::XG* output, const std::string& ref_sample) {
    const gbwt::GBWT& index = *(input->index);
    std::vector<gbwt::size_type> thread_ids = threads_for_sample(index, ref_sample);
    if (thread_ids.empty()) {
        std::cerr << "warning [vg convert]: no threads for reference sample " << ref_sample << " in the GBWT index" << std::endl;
    }
    check_duplicate_contigs(index, thread_ids, ref_sample);

    // Enumerate nodes.
    auto for_each_sequence = [&](const std::function<void(const std::string& seq, const nid_t& node_id)>& lambda) {
        input->for_each_handle([&](const handle_t& handle) {
            lambda(input->get_sequence(handle), input->get_id(handle));
        });
    };

    // Enumerate edges.
    auto for_each_edge = [&](const std::function<void(const nid_t& from_id, const bool& from_rev,
                                                      const nid_t& to_id, const bool& to_rev)>& lambda) {
        input->for_each_edge([&](const edge_t& edge) {
                lambda(input->get_id(edge.first), input->get_is_reverse(edge.first),
                       input->get_id(edge.second), input->get_is_reverse(edge.second));
        });
    };

    // Enumerate path steps.
    auto for_each_path_element = [&](const std::function<void(const std::string& path_name,
                                                              const nid_t& node_id, const bool& is_rev,
                                                              const std::string& cigar, const bool& is_empty, const bool& is_circular)>& lambda) {
        for (gbwt::size_type thread_id : thread_ids) {
            std::string path_name = index.metadata.contig(index.metadata.path(thread_id).contig);
            gbwt::edge_type pos = index.start(gbwt::Path::encode(thread_id, false));
            while (pos.first != gbwt::ENDMARKER) {
                handle_t handle = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
                lambda(path_name, input->get_id(handle), input->get_is_reverse(handle), "", false, false);
                pos = index.LF(pos);
            }
        }
    };

    // Build XG.
    output->from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, false);
}

void add_paths(const gbwtgraph::GBWTGraph* input, MutablePathHandleGraph* output, const std::string& ref_sample) {
    const gbwt::GBWT& index = *(input->index);
    std::vector<gbwt::size_type> thread_ids = threads_for_sample(index, ref_sample);
    if (thread_ids.empty()) {
        std::cerr << "warning [vg convert]: no threads for reference sample " << ref_sample << " in the GBWT index" << std::endl;
    }
    check_duplicate_contigs(index, thread_ids, ref_sample);

    for (gbwt::size_type thread_id : thread_ids) {
        std::string path_name = index.metadata.contig(index.metadata.path(thread_id).contig);
        insert_gbwt_path(*output, index, thread_id, path_name);
    }
}

//------------------------------------------------------------------------------

// Register subcommand
static Subcommand vg_convert("convert", "convert graphs between handle-graph compliant formats as well as GFA", main_convert);
