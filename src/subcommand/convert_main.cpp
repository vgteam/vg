#include "subcommand.hpp"
#include "../vg.hpp"
#include "../utility.hpp"
#include "xg.hpp"
#include "../algorithms/gfa_to_handle.hpp"
#include "../algorithms/find_gbwtgraph.hpp"
#include "../io/save_handle_graph.hpp"
#include "../gfa.hpp"
#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <vg/io/alignment_emitter.hpp>

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/gfa.h>

#include <unistd.h>
#include <getopt.h>

using namespace vg;
using namespace vg::subcommand;
using namespace vg::io;

const string context = "[vg convert]";

//------------------------------------------------------------------------------

// We need a type for describing what kind of input to parse.
enum input_type { input_handlegraph, input_gam, input_gaf, input_gfa, input_gbwtgraph };
const input_type INPUT_DEFAULT = input_handlegraph;

// We also need a type for a tri-state for deciding what kind of GFA output algorithm to use.
enum algorithm_type { algorithm_auto, algorithm_vg, algorithm_gbwtgraph };
const algorithm_type ALGORITHM_DEFAULT = algorithm_auto;

void help_convert(char** argv);
void no_multiple_inputs(input_type input);
// Generate an XG with nodes, edges, and paths from input.
// Promote haplotype-sense paths for the samples in ref_samples to reference sense.
// Promote generic-sense path with the locus hap_locus to haplotype sense.
// Copy across other haplotype-sense paths if unless drop_haplotypes is true.
void graph_to_xg_adjusting_paths(const PathHandleGraph* input, xg::XG* output,
                                 const std::unordered_set<std::string>& ref_samples,
                                 const std::string& hap_locus, const std::string& new_sample, bool drop_haplotypes);
// Copy paths from input to output.
// Promote haplotype-sense paths for the samples in ref_samples to reference sense.
// Promote generic-sense paths with the locus hap_locus to haplotype sense.
// Copy across other haplotype-sense paths if unless drop_haplotypes is true.
void add_and_adjust_paths(const PathHandleGraph* input, MutablePathHandleGraph* output, 
                          const std::unordered_set<std::string>& ref_samples, 
                          const std::string& hap_locus, const std::string& new_sample, bool drop_haplotypes);


//------------------------------------------------------------------------------

int main_convert(int argc, char** argv) {

    string output_format;
    input_type input = INPUT_DEFAULT;
    int64_t input_rgfa_rank = 0;
    string gfa_trans_path;
    string input_aln;
    string gbwt_name;
    unordered_set<string> ref_samples;
    string hap_locus;
    string new_sample;
    bool drop_haplotypes = false;
    set<string> rgfa_paths;
    vector<string> rgfa_prefixes;
    bool rgfa_pline = false;
    bool wline = true;
    algorithm_type gfa_output_algorithm = ALGORITHM_DEFAULT;

    // For GBWTGraph to GFA.
    int num_threads = omp_get_max_threads();
    bool use_translation = true;

    if (argc == 2) {
        help_convert(argv);
        return 1;
    }

    constexpr int OPT_REF_SAMPLE = 1000;
    constexpr int OPT_GBWTGRAPH_ALGORITHM = 1001;
    constexpr int OPT_VG_ALGORITHM = 1002;
    constexpr int OPT_NO_TRANSLATION = 1003;
    constexpr int OPT_HAP_LOCUS = 1004;
    constexpr int OPT_NEW_SAMPLE = 1005;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"gfa-in", no_argument, 0, 'g'},
            {"in-rgfa-rank", required_argument, 0, 'r'},
            {"gbwt-in", required_argument, 0, 'b'},
            {"ref-sample", required_argument, 0, OPT_REF_SAMPLE},
            {"hap-locus", required_argument, 0, OPT_HAP_LOCUS},
            {"new-sample", required_argument, 0, OPT_NEW_SAMPLE},
            {"drop-haplotypes", no_argument, 0, 'H'},
            {"vg-out", no_argument, 0, 'v'},
            {"hash-out", no_argument, 0, 'a'},
            {"packed-out", no_argument, 0, 'p'},
            {"xg-out", no_argument, 0, 'x'},
            {"gfa-out", no_argument, 0, 'f'},
            {"rgfa-path", required_argument, 0, 'P'},
            {"rgfa-prefix", required_argument, 0, 'Q'},
            {"rgfa-pline", no_argument, 0, 'B'},
            {"gfa-trans", required_argument, 0, 'T'},
            {"no-wline", no_argument, 0, 'W'},
            {"gbwtgraph-algorithm", no_argument, 0, OPT_GBWTGRAPH_ALGORITHM},
            {"vg-algorithm", no_argument, 0, OPT_VG_ALGORITHM},
            {"no-translation", no_argument, 0, OPT_NO_TRANSLATION},
            {"gam-to-gaf", required_argument, 0, 'G'},
            {"gaf-to-gam", required_argument, 0, 'F'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}

        };
        int option_index = 0;
        c = getopt_long (argc, argv, "h?gr:b:HvxapfP:Q:BT:WG:F:t:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case '?':
        case 'h':
            help_convert(argv);
            return 1;
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
            gbwt_name = error_if_file_does_not_exist(context, optarg);
            break;
        case OPT_REF_SAMPLE:
            ref_samples.insert(optarg);
            break;
        case OPT_HAP_LOCUS:
            hap_locus = optarg;
            break;
        case OPT_NEW_SAMPLE:
            new_sample = optarg;
            break;
        case 'H':
            drop_haplotypes = true;
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
            gfa_trans_path = error_if_file_cannot_be_written(context, optarg);
            break;
        case 'W':
            wline = false;
            break;
        case OPT_GBWTGRAPH_ALGORITHM:
            gfa_output_algorithm = algorithm_gbwtgraph;
            break;
        case OPT_VG_ALGORITHM:
            gfa_output_algorithm = algorithm_vg;
            break;
        case OPT_NO_TRANSLATION:
            use_translation = false;
            break;
        case 'G':
            no_multiple_inputs(input);
            input = input_gam;
            input_aln = error_if_file_does_not_exist(context, optarg);
            break;
        case 'F':
            no_multiple_inputs(input);
            input = input_gaf;
            input_aln = error_if_file_does_not_exist(context, optarg);
            break;
        case 't':
            omp_set_num_threads(parse_thread_count(context, optarg));
            break;
        default:
            abort();
        }
    }
    
    if (!gfa_trans_path.empty() && input != input_gfa) {
        error_and_exit(context, "-T can only be used with -g");
    }
    if (output_format != "gfa" && (!rgfa_paths.empty() || !rgfa_prefixes.empty() || !wline)) {
        error_and_exit(context, "-P, -Q, and -W can only be used with -f");
    }
    if (gfa_output_algorithm == algorithm_gbwtgraph) {
        if (output_format != "gfa") {
            error_and_exit(context, "Only GFA output format can be used "
                                     "with the GBWTGraph library GFA conversion algorithm");
        }
        if (input == input_gfa) {
            error_and_exit(context, "GFA input cannot be used "
                                    "with the GBWTGraph library GFA conversion algorithm");
        }
        if (!(rgfa_paths.empty() && rgfa_prefixes.empty() && wline)) {
            error_and_exit(context, "GFA output options (-P, -Q, -W) cannot be used "
                                     "with the GBWTGraph library GFA conversion algorithm");
        }
    }
    if (output_format == "gfa" && !ref_samples.empty()) {
        error_and_exit(context, "paths cannot be converted to reference sense when writing GFA output");
    }
    if (output_format == "gfa" && !hap_locus.empty()) {
        error_and_exit(context, "paths cannot be converted to haplotype sense when writing GFA output");
    }
    if (!ref_samples.empty() && !hap_locus.empty()) {
        error_and_exit(context, "cannot convert paths to haplotype and reference sense at the same time");
    }
    if (new_sample.empty() != hap_locus.empty()) {
        error_and_exit(context, "to convert a generic-sense path to a haplotype, "
                                "both --hap-locus and --new-sample are required");
    }
    if (output_format == "vg") {
        emit_warning(context, "vg-protobuf output (-v / --vg-out) is deprecated. "
                              "Please use -p instead");
    }

    
    // with -F or -G we convert an alignment and not a graph
    if (input == input_gam || input == input_gaf) {
        if (!output_format.empty()) {
            error_and_exit(context, "Alignment conversion options (-F and -G) "
                                    "cannot be used with any graph conversion options");
        }

        unique_ptr<HandleGraph> input_graph;
        string input_graph_filename = get_input_file_name(optind, argc, argv);
        input_graph = vg::io::VPKG::load_one<HandleGraph>(input_graph_filename);

        unique_ptr<AlignmentEmitter> emitter = get_non_hts_alignment_emitter(
            "-", (input == input_gam) ? "GAF" : "GAM", {}, get_thread_count(), input_graph.get());
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
        // default to PackedGraph
        output_format = "packed";
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
    } else if (output_format == "gfa") {
        // we need an intermediary for going gfa to gfa, use packed graph
        output_graph = unique_ptr<HandleGraph>(new bdsg::PackedGraph());
    }
    PathHandleGraph* output_path_graph = dynamic_cast<PathHandleGraph*>(output_graph.get());

    unique_ptr<HandleGraph> input_graph;
    unique_ptr<gbwt::GBWT> input_gbwt;
    
    if (input == input_gfa) {
        // we have to check this manually since we're not using the istream-based loading
        // functions in order to be able to use the disk-backed loading algorithm
        if (optind >= argc) {
            error_and_exit(context, "no input graph supplied");
        }
        string input_stream_name = argv[optind];
        if (output_format == "xg") {
            xg::XG* xg_graph = dynamic_cast<xg::XG*>(output_graph.get());
            
            // Need to go through a handle graph
            bdsg::HashGraph intermediate;
            emit_warning(context, "currently cannot convert GFA directly to XG; "
                                  "converting through another format");
            algorithms::gfa_to_path_handle_graph(input_stream_name, &intermediate,
                                                 input_rgfa_rank, gfa_trans_path);
            graph_to_xg_adjusting_paths(&intermediate, xg_graph, ref_samples, hap_locus, new_sample, drop_haplotypes);
        }
        else {
            // If the GFA doesn't have forward references, we can handle it
            // efficiently even if we are streaming it, so we shouldn't warn about
            // "-" input here.
            try {
                if (output_path_graph != nullptr) {
                    MutablePathMutableHandleGraph* mutable_output_graph \
                        = dynamic_cast<MutablePathMutableHandleGraph*>(output_path_graph);
                    assert(mutable_output_graph != nullptr);
                    algorithms::gfa_to_path_handle_graph(input_stream_name, mutable_output_graph,
                                                         input_rgfa_rank, gfa_trans_path);
                }
                else {
                    MutableHandleGraph* mutable_output_graph = dynamic_cast<MutableHandleGraph*>(output_graph.get());
                    assert(mutable_output_graph != nullptr);
                    algorithms::gfa_to_handle_graph(input_stream_name, mutable_output_graph,
                                                    gfa_trans_path);
                }
            } catch (algorithms::GFAFormatError& e) {
                error_and_exit(context, "Input GFA is not acceptable.\n" + string(e.what()));
            } catch (std::ios_base::failure& e) {
                error_and_exit(context, "IO error processing input GFA.\n" + string(e.what()));
            }
        }
    }
    else {
        if (input == input_gbwtgraph) {
            // We need to read the input as a GBWTGraph file and attach it to a GBWT.
            get_input_file(optind, argc, argv, [&](istream& in) {
                input_graph = vg::io::VPKG::load_one<gbwtgraph::GBWTGraph>(in);
            });
            gbwtgraph::GBWTGraph* gbwt_graph = dynamic_cast<gbwtgraph::GBWTGraph*>(input_graph.get());
            if (gbwt_graph == nullptr) {
                error_and_exit(context, "input graph is not a GBWTGraph");
            }
            input_gbwt = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_name);
            gbwt_graph->set_gbwt(*input_gbwt);
        } else if (input == input_handlegraph) {
            string input_graph_filename = get_input_file_name(optind, argc, argv);
            input_graph = vg::io::VPKG::load_one<HandleGraph>(input_graph_filename);
        } else {
            throw std::runtime_error("Unimplemented input type");
        }
        
        PathHandleGraph* input_path_graph = dynamic_cast<PathHandleGraph*>(input_graph.get());

        // Convert HandleGraph to HandleGraph.
        if (output_format != "gfa") {
            // XG output.
            if (output_format == "xg") {
                xg::XG* xg_graph = dynamic_cast<xg::XG*>(output_graph.get());
                if (input_path_graph != nullptr) {
                    // We can convert to XG with paths, which we might adjust
                    graph_to_xg_adjusting_paths(input_path_graph, xg_graph, ref_samples, 
                                               hap_locus, new_sample, drop_haplotypes);
                } else {
                    // No paths, just convert to xg without paths
                    xg_graph->from_handle_graph(*input_graph);
                }
            }
            // PathHandleGraph (possibly with haplotypes) to PathHandleGraph.
            else if (input_path_graph != nullptr && output_path_graph != nullptr) {
                MutablePathMutableHandleGraph* mutable_output_graph \
                    = dynamic_cast<MutablePathMutableHandleGraph*>(output_path_graph);
                assert(mutable_output_graph != nullptr);
                // ID hint (TODO: currently not needed since only odgi used it)
                mutable_output_graph->set_id_increment(input_graph->min_node_id());
                // Copy the graph as-is
                handlealgs::copy_handle_graph(input_graph.get(), mutable_output_graph);
                // Copy the paths across with possibly some rewriting
                add_and_adjust_paths(input_path_graph, mutable_output_graph, ref_samples, 
                                     hap_locus, new_sample, drop_haplotypes);
            }
            // HandleGraph output.
            else {
                if (input_path_graph != nullptr) {
                    emit_warning(context, "output format does not support paths; "
                                          "they are being dropped from the input");
                }
                MutableHandleGraph* mutable_output_graph = dynamic_cast<MutableHandleGraph*>(output_graph.get());
                assert(mutable_output_graph != nullptr);
                // ID hint (TODO: currently not needed since only odgi used it)
                mutable_output_graph->set_id_increment(input_graph->min_node_id());
                handlealgs::copy_handle_graph(input_graph.get(), mutable_output_graph);
            }
        }
    }

    // GFA output.
    if (output_format == "gfa") {
        if (gfa_output_algorithm == algorithm_auto) {
            // Determine algorithm to use.
            if (!rgfa_paths.empty() || !rgfa_prefixes.empty() || rgfa_pline || !wline) {
                // We've asked for special conversion options that only the vg algorithm supports.
                gfa_output_algorithm = algorithm_vg;
            } else if (vg::algorithms::find_gbwtgraph(input_graph.get())) {
                // There's a GBWTGraph available so use that algorithm.
                gfa_output_algorithm = algorithm_gbwtgraph;
            } else {
                // No GBWTGraph is available so use the VG algorithm.
                gfa_output_algorithm = algorithm_vg;
            }
        }
        if (gfa_output_algorithm == algorithm_gbwtgraph) {
            // We need to find a GBWTGraph to use for this
            const gbwtgraph::GBWTGraph* gbwt_graph = vg::algorithms::find_gbwtgraph(input_graph.get());
            if (gbwt_graph == nullptr) {
                error_and_exit(context, "input graph does not have a GBWTGraph, "
                                        "so GBWTGraph library GFA conversion algorithm cannot be used");
            }
            
            gbwtgraph::GFAExtractionParameters parameters;
            parameters.num_threads = num_threads;
            parameters.use_translation = use_translation;
            gbwtgraph::gbwt_to_gfa(*gbwt_graph, std::cout, parameters);
        } else if (gfa_output_algorithm == algorithm_vg) {
            // Use HandleGraph GFA conversion code
            const PathHandleGraph* graph_to_write;
            if (input == input_gfa) {
                graph_to_write = dynamic_cast<const PathHandleGraph*>(output_graph.get());
            } else {
                graph_to_write = dynamic_cast<const PathHandleGraph*>(input_graph.get());
            }
            for (const string& path_name : rgfa_paths) {
                if (!graph_to_write->has_path(path_name)) {
                    error_and_exit(context, "specified path, " + path_name + ", not found in graph");
                }
            }
            if (!rgfa_prefixes.empty()) {
                graph_to_write->for_each_path_matching({}, {}, {}, [&](path_handle_t path_handle) {
                    // Scan for any paths of any sense matching an rGFA prefix.
                    string path_name = graph_to_write->get_path_name(path_handle);
                    for (const string& prefix : rgfa_prefixes) {
                        if (path_name.substr(0, prefix.length()) == prefix) {
                            rgfa_paths.insert(path_name);
                            continue;
                        }
                    }
                });
            }
            graph_to_gfa(graph_to_write, std::cout, rgfa_paths, rgfa_pline, wline);
        } else {
            throw std::runtime_error("Unimplemented GFA output algorithm");
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
         << "  -g, --gfa-in               input in GFA format" << endl
         << "  -r, --in-rgfa-rank N       import rgfa tags with rank <= N as paths [0]" << endl
         << "  -b, --gbwt-in FILE         input graph is a GBWTGraph using the GBWT in FILE" << endl
         << "      --ref-sample STR       change haplotypes for this sample to" << endl
         << "                             reference paths (may repeat)" << endl
         << "      --hap-locus STR        change generic paths with this locus" << endl
         << "                             to haplotype paths (must be used with --new-sample)" << endl
         << "      --new-sample STR       when using --hap-locus, give the new haplotype" << endl
         << "                             this sample name (must be used with --hap-locus)" << endl
         << "gfa input options (use with -g):" << endl
         << "  -T, --gfa-trans FILE       write gfa id conversions to FILE" << endl
         << "output options:" << endl
         << "  -v, --vg-out               output in VG's original Protobuf format" << endl
         << "                             [DEPRECATED: use -p instead]." << endl
         << "  -a, --hash-out             output in HashGraph format" << endl
         << "  -p, --packed-out           output in PackedGraph format (default)" << endl
         << "  -x, --xg-out               output in XG format" << endl
         << "  -f, --gfa-out              output in GFA format" << endl
         << "  -H, --drop-haplotypes      do not include haplotype paths in the output" << endl
         << "                             (useful with GBWTGraph / GBZ inputs)" << endl
         << "gfa output options (use with -f):" << endl
         << "  -P, --rgfa-path STR        write given path as rGFA tags instead of lines" << endl
         << "                             (may repeat, only rank-0 supported)" << endl
         << "  -Q, --rgfa-prefix STR      write paths with this prefix as rGFA tags instead" << endl
         << "                             of lines (may repeat, only rank-0 supported)" << endl
         << "  -B, --rgfa-pline           paths written as rGFA tags also written as lines" << endl
         << "  -W, --no-wline             write all paths as GFA P-lines instead of W-lines." << endl
         << "                             allows handling multiple phase blocks " << endl
         << "                             and subranges used together." << endl
         << "      --gbwtgraph-algorithm  always use the GBWTGraph library GFA algorithm." << endl
         << "                             not compatible with other GFA output options" << endl
         << "                             or non-GBWT graphs." << endl
         << "      --vg-algorithm         always use the VG GFA algorithm. Works with all" << endl
         << "                             options and graph types, but can't preserve" << endl
         << "                             original GFA coordinates" << endl
         << "      --no-translation       when using the GBWTGraph algorithm, convert graph" << endl
         << "                             directly to GFA; do not use the translation" << endl
         << "                             to preserve original coordinates" << endl
         << "alignment options:" << endl
         << "  -G, --gam-to-gaf FILE      convert GAM FILE to GAF" << endl
         << "  -F, --gaf-to-gam FILE      convert GAF FILE to GAM" << endl
         << "general options:" << endl
         << "  -t, --threads N            use N threads [numCPUs]" << endl
         << "  -h, --help                 print this help message to stderr and exit" << endl;
}

void no_multiple_inputs(input_type input) {
    if (input != INPUT_DEFAULT) {
        error_and_exit(context, "cannot combine input types  (GFA, GBWTGraph, GBZ, GAM, GAF)");
    }
}

//------------------------------------------------------------------------------

/// Check to make sure the haplotype paths with sample names in the given set
/// have no more than one phase block per sample/haplotype/contig combination.
/// Also return a map from sample name to set of observed haplotype numbers (so
/// we can include them in reference-sense paths only if needed).
std::unordered_map<std::string, std::unordered_set<int64_t>> check_duplicate_path_names(
    const PathHandleGraph* input, const std::unordered_set<std::string>& ref_samples) {
    // Check to make sure no ref samples have fragmented haplotypes. If they
    // do, we can't drop the phase block and can't change the sense to
    // reference. Store set of phase blocks by sample, haplotype, contig.
    // If we stored just counts, we couldn't handle multiple subranges on the
    // same phase block properly.
    std::unordered_map<std::tuple<std::string, int64_t, std::string>, std::unordered_set<int64_t>> phase_block_sets;
    // Also determine whether to strip the haplotype numbers; if there are
    // multiple haplotypes stored for a sample (and the stored set exceeds size
    // 1) we will keep them.
    std::unordered_map<std::string, std::unordered_set<int64_t>> sample_to_haplotypes;
    if (!ref_samples.empty()) {
        input->for_each_path_matching({PathSense::HAPLOTYPE}, ref_samples, {}, [&](const path_handle_t& path) {
            // For each path in these samples' haplotypes...
            
            auto sample = input->get_sample_name(path);
            auto haplotype = input->get_haplotype(path);
            auto contig = input->get_locus_name(path);
            auto phase_block = input->get_phase_block(path);
            
            // Find the place to remember phase blocks for it
            auto& phase_block_set = phase_block_sets[
                std::tuple<std::string, int64_t, std::string>(sample, haplotype, contig)];
            
            // Insert the phase block
            phase_block_set.insert(phase_block);
            
            if (phase_block_set.size() > 1) {
                // We can't resolve these.
                error_and_exit(context, "multiple phase blocks on sample " + sample +
                                        " haplotype " + std::to_string(haplotype) +
                                        " contig " + contig +
                                        " prevent promoting the sample to a reference");
            }
            
            // Log its haplotypes
            sample_to_haplotypes[sample].insert(haplotype);
        });
    }
    
    return sample_to_haplotypes;
}

void graph_to_xg_adjusting_paths(const PathHandleGraph* input, xg::XG* output, 
                                 const std::unordered_set<std::string>& ref_samples, 
                                 const std::string& hap_locus, const std::string& new_sample, bool drop_haplotypes) {
    // Building an XG uses a slightly different interface, so we duplicate some
    // code from the normal MutablePathMutableHandleGraph build.
    // TODO: Find a way to unify the duplicated code?
    
    // Make sure we can safely promote any haplotypes to reference, and get the
    // information we need to determine if we need to keep haplotype numbers
    // when doing so.
    auto sample_to_haplotypes = check_duplicate_path_names(input, ref_samples);

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
                                                              const std::string& cigar, const bool& is_empty, 
                                                              const bool& is_circular)>& lambda) {
        
        // Define a function to copy over a path.
        // XG constructuon relies on name-encoded path metadata.
        auto copy_path = [&](const path_handle_t& path, const std::string new_name) {
            bool is_circular = input->get_is_circular(path);
            for (handle_t handle : input->scan_path(path)) {
                lambda(new_name, input->get_id(handle), input->get_is_reverse(handle), "", false, is_circular);
            }
            // TODO: Should we preserve empty paths here?
        };
        
        // Copy over the existing generic paths, except ones that get promoted to haplotype paths
        input->for_each_path_matching({PathSense::GENERIC}, {}, {}, [&](const path_handle_t& path) {
            if (hap_locus != input->get_locus_name(path)) {
                copy_path(path, input->get_path_name(path));
            }
        });

        // Copy over the existing reference paths
        input->for_each_path_matching({PathSense::REFERENCE}, {}, {}, [&](const path_handle_t& path) {
            copy_path(path, input->get_path_name(path));
        });
        
        if (!ref_samples.empty()) {
            // Copy all haplotype paths matching the ref samples as reference
            input->for_each_path_matching({PathSense::HAPLOTYPE}, ref_samples, {}, [&](const path_handle_t& path) {
                
                // Compose the new reference-ified metadata
                std::string sample = input->get_sample_name(path);
                std::string locus = input->get_locus_name(path);
                // We should always preserve the haplotype phase number; we
                // will need it if we ever want to go back to haplotype sense.
                int64_t haplotype = input->get_haplotype(path);
                auto subrange = input->get_subrange(path);
                
                // Make a new name with reference-ified metadata.
                // Phase block is safe to discard because we checked for duplicates without it.
                auto new_name = PathMetadata::create_path_name(PathSense::REFERENCE,
                                                               sample,
                                                               locus,
                                                               haplotype,
                                                               PathMetadata::NO_PHASE_BLOCK,
                                                               subrange);
                
                // Copy out to the xg
                copy_path(path, new_name);
            });
        }

        if (!hap_locus.empty()) {
            // Copy the generic paths matching the hap samples as haplotypes
            input->for_each_path_matching({PathSense::GENERIC}, {}, {hap_locus}, [&](const path_handle_t& path) {
                
                // Compose the new haplotype-ified metadata
                std::string sample = new_sample;
                std::string locus = input->get_locus_name(path);
                // We should always preserve the haplotype phase number; we
                // will need it if we ever want to go back to haplotype sense.
                int64_t haplotype = input->get_haplotype(path);
                if (haplotype == -1) {
                    // If there is no haplotype, make it 0
                    haplotype = 0;
                }
                auto phase_block = input->get_phase_block(path);
                auto subrange = input->get_subrange(path);
                
                // Make a new name with haplotype-ified metadata.
                auto new_name = PathMetadata::create_path_name(PathSense::HAPLOTYPE,
                                                               sample,
                                                               locus,
                                                               haplotype,
                                                               phase_block,
                                                               subrange);
                // Copy out to the xg
                copy_path(path, new_name);
            });
        }
        
        if (!drop_haplotypes) {
            // Copy across any other haplotypes.
            input->for_each_path_matching({PathSense::HAPLOTYPE}, {}, {}, [&](const path_handle_t& path) {
                if (ref_samples.count(input->get_sample_name(path))) {
                    // Skip those we already promoted to reference sense
                    return;
                }
                copy_path(path, input->get_path_name(path));
            });
        }
    };

    // Build XG.
    output->from_enumerators(for_each_sequence, for_each_edge, for_each_path_element, false);
}

void add_and_adjust_paths(const PathHandleGraph* input, MutablePathHandleGraph* output, 
                          const std::unordered_set<std::string>& ref_samples, const std::string& hap_locus, 
                          const std::string& new_sample, bool drop_haplotypes) {
    
    // Make sure we aren't working with fragmented haplotypes that can't convert to reference sense.
    auto sample_to_haplotypes = check_duplicate_path_names(input, ref_samples);

    // Copy over the existing generic paths, except ones that get promoted to haplotype paths
    input->for_each_path_matching({PathSense::GENERIC}, {}, {}, [&](const path_handle_t& path) {
        if (hap_locus != input->get_locus_name(path)) {
        handlegraph::algorithms::copy_path(input, path, output);
        }
    });
    
    // Copy all reference paths that exist already
    input->for_each_path_matching({PathSense::REFERENCE}, {}, {}, [&](const path_handle_t& path) {
        handlegraph::algorithms::copy_path(input, path, output);
    });
    
    if (!ref_samples.empty()) {
        // Copy all haplotype paths matching the ref samples as reference
        input->for_each_path_matching({PathSense::HAPLOTYPE}, ref_samples, {}, [&](const path_handle_t& path) {
            
            // Compose the new reference-ified metadata
            std::string sample = input->get_sample_name(path);
            std::string locus = input->get_locus_name(path);
            // We should always preserve the haplotype phase number; we
            // will need it if we ever want to go back to haplotype sense.
            int64_t haplotype = input->get_haplotype(path);
            auto subrange = input->get_subrange(path);
            bool is_circular = input->get_is_circular(path);
            
            // Make a new path with reference-ified metadata.
            // Phase block is safe to discard because we checked for duplicates without it.
            path_handle_t into_path = output->create_path(PathSense::REFERENCE,
                                                          sample,
                                                          locus,
                                                          haplotype,
                                                          PathMetadata::NO_PHASE_BLOCK,
                                                          subrange,
                                                          is_circular);
            
            // Copy across the steps
            handlegraph::algorithms::copy_path(input, path, output, into_path);
        });
    }
    if (!hap_locus.empty()) {
        // Copy all generic paths matching the hap samples as haplotypes
        input->for_each_path_matching({PathSense::GENERIC}, {}, {hap_locus}, [&](const path_handle_t& path) {
            
            // Compose the new reference-ified metadata
            std::string sample = new_sample;
            std::string locus = input->get_locus_name(path);
            // We should always preserve the haplotype phase number; we
            // will need it if we ever want to go back to haplotype sense.
            int64_t haplotype = input->get_haplotype(path);
            if (haplotype == -1) {
                haplotype = 0;
            }
            auto phase_block = input->get_phase_block(path);
            if (phase_block == std::numeric_limits<size_t>::max()) {
                phase_block = 0;
            }
            auto subrange = input->get_subrange(path);
            bool is_circular = input->get_is_circular(path);
            
            // Make a new path with haplotype-ified metadata.
            path_handle_t into_path = output->create_path(PathSense::HAPLOTYPE,
                                                          sample,
                                                          locus,
                                                          haplotype,
                                                          phase_block,
                                                          subrange,
                                                          is_circular);
            
            // Copy across the steps
            handlegraph::algorithms::copy_path(input, path, output, into_path);
        });
    }

    if (!drop_haplotypes) {
        // Copy across any other haplotypes.
        input->for_each_path_matching({PathSense::HAPLOTYPE}, {}, {}, [&](const path_handle_t& path) {
            if (ref_samples.count(input->get_sample_name(path))) {
                // Skip those we already promoted to reference sense
                return;
            }
            handlegraph::algorithms::copy_path(input, path, output);
        });
    }
}

//------------------------------------------------------------------------------

// Register subcommand
static Subcommand vg_convert("convert", "convert graphs between handle-graph compliant formats as well as GFA", main_convert);
