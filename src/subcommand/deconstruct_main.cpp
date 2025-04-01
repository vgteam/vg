/** \file deconstruct_main.cpp
 *
 * Defines the "vg deconstruct" subcommand, which turns graphs back into VCFs.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../deconstructor.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gbwt_helper.hpp"
#include "../gbzgraph.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include <gbwtgraph/utils.h>
#include <gbwtgraph/index.h>

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_deconstruct(char** argv){
    cerr << "usage: " << argv[0] << " deconstruct [options] [-p|-P] <PATH> <GRAPH>" << endl
         << "Outputs VCF records for Snarls present in a graph (relative to a chosen reference path)." << endl
         << "options: " << endl
         << "    -p, --path NAME          A reference path to deconstruct against (multiple allowed)." << endl
         << "    -P, --path-prefix NAME   All paths [excluding GBWT threads / non-reference GBZ paths] beginning with NAME used as reference (multiple allowed)." << endl
         << "                             Other non-ref paths not considered as samples. " << endl
         << "    -r, --snarls FILE        Snarls file (from vg snarls) to avoid recomputing." << endl
         << "    -g, --gbwt FILE          consider alt traversals that correspond to GBWT haplotypes in FILE (not needed for GBZ graph input)." << endl
         << "    -T, --translation FILE   Node ID translation (as created by vg gbwt --translation) to apply to snarl names and AT fields in output" << endl
         << "    -O, --gbz-translation    Use the ID translation from the input gbz to apply snarl names to snarl names and AT fields in output" << endl
         << "    -a, --all-snarls         Process all snarls, including nested snarls (by default only top-level snarls reported)." << endl
         << "    -c, --context-jaccard N  Set context mapping size used to disambiguate alleles at sites with multiple reference traversals (default: 10000)." << endl
         << "    -u, --untangle-travs     Use context mapping to determine the reference-relative positions of each step in allele traversals (AP INFO field)." << endl
         << "    -K, --keep-conflicted    Retain conflicted genotypes in output." << endl
         << "    -S, --strict-conflicts   Drop genotypes when we have more than one haplotype for any given phase (set by default when using GBWT input)." << endl
         << "    -C, --contig-only-ref    Only use the CONTIG name (and not SAMPLE#CONTIG#HAPLOTYPE etc) for the reference if possible (ie there is only one reference sample)." << endl
         << "    -L, --cluster F          Cluster traversals whose (handle) Jaccard coefficient is >= F together (default: 1.0) [experimental]" << endl
         << "    -n, --nested             Write a nested VCF, including special tags. [experimental]" << endl
         << "    -f, --nested-fasta F     Write off-reference FASTA to F (and some indexing information to F.nesting.tsvx)" << endl
         << "    -R, --star-allele        Use *-alleles to denote alleles that span but do not cross the site. Only works with -n" << endl
         << "    -t, --threads N          Use N threads" << endl
         << "    -v, --verbose            Print some status messages" << endl
         << endl;
}

int main_deconstruct(int argc, char** argv){
    if (argc <= 2) {
        help_deconstruct(argv);
        return 1;
    }

    vector<string> refpaths;
    vector<string> refpath_prefixes;
    string graphname;
    string snarl_file_name;
    string gbwt_file_name;
    string translation_file_name;
    bool gbz_translation = false;
    bool path_restricted_traversals = false;
    bool show_progress = false;
    bool all_snarls = false;
    bool keep_conflicted = false;
    bool strict_conflicts = false;
    int context_jaccard_window = 10000;
    bool untangle_traversals = false;
    bool contig_only_ref = false;
    double cluster_threshold = 1.0;
    bool nested = false;
    string nested_fasta_file_name;
    bool star_allele = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"path", required_argument, 0, 'p'},
                {"path-prefix", required_argument, 0, 'P'},
                {"path-sep", required_argument, 0, 'H'},
                {"snarls", required_argument, 0, 'r'},
                {"gbwt", required_argument, 0, 'g'},
                {"translation", required_argument, 0, 'T'},
                {"gbz-translation", no_argument, 0, 'O'},                
                {"path-traversals", no_argument, 0, 'e'},
                {"context-jaccard", required_argument, 0, 'c'},
                {"untangle-travs", no_argument, 0, 'u'},
                {"all-snarls", no_argument, 0, 'a'},
                {"keep-conflicted", no_argument, 0, 'K'},
                {"strict-conflicts", no_argument, 0, 'S'},
                {"contig-only-ref", no_argument, 0, 'C'},
                {"cluster", required_argument, 0, 'L'},
                {"nested", no_argument, 0, 'n'},
                {"nested-fasta", required_argument, 0, 'f'},
                {"start-allele", no_argument, 0, 'R'},
                {"threads", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp:P:H:r:g:T:OeKSCd:c:uaL:nf:Rt:v",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            refpaths.push_back(optarg);
            break;
        case 'P':
            refpath_prefixes.push_back(optarg);
            break;
        case 'H':
            cerr << "Warning [vg deconstruct]: -H is deprecated, and will be ignored" << endl;
            break;
        case 'r':
            snarl_file_name = optarg;
            break;
        case 'g':
            gbwt_file_name = optarg;
            break;
        case 'T':
            translation_file_name = optarg;
            break;
        case 'O':
            gbz_translation = true;
            break;                        
        case 'e':
            cerr << "Warning [vg deconstruct]: -e is deprecated as it's now on default" << endl;
            break;
        case 'd':
            cerr << "Warning [vg deconstruct]: -d is deprecated - ploidy now inferred from haplotypes in path names" << endl;
            break;
            break;            
        case 'c':
            context_jaccard_window = parse<int>(optarg);
            break;
        case 'u':
            untangle_traversals = true;
            break;
        case 'a':
            all_snarls = true;
            break;
        case 'K':
            keep_conflicted = true;
            break;
        case 'S':
            strict_conflicts = true;
            break;
        case 'C':
            contig_only_ref = true;
            break;
        case 'L':
            cluster_threshold = max(0.0, min(1.0, parse<double>(optarg)));
            break;
        case 'n':
            nested = true;
            break;
        case 'f':
            nested_fasta_file_name = optarg;
            break;
        case 'R':
            star_allele = true;
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;
        case 'v':
            show_progress = true;
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

    if (nested == true && contig_only_ref == true) {
        cerr << "Error [vg deconstruct]: -C cannot be used with -n" << endl;
        return 1;
    }
    if (star_allele == true && nested == false) {
        cerr << "Error [vg deconstruct]: -R can only be used with -n" << endl;
        return 1;
    }
    
    // Read the graph
    unique_ptr<PathHandleGraph> path_handle_graph_up;
    unique_ptr<GBZGraph> gbz_graph;
    gbwt::GBWT* gbwt_index = nullptr;
    PathHandleGraph* path_handle_graph = nullptr;
    
    string path_handle_graph_filename = get_input_file_name(optind, argc, argv);    
    auto input = vg::io::VPKG::try_load_first<GBZGraph, PathHandleGraph>(path_handle_graph_filename);
    if (get<0>(input)) {        
        gbz_graph = std::move(get<0>(input));
        path_handle_graph = gbz_graph.get();
        gbwt_index = &gbz_graph->gbz.index;
    } else if (get<1>(input)) {
        path_handle_graph_up = std::move(get<1>(input));
        path_handle_graph = path_handle_graph_up.get();
    } else {
        cerr << "Error [vg deconstruct]: Input graph is not a GBZ or path handle graph" << endl;
        return 1;
    }

    if (!gbz_graph && gbz_translation) {
        cerr << "Error [vg deconstruct]: -O can only be used when input graph is in GBZ format" << endl;
    }

    if (!gbwt_file_name.empty() || gbz_graph) {
        // context jaccard depends on having steps for each alt traversal, which is
        // not something we have on hand when getting traversals from the GBWT/GBZ
        // so we toggle it off in this case to prevent outputting invalid VCFs (GTs go missing)
        // at sites with multiple reference paths
        context_jaccard_window = 0;
    }

    // We might need to apply an overlay to get good path position queries
    bdsg::ReferencePathOverlayHelper overlay_helper;
    
    // Set up to time making the overlay
    clock_t overlay_start_clock = clock();
    std::chrono::time_point<std::chrono::system_clock> overlay_start_time = std::chrono::system_clock::now(); 
    
    // Make the overlay
    PathPositionHandleGraph* graph = overlay_helper.apply(path_handle_graph);
    
    // See how long that took
    clock_t overlay_stop_clock = clock();
    std::chrono::time_point<std::chrono::system_clock> overlay_stop_time = std::chrono::system_clock::now();
    double overlay_cpu_seconds = (overlay_stop_clock - overlay_start_clock) / (double)CLOCKS_PER_SEC;
    std::chrono::duration<double> overlay_seconds = overlay_stop_time - overlay_start_time;
    
    if (show_progress && graph != dynamic_cast<PathPositionHandleGraph*>(path_handle_graph)) {
        std::cerr << "Computed overlay in " << overlay_seconds.count() << " seconds using " << overlay_cpu_seconds << " CPU seconds." << std::endl;
    }

    // Read the GBWT
    unique_ptr<gbwt::GBWT> gbwt_index_up;
    if (!gbwt_file_name.empty()) {
        if (gbwt_index) {
            cerr << "Warning [vg deconstruct]: Using GBWT from -g overrides that in input GBZ (you probably don't want to use -g)" << endl;
        }
        gbwt_index_up = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_file_name);
        if (!gbwt_index_up) {
            cerr << "Error [vg deconstruct]: Unable to load gbwt index file: " << gbwt_file_name << endl;
            return 1;
        }
        gbwt_index = gbwt_index_up.get();
    }
        
    if (!refpaths.empty()) {
        // Check our paths
        for (const string& ref_path : refpaths) {
            if (!graph->has_path(ref_path)) {
                cerr << "error [vg deconstruct]: Reference path \"" << ref_path << "\" not found in graph/gbwt" << endl;
                return 1;
            }
        }
    }
    
    if (refpaths.empty() && refpath_prefixes.empty()) {
        bool found_hap;
        // No paths specified: use them all
        graph->for_each_path_handle([&](path_handle_t path_handle) {
            const string& name = graph->get_path_name(path_handle);
            if (!Paths::is_alt(name) && PathMetadata::parse_sense(name) != PathSense::HAPLOTYPE) {
                refpaths.push_back(name);
            } else {
                found_hap = true;
            }
        });

        if (!found_hap && gbwt_index == nullptr) {
            cerr << "error [vg deconstruct]: All graph paths selected as references (leaving no alts). Please use -P/-p "
                 << "to narrow down the reference to a subset of paths, or GBZ/GBWT input that contains haplotype paths" << endl;
            return 1;
        }        
    }

    // Read the translation
    unique_ptr<unordered_map<nid_t, pair<string, size_t>>> translation;
    if (gbz_graph.get() != nullptr && gbz_translation) {
        // try to get the translation from the graph
        translation = make_unique<unordered_map<nid_t, pair<string, size_t>>>();
        *translation = load_translation_back_map(gbz_graph->gbz.graph);
        if (translation->empty()) {
            // not worth keeping an empty translation
            translation = nullptr;
        }
    }
    if (!translation_file_name.empty()) {
        if (!translation->empty()) {
            cerr << "Warning [vg deconstruct]: Using translation from -T overrides that in input GBZ (you probably don't want to use -T)" << endl;
        }
        ifstream translation_file(translation_file_name.c_str());
        if (!translation_file) {
            cerr << "Error [vg deconstruct]: Unable to load translation file: " << translation_file_name << endl;
            return 1;
        }
        translation = make_unique<unordered_map<nid_t, pair<string, size_t>>>();
        *translation = load_translation_back_map(*graph, translation_file);
    }
    
    // Load or compute the snarls
    unique_ptr<SnarlManager> snarl_manager;    
    if (!snarl_file_name.empty()) {
        ifstream snarl_file(snarl_file_name.c_str());
        if (!snarl_file) {
            cerr << "Error [vg deconstruct]: Unable to load snarls file: " << snarl_file_name << endl;
            return 1;
        }
        if (show_progress) {
            cerr << "Loading snarls" << endl;
        }
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
    } else {
        IntegratedSnarlFinder finder(*graph);
        if (show_progress) {
            cerr << "Finding snarls" << endl;
        }
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
    }
    
    // process the prefixes to find ref paths
    if (!refpath_prefixes.empty()) {
        graph->for_each_path_handle([&](const path_handle_t& path_handle) {
                string path_name = graph->get_path_name(path_handle);
                for (auto& prefix : refpath_prefixes) {                    
                    if (path_name.compare(0, prefix.size(), prefix) == 0) {
                        refpaths.push_back(path_name);
                        break;
                    }
                }
            });
    }

    if (refpaths.empty()) {
        cerr << "Error [vg deconstruct]: No specified reference path or prefix found in graph" << endl;
        return 1;
    }

#ifdef USE_CALLGRIND
    // We want to profile stuff that accesses paths, not the loading.
    CALLGRIND_START_INSTRUMENTATION;
#endif

    // Deconstruct
    Deconstructor dd;
    if (show_progress) {
        cerr << "Deconstructing top-level snarls" << endl;
    }
    dd.set_translation(translation.get());
    dd.set_nested(all_snarls || nested);
    dd.deconstruct(refpaths, graph, snarl_manager.get(),
                   all_snarls,
                   context_jaccard_window,
                   untangle_traversals,
                   keep_conflicted,
                   strict_conflicts,
                   !contig_only_ref,
                   cluster_threshold,
                   gbwt_index,
                   nested,
                   star_allele);

    if (!nested_fasta_file_name.empty()) {
        dd.save_off_ref_sequences(nested_fasta_file_name);
    }
    return 0;
}

// Register subcommand
static Subcommand vg_deconstruct("deconstruct", "create a VCF from variation in the graph", TOOLKIT, main_deconstruct);

