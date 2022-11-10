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
#include <bdsg/overlays/packed_reference_path_overlay.hpp>
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
         << "    -P, --path-prefix NAME   All paths (and/or GBWT threads) beginning with NAME used as reference (multiple allowed)." << endl
         << "                             Other non-ref paths not considered as samples.  When using a GBWT, select only samples with given prefix." << endl
         << "    -H, --path-sep SEP       Obtain alt paths from the set of paths, assuming a path name hierarchy (e.g. SEP='#' and sample#phase#contig)" << endl
         << "    -r, --snarls FILE        Snarls file (from vg snarls) to avoid recomputing." << endl
         << "    -g, --gbwt FILE          only consider alt traversals that correspond to GBWT threads FILE (not needed for GBZ graph input)." << endl
         << "    -T, --translation FILE   Node ID translation (as created by vg gbwt --translation) to apply to snarl names in output" << endl
         << "    -e, --path-traversals    Only consider traversals that correspond to paths in the graph." << endl
         << "    -a, --all-snarls         Process all snarls, including nested snarls (by default only top-level snarls reported)." << endl
         << "    -d, --ploidy N           Expected ploidy.  If more traversals found, they will be flagged as conflicts (default: 2)" << endl
         << "    -c, --context-jaccard N  Set context mapping size used to disambiguate alleles at sites with multiple reference traversals (default: 10000)." << endl
         << "    -u, --untangle-travs     Use context mapping to determine the reference-relative positions of each step in allele traversals (AP INFO field)." << endl
         << "    -K, --keep-conflicted    Retain conflicted genotypes in output." << endl
         << "    -S, --strict-conflicts   Drop genotypes when we have more than one haplotype for any given phase (set by default when using GBWT input)." << endl
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
    bool path_restricted_traversals = false;
    bool show_progress = false;
    int ploidy = 2;
    bool set_ploidy = false;
    bool all_snarls = false;
    bool keep_conflicted = false;
    bool strict_conflicts = false;
    int context_jaccard_window = 10000;
    bool untangle_traversals = false;
    string path_sep;
    
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
                {"path-traversals", no_argument, 0, 'e'},
                {"ploidy", required_argument, 0, 'd'},
                {"context-jaccard", required_argument, 0, 'c'},
                {"untangle-travs", no_argument, 0, 'u'},
                {"all-snarls", no_argument, 0, 'a'},
                {"keep-conflicted", no_argument, 0, 'K'},
                {"strict-conflicts", no_argument, 0, 'S'},
                {"threads", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp:P:H:r:g:T:eKSd:c:uat:v",
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
            path_sep = optarg;
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
        case 'e':
            path_restricted_traversals = true;
            break;
        case 'd':
            ploidy = parse<int>(optarg);
            set_ploidy = true;
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

    if ((!path_sep.empty() || set_ploidy) && !path_restricted_traversals && gbwt_file_name.empty()) {
        cerr << "Error [vg deconstruct]: -H and -d can only be used with -e or -g" << endl;
        return 1;
    }

    if (!gbwt_file_name.empty() && path_restricted_traversals) {
        cerr << "Error [vg deconstruct]: -e cannot be used with -g" << endl;
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

    // We might need to apply an overlay to get good path position queries
    bdsg::OverlayHelper<PathPositionHandleGraph, bdsg::PackedReferencePathOverlay, PathHandleGraph> overlay_helper;
    
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
    
    // Pre-parse some GBWT metadata
    unordered_set<string> gbwt_reference_samples;
    if (gbwt_index) {
        gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(*gbwt_index);
    }
    
    if (!refpaths.empty()) {
        // We need to inventory all the GBWT paths.
        // So we need this precomputed to access them.
        unordered_set<string> gbwt_paths;
        if (gbwt_index) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                // Get the name of this path and put it in our set.
                PathSense sense = gbwtgraph::get_path_sense(*gbwt_index, i, gbwt_reference_samples);
                gbwt_paths.insert(gbwtgraph::compose_path_name(*gbwt_index, i, sense));
            }
        }
        
        // TODO: Should we just make a GBWTGraph?
        
        // Check our paths
        for (const string& ref_path : refpaths) {
            if (!graph->has_path(ref_path) && !gbwt_paths.count(ref_path)) {
                cerr << "error [vg deconstruct]: Reference path \"" << ref_path << "\" not found in graph/gbwt" << endl;
                return 1;
            }
        }
    }
    
    if (refpaths.empty() && refpath_prefixes.empty()) {
        // No paths specified: use them all
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                const string& name = graph->get_path_name(path_handle);
                if (!Paths::is_alt(name) && PathMetadata::parse_sense(name) != PathSense::HAPLOTYPE) {
                    refpaths.push_back(name);
                }
            });
        // Add GBWT threads if no reference paths found or we're running with -a
        if (gbwt_index && (all_snarls || refpaths.empty())) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                refpaths.push_back(compose_short_path_name(*gbwt_index, i));
            }            
        }
    }

    // Read the translation
    unique_ptr<unordered_map<nid_t, pair<string, size_t>>> translation;
    if (gbz_graph.get() != nullptr) {
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

    // We store each sample ploidy specifically, based on the number of named phases
    unordered_map<string, int> sample_ploidy;
    
    // We use this to map, for example, from chromosome to genome (eg S288C.chrXVI --> S288C)
    unordered_map<string, pair<string, int>> alt_path_to_sample_phase;
    
    // process the prefixes
    if (!refpath_prefixes.empty() || !path_sep.empty()) {
        // determine phase identifiers
        map<string, set<string>> sample_seen_phases;
        graph->for_each_path_handle([&](const path_handle_t& path_handle) {
            string path_name = graph->get_path_name(path_handle);
            vector<string> vals = split(path_name, path_sep);
            if (vals.size() > 1) {
                auto& sample_name = vals[0];
                auto& phase_str = vals[1];
                /*
                if (is_number(phase_str)) {
                    seen_phases.insert(phase_str);
                }
                */
                sample_seen_phases[sample_name].insert(phase_str);
            }
        });
        map<string, map<string, int>> sample_phase_name_to_id;
        {
            for (auto& sample : sample_seen_phases) {
                int i = 0;
                for (auto& phase : sample.second) {
                    sample_phase_name_to_id[sample.first][phase] = i++;
                }
            }
        }
        /*
        if (seen_phases.size() > ploidy) {
            cerr << "Error [vg deconstruct]: We saw " << seen_phases.size()
                 << " phases, but ploidy is " << ploidy
                 << ". With -H '" << path_sep << "', the phase identifier in sample'"
                 << path_sep << "'phase1 is \"phase1\"."
                 << " You have too many unique phase identifiers in the input for the given ploidy:" << endl;
            for (auto& phase : seen_phases) {
                cerr << "phase_name: " << phase << " -> phase_id: " << phase_name_to_id[phase] << endl;
            }
            return 1;
        }
        */
        unordered_map<string, set<int>> sample_phases;
        // our phase identifiers now map into a dense range
        graph->for_each_path_handle([&](const path_handle_t& path_handle) {
                string path_name = graph->get_path_name(path_handle);
                // split on our sep
                vector<string> vals = split(path_name, path_sep);
                bool is_ref = false;
                for (auto& prefix : refpath_prefixes) {
                    if (vals[0].compare(0, prefix.size(), prefix) == 0) {
                        refpaths.push_back(path_name);
                        is_ref = true;
                        break;
                    }
                }
                if (!is_ref) {
                    auto& sample_name = vals[0];
                    int phase = 0;
                    if (vals.size() > 1) {
                        //&& is_number(vals[1])) {
                        phase = sample_phase_name_to_id[sample_name][vals[1]];
                    } else {
                        phase = 0;
                    }
                    alt_path_to_sample_phase[path_name] = make_pair(sample_name, phase);
                    sample_phases[sample_name].insert(phase);
                }
            });
        if (gbwt_index) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                std::string path_name = compose_short_path_name(*gbwt_index, i);
                for (auto& prefix : refpath_prefixes) {
                    if (path_name.compare(0, prefix.size(), prefix) == 0) {
                        refpaths.push_back(path_name);
                        break;
                    }
                }
            }
        }
        for (auto& sp : sample_phases) {
            sample_ploidy[sp.first] = sp.second.size();
        }
        if (gbwt_index) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                PathSense sense = gbwtgraph::get_path_sense(*gbwt_index, i, gbwt_reference_samples);
                string sample_name = gbwtgraph::get_path_sample_name(*gbwt_index, i, sense);
                auto phase = gbwtgraph::get_path_haplotype(*gbwt_index, i, sense);
                if (phase == PathMetadata::NO_HAPLOTYPE) {
                    // The Deconstructor defaults this to 0 so we should too.
                    phase = 0;
                }
                alt_path_to_sample_phase[sample_name] = make_pair(sample_name, (int) phase);
            }
        }
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
    dd.set_nested(all_snarls);
    dd.deconstruct(refpaths, graph, snarl_manager.get(), path_restricted_traversals, ploidy,
                   all_snarls,
                   context_jaccard_window,
                   untangle_traversals,
                   keep_conflicted,
                   strict_conflicts,
                   !alt_path_to_sample_phase.empty() ? &alt_path_to_sample_phase : nullptr,
                   &sample_ploidy,
                   gbwt_index);
    return 0;
}

// Register subcommand
static Subcommand vg_deconstruct("deconstruct", "create a VCF from variation in the graph", TOOLKIT, main_deconstruct);

