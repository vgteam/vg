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
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_deconstruct(char** argv){
    cerr << "usage: " << argv[0] << " deconstruct [options] [-p|-P] <PATH> <GRAPH>" << endl
         << "Outputs VCF records for Snarls present in a graph (relative to a chosen reference path)." << endl
         << "options: " << endl
         << "    -p, --path NAME          A reference path to deconstruct against (multiple allowed)." << endl
         << "    -P, --path-prefix NAME   All paths (and/or GBWT threads) beginning with NAME used as reference (multiple allowed)." << endl
        //<< "    -A, --alt-prefix NAME    Non-reference paths (and/or GBWT threads) beginning with NAME get lumped together to same sample in VCF (multiple allowed)." << endl
         << "                             Other non-ref paths not considered as samples.  When using a GBWT, select only samples with given prefix." << endl
         << "    -H, --path-sep SEP       Obtain alt paths from the set of paths, assuming a path name hierarchy (e.g. SEP='#' and sample#phase#contig)" << endl
         << "    -r, --snarls FILE        Snarls file (from vg snarls) to avoid recomputing." << endl
         << "    -g, --gbwt FILE          only consider alt traversals that correspond to GBWT threads FILE." << endl
         << "    -T, --translation FILE   Node ID translation (as created by vg gbwt --translation) to apply to snarl names in output" << endl
         << "    -e, --path-traversals    Only consider traversals that correspond to paths in the graph." << endl
         << "    -a, --all-snarls         Process all snarls, including nested snarls (by default only top-level snarls reported)." << endl
         << "    -d, --ploidy N           Expected ploidy.  If more traversals found, they will be flagged as conflicts (default: 2)" << endl
         << "    -K, --keep-conflicted    Retain conflicted genotypes in output." << endl
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
    //vector<string> altpath_prefixes;
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
    string path_sep;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"path", required_argument, 0, 'p'},
                {"path-prefix", required_argument, 0, 'P'},
                //{"alt-prefix", required_argument, 0, 'A'},
                {"path-sep", required_argument, 0, 'H'},
                {"snarls", required_argument, 0, 'r'},
                {"gbwt", required_argument, 0, 'g'},
                {"translation", required_argument, 0, 'T'},
                {"path-traversals", no_argument, 0, 'e'},
                {"ploidy", required_argument, 0, 'd'},
                {"all-snarls", no_argument, 0, 'a'},
                {"keep-conflicted", no_argument, 0, 'K'},
                {"threads", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp:P:H:r:g:T:eKd:at:v",
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
        case 'a':
            all_snarls = true;
            break;
        case 'K':
            keep_conflicted = true;
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
    unique_ptr<PathHandleGraph> path_handle_graph;
    string path_handle_graph_filename = get_input_file_name(optind, argc, argv);
    path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(path_handle_graph_filename);

    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* graph = overlay_helper.apply(path_handle_graph.get());

    // Read the GBWT
    unique_ptr<gbwt::GBWT> gbwt_index;
    if (!gbwt_file_name.empty()) {
        gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_file_name);
        if (gbwt_index.get() == nullptr) {
            cerr << "Error [vg deconstruct]: Unable to load gbwt index file: " << gbwt_file_name << endl;
            return 1;
        }
    }

    if (!refpaths.empty()) {
        unordered_set<string> gbwt_paths;
        if (gbwt_index.get()) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                gbwt_paths.insert(thread_name(*gbwt_index, i));
            }
        }
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
                if (!Paths::is_alt(name)) {
                    refpaths.push_back(name);
                }
            });
        // Add GBWT threads if no reference paths found or we're running with -a
        if (gbwt_index.get() && (all_snarls || refpaths.empty())) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                refpaths.push_back(thread_name(*gbwt_index, i, true));
            }            
        }
    }

    // Read the translation
    unique_ptr<unordered_map<nid_t, pair<nid_t, size_t>>> translation;
    if (!translation_file_name.empty()) {
        ifstream translation_file(translation_file_name.c_str());
        if (!translation_file) {
            cerr << "Error [vg deconstruct]: Unable to load translation file: " << translation_file_name << endl;
            return 1;
        }
        translation = make_unique<unordered_map<nid_t, pair<nid_t, size_t>>>();
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
        set<string> seen_phases;
        graph->for_each_path_handle([&](const path_handle_t& path_handle) {
            string path_name = graph->get_path_name(path_handle);
            vector<string> vals = split(path_name, path_sep);
            if (vals.size() > 1) {
                auto& phase_str = vals[1];
                if (is_number(phase_str)) {
                    seen_phases.insert(phase_str);
                }
            }
        });
        map<string, int> phase_name_to_id;
        {
            int i = 0;
            for (auto& phase : seen_phases) {
                phase_name_to_id[phase] = i++;
            }
        }
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
                    if (vals.size() > 1
                        && is_number(vals[1])) {
                        phase = phase_name_to_id[vals[1]];
                    } else {
                        phase = 0;
                    }
                    alt_path_to_sample_phase[path_name] = make_pair(sample_name, phase);
                    sample_phases[sample_name].insert(phase);
                }
            });
        if (gbwt_index.get()) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                std::string path_name = thread_name(*gbwt_index, i, true);
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
        if (gbwt_index.get()) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                string sample_name = thread_sample(*gbwt_index.get(), i);
                int phase = thread_phase(*gbwt_index.get(), i);
                alt_path_to_sample_phase[sample_name] = make_pair(sample_name, phase);
            }
        }
    }

    if (refpaths.empty()) {
        cerr << "Error [vg deconstruct]: No specified reference path or prefix found in graph" << endl;
        return 1;
    }

    // Deconstruct
    Deconstructor dd;
    if (show_progress) {
        cerr << "Decsontructing top-level snarls" << endl;
    }
    dd.deconstruct(refpaths, graph, snarl_manager.get(), path_restricted_traversals, ploidy,
                   all_snarls,
                   keep_conflicted,
                   !alt_path_to_sample_phase.empty() ? &alt_path_to_sample_phase : nullptr,
                   &sample_ploidy,
                   gbwt_index.get(), translation.get());
    return 0;
}

// Register subcommand
static Subcommand vg_deconstruct("deconstruct", "create a VCF from variation in the graph", TOOLKIT, main_deconstruct);

