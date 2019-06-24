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
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_deconstruct(char** argv){
    cerr << "usage: " << argv[0] << " deconstruct [options] [-p|-P] <PATH> <my_graph>.vg" << endl
         << "Outputs VCF records for Snarls present in a graph (relative to a chosen reference path)." << endl
         << "options: " << endl
         << "    -p, --path NAME        A reference path to deconstruct against (comma-separated list accepted)." << endl
         << "    -P, --path-prefix NAME All paths beginning with NAME used as reference (comma-separated list accepted)." << endl
         << "    -A, --alt-prefix NAME  Non-reference paths beginning with NAME get lumped together to same sample in VCF (comma-separated list accepted)." << endl
         << "    -r, --snarls FILE      Snarls file (from vg snarls) to avoid recomputing." << endl
         << "    -e, --path-traversals  Only consider traversals that correspond to paths in the grpah." << endl
         << "    -t, --threads N        Use N threads" << endl
         << "    -v, --verbose          Print some status messages" << endl
         << endl;
}

int main_deconstruct(int argc, char** argv){
    if (argc <= 2) {
        help_deconstruct(argv);
        return 1;
    }

    vector<string> refpaths;
    vector<string> refpath_prefixes;
    vector<string> altpath_prefixes;
    string graphname;
    string snarl_file_name;
    bool path_restricted_traversals = false;
    bool show_progress = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"path", required_argument, 0, 'p'},
                {"path-prefix", required_argument, 0, 'P'},
                {"alt-prefix", required_argument, 0, 'A'},
                {"snarls", required_argument, 0, 'r'},
                {"path-traversals", no_argument, 0, 'e'},
                {"threads", required_argument, 0, 't'},
                {"verbose", no_argument, 0, 'v'},
                {0, 0, 0, 0}

            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hp:P:A:r:et:v",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'p':
            refpaths = split(optarg, ",");
            break;
        case 'P':
            refpath_prefixes = split(optarg, ",");
            break;
        case 'A':
            altpath_prefixes = split(optarg, ",");
            break;            
        case 'r':
            snarl_file_name = optarg;
            break;
        case 'e':
            path_restricted_traversals = true;
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

    if (refpaths.empty() && refpath_prefixes.empty()) {
        cerr << "Error [vg deconstruct]: Reference path(s) and/or prefix(es) must be given with -p and/or -P" << endl;
        return 1;
    }

    if (!altpath_prefixes.empty() && !path_restricted_traversals) {
        cerr << "Error [vg decontruct]: -A can only be used with -e" << endl;
    }
    
    string graph_file_name = get_input_file_name(optind, argc, argv);

    vg::VG* graph;
    get_input_file(graph_file_name, [&](istream& in) {
            if (show_progress) {
                cerr << "Loading graph" << endl;
            }
            graph = new VG(in);
    });

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
        CactusSnarlFinder finder(*graph);
        if (show_progress) {
            cerr << "Finding snarls" << endl;
        }
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls())));
    }

    // We use this to map, for example, from chromosome to genome (eg S288C.chrXVI --> S288C)
    unordered_map<string, string> alt_path_to_prefix;
    
    // process the prefixes
    if (!refpath_prefixes.empty()) {
        graph->for_each_path_handle([&](const path_handle_t& path_handle) {
                string path_name = graph->get_path_name(path_handle);
                bool is_ref = false;
                for (auto& prefix : refpath_prefixes) {
                    if (path_name.compare(0, prefix.size(), prefix) == 0) {
                        refpaths.push_back(path_name);
                        is_ref = true;
                        break;
                    }
                }
                if (!is_ref) {
                    for (auto& prefix : altpath_prefixes) {
                        if (path_name.compare(0, prefix.size(), prefix) == 0) {
                            alt_path_to_prefix[path_name] = prefix;
                        }
                    }
                }
            });
    }

    // make sure we have at least one reference
    bool found_refpath = false;
    for (size_t i = 0; i < refpaths.size() && !found_refpath; ++i) {
        found_refpath = found_refpath || graph->has_path(refpaths[i]);
    }

    if (!found_refpath) {
        cerr << "Error [vg deconstruct]: No specified reference path or prefix found in graph" << endl;
        return 1;
    }

    // Deconstruct
    Deconstructor dd;
    if (show_progress) {
        cerr << "Decsontructing top-level snarls" << endl;
    }
    dd.deconstruct(refpaths, graph, snarl_manager.get(), path_restricted_traversals,
                   !alt_path_to_prefix.empty() ? &alt_path_to_prefix : nullptr);
    return 0;
}

// Register subcommand
static Subcommand vg_deconstruct("deconstruct", "convert a graph into VCF relative to a reference", main_deconstruct);

