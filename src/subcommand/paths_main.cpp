/** \file paths_main.cpp
 *
 * Defines the "vg paths" subcommand, which reads paths in the graph.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <algorithm>
#include <iostream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include "../gbwt_helper.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <vg/io/alignment_emitter.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_paths(char** argv) {
    cerr << "usage: " << argv[0] << " paths [options]" << endl
         << "options:" << endl
         << "  input:" << endl
         << "    -x, --xg FILE            use the paths and haplotypes in this graph FILE. Supports GBZ haplotypes." <<endl
         << "                             (Also accepts -v, --vg)" << endl
         << "    -g, --gbwt FILE          use the threads in the GBWT index in FILE" << endl
         << "                             (graph also required for most output options; -g takes priority over -x)" << endl
         << "  output graph (.vg format)" << endl
         << "    -V, --extract-vg         output a path-only graph covering the selected paths" << endl
         << "    -d, --drop-paths         output a graph with the selected paths removed" << endl
         << "    -r, --retain-paths       output a graph with only the selected paths retained" << endl
         << "  output path data:" << endl
         << "    -X, --extract-gam        print (as GAM alignments) the stored paths in the graph" << endl
         << "    -A, --extract-gaf        print (as GAF alignments) the stored paths in the graph" << endl
         << "    -L, --list               print (as a list of names, one per line) the path (or thread) names" << endl
         << "    -E, --lengths            print a list of path names (as with -L) but paired with their lengths" << endl
         << "    -C, --cyclicity          print a list of path names (as with -L) but paired with flag denoting the cyclicity" << endl
         << "    -F, --extract-fasta      print the paths in FASTA format" << endl
         << "    -c, --coverage           print the coverage stats for selected paths (not including cylces)" << endl
         << "  path selection:" << endl
         << "    -p, --paths-file FILE    select the paths named in a file (one per line)" << endl
         << "    -Q, --paths-by STR       select the paths with the given name prefix" << endl
         << "    -S, --sample STR         select the haplotypes or reference paths for this sample" << endl
         << "    -a, --variant-paths      select the variant paths added by 'vg construct -a'" << endl;
}

/// Chunk a path and emit it in Graph messages.
/// Paht must have ranks set.
void chunk_to_emitter(const Path& path, vg::io::ProtobufEmitter<Graph>& graph_emitter) {
    size_t chunk_size = 10000;
                    
    for (size_t start = 0; start < path.mapping_size(); start += chunk_size) {
        // Make sure to chunk.
        // TODO: Can we avoild a copy here somehow?
        Path chunk;
        chunk.set_name(path.name());
        chunk.set_is_circular(path.is_circular());
        
        for (size_t i = 0; i < chunk_size && start + i < path.mapping_size(); i++) {
            // Copy over this batch of mappings
            *chunk.add_mapping() = path.mapping(start + i);
        }
        
        // Emit a graph chunk containing htis part of the path
        Graph g;
        *(g.add_path()) = std::move(chunk);
        graph_emitter.write(std::move(g));
    }
}

int main_paths(int argc, char** argv) {

    if (argc == 2) {
        help_paths(argv);
        return 1;
    }

    bool extract_as_gam = false;
    bool extract_as_gaf = false;
    bool extract_as_vg = false;
    bool list_names = false;
    bool extract_as_fasta = false;
    bool drop_paths = false;
    bool retain_paths = false;
    string graph_file;
    string gbwt_file;
    string path_prefix;
    string sample_name;
    string path_file;
    bool select_alt_paths = false;
    // What kinds of paths are we interested in?
    vector<handlegraph::PathMetadata::Sense> path_senses {handlegraph::PathMetadata::SENSE_REFERENCE,
                                                          handlegraph::PathMetadata::SENSE_GENERIC,
                                                          handlegraph::PathMetadata::SENSE_HAPLOTYPE}
    bool list_lengths = false;
    size_t output_formats = 0, selection_criteria = 0;
    size_t input_formats = 0;
    bool list_cyclicity = false;
    bool coverage = false;
    const size_t coverage_bins = 10;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =

        {
            {"vg", required_argument, 0, 'v'},
            {"xg", required_argument, 0, 'x'},
            {"gbwt", required_argument, 0, 'g'},
            {"extract-vg", no_argument, 0, 'V'},
            {"drop-paths", no_argument, 0, 'd'},
            {"retain-paths", no_argument, 0, 'r'},
            {"extract-gam", no_argument, 0, 'X'},
            {"extract-gaf", no_argument, 0, 'A'},            
            {"list", no_argument, 0, 'L'},
            {"lengths", no_argument, 0, 'E'},
            {"extract-fasta", no_argument, 0, 'F'},
            {"paths-file", required_argument, 0, 'p'},
            {"paths-by", required_argument, 0, 'Q'},
            {"sample", required_argument, 0, 'S'},
            {"variant-paths", no_argument, 0, 'a'},
            {"cyclicity", no_argument, 0, 'C'},
            {"coverage", no_argument, 0, 'c'},            

            // Hidden options for backward compatibility.
            {"threads", no_argument, 0, 'T'},
            {"threads-by", required_argument, 0, 'q'},

            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hLXv:x:g:Q:VEFAS:Tq:drap:Cc",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'v': // Fall through
        case 'x':
            graph_file = optarg;
            ++input_formats;
            break;

        case 'g':
            gbwt_file = optarg;
            ++input_formats;
            break;

        case 'V':
            extract_as_vg = true;
            output_formats++;
            break;
                
        case 'd':
            drop_paths = true;
            output_formats++;
            break;
            
        case 'r':
            retain_paths = true;
            output_formats++;
            break;
                
        case 'X':
            extract_as_gam = true;
            output_formats++;
            break;

        case 'A':
            extract_as_gaf = true;
            output_formats++;
            break;

        case 'L':
            list_names = true;
            output_formats++;
            break;

        case 'E':
            list_names = true;
            list_lengths = true;
            output_formats++;
            break;
                
        case 'F':
            extract_as_fasta = true;
            output_formats++;
            break;
                
        case 'p':
            path_file = optarg;
            selection_criteria++;
            break;

        case 'Q':
            path_prefix = optarg;
            selection_criteria++;
            break;

        case 'S':
            sample_name = optarg;
            // We only care about things with references now.
            path_senses = {handlegraph::PathMetadata::SENSE_REFERENCE, handlegraph::PathMetadata::SENSE_HAPLOTYPE};
            selection_criteria++;
            break;
                
        case 'a':
            select_alt_paths = true;
            selection_criteria++;
            break;

        case 'C':
            list_names = true;          
            list_cyclicity = true;
            output_formats++;
            break;

        case 'c':
            coverage = true;
            output_formats++;
            break;

        case 'T':
            std::cerr << "warning: [vg paths] option --threads is obsolete and unnecessary" << std::endl;
            break;

        case 'q':
            std::cerr << "warning: [vg paths] option --threads-by is deprecated; please use --paths-by" << std::endl;
            path_prefix = optarg;
            selection_criteria++;
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

    if (input_formats != 1 && input_formats != 2) {
        std::cerr << "error: [vg paths] at least one input format (-v, -x, -g) must be specified" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (!gbwt_file.empty()) {
        bool need_graph = (extract_as_gam || extract_as_gaf || extract_as_vg || drop_paths || retain_paths || extract_as_fasta || list_lengths);
        if (need_graph && graph_file.empty()) {
            std::cerr << "error: [vg paths] a graph is needed for extracting threads in -X, -A, -V, -d, -r, -E or -F format" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (!need_graph && !graph_file.empty()) {
            // TODO: This should be an error, but we display a warning instead for backward compatibility.
            //std::cerr << "error: [vg paths] cannot read input from multiple sources" << std::endl;
            //std::exit(EXIT_FAILURE);
            std::cerr << "warning: [vg paths] graph unnecessary for listing GBWT threads" << std::endl;
        }
    } 
    if (output_formats != 1) {
        std::cerr << "error: [vg paths] one output format (-X, -A, -V, -d, -r, -L, -F, -E, -C or -c) must be specified" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (selection_criteria > 1) {
        std::cerr << "error: [vg paths] multiple selection criteria (-Q, -S, -a, -p) cannot be used" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (select_alt_paths && !gbwt_file.empty()) {
        std::cerr << "error: [vg paths] selecting variant allele paths is not compatible with a GBWT index" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if ((drop_paths || retain_paths) && !gbwt_file.empty()) {
        std::cerr << "error: [vg paths] dropping or retaining paths only works on embedded graph paths, not GBWT threads" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (coverage && !gbwt_file.empty()) {
        std::cerr << "error: [vg paths] coverage option -c only works on embedded graph paths, not GBWT threads" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    
    if (select_alt_paths) {
        // alt paths all have a specific prefix
        path_prefix = "_alt_";
        // And are all generic sense.
        path_senses = {handlegraph::PathMetadata::SENSE_GENERIC};
    }
    
    // Load whatever indexes we were given
    // Note: during handlifiction, distinction between -v and -x options disappeared.
    unique_ptr<PathHandleGraph> graph;
    if (!graph_file.empty()) {
        // Load the graph
        graph = vg::io::VPKG::load_one<PathHandleGraph>(graph_file);
    }
    unique_ptr<gbwt::GBWT> gbwt_index;
    if (!gbwt_file.empty()) {
        // We want a gbwt
        
        // Load the GBWT from its container
        gbwt_index = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_file);

        if (gbwt_index.get() == nullptr) {
          // Complain if we couldn't.
          cerr << "error: [vg paths] unable to load gbwt index file" << endl;
          exit(1);
        }
    }
    
    set<string> path_names;
    if (!path_file.empty()) {
        ifstream path_stream(path_file);
        if (!path_stream) {
            cerr << "error: cannot open path name file " << path_file << endl;
            exit(EXIT_FAILURE);
        }
        
        string line;
        while (getline(path_stream, line)) {
            path_names.emplace(move(line));
        }
    }

    // We may need to emit a stream of Alignments
    unique_ptr<vg::io::AlignmentEmitter> aln_emitter; 

    // Or we might need to emit a stream of VG Graph objects
    unique_ptr<vg::io::ProtobufEmitter<Graph>> graph_emitter;
    if (extract_as_gam || extract_as_gaf) {
        // Open up a GAM/GAF output stream
        aln_emitter = vg::io::get_non_hts_alignment_emitter("-", extract_as_gaf ? "GAF" : "GAM", {}, get_thread_count(),
                                                            graph.get());
    } else if (extract_as_vg) {
        // Open up a VG Graph chunk output stream
        graph_emitter = unique_ptr<vg::io::ProtobufEmitter<Graph>>(new vg::io::ProtobufEmitter<Graph>(cout));
    }
    
    if (gbwt_index.get() != nullptr) {

        if (!(gbwt_index->hasMetadata() && gbwt_index->metadata.hasPathNames())) {
            std::cerr << "warning: [vg paths] the GBWT index does not contain thread names" << std::endl;
            std::exit(EXIT_SUCCESS);
        }

        // Select the threads we are interested in.
        std::vector<gbwt::size_type> thread_ids;
        if (!sample_name.empty()) {
            thread_ids = threads_for_sample(*gbwt_index, sample_name);
        } else if(!path_prefix.empty()) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                std::string name = thread_name(*gbwt_index, i);
                if (name.length() >= path_prefix.length() && std::equal(path_prefix.begin(), path_prefix.end(), name.begin())) {
                    thread_ids.push_back(i);
                }
            }
        } else if (!path_file.empty()) {
            // TODO: there doesn't seem to be a look-up by name in the GBWT, so we check all of them
            thread_ids.reserve(path_names.size());
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                std::string name = thread_name(*gbwt_index, i);
                if (path_names.count(name)) {
                    thread_ids.push_back(i);
                }
            }
            if (thread_ids.size() != path_names.size()) {
                std::cerr << "error: [vg paths] could not find all path names from file in GBWT index" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else {
            thread_ids.reserve(gbwt_index->metadata.paths());
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                thread_ids.push_back(i);
            }
        }
        
        // Process the threads.
        for (gbwt::size_type id : thread_ids) {
            std::string name = thread_name(*gbwt_index, id);        

            // We are only interested in the name
            if (list_names && !list_lengths) {
                std::cout << name << endl;
                continue;
            }
            
            // Otherwise we need the actual thread data
            Path path = extract_gbwt_path(*graph, *gbwt_index, id);
            if (extract_as_gam || extract_as_gaf) {
                // Write as an Alignment. Must contain the whole path.
                aln_emitter->emit_singles({alignment_from_path(*graph, path)});
            } else if (extract_as_vg) {
                // Write as a Path in a VG
                chunk_to_emitter(path, *graph_emitter);
            } else if (extract_as_fasta) {
                write_fasta_sequence(name, path_sequence(*graph, path), cout);
            }
            if (list_lengths) {
                cout << path.name() << "\t" << path_to_length(path) << endl;
            }
            if (list_cyclicity) {
                bool cyclic = false;
                unordered_set<pair<nid_t, bool>> visits;
                for (size_t i = 0; i < path.mapping_size() && !cyclic; ++i) {
                    const Mapping& mapping = path.mapping(i);
                    pair<unordered_set<pair<nid_t, bool>>::iterator, bool> ret =
                        visits.insert(make_pair(mapping.position().node_id(), mapping.position().is_reverse()));
                    if (ret.second == false) {
                        cyclic = true;
                    }
                }
                cout << path.name() << "\t" << (cyclic ? "cyclic" : "acyclic") << endl;
            }
        }
    } else if (graph.get() != nullptr) {
        
        // Handle non-thread queries from vg or xg
        
        auto check_path_name = [&](const string& name) {
            // Note: we're using a filter rather than index, so O(#paths in graph).
            if (path_file.empty()) {
                return (path_prefix.empty() ||
                        std::mismatch(name.begin(), name.end(),
                                      path_prefix.begin(), path_prefix.end()).second == path_prefix.end());
            }
            else {
                return bool(path_names.count(name));
            }
        };

        if (drop_paths || retain_paths) {

            vector<string> to_destroy;
            graph->for_each_path_handle([&](const path_handle_t& path_handle) {
                    string name = graph->get_path_name(path_handle);
                    if (check_path_name(name) == drop_paths) {
                        to_destroy.push_back(name);
                }
            });
            for (string& path_name : to_destroy) {
                dynamic_cast<MutablePathMutableHandleGraph*>(graph.get())->destroy_path(graph->get_path_handle(path_name));
            }
            
            // output the graph
            dynamic_cast<SerializableHandleGraph*>(graph.get())->serialize(cout);
        }
        else if (coverage) {
            // for every node, count the number of unique paths.  then add the coverage count to each one
            // (we're doing the whole graph here, which could be inefficient in the case the user is selecting
            //  a small path)
            unordered_map<path_handle_t, vector<int64_t>> coverage_map;
            // big speedup
            unordered_map<path_handle_t, string> path_to_name;
            size_t max_coverage = 0;
            graph->for_each_handle([&](handle_t handle) {
                    vector<step_handle_t> steps = graph->steps_of_handle(handle);
                    unordered_set<string> unique_names;
                    unordered_set<path_handle_t> unique_paths;
                    for (auto step_handle : steps) {
                        path_handle_t step_path_handle = graph->get_path_handle_of_step(step_handle);
                        auto it = path_to_name.find(step_path_handle);
                        if (it == path_to_name.end()) {
                            string step_path_name = graph->get_path_name(step_path_handle);
                            // disregard subpath tags when counting (but not displaying)
                            auto subpath = Paths::parse_subpath_name(step_path_name);
                            string& parsed_name = get<0>(subpath) ? get<1>(subpath) : step_path_name;
                            it = path_to_name.insert(make_pair(step_path_handle, parsed_name)).first;
                        }
                        unique_names.insert(it->second);
                        unique_paths.insert(step_path_handle);
                    }
                    for (auto path : unique_paths) {
                        vector<int64_t>& cov = coverage_map[path];
                        if (cov.size() < unique_paths.size()) {
                            cov.resize(unique_paths.size(), 0);
                        }
                        cov[unique_paths.size() - 1] += graph->get_length(graph->get_handle_of_step(steps[0]));
                        max_coverage = std::max(max_coverage, unique_names.size() - 1);
                    }
                });
            // figure out the bin size
            int64_t bin_size = 1;
            if (max_coverage > coverage_bins) {
                // reserve the first 2 bins for coverage = 0 and 1 no matter 1
                bin_size = (max_coverage - 2) / (coverage_bins - 2);
                if ((max_coverage - 2) % (coverage_bins - 2)) {
                    ++bin_size;
                }
            }
            // compute cumulative coverage
            for (auto& path_cov : coverage_map) {
                int64_t cum_cov = 0;
                vector<int64_t>& cov = path_cov.second;
                cov.resize(max_coverage + 1, 0);
                // bin it up if necessary
                if (cov.size() > coverage_bins) {
                    vector<int64_t> binned_cov(coverage_bins, 0);
                    // reserve the first 2 bins for coverage = 0 and 1 no matter 1
                    binned_cov[0] = cov[0];
                    binned_cov[1] = cov[1];
                    // remaining bins
                    for (size_t bin = 0; bin < coverage_bins - 2; ++bin) {
                        for (size_t i = 0; i < bin_size && (2 + bin * bin_size + i < cov.size()); ++i) {
                            binned_cov[2 + bin] += cov[2 + bin * bin_size + i];
                        }
                    }
                    swap(cov, binned_cov);
                }
                // accumulate
                for (auto cov_it = path_cov.second.rbegin(); cov_it != path_cov.second.rend(); ++cov_it) {
                    cum_cov += *cov_it;
                    *cov_it = cum_cov;
                }
            }
            cout << "PathName";
            for (size_t cov = 0; cov <= min(max_coverage, coverage_bins); ++cov) {
                cout << "\t";
                if (cov < 2 || bin_size == 1) {
                    cout << cov << "-" << cov;
                } else {
                    cout << (2 + (cov - 2) * bin_size) << "-" << (2 + (cov - 2) * bin_size + bin_size - 1);
                } 
            }
            cout << endl;
            graph->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = graph->get_path_name(path_handle);
                    if (check_path_name(path_name)) {
                        cout << path_name;
                        auto& path_cov = coverage_map[path_handle];
                        for (size_t cov = 0; cov < path_cov.size(); ++cov) {
                            cout << "\t" << path_cov[cov];
                        }
                        cout << endl;
                    }
                });
        } else {            
            graph->for_each_path_handle([&](path_handle_t path_handle) {
                string path_name = graph->get_path_name(path_handle);
                if (check_path_name(path_name)) {
                    if (list_names) {
                        cout << path_name;
                        if (list_lengths) {
                            size_t path_length = 0;
                            for (handle_t handle : graph->scan_path(path_handle)) {
                                path_length += graph->get_length(handle);
                            }
                            cout << "\t" << path_length;
                        }
                        if (list_cyclicity) {
                            bool cyclic = false;
                            unordered_set<handle_t> visits;
                            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                                    pair<unordered_set<handle_t>::iterator, bool> ret = visits.insert(graph->get_handle_of_step(step_handle));
                                    if (ret.second == false) {
                                        cyclic = true;
                                    }
                                    return !cyclic;
                                });
                            cout << "\t" << (cyclic ? "cyclic" : "acyclic");
                        }
                        cout << endl;
                    } else {
                        Path path = path_from_path_handle(*graph, path_handle);
                        if (extract_as_gam || extract_as_gaf) {
                            aln_emitter->emit_singles({alignment_from_path(*graph, path)});
                        } else if (extract_as_vg) {
                            chunk_to_emitter(path, *graph_emitter); 
                        } else if (extract_as_fasta) {
                            write_fasta_sequence(path_name, path_sequence(*graph, path), cout);
                        }
                    }
                }
            });
        }
    }
    return 0;

}

// Register subcommand
static Subcommand vg_paths("paths", "traverse paths in the graph", main_paths);

