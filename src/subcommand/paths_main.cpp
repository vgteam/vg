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
#include "../traversal_clusters.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include <vg/io/alignment_emitter.hpp>
#include <gbwtgraph/utils.h>

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
         << "    -n, --normalize-paths    output a graph where all equivalent paths in a site a merged (using selected paths to snap to if possible)" << endl 
         << "  output path data:" << endl
         << "    -X, --extract-gam        print (as GAM alignments) the stored paths in the graph" << endl
         << "    -A, --extract-gaf        print (as GAF alignments) the stored paths in the graph" << endl
         << "    -L, --list               print (as a list of names, one per line) the path (or thread) names" << endl
         << "    -E, --lengths            print a list of path names (as with -L) but paired with their lengths" << endl
         << "    -M, --metadata           print a table of path names and their metadata" << endl
         << "    -C, --cyclicity          print a list of path names (as with -L) but paired with flag denoting the cyclicity" << endl
         << "    -F, --extract-fasta      print the paths in FASTA format" << endl
         << "    -c, --coverage           print the coverage stats for selected paths (not including cylces)" << endl
         << "  path selection:" << endl
         << "    -p, --paths-file FILE    select the paths named in a file (one per line)" << endl
         << "    -Q, --paths-by STR       select the paths with the given name prefix" << endl
         << "    -S, --sample STR         select the haplotypes or reference paths for this sample" << endl
         << "    -a, --variant-paths      select the variant paths added by 'vg construct -a'" << endl
         << "    -G, --generic-paths      select the generic, non-reference, non-haplotype paths" << endl
         << "    -R, --reference-paths    select the reference paths" << endl
         << "    -H, --haplotype-paths    select the haplotype paths paths" << endl;
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

// TODO: promote this to libhandlegraph
unordered_map<PathSense, string> SENSE_TO_STRING {
    {PathSense::REFERENCE, "REFERENCE"},
    {PathSense::GENERIC, "GENERIC"},
    {PathSense::HAPLOTYPE, "HAPLOTYPE"}
};

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
    // Starts empty, but if the options put nothing in it we will add all senses.
    unordered_set<PathSense> path_senses;
    bool list_lengths = false;
    bool list_metadata = false;
    bool list_cyclicity = false;
    size_t output_formats = 0, selection_criteria = 0;
    size_t input_formats = 0;
    bool coverage = false;
    const size_t coverage_bins = 10;
    bool normalize_paths = false;

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
            {"normalize-paths", no_argument, 0, 'n'},
            {"extract-gam", no_argument, 0, 'X'},
            {"extract-gaf", no_argument, 0, 'A'},            
            {"list", no_argument, 0, 'L'},
            {"lengths", no_argument, 0, 'E'},
            {"metadata", no_argument, 0, 'M'},
            {"cyclicity", no_argument, 0, 'C'},
            {"extract-fasta", no_argument, 0, 'F'},
            {"paths-file", required_argument, 0, 'p'},
            {"paths-by", required_argument, 0, 'Q'},
            {"sample", required_argument, 0, 'S'},
            {"variant-paths", no_argument, 0, 'a'},
            {"generic-paths", no_argument, 0, 'G'},
            {"reference-paths", no_argument, 0, 'R'},
            {"haplotype-paths", no_argument, 0, 'H'},
            {"coverage", no_argument, 0, 'c'},            

            // Hidden options for backward compatibility.
            {"threads", no_argument, 0, 'T'},
            {"threads-by", required_argument, 0, 'q'},

            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hLXv:x:g:Q:VEMCFAS:Tq:drnaGRHp:c",
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

        case 'n':
            normalize_paths = true;
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
            
        case 'M':
            list_names = true;          
            list_metadata = true;
            output_formats++;
            break;
            
        case 'C':
            list_names = true;          
            list_cyclicity = true;
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
            selection_criteria++;
            break;
                
        case 'a':
            select_alt_paths = true;
            selection_criteria++;
            break;
            
        case 'G':
            path_senses.insert(PathSense::GENERIC);
            break;

        case 'R':
            path_senses.insert(PathSense::REFERENCE);
            break;

        case 'H':
            path_senses.insert(PathSense::HAPLOTYPE);
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

    if (path_senses.empty()) {
        // No path senses were asked for explicitly.
        // Fill in default ones.
        path_senses = {
            PathSense::REFERENCE,
            PathSense::HAPLOTYPE
        };
        if (sample_name.empty()) {
            // We can support paths with no sample.
            path_senses.insert(PathSense::GENERIC);
        }
    } else {
        // We asked for path senses specifically
        selection_criteria++;
    }

    if (input_formats != 1 && input_formats != 2) {
        std::cerr << "error: [vg paths] at least one input format (-x, -g) must be specified" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (!gbwt_file.empty()) {
        bool need_graph = (extract_as_gam || extract_as_gaf || extract_as_vg || drop_paths || retain_paths || extract_as_fasta || list_lengths);
        if (need_graph && graph_file.empty()) {
            std::cerr << "error: [vg paths] a graph is needed for extracting threads in -X, -A, -V, -d, -r, -n, -E or -F format" << std::endl;
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
        std::cerr << "error: [vg paths] one output format (-X, -A, -V, -d, -r, -n, -L, -F, -E, -C or -c) must be specified" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (selection_criteria > 1) {
        std::cerr << "error: [vg paths] multiple selection criteria (-Q, -S, -a, -G/-R/-H, -p) cannot be used" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (select_alt_paths && !gbwt_file.empty()) {
        std::cerr << "error: [vg paths] selecting variant allele paths is not compatible with a GBWT index" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (list_metadata && !gbwt_file.empty()) {
        std::cerr << "error: [vg paths] listing path metadata is not compatible with a GBWT index" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if ((drop_paths || retain_paths || normalize_paths) && !gbwt_file.empty()) {
        std::cerr << "error: [vg paths] dropping, retaining or normalizing paths only works on embedded graph paths, not GBWT threads" << std::endl;
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
        path_senses = {PathSense::GENERIC};
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
            path_names.emplace(std::move(line));
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
    
    if (gbwt_index) {
        // We want to operate on a GBWT instead of the graph.

        if (!(gbwt_index->hasMetadata() && gbwt_index->metadata.hasPathNames())) {
            std::cerr << "warning: [vg paths] the GBWT index does not contain thread names" << std::endl;
            std::exit(EXIT_SUCCESS);
        }
        
        // Pre-parse some metadata
        auto gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(*gbwt_index);

        // Select the threads we are interested in.
        std::vector<gbwt::size_type> thread_ids;
        if (!sample_name.empty()) {
            thread_ids = threads_for_sample(*gbwt_index, sample_name);
        } else if(!path_prefix.empty()) {
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                PathSense sense = gbwtgraph::get_path_sense(*gbwt_index, i, gbwt_reference_samples);
                std::string name = gbwtgraph::compose_path_name(*gbwt_index, i, sense);
                if (name.length() >= path_prefix.length() && std::equal(path_prefix.begin(), path_prefix.end(), name.begin())) {
                    thread_ids.push_back(i);
                }
            }
        } else if (!path_file.empty()) {
            // TODO: there doesn't seem to be a look-up by name in the GBWT, so we check all of them
            thread_ids.reserve(path_names.size());
            for (size_t i = 0; i < gbwt_index->metadata.paths(); i++) {
                PathSense sense = gbwtgraph::get_path_sense(*gbwt_index, i, gbwt_reference_samples);
                std::string name = gbwtgraph::compose_path_name(*gbwt_index, i, sense);
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
            PathSense sense = gbwtgraph::get_path_sense(*gbwt_index, id, gbwt_reference_samples);
            std::string name = gbwtgraph::compose_path_name(*gbwt_index, id, sense);        

            // We are only interested in the name
            // TODO: do we need to consult list_cyclicity or list_metadata here?
            if (list_names && !list_lengths) {
                std::cout << name << endl;
                continue;
            }
            
            // TODO: implement list_metadata for GBWT threads?
            
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
    } else if (graph) {
        
        // Handle queries from the graph
        
        // Make a helper to loop over the selected paths in the graph
        auto for_each_selected_path = [&](const std::function<void(const path_handle_t)>& iteratee) {
            if (!path_file.empty()) {
                // We only want paths with full names from the file. Just look them all up.
                for (auto& name : path_names) {
                    if (graph->has_path(name)) {
                        auto path_handle = graph->get_path_handle(name);
                        if (path_senses.count(graph->get_sense(path_handle))) {
                            // But only take those with senses we want.
                            iteratee(path_handle);
                        }
                    }
                }
            } else {
                // We have only prefixes or other criteria.
                // We need to do a scan.
                
                // We may restrict to some exact locus names.
                unordered_set<string>* locus_name_filter = nullptr;
                // We may restrict to some exact sample names
                unordered_set<string> sample_name_set;
                unordered_set<string>* sample_name_filter = nullptr;
                if (!sample_name.empty()) {
                    // We only want paths with this exact sample name
                    sample_name_set.insert(sample_name);
                    sample_name_filter = &sample_name_set;
                }
                
                graph->for_each_path_matching(&path_senses, sample_name_filter, locus_name_filter,
                    [&](const path_handle_t& path_handle) {
                
                    // We got a path of appropriate sense, locus, and sample.
                    if (!path_prefix.empty()) {
                        // Filter by name prefix
                        std::string path_name = graph->get_path_name(path_handle);
                        
                        if (std::mismatch(path_name.begin(), path_name.end(),
                                          path_prefix.begin(), path_prefix.end()).second != path_prefix.end()) {
                            // The path does not match the prefix. Skip it.
                            return;
                        }
                    }
                    
                    // It didn't fail a prefix check, so use it.
                    iteratee(path_handle);
                });
            }
            
        };
        // Make a helper to loop over the complement set of un-selected paths in the graph
        auto for_each_unselected_path = [&](const std::function<void(const path_handle_t)>& iteratee) {
            // Get all the selected paths
            unordered_set<path_handle_t> selected;
            for_each_selected_path([&](const path_handle_t& path_handle) {
                selected.insert(path_handle);
            });
            
            unordered_set<PathSense> all_senses {
                PathSense::REFERENCE,
                PathSense::GENERIC,
                PathSense::HAPLOTYPE
            };
            
            graph->for_each_path_matching(&all_senses, nullptr, nullptr, [&](const path_handle_t& path_handle) {
                // And then, for each path of any sense
                if (!selected.count(path_handle)) {
                    // If it isn't selected, yield it.
                    iteratee(path_handle);
                }
            });
            
            
        };
        
        if (drop_paths || retain_paths || normalize_paths) {
            MutablePathMutableHandleGraph* mutable_graph = dynamic_cast<MutablePathMutableHandleGraph*>(graph.get());
            if (!mutable_graph) {
                std::cerr << "error[vg paths]: graph cannot be modified" << std::endl;
                exit(1);
            }
            SerializableHandleGraph* serializable_graph = dynamic_cast<SerializableHandleGraph*>(graph.get());
            if (!serializable_graph) {
                std::cerr << "error[vg paths]: graph cannot be saved after modification" << std::endl;
                exit(1);
            }

            vector<path_handle_t> to_destroy;
            if (drop_paths) {
                for_each_selected_path([&](const path_handle_t& path_handle) {
                    to_destroy.push_back(path_handle);
                });
            } else if (retain_paths) {
                for_each_unselected_path([&](const path_handle_t& path_handle) {
                    to_destroy.push_back(path_handle);
                });
            } else {
                assert(normalize_paths);
                unordered_set<path_handle_t> selected_paths;
                for_each_selected_path([&](const path_handle_t& path_handle) {
                    selected_paths.insert(path_handle);
                });
                merge_equivalent_traversals_in_graph(mutable_graph, selected_paths);
            }
            if (!to_destroy.empty()) {
                mutable_graph->destroy_paths(to_destroy);
            }
            
            // output the graph
            serializable_graph->serialize(cout);
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
                    vector<step_handle_t> steps;
                    for (auto& sense : path_senses) {
                        graph->for_each_step_of_sense(handle, sense, [&](const step_handle_t& step) {
                            // For every step on this handle of any sense we care about, remember it
                            steps.push_back(step);
                        });
                    }
                    unordered_set<string> unique_names;
                    unordered_set<path_handle_t> unique_paths;
                    for (auto step_handle : steps) {
                        path_handle_t step_path_handle = graph->get_path_handle_of_step(step_handle);
                        auto it = path_to_name.find(step_path_handle);
                        if (it == path_to_name.end()) {
                            string step_path_name = graph->get_path_name(step_path_handle);
                            // disregard subpath tags when counting (but not displaying)
                            it = path_to_name.insert(make_pair(step_path_handle, Paths::strip_subrange(step_path_name))).first;
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
            for_each_selected_path([&](path_handle_t path_handle) {
                string path_name = graph->get_path_name(path_handle);
                cout << path_name;
                auto& path_cov = coverage_map[path_handle];
                for (size_t cov = 0; cov < path_cov.size(); ++cov) {
                    cout << "\t" << path_cov[cov];
                }
                cout << endl;
            });
        } else {
            if (list_metadata) {
                // Add a header
                cout << "#NAME";
                if (list_lengths) {
                    cout << "\tLENGTH";
                }
                cout << "\tSENSE";
                cout << "\tSAMPLE";
                cout << "\tHAPLOTYPE";
                cout << "\tLOCUS";
                cout << "\tPHASE_BLOCK";
                cout << "\tSUBRANGE";
                if (list_cyclicity) {
                    cout << "\tCYCLICITY";
                }
                cout << endl;
            }
            for_each_selected_path([&](path_handle_t path_handle) {
                if (list_names) {
                    cout << graph->get_path_name(path_handle);
                    if (list_lengths) {
                        size_t path_length = 0;
                        for (handle_t handle : graph->scan_path(path_handle)) {
                            path_length += graph->get_length(handle);
                        }
                        cout << "\t" << path_length;
                    }
                    if (list_metadata) {
                        // Dump fields for all the metadata
                        cout << "\t" << SENSE_TO_STRING.at(graph->get_sense(path_handle));
                        auto sample = graph->get_sample_name(path_handle);
                        cout << "\t" << (sample == PathMetadata::NO_SAMPLE_NAME ? "NO_SAMPLE_NAME" : sample);
                        auto haplotype = graph->get_haplotype(path_handle);
                        cout << "\t" << (haplotype == PathMetadata::NO_HAPLOTYPE ? "NO_HAPLOTYPE" : std::to_string(haplotype));
                        auto locus = graph->get_locus_name(path_handle);
                        cout << "\t" << (locus == PathMetadata::NO_LOCUS_NAME ? "NO_LOCUS_NAME" : locus);
                        auto phase_block = graph->get_phase_block(path_handle);
                        cout << "\t" << (phase_block == PathMetadata::NO_PHASE_BLOCK ? "NO_PHASE_BLOCK" : std::to_string(phase_block));
                        auto subrange = graph->get_subrange(path_handle);
                        cout << "\t";
                        if (subrange == PathMetadata::NO_SUBRANGE) {
                            cout << "NO_SUBRANGE";
                        } else if (subrange.second == PathMetadata::NO_END_POSITION) {
                            cout << subrange.first;
                        } else {
                            cout << subrange.first << "-" << subrange.second;
                        }
                    }
                    if (list_cyclicity) {
                        bool directed_cyclic = false; // same node visited twice in same orientation
                        bool undirected_cyclic = false; // same not visited twice in any orientation
                        unordered_set<handle_t> visits;
                        graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                                handle_t handle = graph->get_handle_of_step(step_handle);
                                if (visits.count(handle)) {
                                    directed_cyclic = true;
                                    undirected_cyclic = false;
                                } else if (visits.count(graph->flip(handle))) {
                                    undirected_cyclic = true;
                                }
                                visits.insert(handle);
                                return !directed_cyclic || !undirected_cyclic;
                            });
                        cout << "\t" << (directed_cyclic ? "directed-cyclic" : "directed-acyclic")
                             << "\t" << (undirected_cyclic ? "undirected-cyclic" : "undirected-acyclic");
                    }
                    cout << endl;
                } else {
                    Path path = path_from_path_handle(*graph, path_handle);
                    if (extract_as_gam || extract_as_gaf) {
                        aln_emitter->emit_singles({alignment_from_path(*graph, path)});
                    } else if (extract_as_vg) {
                        chunk_to_emitter(path, *graph_emitter); 
                    } else if (extract_as_fasta) {
                        write_fasta_sequence(graph->get_path_name(path_handle), path_sequence(*graph, path), cout);
                    }
                }
            });
        }
    }
    return 0;

}

// Register subcommand
static Subcommand vg_paths("paths", "traverse paths in the graph", main_paths);

