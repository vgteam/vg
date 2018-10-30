/** \file augment_main.cpp
 *
 * Defines the "vg augment" subcommand, which augments a graph using a GAM as a first
 * step before computing sites and calling variants.
 *
 * It presently only supports augmentation via pileup (as originally used by vg call).  
 *
 * vg mod -i is an alternative, but as currently implemented it is too expensive.   Ideally,
 * a more heuristic/fast implementation would get added as an option here.
 *
 * The new vg count structure, provided edge/edit/strand/quality is supported, may be an
 * interesting replacement for the current protobuf pileups which ought to go, one way or the other.
 *
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../option.hpp"

#include "../vg.hpp"
#include "../pileup_augmenter.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

// this used to be pileup_main()
static Pileups* compute_pileups(VG* graph, const string& gam_file_name, int thread_count, int min_quality,
                                int max_mismatches, int window_size, int max_depth, bool use_mapq,
                                bool show_progress);

// this used to be the first half of call_main()
static void augment_with_pileups(PileupAugmenter& augmenter, Pileups& pileups, bool expect_subgraph,
                                 bool show_progress);

void help_augment(char** argv, ConfigurableParser& parser) {
    cerr << "usage: " << argv[0] << " augment [options] <graph.vg> <alignment.gam> > augmented_graph.vg" << endl
         << "Embed GAM alignments into a graph to facilitate variant calling" << endl
         << endl
         << "general options:" << endl
         << "    -a, --augmentation-mode M   augmentation mode.  M = {pileup, direct} [direct]" << endl
         << "    -Z, --translation FILE      save translations from augmented back to base graph to FILE" << endl
         << "    -A, --alignment-out FILE    save augmented GAM reads to FILE" << endl
         << "    -h, --help                  print this help message" << endl
         << "    -p, --progress              show progress" << endl
         << "    -v, --verbose               print information and warnings about vcf generation" << endl
         << "    -t, --threads N             number of threads to use" << endl
         << "pileup options:" << endl
         << "    -P, --pileup FILE           save pileups to FILE" << endl
         << "    -S, --support FILE          save supports to FILE" << endl                
         << "    -g, --min-aug-support N     minimum support to augment graph ["
         << PileupAugmenter::Default_min_aug_support << "]" << endl
         << "    -U, --subgraph              expect a subgraph and ignore extra pileup entries outside it" << endl
         << "    -q, --min-quality N         ignore bases with PHRED quality < N (default=10)" << endl
         << "    -m, --max-mismatches N      ignore bases with > N mismatches within window centered on read (default=1)" << endl
         << "    -w, --window-size N         size of window to apply -m option (default=0)" << endl
         << "    -M, --ignore-mapq           do not combine mapping qualities with base qualities in pileup" << endl;
    
     // Then report more options
     parser.print_help(cerr);
}

int main_augment(int argc, char** argv) {

    // augmentation mode
    string augmentation_mode = "direct";
    
    // load pileupes from here
    string pileup_file_name;

    // minimum support to consider adding a variant to the graph
    int min_aug_support = PileupAugmenter::Default_min_aug_support;
        
    // Should we expect a subgraph and ignore pileups for missing nodes/edges?
    bool expect_subgraph = false;

    // Write the translations (as protobuf) to this path
    string translation_file_name;

    // Write the supports (as protobuf) to this path
    string support_file_name;
    
    // Load in GAM alignments to map over to the augmented graph from here
    string gam_in_file_name;

    // Write the GAM alignments (from gam_in_file_name) projected on the augmented graph here
    string gam_out_file_name;

    // Print some progress messages to screen
    bool show_progress = false;

    // Print verbose message
    bool verbose = false;

    // Number of threads to use (will default to all if not specified)
    int thread_count = 0;

    // Bases wit quality less than 10 will not be added to the pileup
    int min_quality = 10;

    // Bases with more than this many mismatches within the window_size not added
    int max_mismatches = 1;

    // Window size for above (0 effectively turns this check off)
    int window_size = 0;

    // Hack to prevent protobuf messages from getting too big by limiting depth at
    // any given position to max_depth
    int max_depth = 1000;
    
    // Combine MAPQ and PHRED base qualities to determine quality at each position
    // If false, only PHRED base quality will be used. 
    bool use_mapq = true;


    static const struct option long_options[] = {
        // General Options
        {"augmentation-mode", required_argument, 0, 'a'},
        {"translation", required_argument, 0, 'Z'},
        {"alignment-out", required_argument, 0, 'A'}, 
        {"help", no_argument, 0, 'h'},
        {"progress", required_argument, 0, 'p'},
        {"verbose", no_argument, 0, 'v'},
        {"threads", required_argument, 0, 't'},
        // Pileup Options
        {"pileup", required_argument, 0, 'P'},        
        {"support", required_argument, 0, 'S'},
        {"min-quality", required_argument, 0, 'q'},
        {"max-mismatches", required_argument, 0, 'm'},
        {"window-size", required_argument, 0, 'w'},
        {"ignore-mapq", no_argument, 0, 'M'},
        {"min-aug-support", required_argument, 0, 'g'},
        {"subgraph", no_argument, 0, 'U'},
        {0, 0, 0, 0}
    };
    static const char* short_options = "a:Z:A:hpvt:P:S:q:m:w:Mg:U";
    optind = 2; // force optind past command positional arguments

    // This is our command-line parser
    ConfigurableParser parser(short_options, long_options, [&](int c) {
        // Parse all the options we have defined here.
        switch (c)
        {
            // General Options
        case 'a':
            augmentation_mode = optarg;
            break;
        case 'Z':
            translation_file_name = optarg;
            break;
        case 'A':
            gam_out_file_name = optarg;
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_augment(argv, parser);
            exit(1);
            break;
        case 'p':
            show_progress = true;
            break;
        case 'v':
            verbose = true;
            break;            
        case 't':
            thread_count = parse<int>(optarg);
            break;

            // Pileup Options
        case 'P':
            pileup_file_name = optarg;
            break;
        case 'S':
            support_file_name = optarg;
            break;            
        case 'q':
            min_quality = parse<int>(optarg);
            break;
        case 'm':
            max_mismatches = parse<int>(optarg);
            break;
        case 'w':
            window_size = parse<int>(optarg);
            break;
        case 'M':
            use_mapq = false;
            break;            
        case 'g':
            min_aug_support = parse<int>(optarg);
            break;            
        case 'U':
            expect_subgraph = true;
            break;
            
        default:
          abort ();
        }
    });

    // Parse the command line options, updating optind.
    parser.parse(argc, argv);

    if (thread_count != 0) {
        // Use a non-default number of threads
        omp_set_num_threads(thread_count);
    }
    thread_count = get_thread_count();

    // Parse the two positional arguments
    if (optind + 2 > argc) {
        cerr << "[vg augment] error: too few arguments" << endl;
        help_augment(argv, parser);
        return 1;
    }

    string graph_file_name = get_input_file_name(optind, argc, argv);
    gam_in_file_name = get_input_file_name(optind, argc, argv);

    if (gam_in_file_name == "-" && graph_file_name == "-") {
        cerr << "[vg augment] error: graph and gam can't both be from stdin." << endl;
        return 1;
    }
    if (gam_in_file_name == "-" && !gam_out_file_name.empty()) {
        cerr << "[vg augment] error: cannot stream input gam when using -A option (as it requires 2 passes)" << endl;
        return 1;
    }

    if (augmentation_mode != "pileup" && augmentation_mode != "direct") {
        cerr << "[vg augment] error: pileup and direct are currently the only supported augmentation modes (-a)" << endl;
        return 1;
    }

    if (augmentation_mode != "direct" and !gam_out_file_name.empty()) {
        cerr << "[vg augment] error: GAM output only works with \"direct\" augmentation mode" << endl;
        return 1;
    }

    if (augmentation_mode != "pileup" and (!support_file_name.empty() || !pileup_file_name.empty())) {
        cerr << "[vg augment] error: Pileup (-P) and Support (-S) output only work with  \"pileup\" augmentation mode" << endl;
        return 1;
    }
    
    

    // read the graph
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }
    VG* graph;
    get_input_file(graph_file_name, [&](istream& in) {
        graph = new VG(in);
    });
    
    
    Pileups* pileups = nullptr;
    
    if (!pileup_file_name.empty() || augmentation_mode == "pileup") {
        // We will need the computed pileups
        
        // compute the pileups from the graph and gam
        pileups = compute_pileups(graph, gam_in_file_name, thread_count, min_quality, max_mismatches,
                                  window_size, max_depth, use_mapq, show_progress);
    }
        
    if (!pileup_file_name.empty()) {
        // We want to write out pileups.
        if (show_progress) {
            cerr << "Writing pileups" << endl;
        }
        ofstream pileup_file(pileup_file_name);
        if (!pileup_file) {
            cerr << "[vg augment] error: unable to open output pileup file: " << pileup_file_name << endl;
            exit(1);
        }
        pileups->write(pileup_file);
    }

    if (augmentation_mode == "direct") {
        // Augment with the reads
        
        if (!support_file_name.empty()) {
            cerr << "[vg augment] error: support calculation in direct augmentation mode is unimplemented" << endl;
            exit(1);
        }
        
        // We don't need any pileups
        if (pileups != nullptr) {
            delete pileups;
            pileups = nullptr;
        }
    
        // Load all the reads
        vector<Alignment> reads;
        // And pull out their paths
        vector<Path> read_paths;
        get_input_file(gam_in_file_name, [&](istream& alignment_stream) {
            stream::for_each<Alignment>(alignment_stream, [&](Alignment& alignment) {
                // Trim the softclips off of every read
                // Work out were to cut
                int cut_start = softclip_start(alignment);
                int cut_end = softclip_end(alignment);
                // Cut the sequence and quality
                alignment.set_sequence(alignment.sequence().substr(cut_start, alignment.sequence().size() - cut_start - cut_end));
                if (alignment.quality().size() != 0) {
                    alignment.set_quality(alignment.quality().substr(cut_start, alignment.quality().size() - cut_start - cut_end));
                }
                // Trim the path
                *alignment.mutable_path() = trim_hanging_ends(alignment.path());
                
                // Save every read
                reads.push_back(alignment);
                // And the path for the read, separately
                // TODO: Make edit use callbacks or something so it doesn't need a vector of paths necessarily
                read_paths.push_back(alignment.path());
            });
        });
        
        // Augment the graph, rewriting the paths.
        // Don't embed paths or break at ends.
        auto translation = graph->edit(read_paths, false, true, false);
        
        // Write the augmented graph
        if (show_progress) {
            cerr << "Writing augmented graph" << endl;
        }
        graph->serialize_to_ostream(cout);
        
        if (!translation_file_name.empty()) {
            // Write the translations
            if (show_progress) {
                cerr << "Writing translation table" << endl;
            }
            ofstream translation_file(translation_file_name);
            if (!translation_file) {
                cerr << "[vg augment]: Error opening translation file: " << translation_file_name << endl;
                return 1;
            }
            stream::write_buffered(translation_file, translation, 0);
            translation_file.close();
        }        
        
        if (!gam_out_file_name.empty()) {
            // Write out the modified GAM
            
            ofstream gam_out_file(gam_out_file_name);
            if (!gam_out_file) {
                cerr << "[vg augment]: Error opening output GAM file: " << gam_out_file_name << endl;
                return 1;
            }
            
            // We use this buffer and do a buffered write
            vector<Alignment> gam_buffer;
            for (size_t i = 0; i < reads.size(); i++) {
                // Say we are going to write out the alignment
                gam_buffer.push_back(reads[i]);
                
                // Set its path to the corrected embedded path
                *gam_buffer.back().mutable_path() = read_paths[i];
                
                // Write it back out
                stream::write_buffered(gam_out_file, gam_buffer, 100);
            }
            // Flush the buffer
            stream::write_buffered(gam_out_file, gam_buffer, 0);
        }
    } else if (augmentation_mode == "pileup") {
        // We want to augment with pileups
        
        // The PileupAugmenter object will take care of all augmentation
        PileupAugmenter augmenter(graph, PileupAugmenter::Default_default_quality, min_aug_support);    

        // compute the augmented graph from the pileup
        // Note: we can save a fair bit of memory by clearing pileups, and re-reading off of
        //       pileup_file_name
        augment_with_pileups(augmenter, *pileups, expect_subgraph, show_progress);
        delete pileups;
        pileups = nullptr;

        // write the augmented graph
        if (show_progress) {
            cerr << "Writing augmented graph" << endl;
        }
        augmenter.write_augmented_graph(cout, false);

        // write the agumented gam
        if (!gam_out_file_name.empty()) {
            ofstream gam_out_file(gam_out_file_name);
            if (!gam_out_file) {
                cerr << "[vg augment]: Error opening output GAM file: " << gam_out_file_name << endl;
                return 1;
            }
            get_input_file(gam_in_file_name, [&](istream& alignment_stream) {
                    vector<Alignment> gam_buffer;
                    function<void(Alignment&)> lambda = [&gam_out_file, &gam_buffer, &augmenter](Alignment& alignment) {
                        list<mapping_t> aug_path;
                        augmenter.map_path(alignment.path(), aug_path, true);
                        alignment.mutable_path()->clear_mapping();
                        for (auto& aug_mapping : aug_path) {
                            *alignment.mutable_path()->add_mapping() = aug_mapping.to_mapping();
                        }
                        gam_buffer.push_back(alignment);
                        stream::write_buffered(gam_out_file, gam_buffer, 100);
                    };
                    stream::for_each(alignment_stream, lambda);
                    stream::write_buffered(gam_out_file, gam_buffer, 0);
                });
        }

        // write the translation
        if (!translation_file_name.empty()) {
            // write the translations
            if (show_progress) {
                cerr << "Writing translation table" << endl;
            }
            ofstream translation_file(translation_file_name);
            if (!translation_file) {
                cerr << "[vg augment] error: error opening translation file: " << translation_file_name << endl;
                return 1;
            }
            augmenter._augmented_graph.write_translations(translation_file);
            translation_file.close();
        }

        // write the supports
        if (!support_file_name.empty()) {
            // write the supports
            if (show_progress) {
                cerr << "Writing supports" << endl;
            }
            ofstream support_file(support_file_name);
            if (!support_file) {
                cerr << "[vg augment] error: error opening supports file: " << support_file_name << endl;
                return 1;
            }
            augmenter._augmented_graph.write_supports(support_file);
            support_file.close();
        }       
    } else {
        cerr << "[vg augment] error: unrecognized augmentation mode" << endl;
        exit(1);
    }

    if (pileups != nullptr) {
        delete pileups;
        pileups = nullptr;
    }    
    
    delete graph;

    return 0;
}

Pileups* compute_pileups(VG* graph, const string& gam_file_name, int thread_count, int min_quality,
                         int max_mismatches, int window_size, int max_depth, bool use_mapq,
                         bool show_progress) {

    // Make Pileups makers for each thread.
    vector<Pileups*> pileups;
    for (int i = 0; i < thread_count; ++i) {
        pileups.push_back(new Pileups(graph, min_quality, max_mismatches, window_size, max_depth, use_mapq));
    }
    
    // setup alignment stream
    get_input_file(gam_file_name, [&](istream& alignment_stream) {
        // compute the pileups.
        if (show_progress) {
            cerr << "Computing pileups" << endl;
        }
        
        function<void(Alignment&)> lambda = [&pileups, &graph](Alignment& aln) {
            int tid = omp_get_thread_num();
            pileups[tid]->compute_from_alignment(aln);
        };
        stream::for_each_parallel(alignment_stream, lambda);
    });

    // single-threaded (!) merge
    if (show_progress && pileups.size() > 1) {
        cerr << "Merging pileups" << endl;
    }
    for (int i = 1; i < pileups.size(); ++i) {
        pileups[0]->merge(*pileups[i]);
        delete pileups[i];
    }
    return pileups[0];
}

void augment_with_pileups(PileupAugmenter& augmenter, Pileups& pileups, bool expect_subgraph,
                          bool show_progress) {
    
    if (show_progress) {
        cerr << "Computing augmented graph from the pileup" << endl;
    }

    pileups.for_each_node_pileup([&](const NodePileup& node_pileup) {
            if (!augmenter._graph->has_node(node_pileup.node_id())) {
                // This pileup doesn't belong in this graph
                if(!expect_subgraph) {
                    throw runtime_error("Found pileup for nonexistent node " + to_string(node_pileup.node_id()));
                }
                // If that's expected, just skip it
                return;
            }
            // Send approved pileups to the augmenter
            augmenter.call_node_pileup(node_pileup);
            
        });

    pileups.for_each_edge_pileup([&](const EdgePileup& edge_pileup) {
            if (!augmenter._graph->has_edge(edge_pileup.edge())) {
                // This pileup doesn't belong in this graph
                if(!expect_subgraph) {
                    throw runtime_error("Found pileup for nonexistent edge " + pb2json(edge_pileup.edge()));
                }
                // If that's expected, just skip it
                return;
            }
            // Send approved pileups to the augmenter
            augmenter.call_edge_pileup(edge_pileup);            
        });

    // map the edges from original graph
    if (show_progress) {
        cerr << "Mapping edges into augmented graph" << endl;
    }
    augmenter.update_augmented_graph();

    // map the paths from the original graph
    if (show_progress) {
        cerr << "Mapping paths into augmented graph" << endl;
    }
    augmenter.map_paths();
}

// Register subcommand
static Subcommand vg_augment("augment", "augment a graph from an alignment", PIPELINE, 5, main_augment);
