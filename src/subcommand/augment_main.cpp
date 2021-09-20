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
#include "../xg.hpp"
#include "../vg.hpp"
#include "../augment.hpp"
#include "../packer.hpp"
#include "../io/save_handle_graph.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <handlegraph/mutable_path_mutable_handle_graph.hpp>
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"
#include <bdsg/overlays/overlay_helper.hpp>


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_augment(char** argv, ConfigurableParser& parser) {
    cerr << "usage: " << argv[0] << " augment [options] <graph.vg> [alignment.gam] > augmented_graph.vg" << endl
         << "Embed GAM alignments into a graph to facilitate variant calling" << endl
         << endl
         << "general options:" << endl
         << "    -i, --include-paths         merge the paths implied by alignments into the graph" << endl
         << "    -S, --keep-softclips        include softclips from input alignments (they are cut by default)" << endl
         << "    -B, --label-paths           don't augment with alignments, just use them for labeling the graph" << endl
         << "    -Z, --translation FILE      save translations from augmented back to base graph to FILE" << endl
         << "    -A, --alignment-out FILE    save augmented GAM reads to FILE" << endl
         << "    -F, --gaf                   expect (and write) GAF instead of GAM" << endl
         << "    -s, --subgraph              graph is a subgraph of the one used to create GAM. ignore alignments with missing nodes" << endl
         << "    -m, --min-coverage N        minimum coverage of a breakpoint required for it to be added to the graph" << endl
         << "    -c, --expected-cov N        expected coverage.  used only for memory tuning [default : 128]" << endl
         << "    -q, --min-baseq N           ignore edits whose sequence have average base quality < N" << endl
         << "    -Q, --min-mapq N            ignore alignments with mapping quality < N" << endl
         << "    -N, --max-n F               maximum fraction of N bases in an edit for it to be included [default : 0.25]" << endl
         << "    -h, --help                  print this help message" << endl
         << "    -p, --progress              show progress" << endl
         << "    -v, --verbose               print information and warnings about vcf generation" << endl
         << "    -t, --threads N             number of threads (only 1st pass with -m or -q option is multithreaded)" << endl
         << "loci file options:" << endl
         << "    -l, --include-loci FILE     merge all alleles in loci into the graph" << endl       
         << "    -L, --include-gt FILE       merge only the alleles in called genotypes into the graph" << endl;
    
     // Then report more options
     parser.print_help(cerr);
}

int main_augment(int argc, char** argv) {

    // Write the translations (as protobuf) to this path
    string translation_file_name;

    // Include a path in the graph for each GAM
    bool include_paths = false;

    // Include the softclips for each path
    bool include_softclips = false;

    // Just label the paths with the GAM
    bool label_paths = false;

    // Merge alleles from this loci file instead of GAM
    string loci_file;

    // Merge only alleles from called genotypes in the loci file
    bool called_genotypes_only = false;
    
    // Load in GAM alignments to map over to the augmented graph from here
    string gam_in_file_name;

    // Write the GAM alignments (from gam_in_file_name) projected on the augmented graph here
    string gam_out_file_name;

    // Expect given graph to be subgraph of that used to create GAM and not
    // fail when nodes are missing
    bool is_subgraph = false;

    // Min coverage for graph to be broken at a breakpoint
    // Whene non-zero, the Packer will be used to collect breakpoints
    size_t min_coverage = 0;

    // Used to set data_width for Packer
    size_t expected_coverage = 128;

    // Minimum average base quality in an edit's sequence for it to be used
    double min_baseq = 0;

    // Minimum mapping quality of an alignment for it to be used
    double min_mapq = 0;

    // Maximum fraction of Ns
    double max_frac_n = 0.25;

    // GAF format toggle
    string aln_format = "GAM";

    // Print some progress messages to screen
    bool show_progress = false;

    // Print verbose message
    bool verbose = false;

    static const struct option long_options[] = {
        // Deprecated Options
        {"augmentation-mode", required_argument, 0, 'a'},
        // General Options
        {"translation", required_argument, 0, 'Z'},
        {"alignment-out", required_argument, 0, 'A'},
        {"include-paths", no_argument, 0, 'i'},
        {"cut-softclips", no_argument, 0, 'C'},
        {"keep-softclips", no_argument, 0, 'S'},
        {"label-paths", no_argument, 0, 'B'},
        {"subgraph", no_argument, 0, 's'},
        {"min-coverage", required_argument, 0, 'm'},
        {"expected-cov", required_argument, 0, 'c'},
        {"min-baseq", required_argument, 0, 'q'},
        {"min-mapq", required_argument, 0, 'Q'},
        {"max-n", required_argument, 0, 'N'},
        {"gaf", no_argument, 0, 'F'},
        {"help", no_argument, 0, 'h'},
        {"progress", required_argument, 0, 'p'},
        {"verbose", no_argument, 0, 'v'},
        {"threads", required_argument, 0, 't'},
        // Loci Options
        {"include-loci", required_argument, 0, 'l'},
        {"include-gt", required_argument, 0, 'L'},
        {0, 0, 0, 0}
    };
    static const char* short_options = "a:Z:A:iCSBhpvt:l:L:sm:c:q:Q:N:F";
    optind = 2; // force optind past command positional arguments

    // This is our command-line parser
    ConfigurableParser parser(short_options, long_options, [&](int c) {
        // Parse all the options we have defined here.
        switch (c)
        {
            // Deprecated.
        case 'a':
            cerr << "[vg augment] warning: -a / --augmentation-mode option is deprecated" << endl;
            break;
            // General Options
        case 'Z':
            translation_file_name = optarg;
            break;
        case 'A':
            gam_out_file_name = optarg;
            break;
        case 'i':
            include_paths = true;
            break;
        case 'C':
            cerr << "[vg augment] warning: -C / --cut-softclips option is deprecated (now enabled by default)" << endl;
            break;
        case 'S':
            include_softclips = true;
            break;
        case 'B':
            label_paths = true;
            break;
        case 's':
            is_subgraph = true;
            break;
        case 'm':
            min_coverage = parse<size_t>(optarg);
            break;
        case 'c':
            expected_coverage = parse<size_t>(optarg);
            break;
        case 'q':
            min_baseq = parse<double>(optarg);
            break;
        case 'Q':
            min_mapq = parse<double>(optarg);
            break;
        case 'N':
            max_frac_n = parse<double>(optarg);
            break;
        case 'F':
            aln_format = "GAF";
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
        {
            int num_threads = parse<int>(optarg);
            if (num_threads <= 0) {
                cerr << "error:[vg call] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }            
        // Loci Options
        case 'l':
            loci_file = optarg;
            break;
        case 'L':
            loci_file = optarg;
            called_genotypes_only = true;
            break;
            
        default:
          abort ();
        }
    });

    // Parse the command line options, updating optind.
    parser.parse(argc, argv);

    // Parse the two positional arguments
    if (optind + 1 > argc) {
        cerr << "[vg augment] error: too few arguments" << endl;
        help_augment(argv, parser);
        return 1;
    }

    string graph_file_name = get_input_file_name(optind, argc, argv);
    if (optind < argc) {
        gam_in_file_name = get_input_file_name(optind, argc, argv);
    }

    if (gam_in_file_name.empty() && loci_file.empty()) {
        cerr << "[vg augment] error: gam file argument required" << endl;
        return 1;
    }
    if (gam_in_file_name == "-" && graph_file_name == "-") {
        cerr << "[vg augment] error: graph and gam can't both be from stdin." << endl;
        return 1;
    }
    if (label_paths && (!gam_out_file_name.empty() || !translation_file_name.empty())) {
        cerr << "[vg augment] error: Translation (-Z) and GAM (-A) output do not work with \"label-only\" (-B) mode" << endl;
        return 1;
    }
    if (gam_in_file_name == "-" && !label_paths) {
        cerr << "[vg augment] warning: reading the entire GAM from stdin into memory.  it is recommended to pass in"
             << " a filename rather than - so it can be streamed over two passes" << endl;
        if (!gam_out_file_name.empty()) {
            cerr << "             warning: when streaming in a GAM with -A, the output GAM will lose all non-Path related fields from the input" << endl;
        }
    }

    // read the graph
    if (show_progress) {
        cerr << "Reading input graph" << endl;
    }

    // Read the graph
    unique_ptr<MutablePathMutableHandleGraph> graph;
    graph = vg::io::VPKG::load_one<MutablePathMutableHandleGraph>(graph_file_name);
    VG* vg_graph = dynamic_cast<VG*>(graph.get());
    HandleGraph* vectorizable_graph = nullptr;
    unique_ptr<Packer> packer;
    bdsg::VectorizableOverlayHelper overlay_helper;
    
    if (label_paths) {
        // Just add path names with extend()
        function<void(Alignment&)> lambda = [&](Alignment& alignment) {
            if (!include_softclips) {
                softclip_trim(alignment);
            }
            Path simplified_path = simplify(alignment.path());
            *simplified_path.mutable_name() = alignment.name();
            add_path_to_graph(graph.get(), simplified_path);
        };
        if (aln_format == "GAM") {
            get_input_file(gam_in_file_name, [&](istream& alignment_stream) {
                    vg::io::for_each<Alignment>(alignment_stream, lambda);
                });
        } else {
            assert(aln_format == "GAF");
            vg::io::gaf_unpaired_for_each(*graph, gam_in_file_name, lambda);
        }
        if (vg_graph != nullptr) {
            vg_graph->paths.sort_by_mapping_rank();
            vg_graph->paths.rebuild_mapping_aux();
        }
    }
    else {
        // the packer's required for any kind of filtering logic -- so we use it when
        // baseq is present as well, or n-fraction.
        if (min_coverage > 0 || min_baseq || max_frac_n < 1.) {
            vectorizable_graph = dynamic_cast<HandleGraph*>(overlay_helper.apply(graph.get()));
            size_t data_width = Packer::estimate_data_width(expected_coverage);
            size_t bin_count = Packer::estimate_bin_count(get_thread_count());
            packer = make_unique<Packer>(vectorizable_graph, true, false, false, false, 0, bin_count, data_width);
            // makes sure filters are activated. 
            min_coverage = max(size_t(min_coverage), size_t(1));
        }
    
        // Actually do augmentation
        vector<Translation> translation;
        if (!gam_out_file_name.empty()) {
            ofstream gam_out_file(gam_out_file_name);
            if (!gam_out_file) {
                cerr << "[vg augment] error: could not open output GAM file: " << gam_out_file_name << endl;
                return 1;
            }
        }
        if (gam_in_file_name == "-" || !loci_file.empty()) {
            vector<Path> buffer;
            if (gam_in_file_name == "-") {
                // this is usually bad news, but we gave a warning
                function<void(Alignment&)> lambda = [&](Alignment& alignment) {
                    Path& path = *alignment.mutable_path();
                    path.set_name(alignment.name());
                    buffer.push_back(path);
                };
                if (aln_format == "GAM") {
                    get_input_file(gam_in_file_name, [&](istream& alignment_stream) {
                            vg::io::for_each<Alignment>(alignment_stream, lambda);
                        });
                } else {
                    assert(aln_format == "GAF");
                    vg::io::gaf_unpaired_for_each(*graph, gam_in_file_name, lambda);
                }
            } else if (!loci_file.empty()) {
                function<void(Locus&)> lambda = [&graph, &buffer, &called_genotypes_only](Locus& locus) {
                    // if we are only doing called genotypes, record so we can filter alleles
                    set<int> alleles_in_genotype;
                    if (called_genotypes_only) {
                        for (int i = 0; i < locus.genotype_size(); ++i) {
                            for (int j = 0; j < locus.genotype(i).allele_size(); ++j) {
                                alleles_in_genotype.insert(locus.genotype(i).allele(j));
                            }
                        }
                    }
                    for (int i = 0; i < locus.allele_size(); ++i) {
                        // skip alleles not in the genotype if using only called genotypes
                        if (!alleles_in_genotype.empty()) {
                            if (!alleles_in_genotype.count(i)) continue;
                        }
                        Path path = simplify(locus.allele(i));
                        stringstream name;
                        name << locus.name() << ":" << i;
                        path.set_name(name.str());
                        buffer.push_back(path);
                    }
                };
                get_input_file(loci_file, [&](istream& loci_stream) {
                        vg::io::for_each(loci_stream, lambda);
                    });
            }
            
            augment(graph.get(),
                    buffer,
                    aln_format, 
                    translation_file_name.empty() ? nullptr : &translation,
                    gam_out_file_name,
                    include_paths,
                    include_paths,
                    !include_softclips,
                    is_subgraph,
                    min_baseq,
                    min_mapq,
                    packer.get(),
                    min_coverage,
                    max_frac_n);
        } else {
            // much better to stream from a file so we can do two passes without storing in memory
            augment(graph.get(),
                    gam_in_file_name,
                    aln_format,
                    translation_file_name.empty() ? nullptr : &translation,
                    gam_out_file_name,
                    include_paths,
                    include_paths,
                    !include_softclips,
                    is_subgraph,
                    min_baseq,
                    min_mapq,
                    packer.get(),
                    min_coverage,
                    max_frac_n);
        }

        // we don't have a streaming interface for translation:  write the buffer now
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
            vg::io::write_buffered(translation_file, translation, 0);
            translation_file.close();
        }
    } 

    // Serialize the graph using VPKG.
    vg::io::save_handle_graph(graph.get(), cout);
    
    return 0;
}

// Register subcommand
static Subcommand vg_augment("augment", "augment a graph from an alignment", PIPELINE, 8, main_augment);
