/** \file add_main.cpp
 *
 * Defines the "vg add" subcommand, which adds in variation from a VCF to an
 * existing graph.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../variant_adder.hpp"

#include <vg/io/vpkg.hpp>




using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_add(char** argv) {
    cerr << "usage: " << argv[0] << " add [options] old.vg >new.vg" << endl
         << "options:" << endl
         << "    -v, --vcf FILE         add in variants from the given VCF file (may repeat)" << endl
         << "    -n, --rename V=G       rename contig V in the VCFs to contig G in the graph (may repeat)" << endl
         << "    -i, --ignore-missing   ignore contigs in the VCF not found in the graph" << endl
         << "    -r, --variant-range N  range in which to look for nearby variants to make a haplotype" << endl
         << "    -f, --flank-range N    extra flanking sequence to use outside of found variants" << endl
         << "    -p, --progress         show progress" << endl
         << "    -t, --threads N        use N threads (defaults to numCPUs)" << endl;
}

int main_add(int argc, char** argv) {

    if (argc == 2) {
        help_add(argv);
        return 1;
    }

    
    // We can have one or more VCFs
    vector<string> vcf_filenames;
    // And one or more renames
    vector<pair<string, string>> renames;
    bool show_progress = false;
    bool ignore_missing = false;
    int variant_range = -1;
    int flank_range = -1;
    
    // TODO: make variant_adder not hold on to its graph so tightly, so we can
    // set its settings as we parse the options;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"vcf", required_argument, 0, 'v'},
                {"rename", required_argument, 0, 'n'},
                {"ignore-missing", no_argument, 0, 'i'},
                {"variant-range", required_argument, 0, 'r'},
                {"flank-range", required_argument, 0, 'f'},
                {"progress",  no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:n:r:f:ipt:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'v':
            vcf_filenames.push_back(optarg);
            break;
            
        case 'n':
            {
                // Parse the rename old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error:[vg add] could not parse rename " << key_value << endl;
                    exit(1);
                }
                // Parse out the two parts
                string vcf_contig = key_value.substr(0, found);
                string graph_contig = key_value.substr(found + 1);
                // Add the name mapping
                renames.emplace_back(vcf_contig, graph_contig);
            }
            break;

        case 'i':
            ignore_missing = true;
            break;
            
        case 'r':
            variant_range = parse<int>(optarg);
            break;
            
        case 'f':
            flank_range = parse<int>(optarg);
            break;

        case 'p':
            show_progress = true;
            break;
            
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_add(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }
    
    // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
    gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
    
    // Turn on nested parallelism, so we can parallelize over VCFs and over alignment bands
    omp_set_nested(1);
    
    // Open all the VCFs into a list
    vector<unique_ptr<vcflib::VariantCallFile>> vcfs;
    
    for (auto vcf_filename : vcf_filenames) {
        // For each VCF filename
        
        // Open it
        vcfs.emplace_back(new vcflib::VariantCallFile());
        auto& vcf = *vcfs.back();
        vcf.open(vcf_filename);
        if (!vcf.is_open()) {
            cerr << "error:[vg add] could not open " << vcf_filename << endl;
            return 1;
        }
    }
    
    // Load the graph

    unique_ptr<handlegraph::MutablePathDeletableHandleGraph> graph;
    string graph_filename = get_input_file_name(optind, argc, argv);
    graph = vg::io::VPKG::load_one<handlegraph::MutablePathDeletableHandleGraph>(graph_filename);
    
    VG* vg_graph = dynamic_cast<vg::VG*>(graph.get());
    
    // Call this to populate the vg_graph if it isn't populated.
    auto ensure_vg = [&]() -> vg::VG* {
        if (vg_graph == nullptr) {
            // Copy instead.
            vg_graph = new vg::VG();
            handlealgs::copy_path_handle_graph(graph.get(), vg_graph);
            // Give the unique_ptr ownership and delete the graph we loaded.
            graph.reset(vg_graph);
            // Make sure the paths are all synced up
            vg_graph->paths.to_graph(vg_graph->graph);
        }
        return vg_graph;
    };
    
    // TODO: We need to move VariantAdder away from vg::VG eventually.
    // Right now we always need the vg format graph.
    // TODO: deduplicate ensure_vg with other subcommands?
    ensure_vg();
    
    if (vg_graph == nullptr) {
        cerr << "error:[vg add]: Could not load graph" << endl;
        exit(1);
    }
    
    {
        // Clear existing path ranks (since we invalidate them)
        vg_graph->paths.clear_mapping_ranks();
    
        // Make a VariantAdder for the graph
        VariantAdder adder(*vg_graph);
        // Report updates when running interactively
        adder.print_updates = true;
        
        // Set up parameters
        adder.ignore_missing_contigs = ignore_missing;
        if (variant_range != -1) {
            adder.variant_range = variant_range;
        }
        if (flank_range != -1) {
            adder.flank_range = flank_range;
        }
        
        for (auto& rename : renames) {
            // Set up all the VCF contig renames from the command line
            adder.add_name_mapping(rename.first, rename.second);
        }

        #pragma omp parallel for
        for (size_t i = 0; i < vcfs.size(); i++) {
            // For each VCF
            auto& vcf = *vcfs[i];
            
            // Add the variants from the VCF to the graph, at the same
            // time as other VCFs.
            adder.add_variants(&vcf);        
        }
        
        // TODO: should we sort the graph?
        
        // Rebuild all the path ranks and stuff
        vg_graph->paths.rebuild_mapping_aux();
    }
        
    // Output the modified graph
    vg_graph->serialize_to_ostream(std::cout);
    
    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_add("add", "add variants from a VCF to a graph", DEPRECATED, main_add);

