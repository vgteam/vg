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
#include "../vcf_buffer.hpp"
#include "../path_index.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_add(char** argv) {
    cerr << "usage: " << argv[0] << " add [options] old.vg >new.vg" << endl
         << "options:" << endl
         << "    -v, --vcf FILE         add in variants from the given VCF file" << endl
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
    bool show_progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"vcf", required_argument, 0, 'v'},
                {"progress",  no_argument, 0, 'p'},
                {"threads", required_argument, 0, 't'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "v:pt:h?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'v':
            vcf_filenames.push_back(optarg);
            break;

        case 'p':
            show_progress = true;
            break;
            
        case 't':
            omp_set_num_threads(atoi(optarg));
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
    
    // Load the graph
    VG* graph;
    get_input_file(optind, argc, argv, [&](istream& in) {
        graph = new VG(in, show_progress);
    });
    
    if (graph == nullptr) {
        cerr << "error:[vg add]: Could not load graph" << endl;
        exit(1);
    }
    
    {
        
        // We collect all the edits we need to make, and then make them all at
        // once.
        vector<Path> edits_to_make;
        
        // We need indexes of all the paths that variants happen on.
        map<string, PathIndex> indexes;
        
        // We need a function to grab the index for a path
        auto get_path_index = [&](const string& path_name) -> PathIndex& {
            if (!indexes.count(path_name)) {
                // Not already made. Generate it.
                indexes.emplace(piecewise_construct,
                    forward_as_tuple(path_name),
                    forward_as_tuple(*graph, path_name, true));
            }
            return indexes.at(path_name);
        };
        
        for (auto vcf_filename : vcf_filenames) {
            // For each VCF
            
            // Open it
            vcflib::VariantCallFile vcf;
            vcf.open(vcf_filename);
            if (!vcf.is_open()) {
                cerr << "error:[vg add] could not open " << vcf_filename << endl;
                return 1;
            }
            
            // Make a buffer
            VcfBuffer buffer(&vcf);
            buffer.fill_buffer();
            
            // Grab a variant
            auto* variant = buffer.get();
            while(variant != nullptr) {
                // For each variant
            
                // Where is it?
                auto& variant_path_name = variant->sequenceName;
                auto& variant_path_offset = variant->position; // Already made 0-based by the buffer
                
                auto variant_ref_length = variant->ref.size();
                
                if (!graph->paths.has_path(variant_path_name)) {
                    // Explode if the path is not in the vg
                    
                    cerr << "error:[vg add] could not find path" << variant_path_name << " in graph" << endl;
                    return 1;
                }

                // Grab the path index
                auto& index = get_path_index(variant_path_name);
            
                // Extract its left and right context from the appropriate path in the graph
                // On the left we want either 100 bases or all the bases before the first ref base.
                size_t left_context_length = max(min((int64_t)100, (int64_t) variant_path_offset - (int64_t) 1), (int64_t) 0);
                // On the right we want either 100 bases or all the bases after the last ref base.
                size_t right_context_length = min(index.sequence.size() - variant_path_offset - variant_ref_length, (size_t) 100);
            
                string left_context = index.sequence.substr(variant_path_offset - left_context_length, left_context_length);
                string right_context = index.sequence.substr(variant_path_offset + variant_ref_length, right_context_length);
                
                // Find the node that the variant falls on
                NodeSide center = index.at_position(variant_path_offset);
                
                // Extract its graph context for realignment
                VG context;
                graph->nonoverlapping_node_context_without_paths(graph->get_node(center.node), context);
                // TODO: how many nodes should this be?
                // TODO: write/copy the search from xg so we can do this by node length.
                graph->expand_context(context, 10, false);
                
                for (auto& alt : variant->alt) {
                    // For each non-ref alt
                
                    // Paste it in with the left and right context
                    stringstream alt_stream;
                    alt_stream << left_context << alt << right_context;
                    
                    // Align it to the subgraph.
                    
                    // TODO: we want to give a full-length bonus to keep the
                    // ends attached when a variant is near the end of the
                    // reference. But we can't without turning on pinned mode,
                    // which is incorrect because we don't have a read that we
                    // know runs up to the end of the graph.
                    Alignment aln = context.align(alt_stream.str(), 0, false, false, 30);
                    
                    // Add the path to our collection of paths to add
                    edits_to_make.push_back(aln.path());
                }
                
                // Move on to the next variant
                buffer.handle_buffer();
                buffer.fill_buffer();
                variant = buffer.get();
            }
            
            // Then at the end of the VCF, edit the graph
            graph->edit(edits_to_make);
            
        }
    
    }
    
    // Output the modified graph
    graph->serialize_to_ostream(std::cout);
    
    delete graph;

    // NB: If you worry about "still reachable but possibly lost" warnings in valgrind,
    // this would free all the memory used by protobuf:
    //ShutdownProtobufLibrary();

    return 0;
}

// Register subcommand
static Subcommand vg_add("add", "add variants from a VCF to a graph", main_add);

