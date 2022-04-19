// index.cpp: define the "vg index" subcommand, which makes xg, GCSA2, and GBWT indexes

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <random>
#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../haplotype_indexer.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../io/save_handle_graph.hpp"
#include "../stream_index.hpp"
#include "../vg_set.hpp"
#include "../utility.hpp"
#include "../region.hpp"
#include "../snarls.hpp"
#include "../min_distance.hpp"
#include "../source_sink_overlay.hpp"
#include "../gbwt_helper.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gcsa_helper.hpp"

#include <gcsa/algorithms.h>
#include <gbwt/variants.h>
#include <bdsg/overlays/packed_subgraph_overlay.hpp>
#include <handlegraph/algorithms/weakly_connected_components.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_index(char** argv) {
    cerr << "usage: " << argv[0] << " index [options] <graph1.vg> [graph2.vg ...]" << endl
         << "Creates an index on the specified graph or graphs. All graphs indexed must " << endl
         << "already be in a joint ID space." << endl
         << "general options:" << endl
         << "    -b, --temp-dir DIR     use DIR for temporary files" << endl
         << "    -t, --threads N        number of threads to use" << endl
         << "    -p, --progress         show progress" << endl
         << "xg options:" << endl
         << "    -x, --xg-name FILE     use this file to store a succinct, queryable version of the graph(s), or read for GCSA or distance indexing" << endl
         << "    -L, --xg-alts          include alt paths in xg" << endl
         << "gbwt options (more in vg gbwt):" << endl
         << "    -v, --vcf-phasing FILE generate threads from the haplotypes in the VCF file FILE" << endl
         << "    -W, --ignore-missing   don't warn when variants in the VCF are missing from the graph; silently skip them" << endl
         << "    -T, --store-threads    generate threads from the embedded paths" << endl
         << "    -M, --store-gam FILE   generate threads from the alignments in gam FILE (many allowed)" << endl
         << "    -F, --store-gaf FILE   generate threads from the alignments in gaf FILE (many allowed)" << endl
         << "    -G, --gbwt-name FILE   store the threads as GBWT in FILE" << endl
         << "    -z, --actual-phasing   do not make unphased homozygous genotypes phased"<< endl
         << "    -P, --force-phasing    replace unphased genotypes with randomly phased ones" << endl
         << "    -o, --discard-overlaps skip overlapping alternate alleles if the overlap cannot be resolved" << endl
         << "    -B, --batch-size N     number of samples per batch (default 200)" << endl
         << "    -u, --buffer-size N    GBWT construction buffer size in millions of nodes (default 100)" << endl
         << "    -n, --id-interval N    store haplotype ids at one out of N positions (default 1024)" << endl
         << "    -R, --range X..Y       process samples X to Y (inclusive)" << endl
         << "    -r, --rename V=P       rename contig V in the VCFs to path P in the graph (may repeat)" << endl
         << "    --rename-variants      when renaming contigs, find variants in the graph based on the new name" << endl
         << "    -I, --region C:S-E     operate on only the given 1-based region of the given VCF contig (may repeat)" << endl
         << "    -E, --exclude SAMPLE   exclude any samples with the given name from haplotype indexing" << endl
         << "gcsa options:" << endl
         << "    -g, --gcsa-out FILE    output a GCSA2 index to the given file" << endl
         //<< "    -i, --dbg-in FILE      use kmers from FILE instead of input VG (may repeat)" << endl
         << "    -f, --mapping FILE     use this node mapping in GCSA2 construction" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph (default " << gcsa::Key::MAX_LENGTH << ")" << endl
         << "    -X, --doubling-steps N use this number of doubling steps for GCSA2 construction (default " << gcsa::ConstructionParameters::DOUBLING_STEPS << ")" << endl
         << "    -Z, --size-limit N     limit temporary disk space usage to N gigabytes (default " << gcsa::ConstructionParameters::SIZE_LIMIT << ")" << endl
         << "    -V, --verify-index     validate the GCSA2 index using the input kmers (important for testing)" << endl
         << "gam indexing options:" << endl
         << "    -l, --index-sorted-gam input is sorted .gam format alignments, store a GAI index of the sorted GAM in INPUT.gam.gai" << endl
         << "vg in-place indexing options:" << endl
         << "    --index-sorted-vg      input is ID-sorted .vg format graph chunks, store a VGI index of the sorted vg in INPUT.vg.vgi" << endl
         << "snarl distance index options" << endl
         << "    -s  --snarl-name FILE  load snarls from FILE (snarls must include trivial snarls)" << endl
         << "    -j  --dist-name FILE   use this file to store a snarl-based distance index" << endl
         << "    -w  --max_dist N       cap beyond which the maximum distance is no longer accurate. If this is not included or is 0, don't build maximum distance index" << endl;
}

void multiple_thread_sources() {
    std::cerr << "error: [vg index] cannot generate threads from multiple sources (VCF, GAM, GAF, paths)" << std::endl;
    std::cerr << "error: [vg index] GBWT indexes can be built separately and merged with vg gbwt" << std::endl;
    std::exit(EXIT_FAILURE);
}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    #define OPT_BUILD_VGI_INDEX  1000
    #define OPT_RENAME_VARIANTS  1001

    // Which indexes to build.
    bool build_xg = false, build_gbwt = false, build_gcsa = false, build_dist = false;

    // Files we should read.
    string vcf_name, mapping_name;
    vector<string> dbg_names;

    // Files we should write.
    string xg_name, gbwt_name, gcsa_name, dist_name, snarl_name;

    // General
    bool show_progress = false;

    // GBWT
    HaplotypeIndexer haplotype_indexer;
    enum thread_source_type { thread_source_none, thread_source_vcf, thread_source_paths, thread_source_gam, thread_source_gaf };
    thread_source_type thread_source = thread_source_none;
    vector<string> aln_file_names;

    // GCSA
    gcsa::size_type kmer_size = gcsa::Key::MAX_LENGTH;
    gcsa::ConstructionParameters params;
    bool verify_gcsa = false;
    
    // Gam index (GAI)
    bool build_gai_index = false;
    
    // VG in-place index (VGI)
    bool build_vgi_index = false;

    //Distance index
    int cap = -1;
    bool include_maximum = false;

    // Include alt paths in xg
    bool xg_alts = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            // General
            {"temp-dir", required_argument, 0, 'b'},
            {"threads", required_argument, 0, 't'},
            {"progress",  no_argument, 0, 'p'},

            // XG
            {"xg-name", required_argument, 0, 'x'},
            {"thread-db", required_argument, 0, 'F'},
            {"xg-alts", no_argument, 0, 'L'},

            // GBWT
            {"vcf-phasing", required_argument, 0, 'v'},
            {"ignore-missing", no_argument, 0, 'W'},
            {"store-threads", no_argument, 0, 'T'},
            {"store-gam", required_argument, 0, 'M'},
            {"store-gaf", required_argument, 0, 'F'},
            {"gbwt-name", required_argument, 0, 'G'},
            {"actual-phasing", no_argument, 0, 'z'},
            {"force-phasing", no_argument, 0, 'P'},
            {"discard-overlaps", no_argument, 0, 'o'},
            {"batch-size", required_argument, 0, 'B'},
            {"buffer-size", required_argument, 0, 'u'},
            {"id-interval", required_argument, 0, 'n'},
            {"range", required_argument, 0, 'R'},
            {"rename", required_argument, 0, 'r'},
            {"rename-variants", no_argument, 0, OPT_RENAME_VARIANTS},
            {"region", required_argument, 0, 'I'},
            {"exclude", required_argument, 0, 'E'},

            // GCSA
            {"gcsa-out", required_argument, 0, 'g'},
            {"dbg-in", required_argument, 0, 'i'},
            {"mapping", required_argument, 0, 'f'},
            {"kmer-size", required_argument, 0, 'k'},
            {"doubling-steps", required_argument, 0, 'X'},
            {"size-limit", required_argument, 0, 'Z'},
            {"verify-index", no_argument, 0, 'V'},
            
            // GAM index (GAI)
            {"index-sorted-gam", no_argument, 0, 'l'},
            
            // VG in-place index (VGI)
            {"index-sorted-vg", no_argument, 0, OPT_BUILD_VGI_INDEX},

            //Snarl distance index
            {"snarl-name", required_argument, 0, 's'},
            {"dist-name", required_argument, 0, 'j'},
            {"max-dist", required_argument, 0, 'w'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "b:t:px:Lv:WTM:F:G:zPoB:u:n:R:r:I:E:g:i:f:k:X:Z:Vls:j:w:h",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'b':
            temp_file::set_dir(optarg);
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;
        case 'p':
            show_progress = true;
            haplotype_indexer.show_progress = true;
            break;

        // XG
        case 'x':
            build_xg = true;
            xg_name = optarg;
            break;
        case 'L':
            xg_alts = true;
            break;

        // GBWT
        case 'v':
            if (thread_source != thread_source_none) {
                multiple_thread_sources();
            }
            thread_source = thread_source_vcf;
            vcf_name = optarg;
            break;
        case 'W':
            haplotype_indexer.warn_on_missing_variants = false;
            break;
        case 'T':
            if (thread_source != thread_source_none) {
                multiple_thread_sources();
            }
            thread_source = thread_source_paths;
            break;
        case 'M':
            if (thread_source != thread_source_none && thread_source != thread_source_gam) {
                multiple_thread_sources();
            }
            thread_source = thread_source_gam;
            build_gbwt = true;
            aln_file_names.push_back(optarg);
            break;
        case 'F':
            if (thread_source != thread_source_none && thread_source != thread_source_gaf) {
                multiple_thread_sources();
            }
            thread_source = thread_source_gaf;
            build_gbwt = true;
            aln_file_names.push_back(optarg);
            break;
        case 'G':
            build_gbwt = true;
            gbwt_name = optarg;
            break;
        case 'z':
            haplotype_indexer.phase_homozygous = false;
            break;
        case 'P':
            haplotype_indexer.force_phasing = true;
            break;
        case 'o':
            haplotype_indexer.discard_overlaps = true;
            break;
        case 'B':
            haplotype_indexer.samples_in_batch = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'u':
            haplotype_indexer.gbwt_buffer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'n':
            haplotype_indexer.id_interval = parse<size_t>(optarg);
            break;
        case 'R':
            {
                // Parse first..last
                string temp(optarg);
                size_t found = temp.find("..");
                if(found == string::npos || found == 0 || found + 2 == temp.size()) {
                    cerr << "error: [vg index] could not parse range " << temp << endl;
                    exit(1);
                }
                haplotype_indexer.sample_range.first = parse<size_t>(temp.substr(0, found));
                haplotype_indexer.sample_range.second = parse<size_t>(temp.substr(found + 2)) + 1;
            }
            break;
        case 'r':
            {
                // Parse the rename old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error: [vg index] could not parse rename " << key_value << endl;
                    exit(1);
                }
                // Parse out the two parts
                string vcf_contig = key_value.substr(0, found);
                string graph_contig = key_value.substr(found + 1);
                // Add the name mapping
                haplotype_indexer.path_to_vcf[graph_contig] = vcf_contig;
            }
            break;
        case OPT_RENAME_VARIANTS:
            haplotype_indexer.rename_variants = true;
            break;
        case 'I':
            {
                // We want to parse this region specifier
                string region(optarg);
                
                Region parsed;
                parse_region(region, parsed);
                if (parsed.start <= 0 || parsed.end <= 0) {
                    // We need both range bounds, and we can't accept 0 since input is 1-based.
                    cerr << "error: [vg index] could not parse 1-based region " << optarg << endl;
                }
                
                // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                haplotype_indexer.regions[parsed.seq] = make_pair((size_t) (parsed.start - 1), (size_t) parsed.end);
            }
            break;
        case 'E':
            haplotype_indexer.excluded_samples.insert(optarg);
            break;

        // GCSA
        case 'g':
            build_gcsa = true;
            gcsa_name = optarg;
            break;
        case 'i':
            cerr << "warning: -i option is deprecated" << endl;
            dbg_names.push_back(optarg);
            break;
        case 'f':
            mapping_name = optarg;
            break;
        case 'k':
            kmer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'X':
            params.setSteps(parse<size_t>(optarg));
            break;
        case 'Z':
            params.setLimit(parse<size_t>(optarg));
            break;
        case 'V':
            verify_gcsa = true;
            break;
            
        // Gam index (GAI)
        case 'l':
            build_gai_index = true;
            break;
            
        // VGI index
        case OPT_BUILD_VGI_INDEX:
            build_vgi_index = true;
            break;

        //Snarl distance index
        case 's':
            snarl_name = optarg;
            break;
        case 'j':
            build_dist = true;
            dist_name = optarg;
            break;
        case 'w':
            build_dist = true;
            cap = parse<int>(optarg);
            include_maximum = true;
            break;

        case 'h':
        case '?':
            help_index(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    vector<string> file_names;
    while (optind < argc) {
        string file_name = get_input_file_name(optind, argc, argv);
        file_names.push_back(file_name);
    }


    if (xg_name.empty() && gbwt_name.empty() &&
        gcsa_name.empty() && !build_gai_index && !build_vgi_index && dist_name.empty()) {
        cerr << "error: [vg index] index type not specified" << endl;
        return 1;
    }

    if (build_gbwt && thread_source == thread_source_none) {
        cerr << "error: [vg index] cannot build GBWT without threads" << endl;
        return 1;
    }

    if (thread_source != thread_source_none && !build_gbwt) {
        cerr << "error: [vg index] no GBWT output specified for the threads" << endl;
        return 1;
    }

    if (thread_source == thread_source_gam || thread_source == thread_source_gaf) {
        for (const auto& name : aln_file_names) {
            if (name == "-") {
                cerr << "error: [vg index] GAM (-M) and GAF (-F) input files cannot be read from stdin (-)" << endl;
                return 1;
            }
        }
    }

    if (thread_source != thread_source_none && file_names.size() != 1) {
        cerr << "error: [vg index] exactly one graph required for generating threads" << std::endl;
        cerr << "error: [vg index] you may combine the graphs with vg index -x combined.xg --xg-alts" << std::endl;
        return 1;
    }

    if (file_names.size() <= 0 && dbg_names.empty()){
        //cerr << "No graph provided for indexing. Please provide a .vg file or GCSA2-format deBruijn graph to index." << endl;
        //return 1;
    }
    
    if (file_names.size() != 1 && build_gai_index) {
        cerr << "error: [vg index] can only index exactly one sorted GAM file at a time" << endl;
        return 1;
    }
    
    if (file_names.size() != 1 && build_vgi_index) {
        cerr << "error: [vg index] can only index exactly one sorted VG file at a time" << endl;
        return 1;
    }
    
    if (file_names.size() > 1 && build_dist) {
        // Allow zero filenames for the index-from-xg mode
        cerr << "error: [vg index] can only create one distance index at a time" << endl;
        return 1;
    }
    
    if (build_gcsa && kmer_size > gcsa::Key::MAX_LENGTH) {
        cerr << "error: [vg index] GCSA2 cannot index with kmer size greater than " << gcsa::Key::MAX_LENGTH << endl;
        return 1;
    }
    
    if (build_xg && build_gcsa && file_names.empty()) {
        // Really we want to build a GCSA by *reading* and XG
        build_xg = false;
        // We'll continue in the build_gcsa section
        std::cerr << "warning: [vg index] providing input XG with option -x is deprecated" << std::endl;
    }
    if (build_dist && file_names.empty()) {
        //If we want to build the distance index from the xg
        build_xg = false;
        std::cerr << "warning: [vg index] providing input XG with option -x is deprecated" << std::endl;
    }


    // Build XG. Include alt paths in the XG if requested with -L.
    if (build_xg) {
        if (file_names.empty()) {
            // VGset or something segfaults when we feed it no graphs.
            cerr << "error: [vg index] at least one graph is required to build an xg index" << endl;
            return 1;
        }
        if (show_progress) {
            cerr << "Building XG index" << endl;
        }
        xg::XG xg_index;
        VGset graphs(file_names);
        graphs.to_xg(xg_index, (xg_alts ? [](const string&) {return false;} : Paths::is_alt), nullptr);
        if (show_progress) {
            cerr << "Saving XG index to " << xg_name << endl;
        }
        // Save the XG.
        vg::io::save_handle_graph(&xg_index, xg_name);
    }

    // Generate threads
    if (thread_source != thread_source_none) {

        // Load the only input graph.
        unique_ptr<PathHandleGraph> path_handle_graph;
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(file_names[0]);

        std::unique_ptr<gbwt::DynamicGBWT> gbwt_index(nullptr);
        if (thread_source == thread_source_vcf) {
            std::vector<std::string> parse_files = haplotype_indexer.parse_vcf(vcf_name, *path_handle_graph);
            path_handle_graph.reset(); // Save memory by deleting the graph.
            gbwt_index = haplotype_indexer.build_gbwt(parse_files);
        } else if (thread_source == thread_source_paths) {
            gbwt_index = haplotype_indexer.build_gbwt(*path_handle_graph);
        } else if (thread_source == thread_source_gam) {
            gbwt_index = haplotype_indexer.build_gbwt(*path_handle_graph, aln_file_names, "GAM");
        } else if (thread_source == thread_source_gaf) {
            gbwt_index = haplotype_indexer.build_gbwt(*path_handle_graph, aln_file_names, "GAF");
        }
        if (build_gbwt && gbwt_index.get() != nullptr) {
            save_gbwt(*gbwt_index, gbwt_name, show_progress);
        }
    } // End of thread indexing.

    // Build GCSA
    if (build_gcsa) {

        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if (!show_progress) {
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        }

        double start = gcsa::readTimer();

        // Generate temporary kmer files
        bool delete_kmer_files = false;
        if (dbg_names.empty()) {
            if (show_progress) {
                cerr << "Generating kmer files..." << endl;
            }
            
            if (!file_names.empty()) {
                // Get the kmers from a VGset.
                VGset graphs(file_names);
                size_t kmer_bytes = params.getLimitBytes();
                dbg_names = graphs.write_gcsa_kmers_binary(kmer_size, kmer_bytes);
                params.reduceLimit(kmer_bytes);
                delete_kmer_files = true;
            } else if (!xg_name.empty()) {
                // Get the kmers from an XG or other single graph
                
                // Load the graph
                auto single_graph = vg::io::VPKG::load_one<HandleGraph>(xg_name);
                
                auto make_kmers_for_component = [&](const HandleGraph* g) {
                    // Make an overlay on it to add source and sink nodes
                    // TODO: Don't use this directly; unify this code with VGset's code.
                    SourceSinkOverlay overlay(g, kmer_size);
                    
                    // Get the size limit
                    size_t kmer_bytes = params.getLimitBytes();
                    
                    // Write the kmer temp file
                    dbg_names.push_back(write_gcsa_kmers_to_tmpfile(overlay, kmer_size, kmer_bytes,
                        overlay.get_id(overlay.get_source_handle()),
                        overlay.get_id(overlay.get_sink_handle())));
                        
                    // Feed back into the size limit
                    params.reduceLimit(kmer_bytes);
                    delete_kmer_files = true;
                };
                
                if (show_progress) {
                    cerr << "Finding connected components..." << endl;
                }
                
                // Get all the components in the graph, which we can process separately to save memory.
                std::vector<std::unordered_set<nid_t>> components = handlealgs::weakly_connected_components(single_graph.get());
                
                if (components.size() == 1) {
                    // Only one component
                    if (show_progress) {
                        cerr << "Processing single component graph..." << endl;
                    }
                    make_kmers_for_component(single_graph.get());
                } else {
                    for (size_t i = 0; i < components.size(); i++) {
                        // Run separately on each component.
                        // Don't run in parallel or size limit tracking won't work.
                
                        if (show_progress) {
                            cerr << "Selecting component " << i << "/" << components.size() << "..." << endl;
                        }
                        
                        bdsg::PackedSubgraphOverlay component_graph(single_graph.get());
                        for (auto& id : components[i]) {
                            // Add each node to the subgraph.
                            // TODO: use a handle-returning component
                            // finder so we don't need to get_handle here.
                            component_graph.add_node(single_graph->get_handle(id, false));
                        }
                        
                        if (show_progress) {
                            cerr << "Processing component " << i << "/" << components.size() << "..." << endl;
                        }
                        
                        make_kmers_for_component(&component_graph);
                    }
                }
            } else {
                cerr << "error: [vg index] cannot generate GCSA index without either a vg or an xg" << endl;
                exit(1);
            }
        }

        // Build the index
        if (show_progress) {
            cerr << "Building the GCSA2 index..." << endl;
        }
        gcsa::InputGraph input_graph(dbg_names, true, gcsa::Alphabet(), mapping_name);
        gcsa::GCSA gcsa_index(input_graph, params);
        gcsa::LCPArray lcp_array(input_graph, params);
        if (show_progress) {
            double seconds = gcsa::readTimer() - start;
            cerr << "GCSA2 index built in " << seconds << " seconds, "
                 << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
            cerr << "I/O volume: " << gcsa::inGigabytes(gcsa::readVolume()) << " GB read, "
                 << gcsa::inGigabytes(gcsa::writeVolume()) << " GB write" << endl;
        }

        // Save the indexes
        save_gcsa(gcsa_index, gcsa_name, show_progress);
        save_lcp(lcp_array, gcsa_name + ".lcp", show_progress);

        // Verify the index
        if (verify_gcsa) {
            if (show_progress) {
                cerr << "Verifying the index..." << endl;
            }
            if (!gcsa::verifyIndex(gcsa_index, &lcp_array, input_graph)) {
                cerr << "warning: [vg index] GCSA2 index verification failed" << endl;
            }
        }

        // Delete the temporary kmer files
        if (delete_kmer_files) {
            for (auto& filename : dbg_names) {
                temp_file::remove(filename);
            }
        }
    }
    
    if (build_gai_index) {
        // Index a sorted GAM file.
        
        get_input_file(file_names.at(0), [&](istream& in) {
            // Grab the input GAM stream and wrap it in a cursor
            vg::io::ProtobufIterator<Alignment> cursor(in);
            
            // Index the file
            StreamIndex<Alignment> index;
            index.index(cursor);
 
            // Save the GAM index in the appropriate place.
            // TODO: Do we really like this enforced naming convention just beacuse samtools does it?
            ofstream index_out(file_names.at(0) + ".gai");
            if (!index_out.good()) {
                cerr << "error: [vg index] could not open " << file_names.at(0) << ".gai" << endl;
                exit(1);
            }
            index.save(index_out);
        });
    }
    
    if (build_vgi_index) {
        // Index an ID-sorted VG file.
        get_input_file(file_names.at(0), [&](istream& in) {
            // Grab the input VG stream and wrap it in a cursor
            vg::io::ProtobufIterator<Graph> cursor(in);
            
            // Index the file
            StreamIndex<Graph> index;
            index.index(cursor);
 
            // Save the index in the appropriate place.
            // TODO: Do we really like this enforced naming convention just beacuse samtools does it?
            ofstream index_out(file_names.at(0) + ".vgi");
            if (!index_out.good()) {
                cerr << "error: [vg index] could not open " << file_names.at(0) << ".vgi" << endl;
                exit(1);
            }
            index.save(index_out);
        });
        
    }

    //Build new snarl-based minimum distance index
    if (build_dist) {
        if (file_names.empty() && xg_name.empty()) {
            cerr << "error: [vg index] one graph is required to build a distance index" << endl;
            return 1;
        } else if (file_names.size() > 1 || (file_names.size() == 1 && !xg_name.empty())) {
            cerr << "error: [vg index] only one graph at a time can be used to build a distance index" << endl;
        } else if (dist_name.empty()) {
            cerr << "error: [vg index] distance index requires an output file" << endl;
            return 1;
        } else if (snarl_name.empty()) {
            cerr << "error: [vg index] distance index requires a snarl file" << endl;
            return 1;
            
        } else {
            //Get snarl manager
            ifstream snarl_stream(snarl_name);
            if (!snarl_stream) {
                cerr << "error: [vg index] cannot open Snarls file" << endl;
                exit(1);
            }
            SnarlManager* snarl_manager = new SnarlManager(snarl_stream);
            snarl_stream.close();

            //Get graph and build dist index
            if (file_names.empty() && !xg_name.empty()) {
                // We were given a -x specifically to read as XG
                
                auto xg = vg::io::VPKG::load_one<xg::XG>(xg_name);

                // Create the MinimumDistanceIndex
                MinimumDistanceIndex di(xg.get(), snarl_manager);
                // Save the completed DistanceIndex
                vg::io::VPKG::save(di, dist_name);

            } else {
                // May be GBZ or a HandleGraph.
                auto options = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, handlegraph::HandleGraph>(file_names.at(0));
                
                if (get<0>(options)) {
                    // We have a GBZ graph
                    auto& gbz = get<0>(options);
                    
                    // Create the MinimumDistanceIndex
                    MinimumDistanceIndex di(&(gbz->graph), snarl_manager);
                    vg::io::VPKG::save(di, dist_name);
                } else if (get<1>(options)) {
                    // We were given a graph generically
                    auto& graph = get<1>(options);
                    
                    // Create the MinimumDistanceIndex
                    MinimumDistanceIndex di(graph.get(), snarl_manager);
                    vg::io::VPKG::save(di, dist_name);
                } else {
                    cerr << "error: [vg index] input is not a graph or GBZ" << endl;
                    return 1;
                }
            }
        }

    }
    if (show_progress) {
        cerr << "Memory usage: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
    }
    return 0;
}

// Register subcommand
static Subcommand vg_construct("index", "index graphs or alignments for random access or mapping", PIPELINE, 4, main_index);
