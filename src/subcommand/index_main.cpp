// index.cpp: define the "vg index" subcommand, which makes xg, GCSA2, GBWT, and RocksDB indexes

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <random>
#include <string>
#include <vector>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../index.hpp"
#include "../index_manager.hpp"
#include "xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../stream_index.hpp"
#include "../vg_set.hpp"
#include "../utility.hpp"
#include "../region.hpp"
#include "../snarls.hpp"
#include "../min_distance.hpp"
#include "../source_sink_overlay.hpp"
#include "../gbwt_helper.hpp"

#include <gcsa/gcsa.h>
#include <gcsa/algorithms.h>
#include <gbwt/variants.h>

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
         << "gbwt options:" << endl
         << "    -v, --vcf-phasing FILE generate threads from the haplotypes in the VCF file FILE" << endl
         << "    -W, --ignore-missing   don't warn when variants in the VCF are missing from the graph; silently skip them" << endl
         << "    -e, --parse-only FILE  store the VCF parsing with prefix FILE without generating threads" << endl
         << "    -T, --store-threads    generate threads from the embedded paths" << endl
         << "    -M, --store-gam FILE   generate threads from the alignments in FILE (many allowed)" << endl
         << "    -G, --gbwt-name FILE   store the threads as GBWT in FILE" << endl
         << "    -H, --write-haps FILE  store the threads as sequences in FILE" << endl
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
         << "    -i, --dbg-in FILE      use kmers from FILE instead of input VG (may repeat)" << endl
         << "    -f, --mapping FILE     use this node mapping in GCSA2 construction" << endl
         << "    -k, --kmer-size N      index kmers of size N in the graph (default " << gcsa::Key::MAX_LENGTH << ")" << endl
         << "    -X, --doubling-steps N use this number of doubling steps for GCSA2 construction (default " << gcsa::ConstructionParameters::DOUBLING_STEPS << ")" << endl
         << "    -Z, --size-limit N     limit temporary disk space usage to N gigabytes (default " << gcsa::ConstructionParameters::SIZE_LIMIT << ")" << endl
         << "    -V, --verify-index     validate the GCSA2 index using the input kmers (important for testing)" << endl
         << "gam indexing options:" << endl
         << "    -l, --index-sorted-gam input is sorted .gam format alignments, store a GAI index of the sorted GAM in INPUT.gam.gai" << endl
         << "vg in-place indexing options:" << endl
         << "    --index-sorted-vg      input is ID-sorted .vg format graph chunks, store a VGI index of the sorted vg in INPUT.vg.vgi" << endl
         << "rocksdb options:" << endl
         << "    -d, --db-name  <X>     store the RocksDB index in <X>" << endl
         << "    -m, --store-mappings   input is .gam format, store the mappings in alignments by node" << endl
         << "    -a, --store-alignments input is .gam format, store the alignments by node" << endl
         << "    -A, --dump-alignments  graph contains alignments, output them in sorted order" << endl
         << "    -N, --node-alignments  input is (ideally, sorted) .gam format," << endl
         << "                           cross reference nodes by alignment traversals" << endl
         << "    -D, --dump             print the contents of the db to stdout" << endl
         << "    -C, --compact          compact the index into a single level (improves performance)" << endl
         << "snarl distance index options" << endl
         << "    -s  --snarl-name FILE  load snarls from FILE" << endl
         << "    -j  --dist-name FILE   use this file to store a snarl-based distance index" << endl
         << "    -w  --max_dist N       cap beyond which the maximum distance is no longer accurate. If this is not included or is 0, don't build maximum distance index" << endl;
}

int main_index(int argc, char** argv) {

    if (argc == 2) {
        help_index(argv);
        return 1;
    }

    #define OPT_BUILD_VGI_INDEX 1000
    #define OPT_RENAME_VARIANTS 1001

    // Which indexes to build.
    bool build_xg = false, build_gbwt = false, write_threads = false, build_gcsa = false, build_rocksdb = false, build_dist = false;

    // Files we should read.
    string vcf_name, mapping_name;
    vector<string> dbg_names;

    // Files we should write.
    string xg_name, gbwt_name, parse_name, threads_name, gcsa_name, rocksdb_name, dist_name, snarl_name;

    // General
    bool show_progress = false;

    // GBWT
    bool warn_on_missing_variants = true;
    bool index_haplotypes = false, index_paths = false, index_gam = false;
    vector<string> gam_file_names;
    bool phase_homozygous = true, force_phasing = false, discard_overlaps = false;
    size_t samples_in_batch = 200;
    size_t gbwt_buffer_size = gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION; // Millions of nodes.
    size_t id_interval = gbwt::DynamicGBWT::SAMPLE_INTERVAL;
    std::pair<size_t, size_t> sample_range(0, ~(size_t)0); // The semiopen range of samples to process.
    map<string, string> path_to_vcf; // Path name conversion from --rename.
    bool rename_variants = false;
    map<string, pair<size_t, size_t>> regions; // Region restrictions for contigs, in VCF name space, as 0-based exclusive-end ranges.
    unordered_set<string> excluded_samples; // Excluded sample names from --exclude.

    // GCSA
    gcsa::size_type kmer_size = gcsa::Key::MAX_LENGTH;
    gcsa::ConstructionParameters params;
    bool verify_gcsa = false;
    
    // Gam index (GAI)
    bool build_gai_index = false;
    
    // VG in-place index (VGI)
    bool build_vgi_index = false;

    // RocksDB
    bool dump_index = false;
    bool store_alignments = false;
    bool store_node_alignments = false;
    bool store_mappings = false;
    bool dump_alignments = false;

    //Distance index
    int cap = -1;
    bool include_maximum = false;

    // Unused?
    bool compact = false;

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
            {"parse-only", required_argument, 0, 'e'},
            {"store-threads", no_argument, 0, 'T'},
            {"store-gam", required_argument, 0, 'M'},
            {"gbwt-name", required_argument, 0, 'G'},
            {"write-haps", required_argument, 0, 'H'},
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

            // RocksDB
            {"db-name", required_argument, 0, 'd'},
            {"store-mappings", no_argument, 0, 'm'},
            {"store-alignments", no_argument, 0, 'a'},
            {"dump-alignments", no_argument, 0, 'A'},
            {"node-alignments", no_argument, 0, 'N'},
            {"dump", no_argument, 0, 'D'},
            {"compact", no_argument, 0, 'C'},

            //Snarl distance index
            {"snarl-name", required_argument, 0, 's'},
            {"dist-name", required_argument, 0, 'j'},
            {"max-dist", required_argument, 0, 'w'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "b:t:px:Lv:We:TM:G:H:zPoB:u:n:R:r:I:E:g:i:f:k:X:Z:Vld:maANDCs:j:w:h",
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
            index_haplotypes = true;
            build_xg = true;
            vcf_name = optarg;
            break;
        case 'W':
            warn_on_missing_variants = false;
            break;
        case 'e':
            parse_name = optarg;
            break;
        case 'T':
            index_paths = true;
            build_xg = true;
            break;
        case 'M':
            index_gam = true;
            build_gbwt = true;
            gam_file_names.push_back(optarg);
            break;
        case 'G':
            build_gbwt = true;
            gbwt_name = optarg;
            break;
        case 'H':
            write_threads = true;
            threads_name = optarg;
            break;
        case 'z':
            phase_homozygous = false;
            break;
        case 'P':
            force_phasing = true;
            break;
        case 'o':
            discard_overlaps = true;
            break;
        case 'B':
            samples_in_batch = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'u':
            gbwt_buffer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'n':
            id_interval = parse<size_t>(optarg);
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
                sample_range.first = parse<size_t>(temp.substr(0, found));
                sample_range.second = parse<size_t>(temp.substr(found + 2)) + 1;
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
                path_to_vcf[graph_contig] = vcf_contig;
            }
            break;
        case OPT_RENAME_VARIANTS:
            rename_variants = true;
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
                regions[parsed.seq] = make_pair((size_t) (parsed.start - 1), (size_t) parsed.end);
            }
            break;
        case 'E':
            excluded_samples.insert(optarg);
            break;

        // GCSA
        case 'g':
            build_gcsa = true;
            gcsa_name = optarg;
            break;
        case 'i':
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

        // RocksDB
        case 'd':
            build_rocksdb = true;
            rocksdb_name = optarg;
            break;
        case 'm':
            store_mappings = true;
            break;
        case 'a':
            store_alignments = true;
            break;
        case 'A':
            dump_alignments = true;
            break;
        case 'N':
            store_node_alignments = true;
            break;
        case 'D':
            dump_index = true;
            break;

        case 'C':
            compact = true;
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

    if (xg_name.empty() && gbwt_name.empty() && parse_name.empty() && threads_name.empty() &&
        gcsa_name.empty() && rocksdb_name.empty() && !build_gai_index && !build_vgi_index && dist_name.empty()) {
        cerr << "error: [vg index] index type not specified" << endl;
        return 1;
    }

    if ((build_gbwt || write_threads) && !(index_haplotypes || index_paths || index_gam)) {
        cerr << "error: [vg index] cannot build GBWT without threads" << endl;
        return 1;
    }
    
    if (build_gbwt && write_threads) {
        cerr << "error: [vg index] cannot build GBWT and dump threads simultaneously" << endl;
        return 1;
    }

    if (!parse_name.empty() && (!index_haplotypes || index_paths || index_gam)) {
        cerr << "error: [vg index] --parse-only works with --vcf-phasing only" << endl;
        return 1;
    }

    if ((index_haplotypes || index_paths || index_gam) && !(build_gbwt || !parse_name.empty() || write_threads)) {
        cerr << "error: [vg index] no output format specified for the threads" << endl;
        return 1;
    }

    if (index_gam && (index_haplotypes || index_paths)) {
        cerr << "error: [vg index] GAM threads are incompatible with haplotype/path threads" << endl;
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
    }
    if (build_dist && file_names.empty()) {
        //If we want to build the distance index from the xg
        build_xg = false;
    }

    // Build XG
    xg::XG* xg_index = new xg::XG();
    map<string, Path> alt_paths;
    if (build_xg) {
        if (file_names.empty()) {
            // VGset or something segfaults when we feed it no graphs.
            cerr << "error: [vg index] at least one graph is required to build an xg index" << endl;
            return 1;
        }
        VGset graphs(file_names);
        graphs.to_xg(*xg_index, xg_alts ? [](const string&) {return false;} : Paths::is_alt, index_haplotypes ? &alt_paths : nullptr);
        if (index_haplotypes && xg_alts) {
            assert(alt_paths.empty());
            // they weren't filtered in the above: re-extract here
            xg_index->for_each_path_handle([&](path_handle_t path_handle) {
                    string path_name = xg_index->get_path_name(path_handle);
                    if (Paths::is_alt(path_name)) {
                        alt_paths[path_name] = path_from_path_handle(*xg_index, path_handle);
                    }
                });
        }
        if (show_progress) {
            cerr << "Built base XG index" << endl;
        }
    }

#ifdef debug
    cerr << "Alt paths:" << endl;
    for (auto& kv : alt_paths) {
        cerr << kv.first << ": " << kv.second.mapping_size() << " entries" << endl;
    }
#endif

    

    // Generate threads
    if (index_haplotypes || index_paths || index_gam) {
    
        // Make the gadget that will find the paths in the VCF.
        HaplotypeIndexer indexer;
        // Fill in its parameters.
        // TODO: fill these directly from the command line.
        indexer.show_progress = show_progress;
        indexer.warn_on_missing_variants = warn_on_missing_variants;
        indexer.path_to_vcf = path_to_vcf;
        indexer.rename_variants = rename_variants;
        indexer.batch_file_prefix = parse_name;
        indexer.index_paths = index_paths;
        indexer.phase_homozygous = phase_homozygous;
        indexer.force_phasing = force_phasing;
        indexer.discard_overlaps = discard_overlaps;
        indexer.samples_in_batch = samples_in_batch;
        indexer.gbwt_buffer_size = gbwt_buffer_size;
        indexer.id_interval = id_interval;
        indexer.sample_range = sample_range;
        indexer.regions = regions;
        indexer.excluded_samples = excluded_samples;
    
        if (!parse_name.empty()) {
            // Only generate parse files for the VCF. Don't do anything else.
            
            vcflib::VariantCallFile variant_file;
            variant_file.parseSamples = false; // vcflib parsing is very slow if there are many samples.
            variant_file.open(vcf_name);
            if (!variant_file.is_open()) {
                cerr << "error: [vg index] could not open " << vcf_name << endl;
                return 1;
            } else if (show_progress) {
                cerr << "Opened variant file " << vcf_name << endl;
            }
            
            // Process each VCF contig corresponding to an XG path.
            vector<path_handle_t> path_handles;
            // 1st pass: scan for all non-alt paths (they are handled separately)
            xg_index->for_each_path_handle([&](path_handle_t path_handle) {
                    if (!alt_paths.count(xg_index->get_path_name(path_handle))) {
                        path_handles.push_back(path_handle);
                    }
                });

            // Make a place to dump sample names, which we ignore.
            std::vector<std::string> sample_names;

            // Run VCF parsing but do nothing with the generated phasing batches.
            // This will write all the parse files for us.
            indexer.parse_vcf(xg_index, alt_paths, path_handles, variant_file, sample_names,
                [&](size_t contig, const gbwt::VariantPaths& variants, gbwt::PhasingInformation& phasings_batch) {});
                
        } else if (write_threads) {
            // Dump threads instead of building a GBWT
            
            // File to write out haplotypes to.
            gbwt::text_buffer_type binary_file;
            
            indexer.generate_threads(xg_index, alt_paths,
                index_paths, vcf_name, gam_file_names, [&](size_t id_width) {
                
                // We got the ID width
                if (show_progress) { cerr << "Writing the threads to " << threads_name << endl; }
                binary_file = gbwt::text_buffer_type(threads_name, std::ios::out, gbwt::MEGABYTE, id_width);
            }, [&](const gbwt::vector_type& thread, const gbwt::size_type (&thread_name)[4]) {
                // We got a thread!
                // Dump all the thread visits
                for (auto node : thread) { binary_file.push_back(node); }
                binary_file.push_back(gbwt::ENDMARKER);
            });
            binary_file.close();
        } else {
            // Do actual GBWT indexing.
            
            unique_ptr<gbwt::DynamicGBWT> gbwt = indexer.build_gbwt(xg_index, alt_paths, index_paths, vcf_name, gam_file_names);
            
            if (show_progress) {
                cerr << "Saving GBWT to disk..." << endl;
            }
                
            // Save encapsulated in a VPKG
            vg::io::VPKG::save(*gbwt, gbwt_name);
                
        }
    } // End of thread indexing.

    // Save XG
    if (build_xg && !xg_name.empty()) {
        if (show_progress) {
            cerr << "Saving XG index to disk..." << endl;
        }
        // Save encapsulated in a VPKG
        vg::io::VPKG::save(*xg_index, xg_name); 
    }
    delete xg_index; xg_index = nullptr;

    // Build GCSA
    if (build_gcsa) {

        // Configure GCSA2 verbosity so it doesn't spit out loads of extra info
        if (!show_progress) {
            gcsa::Verbosity::set(gcsa::Verbosity::SILENT);
        }

        // Use the same temp directory as VG.
        gcsa::TempFile::setDirectory(temp_file::get_dir());

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
                // Get the kmers from an XG
                
                get_input_file(xg_name, [&](istream& xg_stream) {
                    // Load the XG
                    auto xg = vg::io::VPKG::load_one<xg::XG>(xg_stream);
                
                    // Make an overlay on it to add source and sink nodes
                    // TODO: Don't use this directly; unify this code with VGset's code.
                    SourceSinkOverlay overlay(xg.get(), kmer_size);
                    
                    // Get the size limit
                    size_t kmer_bytes = params.getLimitBytes();
                    
                    // Write just the one kmer temp file
                    dbg_names.push_back(write_gcsa_kmers_to_tmpfile(overlay, kmer_size, kmer_bytes,
                        overlay.get_id(overlay.get_source_handle()),
                        overlay.get_id(overlay.get_sink_handle())));
                        
                    // Feed back into the size limit
                    params.reduceLimit(kmer_bytes);
                    delete_kmer_files = true;
                
                });
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
        if (show_progress) {
            cerr << "Saving the index to disk..." << endl;
        }
        vg::io::VPKG::save(gcsa_index, gcsa_name);
        vg::io::VPKG::save(lcp_array, gcsa_name + ".lcp");

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

    if (build_rocksdb) {

        Index index;

        if (compact) {
            index.open_for_write(rocksdb_name);
            index.compact();
            index.flush();
            index.close();
        }

        if (store_node_alignments && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            int64_t aln_idx = 0;
            function<void(Alignment&)> lambda = [&index,&aln_idx](Alignment& aln) {
                index.cross_alignment(aln_idx++, aln);
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (store_alignments && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                index.put_alignment(aln);
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (dump_alignments) {
            vector<Alignment> output_buf;
            index.open_read_only(rocksdb_name);
            auto lambda = [&output_buf](const Alignment& aln) {
                output_buf.push_back(aln);
                vg::io::write_buffered(cout, output_buf, 100);
            };
            index.for_each_alignment(lambda);
            vg::io::write_buffered(cout, output_buf, 0);
            index.close();
        }

        if (store_mappings && file_names.size() > 0) {
            index.open_for_bulk_load(rocksdb_name);
            function<void(Alignment&)> lambda = [&index](Alignment& aln) {
                const Path& path = aln.path();
                for (int i = 0; i < path.mapping_size(); ++i) {
                    index.put_mapping(path.mapping(i));
                }
            };
            for (auto& file_name : file_names) {
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each(in, lambda);
                });
            }
            index.flush();
            index.close();
        }

        if (dump_index) {
            index.open_read_only(rocksdb_name);
            index.dump(cout);
            index.close();
        }

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
                
                ifstream xg_stream(xg_name);
                auto xg = vg::io::VPKG::load_one<xg::XG>(xg_stream);

                // Create the MinimumDistanceIndex
                MinimumDistanceIndex di(xg.get(), snarl_manager);
                // Save the completed DistanceIndex
                vg::io::VPKG::save(di, dist_name);

            } else {
                // We were given a graph generically
                auto graph = vg::io::VPKG::load_one<handlegraph::HandleGraph>(file_names.at(0));
    
                // Create the MinimumDistanceIndex
                MinimumDistanceIndex di(graph.get(), snarl_manager);
                vg::io::VPKG::save(di, dist_name);
            }
          
            
        }

    }
    if (show_progress) {
        cerr << "Memory usage: " << gcsa::inGigabytes(gcsa::memoryUsage()) << " GB" << endl;
    }
    return 0;
}

// Register subcommand
static Subcommand vg_construct("index", "index graphs or alignments for random access or mapping", PIPELINE, 2, main_index);
