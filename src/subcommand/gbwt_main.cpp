/** \file gbwt_main.cpp
 *
 * Defines the "vg gbwt" subcommand, which wraps up access for commands we'd otherwise find
 * in the gbwt submodule.  */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"
#include "../gbwt_helper.hpp"
#include "../haplotype_indexer.hpp"
#include "../path.hpp"
#include "../region.hpp"

#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/path_cover.h>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <unistd.h>

enum build_mode { build_none, build_vcf, build_paths, build_alignments };
enum merge_mode { merge_none, merge_insert, merge_fast };
enum path_cover_mode { path_cover_none, path_cover_augment, path_cover_local, path_cover_greedy };
enum index_type { index_none, index_compressed, index_dynamic };

void load_gbwt(const std::string& filename, gbwt::GBWT& index, bool show_progress);
void load_gbwt(const std::string& filename, gbwt::DynamicGBWT& index, bool show_progress);

void get_compressed(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress);
void get_dynamic(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress);

void use_or_save(std::unique_ptr<gbwt::DynamicGBWT>& index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, std::vector<std::string>& filenames, size_t i, bool show_progress);

void get_graph(std::unique_ptr<PathHandleGraph>& graph, bool& in_use, const std::string& filename, bool show_progress);
void clear_graph(std::unique_ptr<PathHandleGraph>& graph, bool& in_use);

void help_gbwt(char** argv) {
    std::cerr << "usage: " << argv[0] << " gbwt [options] [args]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Manipulate GBWTs. Input GBWTs are loaded from input args or built in Step 1." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "    -x, --xg-name FILE      read the graph from FILE" << std::endl;
    std::cerr << "    -o, --output FILE       write output GBWT to FILE" << std::endl;
    std::cerr << "    -d, --temp-dir DIR      use directory DIR for temporary files" << std::endl;
    std::cerr << "    -p, --progress          show progress and statistics" << std::endl;
    std::cerr << "    -t, --threads N         use N parallel threads (in -v, -A, or -r; default " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GBWT construction parameters:" << std::endl;
    std::cerr << "    -b, --buffer-size N     GBWT construction buffer size in millions of nodes (default " << (gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION) << ")" << std::endl;
    std::cerr << "    -i, --id-interval N     store path ids at one out of N positions (default " << gbwt::DynamicGBWT::SAMPLE_INTERVAL << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 1: Build input GBWTs (requires -o, -x, and one of { -v, -E, A }):" << std::endl;
    std::cerr << "    -v, --vcf-input         index the haplotypes in the VCF files specified in input args in parallel" << std::endl;
    std::cerr << "                            (implies -f for multiple inputs unless -m is specified)" << std::endl;
    std::cerr << "        --parse-only FILE   store the VCF parses with prefix FILE without GBWTs (skip subsequent steps)" << std::endl;
    std::cerr << "        --ignore-missing    do not warn when variants are missing from the graph" << std::endl;
    std::cerr << "        --actual-phasing    do not interpret unphased homozygous genotypes as phased" << std::endl;
    std::cerr << "        --force-phasing     replace unphased genotypes with randomly phased ones" << std::endl;
    std::cerr << "        --discard-overlaps  skip overlapping alternate alleles if the overlap cannot be resolved" << std::endl;
    std::cerr << "                            instead of creating a phase break" << std::endl;
    std::cerr << "        --batch-size N      index the haplotypes in batches of N samples (default 200)" << std::endl; // FIXME source for the default
    std::cerr << "        --sample-range X-Y  index samples X to Y (inclusive, 0-based)" << std::endl;
    std::cerr << "        --rename V=P        VCF contig V matches path P in the graph (may repeat)" << std::endl;
    std::cerr << "        --rename-variants   the graph stores variants using path names instead of VCF contig names" << std::endl;
    std::cerr << "        --vcf-region C:X-Y  restrict VCF contig to coordinates X to Y (inclusive, 1-based; may repeat)" << std::endl;
    std::cerr << "        --exclude-sample X  do not index the sample with name X (faster than -R; may repeat)" << std::endl;
    std::cerr << "    -E, --index-paths       index the embedded non-alt paths in the graph (no input args)" << std::endl;
    std::cerr << "        --paths-as-samples  each path becomes a sample instead of a contig in the metadata" << std::endl;
    std::cerr << "    -A, --alignment-input   index the alignments in the GAF files specified in input args in parallel" << std::endl;
    std::cerr << "                            (implies -m for multiple inputs unless -f is specified)" << std::endl;
    std::cerr << "        --gam-format        the input files are in GAM format instead of GAF format" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 2: Merge multiple input GBWTs (requires -o; use deps/gbwt/merge_gbwt for more options):" << std::endl;
    std::cerr << "    -m, --merge             use the insertion algorithm" << std::endl;
    std::cerr << "    -f, --fast              fast merging algorithm (node ids must not overlap)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 3: Remove samples (requires -o and one input GBWT):" << std::endl;
    std::cerr << "    -R, --remove-sample X   remove the sample with name X from the index" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 4: Path cover (requires -o, -x, and one of { -a, -l, -P }):" << std::endl;
    std::cerr << "    -a, --augment-gbwt      add a path cover of missing components (one input GBWT)" << std::endl;
    std::cerr << "    -l, --local-haplotypes  sample local haplotypes (one input GBWT)" << std::endl;
    std::cerr << "    -P, --path-cover        build a greedy path cover (no input GBWTs)" << std::endl;
    std::cerr << "    -n, --num-paths N       find N paths per component (default " << gbwtgraph::LOCAL_HAPLOTYPES_DEFAULT_N << " for -l, " << gbwtgraph::PATH_COVER_DEFAULT_N << " otherwise)" << std::endl;
    std::cerr << "    -k, --context-length N  use N-node contexts (default " << gbwtgraph::PATH_COVER_DEFAULT_K << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 5: GBWTGraph construction (requires -x and one input GBWT):" << std::endl;
    std::cerr << "    -g, --graph-name FILE   build GBWTGraph and store it in FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 6: R-index construction (one input GBWT):" << std::endl;
    std::cerr << "    -r, --r-index FILE      build an r-index and store it in FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 7: Metadata (one input GBWT; use deps/gbwt/metadata_tool to modify):" << std::endl;
    std::cerr << "    -M, --metadata          print basic metadata" << std::endl;
    std::cerr << "    -C, --contigs           print the number of contigs" << std::endl;
    std::cerr << "    -H, --haplotypes        print the number of haplotypes" << std::endl;
    std::cerr << "    -S, --samples           print the number of samples" << std::endl;
    std::cerr << "    -L, --list-names        list contig/sample names (use with -C or -S)" << std::endl;
    std::cerr << "    -T, --thread-names      list thread names" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 8: Threads (one input GBWT):" << std::endl;
    std::cerr << "    -c, --count-threads     print the number of threads" << std::endl;
    std::cerr << "    -e, --extract FILE      extract threads in SDSL format to FILE" << std::endl;
    std::cerr << std::endl;
}


int main_gbwt(int argc, char** argv)
{
    if (argc == 2) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Long options in GBWT construction.
    constexpr int OPT_PARSE_ONLY = 1000;
    constexpr int OPT_IGNORE_MISSING = 1001;
    constexpr int OPT_ACTUAL_PHASING = 1002;
    constexpr int OPT_FORCE_PHASING = 1003;
    constexpr int OPT_DISCARD_OVERLAPS = 1004;
    constexpr int OPT_BATCH_SIZE = 1005;
    constexpr int OPT_SAMPLE_RANGE = 1006;
    constexpr int OPT_RENAME = 1007;
    constexpr int OPT_RENAME_VARIANTS = 1008;
    constexpr int OPT_VCF_REGION = 1009;
    constexpr int OPT_EXCLUDE_SAMPLE = 1010;
    constexpr int OPT_PATHS_AS_SAMPLES = 1011;
    constexpr int OPT_GAM_FORMAT = 1012;

    // Requirements and modes.
    bool produces_one_gbwt = false;
    build_mode build = build_none;
    merge_mode merge = merge_none;
    path_cover_mode path_cover = path_cover_none;
    bool metadata_mode = false, thread_mode = false;

    // Input GBWT construction.
    HaplotypeIndexer haplotype_indexer;
    bool gam_format = false;

    // Other parameters and flags.
    bool show_progress = false;
    bool count_threads = false;
    bool metadata = false, contigs = false, haplotypes = false, samples = false, list_names = false, thread_names = false;
    size_t num_paths = gbwtgraph::PATH_COVER_DEFAULT_N, context_length = gbwtgraph::PATH_COVER_DEFAULT_K;
    bool num_paths_set = false;

    // File/sample names.
    std::string gbwt_output, thread_output;
    std::string graph_output, xg_name;
    std::string r_index_name;
    std::string to_remove;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                // General
                { "xg-name", required_argument, 0, 'x' },
                { "output", required_argument, 0, 'o' },
                { "temp-dir", required_argument, 0, 'd' },
                { "progress",  no_argument, 0, 'p' },
                { "threads", required_argument, 0, 't' },

                // GBWT construction parameters
                { "buffer-size", required_argument, 0, 'b' },
                { "id-interval", required_argument, 0, 'i' },

                // Input GBWT construction
                { "vcf-input", no_argument, 0, 'v' },
                { "parse-only", required_argument, 0, OPT_PARSE_ONLY },
                { "ignore-missing", no_argument, 0, OPT_IGNORE_MISSING },
                { "actual-phasing", no_argument, 0, OPT_ACTUAL_PHASING },
                { "force-phasing", no_argument, 0, OPT_FORCE_PHASING },
                { "discard-overlaps", no_argument, 0, OPT_DISCARD_OVERLAPS },
                { "batch-size", required_argument, 0, OPT_BATCH_SIZE },
                { "sample-range", required_argument, 0, OPT_SAMPLE_RANGE },
                { "rename", required_argument, 0, OPT_RENAME },
                { "rename-variants", no_argument, 0, OPT_RENAME_VARIANTS },
                { "vcf-region", required_argument, 0, OPT_VCF_REGION },
                { "exclude-sample", required_argument, 0, OPT_EXCLUDE_SAMPLE },
                { "index-paths", no_argument, 0, 'E' },
                { "paths-as-samples", no_argument, 0, OPT_PATHS_AS_SAMPLES },
                { "alignment-input", no_argument, 0, 'A' },
                { "gam-format", no_argument, 0, OPT_GAM_FORMAT },

                // Merging
                { "merge", no_argument, 0, 'm' },
                { "fast", no_argument, 0, 'f' },

                // Remove sample
                { "remove-sample", required_argument, 0, 'R' },

                // Path cover
                { "augment-gbwt", no_argument, 0, 'a' },
                { "local-haplotypes", no_argument, 0, 'l' },
                { "path-cover", no_argument, 0, 'P' },
                { "num-paths", required_argument, 0, 'n' },
                { "context-length", required_argument, 0, 'k' },

                // GBWTGraph
                { "graph-name", required_argument, 0, 'g' },

                // R-index
                { "r-index", required_argument, 0, 'r' },

                // Metadata
                { "metadata", no_argument, 0, 'M' },
                { "contigs", no_argument, 0, 'C' },
                { "haplotypes", no_argument, 0, 'H' },
                { "samples", no_argument, 0, 'S' },
                { "list-names", no_argument, 0, 'L' },
                { "thread-names", no_argument, 0, 'T' },

                // Threads
                { "count-threads", no_argument, 0, 'c' },
                { "extract", required_argument, 0, 'e' },

                { "help", no_argument, 0, 'h' },
                { 0, 0, 0, 0 }
            };

        int option_index = 0;
        c = getopt_long(argc, argv, "x:o:d:pt:b:i:vEAmfR:alPn:k:g:r:MCHSLTce:h?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'x':
            xg_name = optarg;
            break;
        case 'o':
            gbwt_output = optarg;
            break;
        case 'd':
            temp_file::set_dir(optarg);
            break;
        case 'p':
            show_progress = true;
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        // GBWT construction parameters
        case 'b':
            haplotype_indexer.gbwt_buffer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case 'i':
            haplotype_indexer.id_interval = parse<size_t>(optarg);
            break;

        // Input GBWT construction
        case 'v':
            assert(build == build_none);
            build = build_vcf;
            break;
        case OPT_PARSE_ONLY:
            haplotype_indexer.batch_file_prefix = optarg;
            break;
        case OPT_IGNORE_MISSING:
            haplotype_indexer.warn_on_missing_variants = false;
            break;
        case OPT_ACTUAL_PHASING:
            haplotype_indexer.phase_homozygous = false;
            break;
        case OPT_FORCE_PHASING:
            haplotype_indexer.force_phasing = true;
            break;
        case OPT_DISCARD_OVERLAPS:
            haplotype_indexer.discard_overlaps = true;
            break;
        case OPT_BATCH_SIZE:
            haplotype_indexer.samples_in_batch = std::max(parse<size_t>(optarg), 1ul);
            break;
        case OPT_SAMPLE_RANGE:
            {
                // Parse first-last
                string range(optarg);
                size_t found = range.find("-");
                if(found == std::string::npos || found == 0 || found + 1 == range.size()) {
                    cerr << "error: [vg gbwt] could not parse range " << range << endl;
                    std::exit(EXIT_FAILURE);
                }
                haplotype_indexer.sample_range.first = parse<size_t>(range.substr(0, found));
                haplotype_indexer.sample_range.second = parse<size_t>(range.substr(found + 1)) + 1;
            }
            break;
        case OPT_RENAME:
            {
                // Parse old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error: [vg gbwt] could not parse rename " << key_value << endl;
                    std::exit(EXIT_FAILURE);
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
        case OPT_VCF_REGION:
            {
                // Parse contig:first-last
                std::string region(optarg);
                Region parsed;
                parse_region(region, parsed);
                if (parsed.start <= 0 || parsed.end <= 0) {
                    // We need both range bounds, and we can't accept 0 since input is 1-based.
                    cerr << "error: [vg gbwt] could not parse 1-based region " << region << endl;
                }
                // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                haplotype_indexer.regions[parsed.seq] = std::make_pair((size_t) (parsed.start - 1), (size_t) parsed.end);
            }
            break;
        case OPT_EXCLUDE_SAMPLE:
            haplotype_indexer.excluded_samples.insert(optarg);
            break;
        case 'E':
            assert(build == build_none);
            build = build_paths;
            produces_one_gbwt = true;
            break;
        case OPT_PATHS_AS_SAMPLES:
            haplotype_indexer.paths_as_samples = true;
            break;
        case 'A':
            assert(build == build_none);
            build = build_alignments;
            break;
        case OPT_GAM_FORMAT:
            gam_format = true;
            break;

        // Merging
        case 'm':
            merge = merge_insert;
            produces_one_gbwt = true;
            break;
        case 'f':
            merge = merge_fast;
            produces_one_gbwt = true;
            break;

        // Remove sample
        case 'R':
            to_remove = optarg;
            break;

        // Path cover
        case 'a':
            assert(path_cover == path_cover_none);
            path_cover = path_cover_augment;
            break;
        case 'l':
            assert(path_cover == path_cover_none);
            path_cover = path_cover_local;
            if (!num_paths_set) {
                num_paths = gbwtgraph::LOCAL_HAPLOTYPES_DEFAULT_N;
            }
            break;
        case 'P':
            assert(path_cover == path_cover_none);
            path_cover = path_cover_greedy;
            produces_one_gbwt = true;
            break;
        case 'n':
            num_paths = parse<size_t>(optarg);
            num_paths_set = true;
            break;
        case 'k':
            context_length = parse<size_t>(optarg);
            break;

        // GBWTGraph
        case 'g':
            graph_output = optarg;
            break;

        // Build r-index
        case 'r':
            r_index_name = optarg;
            break;

        // Metadata
        case 'M':
            metadata = true;
            metadata_mode = true;
            break;
        case 'C':
            contigs = true;
            metadata_mode = true;
            break;
        case 'H':
            haplotypes = true;
            metadata_mode = true;
            break;
        case 'S':
            samples = true;
            metadata_mode = true;
            break;
        case 'L':
            list_names = true;
            metadata_mode = true;
            break;
        case 'T':
            thread_names = true;
            metadata_mode = true;
            break;

        // Threads
        case 'c':
            count_threads = true;
            thread_mode = true;
            break;
        case 'e':
            thread_output = optarg;
            thread_mode = true;
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_gbwt(argv);
            std::exit(EXIT_FAILURE);
            break;

        default:
            std::exit(EXIT_FAILURE);
        }
    }
    std::vector<std::string> input_filenames;
    for (int arg = optind; arg < argc; arg++) {
        input_filenames.push_back(argv[arg]);
    }
    std::string gbwt_name = (input_filenames.empty() ? "" : input_filenames.front());

    // Sanity checks.
    if (build != build_none) {
        if (build == build_vcf) {
            if (input_filenames.empty()) {
                std::cerr << "error: [vg gbwt]: GBWT construction from VCF files requires input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (input_filenames.size() > 1 && merge == merge_none && haplotype_indexer.batch_file_prefix.empty()) {
                merge = merge_fast;
            }
        } else if (build == build_alignments) {
            if (input_filenames.empty()) {
                std::cerr << "error: [vg gbwt]: GBWT construction from alignments requires input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (input_filenames.size() > 1 && merge == merge_none) {
                merge = merge_insert;
            }
        } else if (build == build_paths) {
            if (!input_filenames.empty()) {
                std::cerr << "error: [vg gbwt]: GBWT construction from embedded paths does not use input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        if (xg_name.empty() || gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt]: GBWT construction requires a graph and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (merge != merge_none) {
        if (input_filenames.size() < 2 || gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt]: merging requires multiple input GBWTs and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (!to_remove.empty()) {
        if (input_filenames.size() != 1 || gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt]: removing a sample requires one input GBWT and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (path_cover != path_cover_none) {
        if (gbwt_output.empty() || xg_name.empty()) {
            std::cerr << "error: [vg gbwt]: path cover options require output GBWT and a graph" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (path_cover == path_cover_greedy && !input_filenames.empty()) {
            std::cerr << "error: [vg gbwt]: greedy path cover does not use input GBWTs" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if ((path_cover == path_cover_local || path_cover == path_cover_augment) && !(input_filenames.size() == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: path cover options -a and -l require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (num_paths == 0) {
            std::cerr << "error: [vg gbwt] number of paths must be non-zero for path cover" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (context_length < gbwtgraph::PATH_COVER_MIN_K) {
            std::cerr << "error: [vg gbwt] context length must be at least " << gbwtgraph::PATH_COVER_MIN_K << " for path cover" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (!graph_output.empty()) {
        if (xg_name.empty() || !(input_filenames.size() == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: GBWTGraph construction requires a graph and one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (!r_index_name.empty()) {
        if (!(input_filenames.size() == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: r-index construction requires one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (metadata_mode) {
        if (!(input_filenames.size() == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: metadata operations one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (thread_mode) {
        if (!(input_filenames.size() == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: thread operations one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }


    // Let GBWT operate silently and use the same temporary directory as vg.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::TempFile::setDirectory(temp_file::get_dir());

    // This is the data we are using.
    gbwt::GBWT compressed_index;
    gbwt::DynamicGBWT dynamic_index;
    index_type in_use = index_none;
    std::unique_ptr<PathHandleGraph> input_graph;
    bool graph_in_use = false;


    // Input GBWT construction.
    if (build != build_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Building input GBWTs (" << input_filenames.size() << " inputs)" << std::endl;
        }
        get_graph(input_graph, graph_in_use, xg_name, show_progress);
        if (build == build_vcf) {
            if (show_progress) {
                std::cerr << "Input type: VCF" << std::endl;
            }
            // Process each VCF contig corresponding to a non-alt path.
            std::vector<path_handle_t> path_handles;
            input_graph->for_each_path_handle([&](path_handle_t path_handle) {
                if (!Paths::is_alt(input_graph->get_path_name(path_handle))) {
                    path_handles.push_back(path_handle);
                }
            });
            #pragma omp parallel for schedule(static, 1)
            for (size_t i = 0; i < input_filenames.size(); i++) {
                if (show_progress) {
                    #pragma omp critical
                    {
                        std::cerr << "Indexing " << input_filenames[i] << std::endl;
                    }
                }
                if (!haplotype_indexer.batch_file_prefix.empty()) {
                    vcflib::VariantCallFile variant_file;
                    variant_file.parseSamples = false; // vcflib parsing is very slow if there are many samples.
                    variant_file.open(input_filenames[i]);
                    if (!variant_file.is_open()) {
                        cerr << "error: [vg gbwt] could not open VCF file " << input_filenames[i] << endl;
                        std::exit(EXIT_FAILURE);
                    } else if (show_progress) {
                        #pragma omp critical
                        {
                            std::cerr << "Opened variant file " << input_filenames[i] << std::endl;
                        }
                    }
                    // Run VCF parsing but do nothing with the generated phasing batches.
                    std::vector<std::string> sample_names;
                    haplotype_indexer.parse_vcf(input_graph.get(), path_handles, variant_file, sample_names,
                        [&](size_t contig, const gbwt::VariantPaths& variants, gbwt::PhasingInformation& phasings_batch) {},
                        false);
                } else {
                    std::unique_ptr<gbwt::DynamicGBWT> temp = haplotype_indexer.build_gbwt(input_graph.get(), input_filenames[i], false);
                    use_or_save(temp, dynamic_index, in_use, input_filenames, i, show_progress);
                }
            }
        } else if (build == build_paths) {
            if(show_progress) {
                std::cerr << "Input type: embedded paths" << std::endl;
            }
            std::unique_ptr<gbwt::DynamicGBWT> temp = haplotype_indexer.build_gbwt(input_graph.get());
            dynamic_index = std::move(*temp);
            in_use = index_dynamic;
        } else if (build == build_alignments) {
            if (show_progress) {
                std::cerr << "Input type: " << (gam_format ? "GAM" : "GAF") << std::endl;
            }
            #pragma omp parallel for schedule(static, 1)
            for (size_t i = 0; i < input_filenames.size(); i++) {
                if (show_progress) {
                    #pragma omp critical
                    {
                        std::cerr << "Indexing " << input_filenames[i] << std::endl;
                    }
                }
                std::vector<std::string> curr = { input_filenames[i] };
                std::unique_ptr<gbwt::DynamicGBWT> temp = haplotype_indexer.build_gbwt(input_graph.get(), curr, (gam_format ? "GAM" : "GAF"));
                use_or_save(temp, dynamic_index, in_use, input_filenames, i, show_progress);
            }
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWTs built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Merge multiple input GBWTs.
    if (merge != merge_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::string algo_name = (merge == merge_fast ? "fast" : "insertion");
            std::cerr << "Merging " << input_filenames.size() << " input GBWTs (" << algo_name << " algorithm)" << std::endl;
        }
        if (merge == merge_fast) {
            std::vector<gbwt::GBWT> indexes(input_filenames.size());
            for (size_t i = 0; i < input_filenames.size(); i++) {
                load_gbwt(input_filenames[i], indexes[i], show_progress);
            }
            if (show_progress) {
                std::cerr << "Merging the GBWTs" << std::endl;
            }
            compressed_index = gbwt::GBWT(indexes);
            in_use = index_compressed;
        }
        else if (merge == merge_insert) {
            load_gbwt(gbwt_name, dynamic_index, show_progress);
            in_use = index_dynamic;
            for (size_t i = 1; i < input_filenames.size(); i++) {
                gbwt::GBWT next;
                load_gbwt(input_filenames[i], next, show_progress);
                if (next.size() > 2 * dynamic_index.size()) {
                    std::cerr << "warning: [vg gbwt] merging " << input_filenames[i] << " into a substantially smaller index" << std::endl;
                    std::cerr << "warning: [vg gbwt] merging would be faster in another order" << std::endl;
                }
                if (show_progress) {
                    std::cerr << "Inserting " << next.sequences() << " sequences of total length " << next.size() << std::endl;
                }
                dynamic_index.merge(next);
            }
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWTs merged in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Remove a sample from the GBWT.
    if (!to_remove.empty()) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Removing sample " << to_remove << " from the index" << std::endl;
        }
        get_dynamic(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (!(dynamic_index.hasMetadata() && dynamic_index.metadata.hasPathNames() && dynamic_index.metadata.hasSampleNames())) {
            std::cerr << "error: [vg gbwt] the index does not contain metadata with thread and sample names" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        gbwt::size_type sample_id = dynamic_index.metadata.sample(to_remove);
        std::vector<gbwt::size_type> path_ids = dynamic_index.metadata.removeSample(sample_id);
        if (path_ids.empty()) {
            std::cerr << "error: [vg gbwt] the index does not contain sample " << to_remove << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (show_progress) {
            std::cerr << "Removing " << path_ids.size() << " threads" << std::endl;
        }
        size_t foo = dynamic_index.remove(path_ids);
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "Sample removed in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Path cover options.
    if (path_cover != path_cover_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Finding a " << num_paths << "-path cover with context length " << context_length << std::endl;
        }
        get_graph(input_graph, graph_in_use, xg_name, show_progress);
        if (path_cover == path_cover_greedy) {
            if (show_progress) {
                std::cerr << "Algorithm: greedy" << std::endl;
            }
            compressed_index = gbwtgraph::path_cover_gbwt(*input_graph, num_paths, context_length, haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, haplotype_indexer.id_interval, show_progress);
            in_use = index_compressed;
        } else if (path_cover == path_cover_augment) {
            if (show_progress) {
                std::cerr << "Algorithm: augment" << std::endl;
            }
            get_dynamic(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            gbwtgraph::augment_gbwt(*input_graph, dynamic_index, num_paths, context_length, haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, haplotype_indexer.id_interval, show_progress);
        } else {
            if (show_progress) {
                std::cerr << "Algorithm: local haplotypes" << std::endl;
            }
            get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            compressed_index = gbwtgraph::local_haplotypes(*input_graph, compressed_index, num_paths, context_length, haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, haplotype_indexer.id_interval, show_progress);
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "Path cover built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Now we can serialize the GBWT.
    if (!gbwt_output.empty() && in_use != index_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Serializing the GBWT to " << gbwt_output << std::endl;
        }
        if (in_use == index_compressed) {
            vg::io::VPKG::save(compressed_index, gbwt_output);
        } else if (in_use == index_dynamic) {
            vg::io::VPKG::save(dynamic_index, gbwt_output);
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWT serialized in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // GBWTGraph construction.
    if (!graph_output.empty()) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Building GBWTGraph" << std::endl;
        }
        get_graph(input_graph, graph_in_use, xg_name, show_progress);
        get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (show_progress) {
            std::cerr << "Starting the construction" << std::endl;
        }
        gbwtgraph::GBWTGraph graph(compressed_index, *input_graph);
        if (show_progress) {
            std::cerr << "Serializing GBWTGraph to " << graph_output << std::endl;
        }
        vg::io::VPKG::save(graph, graph_output);
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWTGraph built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // We no longer need the graph.
    clear_graph(input_graph, graph_in_use);


    // R-index construction.
    if (!r_index_name.empty()) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Building r-index" << std::endl;
        }
        get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (show_progress) {
            std::cerr << "Starting the construction" << std::endl;
        }
        gbwt::FastLocate r_index(compressed_index);
        if (show_progress) {
            std::cerr << "Serializing the r-index to " << r_index_name << std::endl;
        }
        vg::io::VPKG::save(r_index, r_index_name);
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "R-index built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Metadata options.
    if (metadata_mode) {
        get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
        if (!compressed_index.hasMetadata()) {
            std::cerr << "error: [vg gbwt] the GBWT does not contain metadata" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (metadata) {
            gbwt::operator<<(std::cout, compressed_index.metadata) << std::endl;
        }
        if (contigs) {
            if (list_names) {
                if (compressed_index.metadata.hasContigNames()) {
                    for (size_t i = 0; i < compressed_index.metadata.contigs(); i++) {
                        std::cout << compressed_index.metadata.contig(i) << std::endl;
                    }
                } else {
                    std::cerr << "error: [vg gbwt] the metadata does not contain contig names" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::cout << compressed_index.metadata.contigs() << std::endl;
            }
        }
        if (haplotypes) {
            std::cout << compressed_index.metadata.haplotypes() << std::endl;
        }
        if (samples) {
            if (list_names) {
                if (compressed_index.metadata.hasSampleNames()) {
                    for (size_t i = 0; i < compressed_index.metadata.samples(); i++) {
                        std::cout << compressed_index.metadata.sample(i) << std::endl;
                    }
                } else {
                    std::cerr << "error: [vg gbwt] the metadata does not contain sample names" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::cout << compressed_index.metadata.samples() << std::endl;
            }
        }
        if (thread_names) {
            if (compressed_index.metadata.hasPathNames()) {
                for (size_t i = 0; i < compressed_index.metadata.paths(); i++) {
                    std::cout << thread_name(compressed_index, i) << std::endl;
                }
            } else {
                std::cerr << "error: [vg gbwt] the metadata does not contain thread names" << std::endl;
            }
        }
    }

    // Thread options.
    if (thread_mode) {
        // Extract threads in SDSL format.
        if (!thread_output.empty()) {
            double start = gbwt::readTimer();
            if (show_progress) {
                std::cerr << "Extracting threads to " << thread_output << std::endl;
            }
            get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            if (show_progress) {
                std::cerr << "Starting the extraction" << std::endl;
            }
            gbwt::size_type node_width = gbwt::bit_length(compressed_index.sigma() - 1);
            gbwt::text_buffer_type out(thread_output, std::ios::out, gbwt::MEGABYTE, node_width);
            for (gbwt::size_type id = 0; id < compressed_index.sequences(); id += 2) { // Ignore reverse complements.
                gbwt::vector_type sequence = compressed_index.extract(id);
                for (auto node : sequence) {
                    out.push_back(node);
                }
                out.push_back(gbwt::ENDMARKER);
            }
            out.close();
            if (show_progress) {
                double seconds = gbwt::readTimer() - start;
                std::cerr << "Threads extracted in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
                std::cerr << std::endl;
            }
        }

        // There are two sequences for each thread.
        if (count_threads) {
            get_compressed(compressed_index, dynamic_index, in_use, gbwt_name, show_progress);
            std::cout << (compressed_index.sequences() / 2) << std::endl;
        }
    }

    return 0;
}


//----------------------------------------------------------------------------

// Utility functions

void load_gbwt(const std::string& filename, gbwt::GBWT& index, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading " << filename << std::endl;
    }
    std::unique_ptr<gbwt::GBWT> loaded = vg::io::VPKG::load_one<gbwt::GBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [vg gbwt] could not load compressed GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void load_gbwt(const std::string& filename, gbwt::DynamicGBWT& index, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading " << filename << std::endl;
    }
    std::unique_ptr<gbwt::DynamicGBWT> loaded = vg::io::VPKG::load_one<gbwt::DynamicGBWT>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [vg gbwt] could not load dynamic GBWT " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void get_compressed(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress) {
    if (in_use == index_compressed) {
        return;
    } else if (in_use == index_dynamic) {
        if (show_progress) {
            std::cerr << "Converting dynamic GBWT into compressed GBWT" << std::endl;
        }
        compressed_index = gbwt::GBWT(dynamic_index);
        dynamic_index = gbwt::DynamicGBWT();
        in_use = index_compressed;
    } else {
        load_gbwt(filename, compressed_index, show_progress);
        in_use = index_compressed;
    }
}

void get_dynamic(gbwt::GBWT& compressed_index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, const std::string& filename, bool show_progress) {
    if (in_use == index_dynamic) {
        return;
    } else if (in_use == index_compressed) {
        if (show_progress) {
            std::cerr << "Converting compressed GBWT into dynamic GBWT" << std::endl;
        }
        dynamic_index = gbwt::DynamicGBWT(compressed_index);
        compressed_index = gbwt::GBWT();
        in_use = index_dynamic;
    } else {
        load_gbwt(filename, dynamic_index, show_progress);
        in_use = index_dynamic;
    }
}

void use_or_save(std::unique_ptr<gbwt::DynamicGBWT>& index, gbwt::DynamicGBWT& dynamic_index, index_type& in_use, std::vector<std::string>& filenames, size_t i, bool show_progress) {
    if (filenames.size() == 1) {
        dynamic_index = std::move(*index);
        in_use = index_dynamic;
    } else {
        std::string temp = temp_file::create("gbwt");
        if (show_progress) {
            #pragma omp critical
            {
                std::cerr << "Saving the GBWT of " << filenames[i] << " to " << temp << std::endl;
            }
        }
        vg::io::VPKG::save(*index, temp);
        filenames[i] = temp;
    }
}

void get_graph(std::unique_ptr<PathHandleGraph>& graph, bool& in_use, const std::string& filename, bool show_progress) {
    if (in_use) {
        return;
    } else {
        if (show_progress) {
            std::cerr << "Loading input graph from " << filename << std::endl;
        }
        graph = vg::io::VPKG::load_one<PathHandleGraph>(filename);
        if (graph == nullptr) {
            std::cerr << "error: [vg gbwt] could not load graph " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        in_use = true;
    }
}

void clear_graph(std::unique_ptr<PathHandleGraph>& graph, bool& in_use) {
    if (!in_use) {
        return;
    } else {
        graph.reset();
        in_use = false;
    }
}

//----------------------------------------------------------------------------


// Register subcommand
static Subcommand vg_gbwt("gbwt", "Manipulate GBWTs", main_gbwt);

