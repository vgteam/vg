/** \file gbwt_main.cpp
 *
 * Defines the "vg gbwt" subcommand for building, merging, and manipulating GBWT indexes
 * and GBWTGraphs.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <set>

#include "subcommand.hpp"
#include "../gbwt_helper.hpp"
#include "../haplotype_indexer.hpp"
#include "../path.hpp"
#include "../region.hpp"

#include <vg/io/vpkg.hpp>

#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/path_cover.h>

using namespace vg;

enum build_mode { build_none, build_vcf, build_paths, build_alignments };
enum merge_mode { merge_none, merge_insert, merge_fast, merge_parallel };
enum path_cover_mode { path_cover_none, path_cover_augment, path_cover_local, path_cover_greedy };

void use_or_save(std::unique_ptr<gbwt::DynamicGBWT>& index, GBWTHandler& gbwts, std::vector<std::string>& filenames, size_t i, bool show_progress);

void print_metadata(std::ostream& out, const GBWTHandler& gbwts);

void get_graph(std::unique_ptr<PathHandleGraph>& graph, bool& in_use, const std::string& filename, bool show_progress);
void clear_graph(std::unique_ptr<PathHandleGraph>& graph, bool& in_use);

struct job_type {
    std::string filename;
    std::vector<path_handle_t> paths;
    size_t size;

    // Large jobs first.
    bool operator<(const job_type& another) const {
        return (this->size > another.size);
    }

    typedef std::pair<path_handle_t, size_t> path_type;

    void insert(path_type path) {
        this->paths.push_back(path.first);
        this->size += path.second;
    }
};

std::vector<job_type> determine_jobs(const std::vector<std::string>& vcf_files, std::unique_ptr<PathHandleGraph>& graph, const std::map<std::string, std::string>& path_to_vcf, bool inputs_as_jobs);

size_t default_build_jobs() {
    return std::max(static_cast<size_t>(1), static_cast<size_t>(omp_get_max_threads() / 2));
}

size_t default_merge_jobs() {
    return std::min(static_cast<size_t>(gbwt::MergeParameters::MERGE_JOBS), std::max(static_cast<size_t>(1), static_cast<size_t>(omp_get_max_threads() / 2)));
}

void use_preset(HaplotypeIndexer& haplotype_indexer, std::string preset_name);

void help_gbwt(char** argv) {
    std::cerr << "usage: " << argv[0] << " gbwt [options] [args]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Manipulate GBWTs. Input GBWTs are loaded from input args or built in earlier steps." << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "    -x, --xg-name FILE      read the graph from FILE" << std::endl;
    std::cerr << "    -o, --output FILE       write output GBWT to FILE" << std::endl;
    std::cerr << "    -d, --temp-dir DIR      use directory DIR for temporary files" << std::endl;
    std::cerr << "    -p, --progress          show progress and statistics" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GBWT construction parameters (for steps 1 and 4):" << std::endl;
    std::cerr << "        --buffer-size N     GBWT construction buffer size in millions of nodes (default " << (gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION) << ")" << std::endl;
    std::cerr << "        --id-interval N     store path ids at one out of N positions (default " << gbwt::DynamicGBWT::SAMPLE_INTERVAL << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Search parameters (for -b and -r):" << std::endl;
    std::cerr << "        --num-threads N     use N parallel search threads (default " << omp_get_max_threads() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 1: GBWT construction (requires -o, -x, and one of { -v, -E, A }):" << std::endl;
    std::cerr << "    -v, --vcf-input         index the haplotypes in the VCF files specified in input args in parallel" << std::endl;
    std::cerr << "                            (inputs must be over different contigs; implies -f)" << std::endl;
    std::cerr << "        --preset X          use preset X (available: 1000gp)" << std::endl;
    std::cerr << "        --num-jobs N        use at most N parallel build jobs (default " << default_build_jobs() << ")" << std::endl;
    std::cerr << "        --inputs-as-jobs    create one build job for each input instead of using first-fit heuristic" << std::endl;
    std::cerr << "        --parse-only        store the VCF parses without building GBWTs" << std::endl;
    std::cerr << "                            (use -o for the file name prefix; skips subsequent steps)" << std::endl;
    std::cerr << "        --ignore-missing    do not warn when variants are missing from the graph" << std::endl;
    std::cerr << "        --actual-phasing    do not interpret unphased homozygous genotypes as phased" << std::endl;
    std::cerr << "        --force-phasing     replace unphased genotypes with randomly phased ones" << std::endl;
    std::cerr << "        --discard-overlaps  skip overlapping alternate alleles if the overlap cannot be resolved" << std::endl;
    std::cerr << "                            instead of creating a phase break" << std::endl;
    std::cerr << "        --batch-size N      index the haplotypes in batches of N samples (default 200)" << std::endl; // FIXME source for the default
    std::cerr << "        --sample-range X-Y  index samples X to Y (inclusive, 0-based)" << std::endl;
    std::cerr << "        --rename V=P        VCF contig V matches path P in the graph (may repeat)" << std::endl;
    std::cerr << "        --vcf-variants      variants in the graph use VCF contig names instead of path names" << std::endl;
    std::cerr << "        --vcf-region C:X-Y  restrict VCF contig C to coordinates X to Y (inclusive, 1-based; may repeat)" << std::endl;
    std::cerr << "        --exclude-sample X  do not index the sample with name X (faster than -R; may repeat)" << std::endl;
    std::cerr << "    -E, --index-paths       index the embedded non-alt paths in the graph (no input args)" << std::endl;
    std::cerr << "        --paths-as-samples  each path becomes a sample instead of a contig in the metadata" << std::endl;
    std::cerr << "    -A, --alignment-input   index the alignments in the GAF files specified in input args" << std::endl;
    std::cerr << "        --gam-format        the input files are in GAM format instead of GAF format" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 2: Merge multiple input GBWTs (requires -o):" << std::endl;
    std::cerr << "    -m, --merge             use the insertion algorithm" << std::endl;
    std::cerr << "    -f, --fast              fast merging algorithm (node ids must not overlap)" << std::endl;
    std::cerr << "    -b, --parallel          use the parallel algorithm" << std::endl;
    std::cerr << "        --chunk-size N      search in chunks of N sequences (default " << gbwt::MergeParameters::CHUNK_SIZE << ")" << std::endl;
    std::cerr << "        --pos-buffer N      use N MiB position buffers for each search thread (default " << gbwt::MergeParameters::POS_BUFFER_SIZE << ")" << std::endl;
    std::cerr << "        --thread-buffer N   use N MiB thread buffers for each search thread (default " << gbwt::MergeParameters::THREAD_BUFFER_SIZE << ")" << std::endl;
    std::cerr << "        --merge-buffers N   merge 2^N thread buffers into one file per merge job (default " << gbwt::MergeParameters::MERGE_BUFFERS << ")" << std::endl;
    std::cerr << "        --merge-jobs N      run N parallel merge jobs (default " << default_merge_jobs() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 3: Remove samples (requires -o and one input GBWT):" << std::endl;
    std::cerr << "    -R, --remove-sample X   remove the sample with name X from the index (may repeat)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 4: Path cover GBWT construction (requires -o, -x, and one of { -a, -l, -P }):" << std::endl;
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
    std::cerr << "Step 7: Metadata (one input GBWT):" << std::endl;
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

    // Long options with no corresponding short options.
    constexpr int OPT_BUFFER_SIZE = 1000;
    constexpr int OPT_ID_INTERVAL = 1001;
    constexpr int OPT_NUM_THREADS = 1002;
    constexpr int OPT_PRESET = 1100;
    constexpr int OPT_NUM_JOBS = 1101;
    constexpr int OPT_INPUTS_AS_JOBS = 1102;
    constexpr int OPT_PARSE_ONLY = 1103;
    constexpr int OPT_IGNORE_MISSING = 1104;
    constexpr int OPT_ACTUAL_PHASING = 1105;
    constexpr int OPT_FORCE_PHASING = 1106;
    constexpr int OPT_DISCARD_OVERLAPS = 1107;
    constexpr int OPT_BATCH_SIZE = 1108;
    constexpr int OPT_SAMPLE_RANGE = 1109;
    constexpr int OPT_RENAME = 1110;
    constexpr int OPT_VCF_VARIANTS = 1111;
    constexpr int OPT_VCF_REGION = 1112;
    constexpr int OPT_EXCLUDE_SAMPLE = 1113;
    constexpr int OPT_PATHS_AS_SAMPLES = 1114;
    constexpr int OPT_GAM_FORMAT = 1115;
    constexpr int OPT_CHUNK_SIZE = 1200;
    constexpr int OPT_POS_BUFFER = 1201;
    constexpr int OPT_THREAD_BUFFER = 1202;
    constexpr int OPT_MERGE_BUFFERS = 1203;
    constexpr int OPT_MERGE_JOBS = 1204;

    // Requirements and modes.
    bool produces_one_gbwt = false;
    build_mode build = build_none;
    merge_mode merge = merge_none;
    path_cover_mode path_cover = path_cover_none;
    bool metadata_mode = false, thread_mode = false;

    // Input GBWT construction.
    HaplotypeIndexer haplotype_indexer;
    bool gam_format = false, inputs_as_jobs = false, parse_only = false;
    size_t build_jobs = default_build_jobs();

    // Parallel merging.
    gbwt::MergeParameters merge_parameters;
    merge_parameters.setMergeJobs(default_merge_jobs());

    // Other parameters and flags.
    bool show_progress = false;
    bool count_threads = false;
    bool metadata = false, contigs = false, haplotypes = false, samples = false, list_names = false, thread_names = false;
    size_t num_paths = gbwtgraph::PATH_COVER_DEFAULT_N, context_length = gbwtgraph::PATH_COVER_DEFAULT_K;
    bool num_paths_set = false;
    size_t search_threads = omp_get_max_threads();

    // File/sample names.
    std::string gbwt_output, thread_output;
    std::string graph_output, xg_name;
    std::string r_index_name;
    std::set<std::string> to_remove;

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

                // GBWT construction parameters
                { "buffer-size", required_argument, 0, OPT_BUFFER_SIZE },
                { "id-interval", required_argument, 0, OPT_ID_INTERVAL },

                // Search parameters
                { "num-threads", required_argument, 0, OPT_NUM_THREADS },

                // Input GBWT construction
                { "vcf-input", no_argument, 0, 'v' },
                { "preset", required_argument, 0, OPT_PRESET },
                { "num-jobs", required_argument, 0, OPT_NUM_JOBS },
                { "inputs-as-jobs", no_argument, 0, OPT_INPUTS_AS_JOBS },
                { "parse-only", no_argument, 0, OPT_PARSE_ONLY },
                { "ignore-missing", no_argument, 0, OPT_IGNORE_MISSING },
                { "actual-phasing", no_argument, 0, OPT_ACTUAL_PHASING },
                { "force-phasing", no_argument, 0, OPT_FORCE_PHASING },
                { "discard-overlaps", no_argument, 0, OPT_DISCARD_OVERLAPS },
                { "batch-size", required_argument, 0, OPT_BATCH_SIZE },
                { "sample-range", required_argument, 0, OPT_SAMPLE_RANGE },
                { "rename", required_argument, 0, OPT_RENAME },
                { "vcf-variants", no_argument, 0, OPT_VCF_VARIANTS },
                { "vcf-region", required_argument, 0, OPT_VCF_REGION },
                { "exclude-sample", required_argument, 0, OPT_EXCLUDE_SAMPLE },
                { "index-paths", no_argument, 0, 'E' },
                { "paths-as-samples", no_argument, 0, OPT_PATHS_AS_SAMPLES },
                { "alignment-input", no_argument, 0, 'A' },
                { "gam-format", no_argument, 0, OPT_GAM_FORMAT },

                // Merging
                { "merge", no_argument, 0, 'm' },
                { "fast", no_argument, 0, 'f' },
                { "parallel", no_argument, 0, 'b' },
                { "chunk-size", required_argument, 0, OPT_CHUNK_SIZE },
                { "pos-buffer", required_argument, 0, OPT_POS_BUFFER },
                { "thread-buffer", required_argument, 0, OPT_THREAD_BUFFER },
                { "merge-buffers", required_argument, 0, OPT_MERGE_BUFFERS },
                { "merge-jobs", required_argument, 0, OPT_MERGE_JOBS },

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
        c = getopt_long(argc, argv, "x:o:d:pvEAmfbR:alPn:k:g:r:MCHSLTce:h?", long_options, &option_index);

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
            haplotype_indexer.show_progress = true;
            break;

        // GBWT construction parameters
        case OPT_BUFFER_SIZE:
            haplotype_indexer.gbwt_buffer_size = std::max(parse<size_t>(optarg), 1ul);
            break;
        case OPT_ID_INTERVAL:
            haplotype_indexer.id_interval = parse<size_t>(optarg);
            break;

        // Search parameters
        case OPT_NUM_THREADS:
            search_threads = std::max(parse<size_t>(optarg), 1ul);
            break;

        // Input GBWT construction
        case 'v':
            assert(build == build_none);
            build = build_vcf;
            break;
        case OPT_PRESET:
            use_preset(haplotype_indexer, optarg);
            break;
        case OPT_NUM_JOBS:
            build_jobs = parse<size_t>(optarg);
            break;
        case OPT_INPUTS_AS_JOBS:
            inputs_as_jobs = true;
            break;
        case OPT_PARSE_ONLY:
            parse_only = true;
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
        case OPT_VCF_VARIANTS:
            haplotype_indexer.rename_variants = false;
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
            produces_one_gbwt = true;
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
        case 'b':
            merge = merge_parallel;
            produces_one_gbwt = true;
            break;
        case OPT_CHUNK_SIZE:
            merge_parameters.setChunkSize(parse<size_t>(optarg));
            break;
        case OPT_POS_BUFFER:
            merge_parameters.setPosBufferSize(parse<size_t>(optarg));
            break;
        case OPT_THREAD_BUFFER:
            merge_parameters.setThreadBufferSize(parse<size_t>(optarg));
            break;
        case OPT_MERGE_BUFFERS:
            merge_parameters.setMergeBuffers(parse<size_t>(optarg));
            break;
        case OPT_MERGE_JOBS:
            merge_parameters.setMergeJobs(parse<size_t>(optarg));
            break;

        // Remove sample
        case 'R':
            to_remove.insert(optarg);
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
        if (xg_name.empty() || gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt]: GBWT construction requires a graph and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (build == build_vcf) {
            if (input_filenames.empty()) {
                std::cerr << "error: [vg gbwt]: GBWT construction from VCF files requires input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (parse_only) {
                haplotype_indexer.batch_file_prefix = gbwt_output;
            }
        } else if (build == build_alignments) {
            if (input_filenames.empty()) {
                std::cerr << "error: [vg gbwt]: GBWT construction from alignments requires input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        } else if (build == build_paths) {
            if (!input_filenames.empty()) {
                std::cerr << "error: [vg gbwt]: GBWT construction from embedded paths does not use input args" << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
    if (merge != merge_none) {
        if (input_filenames.size() < 2 || gbwt_output.empty()) {
            std::cerr << "error: [vg gbwt]: merging requires multiple input GBWTs and output GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (!to_remove.empty()) {
        // Ugly hack: We cannot use produces_one_gbwt, because the GBWT may be produced in the next step.
        if (!(input_filenames.size() == 1 || merge != merge_none) || gbwt_output.empty()) {
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
            std::cerr << "error: [vg gbwt]: metadata operations require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    if (thread_mode) {
        if (!(input_filenames.size() == 1 || produces_one_gbwt)) {
            std::cerr << "error: [vg gbwt]: thread operations require one input GBWT" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }


    // Let GBWT operate silently and use the same temporary directory as vg.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    gbwt::TempFile::setDirectory(temp_file::get_dir());

    // This is the data we are using.
    GBWTHandler gbwts;
    gbwts.filename = gbwt_name;
    gbwts.show_progress = show_progress;
    std::unique_ptr<PathHandleGraph> input_graph;
    bool graph_in_use = false;


    // Input GBWT construction.
    if (build != build_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Building input GBWTs" << std::endl;
        }
        get_graph(input_graph, graph_in_use, xg_name, show_progress);
        if (build == build_vcf) {
            if (show_progress) {
                std::cerr << "Input type: VCF" << std::endl;
            }
            omp_set_num_threads(build_jobs);
            // Process each VCF contig corresponding to a non-alt path.
            std::vector<job_type> jobs = determine_jobs(input_filenames, input_graph, haplotype_indexer.path_to_vcf, inputs_as_jobs);
            if (jobs.size() > 1 && merge == merge_none) {
                merge = merge_fast;
            }
            std::vector<std::vector<std::string>> vcf_parses(jobs.size());
            if (show_progress) {
                std::cerr << "Parsing " << jobs.size() << " VCF files using up to " << build_jobs << " parallel jobs" << std::endl;
            }
            #pragma omp parallel for schedule(dynamic, 1)
            for (size_t i = 0; i < jobs.size(); i++) {
                std::string job_name = "Job " + std::to_string(i);
                if (show_progress) {
                    #pragma omp critical
                    {
                        std::cerr << job_name << ": File " << jobs[i].filename << ", paths {";
                        for (path_handle_t handle : jobs[i].paths) {
                            std::cerr << " " << input_graph->get_path_name(handle);
                        }
                        std::cerr << " }" << std::endl;
                    }
                }
                vcf_parses[i] = haplotype_indexer.parse_vcf(jobs[i].filename, *input_graph, jobs[i].paths, job_name);
            }
            // Delete the graph to save memory.
            clear_graph(input_graph, graph_in_use);
            if (!parse_only) {
                std::vector<std::string> gbwt_files(vcf_parses.size(), "");
                if (show_progress) {
                    std::cerr << "Building " << vcf_parses.size() << " GBWTs using up to " << build_jobs << " parallel jobs" << std::endl;
                }
                #pragma omp parallel for schedule(dynamic, 1)
                for (size_t i = 0; i < vcf_parses.size(); i++) {
                    std::string job_name = "Job " + std::to_string(i);
                    std::unique_ptr<gbwt::DynamicGBWT> parsed = haplotype_indexer.build_gbwt(vcf_parses[i], job_name);
                    use_or_save(parsed, gbwts, gbwt_files, i, show_progress);
                }
                if (vcf_parses.size() > 1) {
                    input_filenames = gbwt_files; // Use the temporary GBWTs as inputs.
                }
            }
        } else if (build == build_paths) {
            if(show_progress) {
                std::cerr << "Input type: embedded paths" << std::endl;
            }
            std::unique_ptr<gbwt::DynamicGBWT> temp = haplotype_indexer.build_gbwt(*input_graph);
            gbwts.use(*temp);
        } else if (build == build_alignments) {
            if (show_progress) {
                std::cerr << "Input type: " << (gam_format ? "GAM" : "GAF") << std::endl;
            }
            std::unique_ptr<gbwt::DynamicGBWT> temp = haplotype_indexer.build_gbwt(*input_graph, input_filenames, (gam_format ? "GAM" : "GAF"));
            gbwts.use(*temp);
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWTs built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
        if (parse_only) {
            return 0; // VCF parsing does not produce GBWTs to continue with.
        }
    }


    // Merge multiple input GBWTs.
    if (merge != merge_none) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::string algo_name;
            if (merge == merge_fast) {
                algo_name = "fast";
            } else if (merge == merge_insert) {
                algo_name = "insertion";
            } else if (merge == merge_parallel) {
                algo_name = "parallel";
            }
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
            gbwt::GBWT merged(indexes);
            gbwts.use(merged);
        } else if (merge == merge_insert) {
            gbwts.use_dynamic();
            for (size_t i = 1; i < input_filenames.size(); i++) {
                gbwt::GBWT next;
                load_gbwt(input_filenames[i], next, show_progress);
                if (next.size() > 2 * gbwts.dynamic.size()) {
                    std::cerr << "warning: [vg gbwt] merging " << input_filenames[i] << " into a substantially smaller index" << std::endl;
                    std::cerr << "warning: [vg gbwt] merging would be faster in another order" << std::endl;
                }
                if (show_progress) {
                    std::cerr << "Inserting " << next.sequences() << " sequences of total length " << next.size() << std::endl;
                }
                gbwts.dynamic.merge(next);
            }
        } else if (merge == merge_parallel) {
            gbwts.use_dynamic();
            omp_set_num_threads(search_threads);
            for (size_t i = 1; i < input_filenames.size(); i++) {
                gbwt::DynamicGBWT next;
                load_gbwt(input_filenames[i], next, show_progress);
                if (next.size() > 2 * gbwts.dynamic.size()) {
                    std::cerr << "warning: [vg gbwt] merging " << input_filenames[i] << " into a substantially smaller index" << std::endl;
                    std::cerr << "warning: [vg gbwt] merging would be faster in another order" << std::endl;
                }
                if (show_progress) {
                    std::cerr << "Inserting " << next.sequences() << " sequences of total length " << next.size() << std::endl;
                }
                gbwts.dynamic.merge(next, merge_parameters);
            }
        }
        if (show_progress) {
            print_metadata(std::cerr, gbwts);
            double seconds = gbwt::readTimer() - start;
            std::cerr << "GBWTs merged in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Remove samples from the GBWT.
    if (!to_remove.empty()) {
        double start = gbwt::readTimer();
        if (show_progress) {
            std::cerr << "Removing " << to_remove.size() << " sample(s) from the index" << std::endl;
        }
        gbwts.use_dynamic();
        if (!(gbwts.dynamic.hasMetadata() && gbwts.dynamic.metadata.hasPathNames() && gbwts.dynamic.metadata.hasSampleNames())) {
            std::cerr << "error: [vg gbwt] the index does not contain metadata with thread and sample names" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::set<gbwt::size_type> sample_ids;
        for (const std::string& sample_name : to_remove) {
            gbwt::size_type sample_id = gbwts.dynamic.metadata.sample(sample_name);
            if (sample_id >= gbwts.dynamic.metadata.samples()) {
                std::cerr << "warning: [vg gbwt] the index does not contain sample " << sample_name << std::endl;
            } else {
                sample_ids.insert(sample_id);
            }
        }
        std::vector<gbwt::size_type> path_ids;
        for (gbwt::size_type sample_id : sample_ids) {
            std::vector<gbwt::size_type> current_paths = gbwts.dynamic.metadata.removeSample(sample_id);
            path_ids.insert(path_ids.end(), current_paths.begin(), current_paths.end());
        }
        if (path_ids.empty()) {
            std::cerr << "warning: [vg gbwt] no threads associated with the removed samples" << std::endl;
        } else {
            if (show_progress) {
                std::cerr << "Removing " << path_ids.size() << " threads" << std::endl;
            }
            size_t foo = gbwts.dynamic.remove(path_ids);
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "Samples removed in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
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
            gbwt::GBWT cover = gbwtgraph::path_cover_gbwt(*input_graph, num_paths, context_length, haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, haplotype_indexer.id_interval, show_progress);
            gbwts.use(cover);
        } else if (path_cover == path_cover_augment) {
            if (show_progress) {
                std::cerr << "Algorithm: augment" << std::endl;
            }
            gbwts.use_dynamic();
            gbwtgraph::augment_gbwt(*input_graph, gbwts.dynamic, num_paths, context_length, haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, haplotype_indexer.id_interval, show_progress);
        } else {
            if (show_progress) {
                std::cerr << "Algorithm: local haplotypes" << std::endl;
            }
            gbwts.use_compressed();
            gbwt::GBWT cover = gbwtgraph::local_haplotypes(*input_graph, gbwts.compressed, num_paths, context_length, haplotype_indexer.gbwt_buffer_size * gbwt::MILLION, haplotype_indexer.id_interval, show_progress);
            gbwts.use(cover);
        }
        if (show_progress) {
            double seconds = gbwt::readTimer() - start;
            std::cerr << "Path cover built in " << seconds << " seconds, " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
            std::cerr << std::endl;
        }
    }


    // Now we can serialize the GBWT.
    if (!gbwt_output.empty()) {
        double start = gbwt::readTimer();
        gbwts.serialize(gbwt_output);
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
        gbwts.use_compressed();
        if (show_progress) {
            std::cerr << "Starting the construction" << std::endl;
        }
        gbwtgraph::GBWTGraph graph(gbwts.compressed, *input_graph);
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
        omp_set_num_threads(search_threads);
        gbwts.use_compressed();
        if (show_progress) {
            std::cerr << "Starting the construction" << std::endl;
        }
        gbwt::FastLocate r_index(gbwts.compressed);
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
        gbwts.use_compressed();
        if (!gbwts.compressed.hasMetadata()) {
            std::cerr << "error: [vg gbwt] the GBWT does not contain metadata" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (metadata) {
            print_metadata(std::cout, gbwts);
        }
        if (contigs) {
            if (list_names) {
                if (gbwts.compressed.metadata.hasContigNames()) {
                    for (size_t i = 0; i < gbwts.compressed.metadata.contigs(); i++) {
                        std::cout << gbwts.compressed.metadata.contig(i) << std::endl;
                    }
                } else {
                    std::cerr << "error: [vg gbwt] the metadata does not contain contig names" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::cout << gbwts.compressed.metadata.contigs() << std::endl;
            }
        }
        if (haplotypes) {
            std::cout << gbwts.compressed.metadata.haplotypes() << std::endl;
        }
        if (samples) {
            if (list_names) {
                if (gbwts.compressed.metadata.hasSampleNames()) {
                    for (size_t i = 0; i < gbwts.compressed.metadata.samples(); i++) {
                        std::cout << gbwts.compressed.metadata.sample(i) << std::endl;
                    }
                } else {
                    std::cerr << "error: [vg gbwt] the metadata does not contain sample names" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::cout << gbwts.compressed.metadata.samples() << std::endl;
            }
        }
        if (thread_names) {
            if (gbwts.compressed.metadata.hasPathNames()) {
                for (size_t i = 0; i < gbwts.compressed.metadata.paths(); i++) {
                    std::cout << thread_name(gbwts.compressed, i) << std::endl;
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
            gbwts.use_compressed();
            if (show_progress) {
                std::cerr << "Starting the extraction" << std::endl;
            }
            gbwt::size_type node_width = gbwt::bit_length(gbwts.compressed.sigma() - 1);
            gbwt::text_buffer_type out(thread_output, std::ios::out, gbwt::MEGABYTE, node_width);
            for (gbwt::size_type id = 0; id < gbwts.compressed.sequences(); id += 2) { // Ignore reverse complements.
                gbwt::vector_type sequence = gbwts.compressed.extract(id);
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
            gbwts.use_compressed();
            std::cout << (gbwts.compressed.sequences() / 2) << std::endl;
        }
    }

    return 0;
}


//----------------------------------------------------------------------------

// Utility functions

void print_metadata(std::ostream& out, const GBWTHandler& gbwts) {
    if (gbwts.in_use == GBWTHandler::index_compressed) {
        gbwt::operator<<(out, gbwts.compressed.metadata) << std::endl;
    } else if (gbwts.in_use == GBWTHandler::index_dynamic) {
        gbwt::operator<<(out, gbwts.dynamic.metadata) << std::endl;
    }
}

void use_or_save(std::unique_ptr<gbwt::DynamicGBWT>& index, GBWTHandler& gbwts, std::vector<std::string>& filenames, size_t i, bool show_progress) {
    if (filenames.size() == 1) {
        gbwts.use(*index);
    } else {
        std::string temp = temp_file::create("gbwt-" + std::to_string(i) + "-");
        if (show_progress) {
            #pragma omp critical
            {
                std::cerr << "Job " << i << ": Saving the GBWT to " << temp << std::endl;
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

std::vector<job_type> determine_jobs(const std::vector<std::string>& vcf_files, std::unique_ptr<PathHandleGraph>& graph, const std::map<std::string, std::string>& path_to_vcf, bool inputs_as_jobs) {

    std::vector<job_type> result;

    // Determine the non-alt paths.
    std::vector<job_type::path_type> paths;
    size_t max_length = 0;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
        if (!Paths::is_alt(graph->get_path_name(path_handle))) {
            paths.emplace_back(path_handle, graph->get_step_count(path_handle));
            max_length = std::max(max_length, paths.back().second);
        }
    });
    if (paths.empty()) {
        return result;
    }

    struct vcf_paths {
        size_t file;
        std::vector<job_type::path_type> paths;

        // In descending order by length.
        void sort_paths() {
            std::sort(this->paths.begin(), this->paths.end(), [](job_type::path_type a, job_type::path_type b) -> bool {
                return (a.second > b.second);
            });
        }

        void insert(job_type::path_type path) {
            this->paths.push_back(path);
        }

        bool empty() const {
            return this->paths.empty();
        }
    };

    // Initialize the files.
    std::vector<vcf_paths> paths_by_file;
    for (size_t i = 0; i < vcf_files.size(); i++) {
        paths_by_file.push_back({ i, {} });
    }

    // Determine which VCF file contains each path.
    std::vector<size_t> path_found_in(paths.size(), vcf_files.size());
    for (size_t i = 0; i < vcf_files.size(); i++) {
        std::string filename = vcf_files[i];
        vcflib::VariantCallFile variant_file;
        variant_file.parseSamples = false;
        variant_file.open(filename);
        if (!variant_file.is_open()) {
            std::cerr << "error: [vg gbwt] could not open VCF file " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        for (size_t j = 0; j < paths.size(); j++) {
            std::string contig_name = graph->get_path_name(paths[j].first);
            if (path_to_vcf.find(contig_name) != path_to_vcf.end()) {
                contig_name = path_to_vcf.at(contig_name);
            }
            variant_file.setRegion(contig_name);
            vcflib::Variant var(variant_file);
            if (!(variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == contig_name)) {
                continue;
            }
            if (path_found_in[j] < vcf_files.size()) {
                std::cerr << "error: [vg gbwt] contig " << contig_name << " found in files " << vcf_files[path_found_in[j]] << " and " << filename << std::endl;
                std::exit(EXIT_FAILURE);
            }
            paths_by_file[i].insert(paths[j]);
            path_found_in[j] = i;
        }
    }

    // Special case: Each input file is a single job.
    if (inputs_as_jobs) {
        for (const vcf_paths& curr : paths_by_file) {
            if (curr.empty()) {
                continue;
            }
            job_type job({ vcf_files[curr.file], {}, 0 });
            for (auto path : curr.paths) {
                job.insert(path);
            }
            result.push_back(job);
        }
        return result;
    }

    // First-fit heuristic: Create jobs of size at most max_length from each file.
    for (size_t i = 0; i < paths_by_file.size(); i++) {
        paths_by_file[i].sort_paths();
    }
    for (const vcf_paths& curr : paths_by_file) {
        if (curr.empty()) {
            continue;
        }
        std::vector<job_type> jobs;
        for (auto path : curr.paths) {
            bool inserted = false;
            for (size_t i = 0; i < jobs.size(); i++) {
                if (jobs[i].size + path.second <= max_length) {
                    jobs[i].insert(path);
                    inserted = true;
                    break;
                }
            }
            if (!inserted) {
                jobs.push_back({ vcf_files[curr.file], {}, 0 });
                jobs.back().insert(path);
            }
        }
        result.insert(result.end(), jobs.begin(), jobs.end());
    }

    // Sort the jobs in descending order by size.
    std::sort(result.begin(), result.end());
    return result;
}

//----------------------------------------------------------------------------

void use_preset(HaplotypeIndexer& haplotype_indexer, std::string preset_name) {
    for (char& c : preset_name) {
        c = std::tolower(c);
    }
    if (preset_name == "1000gp") {
        haplotype_indexer.gbwt_buffer_size = 200;
        haplotype_indexer.samples_in_batch = 100;
        haplotype_indexer.force_phasing = true;
        haplotype_indexer.discard_overlaps = true;
    } else {
        std::cerr << "error : [vg gbwt] invalid preset: " << preset_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//----------------------------------------------------------------------------

// Register subcommand
static vg::subcommand::Subcommand vg_gbwt("gbwt", "Manipulate GBWTs", main_gbwt);

