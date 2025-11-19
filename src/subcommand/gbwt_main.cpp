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
#include "../gbwtgraph_helper.hpp"
#include "../haplotype_indexer.hpp"
#include "../path.hpp"
#include "../region.hpp"
#include "../algorithms/find_translation.hpp"

#include <vg/io/vpkg.hpp>

#include <gbwt/fast_locate.h>
#include <gbwtgraph/gfa.h>
#include <gbwtgraph/path_cover.h>

using namespace vg;

struct GBWTConfig {
    // Build mode also defines the type of input args.
    enum build_mode { build_none, build_vcf, build_gfa, build_paths, build_alignments, build_gbz };
    enum merge_mode { merge_none, merge_insert, merge_fast, merge_parallel };
    enum path_cover_mode { path_cover_none, path_cover_augment, path_cover_local, path_cover_greedy };

    // Logger to use when reporting errors/progress
    Logger logger = Logger("vg gbwt");

    // Requirements and modes.
    bool produces_one_gbwt = false; // Steps 1-4 eventually produce one input GBWT regardless of the number of input args.
    build_mode build = build_none;
    merge_mode merge = merge_none;
    path_cover_mode path_cover = path_cover_none;
    bool metadata_mode = false, path_mode = false;

    // Input GBWT construction.
    HaplotypeIndexer haplotype_indexer;
    bool gam_format = false, inputs_as_jobs = false, parse_only = false;
    size_t build_jobs = default_build_jobs();

    // GFA parsing.
    gbwtgraph::GFAParsingParameters gfa_parameters = get_best_gbwtgraph_gfa_parsing_parameters();

    // Parallel merging.
    gbwt::MergeParameters merge_parameters;

    // GBZ construction.
    bool gbz_format = true;

    // Other parameters and flags.
    bool show_progress = false;
    bool count_paths = false;
    bool metadata = false, contigs = false, haplotypes = false, samples = false, list_names = false, path_names = false, tags = false;
    bool include_named_paths = false;
    size_t num_paths = default_num_paths(), context_length = default_context_length();
    bool num_paths_set = false;
    size_t search_threads = omp_get_max_threads();

    // Input data.
    std::vector<std::string> input_filenames;
    std::string gbwt_name; // There is a single input GBWT to load.
    std::string graph_name;
    std::string gbwtgraph_name;

    // File names.
    std::string gbwt_output; // Output GBWT.
    std::string path_output; // Paths in SDSL format.
    std::string graph_output; // Output GBZ.
    std::string segment_translation; // Segment to node translation output.
    std::string r_index_name; // Output r-index.

    // Sample names and metadata
    std::set<std::string> to_remove; // Sample names to remove.
    std::map<std::string, std::string> tags_to_set; // Tag changes to apply to the GBWT
    
    GBWTConfig() {
        this->merge_parameters.setMergeJobs(default_merge_jobs());
    }

    static size_t default_build_jobs() {
        return std::max(static_cast<size_t>(1), static_cast<size_t>(omp_get_max_threads() / 2));
    }

    static constexpr size_t default_num_paths() {
        return gbwtgraph::PATH_COVER_DEFAULT_N;
    }

    static constexpr size_t default_num_paths_local() {
        return gbwtgraph::LOCAL_HAPLOTYPES_DEFAULT_N;
    }

    static constexpr size_t default_context_length() {
        return gbwtgraph::PATH_COVER_DEFAULT_K;
    }

    static size_t default_merge_jobs() {
        return std::min(static_cast<size_t>(gbwt::MergeParameters::MERGE_JOBS), std::max(static_cast<size_t>(1), static_cast<size_t>(omp_get_max_threads() / 2)));
    }

    gbwtgraph::PathCoverParameters path_cover_parameters() const {
        gbwtgraph::PathCoverParameters parameters;
        parameters.num_paths = this->num_paths;
        parameters.context = this->context_length;
        parameters.batch_size = this->haplotype_indexer.gbwt_buffer_size * gbwt::MILLION;
        parameters.sample_interval = this->haplotype_indexer.id_interval;
        parameters.parallel_jobs = this->build_jobs;
        parameters.show_progress = this->show_progress;
        return parameters;
    }
};

struct GraphHandler {
    enum graph_type { graph_none, graph_path, graph_source, graph_gbz };

    std::unique_ptr<PathHandleGraph> path_graph = nullptr;
    std::unique_ptr<gbwtgraph::SequenceSource> sequence_source = nullptr;
    std::unique_ptr<gbwtgraph::GBZ> gbz_graph = nullptr;
    graph_type in_use = graph_none;

    // Returns a pointer to any stored `PathHandleGraph` or loads one according to
    // the config if there is no such graph.
    const PathHandleGraph* get_any_graph(const GBWTConfig& config);

    // Load the `PathHandleGraph` specified in the config and release other graphs.
    // No effect if the handler already contains a `PathHandleGraph`.
    void get_graph(const GBWTConfig& config);

    // Take the ownership of the provided `SequenceSource` and store it in the handler.
    // Releases other graphs.
    void use(std::unique_ptr<gbwtgraph::SequenceSource>& source);

    // Load the GBZ specified in the config, store it in the handler, and tell
    // the GBWTHandler to use the stored GBWT as an external GBWT.
    // NOTE: The graph will become invalid if the GBWT in the GBWTHandler changes.
    void load_gbz(GBWTHandler& gbwts, GBWTConfig& config);

    void clear();

    // If the handler contains a `SequenceSource`, serialize it according to the config.
    void serialize_segment_translation(const GBWTConfig& config) const;
};

//----------------------------------------------------------------------------

GBWTConfig parse_gbwt_config(int argc, char** argv);
void validate_gbwt_config(GBWTConfig& config);

void step_1_build_gbwts(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config);
void step_2_merge_gbwts(GBWTHandler& gbwts, GBWTConfig& config);
void step_3_alter_gbwt(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config);
void step_4_path_cover(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config);
void step_5_gbz(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config);
void step_6_r_index(GBWTHandler& gbwts, GBWTConfig& config);
void step_7_metadata(GBWTHandler& gbwts, GBWTConfig& config);
void step_8_paths(GBWTHandler& gbwts, GBWTConfig& config);

void report_time_memory(const std::string& what, double start_time, const GBWTConfig& config);
void print_metadata(std::ostream& out, const GBWTHandler& gbwts);

//----------------------------------------------------------------------------

int main_gbwt(int argc, char** argv) {
    GBWTConfig config = parse_gbwt_config(argc, argv);
    validate_gbwt_config(config);

    // Let GBWT operate silently.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);

    // This is the data we are using.
    GBWTHandler gbwts;
    gbwts.filename = config.gbwt_name;
    gbwts.show_progress = config.show_progress;
    GraphHandler graphs;

    // Input GBWT construction.
    if (config.build != GBWTConfig::build_none) {
        step_1_build_gbwts(gbwts, graphs, config);
    }

    // Merge multiple input GBWTs.
    if (config.merge != GBWTConfig::merge_none) {
        step_2_merge_gbwts(gbwts, config);
    }

    // Edit the GBWT (remove samples, apply tags).
    if (!config.to_remove.empty() || !config.tags_to_set.empty()) {
        step_3_alter_gbwt(gbwts, graphs, config);
    }

    // Path cover construction.
    if (config.path_cover != GBWTConfig::path_cover_none) {
        step_4_path_cover(gbwts, graphs, config);
    }

    // Now we can serialize the GBWT to a separate file.
    if (!config.gbwt_output.empty() && config.graph_output.empty()) {
        double start = gbwt::readTimer();
        gbwts.serialize(config.gbwt_output);
        report_time_memory("GBWT serialized", start, config);
    }

    // Serialize the segment translation if necessary.
    if (!config.segment_translation.empty()) {
        graphs.serialize_segment_translation(config);
    }

    // GBZ construction and serialization.
    if (!config.graph_output.empty()) {
        step_5_gbz(gbwts, graphs, config);
    }

    // We no longer need the graph, but we may need the GBWT inside a GBZ.
    if (graphs.in_use == GraphHandler::graph_gbz) {
        gbwts.use(graphs.gbz_graph->index);
    }
    graphs.clear();

    // R-index construction.
    if (!config.r_index_name.empty()) {
        step_6_r_index(gbwts, config);
    }

    // Metadata options.
    if (config.metadata_mode) {
        step_7_metadata(gbwts, config);
    }

    // Path options.
    if (config.path_mode) {
        step_8_paths(gbwts, config);
    }

    return 0;
}

//----------------------------------------------------------------------------

void help_gbwt(char** argv) {
    std::cerr << "usage: " << argv[0] << " gbwt [options] [args]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Manipulate GBWTs. Input GBWTs are loaded from input args" << std::endl;
    std::cerr << "or built in earlier steps. See wiki page \"VG GBWT Subcommand\"." << std::endl;
    std::cerr << "The input graph is provided with one of -x, -G, or -Z" << std::endl;
    std::cerr << std::endl;
    std::cerr << "General options:" << std::endl;
    std::cerr << "  -h, --help              print this help message to stderr and exit" << std::endl;
    std::cerr << "  -x, --xg-name FILE      read the graph from FILE" << std::endl;
    std::cerr << "  -o, --output FILE       write output GBWT to FILE" << std::endl;
    std::cerr << "  -d, --temp-dir DIR      use directory DIR for temporary files" << std::endl;
    std::cerr << "  -p, --progress          show progress and statistics" << std::endl;
    std::cerr << std::endl;
    std::cerr << "GBWT construction parameters (for steps 1 and 4):" << std::endl;
    std::cerr << "      --buffer-size N     construction buffer size in millions of nodes"
                                        <<  "[" << (gbwt::DynamicGBWT::INSERT_BATCH_SIZE / gbwt::MILLION) << "]" << std::endl;
    std::cerr << "      --id-interval N     store path IDs at 1/N positions "
                                        << "[" << gbwt::DynamicGBWT::SAMPLE_INTERVAL << "]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Multithreading:" << std::endl;
    std::cerr << "      --num-jobs N        use at most N parallel build jobs" << std::endl;
    std::cerr << "                          (for -v, -G, -A, -l, -P) " 
                                        << "[" << GBWTConfig::default_build_jobs() << "]" << std::endl;
    std::cerr << "      --num-threads N     use N parallel search threads" << std::endl;
    std::cerr << "                          (for -b and -r) [" << omp_get_max_threads() << "]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 1: GBWT construction (requires -o and one of { -v, -G, -Z, -E, A }):" << std::endl;
    std::cerr << "  -v, --vcf-input         index the haplotypes in the VCF files specified in" << std::endl;
    std::cerr << "                          input args in parallel (requires -x, implies -f);" << std::endl;
    std::cerr << "                          (inputs must be over different contigs," << std::endl;
    std::cerr << "                          does not store graph contigs in the GBWT)" << std::endl;
    std::cerr << "      --preset X          use preset X (available: 1000gp)" << std::endl;
    std::cerr << "      --inputs-as-jobs    create one build job for each input" << std::endl;
    std::cerr << "                          instead of using first-fit heuristic" << std::endl;
    std::cerr << "      --parse-only        store the VCF parses without building GBWTs" << std::endl;
    std::cerr << "                          (use -o for file name prefix; skips later steps)" << std::endl;
    std::cerr << "      --ignore-missing    don't warn when variants are missing from the graph" << std::endl;
    std::cerr << "      --actual-phasing    don't treat unphased homozygous genotypes as phased" << std::endl;
    std::cerr << "      --force-phasing     replace unphased genotypes with randomly phased ones" << std::endl;
    std::cerr << "      --discard-overlaps  skip overlapping alternate alleles if the overlap" << std::endl;
    std::cerr << "                          cannot be resolved instead of creating a phase break" << std::endl;
    // TODO source for the default
    std::cerr << "      --batch-size N      index the haplotypes in batches of N samples [200]" << std::endl; 
    std::cerr << "      --sample-range X-Y  index samples X to Y (inclusive, 0-based)" << std::endl;
    std::cerr << "      --rename V=P        VCF contig V matches path P in the graph (may repeat)" << std::endl;
    std::cerr << "      --vcf-variants      variants in graph use VCF contig names, not path names" << std::endl;
    std::cerr << "      --vcf-region C:X-Y  restrict VCF contig C to coordinates X to Y" << std::endl;
    std::cerr << "                          (inclusive, 1-based; may repeat)" << std::endl;
    std::cerr << "      --exclude-sample X  do not index the sample with name X" << std::endl;
    std::cerr << "                          (faster than -R; may repeat)" << std::endl;
    std::cerr << "  -G, --gfa-input         index walks or paths in the GFA file (one input arg)" << std::endl;
    std::cerr << "      --max-node N        chop long segments into nodes of at most N bp" << std::endl;
    std::cerr << "                          (use 0 to disable) [" << gbwtgraph::MAX_NODE_LENGTH << "]" << std::endl;
    std::cerr << "      --path-regex X      parse metadata as haplotypes from path names" << std::endl;
    std::cerr << "                          using regex X instead of vg-parser-compatible rules" << std::endl;
    std::cerr << "      --path-fields X     parse metadata as haplotypes, mapping regex submatches" << std::endl;
    std::cerr << "                          to these fields instead of vg-parser-compatible rules" << std::endl;
    std::cerr << "      --translation FILE  write the segment to node translation table to FILE" << std::endl;
    std::cerr << "  -Z, --gbz-input         extract GBWT & GBWTGraph from GBZ from (one) input arg" << std::endl;
    std::cerr << "  -E, --index-paths       index the embedded non-alt paths in the graph" << std::endl;
    std::cerr << "                          (requires -x, no input args)" << std::endl;
    std::cerr << "  -A, --alignment-input   index the alignments in the GAF files specified" << std::endl;
    std::cerr << "                          in input args (requires -x)" << std::endl;
    std::cerr << "      --gam-format        input files are in GAM format instead of GAF format" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 2: Merge multiple input GBWTs (requires -o):" << std::endl;
    std::cerr << "  -m, --merge             use the insertion algorithm" << std::endl;
    std::cerr << "  -f, --fast              fast merging algorithm (node ids must not overlap)" << std::endl;
    std::cerr << "  -b, --parallel          use the parallel algorithm" << std::endl;
    std::cerr << "      --chunk-size N      search in chunks of N sequences "
                                        << "[" << gbwt::MergeParameters::CHUNK_SIZE << "]" << std::endl;
    std::cerr << "      --pos-buffer N      use N MiB position buffers for each search thread "
                                        << "[" << gbwt::MergeParameters::POS_BUFFER_SIZE << "]" << std::endl;
    std::cerr << "      --thread-buffer N   use N MiB thread buffers for each search thread "
                                        << "[" << gbwt::MergeParameters::THREAD_BUFFER_SIZE << "]" << std::endl;
    std::cerr << "      --merge-buffers N   merge 2^N thread buffers into one file per merge "
                                        << "[" << gbwt::MergeParameters::MERGE_BUFFERS << "]" << std::endl;
    std::cerr << "      --merge-jobs N      run N parallel merge jobs "
                                        << "[" << GBWTConfig::default_merge_jobs() << "]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 3: Alter GBWT (requires -o and one input GBWT):" << std::endl;
    std::cerr << "  -R, --remove-sample X   remove sample X from the index (may repeat)" << std::endl;
    std::cerr << "      --set-tag K=V       set a GBWT tag (may repeat)" << std::endl;
    std::cerr << "      --set-reference X   set sample X as the reference (may repeat)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 4: Path cover GBWT construction " << std::endl;
    std::cerr << "(requires an input graph, -o, and one of { -a, -l, -P }):" << std::endl;
    std::cerr << "  -a, --augment-gbwt      add path cover of missing components (one input GBWT)" << std::endl;
    std::cerr << "  -l, --local-haplotypes  sample local haplotypes (one input GBWT)" << std::endl;
    std::cerr << "  -P, --path-cover        build a greedy path cover (no input GBWTs)" << std::endl;
    std::cerr << "  -n, --num-paths N       find N paths per component"
                                         << "[" << GBWTConfig::default_num_paths_local() << " for -l, " 
                                         << GBWTConfig::default_num_paths() << " otherwise]" << std::endl;
    std::cerr << "  -k, --context-length N  use N-node contexts [" << GBWTConfig::default_context_length() << "]" << std::endl;
    std::cerr << "      --pass-paths        include named graph paths in local haplotype" << std::endl;
    std::cerr << "                          or greedy path cover GBWT" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 5: GBWTGraph construction (requires an input graph and one input GBWT):" << std::endl;
    std::cerr << "  -g, --graph-name FILE   build GBZ graph and store it in FILE" << std::endl;
    std::cerr << "                          (makes -o unnecessary)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 6: R-index construction (one input GBWT):" << std::endl;
    std::cerr << "  -r, --r-index FILE      build an r-index and store it in FILE" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 7: Metadata (one input GBWT):" << std::endl;
    std::cerr << "  -M, --metadata          print basic metadata" << std::endl;
    std::cerr << "  -C, --contigs           print the number of contigs" << std::endl;
    std::cerr << "  -H, --haplotypes        print the number of haplotypes" << std::endl;
    std::cerr << "  -S, --samples           print the number of samples" << std::endl;
    std::cerr << "  -L, --list-names        list contig/sample names (use with -C or -S)" << std::endl;
    std::cerr << "  -T, --path-names        list path names" << std::endl;
    std::cerr << "      --tags              list GBWT tags" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Step 8: Paths (one input GBWT):" << std::endl;
    std::cerr << "  -c, --count-paths       print the number of paths" << std::endl;
    std::cerr << "  -e, --extract FILE      extract paths in SDSL format to FILE" << std::endl;
    std::cerr << std::endl;
}

//----------------------------------------------------------------------------

void use_preset(std::string preset_name, GBWTConfig& config) {
    for (char& c : preset_name) {
        c = std::tolower(c);
    }
    if (preset_name == "1000gp") {
        config.haplotype_indexer.gbwt_buffer_size = 200;
        config.haplotype_indexer.samples_in_batch = 100;
        config.haplotype_indexer.force_phasing = true;
        config.haplotype_indexer.discard_overlaps = true;
    } else {
        config.logger.error() << "unknown preset: " << preset_name << std::endl;
    }
}

void no_multiple_input_types(const GBWTConfig& config) {
    if (config.build != GBWTConfig::build_none) {
        config.logger.error() << "only one input type can be specified for step 1" << std::endl;
    }
}

void no_multiple_cover_types(const GBWTConfig& config) {
    if (config.path_cover != GBWTConfig::path_cover_none)
    {
        config.logger.error() << "only one path cover type can be specified for step 4" << std::endl;
    }
}

void check_tag_validity(const Logger& logger, const std::string& key, const std::string& value, 
                        const std::unordered_set<char>& prohibited, const std::string& description) {
    for (auto& letter : value) {
        if (prohibited.count(letter)) {
            // This letter isn't allowed.
            auto error_msg = logger.error();
            error_msg << "tag \"" << key << "\" contains prohibited character \"" << letter 
                      << "\". It needs to be " << description << " and may not contain any of: ";
            for (auto& c : prohibited) {
                error_msg << " '" << c << "'";
            }
            error_msg << std::endl;
        }
    }
}

GBWTConfig parse_gbwt_config(int argc, char** argv) {
    if (argc == 2) {
        help_gbwt(argv);
        std::exit(EXIT_FAILURE);
    }

    // Long options with no corresponding short options.
    constexpr int OPT_BUFFER_SIZE = 1000;
    constexpr int OPT_ID_INTERVAL = 1001;
    constexpr int OPT_NUM_JOBS = 1002;
    constexpr int OPT_NUM_THREADS = 1003;
    constexpr int OPT_PRESET = 1100;
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
    constexpr int OPT_MAX_NODE = 1114;
    constexpr int OPT_PATH_REGEX = 1115;
    constexpr int OPT_PATH_FIELDS = 1116;
    constexpr int OPT_TRANSLATION = 1117;
    constexpr int OPT_GAM_FORMAT = 1118;
    constexpr int OPT_CHUNK_SIZE = 1200;
    constexpr int OPT_POS_BUFFER = 1201;
    constexpr int OPT_THREAD_BUFFER = 1202;
    constexpr int OPT_MERGE_BUFFERS = 1203;
    constexpr int OPT_MERGE_JOBS = 1204;
    constexpr int OPT_SET_TAG = 1300;
    constexpr int OPT_SET_REFERENCE = 1301;
    constexpr int OPT_PASS_PATHS = 1400;
    constexpr int OPT_TAGS = 1700;

    // Deprecated options.
    constexpr int OPT_THREAD_NAMES = 2000;
    constexpr int OPT_COUNT_THREADS = 2001;
    constexpr int OPT_GBZ_FORMAT = 2002;

    // Make a collection of all the known tags and their descriptions. Use an ordered map so that we can do some typo guessing.
    // Values are description and list of prohibited characters.
    const std::map<std::string, std::pair<std::string, std::unordered_set<char>>> KNOWN_TAGS = {
        {
            gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG,
            {"a space-separated list of PanSN-valid sample/assembly names of references in the graph", {'#'}}
        }
    };

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

        // Multithreading parameters
        { "num-jobs", required_argument, 0, OPT_NUM_JOBS },
        { "num-threads", required_argument, 0, OPT_NUM_THREADS },

        // Input GBWT construction: VCF
        { "vcf-input", no_argument, 0, 'v' },
        { "preset", required_argument, 0, OPT_PRESET },
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

        // Input GBWT construction: GFA
        { "gfa-input", no_argument, 0, 'G' },
        { "max-node", required_argument, 0, OPT_MAX_NODE },
        { "path-regex", required_argument, 0, OPT_PATH_REGEX },
        { "path-fields", required_argument, 0, OPT_PATH_FIELDS },
        { "translation", required_argument, 0, OPT_TRANSLATION },

        // Input GBWT construction: GBZ
        { "gbz-input", no_argument, 0, 'Z' },

        // Input GBWT construction: paths
        { "index-paths", no_argument, 0, 'E' },

        // Input GBWT construction: GAF/GAM
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

        // Alter GBWT
        { "remove-sample", required_argument, 0, 'R' },
        { "set-tag", required_argument, 0, OPT_SET_TAG },
        { "set-reference", required_argument, 0, OPT_SET_REFERENCE },

        // Path cover
        { "augment-gbwt", no_argument, 0, 'a' },
        { "local-haplotypes", no_argument, 0, 'l' },
        { "path-cover", no_argument, 0, 'P' },
        { "num-paths", required_argument, 0, 'n' },
        { "context-length", required_argument, 0, 'k' },
        { "pass-paths", no_argument, 0, OPT_PASS_PATHS },

        // GBZ construction
        { "graph-name", required_argument, 0, 'g' },
        { "gbz-format", no_argument, 0, OPT_GBZ_FORMAT }, // Hidden option enabled by default.

        // R-index
        { "r-index", required_argument, 0, 'r' },

        // Metadata
        { "metadata", no_argument, 0, 'M' },
        { "contigs", no_argument, 0, 'C' },
        { "haplotypes", no_argument, 0, 'H' },
        { "samples", no_argument, 0, 'S' },
        { "list-names", no_argument, 0, 'L' },
        { "path-names", no_argument, 0, 'T' },
        { "thread-names", no_argument, 0, OPT_THREAD_NAMES },
        { "tags", no_argument, 0, OPT_TAGS },

        // Paths
        { "count-paths", no_argument, 0, 'c' },
        { "count-threads", no_argument, 0, OPT_COUNT_THREADS },
        { "extract", required_argument, 0, 'e' },

        { "help", no_argument, 0, 'h' },
        { 0, 0, 0, 0 }
    };

    int c;
    optind = 2; // force optind past command positional argument
    GBWTConfig config;
    while (true) {
        int option_index = 0;
        c = getopt_long(argc, argv, "x:o:d:pvGZEAmfbR:alPn:k:g:r:MCHSLTce:h?", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        // General
        case 'x':
            config.graph_name = require_exists(config.logger, optarg);
            break;
        case 'o':
            config.gbwt_output = ensure_writable(config.logger, optarg);
            break;
        case 'd':
            temp_file::set_dir(optarg);
            break;
        case 'p':
            config.show_progress = true;
            break;

        // GBWT construction parameters
        case OPT_BUFFER_SIZE:
            config.haplotype_indexer.gbwt_buffer_size = std::max(parse<size_t>(optarg), 1ul);
            config.gfa_parameters.automatic_batch_size = false; // User-defined buffer size overrides heuristics.
            break;
        case OPT_ID_INTERVAL:
            config.haplotype_indexer.id_interval = parse<size_t>(optarg);
            break;

        // Multithreading parameters
        case OPT_NUM_JOBS:
            config.build_jobs = parse<size_t>(optarg);
            break;
        case OPT_NUM_THREADS:
            config.search_threads = std::max(parse<size_t>(optarg), 1ul);
            break;

        // Input GBWT construction: VCF
        case 'v':
            no_multiple_input_types(config);
            config.build = GBWTConfig::build_vcf;
            config.produces_one_gbwt = true;
            break;
        case OPT_PRESET:
            use_preset(optarg, config);
            break;
        case OPT_INPUTS_AS_JOBS:
            config.inputs_as_jobs = true;
            break;
        case OPT_PARSE_ONLY:
            config.parse_only = true;
            break;
        case OPT_IGNORE_MISSING:
            config.haplotype_indexer.warn_on_missing_variants = false;
            break;
        case OPT_ACTUAL_PHASING:
            config.haplotype_indexer.phase_homozygous = false;
            break;
        case OPT_FORCE_PHASING:
            config.haplotype_indexer.force_phasing = true;
            break;
        case OPT_DISCARD_OVERLAPS:
            config.haplotype_indexer.discard_overlaps = true;
            break;
        case OPT_BATCH_SIZE:
            config.haplotype_indexer.samples_in_batch = std::max(parse<size_t>(optarg), 1ul);
            break;
        case OPT_SAMPLE_RANGE:
            {
                size_t first, last;
                tie(first, last) = parse_pair<size_t, size_t>(config.logger, optarg, '-', "--sample-range");
                config.haplotype_indexer.sample_range.first = first;
                config.haplotype_indexer.sample_range.second = last + 1;
            }
            break;
        case OPT_RENAME:
            {
                string vcf_contig, graph_contig;
                tie(vcf_contig, graph_contig) = parse_pair(config.logger, optarg, '=', "--rename");
                // Add the name mapping
                config.haplotype_indexer.path_to_vcf[graph_contig] = vcf_contig;
            }
            break;
        case OPT_VCF_VARIANTS:
            config.haplotype_indexer.rename_variants = false;
            break;
        case OPT_VCF_REGION:
            {
                // Parse contig:first-last
                std::string region(optarg);
                Region parsed;
                parse_region(region, parsed);
                if (parsed.start <= 0 || parsed.end <= 0) {
                    // We need both range bounds, and we can't accept 0 since input is 1-based.
                    config.logger.error() << "cannot parse 1-based region " << region << std::endl;
                }
                // Make sure to correct the coordinates to 0-based exclusive-end, from 1-based inclusive-end
                config.haplotype_indexer.regions[parsed.seq] = std::make_pair((size_t) (parsed.start - 1), (size_t) parsed.end);
            }
            break;
        case OPT_EXCLUDE_SAMPLE:
            config.haplotype_indexer.excluded_samples.insert(optarg);
            break;

        // Input GBWT construction: GFA
        case 'G':
            no_multiple_input_types(config);
            config.build = GBWTConfig::build_gfa;
            config.produces_one_gbwt = true;
            break;
        case OPT_MAX_NODE:
            config.gfa_parameters.max_node_length = parse<size_t>(optarg);
            break;
        case OPT_PATH_REGEX:
            if (config.gfa_parameters.path_name_formats.size() != 1) {
                // We need to override the existing rules we set up when we set
                // up the config, and make sure we have a place for the other
                // option.
                config.gfa_parameters.path_name_formats.clear();
                config.gfa_parameters.path_name_formats.emplace_back("", "", PathSense::HAPLOTYPE);
            }
            config.gfa_parameters.path_name_formats.back().regex = optarg;
            break;
        case OPT_PATH_FIELDS:
            if (config.gfa_parameters.path_name_formats.size() != 1) {
                // We need to override the existing rules we set up when we set
                // up the config, and make sure we have a place for the other
                // option.
                config.gfa_parameters.path_name_formats.clear();
                config.gfa_parameters.path_name_formats.emplace_back("", "", PathSense::HAPLOTYPE);
            }
            config.gfa_parameters.path_name_formats.back().fields = optarg;
            break;
        case OPT_TRANSLATION:
            config.segment_translation = ensure_writable(config.logger, optarg);
            break;

        // Input GBWT construction: GBZ
        case 'Z':
            no_multiple_input_types(config);
            config.build = GBWTConfig::build_gbz;
            config.produces_one_gbwt = true;
            break;

        // Input GBWT construction: Paths
        case 'E':
            no_multiple_input_types(config);
            config.build = GBWTConfig::build_paths;
            config.produces_one_gbwt = true;
            break;

        // Input GBWT construction: GAF/GAM
        case 'A':
            no_multiple_input_types(config);
            config.build = GBWTConfig::build_alignments;
            config.produces_one_gbwt = true;
            break;
        case OPT_GAM_FORMAT:
            config.gam_format = true;
            break;

        // Merging
        case 'm':
            config.merge = GBWTConfig::merge_insert;
            config.produces_one_gbwt = true;
            break;
        case 'f':
            config.merge = GBWTConfig::merge_fast;
            config.produces_one_gbwt = true;
            break;
        case 'b':
            config.merge = GBWTConfig::merge_parallel;
            config.produces_one_gbwt = true;
            break;
        case OPT_CHUNK_SIZE:
            config.merge_parameters.setChunkSize(parse<size_t>(optarg));
            break;
        case OPT_POS_BUFFER:
            config.merge_parameters.setPosBufferSize(parse<size_t>(optarg));
            break;
        case OPT_THREAD_BUFFER:
            config.merge_parameters.setThreadBufferSize(parse<size_t>(optarg));
            break;
        case OPT_MERGE_BUFFERS:
            config.merge_parameters.setMergeBuffers(parse<size_t>(optarg));
            break;
        case OPT_MERGE_JOBS:
            config.merge_parameters.setMergeJobs(parse<size_t>(optarg));
            break;

        // Alter GBWT
        case 'R':
            config.to_remove.insert(optarg);
            break;
        case OPT_SET_TAG:
            {
                string tag_name, tag_value;
                tie(tag_name, tag_value) = parse_pair(config.logger, optarg, '=', "--set-tag");
                // See if this tag is known
                auto tag_record = KNOWN_TAGS.lower_bound(tag_name);
                if (tag_record == KNOWN_TAGS.end() && !KNOWN_TAGS.empty()) {
                    // This tag is larger than all known tags. Closest match is last tag.
                    --tag_record;
                }
                if (tag_record != KNOWN_TAGS.end()) {
                    auto& tag_description = tag_record->second.first;
                    auto& tag_prohibited_characters = tag_record->second.second;
                    // Tag is either known, or is unknown but there's a known tag to compare it with.
                    if (tag_name != tag_record->first) {
                        // This is an unknown tag, but we have an idea what it should be.
                        config.logger.warn() << "tag \"" << tag_name 
                                             << "\" is not a tag with a meaning recognized by vg; "
                                             << "maybe you meant \"" << tag_record->first << "\" "
                                             << "which would be " << tag_description << endl;
                    } else {
                        // This is a known tag, so validate it.
                        check_tag_validity(config.logger, tag_name, tag_value, 
                                           tag_prohibited_characters, tag_description);
                    }
                }
                config.tags_to_set.emplace(tag_name, tag_value);
            }
            break;
        case OPT_SET_REFERENCE:
            {
                const std::string& key = gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG;
                auto tag_record = KNOWN_TAGS.find(key);
                std::string sample_name = optarg;
                auto prohibited = tag_record->second.second;
                prohibited.insert(' ');
                check_tag_validity(config.logger, key, sample_name, prohibited, tag_record->second.first);
                auto iter = config.tags_to_set.find(key);
                if (iter != config.tags_to_set.end()) {
                    iter->second += " " + sample_name;
                } else {
                    config.tags_to_set.emplace(key, sample_name);
                }
            }
            break;

        // Path cover
        case 'a':
            no_multiple_cover_types(config);
            config.path_cover = GBWTConfig::path_cover_augment;
            break;
        case 'l':
            no_multiple_cover_types(config);
            config.path_cover = GBWTConfig::path_cover_local;
            if (!config.num_paths_set) {
                config.num_paths = GBWTConfig::default_num_paths_local();
            }
            break;
        case 'P':
            no_multiple_cover_types(config);
            config.path_cover = GBWTConfig::path_cover_greedy;
            config.produces_one_gbwt = true;
            break;
        case 'n':
            config.num_paths = parse<size_t>(optarg);
            config.num_paths_set = true;
            break;
        case 'k':
            config.context_length = parse<size_t>(optarg);
            break;
        case OPT_PASS_PATHS:
            config.include_named_paths = true;
            break;

        // GBWTGraph
        case 'g':
            config.graph_output = ensure_writable(config.logger, optarg);
            break;
        case OPT_GBZ_FORMAT:
            config.gbz_format = true;
            break;

        // Build r-index
        case 'r':
            config.r_index_name = ensure_writable(config.logger, optarg);
            break;

        // Metadata
        case 'M':
            config.metadata = true;
            config.metadata_mode = true;
            break;
        case 'C':
            config.contigs = true;
            config.metadata_mode = true;
            break;
        case 'H':
            config.haplotypes = true;
            config.metadata_mode = true;
            break;
        case 'S':
            config.samples = true;
            config.metadata_mode = true;
            break;
        case 'L':
            config.list_names = true;
            config.metadata_mode = true;
            break;
        case 'T':
            config.path_names = true;
            config.metadata_mode = true;
            break;
        case OPT_THREAD_NAMES:
            config.logger.warn() << "option --thread-names is deprecated; use --path-names instead" << endl;
            config.path_names = true;
            config.metadata_mode = true;
            break;
        case OPT_TAGS:
            config.tags = true;
            config.metadata_mode = true;
            break;

        // Paths
        case 'c':
            config.count_paths = true;
            config.path_mode = true;
            break;
        case OPT_COUNT_THREADS:
            config.logger.warn() << "option --count-threads is deprecated; use --count-paths instead" << endl;
            config.count_paths = true;
            config.path_mode = true;
            break;
        case 'e':
            config.path_output = ensure_writable(config.logger, optarg);
            config.path_mode = true;
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

    // The remaining args are input args.
    for (int arg = optind; arg < argc; arg++) {
        config.input_filenames.push_back(argv[arg]);
    }
    // We can load a single input GBWT if we did not use any build options.
    if (config.input_filenames.size() == 1 && config.build == GBWTConfig::build_none) {
        config.gbwt_name = config.input_filenames.front();
    }

    // Copy information from primary fields to redundant fields.
    config.haplotype_indexer.show_progress = config.show_progress;
    config.gfa_parameters.show_progress = config.show_progress;
    config.gfa_parameters.parallel_jobs = config.build_jobs;
    config.gfa_parameters.batch_size = config.haplotype_indexer.gbwt_buffer_size * gbwt::MILLION;
    config.gfa_parameters.sample_interval = config.haplotype_indexer.id_interval;
    
    return config;
}

//----------------------------------------------------------------------------

void validate_gbwt_config(GBWTConfig& config) {
    // We can either write GBWT in SDSL format to a separate file or as part of a GBZ graph.
    // However, `--parse-only` uses `gbwt_output` for other purposes.
    bool has_gbwt_output =
        (!config.gbwt_output.empty() || (!config.graph_output.empty() && !config.parse_only));

    // We have one input GBWT after steps 1-4.
    bool one_input_gbwt = config.input_filenames.size() == 1 || config.produces_one_gbwt;

    // We can load a PathHandleGraph from a file, get a SequenceSource from parsing GFA, or get a GBWTGraph from GBZ or GG/GBWT.
    bool has_graph_input = (!config.graph_name.empty() || config.build == GBWTConfig::build_gfa 
                          || config.build == GBWTConfig::build_gbz);

    if (config.build == GBWTConfig::build_gbz) {
        // If we "build" a GBWT by loading it from a GBZ, we just need to make
        // sure that we know enough to actually load it.  
        if (!config.graph_name.empty()) {
            config.logger.error() << "GBZ input does not use -x" << std::endl;
        }
        if (config.input_filenames.size() != 1) {
            config.logger.error() << "GBZ input requires one input arg" << std::endl;
        }
    } else if (config.build != GBWTConfig::build_none) {
        if (!has_gbwt_output) {
            // If we build our GBWT by doing anything other than loading it
            // from a GBZ, we need to have somewhere to put it.
            config.logger.error() << "GBWT construction requires output GBWT" << std::endl;
        }
        if (config.build == GBWTConfig::build_vcf) {
            if (config.graph_name.empty() || config.input_filenames.empty()) {
                config.logger.error() << "GBWT construction from VCF requires -x and input args" << std::endl;
            }
            if (config.parse_only) {
                config.haplotype_indexer.batch_file_prefix = config.gbwt_output;
            }
        } else if (config.build == GBWTConfig::build_gfa) {
            if (!config.graph_name.empty()) {
                config.logger.error() << "GBWT construction from GFA does not use -x" << std::endl;
            }
            if (config.input_filenames.size() != 1) {
                config.logger.error() << "GBWT construction from GFA requires one input arg" << std::endl;
            }
        } else if (config.build == GBWTConfig::build_alignments) {
            if (config.graph_name.empty() || config.input_filenames.empty()) {
                config.logger.error() << "GBWT construction from alignments requires -x and input args" << std::endl;
            }
        } else if (config.build == GBWTConfig::build_paths) {
            if (config.graph_name.empty()) {
                config.logger.error() << "GBWT construction from embedded paths requires -x" << std::endl;
            }
            if (!config.input_filenames.empty()) {
                config.logger.error() << "GBWT construction from embedded paths does not use input args" << std::endl;
            }
        }
    }

    if (config.merge != GBWTConfig::merge_none) {
        if (config.input_filenames.size() < 2 || !has_gbwt_output) {
            config.logger.error() << "merging requires at multiple input GBWTs and output GBWT" << std::endl;
        }
    }

    if (!config.to_remove.empty()) {
        if (config.build == GBWTConfig::build_gbz) {
            config.logger.error() << "the GBWT extracted from GBZ cannot have paths modified" << std::endl;
        }
        if (!(config.input_filenames.size() == 1 || config.merge != GBWTConfig::merge_none) || !has_gbwt_output) {
            config.logger.error() << "removing samples requires one input GBWT and output GBWT" << std::endl;
        }
    }
    
    if (!config.tags_to_set.empty()) {
        if (!(config.input_filenames.size() == 1 || config.merge != GBWTConfig::merge_none) || !has_gbwt_output) {
            config.logger.error() << "setting tags requires one input GBWT and output GBWT" << std::endl;
        }
    }

    if (config.path_cover != GBWTConfig::path_cover_none) {
        if (!has_gbwt_output || (config.graph_name.empty() && config.build != GBWTConfig::build_gbz)) {
            // Path cover options needs a graph. We can use the provided graph or the GBZ
            // we took as an input. In the latter case, we know that the corresponding GBWT has not
            // been modified and the graph is hence safe to use.
            config.logger.error() << "path cover options require an input graph and output GBWT" << std::endl;
        }
        if (config.path_cover == GBWTConfig::path_cover_greedy && !config.input_filenames.empty()) {
            config.logger.error() << "greedy path cover does not use input GBWTs" << std::endl;
        }
        if ((config.path_cover == GBWTConfig::path_cover_local || config.path_cover == GBWTConfig::path_cover_augment)
            && !(config.input_filenames.size() == 1 || config.merge != GBWTConfig::merge_none)) {
            config.logger.error() << "path cover options -a and -l require one input GBWT" << std::endl;
        }
        if (config.num_paths == 0) {
            config.logger.error() << "number of paths must be non-zero for path cover" << std::endl;
        }
        if (config.context_length < gbwtgraph::PATH_COVER_MIN_K) {
            config.logger.error() << "context length must be at least " 
                                  << gbwtgraph::PATH_COVER_MIN_K << " for path cover" << std::endl;
        }
    }

    if (!config.segment_translation.empty()) {
        if (config.build != GBWTConfig::build_gfa && config.build != GBWTConfig::build_gbz) {
            config.logger.error() << "node-to-segment translation requires GFA or GBZ input" << std::endl;
        }
    }

    if (!config.graph_output.empty()) {
        if (!has_graph_input || !one_input_gbwt) {
            config.logger.error() << "GBZ construction requires an input graph and one input GBWT" << std::endl;
        }
    }

    if (!config.r_index_name.empty()) {
        if (!one_input_gbwt) {
            config.logger.error() << "r-index construction requires one input GBWT" << std::endl;
        }
    }

    if (config.metadata_mode) {
        if (!one_input_gbwt) {
            config.logger.error() << "metadata operations require one input GBWT" << std::endl;
        }
    }

    if (config.path_mode) {
        if (!one_input_gbwt) {
            config.logger.error() << "path operations require one input GBWT" << std::endl;
        }
    }
    
    for (auto& format : config.gfa_parameters.path_name_formats) {
        // Check the path name format. We don't check the regex syntax here,
        // but we can make sure that we're asking for things consistent with
        // the path sense so that we can make a GBWTGraph out of them later.
        // TODO: Do we still need to let people make GBWTs that can't make a GBWTGraph?
        if (format.regex.empty()) {
            config.logger.error() << "path name format regex is missing" << std::endl;
        }
        switch(format.sense) {
        case PathSense::GENERIC:
            if (format.fields.find("C") == std::string::npos && format.fields.find("c") == std::string::npos) {
                config.logger.error() << "path name fields do not set required contig for regex" 
                                      << format.regex << "and fields " << format.fields << std::endl;
            }
            if (format.fields.find("S") != std::string::npos || format.fields.find("s") != std::string::npos) {
                config.logger.error() << "path name fields set unusable sample for regex" 
                                      << format.regex << "and fields " << format.fields << std::endl;
            }
            if (format.fields.find("H") != std::string::npos || format.fields.find("h") != std::string::npos) {
                config.logger.error() << "path name fields set unusable haplotype for regex" 
                                      << format.regex << "and fields " << format.fields << std::endl;
            }
            break;
        case PathSense::HAPLOTYPE:
            if (format.fields.find("S") == std::string::npos && format.fields.find("s") == std::string::npos) {
                config.logger.error() << "path name fields do not set required sample for regex" 
                                      << format.regex << "and fields " << format.fields << std::endl;
            }
            // Fall-through because haplotypes also need contigs.
        case PathSense::REFERENCE: 
            if (format.fields.find("C") == std::string::npos && format.fields.find("c") == std::string::npos) {
                config.logger.error() << "path name fields do not set required contig for regex" 
                                      << format.regex << "and fields " << format.fields << std::endl;
            }
            break;
        default:
            config.logger.error() << "path name sense is unimplemented " << (int)format.sense << std::endl;
        }
    }
}

//----------------------------------------------------------------------------

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

std::vector<job_type> determine_jobs(std::unique_ptr<PathHandleGraph>& graph, const GBWTConfig& config) {

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
    for (size_t i = 0; i < config.input_filenames.size(); i++) {
        paths_by_file.push_back({ i, {} });
    }

    // Determine which VCF file contains each path.
    std::vector<size_t> path_found_in(paths.size(), config.input_filenames.size());
    for (size_t i = 0; i < config.input_filenames.size(); i++) {
        std::string filename = config.input_filenames[i];
        vcflib::VariantCallFile variant_file;
        variant_file.parseSamples = false;
        variant_file.open(filename);
        if (!variant_file.is_open()) {
            config.logger.error() << "cannot open VCF file " << filename << std::endl;
        }
        for (size_t j = 0; j < paths.size(); j++) {
            std::string contig_name = graph->get_path_name(paths[j].first);
            if (config.haplotype_indexer.path_to_vcf.find(contig_name) != config.haplotype_indexer.path_to_vcf.end()) {
                contig_name = config.haplotype_indexer.path_to_vcf.at(contig_name);
            }
            variant_file.setRegion(contig_name);
            vcflib::Variant var(variant_file);
            if (!(variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == contig_name)) {
                continue;
            }
            if (path_found_in[j] < config.input_filenames.size()) {
                config.logger.error() << "contig " << contig_name << " found in files " 
                                      << config.input_filenames[path_found_in[j]] 
                                      << " and " << filename << std::endl;
            }
            paths_by_file[i].insert(paths[j]);
            path_found_in[j] = i;
        }
    }

    // Special case: Each input file is a single job.
    if (config.inputs_as_jobs) {
        for (const vcf_paths& curr : paths_by_file) {
            if (curr.empty()) {
                continue;
            }
            job_type job({ config.input_filenames[curr.file], {}, 0 });
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
                jobs.push_back({ config.input_filenames[curr.file], {}, 0 });
                jobs.back().insert(path);
            }
        }
        result.insert(result.end(), jobs.begin(), jobs.end());
    }

    // Sort the jobs in descending order by size.
    std::sort(result.begin(), result.end());
    return result;
}

void use_or_save(const Logger& logger, std::unique_ptr<gbwt::DynamicGBWT>& index, GBWTHandler& gbwts,
                 std::vector<std::string>& filenames, size_t i, bool show_progress) {
    if (filenames.size() == 1) {
        gbwts.use(*index);
    } else {
        std::string temp = temp_file::create("gbwt-" + std::to_string(i) + "-");
        if (show_progress) {
            #pragma omp critical
            {
                logger.info() << "Job " << i << ": Saving the GBWT to " << temp << std::endl;
            }
        }
        save_gbwt(*index, temp, false);
        filenames[i] = temp;
    }
}

void step_1_build_gbwts(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        config.logger.info() << "Building input GBWTs" << std::endl;
    }
    gbwts.unbacked(); // We will build a new GBWT.
    if (config.build != GBWTConfig::build_gfa && config.build != GBWTConfig::build_gbz) {
        graphs.get_graph(config);
    }

    if (config.build == GBWTConfig::build_vcf) {
        if (config.show_progress) {
            config.logger.info() << "Input type: VCF" << std::endl;
        }
        omp_set_num_threads(config.build_jobs);
        // Process each VCF contig corresponding to a non-alt path.
        std::vector<job_type> jobs = determine_jobs(graphs.path_graph, config);
        if (jobs.size() > 1 && config.merge == GBWTConfig::merge_none) {
            config.merge = GBWTConfig::merge_fast;
        }
        std::vector<std::vector<std::string>> vcf_parses(jobs.size());
        if (config.show_progress) {
            config.logger.info() << "Parsing " << jobs.size() << " VCF files using up to " 
                                 << config.build_jobs << " parallel jobs" << std::endl;
        }
        #pragma omp parallel for schedule(dynamic, 1)
        for (size_t i = 0; i < jobs.size(); i++) {
            std::string job_name = "Job " + std::to_string(i);
            if (config.show_progress) {
                #pragma omp critical
                {
                    auto job_msg = config.logger.info();
                    job_msg << job_name << ": File " << jobs[i].filename << ", paths {";
                    for (path_handle_t handle : jobs[i].paths) {
                        job_msg << " " << graphs.path_graph->get_path_name(handle);
                    }
                    job_msg << " }" << std::endl;
                }
            }
            vcf_parses[i] = config.haplotype_indexer.parse_vcf(
                jobs[i].filename, *(graphs.path_graph), jobs[i].paths, job_name);
        }
        graphs.clear(); // Delete the graph to save memory.
        if (!config.parse_only) {
            std::vector<std::string> gbwt_files(vcf_parses.size(), "");
            if (config.show_progress) {
                config.logger.info() << "Building " << vcf_parses.size() << " GBWTs using up to " 
                                     << config.build_jobs << " parallel jobs" << std::endl;
            }
            #pragma omp parallel for schedule(dynamic, 1)
            for (size_t i = 0; i < vcf_parses.size(); i++) {
                std::string job_name = "Job " + std::to_string(i);
                std::unique_ptr<gbwt::DynamicGBWT> parsed = config.haplotype_indexer.build_gbwt(vcf_parses[i], job_name);
                use_or_save(config.logger, parsed, gbwts, gbwt_files, i, config.show_progress);
            }
            if (vcf_parses.size() > 1) {
                config.input_filenames = gbwt_files; // Use the temporary GBWTs as inputs.
            }
        }
    } else if (config.build == GBWTConfig::build_gfa) {
        if(config.show_progress) {
            config.logger.info() << "Input type: GFA" << std::endl;
        }
        auto result = gbwtgraph::gfa_to_gbwt(config.input_filenames.front(), config.gfa_parameters);
        if (result.first.get() == nullptr || result.second.get() == nullptr) {
            config.logger.error() << "GBWT construction from GFA failed" << std::endl;
        }
        gbwts.use(*(result.first));
        graphs.use(result.second);
    } else if (config.build == GBWTConfig::build_gbz) {
        if(config.show_progress) {
            config.logger.info() << "Input type: GBZ" << std::endl;
        }
        graphs.load_gbz(gbwts, config);
    } else if (config.build == GBWTConfig::build_paths) {
        if(config.show_progress) {
            config.logger.info() << "Input type: embedded paths" << std::endl;
        }
        std::unique_ptr<gbwt::DynamicGBWT> temp = config.haplotype_indexer.build_gbwt(*(graphs.path_graph));
        gbwts.use(*temp);
    } else if (config.build == GBWTConfig::build_alignments) {
        if (config.show_progress) {
            config.logger.info() << "Input type: " << (config.gam_format ? "GAM" : "GAF") << std::endl;
        }
        std::unique_ptr<gbwt::GBWT> temp =
            config.haplotype_indexer.build_gbwt(*(graphs.path_graph), config.input_filenames, 
                                                (config.gam_format ? "GAM" : "GAF"), config.build_jobs);
        gbwts.use(*temp);
    }

    report_time_memory("GBWTs built", start, config);
    if (config.parse_only) {
        std::exit(EXIT_SUCCESS); // VCF parsing does not produce GBWTs to continue with.
    }
}
//----------------------------------------------------------------------------

void step_2_merge_gbwts(GBWTHandler& gbwts, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        std::string algo_name;
        if (config.merge == GBWTConfig::merge_fast) {
            algo_name = "fast";
        } else if (config.merge == GBWTConfig::merge_insert) {
            algo_name = "insertion";
        } else if (config.merge == GBWTConfig::merge_parallel) {
            algo_name = "parallel";
        }
        config.logger.info() << "Merging " << config.input_filenames.size()
                      << " input GBWTs (" << algo_name << " algorithm)" << std::endl;
    }

    if (config.merge == GBWTConfig::merge_fast) {
        std::vector<gbwt::GBWT> indexes(config.input_filenames.size());
        for (size_t i = 0; i < config.input_filenames.size(); i++) {
            load_gbwt(indexes[i], config.input_filenames[i], config.show_progress);
        }
        if (config.show_progress) {
            config.logger.info() << "Merging the GBWTs" << std::endl;
        }
        gbwt::GBWT merged(indexes);
        gbwts.use(merged);
    } else if (config.merge == GBWTConfig::merge_insert) {
        gbwts.filename = config.input_filenames.front();
        gbwts.use_dynamic();
        for (size_t i = 1; i < config.input_filenames.size(); i++) {
            gbwt::GBWT next;
            load_gbwt(next, config.input_filenames[i], config.show_progress);
            if (next.size() > 2 * gbwts.dynamic.size()) {
                config.logger.warn() << "merging " << config.input_filenames[i] << " into a substantially smaller index; "
                                     << "merging would be faster in another order" << std::endl;
            }
            if (config.show_progress) {
                config.logger.info() << "Inserting " << next.sequences()
                                     << " sequences of total length " << next.size() << std::endl;
            }
            gbwts.dynamic.merge(next);
        }
    } else if (config.merge == GBWTConfig::merge_parallel) {
        gbwts.filename = config.input_filenames.front();
        gbwts.use_dynamic();
        omp_set_num_threads(config.search_threads);
        for (size_t i = 1; i < config.input_filenames.size(); i++) {
            gbwt::DynamicGBWT next;
            load_gbwt(next, config.input_filenames[i], config.show_progress);
            if (next.size() > 2 * gbwts.dynamic.size()) {
                config.logger.warn() << "merging " << config.input_filenames[i] << " into a substantially smaller index; "
                                     << "merging would be faster in another order" << std::endl;
            }
            if (config.show_progress) {
                config.logger.info() << "Inserting " << next.sequences()
                                     << " sequences of total length " << next.size() << std::endl;
            }
            gbwts.dynamic.merge(next, config.merge_parameters);
        }
    }
    gbwts.unbacked(); // We modified the GBWT.

    if (config.show_progress) {
        print_metadata(std::cerr, gbwts);
        report_time_memory("GBWTs merged", start, config);
    }
}

//----------------------------------------------------------------------------

void remove_samples(GBWTHandler& gbwts, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        config.logger.info() << "Removing " << config.to_remove.size()
                             << " sample(s) from the index" << std::endl;
    }

    gbwts.use_dynamic();
    if (!(gbwts.dynamic.hasMetadata() && gbwts.dynamic.metadata.hasPathNames() 
          && gbwts.dynamic.metadata.hasSampleNames())) {
        config.logger.error() << "the index does not contain metadata with path and sample names" << std::endl;
    }

    // Remove the samples one at a time, because old sample/path ids may be invalidated.
    for (const std::string& sample_name : config.to_remove) {
        gbwt::size_type sample_id = gbwts.dynamic.metadata.sample(sample_name);
        if (sample_id >= gbwts.dynamic.metadata.samples()) {
            config.logger.warn() << "the index does not contain sample " << sample_name << std::endl;
            continue;
        }
        std::vector<gbwt::size_type> path_ids = gbwts.dynamic.metadata.removeSample(sample_id);
        if (path_ids.empty()) {
            config.logger.info() << "no paths associated with sample " << sample_name << std::endl;
            continue;
        }
        if (config.show_progress) {
            config.logger.info() << "Removing " << path_ids.size()
                                 << " paths for sample " << sample_name << std::endl;
        }
        gbwts.dynamic.remove(path_ids);
    }
    gbwts.unbacked(); // We modified the GBWT.

    report_time_memory("Samples removed", start, config);
}

void set_tags(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        config.logger.info() << "Setting " << config.tags_to_set.size()
                             << " tags on the GBWT" << std::endl;
    }
    
    gbwts.use_compressed();
    gbwt::GBWT* index = gbwts.get_compressed();
    std::unordered_set<std::string> reference_samples;
    for (auto& kv : config.tags_to_set) {
        if (kv.first == gbwtgraph::REFERENCE_SAMPLE_LIST_GBWT_TAG && graphs.in_use == GraphHandler::graph_gbz) {
            // GBZ has a proper way of setting the reference samples.
            reference_samples = gbwtgraph::parse_reference_samples_tag(kv.second);
        } else {
            index->tags.set(kv.first, kv.second);
        }
    }
    // We modified the GBWT (we assume some tags got set)
    gbwts.unbacked();

    // If we updated reference samples and have an existing GBZ, we need to recache the named paths.
    // This is because we may want to include the reference paths in a path cover GBWT.
    if (!reference_samples.empty() && graphs.in_use == GraphHandler::graph_gbz) {
        graphs.gbz_graph->set_reference_samples(reference_samples);
    }

    report_time_memory("Tags set", start, config);
}

void step_3_alter_gbwt(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    if (!config.to_remove.empty()) {
        remove_samples(gbwts, config);
    }
    if (!config.tags_to_set.empty()) {
        set_tags(gbwts, graphs, config);
    }
}

//----------------------------------------------------------------------------

void step_4_path_cover(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        config.logger.info() << "Finding a " << config.num_paths
                             << "-path cover with context length " << config.context_length << std::endl;
    }

    // Select the appropriate graph.
    const PathHandleGraph* graph = graphs.get_any_graph(config);
    
    // We need to drop paths that are alt allele paths and not pass them
    // through from a graph that has them to the synthesized GBWT.
    std::function<bool(const path_handle_t&)> path_filter = [&graph](const path_handle_t& path) {
        return !Paths::is_alt(graph->get_path_name(path));
    };
    
    if (config.path_cover == GBWTConfig::path_cover_greedy) {
        if (config.show_progress) {
            config.logger.info() << "Algorithm: greedy" << std::endl;
        }
        gbwt::GBWT cover = gbwtgraph::path_cover_gbwt(
            *graph, config.path_cover_parameters(),
            config.include_named_paths, &path_filter
        );
        copy_reference_samples(*graph, cover);
        gbwts.use(cover);
    } else if (config.path_cover == GBWTConfig::path_cover_augment) {
        if (config.show_progress) {
            config.logger.info() << "Algorithm: augment" << std::endl;
        }
        gbwts.use_dynamic();
        gbwtgraph::augment_gbwt(*graph, gbwts.dynamic, config.path_cover_parameters());
    } else {
        if (config.show_progress) {
            config.logger.info() << "Algorithm: local haplotypes" << std::endl;
        }
        gbwts.use_compressed();
        gbwt::GBWT cover = gbwtgraph::local_haplotypes(
            *graph, *gbwts.get_compressed(), config.path_cover_parameters(),
            config.include_named_paths, &path_filter
        );
        copy_reference_samples(*gbwts.get_compressed(), cover);
        gbwts.use(cover);
    }
    gbwts.unbacked(); // We modified the GBWT.

    report_time_memory("Path cover built", start, config);
}

//----------------------------------------------------------------------------

void step_5_gbz(GBWTHandler& gbwts, GraphHandler& graphs, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        config.logger.info() << "Building GBZ" << std::endl;
    }

    gbwts.use_compressed();
    if (graphs.in_use == GraphHandler::graph_source) {
        // We avoid creating a duplicate GBWT by moving it into a std::unique_ptr.
        // Both the GBWT and the SequenceSource get reset during construction.
        std::unique_ptr<gbwt::GBWT> gbwt_ptr = std::make_unique<gbwt::GBWT>(std::move(*gbwts.get_compressed()));
        gbwtgraph::GBZ gbz(gbwt_ptr, graphs.sequence_source);
        graphs.clear(); // We no longer need the SequenceSource.
        save_gbz(gbz, config.graph_output, config.show_progress);
        gbwts.use(gbz.index); // We may need the GBWT later.
    } else if (graphs.in_use == GraphHandler::graph_gbz) {
        if (gbwts.in_use == GBWTHandler::index_external) {
            // We can serialize the GBZ we already have.
            // Step 3 may have changed the tags, but otherwise the GBWT is unchanged.
            save_gbz(*(graphs.gbz_graph), config.graph_output, config.show_progress);
            gbwts.use(graphs.gbz_graph->index); // We may need the GBWT later.
            graphs.clear(); // We no longer need the GBZ.
        } else {
            // We use the subgraph constructor with the new GBWT.
            gbwtgraph::GBZ gbz(std::move(*gbwts.get_compressed()), *(graphs.gbz_graph));
            graphs.clear(); // We no longer need the GBZ.
            save_gbz(gbz, config.graph_output, config.show_progress);
            gbwts.use(gbz.index); // We may need the GBWT later.
        }
    } else {
        graphs.get_graph(config);
        if (config.show_progress) {
            config.logger.info() << "Starting the construction" << std::endl;
        }
        gbwtgraph::GBZ gbz(
            std::move(*gbwts.get_compressed()), *(graphs.path_graph),
            vg::algorithms::find_translation(graphs.path_graph.get())
        );
        graphs.clear(); // We no longer need the graph.
        gbz.compute_pggname(nullptr); // We cannot determine the pggname of a generic HandleGraph efficiently.
        save_gbz(gbz, config.graph_output, config.show_progress);
        gbwts.use(gbz.index); // We may need the GBWT later.
    }

    report_time_memory("GBZ built", start, config);
}

//----------------------------------------------------------------------------

void step_6_r_index(GBWTHandler& gbwts, GBWTConfig& config) {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        config.logger.info() << "Building r-index" << std::endl;
    }

    omp_set_num_threads(config.search_threads);
    gbwts.use_compressed();
    if (config.show_progress) {
        config.logger.info() << "Starting the construction" << std::endl;
    }
    gbwt::FastLocate r_index(*gbwts.get_compressed());
    save_r_index(r_index, config.r_index_name, config.show_progress);

    report_time_memory("R-index built", start, config);
}

//----------------------------------------------------------------------------

void step_7_metadata(GBWTHandler& gbwts, GBWTConfig& config) {
    gbwts.use_compressed();
    gbwt::GBWT* index = gbwts.get_compressed();
    
    // Use this to get the metadata object for the operations that need it, or
    // fail if it's not there.
    auto get_metadata = [&index, &config]() -> const gbwt::Metadata& {
        if (!index->hasMetadata()) {
            config.logger.error() << "the GBWT does not contain metadata" << std::endl;
        }
        return index->metadata;
    };
    
    if (config.metadata) {
        // Make sure the metadata exists.
        get_metadata();
        print_metadata(std::cout, gbwts);
    }

    if (config.contigs) {
        auto& metadata = get_metadata();
        if (config.list_names) {
            if (index->metadata.hasContigNames()) {
                for (size_t i = 0; i < metadata.contigs(); i++) {
                    std::cout << metadata.contig(i) << std::endl;
                }
            } else {
                config.logger.error() << "the metadata does not contain contig names" << std::endl;
            }
        } else {
            std::cout << metadata.contigs() << std::endl;
        }
    }

    if (config.haplotypes) {
        std::cout << get_metadata().haplotypes() << std::endl;
    }

    if (config.samples) {
        auto& metadata = get_metadata();
        if (config.list_names) {
            if (metadata.hasSampleNames()) {
                for (size_t i = 0; i < metadata.samples(); i++) {
                    std::cout << metadata.sample(i) << std::endl;
                }
            } else {
                config.logger.error() << "the metadata does not contain sample names" << std::endl;
            }
        } else {
            std::cout << metadata.samples() << std::endl;
        }
    }

    if (config.path_names) {
        auto& metadata = get_metadata();
        if (metadata.hasPathNames()) {
            // Precompute some metadata
            auto gbwt_reference_samples = gbwtgraph::parse_reference_samples_tag(*index);
            for (size_t i = 0; i < metadata.paths(); i++) {
                PathSense sense = gbwtgraph::get_path_sense(*index, i, gbwt_reference_samples);
                std::cout << gbwtgraph::compose_path_name(*index, i, sense) << std::endl;
            }
        } else {
            config.logger.error() << "the metadata does not contain path names" << std::endl;
        }
    }
    
    if (config.tags) {
        // This only needs GBWT tag metadata.
        // TODO: the gbwt::Tags object doesn't have its own enumeration API.
        // Just reach in and grab them.
        for (auto& kv : index->tags.tags) {
            std::cout << kv.first << "\t" << kv.second << std::endl;
        }
    }
}

//----------------------------------------------------------------------------

void step_8_paths(GBWTHandler& gbwts, GBWTConfig& config) {
    // Extract paths in SDSL format.
    if (!config.path_output.empty()) {
        double start = gbwt::readTimer();
        if (config.show_progress) {
            config.logger.info() << "Extracting paths to " << config.path_output << std::endl;
        }
        gbwts.use_compressed();
        gbwt::GBWT* index = gbwts.get_compressed();
        if (config.show_progress) {
            config.logger.info() << "Starting the extraction" << std::endl;
        }
        gbwt::size_type node_width = gbwt::bit_length(index->sigma() - 1);
        gbwt::text_buffer_type out(config.path_output, std::ios::out, gbwt::MEGABYTE, node_width);
        for (gbwt::size_type id = 0; id < index->sequences(); id += 2) { // Ignore reverse complements.
            gbwt::vector_type sequence = index->extract(id);
            for (auto node : sequence) {
                out.push_back(node);
            }
            out.push_back(gbwt::ENDMARKER);
        }
        out.close();
        report_time_memory("Paths extracted", start, config);
    }

    // There are two sequences for each path.
    if (config.count_paths) {
        gbwts.use_compressed();
        std::cout << (gbwts.get_compressed()->sequences() / 2) << std::endl;
    }
}

//----------------------------------------------------------------------------

const PathHandleGraph* GraphHandler::get_any_graph(const GBWTConfig& config) {
    if (this->in_use == GraphHandler::graph_gbz) {
        return &(this->gbz_graph->graph);
    } else {
        this->get_graph(config);
        return this->path_graph.get();
    }
}

void GraphHandler::get_graph(const GBWTConfig& config) {
    if (this->in_use == graph_path) {
        return;
    } else {
        if (config.show_progress) {
            config.logger.info() << "Loading input graph from " << config.graph_name << std::endl;
        }
        this->clear();
        this->path_graph = vg::io::VPKG::load_one<PathHandleGraph>(config.graph_name);
        if (this->path_graph == nullptr) {
            config.logger.error() << "cannot load graph " << config.graph_name << std::endl;
        }
        this->in_use = graph_path;
    }
}

void GraphHandler::use(std::unique_ptr<gbwtgraph::SequenceSource>& source) {
    this->clear();
    this->sequence_source = std::move(source);
    this->in_use = graph_source;
}

void GraphHandler::load_gbz(GBWTHandler& gbwts, GBWTConfig& config) {
    if (this->in_use == graph_gbz) {
        return;
    } else {
        this->clear();
        this->gbz_graph = std::make_unique<gbwtgraph::GBZ>();
        vg::load_gbz(*this->gbz_graph, config.input_filenames.front(), config.show_progress);
        gbwts.use_external(this->gbz_graph->index);
        this->in_use = graph_gbz;
    }
}

void GraphHandler::clear() {
    this->path_graph.reset();
    this->sequence_source.reset();
    this->gbz_graph.reset();
    this->in_use = graph_none;
}

void GraphHandler::serialize_segment_translation(const GBWTConfig& config) const {
    double start = gbwt::readTimer();
    if (config.show_progress) {
        config.logger.info() << "Serializing segment to node translation to "
                             << config.segment_translation << std::endl;
    }
    std::ofstream out(config.segment_translation, std::ios_base::binary);

    if (this->in_use == graph_source) {
        if (this->sequence_source->uses_translation()) {
            auto& translation = this->sequence_source->segment_translation;
            for (auto iter = translation.begin(); iter != translation.end(); ++iter) {
                out << "T\t" << iter->first << "\t" << iter->second.first;
                for (nid_t i = iter->second.first + 1; i < iter->second.second; i++) {
                    out << "," << i;
                }
                out << "\n";
            }
        }
    } else if (this->in_use == graph_gbz) {
        this->gbz_graph->graph.for_each_segment([&](const std::string& name, std::pair<nid_t, nid_t> nodes) -> bool {
            out << "T\t" << name << "\t" << nodes.first;
            for (nid_t i = nodes.first + 1; i < nodes.second; i++) {
                out << "," << i;
            }
            out << "\n";
            return true;
        });
    }

    out.close();
    report_time_memory("Translation serialized", start, config);
}

//----------------------------------------------------------------------------

void report_time_memory(const std::string& what, double start_time, const GBWTConfig& config) {
    if (config.show_progress) {
        double seconds = gbwt::readTimer() - start_time;
        config.logger.info() << "" << what << " in " << seconds << " seconds, "
                             << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl << std::endl;
    }
}

void print_metadata(std::ostream& out, const GBWTHandler& gbwts) {
    if (gbwts.in_use == GBWTHandler::index_compressed) {
        const gbwt::GBWT* index = gbwts.get_compressed();
        gbwt::operator<<(out, index->metadata) << std::endl;
    } else if (gbwts.in_use == GBWTHandler::index_dynamic) {
        gbwt::operator<<(out, gbwts.dynamic.metadata) << std::endl;
    }
}

//----------------------------------------------------------------------------

// Register subcommand
static vg::subcommand::Subcommand vg_gbwt("gbwt", "build and manipulate GBWT and GBZ files", vg::subcommand::TOOLKIT, main_gbwt);
