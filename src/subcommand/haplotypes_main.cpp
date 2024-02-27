/** \file haplotypes_main.cpp
 *
 * Defines the "vg haplotypes" subcommand, which samples haplotypes by kmer counts in the reads.
 */

#include "subcommand.hpp"

#include "../hash_map.hpp"
#include "../recombinator.hpp"

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>
#include <unordered_map>

#include <getopt.h>
#include <omp.h>

#include <gbwtgraph/index.h>

using namespace vg;

//----------------------------------------------------------------------------

constexpr size_t DEFAULT_MAX_THREADS = 16;

size_t haplotypes_default_threads() {
    size_t threads = omp_get_max_threads();
    threads = std::max(threads, size_t(1));
    return std::min(threads, DEFAULT_MAX_THREADS);
}

constexpr size_t haplotypes_default_k() {
    return Haplotypes::Header::DEFAULT_K;
}

constexpr size_t haplotypes_default_w() {
    return gbwtgraph::Key64::WINDOW_LENGTH;
}

constexpr size_t haplotypes_default_subchain_length() {
    return HaplotypePartitioner::SUBCHAIN_LENGTH;
}

constexpr size_t haplotypes_default_n() {
    return Recombinator::NUM_HAPLOTYPES;
}

constexpr size_t haplotypes_default_candidates() {
    return Recombinator::NUM_CANDIDATES;
}

constexpr size_t haplotypes_default_coverage() {
    return Recombinator::COVERAGE;
}

constexpr double haplotypes_default_discount() {
    return Recombinator::PRESENT_DISCOUNT;
}

constexpr double haplotypes_default_adjustment() {
    return Recombinator::HET_ADJUSTMENT;
}

constexpr double haplotypes_default_absent() {
    return Recombinator::ABSENT_SCORE;
}

struct HaplotypesConfig {
    enum OperatingMode {
        mode_invalid,
        mode_sample_graph,
        mode_preprocess,
        mode_sample_haplotypes,
        mode_map_variants,
        mode_extract,
        mode_classify,
    };

    OperatingMode mode = mode_invalid;
    Haplotypes::Verbosity verbosity = Haplotypes::verbosity_silent;

    // File names.
    std::string graph_name;
    std::string gbz_output, haplotype_output, score_output, kmer_output;
    std::string distance_name, r_index_name;
    std::string haplotype_input, kmer_input, vcf_input;

    // Computational parameters.
    size_t k = haplotypes_default_k(), w = haplotypes_default_w();
    HaplotypePartitioner::Parameters partitioner_parameters;
    Recombinator::Parameters recombinator_parameters;

    // A prefix to add to VCF contig names to get GBWT contig names.
    std::string contig_prefix;

    // For extracting local haplotypes.
    size_t chain_id = std::numeric_limits<size_t>::max();
    size_t subchain_id = std::numeric_limits<size_t>::max();

    // Other parameters.
    size_t threads = haplotypes_default_threads();
    bool validate = false;

    HaplotypesConfig(int argc, char** argv, size_t max_threads);
};

void preprocess_graph(const gbwtgraph::GBZ& gbz, Haplotypes& haplotypes, HaplotypesConfig& config);

void sample_haplotypes(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config);

void map_variants(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config);

void extract_haplotypes(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config);

void classify_kmers(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config);

//----------------------------------------------------------------------------

int main_haplotypes(int argc, char** argv) {
    double start = gbwt::readTimer();
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
    size_t max_threads = omp_get_max_threads();
    omp_set_num_threads(haplotypes_default_threads());

    // Parse the arguments.
    HaplotypesConfig config(argc, argv, max_threads);

    // Load the graph.
    gbwtgraph::GBZ gbz;
    load_gbz(gbz, config.graph_name, config.verbosity >= Haplotypes::verbosity_basic);

    // Generate or load haplotype information.
    Haplotypes haplotypes;
    if (config.mode == HaplotypesConfig::mode_sample_graph || config.mode == HaplotypesConfig::mode_preprocess) {
        preprocess_graph(gbz, haplotypes, config);
    } else {
        if (config.verbosity >= Haplotypes::verbosity_basic) {
            std::cerr << "Loading haplotype information from " << config.haplotype_input << std::endl;
        }
        try {
            sdsl::simple_sds::load_from(haplotypes, config.haplotype_input);
        } catch (const std::runtime_error& e) {
            std::cerr << "error: [vg haplotypes] " << e.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // Save haplotype information if necessary.
    if (!config.haplotype_output.empty()) {
        if (config.verbosity >= Haplotypes::verbosity_basic) {
            std::cerr << "Writing haplotype information to " << config.haplotype_output << std::endl;
        }
        sdsl::simple_sds::serialize_to(haplotypes, config.haplotype_output);
    }

    // Sample the haplotypes.
    if (config.mode == HaplotypesConfig::mode_sample_graph || config.mode == HaplotypesConfig::mode_sample_haplotypes) {
        sample_haplotypes(gbz, haplotypes, config);
    }

    // Map variants to subchains.
    if (config.mode == HaplotypesConfig::mode_map_variants) {
        map_variants(gbz, haplotypes, config);
    }

    // Extract local haplotypes in FASTA format.
    if (config.mode == HaplotypesConfig::mode_extract) {
        extract_haplotypes(gbz, haplotypes, config);
    }

    // Classify kmers.
    if (config.mode == HaplotypesConfig::mode_classify) {
        classify_kmers(gbz, haplotypes, config);
    }

    if (config.verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        double gib = gbwt::inGigabytes(gbwt::memoryUsage());
        std::cerr << "Used " << seconds << " seconds, " << gib << " GiB" << std::endl;
    }
    return 0;
}

static vg::subcommand::Subcommand vg_haplotypes("haplotypes", "haplotype sampling based on kmer counts", vg::subcommand::TOOLKIT, main_haplotypes);

//----------------------------------------------------------------------------

void help_haplotypes(char** argv, bool developer_options) {
    std::string usage = "    " + std::string(argv[0]) + " " + std::string(argv[1]) + " [options] ";
    std::cerr << "Usage:" << std::endl;
    std::cerr << usage << "-k kmers.kff -g output.gbz graph.gbz" << std::endl;
    std::cerr << usage << "-H output.hapl graph.gbz" << std::endl;
    std::cerr << usage << "-i graph.hapl -k kmers.kff -g output.gbz graph.gbz" << std::endl;
    if (developer_options) {
        std::cerr << usage << "-i graph.hapl --vcf-input variants.vcf graph.gbz > output.tsv" << std::endl;
        std::cerr << usage << "-i graph.hapl -k kmers.kff --extract M:N graph.gbz > output.fa" << std::endl;
    }
    std::cerr << std::endl;

    std::cerr << "Haplotype sampling based on kmer counts." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Output files:" << std::endl;
    std::cerr << "    -g, --gbz-output X        write the output GBZ to X" << std::endl;
    std::cerr << "    -H, --haplotype-output X  write haplotype information to X" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Input files:" << std::endl;
    std::cerr << "    -d, --distance-index X    use this distance index (default: <basename>.dist)" << std::endl;
    std::cerr << "    -r, --r-index X           use this r-index (default: <basename>.ri)" << std::endl;
    std::cerr << "    -i, --haplotype-input X   use this haplotype information (default: generate)" << std::endl;
    std::cerr << "    -k, --kmer-input X        use kmer counts from this KFF file (required for --gbz-output)" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options for generating haplotype information:" << std::endl;
    std::cerr << "        --kmer-length N       kmer length for building the minimizer index (default: " << haplotypes_default_k() << ")" << std::endl;
    std::cerr << "        --window-length N     window length for building the minimizer index (default: " << haplotypes_default_w() << ")" << std::endl;
    std::cerr << "        --subchain-length N   target length (in bp) for subchains (default: " << haplotypes_default_subchain_length() << ")" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options for sampling haplotypes:" << std::endl;
    std::cerr << "        --coverage N          kmer coverage in the KFF file (default: estimate)" << std::endl;
    std::cerr << "        --num-haplotypes N    generate N haplotypes (default: " << haplotypes_default_n() << ")" << std::endl;
    std::cerr << "                              sample from N candidates (with --diploid-sampling; default: " << haplotypes_default_candidates() << ")" << std::endl;
    std::cerr << "        --present-discount F  discount scores for present kmers by factor F (default: " << haplotypes_default_discount() << ")" << std::endl;
    std::cerr << "        --het-adjustment F    adjust scores for heterozygous kmers by F (default: " << haplotypes_default_adjustment() << ")" << std::endl;
    std::cerr << "        --absent-score F      score absent kmers -F/+F (default: " << haplotypes_default_absent()  << ")" << std::endl;
    std::cerr << "        --diploid-sampling    choose the best pair from the greedily selected haplotypes" << std::endl;
    std::cerr << "        --include-reference   include named and reference paths in the output" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "    -v, --verbosity N         verbosity level (0 = silent, 1 = basic, 2 = detailed, 3 = debug; default: 0)" << std::endl;
    std::cerr << "    -t, --threads N           approximate number of threads (default: " << haplotypes_default_threads() << " on this system)" << std::endl;
    std::cerr << std::endl;
    if (developer_options) {
        std::cerr << "Developer options:" << std::endl;
        std::cerr << "        --validate            validate the generated information (may be slow)" << std::endl;
        std::cerr << "        --vcf-input X         map the variants in VCF file X to subchains" << std::endl;
        std::cerr << "        --contig-prefix X     a prefix for transforming VCF contig names into GBWT contig names" << std::endl;
        std::cerr << "        --extract M:N         extract haplotypes in chain M, subchain N in FASTA format" << std::endl;
        std::cerr << "        --score-output X      write haplotype scores to X" << std::endl;
        std::cerr << "        --classify X          classify kmers and write output to X" << std::endl;
        std::cerr << std::endl;
    }
}

//----------------------------------------------------------------------------

HaplotypesConfig::HaplotypesConfig(int argc, char** argv, size_t max_threads) {
    constexpr int OPT_KMER_LENGTH = 1200;
    constexpr int OPT_WINDOW_LENGTH = 1201;
    constexpr int OPT_SUBCHAIN_LENGTH = 1202;
    constexpr int OPT_COVERAGE = 1300;
    constexpr int OPT_NUM_HAPLOTYPES = 1301;
    constexpr int OPT_PRESENT_DISCOUNT = 1302;
    constexpr int OPT_HET_ADJUSTMENT = 1303;
    constexpr int OPT_ABSENT_SCORE = 1304;
    constexpr int OPT_DIPLOID_SAMPLING = 1305;
    constexpr int OPT_INCLUDE_REFERENCE = 1306;
    constexpr int OPT_VALIDATE = 1400;
    constexpr int OPT_VCF_INPUT = 1500;
    constexpr int OPT_CONTIG_PREFIX = 1501;
    constexpr int OPT_EXTRACT = 1600;
    constexpr int OPT_SCORE_OUTPUT = 1601;
    constexpr int OPT_CLASSIFY = 1602;

    static struct option long_options[] =
    {
        { "gbz-output", required_argument, 0, 'g' },
        { "haplotype-output", required_argument, 0, 'H' },
        { "distance-index", required_argument, 0, 'd' },
        { "r-index", required_argument, 0, 'r' },
        { "haplotype-input", required_argument, 0, 'i' },
        { "kmer-input", required_argument, 0, 'k' },
        { "kmer-length", required_argument, 0, OPT_KMER_LENGTH },
        { "window-length", required_argument, 0, OPT_WINDOW_LENGTH },
        { "subchain-length", required_argument, 0, OPT_SUBCHAIN_LENGTH },
        { "coverage", required_argument, 0, OPT_COVERAGE },
        { "num-haplotypes", required_argument, 0, OPT_NUM_HAPLOTYPES },
        { "present-discount", required_argument, 0, OPT_PRESENT_DISCOUNT },
        { "het-adjustment", required_argument, 0, OPT_HET_ADJUSTMENT },
        { "absent-score", required_argument, 0, OPT_ABSENT_SCORE },
        { "diploid-sampling", no_argument, 0, OPT_DIPLOID_SAMPLING },
        { "include-reference", no_argument, 0, OPT_INCLUDE_REFERENCE },
        { "verbosity", required_argument, 0, 'v' },
        { "threads", required_argument, 0, 't' },
        { "validate", no_argument, 0,  OPT_VALIDATE },
        { "vcf-input", required_argument, 0, OPT_VCF_INPUT },
        { "contig-prefix", required_argument, 0, OPT_CONTIG_PREFIX },
        { "extract", required_argument, 0, OPT_EXTRACT },
        { "score-output", required_argument, 0, OPT_SCORE_OUTPUT },
        { "classify", required_argument, 0, OPT_CLASSIFY },
        { 0, 0, 0, 0 }
    };

    // Process the arguments.
    int c;
    optind = 2; // force optind past command positional argument
    bool num_haplotypes_set = false;
    while (true) {
        int option_index = 0;
        c = getopt_long(argc, argv, "g:H:d:r:i:k:v:t:h", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'g':
            this->gbz_output = optarg;
            break;
        case 'H':
            this->haplotype_output = optarg;
            break;

        case 'd':
            this->distance_name = optarg;
            break;
        case 'r':
            this->r_index_name = optarg;
            break;
        case 'i':
            this->haplotype_input = optarg;
            break;
        case 'k':
            this->kmer_input = optarg;
            break;

        case OPT_KMER_LENGTH:
            this->k = parse<size_t>(optarg);
            if (this->k == 0 || this->k > gbwtgraph::Key64::KMER_MAX_LENGTH) {
                std::cerr << "error: [vg haplotypes] kmer length must be between 1 and " << gbwtgraph::Key64::KMER_MAX_LENGTH << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
        case OPT_WINDOW_LENGTH:
            this->w = parse<size_t>(optarg);
            if (this->w == 0) {
                std::cerr << "error: [vg haplotypes] window length cannot be 0" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
        case OPT_SUBCHAIN_LENGTH:
            this->partitioner_parameters.subchain_length = parse<size_t>(optarg);
            if (this->partitioner_parameters.subchain_length == 0) {
                std::cerr << "error: [vg haplotypes] subchain length cannot be 0" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
        case OPT_COVERAGE:
            this->recombinator_parameters.coverage = parse<size_t>(optarg);
            break;
        case OPT_NUM_HAPLOTYPES:
            this->recombinator_parameters.num_haplotypes = parse<size_t>(optarg);
            num_haplotypes_set = true;
            if (this->recombinator_parameters.num_haplotypes == 0) {
                std::cerr << "error: [vg haplotypes] number of haplotypes cannot be 0" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
        case OPT_PRESENT_DISCOUNT:
            this->recombinator_parameters.present_discount = parse<double>(optarg);
            if (this->recombinator_parameters.present_discount < 0.0 || this->recombinator_parameters.present_discount > 1.0) {
                std::cerr << "error: [vg haplotypes] present discount must be between 0.0 and 1.0" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
        case OPT_HET_ADJUSTMENT:
            this->recombinator_parameters.het_adjustment = parse<double>(optarg);
            if (this->recombinator_parameters.het_adjustment < 0.0) {
                std::cerr << "error: [vg haplotypes] het adjustment must be non-negative" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
        case OPT_ABSENT_SCORE:
            this->recombinator_parameters.absent_score = parse<double>(optarg);
            if (this->recombinator_parameters.absent_score < 0.0) {
                std::cerr << "error: [vg haplotypes] absent score must be non-negative" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            break;
        case OPT_DIPLOID_SAMPLING:
            this->recombinator_parameters.diploid_sampling = true;
            break;
        case OPT_INCLUDE_REFERENCE:
            this->recombinator_parameters.include_reference = true;
            break;

        case 'v':
            {
                size_t level = parse<size_t>(optarg);
                if (level > Haplotypes::verbosity_debug) {
                    std::cerr << "error: [vg haplotypes] invalid verbosity level: " << level << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                this->verbosity = static_cast<HaplotypePartitioner::Verbosity>(level);
            }
            break;
        case 't':
            this->threads = parse<size_t>(optarg);
            if (this->threads == 0 || this->threads > max_threads) {
                std::cerr << "error: [vg haplotypes] cannot run " << this->threads << " threads in parallel on this system" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            omp_set_num_threads(this->threads);
            break;

        case OPT_VALIDATE:
            this->validate = true;
            break;
        case OPT_VCF_INPUT:
            this->vcf_input = optarg;
            break;
        case OPT_CONTIG_PREFIX:
            this->contig_prefix = optarg;
            break;
        case OPT_EXTRACT:
            {
                std::string arg = optarg;
                size_t offset = arg.find(':');
                if (offset == 0 || offset == std::string::npos || offset + 1 >= arg.length()) {
                    std::cerr << "error: [vg haplotypes] cannot parse chain:subchain from " << arg << std::endl;
                    std::exit(EXIT_FAILURE);
                }
                this->chain_id = parse<size_t>(arg.substr(0, offset));
                this->subchain_id = parse<size_t>(arg.substr(offset + 1));
            }
            break;
        case OPT_SCORE_OUTPUT:
            this->score_output = optarg;
            break;
        case OPT_CLASSIFY:
            this->kmer_output = optarg;
            break;

        case 'h':
        case '?':
            help_haplotypes(argv, true);
            std::exit(EXIT_FAILURE);
        default:
            std::abort();
        }
    }

    // Determine input graph and set operating mode.
    if (optind + 1 != argc) {
        help_haplotypes(argv, false);
        std::exit(EXIT_FAILURE);
    }
    this->graph_name = argv[optind];
    if (this->haplotype_input.empty() && !this->kmer_input.empty() && !this->gbz_output.empty()) {
        this->mode = mode_sample_graph;
    } else if (this->haplotype_input.empty() && !this->haplotype_output.empty()) {
        this->mode = mode_preprocess;
    } else if (!this->haplotype_input.empty() && !this->kmer_input.empty() && !this->gbz_output.empty()) {
        this->mode = mode_sample_haplotypes;
    } else if (!this->haplotype_input.empty() && !this->vcf_input.empty()) {
        this->mode = mode_map_variants;
    } else if (!this->haplotype_input.empty() && !this->kmer_input.empty() &&
        this->chain_id < std::numeric_limits<size_t>::max() && this->subchain_id < std::numeric_limits<size_t>::max()) {
        this->mode = mode_extract;
    } else if (!this->haplotype_input.empty() && !this->kmer_input.empty() && !this->kmer_output.empty()) {
        this->mode = mode_classify;
    }
    if (this->mode == mode_invalid) {
        help_haplotypes(argv, false);
        std::exit(EXIT_FAILURE);
    }

    // Use conditional defaults if the user did not override them.
    if (this->recombinator_parameters.diploid_sampling && !num_haplotypes_set) {
        this->recombinator_parameters.num_haplotypes = haplotypes_default_candidates();
    }
}

//----------------------------------------------------------------------------

void validate_haplotypes(const Haplotypes& haplotypes,
                         const gbwtgraph::GBWTGraph& graph,
                         const gbwt::FastLocate& r_index,
                         const HaplotypePartitioner::minimizer_index_type& minimizer_index,
                         size_t expected_chains,
                         HaplotypePartitioner::Verbosity verbosity);

bool ends_with(const std::string& str, const std::string& suffix) {
    if (str.length() < suffix.length()) {
        return false;
    }
    return (str.substr(str.length() - suffix.length()) == suffix);
}

std::string get_name(const std::string& graph_name, const std::string& extension) {
    size_t length = graph_name.length();
    if (ends_with(graph_name, gbwtgraph::GBZ::EXTENSION)) {
        length -= gbwtgraph::GBZ::EXTENSION.length();
    }
    return graph_name.substr(0, length) + extension;
}

void preprocess_graph(const gbwtgraph::GBZ& gbz, Haplotypes& haplotypes, HaplotypesConfig& config) {
    double start = gbwt::readTimer();
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Generating haplotype information" << std::endl;
    }

    // Distance index.
    if (config.distance_name.empty()) {
        config.distance_name = get_name(config.graph_name, ".dist");
        if (config.verbosity >= Haplotypes::verbosity_basic) {
            std::cerr << "Guessing that distance index is " << config.distance_name << std::endl;
        }
    }
    SnarlDistanceIndex distance_index;
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Loading distance index from " << config.distance_name << std::endl;
    }
    distance_index.deserialize(config.distance_name);
    size_t expected_chains = 0;
    distance_index.for_each_child(distance_index.get_root(), [&](const handlegraph::net_handle_t&) {
        expected_chains++;
    });

    // Minimizer index.
    HaplotypePartitioner::minimizer_index_type minimizer_index(config.k, config.w, false);
    {
        double minimizer = gbwt::readTimer();
        if (config.verbosity >= Haplotypes::verbosity_basic) {
            std::cerr << "Building minimizer index" << std::endl;
        }
        gbwtgraph::index_haplotypes(gbz.graph, minimizer_index);
        if (config.verbosity >= Haplotypes::verbosity_basic) {
            double seconds = gbwt::readTimer() - minimizer;
            std::cerr << "Built the minimizer index in " << seconds << " seconds" << std::endl;
        }
    }

    // R-index.
    if (config.r_index_name.empty()) {
        config.r_index_name = get_name(config.graph_name, gbwt::FastLocate::EXTENSION);
        if (config.verbosity >= Haplotypes::verbosity_basic) {
            std::cerr << "Guessing that r-index is " << config.r_index_name << std::endl;
        }
    }
    gbwt::FastLocate r_index;
    load_r_index(r_index, config.r_index_name, config.verbosity >= Haplotypes::verbosity_basic);
    r_index.setGBWT(gbz.index);

    // Partition the haplotypes.
    HaplotypePartitioner partitioner(gbz, r_index, distance_index, minimizer_index, config.verbosity);
    try {
        haplotypes = partitioner.partition_haplotypes(config.partitioner_parameters);
    }
    catch (const std::runtime_error& e) {
        std::cerr << "error: [vg haplotypes] " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Generated haplotype information in " << seconds << " seconds" << std::endl;
    }

    // Validate the haplotypes.
    if (config.validate) {
        validate_haplotypes(haplotypes, gbz.graph, r_index, minimizer_index, expected_chains, config.verbosity);
    }
}

//----------------------------------------------------------------------------

size_t threads_to_jobs(size_t threads) {
    size_t jobs = std::round(0.85 * threads);
    return std::max(jobs, size_t(1));
}

void validate_subgraph(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph, HaplotypePartitioner::Verbosity verbosity);

void sample_haplotypes(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config) {
    omp_set_num_threads(threads_to_jobs(config.threads));
    Recombinator recombinator(gbz, config.verbosity);
    gbwt::GBWT merged;
    try {
        merged = recombinator.generate_haplotypes(haplotypes, config.kmer_input, config.recombinator_parameters);
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [vg haplotypes] " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    omp_set_num_threads(config.threads); // Restore the number of threads.

    // Build and serialize GBWTGraph.
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Building GBWTGraph" << std::endl;
    }
    double checkpoint = gbwt::readTimer();
    gbwtgraph::GBWTGraph output_graph = gbz.graph.subgraph(merged);
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - checkpoint;
        std::cerr << "Built the GBWTGraph in " << seconds << " seconds" << std::endl;
    }
    save_gbz(merged, output_graph, config.gbz_output, config.verbosity >= Haplotypes::verbosity_basic);

    // Validate the graph.
    if (config.validate) {
        // TODO: How could we validate the haplotypes?
        validate_subgraph(gbz.graph, output_graph, config.verbosity);
    }
}

//----------------------------------------------------------------------------

gbwt::size_type path_for_contig(const gbwtgraph::GBZ& gbz, gbwt::size_type contig_id, const std::string& contig_name) {
    gbwt::size_type path_id = gbz.index.metadata.paths();
    size_t found_paths = 0;
    for (size_t i = 0; i < gbz.graph.named_paths.size(); i++) {
        gbwt::size_type candidate = gbz.graph.named_paths[i].id;
        if (gbz.index.metadata.path(candidate).contig == contig_id) {
            path_id = candidate;
            found_paths++;
        }
    }
    if (found_paths != 1) {
        std::cerr << "error: [vg haplotypes] found " << found_paths << " named/reference paths for contig " << contig_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return path_id;
}

std::pair<gbwt::size_type, size_t> seq_chain_for_path(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, gbwt::size_type path_id, const std::string& contig_name) {
    gbwt::size_type sequence_id = gbwt::Path::encode(path_id, false);
    gbwt::size_type reverse_id = gbwt::Path::encode(path_id, true);
    size_t found_chains = 0;
    std::pair<gbwt::size_type, size_t> result(gbwt::invalid_sequence(), haplotypes.components());
    for (size_t chain_id = 0; chain_id < haplotypes.components(); chain_id++) {
        const Haplotypes::Subchain& subchain = haplotypes.chains[chain_id].subchains.front();
        for (size_t i = 0; i < subchain.sequences.size(); i++) {
            if (subchain.sequences[i].first == sequence_id) {
                result.first = sequence_id;
                result.second = chain_id;
                found_chains++;
                break;
            }
            if (subchain.sequences[i].first == reverse_id) {
                result.first = reverse_id;
                result.second = chain_id;
                found_chains++;
                break;
            }
        }
    }
    if (found_chains != 1) {
        std::cerr << "error: [vg haplotypes] found " << found_chains << " top-level chains for contig " << contig_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return result;
}

struct ReferenceInterval {
    enum order { before, overlap, after };

    Haplotypes::Subchain::subchain_t type;

    size_t id;

    // Semiopen range of reference positions for the internal parts of the subchain.
    size_t start, end;

    // Where is this interval relative to the specified interval?
    order compare(std::pair<size_t, size_t> interval) {
        if (this->end <= interval.first) {
            return before;
        } else if (this->start >= interval.second) {
            return after;
        } else {
            return overlap;
        }
    }

    size_t length() const {
        return this->end - this->start;
    }

    std::string to_string() const {
        std::string result;
        switch (this->type) {
            case Haplotypes::Subchain::normal:
                result.push_back('N');
                break;
            case Haplotypes::Subchain::prefix:
                result.push_back('P');
                break;
            case Haplotypes::Subchain::suffix:
                result.push_back('S');
                break;
            case Haplotypes::Subchain::full_haplotype:
                result.push_back('F');
                break;
        }
        result += std::to_string(this->id) + "(" + std::to_string(this->start) + ".." + std::to_string(this->end) + ")";
        return result;
    }
};

std::vector<ReferenceInterval> subchain_intervals(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, gbwt::size_type sequence_id, size_t chain_id, bool reverse) {
    gbwt::size_type actual_sequence_id = (reverse ? gbwt::Path::reverse(sequence_id) : sequence_id);
    gbwt::vector_type path = gbz.index.extract(actual_sequence_id);
    size_t total_length = 0;
    for (auto gbwt_node : path) {
        total_length += gbz.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(gbwt_node));
    }

    const Haplotypes::TopLevelChain& chain = haplotypes.chains[chain_id];
    std::vector<ReferenceInterval> result;
    size_t seq_offset = 0, node_offset = 0;
    for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
        size_t actual_subchain_id;
        Haplotypes::Subchain subchain;
        if (reverse) {
            actual_subchain_id = chain.subchains.size() - 1 - subchain_id;
            switch (chain.subchains[actual_subchain_id].type) {
                case Haplotypes::Subchain::prefix:
                    subchain.type = Haplotypes::Subchain::suffix;
                    break;
                case Haplotypes::Subchain::suffix:
                    subchain.type = Haplotypes::Subchain::prefix;
                    break;
                default:
                    subchain.type = chain.subchains[actual_subchain_id].type;
                    break;
            }
            subchain.start = gbwt::Node::reverse(chain.subchains[actual_subchain_id].end);
            subchain.end = gbwt::Node::reverse(chain.subchains[actual_subchain_id].start);
        } else {
            actual_subchain_id = subchain_id;
            subchain.type = chain.subchains[actual_subchain_id].type;
            subchain.start = chain.subchains[actual_subchain_id].start;
            subchain.end = chain.subchains[actual_subchain_id].end;
        }
        ReferenceInterval interval { subchain.type, actual_subchain_id, 0, total_length };
        if (subchain.has_start()) {
            while (node_offset < path.size() && path[node_offset] != subchain.start) {
                seq_offset += gbz.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(path[node_offset]));
                node_offset++;
            }
            if (node_offset < path.size()) {
                seq_offset += gbz.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(path[node_offset]));
                node_offset++;
            }
            interval.start = seq_offset;
        } else if (subchain.type == Haplotypes::Subchain::prefix) {
            // If a prefix follows a suffix, they cover the same interval.
            interval.start = result.back().start;
        }
        if (subchain.has_end()) {
            while (node_offset < path.size() && path[node_offset] != subchain.end) {
                seq_offset += gbz.graph.get_length(gbwtgraph::GBWTGraph::node_to_handle(path[node_offset]));
                node_offset++;
            }
            interval.end = seq_offset;
            // If a prefix follows a suffix, they cover the same interval.
            if (subchain.type == Haplotypes::Subchain::prefix) {
                result.back().end = interval.end;
            }
        }
        result.push_back(interval);
    }

    return result;
}

void map_variants(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config) {
    if (!gbz.index.metadata.hasContigNames()) {
        std::cerr << "error: [vg haplotypes] cannot map variant positions without contig names in the GBWT index" << std::endl;
    }

    // Read variants from the VCF file.
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Reading VCF file " << config.vcf_input << std::endl;
    }
    vcflib::VariantCallFile variant_file;
    variant_file.parseSamples = false; // Just in case there are many samples.
    std::string temp_filename = config.vcf_input;
    variant_file.open(temp_filename);
    if (!variant_file.is_open()) {
        std::cerr << "error: [vg haplotypes] cannot open VCF file " << config.vcf_input << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::unordered_map<std::string, size_t> contig_to_offset; // VCF contig name to offset in `variant positions`.
    std::vector<std::vector<std::pair<size_t, size_t>>> variant_positions; // Semiopen 0-based ranges of sequence positions.
    vcflib::Variant var(variant_file);
    size_t total_variants = 0;
    while (variant_file.is_open() && variant_file.getNextVariant(var)) {
        size_t offset;
        auto iter = contig_to_offset.find(var.sequenceName);
        if (iter == contig_to_offset.end()) {
            offset = variant_positions.size();
            contig_to_offset[var.sequenceName] = offset;
            variant_positions.push_back({});
        } else {
            offset = iter->second;
        }
        size_t start = var.zeroBasedPosition();
        variant_positions[offset].push_back({ start, start + var.ref.length() });
        total_variants++;
    }
    for (auto& positions : variant_positions) {
        std::sort(positions.begin(), positions.end());
    }
    if (config.verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "Read " << total_variants << " variants over " << variant_positions.size() << " contigs" << std::endl;
    }

    // Map VCF contig names to GBWT sequence ids for named/reference paths and top-level chain.
    std::vector<std::string> contig_names(contig_to_offset.size(), "");
    std::vector<std::pair<gbwt::size_type, size_t>> offset_to_seq_chain(contig_to_offset.size(), { gbwt::invalid_sequence(), haplotypes.components() });
    for (auto iter = contig_to_offset.begin(); iter != contig_to_offset.end(); ++iter) {
        std::string contig_name = config.contig_prefix + iter->first;
        gbwt::size_type contig_id = gbz.index.metadata.contig(contig_name);
        if (contig_id >= gbz.index.metadata.contigs()) {
            std::cerr << "error: [vg haplotypes] no contig " << contig_name << " in the GBWT index" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        contig_names[iter->second] = contig_name;
        gbwt::size_type path_id = path_for_contig(gbz, contig_id, contig_name);
        std::pair<gbwt::size_type, size_t> seq_chain = seq_chain_for_path(gbz, haplotypes, path_id, contig_name);
        offset_to_seq_chain[iter->second] = seq_chain;
        if (config.verbosity >= Haplotypes::verbosity_debug) {
            std::cerr << "VCF contig " << iter->first << ", GBWT contig " << contig_name
                << ": contig id " << contig_id
                << ", path id " << path_id
                << ", reverse " << gbwt::Path::is_reverse(seq_chain.first)
                << ", chain " << seq_chain.second << std::endl;
        }
    }

    // Output (contig[interval], top-level chain, subchains, subchain lengths)
    for (auto iter = contig_to_offset.begin(); iter != contig_to_offset.end(); ++iter) {
        std::string contig_name = config.contig_prefix + iter->first;
        size_t offset = iter->second;
        gbwt::size_type sequence_id = offset_to_seq_chain[offset].first;
        gbwt::size_type chain_id = offset_to_seq_chain[offset].second;
        auto ref_intervals = subchain_intervals(gbz, haplotypes, sequence_id, chain_id, gbwt::Path::is_reverse(sequence_id));
        for (auto interval : variant_positions[offset]) {
            size_t low = 0, high = ref_intervals.size();
            bool found = false;
            while (!found && low < high) {
                size_t mid = low + (high - low) / 2;
                switch (ref_intervals[mid].compare(interval)) {
                    case ReferenceInterval::before:
                        low = mid + 1;
                        break;
                    case ReferenceInterval::overlap:
                        low = mid;
                        while (low > 0 && ref_intervals[low - 1].compare(interval) == ReferenceInterval::overlap) {
                            low--;
                        }
                        high = mid + 1;
                        while (high < ref_intervals.size() && ref_intervals[high].compare(interval) == ReferenceInterval::overlap) {
                            high++;
                        }
                        found = true;
                        break;
                    case ReferenceInterval::after:
                        high = mid;
                        break;
                }
            }
            std::cout << iter->first << "[" << interval.first << ".." << interval.second << "]\t" << chain_id << "\t";
            if (low >= high) {
                if (low > 0) {
                    std::cout << ref_intervals[low - 1].to_string();
                }
                std::cout << "..";
                if (low < ref_intervals.size()) {
                    std::cout << ref_intervals[low].to_string();
                }
            } else {
                for (size_t i = low; i < high; i++) {
                    if (i > low) {
                        std::cout << ",";
                    }
                    std::cout << ref_intervals[i].to_string();
                }
            }
            std::cout << "\t";
            for (size_t i = low; i < high; i++) {
                if (i > low) {
                    std::cout << ",";
                }
                std::cout << ref_intervals[i].length();
            }
            std::cout << std::endl;
        }
    }
}
//----------------------------------------------------------------------------

void extract_haplotypes(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config) {
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Extracting haplotypes from chain " << config.chain_id << ", subchain " << config.subchain_id << std::endl;
    }

    Recombinator recombinator(gbz, config.verbosity);
    std::vector<Recombinator::LocalHaplotype> result;
    try {
        result = recombinator.extract_sequences(
            haplotypes, config.kmer_input,
            config.chain_id, config.subchain_id, config.recombinator_parameters
        );
    } catch (const std::runtime_error& e) {
        std::cerr << "error: [vg haplotypes] " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (config.verbosity >= Haplotypes::verbosity_detailed) {
        std::cerr << "Found " << result.size() << " haplotypes" << std::endl;
    }
    for (auto& sequence : result) {
        write_fasta_sequence(sequence.name, sequence.sequence, std::cout);
    }

    if (!config.score_output.empty()) {
        std::ofstream out(config.score_output, std::ios_base::binary);
        if (!out) {
            std::cerr << "error: [vg haplotypes] cannot open score file " << config.score_output << " for writing" << std::endl;
            return;
        }
        for (auto& sequence : result) {
            out << sequence.name;
            for (size_t i = 0; i < config.recombinator_parameters.num_haplotypes; i++) {
                if (i < sequence.scores.size()) {
                    out << "\t" << sequence.scores[i].first << "\t" << sequence.scores[i].second;
                } else {
                    out << "\t-\t-";
                }
            }
            out << "\n";
        }
    }
}

//----------------------------------------------------------------------------

void classify_kmers(const gbwtgraph::GBZ& gbz, const Haplotypes& haplotypes, const HaplotypesConfig& config) {
    if (config.verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Classifying kmers" << std::endl;
    }
    Recombinator recombinator(gbz, config.verbosity);
    std::vector<char> classifications = recombinator.classify_kmers(haplotypes, config.kmer_input, config.recombinator_parameters);

    if (config.verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Writing " << classifications.size() << " classifications to " << config.kmer_output << std::endl;
    }
    std::ofstream out(config.kmer_output, std::ios_base::binary);
    if (!out) {
        std::cerr << "error: [vg haplotypes] cannot open kmer classification file " << config.kmer_output << " for writing" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // Multi-gigabyte writes do not work in every environment, but we assume that
    // there are only at most a few hundred million kmers.
    out.write(classifications.data(), classifications.size());
}

//----------------------------------------------------------------------------

void validate_error(const std::string& header, const std::string& message) {
    std::cerr << "error: [vg haplotypes] ";
    if (!header.empty()) {
        std::cerr << header << ": ";
    }
    std::cerr << message << std::endl;
    std::exit(EXIT_FAILURE);
}

template<typename T>
std::string expected_got(T expected, T got) {
    return "expected " + std::to_string(expected) + ", got " + std::to_string(got);
}

template<typename T>
std::string pair_to_string(std::pair<T, T> value) {
    return "(" + std::to_string(value.first) + ", " + std::to_string(value.second) + ")";
}

void validate_error_chain(size_t chain_id, const std::string& message) {
    validate_error("chain " + std::to_string(chain_id), message);
}

void validate_error_subchain(size_t chain_id, size_t subchain_id, const std::string& message) {
    validate_error("chain " + std::to_string(chain_id) + ", subchain " + std::to_string(subchain_id), message);
}

void validate_error_sequence(size_t chain_id, size_t subchain_id, size_t sequence_id, const std::string& message) {
    std::string header = "chain " + std::to_string(chain_id) + ", subchain " + std::to_string(subchain_id) + ", sequence " + std::to_string(sequence_id);
    validate_error(header, message);
}

std::string validate_unary_path(const HandleGraph& graph, handle_t from, handle_t to) {
    hash_set<handle_t> visited;
    handle_t curr = from;
    while (curr != to) {
        if (visited.find(curr) != visited.end()) {
            return "incoming path contains a cycle";
        }
        visited.insert(curr);
        handle_t successor = empty_gbwtgraph_handle();
        size_t successors = 0;
        graph.follow_edges(curr, false, [&](const handle_t& next) {
            successor = next;
            successors++;
        });
        if (successors != 1) {
            return "incoming path is not unary";
        }
        curr = successor;
    }
    return "";
}

// Returns true if the path from (start, offset) reaches end without revisiting start.
bool trace_path(const gbwt::GBWT& index, gbwt::node_type start, gbwt::size_type offset, gbwt::node_type end) {
    gbwt::edge_type pos(start, offset);
    while (pos.first != end) {
        pos = index.LF(pos);
        if (pos.first == gbwt::ENDMARKER || pos.first == start) {
            return false;
        }
    }
    return true;
}

// Returns the given haplotype over the given subchain.
std::string get_haplotype(const gbwtgraph::GBWTGraph& graph, Haplotypes::sequence_type sequence,
                          gbwt::node_type from, gbwt::node_type to, size_t k) {
    std::string haplotype;
    gbwt::edge_type pos;

    // Initial node with three cases (from start, suffix of a long `from`, short `from`).
    if (from == gbwt::ENDMARKER) {
        pos = graph.index->start(sequence.first);
        gbwtgraph::view_type view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(pos.first));
        haplotype.append(view.first, view.second);
    } else {
        pos = gbwt::edge_type(from, sequence.second);
        gbwtgraph::view_type view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(pos.first));
        if (view.second >= k) {
            haplotype.append(view.first + view.second - (k - 1), k - 1);
        } else {
            haplotype.append(view.first, view.second);
        }
    }

    while (true) {
        pos = graph.index->LF(pos);
        if (pos.first == gbwt::ENDMARKER) {
            break;
        }
        gbwtgraph::view_type view = graph.get_sequence_view(gbwtgraph::GBWTGraph::node_to_handle(pos.first));
        if (pos.first == to) {
            haplotype.append(view.first, std::min(view.second, k - 1));
            break;
        } else {
            haplotype.append(view.first, view.second);
        }
    }

    return haplotype;
}

void validate_chain(const Haplotypes::TopLevelChain& chain,
                    const gbwtgraph::GBWTGraph& graph,
                    const gbwt::FastLocate& r_index,
                    const HaplotypePartitioner::minimizer_index_type& minimizer_index,
                    size_t chain_id,
                    HaplotypePartitioner::Verbosity verbosity) {
    if (chain.offset != chain_id) {
        validate_error_chain(chain_id, "stored id is " + std::to_string(chain.offset));
    }
    if (chain.subchains.empty()) {
        validate_error_chain(chain_id, "the chain is empty");
    }

    const Haplotypes::Subchain* prev = nullptr;
    for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
        const Haplotypes::Subchain& subchain = chain.subchains[subchain_id];

        // Check that the subchain is of an appropriate type.
        switch (subchain.type) {
        case Haplotypes::Subchain::normal:
            break;
        case Haplotypes::Subchain::prefix:
            if (subchain_id > 0 && prev->type != Haplotypes::Subchain::suffix) {
                validate_error_subchain(chain_id, subchain_id, "a prefix inside a fragment");
            }
            break;
        case Haplotypes::Subchain::suffix:
            break;
        case Haplotypes::Subchain::full_haplotype:
            if (chain.subchains.size() != 1) {
                validate_error_subchain(chain_id, subchain_id, "full haplotypes in a nontrivial chain");
            }
            break;
        }

        // Check that the boundary nodes have been defined.
        if (subchain.has_start() && subchain.start == gbwt::ENDMARKER) {
            validate_error_subchain(chain_id, subchain_id, "missing start node");
        }
        if (subchain.has_end() && subchain.end == gbwt::ENDMARKER) {
            validate_error_subchain(chain_id, subchain_id, "missing end node");
        }

        // Check that the kmer presence bitvector is of appropriate length.
        size_t total_kmers = subchain.sequences.size() * subchain.kmers.size();
        if (subchain.kmers_present.size() != total_kmers) {
            std::string message = expected_got<size_t>(total_kmers, subchain.kmers_present.size()) + " kmer occurrences";
            validate_error_subchain(chain_id, subchain_id, message);
        }

        // Check that there is a unary path from the previous subchain if the
        // appropriate boundary nodes are present.
        if (subchain_id > 0 && prev->has_end() && subchain.has_start()) {
            std::string message = validate_unary_path(graph, gbwtgraph::GBWTGraph::node_to_handle(prev->end), gbwtgraph::GBWTGraph::node_to_handle(subchain.start));
            if (!message.empty()) {
                validate_error_subchain(chain_id, subchain_id, message);
            }
        }

        // Sequences: normal subchains.
        if (subchain.type == Haplotypes::Subchain::normal) {
            std::vector<gbwt::size_type> da = r_index.decompressDA(subchain.start);
            hash_set<Haplotypes::sequence_type> selected;
            for (size_t i = 0; i < da.size(); i++) {
                if (trace_path(*(graph.index), subchain.start, i, subchain.end)) {
                    selected.insert(Haplotypes::sequence_type(da[i], i));
                }
            }
            if (subchain.sequences.size() != selected.size()) {
                std::string message = expected_got(selected.size(), subchain.sequences.size()) + " sequences (normal)";
                validate_error_subchain(chain_id, subchain_id, message);
            }
            for (size_t i = 0; i < subchain.sequences.size(); i++) {
                if (selected.find(subchain.sequences[i]) == selected.end()) {
                    std::string message = "invalid value " + pair_to_string(subchain.sequences[i]);
                    validate_error_sequence(chain_id, subchain_id, i, message);
                }
            }
        }

        // Sequences: prefixes and suffixes.
        if (subchain.type == Haplotypes::Subchain::prefix || subchain.type == Haplotypes::Subchain::suffix) {
            gbwt::node_type node = (subchain.has_start() ? subchain.start : subchain.end);
            std::vector<gbwt::size_type> da = r_index.decompressDA(node);
            if (subchain.sequences.size() != da.size()) {
                std::string message = expected_got(da.size(), subchain.sequences.size()) + " sequences (prefix / suffix)";
                validate_error_subchain(chain_id, subchain_id, message);
            }
            hash_set<Haplotypes::sequence_type> truth;
            for (size_t i = 0; i < da.size(); i++) {
                truth.insert({ da[i], i });
            }
            for (size_t i = 0; i < subchain.sequences.size(); i++) {
                if (truth.find(subchain.sequences[i]) == truth.end()) {
                    std::string message = "invalid value " + pair_to_string(subchain.sequences[i]);
                    validate_error_sequence(chain_id, subchain_id, i, message);
                }
            }
        }

        // Sequences: full haplotypes.
        if (subchain.type == Haplotypes::Subchain::full_haplotype) {
            if (subchain.sequences.empty()) {
                validate_error_subchain(chain_id, subchain_id, "full haplotypes without sequences");
            }
        }

        // Kmers.
        if (subchain.type != Haplotypes::Subchain::full_haplotype) {
            hash_set<Haplotypes::Subchain::kmer_type> all_kmers;
            for (size_t i = 0; i < subchain.kmers.size(); i++) {
                all_kmers.insert(subchain.kmers[i].first);
            }
            if (all_kmers.size() != subchain.kmers.size()) {
                std::string message = expected_got(subchain.kmers.size(), all_kmers.size()) + " kmers";
                validate_error_subchain(chain_id, subchain_id, message);
            }
            hash_map<Haplotypes::Subchain::kmer_type, size_t> used_kmers; // (kmer used in haplotypes, number of sequences that contain it)
            hash_map<Haplotypes::Subchain::kmer_type, size_t> missing_kmers; // (kmer not used in haplotypes, number of sequences that contain it)
            for (size_t i = 0; i < subchain.sequences.size(); i++) {
                std::string haplotype = get_haplotype(graph, subchain.sequences[i], subchain.start, subchain.end, minimizer_index.k());
                auto minimizers = minimizer_index.minimizers(haplotype);
                hash_map<Haplotypes::Subchain::kmer_type, bool> unique_minimizers; // (kmer, used in the sequence)
                for (auto& minimizer : minimizers) {
                    if (minimizer_index.count(minimizer) == 1) {
                        unique_minimizers[minimizer.key.get_key()] = false;
                    }
                }
                for (size_t j = 0, offset = i * subchain.kmers.size(); j < subchain.kmers.size(); j++, offset++) {
                    if (subchain.kmers_present[offset]) {
                        auto iter = unique_minimizers.find(subchain.kmers[j].first);
                        if (iter == unique_minimizers.end()) {
                            std::string message = "kmer " + std::to_string(j) + " not present in the haplotype";
                            validate_error_sequence(chain_id, subchain_id, i, message);
                        }
                        used_kmers[subchain.kmers[j].first]++;
                        iter->second = true;
                    } else {
                        if (unique_minimizers.find(subchain.kmers[j].first) != unique_minimizers.end()) {
                            std::string message = "kmer " + std::to_string(j) + " is present in the haplotype";
                            validate_error_sequence(chain_id, subchain_id, i, message);
                        }
                    }
                }
                for (auto iter = unique_minimizers.begin(); iter != unique_minimizers.end(); ++iter) {
                    if (!iter->second) {
                        missing_kmers[iter->first]++;
                    }
                }
            }
            size_t invalid_count = 0;
            for (size_t kmer_id = 0; kmer_id < subchain.kmers.size(); kmer_id++) {
                size_t count = 0;
                auto iter = used_kmers.find(subchain.kmers[kmer_id].first);
                if (iter == used_kmers.end() || iter->second != subchain.kmers[kmer_id].second) {
                    invalid_count++;
                }
            }
            if (invalid_count > 0) {
                std::string message = "invalid occurrence count for "+ std::to_string(invalid_count) + " kmers";
                validate_error_subchain(chain_id, subchain_id, message);
            }
            size_t missing_informative_kmers = 0;
            for (auto iter = missing_kmers.begin(); iter != missing_kmers.end(); ++iter) {
                if (iter->second < subchain.sequences.size()) {
                    missing_informative_kmers++;
                }
            }
            if (missing_informative_kmers > 0) {
                std::string message = "missing " + std::to_string(missing_informative_kmers) + " informative kmers";
                validate_error_subchain(chain_id, subchain_id, message);
            }
        }

        prev = &subchain;
    }
}

std::string subchain_to_string(size_t chain_id, size_t subchain_id, const Haplotypes::Subchain& subchain) {
    return "chain " + std::to_string(chain_id) + ", subchain " + std::to_string(subchain_id) + " (" + subchain.to_string() + ")";
}

void validate_haplotypes(const Haplotypes& haplotypes,
                         const gbwtgraph::GBWTGraph& graph,
                         const gbwt::FastLocate& r_index,
                         const HaplotypePartitioner::minimizer_index_type& minimizer_index,
                         size_t expected_chains,
                         HaplotypePartitioner::Verbosity verbosity) {
    if (verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Validating the haplotype information" << std::endl;
    }
    double start = gbwt::readTimer();

    // Header information.
    if (haplotypes.k() != minimizer_index.k()) {
        validate_error("k-mer length", expected_got(minimizer_index.k(), haplotypes.k()));
    }
    if (haplotypes.components() != expected_chains) {
        validate_error("graph components", expected_got(expected_chains, haplotypes.components()));
    }
    if (haplotypes.components() != haplotypes.chains.size()) {
        validate_error("top-level chains", expected_got(haplotypes.components(), haplotypes.chains.size()));
    }
    std::vector<size_t> chains_per_job(haplotypes.jobs(), 0);
    for (size_t chain = 0; chain < haplotypes.components(); chain++) {
        size_t job_id = haplotypes.chains[chain].job_id;
        if (job_id >= haplotypes.jobs()) {
            validate_error_chain(chain, "job id " + std::to_string(job_id) + " >= " + std::to_string(haplotypes.jobs()));
        }
        chains_per_job[job_id]++;
    }
    for (size_t job_id = 0; job_id < chains_per_job.size(); job_id++) {
        if (chains_per_job[job_id] == 0) {
            validate_error("", "job " + std::to_string(job_id) + " is empty");
        }
    }

    // Cached paths.
    if (haplotypes.jobs_for_cached_paths.size() != graph.named_paths.size()) {
        validate_error("cached paths", expected_got(graph.named_paths.size(), haplotypes.jobs_for_cached_paths.size()));
    }

    // Haplotype information is valid
    if (verbosity >= HaplotypePartitioner::Verbosity::verbosity_detailed) {
        std::cerr << "Validating subchains, sequences, and kmers" << std::endl;
    }
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t chain = 0; chain < haplotypes.components(); chain++) {
        validate_chain(haplotypes.chains[chain], graph, r_index, minimizer_index, chain, verbosity);
    }

    // Kmers are globally unique.
    if (verbosity >= HaplotypePartitioner::Verbosity::verbosity_detailed) {
        std::cerr << "Validating kmer specificity" << std::endl;
    }
    hash_map<Haplotypes::Subchain::kmer_type, std::pair<size_t, size_t>> kmers;
    for (size_t chain_id = 0; chain_id < haplotypes.components(); chain_id++) {
        const Haplotypes::TopLevelChain& chain = haplotypes.chains[chain_id];
        for (size_t subchain_id = 0; subchain_id < chain.subchains.size(); subchain_id++) {
            const Haplotypes::Subchain& subchain = chain.subchains[subchain_id];
            for (size_t i = 0; i < subchain.kmers.size(); i++) {
                auto iter = kmers.find(subchain.kmers[i].first);
                if (iter != kmers.end()) {
                    const Haplotypes::Subchain& prev = haplotypes.chains[iter->second.first].subchains[iter->second.second];
                    if (chain_id == iter->second.first && subchain_id == iter->second.second + 1 && subchain.type == Haplotypes::Subchain::prefix && prev.type == Haplotypes::Subchain::suffix) {
                        // A prefix subchain may overlap the preceding suffix subchain and
                        // contain the same kmers.
                    } else {
                        std::string message = subchain.to_string() + ": kmer " + std::to_string(i) + " also found in " + subchain_to_string(iter->second.first, iter->second.second, prev);
                        validate_error_subchain(chain_id, subchain_id, message);
                    }
                }
                kmers[subchain.kmers[i].first] = { chain_id, subchain_id };
            }
        }
    }
    kmers.clear();

    if (verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Validated the haplotype information in " << seconds << " seconds" << std::endl;
    }
}

//----------------------------------------------------------------------------

void validate_nodes(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph) {
    nid_t last_node = 0;
    bool nodes_ok = subgraph.for_each_handle([&](const handle_t& handle) -> bool {
        last_node = subgraph.get_id(handle);
        return graph.has_node(last_node);
    });
    if (!nodes_ok) {
        validate_error("", "invalid node " + std::to_string(last_node));
    }
}

void validate_edges(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph) {
    edge_t last_edge(gbwtgraph::GBWTGraph::node_to_handle(0), gbwtgraph::GBWTGraph::node_to_handle(0));
    bool edges_ok = subgraph.for_each_edge([&](const edge_t& edge) -> bool {
        last_edge = edge;
        return graph.has_edge(edge.first, edge.second);
    });
    if (!edges_ok) {
        validate_error("", "invalid edge " + to_string_gbwtgraph(last_edge.first) + " to " + to_string_gbwtgraph(last_edge.second));
    }
}

void validate_subgraph(const gbwtgraph::GBWTGraph& graph, const gbwtgraph::GBWTGraph& subgraph, HaplotypePartitioner::Verbosity verbosity) {
    if (verbosity >= Haplotypes::verbosity_basic) {
        std::cerr << "Validating the output subgraph" << std::endl;
    }
    double start = gbwt::readTimer();

    std::thread nodes(validate_nodes, std::cref(graph), std::cref(subgraph));
    std::thread edges(validate_edges, std::cref(graph), std::cref(subgraph));
    nodes.join();
    edges.join();

    if (verbosity >= Haplotypes::verbosity_basic) {
        double seconds = gbwt::readTimer() - start;
        std::cerr << "Validated the subgraph in " << seconds << " seconds" << std::endl;
    }
}

//----------------------------------------------------------------------------

