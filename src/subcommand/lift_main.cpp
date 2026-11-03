/**
 * \file lift_main.cpp: lift features on a cover or on paths to the graph, the gref and the paths
 */

#include <getopt.h>
#include <omp.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <unordered_set>

#include <htslib/hts.h>
#include <htslib/kstring.h>

#include "subcommand.hpp"
#include "../feature_lift.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gref.hpp"
#include "../haplotype_depth.hpp"
#include "../nested_variants.hpp"
#include "../utility.hpp"
#include "../vcf_alleles.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_lift(char** argv) {
    cerr << "usage: " << argv[0] << " lift [options] graph.gbz features.bed > features.gaf" << endl
         << "       " << argv[0] << " lift [options] graph.gbz variants.vcf > variants.lifted.vcf" << endl
         << "       " << argv[0] << " lift --depth [options] graph.gbz [regions.bed] > depth.bedGraph" << endl
         << "Lift features to the graph and the gref, and to each path that visits them." << endl
         << "A feature can be on a cover sequence (a gref path of the GBZ or a tig of" << endl
         << "vg unitigs) or on any path in the GBZ, named without a subrange or phase" << endl
         << "block and in coordinates of the full path (for haplotypes, the assembly" << endl
         << "contig). A feature on several sequences at the same graph place is reported" << endl
         << "once." << endl
         << endl
         << "Each feature gets an ID (f1, f2, ...) and a GAF record with that query name," << endl
         << "its graph walk, and the tags sp:i: (visits of the whole walk by paths other" << endl
         << "than gref paths), gr:Z: (gref placement, see vg unitigs), pc:Z: (path" << endl
         << "coordinate) and an:Z: (BED name)." << endl
         << endl
         << "The path coordinate names a feature exactly, as PATH:START-END:STRAND on the" << endl
         << "first path, in priority order, that carries the whole feature: the gref path" << endl
         << "it lies along, else a reference or assembly path. Paths keep their" << endl
         << "coordinates in clipped, filtered and sampled graphs, and no two places in the" << endl
         << "graph share one. Naming a feature off the gref takes a pass over all paths." << endl
         << endl
         << "A VCF of variants in the graph gives a feature for each ALT allele: its walk" << endl
         << "without the first and last nodes, which lie on the reference on either side." << endl
         << "The walk comes from INFO/AT (vg deconstruct) or, for a record with one ALT" << endl
         << "allele and no AT, from a PanGenie variant ID (CHROM-POS-TYPE-WALK-LENGTH)." << endl
         << "The VCF is written out with each allele's gref placement in INFO/GR, its" << endl
         << "segments separated by |, and its path coordinate in INFO/PC. A deletion has" << endl
         << "no nodes left to place and gets neither. The VCF must be a file, plain or" << endl
         << "compressed. Its features are numbered in file order and named by the ID of" << endl
         << "their record, with :N added for the Nth ALT allele of a record with several." << endl
         << endl
         << "Where an allele runs along a gref alt contig and differs from it, such as an" << endl
         << "insertion with SNPs in it, each difference becomes a biallelic record on the" << endl
         << "contig, written after the calls as vg deconstruct writes nested records. A" << endl
         << "haplotype gets 1 if its called allele has the variant, 0 if the allele reads" << endl
         << "that part of the contig without it, and . otherwise, so allele frequencies" << endl
         << "are among the carriers of the contig. Its GQ is the lowest GQ of the calls" << endl
         << "that give the genotype, if they have GQ. With --gref-segs, each such record" << endl
         << "gets INFO/PS, the top-level snarl of its contig. A VCF that already has" << endl
         << "records on alt contigs gets none added." << endl
         << endl
         << "With --depth, write how many haplotypes read each base of the gref paths, or" << endl
         << "of the paths chosen with -Q, as bedGraph: the path, the 0-based start and end" << endl
         << "of a run of bases of equal depth, and the depth, which is the number of" << endl
         << "distinct haplotypes (sample and haplotype number) that read them. Runs of" << endl
         << "depth 0 are written too. Haplotype paths count, and reference paths with" << endl
         << "--references; gref paths and generic paths never do. A regions BED, on paths" << endl
         << "named and placed as features are, limits the count to the bases inside it," << endl
         << "and a path without regions counts nothing. Given one assembly's QC flags as" << endl
         << "regions, the depth is the number of haplotypes flagged there; given the" << endl
         << "regions each assembly's QC could assess, it is the number assessed. Unlike" << endl
         << "vg depth -P, this counts haplotypes, not path names." << endl
         << endl
         << "options:" << endl
         << "  -u, --tigs FILE         GAF of the tigs, from vg unitigs -g" << endl
         << "  -c, --canonical FILE    write each feature once to FILE as BED at its path" << endl
         << "                          coordinate, in priority order, with the columns of -a" << endl
         << "  -P, --priority PREFIX   name features off the gref on paths whose name starts" << endl
         << "                          with PREFIX first (may repeat, first given first);" << endl
         << "                          then come reference paths, then the rest by name" << endl
         << "  -a, --assemblies FILE   write each place where a path visits a feature to FILE" << endl
         << "                          as BED on the full path: the feature's name, score" << endl
         << "                          and strand (., 0 and . if it has none), any further" << endl
         << "                          BED columns, and then its ID" << endl
         << "  -m, --min-match F       write a place wherever a path visits at least F of the" << endl
         << "                          feature's bases in order (liftover), and end each row" << endl
         << "                          with that fraction. Without -m, a place must hold the" << endl
         << "                          feature's whole walk. A gap on the path longer than" << endl
         << "                          the feature splits a place." << endl
         << "  -b, --blocks            write places as BED12, with a block for each stretch" << endl
         << "                          that holds feature bases" << endl
         << "      --gref-segs FILE    with a VCF, take the top-level snarl of each alt" << endl
         << "                          contig for INFO/PS from FILE, the table of gref" << endl
         << "                          fragments from vg paths --gref-segs" << endl
         << "  -Q, --paths-by PREFIX   only write places (with --depth, rows) on paths whose" << endl
         << "                          name starts with PREFIX (may repeat). With --depth," << endl
         << "                          the default is the gref paths." << endl
         << "      --depth             write the depth of haplotypes along paths as bedGraph" << endl
         << "                          (see above); takes no -u, -c, -P, -a, -m or -b" << endl
         << "      --visits            with --depth, count the visits of haplotype paths" << endl
         << "                          instead of distinct haplotypes" << endl
         << "      --references        with --depth, count reference paths as haplotypes" << endl
         << "  -t, --threads N         number of threads to use [all available]" << endl
         << "  -p, --progress          show progress" << endl
         << "  -h, --help              print this help message to stderr and exit" << endl;
}

/// The options and arguments of vg lift.
struct LiftConfig {
    string graph_file;
    /// The BED or VCF of features, or with --depth the regions BED, if any.
    string input_file;
    /// True if the features are the ALT alleles of a VCF.
    bool vcf_input = false;
    string tigs_file;
    string canonical_file;
    vector<string> priority_prefixes;
    string assemblies_file;
    vector<string> path_prefixes;
    double min_match = 0;
    bool blocks = false;
    bool progress = false;
    bool depth_mode = false;
    HaplotypeDepth::Options depth_options;
    string gref_segs_file;

    /**
     * Read the command line.
     *
     * Print the help and exit on -h or an unknown option. Report an error on
     * options that do not go together, and on missing or extra arguments.
     */
    LiftConfig(int argc, char** argv, const Logger& logger);

private:
    /// Read the options, leaving optind at the first argument.
    void read_options(int argc, char** argv, const Logger& logger);
    /// Return true if an option for lifting features is set.
    bool has_lift_options() const;
};

LiftConfig::LiftConfig(int argc, char** argv, const Logger& logger) {
    read_options(argc, argv, logger);
    if (!depth_mode && (depth_options.count_visits || depth_options.count_references)) {
        logger.error() << "--visits and --references need --depth" << endl;
    }
    if (depth_mode) {
        if (has_lift_options()) {
            logger.error() << "--depth takes no -u, -c, -P, -a, -m, -b or --gref-segs" << endl;
        }
        if (optind + 1 != argc && optind + 2 != argc) {
            logger.error() << "expected a GBZ file and optionally a regions BED" << endl;
        }
        graph_file = require_exists(logger, argv[optind]);
        input_file = optind + 2 == argc ? argv[optind + 1] : "";
        return;
    }
    if (optind + 2 != argc) {
        logger.error() << "expected a GBZ file and a BED or VCF file" << endl;
    }
    graph_file = require_exists(logger, argv[optind]);
    input_file = argv[optind + 1];
    vcf_input = input_file != "-" && is_variant_file(input_file);
    if (!vcf_input && !gref_segs_file.empty()) {
        logger.error() << "--gref-segs needs a VCF" << endl;
    }
}

void LiftConfig::read_options(int argc, char** argv, const Logger& logger) {
    constexpr int OPT_DEPTH = 1001;
    constexpr int OPT_VISITS = 1002;
    constexpr int OPT_REFERENCES = 1003;
    constexpr int OPT_GREF_SEGS = 1004;

    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
            {"depth", no_argument, nullptr, OPT_DEPTH},
            {"visits", no_argument, nullptr, OPT_VISITS},
            {"references", no_argument, nullptr, OPT_REFERENCES},
            {"gref-segs", required_argument, nullptr, OPT_GREF_SEGS},
            {"tigs", required_argument, nullptr, 'u'},
            {"canonical", required_argument, nullptr, 'c'},
            {"priority", required_argument, nullptr, 'P'},
            {"assemblies", required_argument, nullptr, 'a'},
            {"min-match", required_argument, nullptr, 'm'},
            {"blocks", no_argument, nullptr, 'b'},
            {"paths-by", required_argument, nullptr, 'Q'},
            {"threads", required_argument, nullptr, 't'},
            {"progress", no_argument, nullptr, 'p'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}
        };
        int option_index = 0;
        int c = getopt_long(argc, argv, "u:c:P:a:m:bQ:t:ph?", long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
        case 'u':
            tigs_file = require_exists(logger, optarg);
            break;
        case 'c':
            canonical_file = ensure_writable(logger, optarg);
            break;
        case 'P':
            priority_prefixes.emplace_back(optarg);
            break;
        case 'a':
            assemblies_file = ensure_writable(logger, optarg);
            break;
        case 'm':
            min_match = parse<double>(optarg);
            if (min_match <= 0 || min_match > 1) {
                logger.error() << "-m/--min-match must be in (0, 1]" << endl;
            }
            break;
        case 'b':
            blocks = true;
            break;
        case 'Q':
            path_prefixes.emplace_back(optarg);
            break;
        case 't':
            set_thread_count(logger, optarg);
            break;
        case 'p':
            progress = true;
            break;
        case OPT_DEPTH:
            depth_mode = true;
            break;
        case OPT_VISITS:
            depth_options.count_visits = true;
            break;
        case OPT_REFERENCES:
            depth_options.count_references = true;
            break;
        case OPT_GREF_SEGS:
            gref_segs_file = require_exists(logger, optarg);
            break;
        case 'h':
        case '?':
        default:
            help_lift(argv);
            std::exit(EXIT_FAILURE);
        }
    }
}

bool LiftConfig::has_lift_options() const {
    return !tigs_file.empty() || !canonical_file.empty() || !priority_prefixes.empty() || !assemblies_file.empty() ||
           min_match > 0 || blocks || !gref_segs_file.empty();
}

//------------------------------------------------------------------------------

/// Return true if the name starts with one of the prefixes.
static bool starts_with_any(string_view name, const vector<string>& prefixes) {
    return std::any_of(prefixes.begin(), prefixes.end(), [&](const string& prefix) {
        return name.starts_with(prefix);
    });
}

/// Return the number of features without a path coordinate.
static size_t count_unnamed(const vector<GraphFeature>& features) {
    return std::count_if(features.begin(), features.end(), [](const GraphFeature& feature) {
        return !feature.coordinate.has_value();
    });
}

/// Return the ID of the feature with the given index.
static string feature_id(size_t i) {
    return "f" + to_string(i + 1);
}

/**
 * Return the contents of a text file, plain or compressed, as a stream.
 *
 * @throws std::runtime_error if the file cannot be read.
 */
static istringstream read_text_file(const string& path) {
    htsFile* file = hts_open(path.c_str(), "r");
    if (file == nullptr) {
        throw runtime_error("cannot open " + path);
    }
    string text;
    kstring_t line = KS_INITIALIZE;
    int status;
    while ((status = hts_getline(file, '\n', &line)) >= 0) {
        text.append(line.s, line.l).push_back('\n');
    }
    ks_free(&line);
    hts_close(file);
    if (status < -1) {
        throw runtime_error("cannot read " + path);
    }
    return istringstream(text);
}

//------------------------------------------------------------------------------

/// Read the features of a BED file, and warn that BED12 blocks are not lifted.
static vector<BedFeature> read_bed_file(const string& features_file, const Logger& logger, bool progress) {
    vector<BedFeature> bed_features;
    get_input_file(features_file, [&](istream& in) {
        try {
            bed_features = read_bed_features(in);
        } catch (const runtime_error& e) {
            logger.error() << e.what() << endl;
        }
    });
    if (progress) {
        cerr << "[vg lift] " << bed_features.size() << " features" << endl;
    }
    if (std::any_of(bed_features.begin(), bed_features.end(), [](const BedFeature& feature) {
            return feature.fields.size() >= 9;
        })) {
        logger.warn() << "BED12 blocks are not lifted; columns from the 7th on are copied as they are" << endl;
    }
    return bed_features;
}

/// Report an error for a feature that ends past the end of its sequence, or
/// that is on a sequence that is not in the cover or the path directory.
static void check_sequences(const vector<BedFeature>& bed_features, const Cover& cover, const PathDirectory& paths,
                            const LiftConfig& config, const Logger& logger) {
    for (const BedFeature& feature : bed_features) {
        optional<size_t> length = cover.sequence_length(feature.sequence);
        if (length.has_value() && feature.end > *length) {
            logger.error() << "feature " << feature.sequence << ":" << feature.start << "-" << feature.end
                           << " ends past the end of the sequence (" << *length << " bp)" << endl;
        }
        if (!length.has_value() && !paths.has_path(feature.sequence)) {
            logger.error() << "feature on " << feature.sequence << ", which is not a gref path or a path in "
                           << config.graph_file << ", or a tig"
                           << (config.tigs_file.empty() ? " (use -u)" : " in " + config.tigs_file) << endl;
        }
    }
}

/**
 * Lift each feature to the graph, from its sequence in the cover or its path
 * in the directory, in order.
 *
 * A feature on a path lifts only if one fragment of the path holds it. Leave
 * out the features that no fragment holds, and count them in skipped.
 */
static vector<GraphFeature> lift_to_graph(vector<BedFeature> bed_features, const Cover& cover,
                                          const PathDirectory& paths, const gbwtgraph::GBWTGraph& graph,
                                          const GrefIndex* gref, size_t& skipped) {
    vector<GraphFeature> features(bed_features.size());
    vector<bool> lifted(bed_features.size(), true);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < bed_features.size(); i++) {
        BedFeature& feature = bed_features[i];
        optional<WalkInterval> interval;
        if (cover.sequence_length(feature.sequence).has_value()) {
            interval = cover.interval(feature.sequence, feature.start, feature.end);
        } else {
            interval = paths.interval(feature.sequence, feature.start, feature.end);
        }
        if (interval.has_value()) {
            features[i] = make_graph_feature(graph, gref, std::move(*interval), std::move(feature.fields));
        } else {
            lifted[i] = false;
        }
    }
    skipped = 0;
    for (size_t i = 0; i < features.size(); i++) {
        if (!lifted[i]) {
            skipped++;
        } else if (skipped > 0) {
            features[i - skipped] = std::move(features[i]);
        }
    }
    features.resize(features.size() - skipped);
    return features;
}

/**
 * Read the BED features and lift each to the graph, keeping each feature once.
 *
 * A feature is on a gref path, on a tig in the tigs file, or on a path in the
 * GBZ. Features on a path outside its fragments in the graph are skipped.
 */
static vector<GraphFeature> lift_bed_features(const LiftConfig& config, const gbwtgraph::GBZ& gbz,
                                              const GrefIndex* gref_pointer, const Logger& logger) {
    Cover cover(gbz.graph, gref_pointer);
    vector<BedFeature> bed_features = read_bed_file(config.input_file, logger, config.progress);

    // Find the sequence of each feature: a gref path, a tig, or a path in the
    // GBZ. Load only the tigs and paths that features are on.
    unordered_set<string> wanted;
    for (const BedFeature& feature : bed_features) {
        if (!cover.is_gref_path(feature.sequence)) {
            wanted.insert(feature.sequence);
        }
    }
    if (!wanted.empty() && !config.tigs_file.empty()) {
        get_input_file(config.tigs_file, [&](istream& in) {
            try {
                cover.add_tigs(in, [&](const string& name) { return wanted.contains(name); });
            } catch (const runtime_error& e) {
                logger.error() << e.what() << endl;
            }
        });
    }
    PathDirectory paths(gbz, [&](const string& name) {
        return wanted.contains(name) && !cover.sequence_length(name).has_value();
    });
    check_sequences(bed_features, cover, paths, config, logger);

    size_t skipped = 0;
    vector<GraphFeature> features = lift_to_graph(std::move(bed_features), cover, paths, gbz.graph, gref_pointer,
                                                  skipped);
    if (skipped > 0) {
        logger.warn() << skipped << " features on paths are not within one fragment of the path in "
                      << config.graph_file << " and were skipped" << endl;
    }
    sort_and_merge_features(features);
    if (config.progress) {
        cerr << "[vg lift] " << features.size() << " distinct features in the graph" << endl;
    }
    return features;
}

//------------------------------------------------------------------------------

/// The feature index of an ALT allele without a feature. The VCF gives no walk
/// for it, or its walk has no nodes besides the first and last.
constexpr size_t NO_FEATURE = numeric_limits<size_t>::max();

/// Return the INFO/GR and INFO/PC values of a feature, "." where it has none.
static pair<string, string> lifted_values(const GraphFeature& feature) {
    pair<string, string> values(".", ".");
    if (!feature.placement.text.empty()) {
        // Commas separate the values of the alleles.
        values.first = feature.placement.text;
        std::replace(values.first.begin(), values.first.end(), ',', '|');
    }
    if (feature.coordinate) {
        values.second = feature.coordinate->text();
    }
    return values;
}

/**
 * A VCF of graph variants whose ALT alleles are features to lift.
 *
 * It adds a feature for each ALT allele with bases to place. It also adds one
 * for each variant that the alleles carry along gref alt contigs, unless the
 * VCF already has records on the contigs. It writes the VCF back with the gref
 * placement and path coordinate of each allele, followed by the records of
 * those nested variants.
 */
class LiftedVcf {
public:
    /// Read the VCF, and add the features of its ALT alleles and of their
    /// nested variants to the features.
    LiftedVcf(const string& vcf_file, const gbwtgraph::GBWTGraph& graph, const GrefIndex& gref,
              vector<GraphFeature>& features, const Logger& logger, bool progress);

    /**
     * Write the VCF to standard output with the tags of its ALT alleles,
     * followed by the records of the nested variants.
     *
     * @param features The features, named and placed.
     * @param top_level_snarls The top-level snarl of each alt contig, by its
     * name in the VCF, for the PS of the nested records.
     */
    void write(const vector<GraphFeature>& features, const unordered_map<string, string>& top_level_snarls,
               const Logger& logger) const;

private:
    /**
     * Add a feature for each ALT allele whose walk has nodes besides the
     * first and last, in file order.
     *
     * Each feature is stranded and reads along its allele, so that its
     * placement and path coordinate say which way the allele reads its paths.
     *
     * @returns For each ALT allele, the index of its feature, or NO_FEATURE.
     */
    vector<size_t> add_allele_features(vector<GraphFeature>& features, const GrefIndex* gref, const Logger& logger,
                                       bool progress);
    /// Return true if a record is on a gref alt contig. Such a VCF already
    /// has its nested variants, as vg deconstruct or vg lift wrote them.
    bool has_contig_records(const ContigDecomposer& decomposer, const GrefIndex& gref) const;
    /**
     * Add a feature for the ALT allele of each nested variant, its walk
     * without the first and last nodes, if that has any nodes.
     *
     * @returns For each variant, the index of its feature, or NO_FEATURE.
     */
    vector<size_t> add_nested_features(vector<GraphFeature>& features, const GrefIndex* gref) const;
    /// Return INFO/GR, the gref placement, and INFO/PC, the path coordinate,
    /// for each ALT allele.
    vector<VcfAlleleInfo> lifted_tags(const vector<GraphFeature>& features) const;
    /// Return the INFO/GR and INFO/PC entries of the record of each nested
    /// variant, each starting with ";", or empty if the variant has neither.
    vector<string> nested_tags(const vector<GraphFeature>& features) const;

    string vcf_file;
    const gbwtgraph::GBWTGraph& graph;
    VcfCalls calls;
    /// The variants that the ALT alleles carry along gref alt contigs.
    optional<NestedVariantCatalog> catalog;
    /// The index of the feature of each ALT allele, or NO_FEATURE.
    vector<size_t> allele_features;
    /// The index of the feature of each nested variant, or NO_FEATURE.
    vector<size_t> nested_features;
};

LiftedVcf::LiftedVcf(const string& vcf_file, const gbwtgraph::GBWTGraph& graph, const GrefIndex& gref,
                     vector<GraphFeature>& features, const Logger& logger, bool progress) :
    vcf_file(vcf_file), graph(graph) {
    try {
        calls = read_vcf_calls(vcf_file, graph, vg::get_thread_count());
    } catch (const runtime_error& e) {
        logger.error() << e.what() << endl;
    }
    const GrefIndex* gref_pointer = gref.path_count() > 0 ? &gref : nullptr;
    allele_features = add_allele_features(features, gref_pointer, logger, progress);
    ContigDecomposer decomposer(graph, gref);
    if (has_contig_records(decomposer, gref)) {
        logger.warn() << vcf_file << " already has records on gref alt contigs, so no nested variants "
                      << "are added" << endl;
        return;
    }
    if (gref_pointer == nullptr) {
        return;
    }
    catalog.emplace(calls, decomposer, gref, graph);
    nested_features = add_nested_features(features, gref_pointer);
    if (progress) {
        const auto& contigs = catalog->contigs();
        size_t contigs_with_variants = count_if(contigs.begin(), contigs.end(),
                                                [](const NestedVariantCatalog::Contig& contig) {
                                                    return !contig.variants.empty();
                                                });
        cerr << "[vg lift] " << catalog->variants().size() << " nested variants on " << contigs_with_variants
             << " gref alt contigs, of " << contigs.size() << " that ALT alleles read" << endl;
    }
}

vector<size_t> LiftedVcf::add_allele_features(vector<GraphFeature>& features, const GrefIndex* gref,
                                              const Logger& logger, bool progress) {
    vector<VcfAllele>& alleles = calls.alleles;
    vector result(alleles.size(), NO_FEATURE);
    size_t first_feature = features.size();
    size_t without_walk = 0;
    size_t foreign = 0;
    string foreign_example;
    for (size_t i = 0; i < alleles.size(); i++) {
        if (alleles[i].walk.empty()) {
            without_walk++;
        }
        if (!alleles[i].holds_alt) {
            if (foreign == 0) {
                foreign_example = alleles[i].name;
            }
            foreign++;
        }
        if (alleles[i].interval) {
            result[i] = features.size();
            features.emplace_back();
        }
    }
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < alleles.size(); i++) {
        if (result[i] != NO_FEATURE) {
            features[result[i]] = make_graph_feature(graph, gref, std::move(*alleles[i].interval),
                                                     { alleles[i].name, "0", "+" },
                                                     alleles[i].walk.front(), alleles[i].walk.back());
        }
    }
    if (progress) {
        cerr << "[vg lift] " << alleles.size() << " ALT alleles, " << (features.size() - first_feature)
             << " with bases to place" << endl;
    }
    if (without_walk > 0) {
        logger.warn() << without_walk << " ALT alleles have no walk in INFO/AT or in a PanGenie variant ID" << endl;
    }
    if (foreign > 0) {
        logger.warn() << foreign << " ALT alleles, such as " << foreign_example
                      << ", are not in the sequences of their walks; is the VCF from this graph?" << endl;
    }
    return result;
}

bool LiftedVcf::has_contig_records(const ContigDecomposer& decomposer, const GrefIndex& gref) const {
    unordered_set<string> contig_names;
    for (uint32_t path = 0; path < gref.path_count(); path++) {
        if (decomposer.is_contig(path)) {
            contig_names.insert(PathMetadata::parse_locus_name(gref.path_name(path)));
        }
    }
    return std::any_of(calls.sites.begin(), calls.sites.end(), [&](const VcfSite& site) {
        return contig_names.contains(site.chrom);
    });
}

vector<size_t> LiftedVcf::add_nested_features(vector<GraphFeature>& features, const GrefIndex* gref) const {
    const vector<NestedVariant>& variants = catalog->variants();
    vector result(variants.size(), NO_FEATURE);
    for (size_t i = 0; i < variants.size(); i++) {
        if (variants[i].alt_walk.size() > 2) {
            result[i] = features.size();
            features.emplace_back();
        }
    }
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < variants.size(); i++) {
        if (result[i] != NO_FEATURE) {
            const vector<handle_t>& walk = variants[i].alt_walk;
            WalkInterval interval;
            interval.walk.assign(walk.begin() + 1, walk.end() - 1);
            interval.length = interval.walk_length(graph);
            features[result[i]] = make_graph_feature(graph, gref, std::move(interval), { catalog->id(i), "0", "+" },
                                                     walk.front(), walk.back());
        }
    }
    return result;
}

vector<VcfAlleleInfo> LiftedVcf::lifted_tags(const vector<GraphFeature>& features) const {
    VcfAlleleInfo placements { "GR",
        "##INFO=<ID=GR,Number=A,Type=String,Description=\"Gref placement of the ALT allele, its walk without the "
        "first and last nodes: PATH:START-END:STRAND:CIGAR segments, 1-based and inclusive, separated by |, or *:NI "
        "if none of its N bases is on the gref (vg lift)\">",
        vector<string>(allele_features.size(), ".") };
    VcfAlleleInfo coordinates { "PC",
        "##INFO=<ID=PC,Number=A,Type=String,Description=\"Path coordinate of the ALT allele, its walk without the "
        "first and last nodes: PATH:START-END:STRAND, 1-based and inclusive, on the first path in priority order "
        "that carries it whole (vg lift)\">",
        vector<string>(allele_features.size(), ".") };
    for (size_t i = 0; i < allele_features.size(); i++) {
        if (allele_features[i] != NO_FEATURE) {
            tie(placements.values[i], coordinates.values[i]) = lifted_values(features[allele_features[i]]);
        }
    }
    return { placements, coordinates };
}

vector<string> LiftedVcf::nested_tags(const vector<GraphFeature>& features) const {
    vector<string> tags(catalog->variants().size());
    for (size_t i = 0; i < tags.size(); i++) {
        if (nested_features[i] != NO_FEATURE) {
            auto [placement, coordinate] = lifted_values(features[nested_features[i]]);
            tags[i] = (placement == "." ? "" : ";GR=" + placement) + (coordinate == "." ? "" : ";PC=" + coordinate);
        }
    }
    return tags;
}

void LiftedVcf::write(const vector<GraphFeature>& features, const unordered_map<string, string>& top_level_snarls,
                      const Logger& logger) const {
    // The genotyper writes the nested records to a temporary file as it goes,
    // and they follow the calls.
    string nested_file = temp_file::create("vg-lift-nested");
    ofstream nested_out(nested_file);
    optional<NestedGenotyper> genotyper;
    vector<string> nested_header;
    if (catalog.has_value() && !catalog->variants().empty()) {
        genotyper.emplace(*catalog, calls, nested_tags(features), top_level_snarls, graph, nested_out);
        nested_header = genotyper->header_lines();
    }
    try {
        write_vcf_with_allele_info(vcf_file, cout, lifted_tags(features), nested_header,
                                   [&](size_t record, const VcfCursor& cursor) {
            if (genotyper) {
                genotyper->add_call(record, cursor);
            }
        }, vg::get_thread_count());
    } catch (const runtime_error& e) {
        logger.error() << e.what() << endl;
    }
    nested_out.close();
    ifstream nested_in(nested_file);
    if (nested_in.peek() != ifstream::traits_type::eof()) {
        cout << nested_in.rdbuf();
    }
    nested_in.close();
    temp_file::remove(nested_file);
    if (genotyper && genotyper->conflicts() > 0) {
        logger.warn() << genotyper->conflicts() << " nested genotypes are missing because the sample carries more "
                      << "walks there than it has haplotypes" << endl;
    }
}

/**
 * Return the top-level snarl of each alt contig, by its name in the VCF, from
 * a table of gref fragments, or none without a table.
 */
static unordered_map<string, string> read_contig_snarls(const string& gref_segs_file, const Logger& logger) {
    unordered_map<string, string> top_level_snarls;
    if (gref_segs_file.empty()) {
        return top_level_snarls;
    }
    try {
        istringstream segs = read_text_file(gref_segs_file);
        for (auto& [contig, snarl] : GrefCover::read_top_level_snarls(segs)) {
            top_level_snarls.try_emplace(PathMetadata::parse_locus_name(contig), std::move(snarl));
        }
    } catch (const runtime_error& e) {
        logger.error() << gref_segs_file << ": " << e.what() << endl;
    }
    return top_level_snarls;
}

//------------------------------------------------------------------------------

/**
 * Write BED columns 4 to 6 (name, score and strand) for a place of the
 * feature, each after a tab. The place reads the feature in reverse if
 * is_reverse is set.
 */
static void write_name_score_strand(ostream& out, const GraphFeature& feature, bool is_reverse) {
    const vector<string>& fields = feature.fields;
    out << '\t' << (!fields.empty() ? fields[0] : ".") << '\t' << (fields.size() > 1 ? fields[1] : "0") << '\t';
    if (feature.is_stranded()) {
        out << (is_reverse ? '-' : '+');
    } else {
        out << '.';
    }
}

/// Write the BED columns of the feature from the 7th on, then the ID of the
/// feature, which has index i, each after a tab.
static void write_rest(ostream& out, const GraphFeature& feature, size_t i) {
    const vector<string>& fields = feature.fields;
    for (size_t j = BED_STRAND_FIELD + 1; j < fields.size(); j++) {
        out << '\t' << fields[j];
    }
    out << '\t' << feature_id(i);
}

/// Write BED columns 7 to 12 for a place, each after a tab. The thick part is
/// the whole place, the color is 0, and there is a block for each stretch that
/// holds feature bases.
static void write_blocks(ostream& out, const PathOccurrence& place) {
    out << '\t' << place.start << '\t' << place.end << "\t0\t" << place.blocks.size() << '\t';
    for (const auto& [block_start, block_end] : place.blocks) {
        out << (block_end - block_start) << ',';
    }
    out << '\t';
    for (const auto& [block_start, block_end] : place.blocks) {
        out << (block_start - place.start) << ',';
    }
}

/**
 * Write the BED row of a place on the named path, as -a writes it.
 *
 * With blocks, the row is BED12. With a min_match, it ends with the fraction
 * of the feature's bases that the place holds.
 */
static void write_place(ostream& out, const string& path, const PathOccurrence& place, const GraphFeature& feature,
                        bool blocks, double min_match) {
    out << path << '\t' << place.start << '\t' << place.end;
    write_name_score_strand(out, feature, place.is_reverse);
    if (blocks) {
        write_blocks(out, place);
    }
    write_rest(out, feature, place.feature);
    if (min_match > 0) {
        out << '\t' << static_cast<double>(place.matched) / static_cast<double>(feature.interval.length);
    }
    out << '\n';
}

/// Name the feature on the path at a place, if the place holds the whole
/// feature in one block and the path comes before the feature's naming path,
/// if any, in the priority order.
static void name_by_place(GraphFeature& feature, const string& path, size_t rank, const PathOccurrence& place) {
    // A place that holds the whole feature in one block is a path that carries it.
    if (size_t length = feature.interval.length;
        place.matched != length || place.blocks.size() != 1 || place.end - place.start != length) {
        return;
    }
    PathCoordinate here { path, rank, place.start, place.end, place.is_reverse };
    if (!feature.coordinate || here < *feature.coordinate) {
        feature.coordinate = std::move(here);
    }
}

/**
 * Name the features off the gref by the paths that carry them, and write the
 * places of the features on the paths chosen with -Q for -a.
 *
 * One pass over the paths other than gref paths does both. Naming uses every
 * path, whatever -Q says.
 */
static void name_and_place_features(const gbwtgraph::GBZ& gbz, vector<GraphFeature>& features,
                                    const LiftConfig& config, size_t gref_path_count, const Logger& logger) {
    PathPriority priority(gbz, config.priority_prefixes, GrefCover::gref_prefix, gref_path_count);
    ofstream out;
    if (!config.assemblies_file.empty()) {
        out.open(config.assemblies_file);
        if (!out) {
            logger.error() << "could not open " << config.assemblies_file << " for writing" << endl;
        }
    }
    size_t written = 0;
    auto visit = [&](const string& name, const vector<PathOccurrence>& places) {
        bool write = out.is_open() && (config.path_prefixes.empty() || starts_with_any(name, config.path_prefixes));
        size_t rank = priority.rank(name);
        for (const PathOccurrence& place : places) {
            name_by_place(features[place.feature], name, rank, place);
            if (write) {
                write_place(out, name, place, features[place.feature], config.blocks, config.min_match);
                written++;
            }
        }
    };
    auto not_gref = [&](string_view name) {
        return !name.starts_with(GrefCover::gref_prefix);
    };
    if (config.min_match > 0) {
        for_each_path_liftover(gbz, features, config.min_match, not_gref, visit);
    } else {
        for_each_path_occurrences(gbz, features, not_gref, visit);
    }
    if (config.progress && out.is_open()) {
        cerr << "[vg lift] " << written << " places on paths" << endl;
    }
}

/// Write a GAF record for each feature to standard output, with its tags.
static void write_feature_gaf(const vector<GraphFeature>& features, const gbwtgraph::GBWTGraph& graph) {
    for (size_t i = 0; i < features.size(); i++) {
        const GraphFeature& feature = features[i];
        write_gaf_columns(cout, graph, feature_id(i), feature.interval);
        cout << "\tsp:i:" << feature.support;
        if (!feature.placement.text.empty()) {
            cout << "\tgr:Z:" << feature.placement.text;
        }
        if (feature.coordinate) {
            cout << "\tpc:Z:" << feature.coordinate->text();
        }
        if (!feature.fields.empty()) {
            cout << "\tan:Z:" << feature.fields.front();
        }
        cout << '\n';
    }
}

/// Write each feature with a path coordinate once, as a BED row at that
/// coordinate, in priority order.
static void write_canonical(const vector<GraphFeature>& features, const string& canonical_file,
                            const Logger& logger) {
    ofstream out(canonical_file);
    if (!out) {
        logger.error() << "could not open " << canonical_file << " for writing" << endl;
    }
    vector<size_t> order;
    for (size_t i = 0; i < features.size(); i++) {
        if (features[i].coordinate) {
            order.push_back(i);
        }
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return *features[a].coordinate < *features[b].coordinate;
    });
    for (size_t i : order) {
        const PathCoordinate& place = *features[i].coordinate;
        out << place.path << '\t' << place.start << '\t' << place.end;
        write_name_score_strand(out, features[i], place.is_reverse);
        write_rest(out, features[i], i);
        out << '\n';
    }
}

//------------------------------------------------------------------------------

/// Load the regions in a BED file, or none if the file name is empty.
static unique_ptr<PathRegions> read_regions(const string& regions_file, const Logger& logger) {
    unique_ptr<PathRegions> regions;
    if (regions_file.empty()) {
        return regions;
    }
    get_input_file(regions_file, [&](istream& in) {
        try {
            regions = make_unique<PathRegions>(in);
        } catch (const runtime_error& e) {
            logger.error() << e.what() << endl;
        }
    });
    return regions;
}

/// A path to write depth along.
struct OutputPath {
    /// The name of the full path, as for_each_path_start() gives it.
    string name;
    /// Where the path starts on the full path.
    uint64_t start;
    gbwt::size_type id;
};

/**
 * Return the paths to write depth along, sorted by name and then along each
 * full path. They are the gref paths, or those whose names start with one of
 * the prefixes.
 */
static vector<OutputPath> output_paths(const gbwtgraph::GBZ& gbz, const vector<string>& path_prefixes) {
    vector<OutputPath> outputs;
    for_each_path_start(gbz.index, [&](gbwt::size_type path_id, const string& name, uint64_t start) {
        if (path_prefixes.empty() ? GrefCover::is_gref_derived(name) : starts_with_any(name, path_prefixes)) {
            outputs.emplace_back(name, start, path_id);
        }
    });
    std::sort(outputs.begin(), outputs.end(), [](const OutputPath& a, const OutputPath& b) {
        return std::tie(a.name, a.start) < std::tie(b.name, b.start);
    });
    return outputs;
}

/**
 * Write the depth of haplotypes along paths as bedGraph on standard output.
 *
 * The paths are the gref paths, or those whose names start with one of the
 * prefixes, in name order and then in order along each full path.
 *
 * @param regions_file A BED of the regions to count in, or empty to count
 * every base.
 */
static void write_depth(const gbwtgraph::GBZ& gbz, const string& graph_file, const string& regions_file,
                        HaplotypeDepth::Options options, const vector<string>& path_prefixes, const Logger& logger,
                        bool progress) {
    unique_ptr<PathRegions> regions = read_regions(regions_file, logger);
    options.regions = regions.get();
    HaplotypeDepth depth(gbz, options);
    if (progress) {
        cerr << "[vg lift] counted " << depth.haplotype_count() << " haplotypes";
        if (regions) {
            cerr << ", " << depth.haplotypes_without_regions() << " of them without regions";
        }
        cerr << endl;
    }
    if (const vector<string>& unmatched = depth.uncounted_region_paths(); !unmatched.empty()) {
        logger.warn() << unmatched.size() << " paths in " << regions_file << ", such as " << unmatched.front()
                      << ", are not counted paths in " << graph_file << ", so their regions are ignored" << endl;
    }

    vector<OutputPath> outputs = output_paths(gbz, path_prefixes);
    if (outputs.empty()) {
        logger.warn() << "no " << (path_prefixes.empty() ? "gref paths" : "paths chosen with -Q") << " in "
                      << graph_file << ", so no depth is written" << endl;
    }
    size_t rows = 0;
    const gbwtgraph::GBWTGraph& graph = gbz.graph;
    // The nodes that the output paths visit, by node ID minus the smallest.
    nid_t min_id = graph.get_node_count() > 0 ? graph.min_node_id() : 0;
    vector<bool> on_output(graph.get_node_count() > 0 ? graph.max_node_id() - min_id + 1 : 0, false);
    for (const OutputPath& path : outputs) {
        PathSteps steps(graph, path.id);
        for (const PathSteps::Step& step : steps.steps()) {
            on_output[graph.get_id(step.handle) - min_id] = true;
        }
        depth.for_each_run(steps, path.start, [&](uint64_t start, uint64_t end, uint32_t value) {
            cout << path.name << '\t' << start << '\t' << end << '\t' << value << '\n';
            rows++;
        });
    }
    cout.flush();
    if (progress) {
        cerr << "[vg lift] " << rows << " depth rows on " << outputs.size() << " paths" << endl;
    }

    // Say how much of what the haplotypes read the rows leave out, such as
    // haplotype sequence in gref fragments too short to keep.
    uint64_t read = 0;
    uint64_t left_out = 0;
    graph.for_each_handle([&](const handle_t& handle) {
        nid_t id = graph.get_id(handle);
        uint64_t bases = depth.bases_read(id, static_cast<uint32_t>(graph.get_length(handle)));
        read += bases;
        if (!on_output[id - min_id]) {
            left_out += bases;
        }
    });
    if (left_out > 0) {
        logger.info() << left_out << " of the " << read << " " << (options.count_visits ? "visited" : "haplotype")
                      << " bases counted are on nodes that no output path visits, so the rows leave them out"
                      << endl;
    }
}

//------------------------------------------------------------------------------

int main_lift(int argc, char** argv) {
    Logger logger("vg lift");
    if (argc == 2) {
        help_lift(argv);
        return 1;
    }
    LiftConfig config(argc, argv, logger);
    if (config.depth_mode) {
        gbwtgraph::GBZ gbz;
        load_gbz(gbz, config.graph_file, config.progress);
        write_depth(gbz, config.graph_file, config.input_file, config.depth_options, config.path_prefixes, logger,
                    config.progress);
        return 0;
    }
    unordered_map<string, string> top_level_snarls = read_contig_snarls(config.gref_segs_file, logger);

    gbwtgraph::GBZ gbz;
    load_gbz(gbz, config.graph_file, config.progress);
    const gbwtgraph::GBWTGraph& graph = gbz.graph;
    if (!config.tigs_file.empty()) {
        // Tigs from another graph would name the wrong walks.
        require_compatible_reference(config.tigs_file, nullptr, &gbz, false);
    }
    GrefIndex gref(gbz, GrefCover::gref_prefix);
    const GrefIndex* gref_pointer = gref.path_count() > 0 ? &gref : nullptr;
    if (config.progress) {
        cerr << "[vg lift] gref: " << gref.path_count() << " paths, " << gref.total_length() << " bp" << endl;
    }

    vector<GraphFeature> features;
    optional<LiftedVcf> lifted_vcf;
    if (config.vcf_input) {
        lifted_vcf.emplace(config.input_file, graph, gref, features, logger, config.progress);
    } else {
        features = lift_bed_features(config, gbz, gref_pointer, logger);
    }
    if (count_unnamed(features) > 0 || !config.assemblies_file.empty()) {
        name_and_place_features(gbz, features, config, gref.path_count(), logger);
    }
    if (lifted_vcf.has_value()) {
        lifted_vcf->write(features, top_level_snarls, logger);
    } else {
        write_feature_gaf(features, graph);
    }
    cout.flush();

    if (!config.canonical_file.empty()) {
        write_canonical(features, config.canonical_file, logger);
    }
    if (size_t nameless = count_unnamed(features); nameless > 0) {
        logger.warn() << nameless << " features are not carried whole by any path, so they have no path coordinate"
                      << endl;
    }
    return 0;
}

// Register subcommand
static Subcommand vg_lift("lift", "lift features on a cover or on paths to the graph, gref and assemblies", WIDGET, main_lift);
