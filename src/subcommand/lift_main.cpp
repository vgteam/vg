/**
 * \file lift_main.cpp: lift features on a cover or on paths to the graph, the gref and the paths
 */

#include <getopt.h>
#include <omp.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_set>

#include "subcommand.hpp"
#include "../feature_lift.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gref.hpp"
#include "../nested_variants.hpp"
#include "../utility.hpp"
#include "../vcf_alleles.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_lift(char** argv) {
    cerr << "usage: " << argv[0] << " lift [options] graph.gbz features.bed > features.gaf" << endl
         << "       " << argv[0] << " lift [options] graph.gbz variants.vcf > variants.lifted.vcf" << endl
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
         << "are among the carriers of the contig. A VCF that already has records on alt" << endl
         << "contigs gets none added." << endl
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
         << "  -Q, --paths-by PREFIX   only write places on paths whose name starts with" << endl
         << "                          PREFIX (may repeat)" << endl
         << "  -t, --threads N         number of threads to use [all available]" << endl
         << "  -p, --progress          show progress" << endl
         << "  -h, --help              print this help message to stderr and exit" << endl;
}

/**
 * Read the BED features and lift each to the graph, keeping each feature once.
 *
 * A feature is on a gref path, on a tig in the tigs file, or on a path in the
 * GBZ. Features on a path outside its fragments in the graph are skipped.
 */
static vector<GraphFeature> lift_bed_features(const string& features_file, const string& tigs_file,
                                              const string& graph_file, const gbwtgraph::GBZ& gbz,
                                              const GrefIndex* gref_pointer, const Logger& logger, bool progress) {
    const gbwtgraph::GBWTGraph& graph = gbz.graph;
    Cover cover(graph, gref_pointer);

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
    for (const BedFeature& feature : bed_features) {
        if (feature.fields.size() >= 9) {
            logger.warn() << "BED12 blocks are not lifted; columns from the 7th on are copied as they are" << endl;
            break;
        }
    }

    // Find the sequence of each feature: a gref path, a tig, or a path in the
    // GBZ. Load only the tigs and paths that features are on.
    unordered_set<string> wanted;
    for (const BedFeature& feature : bed_features) {
        if (!cover.is_gref_path(feature.sequence)) {
            wanted.insert(feature.sequence);
        }
    }
    if (!wanted.empty() && !tigs_file.empty()) {
        get_input_file(tigs_file, [&](istream& in) {
            try {
                cover.add_tigs(in, [&](const string& name) { return wanted.count(name) > 0; });
            } catch (const runtime_error& e) {
                logger.error() << e.what() << endl;
            }
        });
    }
    PathDirectory paths(gbz, [&](const string& name) {
        return wanted.count(name) > 0 && !cover.sequence_length(name);
    });
    for (const BedFeature& feature : bed_features) {
        optional<size_t> length = cover.sequence_length(feature.sequence);
        if (length && feature.end > *length) {
            logger.error() << "feature " << feature.sequence << ":" << feature.start << "-" << feature.end
                           << " ends past the end of the sequence (" << *length << " bp)" << endl;
        }
        if (!length && !paths.has_path(feature.sequence)) {
            logger.error() << "feature on " << feature.sequence << ", which is not a gref path or a path in "
                           << graph_file << ", or a tig" << (tigs_file.empty() ? " (use -u)" : " in " + tigs_file) << endl;
        }
    }

    // A feature on a path lifts only if one fragment of the path holds it.
    vector<GraphFeature> features(bed_features.size());
    vector<bool> lifted(bed_features.size(), true);
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < bed_features.size(); i++) {
        BedFeature& feature = bed_features[i];
        optional<WalkInterval> interval;
        if (cover.sequence_length(feature.sequence)) {
            interval = cover.interval(feature.sequence, feature.start, feature.end);
        } else {
            interval = paths.interval(feature.sequence, feature.start, feature.end);
        }
        if (interval) {
            features[i] = make_graph_feature(graph, gref_pointer, std::move(*interval), std::move(feature.fields));
        } else {
            lifted[i] = false;
        }
    }
    size_t outside = 0;
    for (size_t i = 0; i < features.size(); i++) {
        if (!lifted[i]) {
            outside++;
        } else if (outside > 0) {
            features[i - outside] = std::move(features[i]);
        }
    }
    features.resize(features.size() - outside);
    if (outside > 0) {
        logger.warn() << outside << " features on paths are not within one fragment of the path in "
                      << graph_file << " and were skipped" << endl;
    }
    bed_features = vector<BedFeature>();
    sort_and_merge_features(features);
    if (progress) {
        cerr << "[vg lift] " << features.size() << " distinct features in the graph" << endl;
    }
    return features;
}

/// The feature index of an ALT allele without a feature. The VCF gives no walk
/// for it, or its walk has no nodes besides the first and last.
constexpr size_t NO_FEATURE = numeric_limits<size_t>::max();

/**
 * Add a feature for each ALT allele of the VCF whose walk has nodes besides
 * the first and last, in file order.
 *
 * Each feature is stranded and reads along its allele, so that its placement
 * and path coordinate say which way the allele reads its paths.
 *
 * @returns For each ALT allele, the index of its feature, or NO_FEATURE.
 */
static vector<size_t> add_allele_features(VcfCalls& calls, const gbwtgraph::GBWTGraph& graph, const GrefIndex* gref,
                                          vector<GraphFeature>& features, const Logger& logger, bool progress) {
    vector<VcfAllele>& alleles = calls.alleles;
    vector<size_t> allele_features(alleles.size(), NO_FEATURE);
    size_t first_feature = features.size();
    size_t without_walk = 0;
    size_t foreign = 0;
    string foreign_example;
    for (size_t i = 0; i < alleles.size(); i++) {
        without_walk += alleles[i].walk.empty();
        if (!alleles[i].holds_alt && foreign++ == 0) {
            foreign_example = alleles[i].name;
        }
        if (alleles[i].interval) {
            allele_features[i] = features.size();
            features.emplace_back();
        }
    }
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < alleles.size(); i++) {
        if (allele_features[i] != NO_FEATURE) {
            features[allele_features[i]] = make_graph_feature(graph, gref, std::move(*alleles[i].interval),
                                                              { alleles[i].name, "0", "+" });
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
    return allele_features;
}

/**
 * Return true if a record of the calls is on a gref alt contig. Such a VCF
 * already has its nested variants, as vg deconstruct or vg lift wrote them.
 */
static bool has_contig_records(const VcfCalls& calls, const ContigDecomposer& decomposer, const GrefIndex& gref) {
    unordered_set<string> contig_names;
    for (uint32_t path = 0; path < gref.path_count(); path++) {
        if (decomposer.is_contig(path)) {
            contig_names.insert(PathMetadata::parse_locus_name(gref.path_name(path)));
        }
    }
    return std::any_of(calls.sites.begin(), calls.sites.end(), [&](const VcfSite& site) {
        return contig_names.count(site.chrom) > 0;
    });
}

/**
 * Add a feature for the ALT allele of each nested variant, its walk without
 * the first and last nodes, if that has any nodes.
 *
 * @returns For each variant, the index of its feature, or NO_FEATURE.
 */
static vector<size_t> add_nested_features(const NestedVariantCatalog& catalog, const gbwtgraph::GBWTGraph& graph,
                                          const GrefIndex* gref, vector<GraphFeature>& features) {
    const vector<NestedVariant>& variants = catalog.variants();
    vector<size_t> nested_features(variants.size(), NO_FEATURE);
    for (size_t i = 0; i < variants.size(); i++) {
        if (variants[i].alt_walk.size() > 2) {
            nested_features[i] = features.size();
            features.emplace_back();
        }
    }
    #pragma omp parallel for schedule(dynamic, 1024)
    for (size_t i = 0; i < variants.size(); i++) {
        if (nested_features[i] != NO_FEATURE) {
            const vector<handle_t>& walk = variants[i].alt_walk;
            WalkInterval interval;
            interval.walk.assign(walk.begin() + 1, walk.end() - 1);
            interval.length = interval.walk_length(graph);
            features[nested_features[i]] = make_graph_feature(graph, gref, std::move(interval),
                                                              { catalog.id(i), "0", "+" });
        }
    }
    return nested_features;
}

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
 * Return INFO/GR, the gref placement, and INFO/PC, the path coordinate, for
 * each ALT allele of a VCF.
 */
static vector<VcfAlleleInfo> lifted_tags(const vector<GraphFeature>& features, const vector<size_t>& allele_features) {
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

/**
 * Return the INFO/GR and INFO/PC entries of the record of each nested variant,
 * each starting with ";", or empty if the variant has neither.
 */
static vector<string> nested_tags(const NestedVariantCatalog& catalog, const vector<GraphFeature>& features,
                                  const vector<size_t>& nested_features) {
    vector<string> tags(catalog.variants().size());
    for (size_t i = 0; i < tags.size(); i++) {
        if (nested_features[i] != NO_FEATURE) {
            auto [placement, coordinate] = lifted_values(features[nested_features[i]]);
            tags[i] = (placement == "." ? "" : ";GR=" + placement) + (coordinate == "." ? "" : ";PC=" + coordinate);
        }
    }
    return tags;
}

/**
 * Write the VCF to standard output with the tags of its ALT alleles, followed
 * by the records of the nested variants in the catalog, if there is one.
 *
 * @param nested_info Extra INFO for the record of each nested variant.
 */
static void write_lifted_vcf(const string& vcf_file, const VcfCalls& calls, const vector<VcfAlleleInfo>& tags,
                             const NestedVariantCatalog* catalog, vector<string> nested_info,
                             const gbwtgraph::GBWTGraph& graph, const Logger& logger) {
    // The genotyper writes the nested records to a temporary file as it goes,
    // and they follow the calls.
    string nested_file = temp_file::create("vg-lift-nested");
    ofstream nested_out(nested_file);
    optional<NestedGenotyper> genotyper;
    vector<string> nested_header;
    if (catalog != nullptr && !catalog->variants().empty()) {
        genotyper.emplace(*catalog, calls, std::move(nested_info), graph, nested_out);
        nested_header = genotyper->header_lines();
    }
    try {
        write_vcf_with_allele_info(vcf_file, cout, tags, nested_header, [&](size_t record, const VcfCursor& vcf) {
            if (genotyper) {
                genotyper->add_call(record, vcf);
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

int main_lift(int argc, char** argv) {
    Logger logger("vg lift");
    if (argc == 2) {
        help_lift(argv);
        return 1;
    }

    string tigs_file;
    string canonical_file;
    vector<string> priority_prefixes;
    string assemblies_file;
    vector<string> path_prefixes;
    double min_match = 0;
    bool blocks = false;
    bool progress = false;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
            {"tigs", required_argument, 0, 'u'},
            {"canonical", required_argument, 0, 'c'},
            {"priority", required_argument, 0, 'P'},
            {"assemblies", required_argument, 0, 'a'},
            {"min-match", required_argument, 0, 'm'},
            {"blocks", no_argument, 0, 'b'},
            {"paths-by", required_argument, 0, 'Q'},
            {"threads", required_argument, 0, 't'},
            {"progress", no_argument, 0, 'p'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "u:c:P:a:m:bQ:t:ph?", long_options, &option_index);
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
            priority_prefixes.push_back(optarg);
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
            path_prefixes.push_back(optarg);
            break;
        case 't':
            set_thread_count(logger, optarg);
            break;
        case 'p':
            progress = true;
            break;
        case 'h':
        case '?':
        default:
            help_lift(argv);
            return 1;
        }
    }
    if (optind + 2 != argc) {
        logger.error() << "expected a GBZ file and a BED or VCF file" << endl;
    }
    string graph_file = require_exists(logger, argv[optind]);
    string features_file = argv[optind + 1];
    bool vcf_input = features_file != "-" && is_variant_file(features_file);

    gbwtgraph::GBZ gbz;
    load_gbz(gbz, graph_file, progress);
    const gbwtgraph::GBWTGraph& graph = gbz.graph;

    if (!tigs_file.empty()) {
        // Tigs from another graph would name the wrong walks.
        require_compatible_reference(tigs_file, nullptr, &gbz, false);
    }

    GrefIndex gref(gbz, GrefCover::gref_prefix);
    const GrefIndex* gref_pointer = gref.path_count() > 0 ? &gref : nullptr;
    if (progress) {
        cerr << "[vg lift] gref: " << gref.path_count() << " paths, " << gref.total_length() << " bp" << endl;
    }

    vector<GraphFeature> features;
    // For a VCF, its calls, the variants that their ALT alleles carry along
    // gref alt contigs, and the index of the feature of each ALT allele and of
    // each such variant.
    VcfCalls calls;
    optional<NestedVariantCatalog> catalog;
    vector<size_t> allele_features;
    vector<size_t> nested_features;
    if (vcf_input) {
        try {
            calls = read_vcf_calls(features_file, graph, vg::get_thread_count());
        } catch (const runtime_error& e) {
            logger.error() << e.what() << endl;
        }
        allele_features = add_allele_features(calls, graph, gref_pointer, features, logger, progress);
        ContigDecomposer decomposer(graph, gref);
        if (has_contig_records(calls, decomposer, gref)) {
            logger.warn() << features_file << " already has records on gref alt contigs, so no nested variants "
                          << "are added" << endl;
        } else if (gref_pointer != nullptr) {
            catalog.emplace(calls, decomposer, gref, graph);
            nested_features = add_nested_features(*catalog, graph, gref_pointer, features);
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
    } else {
        features = lift_bed_features(features_file, tigs_file, graph_file, gbz, gref_pointer, logger, progress);
    }

    auto feature_id = [](size_t i) {
        return "f" + to_string(i + 1);
    };
    // Write BED columns 4 to 6 (name, score and strand) for a place of the
    // feature. The place reads the feature in reverse if is_reverse is set.
    auto write_name_score_strand = [](ostream& out, const GraphFeature& feature, bool is_reverse) {
        const vector<string>& fields = feature.fields;
        out << '\t' << (fields.size() > 0 ? fields[0] : ".") << '\t' << (fields.size() > 1 ? fields[1] : "0") << '\t';
        if (feature.is_stranded()) {
            out << (is_reverse ? '-' : '+');
        } else {
            out << '.';
        }
    };
    // Write the feature's BED columns from the 7th on, then its ID.
    auto write_rest = [&](ostream& out, size_t i) {
        const vector<string>& fields = features[i].fields;
        for (size_t j = BED_STRAND_FIELD + 1; j < fields.size(); j++) {
            out << '\t' << fields[j];
        }
        out << '\t' << feature_id(i);
    };

    // One pass over the paths other than gref paths finds the places for -a and
    // names the features off the gref. Naming uses every path, whatever -Q says.
    size_t unnamed = 0;
    for (const GraphFeature& feature : features) {
        unnamed += !feature.coordinate;
    }
    if (unnamed > 0 || !assemblies_file.empty()) {
        PathPriority priority(gbz, priority_prefixes, GrefCover::gref_prefix, gref.path_count());
        ofstream out;
        if (!assemblies_file.empty()) {
            out.open(assemblies_file);
            if (!out) {
                logger.error() << "could not open " << assemblies_file << " for writing" << endl;
            }
        }
        auto is_gref = [&](const string& name) {
            return name.compare(0, GrefCover::gref_prefix.size(), GrefCover::gref_prefix) == 0;
        };
        auto selected = [&](const string& name) {
            if (path_prefixes.empty()) {
                return true;
            }
            for (const string& prefix : path_prefixes) {
                if (name.compare(0, prefix.size(), prefix) == 0) {
                    return true;
                }
            }
            return false;
        };
        size_t written = 0;
        auto visit = [&](const string& name, const vector<PathOccurrence>& places) {
            bool write = out.is_open() && selected(name);
            size_t rank = priority.rank(name);
            for (const PathOccurrence& place : places) {
                GraphFeature& feature = features[place.feature];
                size_t length = feature.interval.length;
                // A place that holds the whole feature in one block is a path that carries it.
                if (place.matched == length && place.blocks.size() == 1 && place.end - place.start == length) {
                    PathCoordinate here { name, rank, place.start, place.end, place.is_reverse };
                    if (!feature.coordinate || here < *feature.coordinate) {
                        feature.coordinate = std::move(here);
                    }
                }
                if (!write) {
                    continue;
                }
                out << name << '\t' << place.start << '\t' << place.end;
                write_name_score_strand(out, feature, place.is_reverse);
                if (blocks) {
                    out << '\t' << place.start << '\t' << place.end << "\t0\t" << place.blocks.size() << '\t';
                    for (const auto& block : place.blocks) {
                        out << (block.second - block.first) << ',';
                    }
                    out << '\t';
                    for (const auto& block : place.blocks) {
                        out << (block.first - place.start) << ',';
                    }
                }
                write_rest(out, place.feature);
                if (min_match > 0) {
                    out << '\t' << double(place.matched) / length;
                }
                out << '\n';
                written++;
            }
        };
        auto not_gref = [&](const string& name) {
            return !is_gref(name);
        };
        if (min_match > 0) {
            for_each_path_liftover(gbz, features, min_match, not_gref, visit);
        } else {
            for_each_path_occurrences(gbz, features, not_gref, visit);
        }
        if (progress && out.is_open()) {
            cerr << "[vg lift] " << written << " places on paths" << endl;
        }
    }

    if (vcf_input) {
        vector<string> nested_info;
        if (catalog) {
            nested_info = nested_tags(*catalog, features, nested_features);
        }
        write_lifted_vcf(features_file, calls, lifted_tags(features, allele_features), catalog ? &*catalog : nullptr,
                         std::move(nested_info), graph, logger);
    } else {
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
    cout.flush();

    if (!canonical_file.empty()) {
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
            write_rest(out, i);
            out << '\n';
        }
    }
    size_t nameless = 0;
    for (const GraphFeature& feature : features) {
        nameless += !feature.coordinate;
    }
    if (nameless > 0) {
        logger.warn() << nameless << " features are not carried whole by any path, so they have no path coordinate"
                      << endl;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_lift("lift", "lift features on a cover or on paths to the graph, gref and assemblies", WIDGET, main_lift);
