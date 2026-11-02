/** \file call_main.cpp
 *
 * Defines the "vg call" subcommand, which calls variation from an augmented graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <regex>
#include <algorithm>
#include <list>
#include <fstream>
#include <sstream>

#include "subcommand.hpp"
#include "../path.hpp"
#include "../graph_caller.hpp"
#include "../flow_caller.hpp"
#include "../legacy_caller.hpp"
#include "../../artifacts/multipass_caller.hpp"
#include "../snarl_manager_decomposition.hpp"
#include "../nested_flow_caller.hpp"
#include "../vcf_genotyper.hpp"
#include "../sorted_vcf_writer.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../xg.hpp"
#include "../gbzgraph.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gref.hpp"
#include "../traversal_clusters.hpp"
#include "../../artifacts/read_likelihood_caller.hpp"
#include "../../artifacts/site_read_source.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

const string DEFAULT_SAMPLE_NAME = "SAMPLE";

// The default --read-window, in node IDs, of each read source that fetches reads by node-ID
// window. The window size changes the order in which a site sees its reads, and so the last
// digits of likelihoods summed over them.

/// A GAM index query is a seek in an open file, so its window is narrow, to avoid fetching reads
/// that no site needs.
const size_t DEFAULT_GAM_INDEX_WINDOW = 256;
/// A GAF-Base query decodes whole blocks of alignments in a gaf-base process, so its window is
/// wide, to serve many sites per query. It is also much wider than the node-ID span of a long read,
/// so that few reads cross a window boundary and are fetched twice.
const size_t DEFAULT_GAF_BASE_WINDOW = 16384;
/// An indexed GAF query is a seek, but its window is as wide as that of GAF-Base, so that the two
/// fetch the same windows, and few long reads cross a window boundary and are decoded twice.
const size_t DEFAULT_GAF_INDEX_WINDOW = 16384;

/// Count the haplotypes of the HAPLOTYPE-sense paths of a graph.
///
/// A haplotype is a sample name and haplotype number, so a haplotype stored as several paths
/// counts once. The count bounds the number of non-reference alleles that haplotype-based allele
/// enumeration can offer at a site.
static size_t count_panel_haplotypes(const PathHandleGraph& graph) {
    set<pair<string, size_t>> haplotypes;
    graph.for_each_path_of_sense(PathSense::HAPLOTYPE, [&](const path_handle_t& path) {
        haplotypes.emplace(graph.get_sample_name(path), graph.get_haplotype(path));
    });
    return haplotypes.size();
}

/// Join items into one list, separated by ", ".
static string join_with_commas(const vector<string>& items) {
    stringstream joined;
    for (size_t i = 0; i < items.size(); ++i) {
        joined << (i ? ", " : "") << items[i];
    }
    return joined.str();
}

/// The long options that configure one subsystem of vg call, and the option that turns the
/// subsystem on. Options of a subsystem that is off would have no effect, so we refuse them.
struct OptionGroup {
    /// The option that turns the subsystem on, or null for options that always apply.
    const char* switch_option;
    /// The name of the subsystem that contains this one, or null. A subsystem can only be on when
    /// the one containing it is.
    const char* within;
    /// Whether the subsystem is on, once all the options are parsed.
    function<bool()> on;
    vector<struct option> options;
};

/**
 * Refuse the options given for a subsystem that is off.
 *
 * An option is refused with the outermost subsystem that contains it and is off, since that is
 * the one to turn on first. Only the first such subsystem is reported, since logger.error() exits.
 *
 * @param groups The option group of each subsystem, by the subsystem's name.
 * @param options The option numbers of the options given.
 * @param option_name The long name of each option, by its option number.
 * @param option_subsystem The subsystem each option configures, by its option number.
 */
static void refuse_inactive_subsystem_options(const Logger& logger,
                                              const map<string, OptionGroup>& groups,
                                              const unordered_set<int>& options,
                                              const unordered_map<int, string>& option_name,
                                              const unordered_map<int, string>& option_subsystem) {
    // The options to refuse, by the subsystem they are refused with.
    map<string, vector<string>> refused;
    for (int number : options) {
        const char* outermost_off = nullptr;
        for (const char* subsystem = option_subsystem.at(number).c_str(); subsystem != nullptr;
             subsystem = groups.at(subsystem).within) {
            if (!groups.at(subsystem).on()) {
                outermost_off = subsystem;
            }
        }
        if (outermost_off != nullptr) {
            refused[outermost_off].push_back("--" + option_name.at(number));
        }
    }
    for (auto& [subsystem, names] : refused) {
        std::ranges::sort(names);
        logger.error() << join_with_commas(names)
                       << (names.size() == 1 ? " only applies to " : " only apply to ") << subsystem
                       << ", which " << groups.at(subsystem).switch_option << " turns on" << endl;
    }
}

/// The reads given to the read-likelihood genotyper, and how to fetch them. A file name is empty
/// when the file was not given.
struct ReadInputOptions {
    /// Reads in GAM (--gam).
    string gam;
    /// The index of the GAM (--gam-index).
    string gam_index;
    /// Reads in GAF (--gaf-reads).
    string gaf;
    /// The tabix index of the sorted GAF (--gaf-index).
    string gaf_index;
    /// A GAF-Base database of reads (--gaf-base).
    string gaf_base;
    /// The gaf-base program, which answers GAF-Base queries (--gaf-base-binary).
    string gaf_base_binary = "gaf-base";
    /// The width, in node IDs, of the windows that reads are fetched by (--read-window), or 0 for
    /// the default of the read source.
    size_t window_size = 0;

    /// Get the name of the file the reads come from.
    const string& reads_file() const {
        if (!gaf_base.empty()) {
            return gaf_base;
        }
        if (!gaf.empty()) {
            return gaf;
        }
        return gam;
    }
};

/// Return the default --read-window of the read source that `reads` selects, or 0 for a source
/// that loads every read at once.
static size_t default_read_window(const ReadInputOptions& reads) {
    if (!reads.gaf_base.empty()) {
        return DEFAULT_GAF_BASE_WINDOW;
    }
    if (!reads.gaf_index.empty()) {
        return DEFAULT_GAF_INDEX_WINDOW;
    }
    if (!reads.gam_index.empty()) {
        return DEFAULT_GAM_INDEX_WINDOW;
    }
    return 0;
}

/**
 * Open the read source that `reads` selects.
 *
 * GAF-Base, an indexed GAF and an indexed GAM fetch reads as the sites need them, in windows of
 * `reads.window_size` node IDs. Otherwise every read is loaded into memory now. A source that
 * cannot be set up is a user error, which exits.
 */
static unique_ptr<SiteReadSource> open_site_read_source(const PathHandleGraph& graph,
                                                        const ReadInputOptions& reads,
                                                        const SiteReadFilter& read_filter,
                                                        bool show_progress, const Logger& logger) {
    // Each windowed source keeps four windows per thread, in one cache that all threads share.
    const size_t cached_windows = 4 * (size_t)vg::get_thread_count();
    if (!reads.gaf_base.empty()) {
        // Extra arguments for each gaf-base process, such as the sizes of its caches, separated by
        // spaces. They come from an environment variable until the right settings are known.
        vector<string> batch_arguments;
        if (const char* extra = getenv("VG_GAFBASE_BATCH_ARGS"); extra != nullptr) {
            istringstream words(extra);
            for (string word; words >> word;) {
                batch_arguments.push_back(word);
            }
        }
        unique_ptr<GafBaseSiteReadSource> gaf_base_source;
        // Check the setup now, so that a missing binary or unreadable database is reported as a
        // user error rather than as a crash.
        try {
            gaf_base_source = make_unique<GafBaseSiteReadSource>(graph, reads.gaf_base, read_filter,
                                                                 reads.window_size, cached_windows,
                                                                 reads.gaf_base_binary,
                                                                 batch_arguments);
            gaf_base_source->check_setup();
        } catch (const std::exception& e) {
            logger.error() << e.what() << endl;
        }
        if (show_progress) {
            string extra;
            for (const string& argument : batch_arguments) {
                extra += " " + argument;
            }
            logger.info() << "Using GAF-Base " << reads.gaf_base << " via '"
                          << reads.gaf_base_binary << " batch" << extra << "', one process per thread"
                          << endl;
            logger.info() << "GAF-Base: database opened by gaf-base "
                          << (gaf_base_source->immutable_database()
                              ? "as a SQLite immutable URI, without file locking"
                              : "as a plain path, with SQLite file locking (VG_GAFBASE_LOCKING=1)")
                          << endl;
        }
        return gaf_base_source;
    }
    if (!reads.gaf_index.empty()) {
        unique_ptr<SiteReadSource> source;
        try {
            source = make_unique<TabixGafSiteReadSource>(graph, reads.gaf, reads.gaf_index,
                                                         read_filter, reads.window_size,
                                                         cached_windows);
        } catch (const std::exception& e) {
            logger.error() << e.what() << endl;
        }
        if (show_progress) {
            logger.info() << "Using indexed GAF " << reads.gaf
                          << " with index " << reads.gaf_index << endl;
        }
        return source;
    }
    if (!reads.gam_index.empty()) {
        auto source = make_unique<IndexedGamSiteReadSource>(reads.gam, reads.gam_index, read_filter,
                                                            reads.window_size, cached_windows);
        if (show_progress) {
            logger.info() << "Using indexed GAM " << reads.gam
                          << " with index " << reads.gam_index << endl;
        }
        return source;
    }
    auto in_memory_source = make_unique<InMemorySiteReadSource>();
    if (!reads.gam.empty()) {
        in_memory_source->load_gam(reads.gam, read_filter);
    } else {
        in_memory_source->load_gaf(graph, reads.gaf, read_filter);
    }
    if (show_progress) {
        logger.info() << "Loaded " << in_memory_source->get_read_count()
                      << " reads (" << in_memory_source->get_filtered_count()
                      << " filtered out)" << endl;
    }
    return in_memory_source;
}

/// Log how the fetches and the cache of a windowed read source performed.
static void log_windowed_read_source_stats(const Logger& logger,
                                           const WindowedSiteReadSource& windowed, bool off_panel) {
    const WindowedSiteReadSource::FetchStats stats = windowed.fetch_stats();
    const size_t hits = stats.cache_hits;
    const size_t misses = stats.cache_misses;
    const size_t total = hits + misses;
    const auto* gaf_base = dynamic_cast<const GafBaseSiteReadSource*>(&windowed);
    const auto* tabix_gaf = dynamic_cast<const TabixGafSiteReadSource*>(&windowed);
    const char* label = "Indexed GAM: ";
    if (gaf_base != nullptr) {
        label = "GAF-Base: ";
    } else if (tabix_gaf != nullptr) {
        label = "Indexed GAF: ";
    }
    logger.info() << label
                  << stats.fetched << " reads fetched, "
                  << hits << "/" << total << " site queries served from cache"
                  << (total > 0 ? " (" + std::to_string((int)(100.0 * hits / total)) + "%)" : "")
                  << endl;
    // Reads examined by site queries, against reads delivered to sites.
    const size_t seen = stats.scanned;
    const size_t used = stats.delivered;
    logger.info() << label
                  << seen << " read candidates examined, " << used << " delivered"
                  << (seen > 0 ? " (" + std::to_string((int)(100.0 * used / seen)) + "%)" : "")
                  << endl;
    // A site across a few windows is served from those windows. A wider one is fetched uncached,
    // by its exact node ranges.
    logger.info() << label
                  << stats.wide_from_windows << " site queries across windows "
                  << "answered from the cached windows" << endl;
    logger.info() << label
                  << stats.straddles << " site queries too wide for a "
                  << "window, fetched uncached over " << stats.straddle_wanted
                  << " node IDs (spanning " << stats.straddle_nodes << ")"
                  << endl;
    logger.info() << label
                  << stats.whole_fetches << " windows fetched whole, holding "
                  << stats.whole_fetch_reads << " reads" << endl;
    if (off_panel) {
        logger.info() << label << stats.refetches
                      << " of the whole fetches were of a window the cache had dropped, and "
                      << stats.edit_fetches
                      << " only to count a window's off-panel candidates; "
                      << stats.held_queries
                      << " site and read-start queries were answered from reads held for edit sites"
                      << endl;
    }
    if (gaf_base != nullptr) {
        auto totals = gaf_base->get_query_totals();
        logger.info() << "GAF-Base: " << totals.processes << " gaf-base processes answered "
                      << totals.queries << " queries; " << gaf_base->get_skipped_count()
                      << " records lay between a fetch's ranges without touching its nodes; "
                      << gaf_base->get_duplicate_count() << " duplicate reads dropped" << endl;
        // Where the time went. The processes' CPU and memory are known once they have stopped.
        constexpr double GIB = 1024.0 * 1024.0 * 1024.0;
        logger.info() << "GAF-Base: threads waited " << totals.wait_s / 3600 << " h for answers and "
                      << "spent " << totals.parse_s / 3600 << " h parsing their "
                      << totals.gaf_bytes / 1e9 << " GB of GAF; the processes used "
                      << totals.child_user_s / 3600 << " h of CPU in their own code and "
                      << totals.child_system_s / 3600 << " h in the kernel, and at most "
                      << totals.child_max_rss / GIB << " GiB of memory each ("
                      << totals.child_rss_sum / GIB << " GiB added up)" << endl;
    }
    if (tabix_gaf != nullptr) {
        logger.info() << "Indexed GAF: " << tabix_gaf->get_query_count() << " index lookups; "
                      << tabix_gaf->get_skipped_count() << " records overlapped a fetch's node "
                      << "interval without touching its nodes; "
                      << tabix_gaf->get_duplicate_count() << " duplicate reads dropped" << endl;
        logger.info() << "Indexed GAF: of those records, " << tabix_gaf->get_unparsed_count()
                      << " were dropped from their path text without being parsed; threads spent "
                      << tabix_gaf->get_read_seconds() / 3600 << " h reading "
                      << tabix_gaf->get_gaf_bytes() / 1e9 << " GB of GAF through the index and "
                      << tabix_gaf->get_parse_seconds() / 3600 << " h parsing it" << endl;
    }
}

/// Return the gRef path with a step on `handle`, or "" if there is none.
static string gref_path_on(const PathHandleGraph& graph, handle_t handle) {
    string owner;
    graph.for_each_step_on_handle(handle, [&](const step_handle_t& step) {
        if (const string n = graph.get_path_name(graph.get_path_handle_of_step(step));
            GrefCover::is_gref_derived(n)) {
            owner = n;
            return false;
        }
        return true;
    });
    return owner;
}

/**
 * Return a gRef path, other than `path_name` itself, that has a node next to an end of path
 * `path_name`, or "" if there is none.
 *
 * For a fragment, either end will do. A fragment is a maximal run of uncovered nodes, so its
 * neighbours at both ends are boundary nodes of the same enclosing snarl.
 */
static string gref_neighbor_of(const PathHandleGraph& graph, const string& path_name) {
    string neighbor;
    if (!graph.has_path(path_name)) {
        return neighbor;
    }
    const path_handle_t p = graph.get_path_handle(path_name);
    if (graph.is_empty(p)) {
        return neighbor;
    }
    for (bool left : {true, false}) {
        const handle_t end = graph.get_handle_of_step(left ? graph.path_begin(p) : graph.path_back(p));
        graph.follow_edges(end, left, [&](const handle_t& next) {
            if (const string owner = gref_path_on(graph, next); !owner.empty() && owner != path_name) {
                neighbor = owner;
                return false;
            }
            return true;
        });
        if (!neighbor.empty()) {
            break;
        }
    }
    return neighbor;
}

/**
 * Return how many levels of non-reference sequence a path lies in.
 *
 * A base contig, or the gRef copy of one, is level 0. A fragment is one level deeper than the gRef
 * path next to it, since the paths of a gRef cover share no nodes. A fragment with no gRef path
 * next to it is level 1, as when the fragment next to it was too short to be written.
 *
 * @param by_path The levels found so far, by path name.
 * @param in_progress The paths whose levels are being found.
 */
static int gref_level_of(const PathHandleGraph& graph, const string& path_name,
                         map<string, int>& by_path, unordered_set<string>& in_progress) {
    // A clipped fragment is a subpath, so its name ends in [offset].
    if (!GrefCover::is_gref_name(Paths::strip_subrange(path_name))) {
        return 0;
    }
    if (auto seen = by_path.find(path_name); seen != by_path.end()) {
        return seen->second;
    }
    if (!in_progress.insert(path_name).second) {
        // A cover is a tree, so this means a broken cover. Stop rather than recurse forever.
        return 1;
    }
    const string neighbor = gref_neighbor_of(graph, path_name);
    const int level = neighbor.empty() ? 1 : gref_level_of(graph, neighbor, by_path, in_progress) + 1;
    in_progress.erase(path_name);
    by_path[path_name] = level;
    return level;
}

/// Return the gRef nesting level of each reference path above level 0, by the contig name the
/// VCF uses for it.
static map<string, int> gref_contig_levels(const PathHandleGraph& graph,
                                           const vector<string>& ref_paths) {
    map<string, int> gref_levels;
    map<string, int> by_path;
    unordered_set<string> in_progress;
    for (const string& ref_path : ref_paths) {
        const int level = gref_level_of(graph, ref_path, by_path, in_progress);
        if (level > 0) {
            const string locus = PathMetadata::parse_locus_name(ref_path);
            gref_levels[locus != PathMetadata::NO_LOCUS_NAME ? locus : ref_path] = level;
        }
    }
    return gref_levels;
}

/// A haplotype of the panel: a GBWT sample number and a phase.
using PanelHaplotype = pair<size_t, size_t>;

/// Return whether the linkage panel leaves out a GBWT path, because it is part of a gRef cover.
///
/// The fragment paths of a cover all share one sample and phase, so they would make one haplotype
/// stitched together from many donors. The gRef copy of a base reference path is left out only if
/// the sample of the base path is in `base_samples`, so that the reference is in the panel once.
static bool is_left_out_of_panel(const gbwt::Metadata& meta, const gbwt::PathName& name,
                                 const unordered_set<string>& base_samples) {
    const string path_sample = (size_t)name.sample < meta.sample_names.size()
                                   ? meta.sample(name.sample) : string();
    const string path_contig = (size_t)name.contig < meta.contig_names.size()
                                   ? meta.contig(name.contig) : string();
    const bool gref_fragment = GrefCover::is_gref_derived(path_sample)
                               && GrefCover::is_gref_name(path_contig);
    const bool gref_shadowing_base =
        GrefCover::is_gref_derived(path_sample) && !gref_fragment
        && base_samples.contains(path_sample.substr(GrefCover::gref_prefix.size()));
    return gref_fragment || gref_shadowing_base;
}

/**
 * Number the haplotypes of the panel in a GBWT, and give each GBWT sequence the number of its
 * haplotype.
 *
 * A haplotype may be stored as several paths, one per contig or fragment, each with a sequence per
 * orientation. They all get the number of the haplotype, so a haplotype counts once. The sequences
 * of paths left out of the panel get LinkageModel::NO_HAPLOTYPE.
 *
 * @param sequence_to_haplotype Filled with the haplotype number of each GBWT sequence.
 * @return The number of each haplotype.
 */
static map<PanelHaplotype, size_t> number_panel_haplotypes(const gbwt::GBWT& gbwt_index,
                                                           vector<size_t>& sequence_to_haplotype) {
    const gbwt::Metadata& meta = gbwt_index.metadata;
    // The samples that are not gRef-derived.
    unordered_set<string> base_samples;
    for (gbwt::size_type i = 0; i < meta.sample_names.size(); ++i) {
        const string s = meta.sample(i);
        if (!GrefCover::is_gref_derived(s)) {
            base_samples.insert(s);
        }
    }
    map<PanelHaplotype, size_t> hap_index;
    sequence_to_haplotype.assign(gbwt_index.sequences(), 0);
    for (gbwt::size_type path = 0; path < meta.paths(); ++path) {
        const gbwt::PathName& name = meta.path(path);
        size_t index = LinkageModel::NO_HAPLOTYPE;
        if (!is_left_out_of_panel(meta, name, base_samples)) {
            // try_emplace numbers a new haplotype and finds an old one.
            index = hap_index.try_emplace(make_pair((size_t)name.sample, (size_t)name.phase),
                                          hap_index.size()).first->second;
        }
        for (gbwt::size_type orientation = 0; orientation < 2; ++orientation) {
            gbwt::size_type seq = gbwt::Path::encode(path, orientation);
            if (seq < sequence_to_haplotype.size()) {
                sequence_to_haplotype[seq] = index;
            }
        }
    }
    return hap_index;
}

/// Return the name of each panel haplotype, as "sample#phase", by its number.
static vector<string> panel_haplotype_names(const gbwt::Metadata& meta,
                                            const map<PanelHaplotype, size_t>& hap_index) {
    vector<string> hap_names(hap_index.size());
    for (const auto& [haplotype, index] : hap_index) {
        const auto& [sample_number, phase] = haplotype;
        string sample = sample_number < meta.sample_names.size()
                            ? meta.sample(sample_number)
                            : string("sample") + std::to_string(sample_number);
        hap_names[index] = sample + "#" + std::to_string(phase);
    }
    return hap_names;
}

/// One option that a --preset sets.
struct PresetSetting {
    /// The long name of the option, without "--".
    const char* option;
    /// The option's argument, or null for a switch.
    const char* value;
};

/// The settings of each --preset, by preset name.
///
/// help_call lists these and main_call applies them, so the help lists exactly the settings
/// applied. The values for `ont` were fitted to Oxford Nanopore reads, and those for `hifi` to
/// PacBio HiFi reads (HG002, Revio).
static const vector<pair<string, vector<PresetSetting>>> PRESETS = {
    {"ont", {{"read-phasing", nullptr}, {"regenotype", nullptr}, {"gap-open", "1"},
             {"gap-extend", "1"}, {"mismap-min", "0.05"}, {"read-min-mapq", "5"},
             {"insertion-nats", "0.9"}, {"hp-prior", "20"},
             {"off-panel-indels", "7,0.25,5"}}},
    {"hifi", {{"read-phasing", nullptr}, {"regenotype", nullptr}, {"gap-open", "1"},
              {"gap-extend", "1"}, {"mismap-min", "0.05"}, {"read-min-mapq", "5"},
              {"insertion-nats", "-1.8"}, {"hp-prior", "0"}, {"linkage-weight", "1"},
              {"off-panel-indels", "7,0.25,5"}}},
};

/// Return the settings of each preset for the --preset help, wrapped at the description column of
/// the help.
static string preset_help_lines() {
    const string indent(28, ' ');
    const size_t width = 80;
    stringstream out;
    for (const auto& [name, settings] : PRESETS) {
        string line = indent + name + ":";
        for (const PresetSetting& setting : settings) {
            string item = string("--") + setting.option;
            if (setting.value != nullptr) {
                item += string(" ") + setting.value;
            }
            if (line.size() + 1 + item.size() > width) {
                out << line << endl;
                line = indent + item;
            } else {
                line += " " + item;
            }
        }
        out << line << endl;
    }
    return out.str();
}

void help_call(char** argv) {
    cerr << "usage: " << argv[0] << " call [options] <graph> > output.vcf" << endl
         << "Call variants or genotype known variants" << endl
         << endl
         << "support calling options:" << endl
         << "  -k, --pack FILE           supports created from vg pack for given input graph" << endl
         << "  -m, --min-support M,N     min allele (M) and site (N) support to call [2,4]" << endl
         << "  -e, --baseline-error X,Y  baseline error rates for Poisson model for small (X)" << endl
         << "                            and large (Y) variants [0.005,0.01]" << endl
         << "  -B, --bias-mode           use old ratio-based genotyping algorithm" << endl
         << "                            as opposed to probablistic model" << endl
         << "  -b, --het-bias M,N        homozygous alt/ref allele must have >= M/N times" << endl
         << "                            more support than the next best allele [6,6]" << endl
         << "read-likelihood calling options (all need --read-likelihood):" << endl
         << "      --read-likelihood     genotype from the likelihood of the reads under" << endl
         << "                            each genotype, instead of from pack support" << endl
         << "      --preset NAME         set options to values suited to a read type;" << endl
         << "                            options given explicitly keep their values:" << endl
         << preset_help_lines()
         << "      --enumerate-support   take candidate alleles from read support, which" << endl
         << "                            needs -k, rather than from the GBZ haplotypes" << endl
         << "  read input (give one of --gam, --gaf-reads or --gaf-base; for many reads," << endl
         << "  such as a whole genome, --gaf-reads with --gaf-index is recommended):" << endl
         << "      --gam FILE            read alignments in GAM format" << endl
         << "      --gaf-reads FILE      read alignments in GAF format, all loaded into" << endl
         << "                            memory unless --gaf-index is given" << endl
         << "      --gam-index FILE      index of --gam (from vg gamsort -i), to fetch" << endl
         << "                            reads as sites need them instead of loading them" << endl
         << "                            all" << endl
         << "      --gaf-index FILE      tabix index of --gaf-reads, to fetch reads from the" << endl
         << "                            GAF file itself as sites need them, with no" << endl
         << "                            GAF-Base database. --gaf-reads must be sorted" << endl
         << "                            (vg gamsort -G), compressed with bgzip and" << endl
         << "                            indexed with tabix -p gaf" << endl
         << "      --gaf-base FILE       GAF-Base database of read alignments (gaf-base" << endl
         << "                            construct, from a sorted GAF), queried as sites" << endl
         << "                            need them through gaf-base batch, which must be" << endl
         << "                            on the PATH" << endl
         << "      --gaf-base-binary P   gaf-base executable to run [gaf-base]" << endl
         << "      --read-window N       node-ID window for indexed read fetches" << endl
         << "                            [16384 for --gaf-base and --gaf-index, 256 for" << endl
         << "                            --gam-index]" << endl
         << "      --read-min-mapq N     ignore reads with MAPQ below N [0]" << endl
         << "  nested snarls:" << endl
         << "      --nested              genotype nested snarls in records of their own," << endl
         << "                            not as part of the enclosing snarl's alleles [on]" << endl
         << "      --no-nested           report variation inside nested snarls as part of" << endl
         << "                            the enclosing snarl's alleles" << endl
         << "      --atomize-blocks      write one record per separate difference between" << endl
         << "                            a called allele and the reference, rather than one" << endl
         << "                            per snarl [on]" << endl
         << "      --no-atomize-blocks   write one record per snarl" << endl
         << "  read scoring:" << endl
         << "      --gap-open N          gap-open penalty for scoring reads [6]" << endl
         << "      --gap-extend N        gap-extension penalty for scoring reads [1]" << endl
         << "      --insertion-nats X    add X nats to a read's log-likelihood for each gap" << endl
         << "                            where the read has bases the allele lacks [0]" << endl
         << "      --optimal-pairing     pair a read's node visits with an allele's" << endl
         << "                            optimally, rather than greedily" << endl
         << "      --no-optimal-pairing  pair them greedily (default)" << endl
         << "      --no-mismap-term      do not model the chance that a read is mismapped" << endl
         << "      --mismap-max P        cap on the mismapping probability derived from a" << endl
         << "                            read's MAPQ [0.95]" << endl
         << "      --mismap-min P        floor on any read's mismapping probability [0.02]" << endl
         << "      --depth-term W        weight of the term ln P(read count | genotype) in" << endl
         << "                            the likelihood; 0 disables it [0.1]" << endl
         << "      --depth-count-raw     count each read as one in the read-count term," << endl
         << "                            instead of as its chance of being correctly mapped" << endl
         << "  linkage between sites:" << endl
         << "      --linkage-weight W    weight of the linkage model, which favours runs of" << endl
         << "                            genotypes that the -z/-g haplotypes carry; 0" << endl
         << "                            disables it [2]" << endl
         << "      --linkage-scale N     distance scale of linkage decay, in bp [10000]" << endl
         << "      --linkage-prior F     exponent on the haplotypes' allele-frequency prior" << endl
         << "                            [5]" << endl
         << "      --linkage-mutation P  chance that a strand carries another allele than the" << endl
         << "                            haplotype it copies, so that an allele no haplotype" << endl
         << "                            carries is a mutation on one ["
         << LinkageModel::Params().mutation << "]" << endl
         << "      --hp-prior F          exponent used instead of --linkage-prior where an" << endl
         << "                            allele differs from the reference only in the" << endl
         << "                            length of a homopolymer run; 0 disables it [0]" << endl
         << "      --hp-prior-run N      shortest homopolymer run for --hp-prior [11]" << endl
         << "      --no-phased           write unphased genotypes, without FORMAT/PS" << endl
         << "      --phased              write phased genotypes, failing if the linkage" << endl
         << "                            model cannot run [on when the linkage model runs]" << endl
         << "  off-panel variants:" << endl
         << "      --off-panel           also call SNVs and indels of 1-49 bp that no graph" << endl
         << "                            walk spells, found as edits recurring in the reads" << endl
         << "                            [on]" << endl
         << "      --no-off-panel        do not call off-panel variants" << endl
         << "      --off-panel-dump F    write each off-panel candidate, with its read" << endl
         << "                            counts and what became of it, to F as TSV" << endl
         << "      --off-panel-bar B     call a candidate whose ALT is on >= K fragments" << endl
         << "                            and a fraction >= F of them, at GQ >= G and a" << endl
         << "                            depth ratio (DR) <= D (D 0: no limit), with B" << endl
         << "                            given as K,F,G[,D] [7,0.25,0,3]" << endl
         << "      --off-panel-indels B  call an indel candidate whose ALT is on >= K" << endl
         << "                            fragments and a fraction >= F of them, unless its" << endl
         << "                            tract starts with a run of >= H of one base (H 0:" << endl
         << "                            no limit; 1: no indels), with B given as K,F,H" << endl
         << "                            [7,0.25,0]" << endl
         << "      --off-panel-svs       also call insertions and deletions of 50 bp or" << endl
         << "                            more, from long reads' indels (the medoid read's" << endl
         << "                            insertion as the ALT) and short reads' clips" << endl
         << "                            (deletions only); one crossing a site is refused" << endl
         << "                            [off]" << endl
         << "  read-based phasing (reliable heterozygous sites are joined into a phase" << endl
         << "  chain by the reads they share; other sites are phased from the chain):" << endl
         << "      --read-phasing        phase heterozygous sites with the reads that span" << endl
         << "                            them, not only with the haplotypes" << endl
         << "      --no-read-phasing     phase with the haplotypes only (default)" << endl
         << "      --phase-min-q N       minimum mean read confidence (phred) for a site to" << endl
         << "                            join the phase chain" << endl
         << "                            [9.5, or 8.5 under --optimal-pairing]" << endl
         << "      --phase-break N       break the phase chain where the reads linking two" << endl
         << "                            neighbouring sites give under N log10 units of" << endl
         << "                            evidence [20]" << endl
         << "      --phase-relink N      sites of the chain used on each side of a break to" << endl
         << "                            decide how to rejoin it [10]" << endl
         << "      --phase-hang N        sites of the chain used to phase a site outside" << endl
         << "                            it: the nearest N/2+1 on each side [4]" << endl
         << "      --phase-prior N       weight, in log10 units, of the haplotypes' phase" << endl
         << "                            when phasing a site that is not in the chain [3]" << endl
         << "      --phase-cap N         maximum evidence, in log10 units, from the reads" << endl
         << "                            linking two sites; 0 for no maximum [0]" << endl
         << "      --phase-coherence F   minimum fraction of a site's reads that agree with" << endl
         << "                            their phase at other sites, for the site to stay" << endl
         << "                            in the phase chain; 0 disables the check [0.70]" << endl
         << "      --phase-coh-rounds N  maximum rounds of removing incoherent sites from" << endl
         << "                            the phase chain [2]" << endl
         << "  re-genotyping (a read's strand log-odds, in nats, say which haplotype its" << endl
         << "  other heterozygous sites place it on):" << endl
         << "      --regenotype          re-genotype sites after phasing, weighting each" << endl
         << "                            read toward the haplotype its other sites put it" << endl
         << "                            on (needs --read-phasing)" << endl
         << "      --no-regenotype       do not re-genotype sites (default)" << endl
         << "      --regeno-temper N     scale on each read's strand log-odds; 0 leaves the" << endl
         << "                            genotypes unchanged [fitted to the data]" << endl
         << "      --regeno-passes N     rounds to run, each with one linkage pass," << endl
         << "                            counting the first, which comes before any" << endl
         << "                            re-genotyping; stops early if genotypes stop" << endl
         << "                            changing. 1 only reports what re-genotyping" << endl
         << "                            would change [2]" << endl
         << "      --regeno-ceiling N    scale on how far a read's haplotype probability" << endl
         << "                            may move from 1/2; 1 means no limit [1]" << endl
         << "      --no-regeno-haploid   at a nested site that only one haplotype carries," << endl
         << "                            do not down-weight reads whose strand log-odds" << endl
         << "                            place them on the other haplotype" << endl
         << "  anchors:" << endl
         << "      --anchors-out FILE    write assembly anchors to FILE: the reads at each" << endl
         << "                            genotyped snarl's boundaries, grouped by the" << endl
         << "                            called allele they fit" << endl
         << "      --anchors-het-only    write anchors only at heterozygous sites" << endl
         << "      --anchors-leaf-only   write anchors only at leaf snarls" << endl
         << "      --anchors-reads N     minimum reads per anchor [2]" << endl
         << "      --anchors-min-gqn F   minimum site GQN for a site's anchors [0]" << endl
         << "      --anchors-min-q F     minimum confidence (phred) of a read's placement" << endl
         << "                            [0]" << endl
         << "      --anchors-keep-off-call" << endl
         << "                            keep reads whose best-fitting allele was not" << endl
         << "                            called" << endl
         << "      --anchors-end-new N   write a snarl's end anchor only if it has at least" << endl
         << "                            N reads that its start anchor lacks [0]" << endl
         << "      --anchors-hom-split   with --read-phasing, split the reads at a" << endl
         << "                            homozygous site into two haplotypes by their" << endl
         << "                            strand log-odds" << endl
         << "      --split-min-q N       minimum |strand log-odds|, tempered, for a read" << endl
         << "                            to count as confidently placed when splitting [0.5]" << endl
         << "      --split-min-side N    confidently placed reads needed on each haplotype" << endl
         << "                            to split a homozygous site [10]" << endl
         << "      --no-anchors-phase-hets" << endl
         << "                            place reads at a heterozygous site by allele" << endl
         << "                            match alone, not also by strand log-odds" << endl
         << "      --no-off-ref-nesting  do not genotype chains off the reference just to" << endl
         << "                            write their anchors" << endl
         << "  mosaic:" << endl
         << "      --mosaic-out FILE     write the genome as a mosaic of the -z/-g" << endl
         << "                            haplotypes to FILE (implies --phased)" << endl
         << "      --mosaic-break-unexplained" << endl
         << "                            end a mosaic segment where no haplotype explains" << endl
         << "                            a stretch, instead of continuing the flanking one" << endl
         << "      --no-mosaic-nested    ignore nested sites, following each enclosing" << endl
         << "                            site's haplotype through its nested snarls" << endl
         << "      --no-mosaic-patch-gaps" << endl
         << "                            leave a gap where no haplotype can be followed," << endl
         << "                            instead of filling it with the reference" << endl
         << "  quality reporting (does not change genotypes):" << endl
         << "      --no-share-quality    do not scale GQ by the fraction of reads the" << endl
         << "                            called genotype explains (GQI is never scaled)" << endl
         << "      --depth-quality A     scale GQ by exp(-A * |ln DR|) at records whose" << endl
         << "                            alleles differ in length by 50 bp or more [0]" << endl
         << "      --min-confidence X    set FILTER to lowconf on records with GQN below X;" << endl
         << "                            0 disables it [0]" << endl
         << "  debugging and evaluation:" << endl
         << "      --dump-likelihoods F  write the per-site read/allele matrix to F as TSV" << endl
         << "      --flat-mixture        weight a genotype's haplotypes equally, instead of" << endl
         << "                            by the reads each is expected to produce" << endl
         << "      --regeno-ledger FILE  write one line for each site that re-genotyping" << endl
         << "                            changes" << endl
         << "      --regeno-shuffle      randomize the sign of each read's strand log-odds" << endl
         << "                            before re-genotyping" << endl
         << "      --anchors-strict-hets" << endl
         << "                            place reads at a heterozygous site by the sign of" << endl
         << "                            their strand log-odds alone, ignoring allele match" << endl
         << "GAF options:" << endl
         << "  -G, --gaf                 output GAF genotypes instead of VCF" << endl
         << "  -T, --traversals          output all candidate traversals in GAF" << endl
         << "                            without doing any genotyping" << endl
         << "  -M, --trav-padding N      extend each flank of traversals (from -T) with" << endl
         << "                            reference path by N bases if possible" << endl
         << "general options:" << endl
         << "  -v, --vcf FILE            VCF file to genotype (must have been used" << endl
         << "                            to construct input graph with -a)" << endl
         << "  -a, --genotype-snarls     genotype every snarl, including reference calls" << endl
         << "                            (use to compare multiple samples)" << endl
         << "  -A, --all-snarls          genotype every snarl, nested ones included, each" << endl
         << "                            independently of the snarl it is nested in" << endl
         << "      --max-snarl-edges N   call a snarl's children instead of the snarl when" << endl
         << "                            it has more than N edges; 0 for no limit" << endl
         << "                            [0 with --read-likelihood and -z/-g, else 10000]" << endl
         << "  -c, --min-length N        genotype only snarls with" << endl
         << "                            at least one traversal of length >= N" << endl
         << "  -C, --max-length N        genotype only snarls where" << endl
         << "                            all traversals have length <= N" << endl
         << "  -f, --ref-fasta FILE      reference FASTA" << endl
         << "                            (required if VCF has symbolic deletions/inversions)" << endl
         << "  -i, --ins-fasta FILE      insertions (required if VCF has symbolic insertions)" << endl
         << "  -s, --sample NAME         sample name [" << DEFAULT_SAMPLE_NAME << "]" << endl
         << "  -r, --snarls FILE         snarls (from vg snarls) to avoid recomputing." << endl
         << "  -g, --gbwt FILE           only call genotypes present in given GBWT index" << endl
         << "  -z, --gbz                 only call genotypes present in GBZ index" << endl
         << "                            (applies only if input graph is GBZ; on by" << endl
         << "                            default with --read-likelihood)" << endl
         << "  -N, --translation FILE    node ID translation (from vg gbwt --translation)" << endl
         << "                            to apply to snarl names in output" << endl
         << "  -O, --gbz-translation     use the ID translation from the input GBZ to" << endl
         << "                            apply snarl names to snarl names/AT fields in output" << endl
         << "  -p, --ref-path NAME       reference path to call on (may repeat; default all)" << endl
         << "  -P, --path-prefix NAME    call on all paths with this prefix (may repeat)" << endl
         << "  -S, --ref-sample NAME     call on all paths with this sample" << endl
         << "                            (cannot use with -p or -P)" << endl
         << "  -o, --ref-offset N        offset in reference path (may repeat; 1 per path)" << endl
         << "  -l, --ref-length N        override reference length for output VCF contig" << endl
         << "  -d, --ploidy N            ploidy wherever --ploidy-bed, --ploidy-preset and" << endl
         << "                            -R give none {1, 2} [2]" << endl
         << "      --ploidy-bed FILE     the sample's ploidy by region, overriding -d and" << endl
         << "                            -R, as a BED of CHROM START END PLOIDY with the" << endl
         << "                            VCF contig name as CHROM, intervals that do not" << endl
         << "                            overlap, and ploidy 0 for no calls" << endl
         << "      --ploidy-preset NAME  the sample's ploidy by region from a BED built" << endl
         << "                            into vg for one assembly and karyotype:" << endl
         << "                            chm13v2.XY, chm13v2.XX, grch38.XY, grch38.XX," << endl
         << "                            grch37.ucsc.XY, grch37.ucsc.XX, grch37.b37.XY or" << endl
         << "                            grch37.b37.XX (cannot use with --ploidy-bed)" << endl
         << "  -R, --ploidy-regex RULES  comma-separated REGEX:PLOIDY rules, each giving" << endl
         << "                            the ploidy of every reference path whose whole" << endl
         << "                            name REGEX matches. The first matching rule wins," << endl
         << "                            and unmatched paths get the ploidy from -d." << endl
         << "      --top-down            genotype nested snarls after their parents, each" << endl
         << "                            child's candidate alleles taken from its parent's" << endl
         << "                            genotype" << endl
         << "      --bottom-up           genotype nested snarls before their parents" << endl
         << "  -I, --chains              call chains instead of snarls (experimental)" << endl
         << "  -L, --cluster F           merge called alt alleles whose length-weighted" << endl
         << "                            similarity is >= F, so 1/2 of two effectively" << endl
         << "                            identical alleles becomes 1/1 [1.0; experimental]" << endl
         << "      --cluster-min-len N   apply -L only at sites whose longest allele, less" << endl
         << "                            the prefix and suffix all alleles share, is at" << endl
         << "                            least N bp; 0 for every site [50]" << endl
         << "  -Y, --star-allele         use * alleles for spanning haplotypes" << endl
         << "                            (requires --top-down)" << endl
         << "      --progress            show progress" << endl
         << "  -t, --threads N           number of threads to use" << endl
         << "  -h, --help                print this help message to stderr and exit" << endl;
}

int main_call(int argc, char** argv) {
    Logger logger("vg call");

    string pack_filename;
    string vcf_filename;
    string sample_name = DEFAULT_SAMPLE_NAME;
    string snarl_filename;
    string gbwt_filename;
    bool   gbz_paths = false;
    bool   enumerate_support = false;
    // Phasing comes from the linkage model. Where the model cannot run, an explicit --phased is
    // an error, and the default is turned off.
    bool   phased_output = true;
    string mosaic_out;
    string anchors_out;
    // Whether to genotype off-reference chains: nested chains that the reference paths do not
    // pass through. Their variants have no position on the reference, so they only get VCF records
    // when a gref fragment path (see gref.hpp) gives them a contig of their own. Selecting such a
    // path as a reference turns this on, as does --anchors-out. VG_CALL_OFF_REF_NESTING turns it
    // on without either, for testing.
    bool off_ref_nesting = getenv("VG_CALL_OFF_REF_NESTING") != nullptr;
    AnchorParams anchor_params;
    /// Counts of how the anchor code placed reads in this run.
    AnchorCounters anchor_run_counters;
    // Mosaic settings. Where no haplotype can be followed, fill the gap with the reference.
    bool mosaic_patch_gaps = true;
    // Whether a strand can switch haplotypes inside nested snarls. If not, it follows its
    // enclosing site's haplotype through them.
    bool mosaic_keep_nested = true;
    // Where no haplotype explains a stretch, continue the flanking haplotype through it.
    bool mosaic_connect_unexplained = true;
    string translation_file_name;
    bool   gbz_translation = false;
    string ref_fasta_filename;
    string ins_fasta_filename;
    vector<string> ref_paths;
    vector<string> ref_path_prefixes;
    string ref_sample;
    vector<size_t> ref_path_offsets;
    vector<size_t> ref_path_lengths;
    string min_support_string;
    string baseline_error_string;
    string bias_string;
    // require at least some support for all breakpoint edges
    // inceases sv precision, but at some recall cost.
    // think this is worth leaving on by default and not adding an option (famouse last words)
    bool expect_bp_edges = true;
    bool ratio_caller = false;
    bool legacy = false;
    int ploidy = 2;
    // copied over from vg sim
    std::vector<std::pair<std::regex, size_t>> ploidy_rules;
    string ploidy_bed_filename;
    optional<string> ploidy_preset;

    bool traversals_only = false;
    bool gaf_output = false;
    size_t trav_padding = 0;
    bool genotype_snarls = false;
    bool top_down = false;
    bool bottom_up = false;
    // Both on by default, and turned off below where they cannot apply.
    bool atomize_blocks = true;
    bool nested_calling = true;
    bool call_chains = false;
    bool all_snarls = false;
    size_t min_allele_len = 0;
    size_t max_allele_len = numeric_limits<size_t>::max();
    bool show_progress = false;

    // Nested calling option (for use with --top-down)
    bool star_allele = false;

    // -L: after genotyping, merge called ALT alleles whose similarity is at least this.
    double cluster_threshold = 1.0;
    // --cluster-min-len: the shortest site at which -L merges alleles. The similarity counts
    // sequence that all alleles share, so at a short site almost any two alleles look alike. The
    // default is the usual size of a structural variant.
    int64_t cluster_min_allele_len = 50;
    bool cluster_min_len_set = false;

    // Options for the read-likelihood genotyper (--read-likelihood).
    bool read_likelihood = false;
    ReadInputOptions reads;
    string dump_likelihoods_filename;
    bool no_mismap_term = false;
    bool no_share_quality = false;
    double depth_quality = 0.0;
    // Gap penalties for scoring reads against alleles.
    int gap_open = default_gap_open;
    int gap_extend = default_gap_extension;
    // Phase heterozygous sites with the reads, as well as with the haplotypes.
    bool read_phasing = false;
    bool regenotype = false;
    RegenotypeParams regenotype_params;
    // The number of linkage passes under --regenotype. The first pass chooses the genotypes. Each
    // later pass applies a correction from the read phasing and chooses them again.
    size_t regenotype_passes = 2;
    string regenotype_ledger;
    // Off-panel detection, and where to write its candidates.
    bool off_panel = true;
    string off_panel_dump;
    EditCallBar off_panel_bar;
    bool off_panel_svs = false;
    ReadPhasingParams read_phasing_params;
    // The preset, which sets the options the command line does not give.
    string preset;
    double insertion_gap_nats = 0.0;
    bool optimal_pairing = false;
    double min_confidence = 0.0;
    double linkage_weight = 2.0;
    double linkage_scale = 10000.0;
    double linkage_freq_prior = 5.0;
    double linkage_mutation = LinkageModel::Params().mutation;
    double hp_prior = 0.0;
    int hp_prior_run = 11;
    bool flat_mixture = false;
    double depth_weight = 0.1;
    bool depth_count_raw = false;
    double max_mismap_prob = 0.95;
    double min_mismap_prob = 0.02;
    // 0 keeps reads from aligners that give every read MAPQ 0.
    int read_min_mapq = 0;

    // constants
    const size_t avg_trav_threshold = 50;
    const size_t avg_node_threshold = 50;
    const size_t min_depth_bin_width = 50;
    const size_t max_depth_bin_width = 50000000;
    const double depth_scale_fac = 1.5;
    // Set after parsing, since it depends on -T.
    size_t max_yens_traversals = 50;
    // If not given, set after parsing, since the default depends on the traversal finder.
    size_t max_snarl_edges_opt = 0;
    // used to merge up snarls from chains when generating traversals
    const size_t max_chain_edges = 1000; 
    const size_t max_chain_trivial_travs = 5;
    constexpr int OPT_PROGRESS = 1000;
    constexpr int OPT_CLUSTER_MIN_LEN = 1002;
    constexpr int OPT_CLUSTER_POST = 1003;
    constexpr int OPT_LEGACY = 1004;
    constexpr int OPT_BOTTOM_UP = 1005;
    constexpr int OPT_ATOMIZE_BLOCKS = 1047;
    constexpr int OPT_NO_ATOMIZE_BLOCKS = 1048;
    constexpr int OPT_TOP_DOWN = 1006;
    constexpr int OPT_READ_LIKELIHOOD = 1007;
    constexpr int OPT_GAM = 1008;
    constexpr int OPT_GAF = 1009;
    constexpr int OPT_DUMP_LIKELIHOODS = 1010;
    constexpr int OPT_NO_MISMAP_TERM = 1011;
    constexpr int OPT_READ_MIN_MAPQ = 1012;
    constexpr int OPT_GAM_INDEX = 1013;
    constexpr int OPT_READ_WINDOW = 1014;
    constexpr int OPT_GAF_BASE = 1015;
    constexpr int OPT_GAF_BASE_BINARY = 1017;
    constexpr int OPT_GAF_INDEX = 1112;
    constexpr int OPT_MISMAP_MAX = 1019;
    constexpr int OPT_MISMAP_MIN = 1020;
    constexpr int OPT_INSERTION_GAP_NATS = 1091;
    constexpr int OPT_OPTIMAL_PAIRING = 1092;
    constexpr int OPT_ANCHORS_HOM_SPLIT = 1094;
    constexpr int OPT_ANCHORS_PHASE_HETS = 1102;
    constexpr int OPT_NO_ANCHORS_PHASE_HETS = 1103;
    constexpr int OPT_ANCHORS_STRICT_HETS = 1104;
    constexpr int OPT_ANCHORS_PHASE_MIN = 1095;
    constexpr int OPT_ANCHORS_PHASE_MIN_SIDE = 1096;
    constexpr int OPT_NO_OPTIMAL_PAIRING = 1093;
    constexpr int OPT_NO_SHARE_QUALITY = 1021;
    constexpr int OPT_FLAT_MIXTURE = 1023;
    constexpr int OPT_MAX_SNARL_EDGES = 1109;
    constexpr int OPT_HP_PRIOR = 1110;
    constexpr int OPT_HP_PRIOR_RUN = 1111;
    constexpr int OPT_LINKAGE_MUTATION = 1113;
    constexpr int OPT_OFF_PANEL = 1114;
    constexpr int OPT_OFF_PANEL_DUMP = 1115;
    constexpr int OPT_OFF_PANEL_BAR = 1116;
    constexpr int OPT_OFF_PANEL_INDELS = 1117;
    constexpr int OPT_OFF_PANEL_SVS = 1118;
    constexpr int OPT_NO_OFF_PANEL = 1119;
    constexpr int OPT_DEPTH_TERM = 1025;
    constexpr int OPT_DEPTH_COUNT_RAW = 1026;
    constexpr int OPT_DEPTH_QUALITY = 1027;
    constexpr int OPT_PRESET = 1072;
    constexpr int OPT_GAP_OPEN = 1070;
    constexpr int OPT_GAP_EXTEND = 1071;
    constexpr int OPT_READ_PHASING = 1074;
    constexpr int OPT_NO_READ_PHASING = 1081;
    constexpr int OPT_PHASE_MIN_Q = 1075;
    constexpr int OPT_PHASE_BREAK = 1076;
    constexpr int OPT_PHASE_RELINK = 1077;
    constexpr int OPT_PHASE_HANG = 1078;
    constexpr int OPT_PHASE_PRIOR = 1079;
    constexpr int OPT_PHASE_CAP = 1080;
    constexpr int OPT_PHASE_COHERENCE = 1107;
    constexpr int OPT_PHASE_COH_ROUNDS = 1108;
    constexpr int OPT_REGENOTYPE = 1082;
    constexpr int OPT_NO_REGENOTYPE = 1087;
    constexpr int OPT_REGENO_CEILING = 1088;
    constexpr int OPT_REGENO_HAPLOID = 1089;
    constexpr int OPT_NO_REGENO_HAPLOID = 1090;
    constexpr int OPT_REGENO_TEMPER = 1083;
    constexpr int OPT_REGENO_PASSES = 1084;
    constexpr int OPT_REGENO_SHUFFLE = 1085;
    constexpr int OPT_REGENO_LEDGER = 1086;
    constexpr int OPT_MIN_CONFIDENCE = 1042;
    constexpr int OPT_PLOIDY_BED = 1043;
    constexpr int OPT_PLOIDY_PRESET = 1120;
    constexpr int OPT_NESTED = 1044;
    constexpr int OPT_NO_NESTED = 1045;
    constexpr int OPT_NO_OFF_REF_NESTING = 1101;
    constexpr int OPT_NO_PHASED = 1046;
    constexpr int OPT_LINKAGE_WEIGHT = 1028;
    constexpr int OPT_LINKAGE_SCALE = 1030;
    constexpr int OPT_LINKAGE_FREQ_PRIOR = 1031;
    constexpr int OPT_ENUMERATE_SUPPORT = 1032;
    constexpr int OPT_PHASED = 1033;
    constexpr int OPT_MOSAIC_OUT = 1034;
    constexpr int OPT_ANCHORS_OUT = 1060;
    constexpr int OPT_ANCHORS_HET_ONLY = 1061;
    constexpr int OPT_ANCHORS_LEAF_ONLY = 1062;
    constexpr int OPT_ANCHORS_MIN_READS = 1063;
    constexpr int OPT_ANCHORS_MIN_GQN = 1064;
    constexpr int OPT_ANCHORS_MIN_SCORE = 1065;
    constexpr int OPT_ANCHORS_KEEP_OFF_CALL = 1066;
    constexpr int OPT_ANCHORS_END_PIN_MIN_NEW = 1067;
    constexpr int OPT_NO_MOSAIC_PATCH = 1053;
    constexpr int OPT_MOSAIC_PATCH = 1054;
    constexpr int OPT_NO_MOSAIC_NESTED = 1055;
    constexpr int OPT_MOSAIC_BREAK_UNEXPL = 1056;
    optind = 2; // force optind past command positional argument
    // The long options, grouped by the subsystem they configure, by its name. Every subsystem with
    // a switch is part of the read-likelihood genotyper.
    map<string, OptionGroup> long_option_groups = {
        {"vg call", {nullptr, nullptr,
         [] { return true; }, {
            {"pack", required_argument, nullptr, 'k'},
            {"bias-mode", no_argument, nullptr, 'B'},
            {"baseline-error", required_argument, nullptr, 'e'},
            {"het-bias", required_argument, nullptr, 'b'},
            {"min-support", required_argument, nullptr, 'm'},
            {"vcf", required_argument, nullptr, 'v'},
            {"genotype-snarls", no_argument, nullptr, 'a'},
            {"all-snarls", no_argument, nullptr, 'A'},
            {"min-length", required_argument, nullptr, 'c'},
            {"max-length", required_argument, nullptr, 'C'},
            {"ref-fasta", required_argument, nullptr, 'f'},
            {"ins-fasta", required_argument, nullptr, 'i'},
            {"sample", required_argument, nullptr, 's'},
            {"snarls", required_argument, nullptr, 'r'},
            {"gbwt", required_argument, nullptr, 'g'},
            {"gbz", no_argument, nullptr, 'z'},
            {"translation", required_argument, nullptr, 'N'},
            {"gbz-translation", no_argument, nullptr, 'O'},
            {"ref-path", required_argument, nullptr, 'p'},
            {"path-prefix", required_argument, nullptr, 'P'},
            {"ref-sample", required_argument, nullptr, 'S'},
            {"ref-offset", required_argument, nullptr, 'o'},
            {"ref-length", required_argument, nullptr, 'l'},
            {"ploidy", required_argument, nullptr, 'd'},
            {"ploidy-regex", required_argument, nullptr, 'R'},
            {"ploidy-bed", required_argument, nullptr, OPT_PLOIDY_BED},
            {"ploidy-preset", required_argument, nullptr, OPT_PLOIDY_PRESET},
            {"gaf", no_argument, nullptr, 'G'},
            {"traversals", no_argument, nullptr, 'T'},
            {"trav-padding", required_argument, nullptr, 'M'},
            {"legacy", no_argument, nullptr, OPT_LEGACY},
            {"top-down", no_argument, nullptr, OPT_TOP_DOWN},
            {"bottom-up", no_argument, nullptr, OPT_BOTTOM_UP},
            {"read-likelihood", no_argument, nullptr, OPT_READ_LIKELIHOOD},
            {"max-snarl-edges", required_argument, nullptr, OPT_MAX_SNARL_EDGES},
            {"chains", no_argument, nullptr, 'I'},
            {"cluster", required_argument, nullptr, 'L'},
            {"cluster-min-len", required_argument, nullptr, OPT_CLUSTER_MIN_LEN},
            // Accepted and ignored, and left out of the help, so that command lines that still
            // pass it keep working.
            {"cluster-post", no_argument, nullptr, OPT_CLUSTER_POST},
            {"star-allele", no_argument, nullptr, 'Y'},
            {"threads", required_argument, nullptr, 't'},
            {"progress", no_argument, nullptr, OPT_PROGRESS},
            {"help", no_argument, nullptr, 'h'},
        }}},
        {"the read-likelihood genotyper", {"--read-likelihood", nullptr,
         [&] { return read_likelihood; }, {
            {"no-phased", no_argument, nullptr, OPT_NO_PHASED},
            {"gam", required_argument, nullptr, OPT_GAM},
            {"gaf-reads", required_argument, nullptr, OPT_GAF},
            {"dump-likelihoods", required_argument, nullptr, OPT_DUMP_LIKELIHOODS},
            {"no-mismap-term", no_argument, nullptr, OPT_NO_MISMAP_TERM},
            {"mismap-max", required_argument, nullptr, OPT_MISMAP_MAX},
            {"mismap-min", required_argument, nullptr, OPT_MISMAP_MIN},
            {"insertion-nats", required_argument, nullptr, OPT_INSERTION_GAP_NATS},
            {"optimal-pairing", no_argument, nullptr, OPT_OPTIMAL_PAIRING},
            {"no-optimal-pairing", no_argument, nullptr, OPT_NO_OPTIMAL_PAIRING},
            {"no-share-quality", no_argument, nullptr, OPT_NO_SHARE_QUALITY},
            {"flat-mixture", no_argument, nullptr, OPT_FLAT_MIXTURE},
            {"depth-term", required_argument, nullptr, OPT_DEPTH_TERM},
            {"depth-count-raw", no_argument, nullptr, OPT_DEPTH_COUNT_RAW},
            {"depth-quality", required_argument, nullptr, OPT_DEPTH_QUALITY},
            {"preset", required_argument, nullptr, OPT_PRESET},
            {"gap-open", required_argument, nullptr, OPT_GAP_OPEN},
            {"gap-extend", required_argument, nullptr, OPT_GAP_EXTEND},
            {"read-phasing", no_argument, nullptr, OPT_READ_PHASING},
            {"no-read-phasing", no_argument, nullptr, OPT_NO_READ_PHASING},
            {"phase-min-q", required_argument, nullptr, OPT_PHASE_MIN_Q},
            {"phase-break", required_argument, nullptr, OPT_PHASE_BREAK},
            {"phase-relink", required_argument, nullptr, OPT_PHASE_RELINK},
            {"phase-hang", required_argument, nullptr, OPT_PHASE_HANG},
            {"phase-prior", required_argument, nullptr, OPT_PHASE_PRIOR},
            {"phase-cap", required_argument, nullptr, OPT_PHASE_CAP},
            {"phase-coherence", required_argument, nullptr, OPT_PHASE_COHERENCE},
            {"phase-coh-rounds", required_argument, nullptr, OPT_PHASE_COH_ROUNDS},
            {"regenotype", no_argument, nullptr, OPT_REGENOTYPE},
            {"no-regenotype", no_argument, nullptr, OPT_NO_REGENOTYPE},
            {"min-confidence", required_argument, nullptr, OPT_MIN_CONFIDENCE},
            {"linkage-weight", required_argument, nullptr, OPT_LINKAGE_WEIGHT},
            {"linkage-scale", required_argument, nullptr, OPT_LINKAGE_SCALE},
            {"linkage-prior", required_argument, nullptr, OPT_LINKAGE_FREQ_PRIOR},
            {"linkage-mutation", required_argument, nullptr, OPT_LINKAGE_MUTATION},
            {"off-panel", no_argument, nullptr, OPT_OFF_PANEL},
            {"no-off-panel", no_argument, nullptr, OPT_NO_OFF_PANEL},
            {"hp-prior", required_argument, nullptr, OPT_HP_PRIOR},
            {"hp-prior-run", required_argument, nullptr, OPT_HP_PRIOR_RUN},
            {"enumerate-support", no_argument, nullptr, OPT_ENUMERATE_SUPPORT},
            {"nested", no_argument, nullptr, OPT_NESTED},
            {"no-nested", no_argument, nullptr, OPT_NO_NESTED},
            {"atomize-blocks", no_argument, nullptr, OPT_ATOMIZE_BLOCKS},
            {"no-atomize-blocks", no_argument, nullptr, OPT_NO_ATOMIZE_BLOCKS},
            {"phased", no_argument, nullptr, OPT_PHASED},
            {"mosaic-out", required_argument, nullptr, OPT_MOSAIC_OUT},
            {"anchors-out", required_argument, nullptr, OPT_ANCHORS_OUT},
            {"read-min-mapq", required_argument, nullptr, OPT_READ_MIN_MAPQ},
            {"gam-index", required_argument, nullptr, OPT_GAM_INDEX},
            {"gaf-index", required_argument, nullptr, OPT_GAF_INDEX},
            {"gaf-base", required_argument, nullptr, OPT_GAF_BASE},
            {"gaf-base-binary", required_argument, nullptr, OPT_GAF_BASE_BINARY},
            {"read-window", required_argument, nullptr, OPT_READ_WINDOW},
        }}},
        {"off-panel variant calling", {"--off-panel", "the read-likelihood genotyper",
         [&] { return off_panel; }, {
            {"off-panel-dump", required_argument, nullptr, OPT_OFF_PANEL_DUMP},
            {"off-panel-bar", required_argument, nullptr, OPT_OFF_PANEL_BAR},
            {"off-panel-indels", required_argument, nullptr, OPT_OFF_PANEL_INDELS},
            {"off-panel-svs", no_argument, nullptr, OPT_OFF_PANEL_SVS},
        }}},
        {"the anchor file", {"--anchors-out", "the read-likelihood genotyper",
         [&] { return !anchors_out.empty(); }, {
            {"no-off-ref-nesting", no_argument, nullptr, OPT_NO_OFF_REF_NESTING},
            {"anchors-hom-split", no_argument, nullptr, OPT_ANCHORS_HOM_SPLIT},
            {"anchors-phase-hets", no_argument, nullptr, OPT_ANCHORS_PHASE_HETS},
            {"no-anchors-phase-hets", no_argument, nullptr, OPT_NO_ANCHORS_PHASE_HETS},
            {"anchors-strict-hets", no_argument, nullptr, OPT_ANCHORS_STRICT_HETS},
            {"split-min-q", required_argument, nullptr, OPT_ANCHORS_PHASE_MIN},
            {"split-min-side", required_argument, nullptr, OPT_ANCHORS_PHASE_MIN_SIDE},
            {"anchors-het-only", no_argument, nullptr, OPT_ANCHORS_HET_ONLY},
            {"anchors-leaf-only", no_argument, nullptr, OPT_ANCHORS_LEAF_ONLY},
            {"anchors-reads", required_argument, nullptr, OPT_ANCHORS_MIN_READS},
            {"anchors-min-gqn", required_argument, nullptr, OPT_ANCHORS_MIN_GQN},
            {"anchors-min-q", required_argument, nullptr, OPT_ANCHORS_MIN_SCORE},
            {"anchors-keep-off-call", no_argument, nullptr, OPT_ANCHORS_KEEP_OFF_CALL},
            {"anchors-end-new", required_argument, nullptr, OPT_ANCHORS_END_PIN_MIN_NEW},
        }}},
        {"the mosaic", {"--mosaic-out", "the read-likelihood genotyper",
         [&] { return !mosaic_out.empty(); }, {
            {"mosaic-patch-gaps", no_argument, nullptr, OPT_MOSAIC_PATCH},
            {"no-mosaic-patch-gaps", no_argument, nullptr, OPT_NO_MOSAIC_PATCH},
            {"no-mosaic-nested", no_argument, nullptr, OPT_NO_MOSAIC_NESTED},
            {"mosaic-break-unexplained", no_argument, nullptr, OPT_MOSAIC_BREAK_UNEXPL},
        }}},
        // Re-genotyping needs read phasing, so it is off without it, even where a preset turned it
        // on.
        {"re-genotyping", {"--regenotype", "the read-likelihood genotyper",
         [&] { return regenotype && read_phasing; }, {
            {"regeno-ceiling", required_argument, nullptr, OPT_REGENO_CEILING},
            {"regeno-haploid", no_argument, nullptr, OPT_REGENO_HAPLOID},
            {"no-regeno-haploid", no_argument, nullptr, OPT_NO_REGENO_HAPLOID},
            {"regeno-temper", required_argument, nullptr, OPT_REGENO_TEMPER},
            {"regeno-passes", required_argument, nullptr, OPT_REGENO_PASSES},
            {"regeno-shuffle", no_argument, nullptr, OPT_REGENO_SHUFFLE},
            {"regeno-ledger", required_argument, nullptr, OPT_REGENO_LEDGER},
        }}},
    };

    // getopt_long takes a single array ending in an all-zero entry. Each option has a unique
    // option number, its `val`, which getopt_long returns for it. We record each option's name and
    // subsystem by its number, and its number by its name.
    vector<struct option> long_options;
    unordered_map<int, string> option_name;
    unordered_map<int, string> option_subsystem;
    unordered_map<string, int> option_number;
    for (const auto& [subsystem, group] : long_option_groups) {
        for (const struct option& o : group.options) {
            long_options.emplace_back(o);
            option_name[o.val] = o.name;
            option_subsystem[o.val] = subsystem;
            option_number[o.name] = o.val;
        }
    }
    long_options.push_back({nullptr, 0, nullptr, 0});

    const char* short_options = "k:Be:b:m:v:aAc:C:f:i:s:r:g:zN:Op:P:S:o:l:d:R:GTM:IL:Yt:h?";

    // Set the option with option number `c`, with its argument `optarg`, or null for a switch. The
    // argument has the name of getopt's, so that the cases read as in other subcommands.
    auto apply_option = [&](int c, const char* optarg) {
        switch (c)
        {
        case 'k':
            pack_filename = require_exists(logger, optarg);
            break;
        case 'B':
            ratio_caller = true;
            break;
        case 'b':
            bias_string = optarg;
            break;
        case 'm':
            min_support_string = optarg;
            break;
        case 'e':
            baseline_error_string = optarg;
            break;            
        case 'v':
            vcf_filename = require_exists(logger, optarg);
            break;
        case 'a':
            genotype_snarls = true;
            break;
        case 'A':
            all_snarls = true;
            break;
        case 'c':
            min_allele_len = parse<size_t>(optarg);
            break;
        case 'C':
            max_allele_len = parse<size_t>(optarg);
            break;
        case 'f':
            ref_fasta_filename = require_exists(logger, optarg);
            break;
        case 'i':
            ins_fasta_filename = require_exists(logger, optarg);
            break;
        case 's':
            sample_name = optarg;
            break;
        case 'r':
            snarl_filename = require_exists(logger, optarg);
            break;
        case 'g':
            gbwt_filename = require_exists(logger, optarg);
            break;
        case 'z':
            gbz_paths = true;
            break;
        case OPT_ENUMERATE_SUPPORT:
            enumerate_support = true;
            break;
        case OPT_PHASED:
            phased_output = true;
            break;
        case OPT_MOSAIC_BREAK_UNEXPL:
            mosaic_connect_unexplained = false;
            break;
        case OPT_NO_MOSAIC_NESTED:
            mosaic_keep_nested = false;
            break;
        case OPT_MOSAIC_PATCH:
            mosaic_patch_gaps = true;
            break;
        case OPT_NO_MOSAIC_PATCH:
            mosaic_patch_gaps = false;
            break;
        case OPT_ANCHORS_OUT:
            anchors_out = optarg;
            anchor_params.enabled = true;
            break;
        case OPT_ANCHORS_HET_ONLY:
            anchor_params.het_only = true;
            break;
        case OPT_ANCHORS_LEAF_ONLY:
            anchor_params.leaf_only = true;
            break;
        case OPT_ANCHORS_MIN_READS:
            anchor_params.min_reads = parse<size_t>(optarg);
            break;
        case OPT_ANCHORS_MIN_GQN:
            anchor_params.min_gqn = parse<double>(optarg);
            break;
        case OPT_ANCHORS_MIN_SCORE:
            anchor_params.min_read_score = parse<double>(optarg);
            break;
        case OPT_ANCHORS_END_PIN_MIN_NEW:
            anchor_params.end_pin_min_new = parse<size_t>(optarg);
            break;
        case OPT_ANCHORS_KEEP_OFF_CALL:
            anchor_params.keep_off_call = true;
            break;
        case OPT_MOSAIC_OUT:
            mosaic_out = optarg;
            break;
        case 'N':
            translation_file_name = require_exists(logger, optarg);
            break;
        case 'O':
            gbz_translation = true;
            break;            
        case 'p':
            ref_paths.push_back(optarg);
            break;
        case 'P':
            ref_path_prefixes.push_back(optarg);
            break;
        case 'S':
            ref_sample = optarg;
            break;            
        case 'o':
            ref_path_offsets.push_back(parse<int>(optarg));
            break;
        case 'l':
            ref_path_lengths.push_back(parse<int>(optarg));
            break;
        case 'd':
            ploidy = parse<int>(optarg);
            break;
        case 'R':
            for (auto& rule : split_delims(optarg, ",")) {
                // For each comma-separated rule
                auto parts = split_delims(rule, ":");
                if (parts.size() != 2) {
                    logger.error() << "ploidy rules must be REGEX:PLOIDY" << endl;
                }
                try {
                    // Parse the regex
                    std::regex match(parts[0]);
                    size_t weight = parse<size_t>(parts[1]);
                    // Save the rule.  The {1,2} restriction the callers impose is checked where the
                    // rule is applied, not here: a rule that matches none of the called contigs
                    // never reaches a caller, and rejecting it would break command lines that work.
                    ploidy_rules.emplace_back(match, weight);
                } catch (const std::regex_error& e) {
                    // This is not a good regex
                    logger.error() << "unacceptable regular expression \""
                                   << parts[0] << "\": " << e.what() << endl;
                }
            }
            break;            
        case 'G':
            gaf_output = true;
            break;
        case 'T':
            traversals_only = true;
            gaf_output = true;
            break;
        case 'M':
            trav_padding = parse<size_t>(optarg);
            break;
        case 'L':
            cluster_threshold = parse<double>(optarg);
            break;
        case OPT_TOP_DOWN:
            top_down = true;
            break;
        case OPT_ATOMIZE_BLOCKS:
            atomize_blocks = true;
            break;
        case OPT_NO_ATOMIZE_BLOCKS:
            atomize_blocks = false;
            break;
        case OPT_BOTTOM_UP:
            bottom_up = true;
            break;
        case OPT_READ_LIKELIHOOD:
            read_likelihood = true;
            break;
        case OPT_GAM:
            reads.gam = require_exists(logger, optarg);
            break;
        case OPT_GAF:
            reads.gaf = require_exists(logger, optarg);
            break;
        case OPT_DUMP_LIKELIHOODS:
            dump_likelihoods_filename = ensure_writable(logger, optarg);
            break;
        case OPT_NO_MISMAP_TERM:
            no_mismap_term = true;
            break;
        case OPT_FLAT_MIXTURE:
            flat_mixture = true;
            break;
        case OPT_MAX_SNARL_EDGES:
            max_snarl_edges_opt = parse<size_t>(optarg);
            break;
        case OPT_DEPTH_TERM:
            depth_weight = parse<double>(optarg);
            break;
        case OPT_DEPTH_COUNT_RAW:
            depth_count_raw = true;
            break;
        case OPT_PRESET:
            preset = optarg;
            break;
        case OPT_GAP_OPEN:
            gap_open = parse<int>(optarg);
            if (gap_open < 1 || gap_open > 127) {
                logger.error() << "--gap-open must be between 1 and 127" << endl;
            }
            break;
        case OPT_READ_PHASING:
            read_phasing = true;
            break;
        case OPT_NO_READ_PHASING:
            read_phasing = false;
            break;
        case OPT_PHASE_MIN_Q:
            read_phasing_params.reliability = parse<double>(optarg);
            break;
        case OPT_PHASE_BREAK:
            read_phasing_params.break_threshold = parse<double>(optarg);
            break;
        case OPT_PHASE_RELINK:
            read_phasing_params.relink = parse<size_t>(optarg);
            break;
        case OPT_PHASE_HANG:
            read_phasing_params.hang = parse<size_t>(optarg);
            break;
        case OPT_PHASE_PRIOR:
            read_phasing_params.panel_weight = parse<double>(optarg);
            break;
        case OPT_PHASE_CAP:
            read_phasing_params.cap = parse<double>(optarg);
            break;
        case OPT_PHASE_COHERENCE:
            read_phasing_params.coherence_min = parse<double>(optarg);
            if (read_phasing_params.coherence_min < 0.0
                || read_phasing_params.coherence_min > 1.0) {
                logger.error() << "--phase-coherence is a fraction in [0,1]" << endl;
            }
            break;
        case OPT_PHASE_COH_ROUNDS:
            read_phasing_params.coherence_rounds = parse<size_t>(optarg);
            if (read_phasing_params.coherence_rounds < 1) {
                logger.error() << "--phase-coh-rounds must be >= 1" << endl;
            }
            break;
        case OPT_REGENOTYPE:
            regenotype = true;
            break;
        case OPT_NO_REGENOTYPE:
            regenotype = false;
            break;
        case OPT_REGENO_TEMPER:
            regenotype_params.temper = parse<double>(optarg);
            if (regenotype_params.temper < 0) {
                logger.error() << "--regeno-temper must be >= 0 (0 leaves the genotypes unchanged)"
                               << endl;
            }
            break;
        case OPT_REGENO_PASSES:
            regenotype_passes = parse<size_t>(optarg);
            if (regenotype_passes < 1 || regenotype_passes > 20) {
                logger.error() << "--regeno-passes must be between 1 and 20" << endl;
            }
            break;
        case OPT_REGENO_SHUFFLE:
            regenotype_params.shuffle = true;
            break;
        case OPT_REGENO_CEILING:
            regenotype_params.ceiling = parse<double>(optarg);
            if (regenotype_params.ceiling <= 0.0 || regenotype_params.ceiling > 1.0) {
                logger.error() << "--regeno-ceiling must be in (0, 1]; 1 means no limit" << endl;
            }
            break;
        case OPT_REGENO_HAPLOID:
            regenotype_params.haploid_include = true;
            break;
        case OPT_NO_REGENO_HAPLOID:
            regenotype_params.haploid_include = false;
            break;
        case OPT_REGENO_LEDGER:
            regenotype_ledger = optarg;
            break;
        case OPT_OFF_PANEL:
            off_panel = true;
            break;
        case OPT_NO_OFF_PANEL:
            off_panel = false;
            break;
        case OPT_OFF_PANEL_DUMP:
            off_panel_dump = optarg;
            break;
        case OPT_OFF_PANEL_SVS:
            off_panel_svs = true;
            break;
        case OPT_OFF_PANEL_BAR: {
            vector<string> parts = split_delims(optarg, ",");
            if (parts.size() != 3 && parts.size() != 4) {
                logger.error() << "--off-panel-bar takes K,F,G[,D]" << endl;
            }
            off_panel_bar.min_alt_fragments = parse<size_t>(parts[0]);
            off_panel_bar.min_fraction = parse<double>(parts[1]);
            off_panel_bar.min_gq = parse<double>(parts[2]);
            if (parts.size() == 4) {
                off_panel_bar.max_depth_ratio = parse<double>(parts[3]);
            }
            break;
        }
        case OPT_OFF_PANEL_INDELS: {
            vector<string> parts = split_delims(optarg, ",");
            if (parts.size() != 3) {
                logger.error() << "--off-panel-indels takes K,F,H" << endl;
            }
            off_panel_bar.indel_min_alt_fragments = parse<size_t>(parts[0]);
            off_panel_bar.indel_min_fraction = parse<double>(parts[1]);
            off_panel_bar.indel_max_homopolymer = parse<uint32_t>(parts[2]);
            break;
        }
        case OPT_GAP_EXTEND:
            gap_extend = parse<int>(optarg);
            if (gap_extend < 1 || gap_extend > 127) {
                logger.error() << "--gap-extend must be between 1 and 127" << endl;
            }
            break;
        case OPT_DEPTH_QUALITY:
            depth_quality = parse<double>(optarg);
            break;
        case OPT_MIN_CONFIDENCE:
            min_confidence = parse<double>(optarg);
            break;
        case OPT_PLOIDY_BED:
            ploidy_bed_filename = require_exists(logger, optarg);
            break;
        case OPT_PLOIDY_PRESET:
            ploidy_preset = optarg;
            break;
        case OPT_NESTED:
            nested_calling = true;
            break;
        case OPT_NO_OFF_REF_NESTING:
            off_ref_nesting = false;
            break;
        case OPT_NO_NESTED:
            nested_calling = false;

            break;
        case OPT_NO_PHASED:
            phased_output = false;
            break;
        case OPT_LINKAGE_WEIGHT:
            linkage_weight = parse<double>(optarg);
            break;
        case OPT_LINKAGE_SCALE:
            linkage_scale = parse<double>(optarg);
            break;
        case OPT_LINKAGE_FREQ_PRIOR:
            linkage_freq_prior = parse<double>(optarg);
            break;
        case OPT_LINKAGE_MUTATION:
            linkage_mutation = parse<double>(optarg);
            break;
        case OPT_HP_PRIOR:
            hp_prior = parse<double>(optarg);
            break;
        case OPT_HP_PRIOR_RUN:
            hp_prior_run = parse<int>(optarg);
            break;
        case OPT_NO_SHARE_QUALITY:
            no_share_quality = true;
            break;
        case OPT_MISMAP_MAX:
            max_mismap_prob = parse<double>(optarg);
            break;
        case OPT_MISMAP_MIN:
            min_mismap_prob = parse<double>(optarg);
            break;
        case OPT_INSERTION_GAP_NATS:
            insertion_gap_nats = parse<double>(optarg);
            break;
        case OPT_OPTIMAL_PAIRING:
            optimal_pairing = true;
            break;
        case OPT_ANCHORS_HOM_SPLIT:
            anchor_params.hom_split = true;
            break;
        case OPT_ANCHORS_PHASE_HETS:
            anchor_params.phase_hets = true;
            break;
        case OPT_NO_ANCHORS_PHASE_HETS:
            anchor_params.phase_hets = false;
            break;
        case OPT_ANCHORS_STRICT_HETS:
            anchor_params.strict_hets = true;
            break;
        case OPT_ANCHORS_PHASE_MIN:
            anchor_params.phase_min = parse<double>(optarg);
            if (anchor_params.phase_min < 0) {
                logger.error() << "--split-min-q must be at least 0" << endl;
            }
            break;
        case OPT_ANCHORS_PHASE_MIN_SIDE:
            anchor_params.phase_min_side = parse<size_t>(optarg);
            if (anchor_params.phase_min_side < 1) {
                logger.error() << "--split-min-side must be at least 1" << endl;
            }
            break;
        case OPT_NO_OPTIMAL_PAIRING:
            optimal_pairing = false;
            break;
        case OPT_READ_MIN_MAPQ:
            read_min_mapq = parse<int>(optarg);
            break;
        case OPT_GAM_INDEX:
            reads.gam_index = require_exists(logger, optarg);
            break;
        case OPT_GAF_INDEX:
            reads.gaf_index = require_exists(logger, optarg);
            break;
        case OPT_GAF_BASE:
            reads.gaf_base = require_exists(logger, optarg);
            break;
        case OPT_GAF_BASE_BINARY:
            reads.gaf_base_binary = optarg;
            break;
        case OPT_READ_WINDOW:
            reads.window_size = parse<size_t>(optarg);
            break;
        case 'I':
            call_chains = true;
            break;
        case OPT_CLUSTER_POST:
            logger.warn() << "--cluster-post is deprecated and ignored: -L/--cluster now always "
                          << "merges after genotyping" << endl;
            break;
        case OPT_CLUSTER_MIN_LEN:
            cluster_min_allele_len = parse<int64_t>(optarg);
            cluster_min_len_set = true;
            if (cluster_min_allele_len < 0) {
                logger.error() << "--cluster-min-len must be >= 0" << endl;
            }
            break;
        case 'Y':
            star_allele = true;
            break;
        case OPT_LEGACY:
            legacy = true;
            break;
        case OPT_PROGRESS:
            show_progress = true;
            break;
        case 't':
            set_thread_count(logger, optarg);
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_call(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    };

    // The option number of each option given on the command line.
    unordered_set<int> options_seen;
    while (true) {
        int option_index = 0;
        int c = getopt_long(argc, argv, short_options, long_options.data(), &option_index);
        if (c == -1) {
            break;
        }
        options_seen.insert(c);
        apply_option(c, optarg);
    }

    // Return whether the command line gave the option `name`, or its "no-" form, which sets the
    // same thing. A setting the command line did not give has its default, or the preset's value.
    auto is_explicit = [&](const string& name) {
        auto given = [&](const string& option) {
            auto number = option_number.find(option);
            return number != option_number.end() && options_seen.contains(number->second);
        };
        return given(name) || given("no-" + name);
    };

    // Apply the preset to the settings the command line did not give.
    if (!preset.empty()) {
        auto chosen = std::ranges::find_if(PRESETS, [&](const pair<string, vector<PresetSetting>>& p) {
            return p.first == preset;
        });
        if (chosen == PRESETS.end()) {
            vector<string> names;
            for (const auto& known : PRESETS) {
                names.push_back(known.first);
            }
            logger.error() << "unknown --preset '" << preset << "'; known presets: "
                           << join_with_commas(names) << endl;
        }
        for (const PresetSetting& setting : chosen->second) {
            if (!is_explicit(setting.option)) {
                apply_option(option_number.at(setting.option), setting.value);
            }
        }
    }

    if (argc <= 2) {
        help_call(argv);
        return 1;
    }

    // parse the supports (stick together to keep number of options down)
    vector<string> support_toks = split_delims(min_support_string, ",");
    double min_allele_support = -1;
    double min_site_support = -1;
    if (support_toks.size() >= 1) {
        min_allele_support = parse<double>(support_toks[0]);
        min_site_support = min_allele_support;
    }
    if (support_toks.size() == 2) {
        min_site_support = parse<double>(support_toks[1]);
    } else if (support_toks.size() > 2) {
        logger.error() << "-m option expects at most two comma separated numbers M,N" << endl;
    }
    // parse the biases
    vector<string> bias_toks = split_delims(bias_string, ",");
    double het_bias = -1;
    double ref_het_bias = -1;
    if (bias_toks.size() >= 1) {
        het_bias = parse<double>(bias_toks[0]);
        ref_het_bias = het_bias;
    }
    if (bias_toks.size() == 2) {
        ref_het_bias = parse<double>(bias_toks[1]);
    } else if (bias_toks.size() > 2) {
        logger.error() << "-b option expects at most two comma separated numbers M,N" << endl;
    }
    // parse the baseline errors (defaults are in snarl_caller.hpp)
    vector<string> error_toks = split_delims(baseline_error_string, ",");
    double baseline_error_large = -1;
    double baseline_error_small = -1;
    if (error_toks.size() == 2) {
        baseline_error_small = parse<double>(error_toks[0]);
        baseline_error_large = parse<double>(error_toks[1]);
        if (baseline_error_small > baseline_error_large) {
            logger.warn() << "with baseline error -e X,Y option, "
                          << "small variant error (X) normally less than large (Y)" << endl;
        }
    } else if (error_toks.size() != 0) {
        logger.error() << "-e option expects exactly two comma-separated numbers X,Y" << endl;
    }

    if (trav_padding > 0 && traversals_only == false) {
        logger.error() << "-M option can only be used in conjunction with -T" << endl;
    }

    // With -T the candidate traversals are themselves the output, and none of them has to be
    // genotyped, so we can afford to look for more of them.
    max_yens_traversals = traversals_only ? 100 : 50;

    if (!(linkage_mutation > 0.0 && linkage_mutation <= 0.5)) {
        cerr << "error [vg call]: --linkage-mutation takes a probability above 0 and at most 0.5"
             << endl;
        return 1;
    }
    if (hp_prior < 0.0 || hp_prior_run < 1) {
        cerr << "error [vg call]: --hp-prior takes a value >= 0, and --hp-prior-run a run of at least 1"
             << endl;
        return 1;
    }
    // --atomize-blocks is on by default, so where it cannot apply we turn it off, and only an
    // explicit --atomize-blocks is an error. These checks depend only on the options, so they run
    // before the graph is loaded.
    if (atomize_blocks && genotype_snarls) {
        // -a writes the same records, one per snarl, whatever the sample, but block emission
        // writes one record per difference in the called haplotypes.
        if (is_explicit("atomize-blocks")) {
            logger.error() << "--atomize-blocks cannot be combined with -a/--genotype-snarls: "
                           << "a block list depends on the called haplotypes, so the record set "
                           << "would stop being sample-independent" << endl;
        }
        atomize_blocks = false;
    }
    if (atomize_blocks && top_down) {
        if (is_explicit("atomize-blocks")) {
            logger.error() << "--atomize-blocks cannot be combined with --top-down" << endl;
        }
        atomize_blocks = false;
    }

    if (!vcf_filename.empty() && genotype_snarls) {
        logger.error() << "-v and -a options cannot be used together" << endl;
    }

    if ((min_allele_len > 0 || max_allele_len < numeric_limits<size_t>::max())
        && (legacy || !vcf_filename.empty() || bottom_up)) {
        logger.error() << "-c/-C not supported with -v, --legacy, or --bottom-up" << endl;
    }
    if (!ref_paths.empty() && !ref_sample.empty()) {
        logger.error() << "-S cannot be used with -p" << endl;
    }
    if (!ref_path_prefixes.empty() && !ref_sample.empty()) {
        logger.error() << "-S cannot be used with -P" << endl;
    }
    if (!ref_path_prefixes.empty() && !ref_paths.empty()) {
        logger.error() << "-P cannot be used with -p" << endl;
    }
    if (!ploidy_bed_filename.empty() && ploidy_preset.has_value()) {
        logger.error() << "--ploidy-bed cannot be used with --ploidy-preset" << endl;
    }

    // The sample's ploidy by region, if --ploidy-preset or --ploidy-bed gives it. We read it before
    // loading the graph, so that a misspelled preset name or a malformed BED fails quickly.
    optional<PloidyRegions> ploidy_regions;
    if (ploidy_preset.has_value()) {
        ploidy_regions = PloidyRegions::preset(*ploidy_preset);
    } else if (!ploidy_bed_filename.empty()) {
        ploidy_regions = PloidyRegions(ploidy_bed_filename);
    }

    // Check -L/--cluster before loading the graph, so that a mistake fails quickly. We reject a
    // threshold outside [0, 1] rather than clamping it, since "-L 5" is probably a typo for
    // "-L 0.5".
    if (cluster_threshold < 0.0 || cluster_threshold > 1.0 || std::isnan(cluster_threshold)) {
        logger.error() << "-L/--cluster threshold must be in range [0.0, 1.0]" << endl;
    }
    // Warn only when --cluster-min-len was given: its default is nonzero, so otherwise every run
    // without -L would warn.
    if (cluster_min_len_set && cluster_min_allele_len > 0 && cluster_threshold >= 1.0) {
        logger.warn() << "--cluster-min-len has no effect without -L (cluster threshold < 1.0)" << endl;
    }
    // -L merges alleles as a SiteEmitter builds each record, and neither the VCF genotyper (-v) nor
    // GAF output builds records with one.
    if (cluster_threshold < 1.0 && !vcf_filename.empty()) {
        logger.error() << "-L/--cluster cannot be used when genotyping a VCF (-v)" << endl;
    }
    if (cluster_threshold < 1.0 && (gaf_output || traversals_only)) {
        logger.error() << "-L/--cluster cannot be used with GAF output (-G/-T)" << endl;
    }
    // The ratio caller's QUAL and lowxadl filter describe a heterozygous genotype that the merge
    // would remove.
    if (cluster_threshold < 1.0 && ratio_caller) {
        logger.error() << "-L/--cluster cannot be used with the ratio caller (-B)" << endl;
    }
    // --bottom-up's NestedFlowCaller writes a child snarl into its parent's allele as a Visit to the
    // snarl rather than to a node. The merge cannot compare such alleles, so it would do nothing at
    // nested sites.
    if (cluster_threshold < 1.0 && bottom_up) {
        logger.error() << "-L/--cluster cannot be used with --bottom-up mode" << endl;
    }
    // -Y writes "*" in a child record when a deletion in the parent's allele spans the child. If the
    // merge removes that parent allele, the "*" refers to an allele that is not in the file. Without
    // -Y, a merged parent record may disagree with its children's records, which we allow: the
    // parent gives an approximate view of a large variant and the children the exact one, and
    // INFO/MAT records the merge.
    if (cluster_threshold < 1.0 && star_allele) {
        logger.error() << "-L/--cluster cannot be used with -Y/--star-allele" << endl;
    }
    // -L is not supported by LegacyCaller, which has its own traversal finder and support model.
    if (cluster_threshold < 1.0 && legacy) {
        logger.error() << "-L/--cluster cannot be used with the legacy caller (--legacy)" << endl;
    }
    // A merge needs two different called ALT alleles, which a haploid genotype cannot have. We warn
    // only for -d 1, since we cannot tell here which contigs -R/--ploidy-regex makes haploid.
    if (cluster_threshold < 1.0 && ploidy == 1) {
        logger.warn() << "-L/--cluster has no effect at ploidy 1 (-d 1)" << endl;
    }
    // The GAF writers (-G/-T) cannot write NestedFlowCaller's Visits to snarls either.
    if (bottom_up && (gaf_output || traversals_only)) {
        logger.error() << "--bottom-up cannot be used with GAF output (-G/-T)" << endl;
    }

    // Read the graph
    unique_ptr<PathHandleGraph> path_handle_graph;
    unique_ptr<GBZGraph> gbz_graph;
    gbwt::GBWT* gbwt_index = nullptr;
    PathHandleGraph* graph = nullptr;
    string graph_filename = get_input_file_name(optind, argc, argv);
    if (show_progress) logger.info() << "Loading graph " << graph_filename << endl;
    auto input = vg::io::VPKG::try_load_first<GBZGraph, PathHandleGraph>(graph_filename);
    if (show_progress) logger.info() << "Loaded graph" << endl;
    if (get<0>(input)) {        
        gbz_graph = std::move(get<0>(input));
        graph = gbz_graph.get();
        if (show_progress) logger.info() << "GBZ input detected" << endl;
        if (gbz_paths) {
            if (show_progress) logger.info() << "Restricting search to GBZ haplotypes" << endl;
            gbwt_index = &gbz_graph->gbz.index;
        } else if (!read_likelihood) {
            // Under --read-likelihood, haplotype enumeration is chosen automatically below.
            logger.info() << "You can restrict the search to GBZ haplotypes, "
                          << "often to the benefict of speed and accuracy, with the -z option" << endl;
        }
    } else if (get<1>(input)) {
        path_handle_graph = std::move(get<1>(input));
        graph = path_handle_graph.get();
    } else {
        logger.error() << "Input graph is not a GBZ or path handle graph" << endl;
    }
    if (gbz_paths && !gbz_graph) {
        logger.error() << "-z can only be used when input graph is in GBZ format" << endl;
    }
    if (gbz_translation && !gbz_graph) {
        logger.error() << "-O can only be used when input graph is in GBZ format" << endl;
    }
    
    // Read the translation
    unique_ptr<unordered_map<nid_t, pair<string, size_t>>> translation;
    if (gbz_graph.get() != nullptr && gbz_translation) {
        // try to get the translation from the graph
        translation = make_unique<unordered_map<nid_t, pair<string, size_t>>>();
        *translation = load_translation_back_map(gbz_graph->gbz.graph);
        if (translation->empty()) {
            // not worth keeping an empty translation
            translation = nullptr;
        }
    }
    if (!translation_file_name.empty()) {
        if (!translation->empty()) {
            logger.warn() << "Using translation from -N overrides that in input GBZ "
                          << "(you probably don't want to use -N)" << endl;
        }        
        ifstream translation_file(translation_file_name.c_str());
        translation = make_unique<unordered_map<nid_t, pair<string, size_t>>>();
        *translation = load_translation_back_map(*graph, translation_file);
    }    
    
    // Apply overlays as necessary
    bool need_path_positions = vcf_filename.empty();
    bool need_vectorizable = !pack_filename.empty();
    // When not using GBWT/GBZ, embedded HAPLOTYPE paths are the sample alleles
    bool embedded_haplotype_paths = gbwt_filename.empty() && !gbz_graph;
    bdsg::ReferencePathOverlayHelper pp_overlay_helper;
    bdsg::ReferencePathVectorizableOverlayHelper ppv_overlay_helper;
    bdsg::PathVectorizableOverlayHelper pv_overlay_helper;
    if (show_progress) {
        logger.info() << "Applying overlays if necessary (i.e. input not in XG format)" << endl;
    }
    if (need_path_positions && need_vectorizable) {
        graph = dynamic_cast<PathHandleGraph*>(ppv_overlay_helper.apply(graph, embedded_haplotype_paths));
    } else if (need_path_positions && !need_vectorizable) {
        graph = dynamic_cast<PathHandleGraph*>(pp_overlay_helper.apply(graph, embedded_haplotype_paths));
    } else if (!need_path_positions && need_vectorizable) {
        graph = dynamic_cast<PathHandleGraph*>(pv_overlay_helper.apply(graph));
    }
    if (show_progress) logger.info() << "Applied overlays" << endl;
    
    // Check our offsets
    if (ref_path_offsets.size() != 0 && ref_path_offsets.size() != ref_paths.size()) {
        logger.error() << "when using -o, the same number of paths must be given with -p" << endl;
    }
    if (!ref_path_offsets.empty() && !vcf_filename.empty()) {
        logger.error() << "-o cannot be used with -v" << endl;
    }
    // Check our ref lengths
    if (ref_path_lengths.size() != 0 && ref_path_lengths.size() != ref_paths.size()) {
        logger.error() << "when using -l, the same number of paths must be given with -p" << endl;
    }
    // Check bias option
    if (!bias_string.empty() && !ratio_caller) {
        logger.error() << "-b can only be used with -B" << endl;
    }
    // Check ploidy option
    if (ploidy < 1 || ploidy > 2) {
        logger.error() << "ploidy (-d) must be either 1 or 2" << endl;
    }
    if (ratio_caller == true && ploidy != 2) {
        logger.error() << "ploidy (-d) must be 2 when using ratio caller (-B)" << endl;
    }
    if (legacy == true && ploidy != 2) {
        logger.error() << "ploidy (-d) must be 2 when using legacy caller (--legacy)" << endl;
    }
    if (!vcf_filename.empty() && !gbwt_filename.empty()) {
        logger.error() << "GBWT (-g) cannot be used when genotyping VCF (-v)" << endl;
    }
    if (legacy == true && !gbwt_filename.empty()) {
        logger.error() << "GBWT (-g) cannot be used with legacy caller (--legacy)" << endl;
    }
    if (gbz_paths && !gbwt_filename.empty()) {
        logger.error() << "GBWT (-g) cannot be used with GBZ graph (-z): choose one or the other" << endl;
    }

    // An index is useless without the reads it indexes, whatever else is given.
    if (!reads.gam_index.empty() && reads.gam.empty()) {
        logger.error() << "--gam-index requires --gam" << endl;
    }
    if (!reads.gaf_index.empty() && reads.gaf.empty()) {
        logger.error() << "--gaf-index requires --gaf-reads" << endl;
    }

    // Refuse options for a subsystem that is off.
    refuse_inactive_subsystem_options(logger, long_option_groups, options_seen, option_name,
                                      option_subsystem);

    // Refuse or turn off the options that start from the phase the linkage model gives each site:
    // read phasing, re-genotyping and hom splitting. An option given explicitly is an error, and
    // one a preset turned on is turned off. `why` ends the error message, saying why there is no
    // phase.
    auto refuse_phase_dependents = [&](const string& why) {
        vector<string> offenders;
        if (read_phasing && is_explicit("read-phasing")) {
            offenders.emplace_back("--read-phasing");
        }
        if (regenotype && is_explicit("regenotype")) {
            offenders.emplace_back("--regenotype");
        }
        if (anchor_params.hom_split) {
            offenders.emplace_back("--anchors-hom-split");
        }
        if (!offenders.empty()) {
            logger.error() << join_with_commas(offenders) << (offenders.size() == 1 ? " starts" : " start")
                           << " from the phase the linkage model gives each site, " << why
                           << endl;
        }
        read_phasing = false;
        regenotype = false;
    };
    if (is_explicit("phased") && !phased_output) {
        refuse_phase_dependents("which --no-phased turns off");
    }

    // --read-likelihood needs exactly one read source, and cannot be combined with the ratio
    // support caller or with the calling modes the read-likelihood caller does not have.
    if (read_likelihood) {
        int read_source_count = (reads.gam.empty() ? 0 : 1) + (reads.gaf.empty() ? 0 : 1) +
                                (reads.gaf_base.empty() ? 0 : 1);
        if (read_source_count == 0) {
            logger.error() << "--read-likelihood requires reads: pass --gam, --gaf-reads, "
                           << "or --gaf-base" << endl;
        }
        if (read_source_count > 1) {
            logger.error() << "--gam, --gaf-reads, and --gaf-base are mutually exclusive" << endl;
        }
        if (ratio_caller) {
            logger.error() << "--read-likelihood and -B/--bias-mode are mutually exclusive" << endl;
        }
        vector<string> modes;
        if (legacy) {
            modes.emplace_back("--legacy");
        }
        if (bottom_up) {
            modes.emplace_back("--bottom-up");
        }
        if (call_chains) {
            modes.emplace_back("-I/--chains");
        }
        if (gaf_output && !traversals_only) {
            modes.emplace_back("-G/--gaf");
        }
        if (traversals_only) {
            modes.emplace_back("-T/--traversals");
        }
        if (!vcf_filename.empty()) {
            modes.emplace_back("-v/--vcf");
        }
        if (!modes.empty()) {
            logger.error() << "--read-likelihood cannot be used with " << join_with_commas(modes) << endl;
        }
    }

    // Under --read-likelihood on a GBZ, take candidate alleles from the GBZ's haplotypes by default,
    // as -z does, rather than from read support. This needs no pack file, but it can only offer
    // alleles that some haplotype carries. --enumerate-support turns it off.
    if (read_likelihood && !gbz_paths && !enumerate_support && gbz_graph &&
        gbwt_filename.empty()) {
        // A GBZ may carry only reference paths, which would offer nothing but the reference
        // allele, so we only choose this automatically with at least two haplotypes. An
        // explicit -z uses the haplotypes regardless.
        size_t panel = count_panel_haplotypes(*gbz_graph);
        if (panel >= 2) {
            gbz_paths = true;
            gbwt_index = &gbz_graph->gbz.index;
            if (show_progress) {
                logger.info() << "Enumerating alleles from the " << panel
                              << " GBZ panel haplotypes (default under --read-likelihood; "
                              << "--enumerate-support to enumerate from read support instead)"
                              << endl;
            }
            if (!pack_filename.empty()) {
                // Only support-based enumeration uses the pack file.
                logger.warn() << "-k/--pack is unused when alleles come from the haplotype "
                              << "panel; pass --enumerate-support to enumerate from read "
                              << "support and use it" << endl;
            }
        } else {
            logger.warn() << "GBZ carries " << panel << " panel haplotype(s), too few to "
                             << "enumerate alleles from; falling back to support-based "
                             << "enumeration, which needs -k/--pack" << endl;
        }
    }
    if (enumerate_support && is_explicit("gbz")) {
        logger.error() << "--enumerate-support and -z/--gbz select different allele "
                       << "enumeration strategies: choose one or the other" << endl;
    }
    if (enumerate_support && !gbwt_filename.empty()) {
        logger.error() << "--enumerate-support and -g/--gbwt select different allele "
                       << "enumeration strategies: choose one or the other" << endl;
    }

    // The linkage model needs --read-likelihood, haplotypes to enumerate alleles from (-z or -g)
    // and a positive --linkage-weight. We check this before the likelihood calculator is built,
    // so that the calculator does not collect read-phasing evidence that nothing will use. A GBWT
    // with too few haplotypes is found only once it is loaded, below.
    if (const bool haplotypes_given = gbwt_index != nullptr || !gbwt_filename.empty();
        !(read_likelihood && haplotypes_given && linkage_weight > 0.0)) {
        refuse_phase_dependents("and the linkage model needs --read-likelihood, haplotype "
                                "enumeration (-z or -g) and a --linkage-weight above 0");
    }

    if (max_mismap_prob <= 0.0 || max_mismap_prob >= 1.0) {
        logger.error() << "--mismap-max must be in (0, 1)" << endl;
    }
    if (min_mismap_prob <= 0.0 || min_mismap_prob > max_mismap_prob) {
        logger.error() << "--mismap-min must be in (0, --mismap-max]" << endl;
    }
    // --optimal-pairing also finds good pairings for the alleles a read did not come from, which
    // lowers reads' confidence at heterozygous sites, so --phase-min-q has a lower default with it.
    if (optimal_pairing && !is_explicit("phase-min-q")) {
        read_phasing_params.reliability = 8.5;
    }
    // Under read phasing, refuse a --phase-min-q above the heterozygous score ceiling. The ceiling
    // is the confidence of a read that fits one of two equal-length alleles perfectly and the other
    // not at all, phred(e / (e + (1 - e) / 2)) for e = --mismap-min. Sites whose alleles are of
    // similar length cannot reach it, so read phasing would do almost nothing.
    if (read_phasing) {
        const double het_ceiling =
            -10.0 * log10(min_mismap_prob / (min_mismap_prob + (1.0 - min_mismap_prob) / 2.0));
        if (read_phasing_params.reliability > het_ceiling) {
            logger.error() << "--phase-min-q " << read_phasing_params.reliability
                           << " is above the heterozygous score ceiling of " << het_ceiling
                           << " implied by --mismap-min " << min_mismap_prob
                           << "; no site could be reliable" << endl;
        }
    }

    // Validation: -A, --top-down, and --bottom-up are mutually exclusive
    if (int nested_mode_count = (all_snarls ? 1 : 0) + (top_down ? 1 : 0) + (bottom_up ? 1 : 0);
        nested_mode_count > 1) {
        logger.error() << "-A, --top-down, and --bottom-up are mutually exclusive" << endl;
    }

    // Validation for nested calling options
    if (star_allele && !top_down) {
        logger.error() << "-Y/--star-allele requires --top-down mode" << endl;
    }

    // Validation for bottom-up mode
    if (bottom_up && star_allele) {
        logger.error() << "-Y/--star-allele cannot be used with --bottom-up mode" << endl;
    }

    // in order to add subpath support, we let all ref_paths be subpaths and then convert coordinates
    // on VCF export.  the exception is writing the header where we need base paths. we keep
    // track of them the best we can here (just for writing the ##contigs)
    unordered_map<string, size_t> basepath_length_map;

    // call doesn't always require path positions .. .don't change that now
    function<size_t(path_handle_t)> compute_path_length = [&] (path_handle_t path_handle) {
        PathPositionHandleGraph* pp_graph = dynamic_cast<PathPositionHandleGraph*>(graph);
        if (pp_graph) {
            return pp_graph->get_path_length(path_handle);
        } else {
            size_t len = 0;
            graph->for_each_step_in_path(path_handle, [&] (step_handle_t step) {
                    len += graph->get_length(graph->get_handle_of_step(step));
                });
            return len;
        }
    };

    // Add each base path name that -P or -S selects to ref_paths, once. The paths of the graph are
    // matched to ref_paths by base name below, so the name of a subpath, such as NAME[offset] for a
    // gref fragment that vg haplotypes clipped, would never match.
    unordered_set<string> selected_base_paths(ref_paths.begin(), ref_paths.end());
    auto select_base_path = [&](const string& path_name) {
        string base_name = Paths::strip_subrange(path_name);
        if (selected_base_paths.insert(base_name).second) {
            ref_paths.push_back(base_name);
        }
    };

    // Prefix given: find all paths matching it
    if (!ref_path_prefixes.empty()) {
        graph->for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC, PathSense::HAPLOTYPE}, [&](const path_handle_t& path_handle) {
            string path_name = graph->get_path_name(path_handle);
            // Never include alt paths in reference paths
            if (Paths::is_alt(path_name)) {
                return;
            }
            for (const auto& prefix : ref_path_prefixes) {
                if (path_name.compare(0, prefix.size(), prefix) == 0) {
                    select_base_path(path_name);
                    break;
                }
            }
        });
        if (ref_paths.empty()) {
            logger.error() << "No non-alt paths found matching prefix(es) (see vg paths --list)" << endl;
        }
    }

    // Sample given: find all non-alt paths matching it
    if (!ref_sample.empty()) {
        graph->for_each_path_of_sample({ref_sample}, [&](path_handle_t path_handle) {
            const string& name = graph->get_path_name(path_handle);
            if (!Paths::is_alt(name)) {
                select_base_path(name);
            }
        });
        if (ref_paths.empty()) {
            logger.error() << "No REFERENCE or HAPLOTYPE paths for sample \"" << ref_sample << "\" found.\n"
                           << "Use vg paths -M to check which paths exist in this graph\n" 
                           << "Also see: https://github.com/vgteam/vg/wiki/Changing-References" << endl;
        }
    }

    // No paths specified: use all reference/generic
    if (ref_paths.empty()) {
        unordered_set<string> ref_sample_names;
        graph->for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC}, [&](path_handle_t path_handle) {
                const string& name = graph->get_path_name(path_handle);
                if (!Paths::is_alt(name)) {
                    string sample_name = graph->get_sample_name(path_handle);                   
                    ref_paths.push_back(name);
                    // keep track of length best we can using maximum coordinate in event of subpaths
                    
                    // TODO: We can get the subrange from the graph but not
                    // the base path name yet, so we do this from the path
                    // name.
                    subrange_t subrange;
                    string base_name = Paths::strip_subrange(name, &subrange);
                    size_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
                    size_t& cur_len = basepath_length_map[base_name];
                    cur_len = max(cur_len, compute_path_length(path_handle) + offset);
                    if (sample_name != PathMetadata::NO_SAMPLE_NAME) {
                        ref_sample_names.insert(sample_name);
                    }
                }
            });
        if (ref_sample_names.size() > 1) {
            auto err_msg = logger.error();
            err_msg << "Multiple reference samples detected: [";
            size_t count = 0;
            for (const string& n : ref_sample_names) {                
                err_msg << n;
                if (++count >= std::min(ref_sample_names.size(), (size_t)5)) {
                    if (ref_sample_names.size() > 5) {
                        err_msg << ", ...";
                    }
                    break;
                } else {
                    err_msg << ", ";
                }
            }
            err_msg << "]. Please use -S to specify a single reference sample "
                    << "or use -p to specify reference paths" << endl;
        }                
    } else {
        // if paths are given, we convert them to subpaths so that ref paths list corresponds
        // to path names in graph.  subpath handling will only be handled when writing the vcf
        // (this way, we add subpath support without changing anything in between)
        vector<string> ref_subpaths;
        unordered_map<string, bool> ref_path_set;
        for (const string& ref_path : ref_paths) {
            ref_path_set[ref_path] = false;
        }
        graph->for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC, PathSense::HAPLOTYPE}, [&](path_handle_t path_handle) {
                const string& name = graph->get_path_name(path_handle);
                subrange_t subrange;
                string base_name = Paths::strip_subrange(name, &subrange);
                size_t offset = subrange == PathMetadata::NO_SUBRANGE ? 0 : subrange.first;
                if (ref_path_set.count(base_name)) {
                    ref_subpaths.push_back(name);
                    // keep track of length best we can
                    if (ref_path_lengths.empty()) {
                        size_t& cur_len = basepath_length_map[base_name];
                        cur_len = max(cur_len, compute_path_length(path_handle) + offset);
                    }
                    ref_path_set[base_name] = true;
                }
            });

        // if we have reference lengths, great!
        // this will be the only way to get a correct header in the presence of supbpaths
        if (!ref_path_lengths.empty()) {
            assert(ref_path_lengths.size() == ref_paths.size());
            for (size_t i = 0; i < ref_paths.size(); ++i) {
                basepath_length_map[ref_paths[i]] = ref_path_lengths[i];
            }
        }

        // Check our paths
        for (const auto& ref_path_used : ref_path_set) {
            if (!ref_path_used.second) {
                logger.error() << "Path \"" << ref_path_used.first 
                               << "\" not found in graph as a non-alt sense (see vg paths -M)\n"
                               << "Also see: https://github.com/vgteam/vg/wiki/Changing-References" << endl;
            }
        }
        
        swap(ref_paths, ref_subpaths);
    }

    // make sure we have some ref paths
    if (ref_paths.empty()) {
        logger.error() << "No reference paths found. "
                       << "Paths must be REFERENCE or GENERIC sense (see vg paths -M)\n"
                       << "Alternatively, use --ref-path, --path-prefix, or --ref-sample to force a HAPLOTYPE path to be treated as a reference\n"
                       << "Also see: https://github.com/vgteam/vg/wiki/Changing-References" << endl;
    }

    for (const string& ref_path : ref_paths) {
        // Only a fragment counts. The gref copy of a base contig is the reference under another
        // name, and gives no chain a contig of its own. A clipped fragment is a subpath, so its
        // name ends in [offset].
        if (GrefCover::is_gref_name(Paths::strip_subrange(ref_path))) {
            off_ref_nesting = true;
            if (show_progress) {
                logger.info() << "gref reference selected: descending into chains the reference "
                              << "does not cross, and reporting them against their gref contig"
                              << endl;
            }
            break;
        }
    }

    // Off-reference chains have anchors even without VCF records, so --anchors-out turns their
    // genotyping on unless --no-off-ref-nesting is given.
    if (!anchors_out.empty() && !is_explicit("off-ref-nesting")) {
        off_ref_nesting = true;
        if (show_progress) {
            logger.info() << "anchors requested: descending into chains the reference does not "
                          << "cross, which have no VCF record but do have anchors" << endl;
        }
    }

    // The gRef nesting level of each selected gRef contig, for INFO/CH, by the contig name the VCF
    // uses.
    map<string, int> gref_levels;
    if (off_ref_nesting && graph != nullptr) {
        gref_levels = gref_contig_levels(*graph, ref_paths);
        if (show_progress && !gref_levels.empty()) {
            map<int, size_t> hist;
            for (const auto& [contig, level] : gref_levels) {
                ++hist[level];
            }
            auto& l = logger.info() << "gref nesting levels:";
            for (const auto& [level, count] : hist) {
                l << " " << level << ":" << count;
            }
            l << endl;
        }
    }

    // build table of ploidys
    vector<int> ref_path_ploidies;
    // Paths which aren't REFERENCE/GENERIC sense that we want to call against
    unordered_set<string> pretend_ref_paths;
    for (const string& ref_path : ref_paths) {
        int path_ploidy = ploidy;
        for (auto& rule : ploidy_rules) {
            if (std::regex_match(ref_path, rule.first)) {
                path_ploidy = rule.second;
                break;
            }
        }
        // The callers only implement ploidy 1 and 2. We check -d for this above, and each -R rule
        // here, since an unsupported ploidy would abort in the caller.
        if (path_ploidy != 1 && path_ploidy != 2) {
            logger.error() << "ploidy " << path_ploidy << " assigned to path \"" << ref_path
                           << "\" by -R/--ploidy-regex must be 1 or 2" << endl;
        }
        ref_path_ploidies.push_back(path_ploidy);

        if (graph->get_sense(graph->get_path_handle(ref_path)) == PathSense::HAPLOTYPE) {
            pretend_ref_paths.emplace(ref_path);
        }
    }

    // Use an overlay so that all ref paths are treated as refs
    bdsg::ReferencePathOverlayHelper overlay_helper;
    if (!pretend_ref_paths.empty()) {
        if (show_progress) logger.info() << "Applying overlay to treat HAPLOTYPE paths as REFERENCE" << endl;
        graph = overlay_helper.apply(graph, pretend_ref_paths);
    }

    // Load or compute the snarls
    unique_ptr<SnarlManager> snarl_manager;    
    if (!snarl_filename.empty()) {
        ifstream snarl_file(snarl_filename.c_str());
        if (show_progress) logger.info() << "Loading snarls from " << snarl_filename << endl;
        snarl_manager = vg::io::VPKG::load_one<SnarlManager>(snarl_file);
        if (show_progress) logger.info() << "Loaded snarls" << endl;
    } else {
        if (show_progress) logger.info() << "Computing snarls" << endl;
        std::unordered_map<nid_t, size_t> extra_node_weight;
        constexpr size_t EXTRA_WEIGHT = 10000000000;
        for (const string& refpath_name : ref_paths) {
            // Skip gref fragments, clipped or not (they shouldn't influence snarl decomposition)
            if (GrefCover::is_gref_name(Paths::strip_subrange(refpath_name))) {
                continue;
            }
            path_handle_t refpath_handle = graph->get_path_handle(refpath_name);
            extra_node_weight[graph->get_id(graph->get_handle_of_step(graph->path_begin(refpath_handle)))] += EXTRA_WEIGHT;
            extra_node_weight[graph->get_id(graph->get_handle_of_step(graph->path_back(refpath_handle)))] += EXTRA_WEIGHT;
        }        
        IntegratedSnarlFinder finder(*graph, extra_node_weight);
        // The decomposition happens in find_snarls_parallel(), so we report after it.
        snarl_manager = unique_ptr<SnarlManager>(new SnarlManager(std::move(finder.find_snarls_parallel())));
        if (show_progress) logger.info() << "Computed snarls" << endl;
    }
    
    // Make a Packed Support Caller
    unique_ptr<SnarlCaller> snarl_caller;
    // The read-likelihood genotyper.
    unique_ptr<ReadLikelihoodSnarlCaller> read_likelihood_caller;
    vg::algorithms::BinnedDepthIndex depth_index;

    unique_ptr<Packer> packer;
    unique_ptr<TraversalSupportFinder> support_finder;
    // Only used by --read-likelihood, but declared out here so they outlive the
    // caller, which holds references to them.
    unique_ptr<SiteReadSource> read_source;
    unique_ptr<EditAlignmentScorer> qual_scorer;
    unique_ptr<EditAlignmentScorer> plain_scorer;
    unique_ptr<AlleleLikelihoodCalculator> likelihood_calculator;
    unique_ptr<ofstream> likelihood_dump;
    // The read-likelihood genotyper can run without a pack file if allele enumeration needs no
    // support either. GBWTTraversalFinder enumerates from haplotypes, while FlowTraversalFinder
    // needs node and edge support.
    bool gbwt_enumeration = !gbwt_filename.empty() || gbz_paths;
    bool support_free = read_likelihood && pack_filename.empty();

    if (!pack_filename.empty() || support_free) {
        if (!support_free) {
            // Load our packed supports (they must have come from vg pack on graph)
            packer = unique_ptr<Packer>(new Packer(graph));
            if (show_progress) logger.info() << "Loading pack file " << pack_filename << endl;
            packer->load_from_file(pack_filename);
            if (show_progress) logger.info() << "Loaded pack file" << endl;
            if (bottom_up) {
                // Make a nested packed traversal support finder (required by NestedFlowCaller)
                support_finder.reset(new NestedCachedPackedTraversalSupportFinder(*packer, *snarl_manager));
            } else {
                // Make a packed traversal support finder (using cached version important for
                // poisson caller)
                support_finder.reset(new CachedPackedTraversalSupportFinder(*packer, *snarl_manager));
            }
        }
                
        if (support_finder) {
            // need to use average support when genotyping as small differences in between sample and graph
            // will lead to spots with 0-support, espeically in and around SVs. 
            support_finder->set_support_switch_threshold(avg_trav_threshold, avg_node_threshold);

            // upweight breakpoint edges even when taking average support otherwise
            support_finder->set_min_bp_edge_override(expect_bp_edges);
        }

        // todo: toggle between min / average (or thresholds) via command line
        
        SupportBasedSnarlCaller* packed_caller = nullptr;

        if (read_likelihood) {
            // Read-likelihood genotyping. Support from a pack file, if given, is only used to
            // enumerate candidate alleles.
            if (show_progress) logger.info() << "Loading reads for read-level genotyping" << endl;

            SiteReadFilter read_filter;
            read_filter.min_mapq = read_min_mapq;

            if (reads.window_size == 0) {
                reads.window_size = default_read_window(reads);
            }
            read_source = open_site_read_source(*graph, reads, read_filter, show_progress, logger);

            // Two scorers: quality-adjusted for reads that have base qualities, and plain for
            // reads that do not. We pick one per read, so we neither invent qualities for a read
            // nor ignore the ones it has.
            qual_scorer.reset(new QualAdjAlignmentScorer(default_score_matrix,
                                                        (int8_t)gap_open,
                                                        (int8_t)gap_extend));
            plain_scorer.reset(new MatrixAlignmentScorer(default_score_matrix,
                                                        (int8_t)gap_open,
                                                        (int8_t)gap_extend));

            AlleleLikelihoodParams likelihood_params;
            likelihood_params.use_mismap_term = !no_mismap_term;
            likelihood_params.length_weighted_mixture = !flat_mixture;
            likelihood_params.depth_weight = depth_weight;
            likelihood_params.depth_effective_reads = !depth_count_raw;
            likelihood_params.max_mismap_prob = max_mismap_prob;
            likelihood_params.min_mismap_prob = min_mismap_prob;
            likelihood_params.insertion_gap_nats = insertion_gap_nats;
            likelihood_params.optimal_pairing = optimal_pairing;
            likelihood_params.collect_anchors = anchor_params.enabled;
            // The likelihood calculator and the anchor writer count into the same counters.
            likelihood_params.anchor_counters = &anchor_run_counters;
            anchor_params.counters = &anchor_run_counters;
            likelihood_params.collect_read_phasing = read_phasing;

            auto owned_calculator = make_unique<GraphAlignedAlleleLikelihoodCalculator>(
                *graph, *read_source, *qual_scorer, *plain_scorer,
                likelihood_params);
            GraphAlignedAlleleLikelihoodCalculator* graph_calculator = owned_calculator.get();
            likelihood_calculator = std::move(owned_calculator);
            // Place the depth-rate windows on the reference paths, so that they do not depend
            // on the node numbering. A graph without positions keeps node-ID windows.
            if (const auto* position_graph = dynamic_cast<const PathPositionHandleGraph*>(graph)) {
                vector<path_handle_t> rate_paths;
                for (const string& ref_path : ref_paths) {
                    if (position_graph->has_path(ref_path)) {
                        rate_paths.push_back(position_graph->get_path_handle(ref_path));
                    }
                }
                graph_calculator->set_rate_reference(position_graph, rate_paths);
            }

            auto rl_caller = new ReadLikelihoodSnarlCaller(*graph, *likelihood_calculator);

            if (!dump_likelihoods_filename.empty()) {
                likelihood_dump = make_unique<ofstream>(dump_likelihoods_filename);
                if (!(*likelihood_dump)) {
                    logger.error() << "could not open " << dump_likelihoods_filename
                                   << " for writing" << endl;
                }
                rl_caller->set_likelihood_dump(likelihood_dump.get());
            }

            rl_caller->set_share_discount(!no_share_quality);
            rl_caller->set_depth_quality(depth_quality);
            rl_caller->set_min_confidence(min_confidence);

            read_likelihood_caller.reset(rl_caller);
        } else if (ratio_caller == false) {
            // Make a depth index
            if (show_progress) logger.info() << "Computing coverage statistics" << endl;
            depth_index = vg::algorithms::binned_packed_depth_index(*packer, ref_paths, min_depth_bin_width, max_depth_bin_width,
                                                                depth_scale_fac, 0, true, true);
            if (show_progress) logger.info() << "Computed coverage statistics" << endl;
            // Make a new-stype probablistic caller
            auto poisson_caller = new PoissonSupportSnarlCaller(*graph, *snarl_manager, *support_finder, depth_index,
                                                                //todo: qualities need to be used better in conjunction with
                                                                //expected depth.
                                                                //packer->has_qualities());
                                                                false);

            // Pass the errors through
            poisson_caller->set_baseline_error(baseline_error_small, baseline_error_large);
                
            packed_caller = poisson_caller;
        } else {
            // Make an old-style ratio support caller
            auto ratio_caller = new RatioSupportSnarlCaller(*graph, *snarl_manager, *support_finder);
            if (het_bias >= 0) {
                ratio_caller->set_het_bias(het_bias, ref_het_bias);
            }
            packed_caller = ratio_caller;
        }
        if (packed_caller != nullptr) {
            if (min_allele_support >= 0) {
                packed_caller->set_min_supports(min_allele_support, min_allele_support, min_site_support);
            }
            snarl_caller = unique_ptr<SnarlCaller>(packed_caller);
        }
    }

    if (!snarl_caller && !read_likelihood_caller) {
        logger.error() << "pack file (-k) is required" << endl;
    }

    // Running without a pack file is only sound where nothing consults support.
    if (support_free && !gbwt_enumeration) {
        logger.error() << "--read-likelihood without -k/--pack requires haplotype-based allele "
                       << "enumeration (-g/--gbwt or -z/--gbz); support-based enumeration needs "
                       << "a pack file" << endl;
    }

    unique_ptr<AlignmentEmitter> alignment_emitter;
    if (gaf_output) {
        alignment_emitter = vg::io::get_non_hts_alignment_emitter("-", "GAF", {}, vg::get_thread_count(), graph);
        // TODO: There should be a general function for emitting headers. See giraffe_main.cpp.
        io::GafAlignmentEmitter* gaf_emitter = dynamic_cast<io::GafAlignmentEmitter*>(alignment_emitter.get());
        if (gbz_graph.get() != nullptr && gaf_emitter != nullptr) {
            gbwtgraph::GraphName graph_name = gbz_graph->gbz.graph_name();
            std::vector<std::string> header_lines = graph_name.gaf_header_lines();
            gaf_emitter->emit_header_lines(header_lines);
        }
    }

    // Every run without --read-likelihood calls through a GraphCaller, and every run with it
    // through a MultiPassCaller, which writes its records through `multipass_output`.
    unique_ptr<GraphCaller> graph_caller;
    unique_ptr<SnarlManagerDecomposition> multipass_sites;
    unique_ptr<SortedVCFWriter> multipass_output;
    unique_ptr<multipass::MultiPassCaller> multipass_caller;
    unique_ptr<TraversalFinder> traversal_finder;
    unique_ptr<gbwt::GBWT> gbwt_index_up;

    vcflib::VariantCallFile variant_file;
    unique_ptr<FastaReference> ref_fasta;
    unique_ptr<FastaReference> ins_fasta;
    if (!vcf_filename.empty()) {
        // Genotype the VCF
        variant_file.parseSamples = false;
        variant_file.open(vcf_filename);
        if (!variant_file.is_open()) {
            logger.error() << "could not open " << vcf_filename << endl;
        }

        // load up the fasta
        if (!ref_fasta_filename.empty()) {
            ref_fasta = unique_ptr<FastaReference>(new FastaReference);
            ref_fasta->open(ref_fasta_filename);
        }
        if (!ins_fasta_filename.empty()) {
            ins_fasta = unique_ptr<FastaReference>(new FastaReference);
            ins_fasta->open(ins_fasta_filename);
        }
        
        VCFGenotyper* vcf_genotyper = new VCFGenotyper(*graph, *snarl_caller,
                                                       *snarl_manager, variant_file,
                                                       sample_name, ref_paths, ref_path_ploidies,
                                                       ref_fasta.get(),
                                                       ins_fasta.get(),
                                                       alignment_emitter.get(),
                                                       traversals_only,
                                                       gaf_output,
                                                       trav_padding);
        graph_caller = unique_ptr<GraphCaller>(vcf_genotyper);
    } else if (legacy) {
        // de-novo caller (port of the old vg call code, which requires a support based caller)
        LegacyCaller* legacy_caller = new LegacyCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                                       *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                                       *snarl_manager,
                                                       sample_name, ref_paths, ref_path_offsets, ref_path_ploidies);
        graph_caller = unique_ptr<GraphCaller>(legacy_caller);
    } else {
        // flow caller can take any kind of traversal finder.  two are supported for now:
        
        if (!gbwt_filename.empty() || gbz_paths) {
            // GBWT traversals
            if (!gbz_paths) {
                gbwt_index_up = vg::io::VPKG::load_one<gbwt::GBWT>(gbwt_filename);
                gbwt_index = gbwt_index_up.get();
                if (gbwt_index == nullptr) {
                    logger.error() << "unable to load GBWT index from file: " << gbwt_filename << endl;
                }
            }
            GBWTTraversalFinder* gbwt_traversal_finder = new GBWTTraversalFinder(*graph, *gbwt_index);
            traversal_finder = unique_ptr<TraversalFinder>(gbwt_traversal_finder);
        } else {
            // Flow traversals (Yen's algorithm)
            if (!support_finder) {
                logger.error() << "finding candidate alleles by flow needs read support from a pack "
                               << "file (-k)" << endl;
            }
            
            // todo: do we ever want to toggle in min-support?
            function<double(handle_t)> node_support = [&] (handle_t h) {
                return support_finder->support_val(support_finder->get_avg_node_support(graph->get_id(h)));
            };
            
            function<double(edge_t)> edge_support = [&] (edge_t e) {
                return support_finder->support_val(support_finder->get_edge_support(e));
            };

            // create the flow traversal finder
            // A site longer than the average support threshold is searched by average flow too.
            traversal_finder = make_unique<FlowTraversalFinder>(*graph, max_yens_traversals, node_support, edge_support,
                                                                max_allele_len, avg_trav_threshold);
        }

        if (read_likelihood) {
            // The caller reads the sites through the SnarlDecomposition interface. The adapter over
            // the loaded snarls gives it the same sites as every other caller.
            try {
                multipass_sites.reset(new SnarlManagerDecomposition(*snarl_manager, *graph));
            } catch (const std::invalid_argument&) {
                // The view throws only on a unary snarl, which only the cactus snarl finder makes.
                logger.error() << "the snarls include a unary snarl, which --read-likelihood cannot "
                               << "call; find the snarls with vg snarls -A integrated, which makes none"
                               << endl;
            }
            multipass_output.reset(new SortedVCFWriter(sample_name));
            multipass_caller.reset(new multipass::MultiPassCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                                       *read_likelihood_caller, *multipass_sites,
                                                       *multipass_output, sample_name,
                                                       dynamic_cast<SiteTraversalFinder&>(*traversal_finder),
                                                       multipass::ReferencePathTable(
                                                           ref_paths, ref_path_offsets,
                                                           ref_path_ploidies),
                                                       genotype_snarls,
                                                       make_pair(min_allele_len, max_allele_len),
                                                       top_down, star_allele));
        } else if (top_down) {
            // Use FlowCaller with nested mode enabled (top-down genotype propagation)
            graph_caller.reset(new FlowCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                              *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                              *snarl_manager,
                                              sample_name, *traversal_finder, ref_paths, ref_path_offsets,
                                              ref_path_ploidies,
                                              alignment_emitter.get(),
                                              traversals_only,
                                              gaf_output,
                                              trav_padding,
                                              genotype_snarls,
                                              make_pair(min_allele_len, max_allele_len),
                                              true,  // nested mode enabled
                                              star_allele));
        } else if (bottom_up) {
            // Use NestedFlowCaller (bottom-up snarl merging, original nested algorithm)
            graph_caller.reset(new NestedFlowCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                                    *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                                    *snarl_manager,
                                                    sample_name, *traversal_finder, ref_paths, ref_path_offsets,
                                                    ref_path_ploidies,
                                                    alignment_emitter.get(),
                                                    traversals_only,
                                                    gaf_output,
                                                    trav_padding,
                                                    genotype_snarls));
        } else {
            graph_caller.reset(new FlowCaller(*dynamic_cast<PathPositionHandleGraph*>(graph),
                                              *dynamic_cast<SupportBasedSnarlCaller*>(snarl_caller.get()),
                                              *snarl_manager,
                                              sample_name, *traversal_finder, ref_paths, ref_path_offsets,
                                              ref_path_ploidies,
                                              alignment_emitter.get(),
                                              traversals_only,
                                              gaf_output,
                                              trav_padding,
                                              genotype_snarls,
                                              make_pair(min_allele_len, max_allele_len)));
        }
    }

    // By default, snarls have no size cap when the read-likelihood genotyper takes alleles from
    // haplotypes. Everywhere else they are capped at 10000 edges, which keeps Yen's traversal
    // search off very large snarls.
    if (!is_explicit("max-snarl-edges")) {
        max_snarl_edges_opt = (read_likelihood && gbwt_enumeration) ? 0 : 10000;
    }
    // The constructors do not take the cap, so we set it on whichever caller was built.
    if (auto* flow_caller = dynamic_cast<FlowCaller*>(graph_caller.get())) {
        flow_caller->set_max_snarl_edges(max_snarl_edges_opt);
    }
    if (multipass_caller != nullptr) {
        multipass_caller->set_max_snarl_edges(max_snarl_edges_opt);
    }

    // Where the records go: the writer of whichever caller was built.
    SortedVCFWriter* vcf_out = nullptr;
    if (multipass_caller != nullptr) {
        vcf_out = multipass_output.get();
    } else if (auto* flow_caller = dynamic_cast<FlowCaller*>(graph_caller.get())) {
        vcf_out = &flow_caller->vcf_writer;
    } else if (auto* nested_caller = dynamic_cast<NestedFlowCaller*>(graph_caller.get())) {
        vcf_out = &nested_caller->vcf_writer;
    } else if (auto* legacy_caller = dynamic_cast<LegacyCaller*>(graph_caller.get())) {
        vcf_out = &legacy_caller->vcf_writer;
    } else if (auto* vcf_genotyper = dynamic_cast<VCFGenotyper*>(graph_caller.get())) {
        vcf_out = &vcf_genotyper->vcf_writer;
    }

    // Give the per-region ploidy, if any, to the caller.
    if (ploidy_regions.has_value()) {
        if (multipass_caller != nullptr) {
            multipass_caller->set_ploidy_regions(std::move(*ploidy_regions));
        } else {
            graph_caller->set_ploidy_regions(std::move(*ploidy_regions));
        }
    }

    // Nested calling: a called traversal that takes the reference's route through a snarl,
    // differing only inside nested snarls, is called as the reference allele, and the differences
    // are called at the nested snarls. It is part of the read-likelihood genotyper, and on by
    // default there.
    if (!read_likelihood) {
        nested_calling = false;
    }
    // An explicit --regenotype without read phasing is an error. One set by a preset is turned
    // off, as with --preset ont --no-read-phasing.
    if (regenotype && !read_phasing) {
        if (is_explicit("regenotype")) {
            cerr << "error [vg call]: --regenotype needs --read-phasing, which gives each read"
                 << " the strand log-odds it uses" << endl;
            return 1;
        }
        regenotype = false;
    }
    // --top-down gives each child snarl candidate traversals derived from its parent's called
    // genotype, so re-genotyping, which can change that genotype afterwards, would leave the child
    // genotyped against the wrong alleles. An explicit --regenotype with it is an error. One set by
    // a preset is turned off.
    if (regenotype && top_down && !is_explicit("regenotype")) {
        regenotype = false;
    }
    if (regenotype && top_down) {
        cerr << "error [vg call]: --regenotype cannot be combined with --top-down; --top-down"
             << " derives each child's candidate traversals from its parent's called genotype, so"
             << " re-genotyping a parent would leave its children genotyped against alleles it no"
             << " longer carries" << endl;
        return 1;
    }
    // Splitting homozygous sites, placing heterozygous reads by phase, and re-genotyping all use
    // each read's strand log-odds, which read phasing computes.
    if (anchor_params.hom_split && !read_phasing) {
        cerr << "error [vg call]: --anchors-hom-split needs --read-phasing, which gives each read"
             << " the strand log-odds the split uses" << endl;
        return 1;
    }
    // Placing heterozygous reads by phase is on by default, and has no effect without read
    // phasing, so we only refuse it when it was asked for.
    if ((anchor_params.phase_hets || anchor_params.strict_hets) && (is_explicit("anchors-phase-hets") || is_explicit("anchors-strict-hets"))
        && !read_phasing) {
        cerr << "error [vg call]: --anchors-phase-hets / --anchors-strict-hets need --read-phasing,"
             << " which gives each read the strand log-odds they use" << endl;
        return 1;
    }
    // -A already calls every nested snarl as a snarl of its own, so with nested calling each would
    // be called twice.
    if (nested_calling && all_snarls) {
        if (is_explicit("nested")) {
            cerr << "error [vg call]: -A/--all-snarls calls every snarl independently while"
                 << " --nested calls children through their parent; they are alternatives, not a"
                 << " combination" << endl;
            return 1;
        }
        nested_calling = false;
    }
    if (multipass_caller != nullptr) {
        // The linkage model can change a record's GQN, so the caller re-applies the lowconf
        // filter with the same threshold.
        multipass_caller->set_linkage_min_confidence(min_confidence);
    }
    // --atomize-blocks splits up the records of nested calling, so it needs nested calling. The
    // checks that depend only on the options were made before the graph was loaded.
    if (atomize_blocks) {
        if (!nested_calling) {
            if (is_explicit("atomize-blocks")) {
                cerr << "error [vg call]: --atomize-blocks needs nested calling (--nested), which is"
                     << " on by default under --read-likelihood" << endl;
                return 1;
            }
            atomize_blocks = false;
        }
        if (atomize_blocks) {
            multipass_caller->set_atomize_blocks(true);
        }
    }

    if (nested_calling) {
        multipass_caller->set_nested_calling(true);
    }

    // Owned here, since the caller holds only a pointer to it.
    unique_ptr<LinkageCollector> linkage_collector;
    vector<size_t> linkage_sequence_to_haplotype;

    // Write the nesting INFO tags (such as LV and PS) with -A, --top-down or --bottom-up, and
    // when off-reference chains are genotyped. A record on a gref fragment contig lies inside a
    // non-reference allele, and only these tags say so.
    bool nesting_tags = all_snarls || top_down || bottom_up || off_ref_nesting;
    // The snarl tree the nesting tags are read from, for a caller other than the read-likelihood
    // one, which has its own. It is made now so that a unary snarl is refused before calling.
    unique_ptr<SnarlManagerDecomposition> nesting_sites;
    if (nesting_tags && !gaf_output && multipass_caller == nullptr) {
        try {
            nesting_sites.reset(new SnarlManagerDecomposition(*snarl_manager, *graph));
        } catch (const std::invalid_argument&) {
            // The view throws only on a unary snarl, which only the cactus snarl finder makes.
            logger.error() << "the snarls include a unary snarl, which the nesting INFO tags "
                           << "cannot place; find the snarls with vg snarls -A integrated, which "
                           << "makes none" << endl;
        }
    }

    if (!gaf_output) {
        // Init The VCF
        assert(vcf_out != nullptr);
        vcf_out->set_nested(nesting_tags);
        if (multipass_caller != nullptr) {
            multipass_caller->set_off_reference_nesting(off_ref_nesting);
        }
        vcf_out->set_gref_levels(std::move(gref_levels));
        vcf_out->set_translation(translation.get());

        // The linkage model compares genotypes with the haplotypes (-z or -g), and uses the
        // read-likelihood genotyper's likelihoods. If either is missing, the default weight is set
        // to 0, and an explicit weight is an error. This is decided first because phasing, the mosaic
        // and nested calling depend on it.
        if (linkage_weight > 0.0 && !(gbwt_index != nullptr && read_likelihood)) {
            if (is_explicit("linkage-weight")) {
                cerr << "error [vg call]: --linkage-weight needs haplotype enumeration (-z or -g) "
                     << "and --read-likelihood" << endl;
                return 1;
            }
            linkage_weight = 0.0;
        }
        if (!anchors_out.empty() && !read_likelihood) {
            // Anchors need the read-likelihood genotyper's per-read allele likelihoods.
            logger.error() << "--anchors-out needs --read-likelihood, which computes the per-read "
                           << "allele likelihoods the anchors are built from" << endl;
        }
        if (!mosaic_out.empty() && linkage_weight <= 0.0) {
            // The mosaic is built from the linkage model's phasing.
            logger.error() << "--mosaic-out needs the linkage model, which needs haplotype "
                           << "enumeration (-z or -g) and --read-likelihood" << endl;
        }
        if (!mosaic_out.empty() && is_explicit("phased") && !phased_output) {
            // The mosaic is written from the phasing, so it cannot be written without it.
            logger.error() << "--mosaic-out cannot be combined with --no-phased, because the mosaic "
                           << "is written from the phasing" << endl;
        }
        if (phased_output && linkage_weight <= 0.0) {
            // Phasing comes from the linkage model. Without it, an explicit --phased is an error,
            // and the default is turned off.
            if (is_explicit("phased")) {
                logger.error() << "--phased needs the linkage model, which needs haplotype "
                               << "enumeration (-z or -g) and --read-likelihood" << endl;
            }
            phased_output = false;
        }
        if (nested_calling && linkage_weight > 0.0 && !phased_output) {
            // Where the linkage model runs, a nested site's ploidy and strand come from its
            // parent's phased genotype, so nested calling needs phasing. An explicit --nested is
            // an error. The default is turned off.
            if (is_explicit("nested")) {
                cerr << "error [vg call]: --nested needs --phased where the linkage model runs, "
                     << "because a nested site's ploidy and strand come from its parent's phased "
                     << "genotype" << endl;
                return 1;
            }
            // Nested calling was turned on in the caller above, so it is turned off there too.
            nested_calling = false;
            multipass_caller->set_nested_calling(false);
            if (show_progress) {
                logger.info() << "Nested calling is off under --no-phased: a nested site's ploidy "
                              << "and strand come from its parent's phased genotype" << endl;
            }
        }
        if (linkage_weight > 0.0) {
            const map<PanelHaplotype, size_t> hap_index =
                number_panel_haplotypes(*gbwt_index, linkage_sequence_to_haplotype);
            // The model's states are pairs of haplotypes, so it needs at least two.
            if (hap_index.size() < 2) {
                cerr << "warning [vg call]: linkage disabled -- the GBWT carries "
                     << hap_index.size() << " haplotype(s) and the model needs at least 2" << endl;
                linkage_weight = 0.0;
                refuse_phase_dependents("and the linkage model needs at least 2 haplotypes");
            } else {
                // The default weight suits panels of tens of haplotypes, and does little on very
                // small ones.
                if (hap_index.size() < 8) {
                    cerr << "warning [vg call]: linkage on a " << hap_index.size()
                         << "-haplotype panel; the default --linkage-weight suits panels of tens of"
                         << " haplotypes and does little on small ones (consider --linkage-weight 0)"
                         << endl;
                }
                LinkageModel::Params linkage_params;
                linkage_params.weight = linkage_weight;
                linkage_params.scale = linkage_scale;
                linkage_params.freq_prior = linkage_freq_prior;
                linkage_params.mutation = linkage_mutation;
                linkage_params.hp_prior = hp_prior;
                linkage_params.hp_prior_run = (size_t)hp_prior_run;
                linkage_collector.reset(new LinkageCollector(linkage_params, hap_index.size()));
                multipass_caller->set_linkage(linkage_collector.get(), gbwt_index,
                                              &linkage_sequence_to_haplotype);
                multipass_caller->set_emit_phasing(phased_output);
                vector<string> hap_names = panel_haplotype_names(gbwt_index->metadata, hap_index);
                // The mosaic's rows name only the contig, so it also records the full names of
                // the reference paths.
                multipass_caller->set_mosaic_out(multipass::MosaicParams{
                    .path = mosaic_out,
                    .graph_name = graph_filename,
                    .haplotype_names = hap_names,
                    .reference_paths = ref_paths,
                    .patch_gaps = mosaic_patch_gaps,
                    .keep_nested = mosaic_keep_nested,
                    .connect_unexplained = mosaic_connect_unexplained,
                });
            }
            if (show_progress) {
                logger.info() << "Linkage: " << hap_index.size() << " panel haplotypes over "
                              << gbwt_index->sequences() << " GBWT sequences" << endl;
            }
        }
        if (!anchors_out.empty()) {
            // The anchor file names the read file, since its read offsets refer to those reads.
            multipass_caller->set_anchors_out(anchors_out, anchor_params, graph_filename,
                                              reads.reads_file(), min_mismap_prob);
        }
        if (multipass_caller != nullptr) {
            multipass_caller->set_read_phasing(read_phasing, read_phasing_params);
            if (off_panel) {
                // The reads are counted as they are fetched, from the first site on.
                EditCountParams count_params;
                count_params.large_indels = off_panel_svs;
                count_params.bar = &off_panel_bar;
                read_source->count_edits(graph, count_params);
                multipass_caller->set_off_panel(read_source.get(), off_panel_dump, off_panel_bar);
            }
            multipass_caller->set_read_source(read_source.get());
            multipass_caller->set_regenotype(regenotype, regenotype_params, regenotype_passes,
                                             regenotype_ledger);
        }
        // Merge similar called ALT alleles (-L/--cluster) in the records of each caller that
        // builds its records with a SiteEmitter.
        SiteEmitter* site_emitter = nullptr;
        if (multipass_caller != nullptr) {
            site_emitter = &multipass_caller->emitter();
        } else if (auto* flow_caller = dynamic_cast<FlowCaller*>(graph_caller.get())) {
            site_emitter = &flow_caller->site_emitter;
        } else if (auto* nested_caller = dynamic_cast<NestedFlowCaller*>(graph_caller.get())) {
            site_emitter = &nested_caller->site_emitter;
        } else if (auto* legacy_caller = dynamic_cast<LegacyCaller*>(graph_caller.get())) {
            site_emitter = &legacy_caller->site_emitter;
        }
        if (site_emitter != nullptr) {
            site_emitter->allele_merge_threshold = cluster_threshold;
            site_emitter->allele_merge_min_len = cluster_min_allele_len;
        }
        // The fields of the read-likelihood caller come before MAT, the field that records merged
        // alleles, and the genotyper's fields come after it.
        if (multipass_caller != nullptr) {
            multipass_caller->add_header_lines();
        }
        if (cluster_threshold < 1.0) {
            vcf_out->add_info_line(
                "##INFO=<ID=MAT,Number=.,Type=String,Description=\"Merged Allele Traversal: "
                "ALT alleles merged after genotyping by -L/--cluster, as OLD>NEW:SIMILARITY using "
                "pre-merge allele numbers. AD and GL are folded onto the surviving allele and MAD "
                "is recomputed; DP, QUAL, GQ, GP and FILTER are as computed over the pre-merge allele "
                "set. In a nested run this record gives the collapsed view of the site and its "
                "child records the precise one, so they disagree by design.\">");
        }
        if (multipass_caller != nullptr) {
            read_likelihood_caller->add_header_lines(*vcf_out);
        } else {
            snarl_caller->add_header_lines(*vcf_out);
        }
        // The header declares the base paths found above, by the names the records have in their
        // CHROM column, rather than any path of the graph.
        auto* vcf_genotyper = dynamic_cast<VCFGenotyper*>(graph_caller.get());
        vector<string> contig_names;
        vector<size_t> contig_lengths;
        for (const auto& [path_name, length] : basepath_length_map) {
            // Use the locus name (eg "chrI") to match what SiteEmitter::emit()
            // writes to the CHROM column in data lines
            string contig_name = PathMetadata::parse_locus_name(path_name);
            contig_names.push_back(contig_name != PathMetadata::NO_LOCUS_NAME ? contig_name
                                                                             : path_name);
            contig_lengths.push_back(length);
        }
        if (vcf_genotyper != nullptr) {
            // The VCF being genotyped gives the lengths of its contigs.
            contig_lengths = vcf_genotyper->contig_lengths(contig_names);
        }
        vector<SortedVCFWriter::ReferenceContig> contigs;
        for (size_t i = 0; i < contig_names.size(); ++i) {
            contigs.push_back({contig_names[i], contig_lengths[i]});
        }
        vcf_out->set_contigs(std::move(contigs), true);
    }

    if (multipass_caller != nullptr) {
        multipass_caller->set_show_progress(show_progress);
    } else {
        graph_caller->set_show_progress(show_progress);
    }
    
    // Call the graph
    // Determine recursion strategy based on mode:
    // - top_down: the caller handles recursion internally, so RecurseNever
    // - all_snarls (-A): visit every snarl independently, so RecurseAlways
    // - default: only recurse into children of failed snarls, so RecurseOnFail
    GraphCaller::RecurseType recurse_type;
    if (top_down) {
        recurse_type = GraphCaller::RecurseNever;
    } else if (all_snarls) {
        recurse_type = GraphCaller::RecurseAlways;
    } else {
        recurse_type = GraphCaller::RecurseOnFail;
    }

    // A read source that fetches reads by node-ID window works best when snarls are called in
    // node-ID order, batched by its window, so that each fetched window serves many sites in a row.
    if (dynamic_cast<WindowedSiteReadSource*>(read_source.get()) != nullptr) {
        // Only a run with --read-likelihood has a read source, and it always has multipass_caller.
        multipass_caller->set_snarl_batching(reads.window_size);
        if (show_progress) {
            logger.info() << "Visiting snarls in node-ID order, window " << reads.window_size << endl;
        }
    }

    if (multipass_caller != nullptr) {
        // Run every pass of the multipass caller.
        if (show_progress) logger.info() << "Calling top-level snarls" << endl;
        multipass_caller->call(recurse_type, [&]() {
            // The passes after the first use the read likelihoods that the first pass stored, and
            // read no more reads. So we free the read source's cached windows, which use the most
            // memory, before those passes run.
            if (auto* windowed = dynamic_cast<WindowedSiteReadSource*>(read_source.get())) {
                size_t freed = windowed->drop_cached_windows();
                if (show_progress) {
                    logger.info() << "Freed the read cache: " << freed << " reads" << endl;
                }
            }
            // The gaf-base processes hold caches of their own, and no more reads are needed.
            if (auto* gaf_base = dynamic_cast<GafBaseSiteReadSource*>(read_source.get())) {
                gaf_base->stop_processes();
                if (show_progress) {
                    logger.info() << "Stopped the gaf-base processes" << endl;
                }
            }
        });
    } else if (!call_chains) {
        // Call each snarl
        if (show_progress) logger.info() << "Calling top-level snarls" << endl;
        graph_caller->call_top_level_snarls(*graph, recurse_type);
    } else {
        // Attempt to call chains instead of snarls so that the output traversals are longer
        // Todo: this could probably help in some cases when making VCFs too
        if (show_progress) logger.info() << "Calling top-level chains" << endl;
        graph_caller->call_top_level_chains(*graph, max_chain_edges, max_chain_trivial_travs, recurse_type);
    }

    if (show_progress) {
        if (const auto* windowed = dynamic_cast<const WindowedSiteReadSource*>(read_source.get())) {
            log_windowed_read_source_stats(logger, *windowed, off_panel);
        }
    }

    if (show_progress) {
        const auto* graph_calculator =
            dynamic_cast<const GraphAlignedAlleleLikelihoodCalculator*>(likelihood_calculator.get());
        if (graph_calculator != nullptr && graph_calculator->rate_id_fallbacks() > 0) {
            logger.info() << graph_calculator->rate_id_fallbacks() << " sites had no reference "
                          << "position for their depth-rate window and used a node-ID window"
                          << endl;
        }
    }

    if (show_progress) logger.info() << "Calling complete" << endl;

    if (!gaf_output) {
        // Output VCF
        assert(vcf_out != nullptr);
        if (show_progress) logger.info() << "Writing VCF Variants" << endl;
        if (multipass_caller != nullptr) {
            // The read-likelihood caller places each record among the sites it called, which
            // include its edit sites.
            vcf_out->write_variants(
                cout, &multipass_caller->site_decomposition(), graph,
                [&](const net_handle_t& site) {
                    return multipass_caller->record_id_of_site(site);
                });
        } else {
            vcf_out->write_variants(cout, nesting_sites.get(), graph);
        }
        if (show_progress) logger.info() << "VCF complete" << endl;        
    }

    // Everything is written. Close what is still open and exit without freeing what calling built,
    // which would take minutes on a whole genome.
    if (likelihood_dump) {
        likelihood_dump->close();
    }
    // The GAF emitter (-G) holds its last batch of records until it is destroyed.
    alignment_emitter.reset();
    cout.flush();
    cerr.flush();
    exit(EXIT_SUCCESS);
}

// Register subcommand
static Subcommand vg_call("call", "call or genotype VCF variants", PIPELINE, 10, main_call);

