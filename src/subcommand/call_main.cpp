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

#include "subcommand.hpp"
#include "../path.hpp"
#include "../graph_caller.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../xg.hpp"
#include "../gbzgraph.hpp"
#include "../gbwtgraph_helper.hpp"
#include "../gref.hpp"
#include "../traversal_clusters.hpp"
#include "../read_likelihood_caller.hpp"
#include "../site_read_source.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

const string DEFAULT_SAMPLE_NAME = "SAMPLE";

// The default --read-window, in node IDs, for the two read sources that fetch reads by node-ID
// window. The window size can change the order in which a site sees its reads, and so the last
// digits of likelihoods summed over them.

/// A GAM index query is a seek in an open file, so its window is narrow, to avoid fetching reads
/// that no site needs.
const size_t DEFAULT_GAM_INDEX_WINDOW = 256;
/// A GAF-Base query runs a gbz-base subprocess, so its window is wide enough to serve many sites
/// per query, and much wider than a long read's node-ID span, so that few reads are fetched twice
/// by crossing a window boundary.
const size_t DEFAULT_GAF_BASE_WINDOW = 16384;
/// An indexed GAF query is a seek, but its window is as wide as GAF-Base's, so that the two fetch
/// the same windows and few long reads are decoded twice by crossing a window boundary.
const size_t DEFAULT_GAF_INDEX_WINDOW = 16384;

/// Count the haplotypes that a graph's HAPLOTYPE-sense paths belong to. A haplotype is identified
/// by its sample name and haplotype number, so one stored as several paths (one per contig or
/// fragment) counts once. Reference and generic paths are not counted.
///
/// This is an upper bound on the number of distinct non-reference alleles that haplotype-based
/// allele enumeration can offer at a site.
static size_t count_panel_haplotypes(const PathHandleGraph& graph) {
    set<pair<string, size_t>> haplotypes;
    graph.for_each_path_of_sense(PathSense::HAPLOTYPE, [&](const path_handle_t& path) {
        haplotypes.emplace(graph.get_sample_name(path), graph.get_haplotype(path));
    });
    return haplotypes.size();
}

/// One option that a --preset sets, written as on the command line. A switch has an empty value.
struct PresetSetting {
    const char* option;
    const char* value;
};

/// The settings of each --preset. help_call lists them, and main_call applies each one whose option
/// was not given explicitly, so the help cannot list a setting that is not applied.
static const vector<pair<string, vector<PresetSetting>>> PRESETS = {
    {"ont", {{"--read-phasing", ""}, {"--regenotype", ""}, {"--gap-open", "1"},
             {"--gap-extend", "1"}, {"--mismap-min", "0.05"}, {"--read-min-mapq", "5"},
             {"--insertion-nats", "0.9"}, {"--hp-prior", "20"}}},
};

/// Each preset's settings for the --preset help, wrapped at the help's description column.
static string preset_help_lines() {
    const string indent(28, ' ');
    const size_t width = 80;
    stringstream out;
    for (const auto& preset : PRESETS) {
        string line = indent + preset.first + ":";
        for (const PresetSetting& setting : preset.second) {
            string item = setting.option;
            if (*setting.value != '\0') {
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
         << "                            GAF-Base database or gbz-base. --gaf-reads must be" << endl
         << "                            sorted (vg gamsort -G), compressed with bgzip and" << endl
         << "                            indexed with tabix -p gaf" << endl
         << "      --gaf-base FILE       GAF-Base database of read alignments, queried as" << endl
         << "                            sites need them by running gbz-base, which must be" << endl
         << "                            on the PATH" << endl
         << "      --gbz-base FILE       graph for --gaf-base queries, as GBZ-Base or GBZ" << endl
         << "                            [the input graph]" << endl
         << "      --gaf-base-binary P   gbz-base executable to run [gbz-base]" << endl
         << "      --read-window N       node-ID window for indexed read fetches" << endl
         << "                            [16384 for --gaf-base and --gaf-index, 256 for" << endl
         << "                            --gam-index]" << endl
         << "      --read-min-mapq N     ignore reads with MAPQ below N [0]" << endl
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
         << "      --hp-prior F          exponent used instead of --linkage-prior where an" << endl
         << "                            allele differs from the reference only in the" << endl
         << "                            length of a homopolymer run; 0 disables it [0]" << endl
         << "      --hp-prior-run N      shortest homopolymer run for --hp-prior [11]" << endl
         << "      --no-phased           write unphased genotypes, without FORMAT/PS" << endl
         << "      --phased              write phased genotypes, failing if the linkage" << endl
         << "                            model cannot run [on when the linkage model runs]" << endl
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
         << "  -d, --ploidy N            ploidy of sample. {1, 2} [2]" << endl
         << "      --no-nested           report variation inside nested snarls as part of" << endl
         << "                            the enclosing snarl's alleles" << endl
         << "      --nested              genotype nested snarls in records of their own," << endl
         << "                            not as part of the enclosing snarl's alleles" << endl
         << "                            [on with --read-likelihood]" << endl
         << "      --atomize-blocks      write one record per separate difference between" << endl
         << "                            a called allele and the reference, rather than one" << endl
         << "                            per snarl [on with --read-likelihood]" << endl
         << "      --no-atomize-blocks   write one record per snarl" << endl
         << "      --ploidy-bed FILE     BED of CHROM START END PLOIDY giving the ploidy of" << endl
         << "                            each region, overriding -d and -R; CHROM is the" << endl
         << "                            VCF contig name, and intervals must not overlap" << endl
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
    bool   gbz_paths_explicit = false;
    bool   enumerate_support = false;
    // Phasing comes from the linkage model; see the check where that model is set up.
    bool   phased_output = true;
    bool   phased_explicit = false;
    string mosaic_out;
    string anchors_out;
    bool no_off_ref_nesting = false;
    AnchorParams anchor_params;
    /// Counts of how the anchor code placed reads in this run.
    AnchorCounters anchor_run_counters;
    // Mosaic settings. Where no haplotype can be followed, fill the gap with the reference.
    bool mosaic_patch_gaps = true;
    // Include nested sites; if false, each strand follows its enclosing site's haplotype through
    // nested snarls instead of switching haplotypes inside them.
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

    bool traversals_only = false;
    bool gaf_output = false;
    size_t trav_padding = 0;
    bool genotype_snarls = false;
    bool top_down = false;
    bool bottom_up = false;
    // Both on by default, and turned off below where they cannot apply.
    bool atomize_blocks = true;
    bool atomize_explicit = false;
    bool nested_calling = true;
    bool nested_explicit = false;
    bool call_chains = false;
    bool all_snarls = false;
    size_t min_allele_len = 0;
    size_t max_allele_len = numeric_limits<size_t>::max();
    bool show_progress = false;

    // Nested calling option (for use with --top-down)
    bool star_allele = false;

    // -L: after genotyping, merge called ALT alleles whose similarity is at least this.
    double cluster_threshold = 1.0;
    // --cluster-min-len: the similarity counts sequence that all alleles share, so at a short site
    // almost any two alleles look alike, so we merge only at sites at least this long. The default
    // is the usual structural-variant size.
    int64_t cluster_min_allele_len = 50;
    bool cluster_min_len_set = false;

    // Options for the read-likelihood genotyper (--read-likelihood).
    bool read_likelihood = false;
    string gam_filename;
    string gaf_filename;
    string dump_likelihoods_filename;
    string gam_index_filename;
    string gaf_index_filename;
    string gaf_base_filename;
    string gbz_base_filename;
    string gaf_base_binary = "gbz-base";
    // 0 means the read source's default; see the window defaults at the top of the file.
    size_t read_window_size = 0;
    bool no_mismap_term = false;
    bool no_share_quality = false;
    double depth_quality = 0.0;
    // Gap penalties for scoring reads against alleles.
    int gap_open = default_gap_open;
    int gap_extend = default_gap_extension;
    // Phase heterozygous sites with the reads, as well as with the haplotypes.
    bool read_phasing = false;
    bool read_phasing_explicit = false;
    bool regenotype = false;
    bool regenotype_explicit = false;
    bool phase_hets_explicit = false;
    RegenotypeParams regenotype_params;
    // Linkage passes under --regenotype, one per round. The first chooses the genotypes; each
    // later one follows a correction from the read phase and chooses them again.
    size_t regenotype_passes = 2;
    string regenotype_ledger;
    ReadPhasingParams read_phasing_params;
    // The preset is applied after all options are parsed. The *_explicit flags record which of its
    // settings were given explicitly, and those keep their given values.
    string preset;
    double insertion_gap_nats = 0.0;
    bool insertion_nats_explicit = false;
    bool optimal_pairing = false;
    bool optimal_pairing_explicit = false;
    bool phase_min_q_explicit = false;
    bool gap_open_explicit = false, gap_extend_explicit = false, mismap_min_explicit = false;
    bool read_min_mapq_explicit = false;
    double min_confidence = 0.0;
    double linkage_weight = 2.0;
    /// Whether --linkage-weight was given. Where the linkage model cannot run, an explicit weight
    /// is an error, while the default weight is set to 0.
    bool linkage_weight_explicit = false;
    double linkage_scale = 10000.0;
    double linkage_freq_prior = 5.0;
    double hp_prior = 0.0;
    bool hp_prior_explicit = false;
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
    bool max_snarl_edges_explicit = false;
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
    constexpr int OPT_GBZ_BASE = 1016;
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
    int c;
    optind = 2; // force optind past command positional argument
    // The long options, grouped by the subsystem they configure. The options of a subsystem other
    // than "core" only have an effect when that subsystem is turned on, and we refuse them otherwise
    // (see below). Every such subsystem is part of the read-likelihood genotyper.
    const map<string, vector<struct option>> long_options_by_subsystem = {
        {"core", {
            {"pack", required_argument, 0, 'k'},
            {"bias-mode", no_argument, 0, 'B'},
            {"baseline-error", required_argument, 0, 'e'},
            {"het-bias", required_argument, 0, 'b'},
            {"min-support", required_argument, 0, 'm'},
            {"vcf", required_argument, 0, 'v'},
            {"genotype-snarls", no_argument, 0, 'a'},
            {"all-snarls", no_argument, 0, 'A'},
            {"min-length", required_argument, 0, 'c'},
            {"max-length", required_argument, 0, 'C'},
            {"ref-fasta", required_argument, 0, 'f'},
            {"ins-fasta", required_argument, 0, 'i'},
            {"sample", required_argument, 0, 's'},
            {"snarls", required_argument, 0, 'r'},
            {"gbwt", required_argument, 0, 'g'},
            {"gbz", no_argument, 0, 'z'},
            {"translation", required_argument, 0, 'N'},
            {"gbz-translation", no_argument, 0, 'O'},
            {"ref-path", required_argument, 0, 'p'},
            {"path-prefix", required_argument, 0, 'P'},
            {"ref-sample", required_argument, 0, 'S'},
            {"ref-offset", required_argument, 0, 'o'},
            {"ref-length", required_argument, 0, 'l'},
            {"ploidy", required_argument, 0, 'd'},
            {"ploidy-regex", required_argument, 0, 'R'},
            {"ploidy-bed", required_argument, 0, OPT_PLOIDY_BED},
            {"nested", no_argument, 0, OPT_NESTED},
            {"no-nested", no_argument, 0, OPT_NO_NESTED},
            {"gaf", no_argument, 0, 'G'},
            {"traversals", no_argument, 0, 'T'},
            {"trav-padding", required_argument, 0, 'M'},
            {"legacy", no_argument, 0, OPT_LEGACY},
            {"top-down", no_argument, 0, OPT_TOP_DOWN},
            {"bottom-up", no_argument, 0, OPT_BOTTOM_UP},
            {"atomize-blocks", no_argument, 0, OPT_ATOMIZE_BLOCKS},
            {"no-atomize-blocks", no_argument, 0, OPT_NO_ATOMIZE_BLOCKS},
            {"read-likelihood", no_argument, 0, OPT_READ_LIKELIHOOD},
            {"max-snarl-edges", required_argument, 0, OPT_MAX_SNARL_EDGES},
            {"chains", no_argument, 0, 'I'},
            {"cluster", required_argument, 0, 'L'},
            {"cluster-min-len", required_argument, 0, OPT_CLUSTER_MIN_LEN},
            // Accepted and ignored, and left out of the help, so that command lines that still
            // pass it keep working.
            {"cluster-post", no_argument, 0, OPT_CLUSTER_POST},
            {"star-allele", no_argument, 0, 'Y'},
            {"threads", required_argument, 0, 't'},
            {"progress", no_argument, 0, OPT_PROGRESS},
            {"help", no_argument, 0, 'h'},
        }},
        {"read-likelihood", {  // need --read-likelihood
            {"no-phased", no_argument, 0, OPT_NO_PHASED},
            {"gam", required_argument, 0, OPT_GAM},
            {"gaf-reads", required_argument, 0, OPT_GAF},
            {"dump-likelihoods", required_argument, 0, OPT_DUMP_LIKELIHOODS},
            {"no-mismap-term", no_argument, 0, OPT_NO_MISMAP_TERM},
            {"mismap-max", required_argument, 0, OPT_MISMAP_MAX},
            {"mismap-min", required_argument, 0, OPT_MISMAP_MIN},
            {"insertion-nats", required_argument, 0, OPT_INSERTION_GAP_NATS},
            {"optimal-pairing", no_argument, 0, OPT_OPTIMAL_PAIRING},
            {"no-optimal-pairing", no_argument, 0, OPT_NO_OPTIMAL_PAIRING},
            {"no-share-quality", no_argument, 0, OPT_NO_SHARE_QUALITY},
            {"flat-mixture", no_argument, 0, OPT_FLAT_MIXTURE},
            {"depth-term", required_argument, 0, OPT_DEPTH_TERM},
            {"depth-count-raw", no_argument, 0, OPT_DEPTH_COUNT_RAW},
            {"depth-quality", required_argument, 0, OPT_DEPTH_QUALITY},
            {"preset", required_argument, 0, OPT_PRESET},
            {"gap-open", required_argument, 0, OPT_GAP_OPEN},
            {"gap-extend", required_argument, 0, OPT_GAP_EXTEND},
            {"read-phasing", no_argument, 0, OPT_READ_PHASING},
            {"no-read-phasing", no_argument, 0, OPT_NO_READ_PHASING},
            {"phase-min-q", required_argument, 0, OPT_PHASE_MIN_Q},
            {"phase-break", required_argument, 0, OPT_PHASE_BREAK},
            {"phase-relink", required_argument, 0, OPT_PHASE_RELINK},
            {"phase-hang", required_argument, 0, OPT_PHASE_HANG},
            {"phase-prior", required_argument, 0, OPT_PHASE_PRIOR},
            {"phase-cap", required_argument, 0, OPT_PHASE_CAP},
            {"phase-coherence", required_argument, 0, OPT_PHASE_COHERENCE},
            {"phase-coh-rounds", required_argument, 0, OPT_PHASE_COH_ROUNDS},
            {"regenotype", no_argument, 0, OPT_REGENOTYPE},
            {"no-regenotype", no_argument, 0, OPT_NO_REGENOTYPE},
            {"min-confidence", required_argument, 0, OPT_MIN_CONFIDENCE},
            {"linkage-weight", required_argument, 0, OPT_LINKAGE_WEIGHT},
            {"linkage-scale", required_argument, 0, OPT_LINKAGE_SCALE},
            {"linkage-prior", required_argument, 0, OPT_LINKAGE_FREQ_PRIOR},
            {"hp-prior", required_argument, 0, OPT_HP_PRIOR},
            {"hp-prior-run", required_argument, 0, OPT_HP_PRIOR_RUN},
            {"enumerate-support", no_argument, 0, OPT_ENUMERATE_SUPPORT},
            {"phased", no_argument, 0, OPT_PHASED},
            {"mosaic-out", required_argument, 0, OPT_MOSAIC_OUT},
            {"anchors-out", required_argument, 0, OPT_ANCHORS_OUT},
            {"read-min-mapq", required_argument, 0, OPT_READ_MIN_MAPQ},
            {"gam-index", required_argument, 0, OPT_GAM_INDEX},
            {"gaf-index", required_argument, 0, OPT_GAF_INDEX},
            {"gaf-base", required_argument, 0, OPT_GAF_BASE},
            {"gbz-base", required_argument, 0, OPT_GBZ_BASE},
            {"gaf-base-binary", required_argument, 0, OPT_GAF_BASE_BINARY},
            {"read-window", required_argument, 0, OPT_READ_WINDOW},
        }},
        {"anchors", {  // need --anchors-out
            {"no-off-ref-nesting", no_argument, 0, OPT_NO_OFF_REF_NESTING},
            {"anchors-hom-split", no_argument, 0, OPT_ANCHORS_HOM_SPLIT},
            {"anchors-phase-hets", no_argument, 0, OPT_ANCHORS_PHASE_HETS},
            {"no-anchors-phase-hets", no_argument, 0, OPT_NO_ANCHORS_PHASE_HETS},
            {"anchors-strict-hets", no_argument, 0, OPT_ANCHORS_STRICT_HETS},
            {"split-min-q", required_argument, 0, OPT_ANCHORS_PHASE_MIN},
            {"split-min-side", required_argument, 0, OPT_ANCHORS_PHASE_MIN_SIDE},
            {"anchors-het-only", no_argument, 0, OPT_ANCHORS_HET_ONLY},
            {"anchors-leaf-only", no_argument, 0, OPT_ANCHORS_LEAF_ONLY},
            {"anchors-reads", required_argument, 0, OPT_ANCHORS_MIN_READS},
            {"anchors-min-gqn", required_argument, 0, OPT_ANCHORS_MIN_GQN},
            {"anchors-min-q", required_argument, 0, OPT_ANCHORS_MIN_SCORE},
            {"anchors-keep-off-call", no_argument, 0, OPT_ANCHORS_KEEP_OFF_CALL},
            {"anchors-end-new", required_argument, 0, OPT_ANCHORS_END_PIN_MIN_NEW},
        }},
        {"mosaic", {  // need --mosaic-out
            {"mosaic-patch-gaps", no_argument, 0, OPT_MOSAIC_PATCH},
            {"no-mosaic-patch-gaps", no_argument, 0, OPT_NO_MOSAIC_PATCH},
            {"no-mosaic-nested", no_argument, 0, OPT_NO_MOSAIC_NESTED},
            {"mosaic-break-unexplained", no_argument, 0, OPT_MOSAIC_BREAK_UNEXPL},
        }},
        {"regenotype", {  // need --regenotype
            {"regeno-ceiling", required_argument, 0, OPT_REGENO_CEILING},
            {"regeno-haploid", no_argument, 0, OPT_REGENO_HAPLOID},
            {"no-regeno-haploid", no_argument, 0, OPT_NO_REGENO_HAPLOID},
            {"regeno-temper", required_argument, 0, OPT_REGENO_TEMPER},
            {"regeno-passes", required_argument, 0, OPT_REGENO_PASSES},
            {"regeno-shuffle", no_argument, 0, OPT_REGENO_SHUFFLE},
            {"regeno-ledger", required_argument, 0, OPT_REGENO_LEDGER},
        }},
    };

    // getopt_long takes a single array ending in an all-zero entry. We also record each option's
    // name and subsystem by its `val`, which is unique.
    vector<struct option> long_options;
    unordered_map<int, string> option_name;
    unordered_map<int, string> option_subsystem;
    for (const auto& subsystem_and_options : long_options_by_subsystem) {
        for (const struct option& o : subsystem_and_options.second) {
            long_options.push_back(o);
            option_name[o.val] = o.name;
            option_subsystem[o.val] = subsystem_and_options.first;
        }
    }
    long_options.push_back({0, 0, 0, 0});

    // The options given, by `val`, in the order first given, for the subsystem check below.
    vector<int> options_seen;

    while (true) {


        int option_index = 0;

        c = getopt_long (argc, argv, "k:Be:b:m:v:aAc:C:f:i:s:r:g:zN:Op:P:S:o:l:d:R:GTM:IL:Yt:h?",
                         long_options.data(), &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        if (std::find(options_seen.begin(), options_seen.end(), c) == options_seen.end()) {
            options_seen.push_back(c);
        }

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
            gbz_paths_explicit = true;
            break;
        case OPT_ENUMERATE_SUPPORT:
            enumerate_support = true;
            break;
        case OPT_PHASED:
            phased_output = true;
            phased_explicit = true;
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
            atomize_explicit = true;
            break;
        case OPT_NO_ATOMIZE_BLOCKS:
            atomize_blocks = false;
            atomize_explicit = true;
            break;
        case OPT_BOTTOM_UP:
            bottom_up = true;
            break;
        case OPT_READ_LIKELIHOOD:
            read_likelihood = true;
            break;
        case OPT_GAM:
            gam_filename = require_exists(logger, optarg);
            break;
        case OPT_GAF:
            gaf_filename = require_exists(logger, optarg);
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
            max_snarl_edges_explicit = true;
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
            gap_open_explicit = true;
            gap_open = parse<int>(optarg);
            if (gap_open < 1 || gap_open > 127) {
                cerr << "error [vg call]: --gap-open must be between 1 and 127" << endl;
                return 1;
            }
            break;
        case OPT_READ_PHASING:
            read_phasing = true;
            read_phasing_explicit = true;
            break;
        case OPT_NO_READ_PHASING:
            read_phasing = false;
            read_phasing_explicit = true;
            break;
        case OPT_PHASE_MIN_Q:
            read_phasing_params.reliability = parse<double>(optarg);
            phase_min_q_explicit = true;
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
                cerr << "error [vg call]: --phase-coherence is a fraction in [0,1]" << endl;
                return 1;
            }
            break;
        case OPT_PHASE_COH_ROUNDS:
            read_phasing_params.coherence_rounds = parse<size_t>(optarg);
            if (read_phasing_params.coherence_rounds < 1) {
                cerr << "error [vg call]: --phase-coh-rounds must be >= 1" << endl;
                return 1;
            }
            break;
        case OPT_REGENOTYPE:
            regenotype = true;
            regenotype_explicit = true;
            break;
        case OPT_NO_REGENOTYPE:
            regenotype = false;
            regenotype_explicit = true;
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
        case OPT_GAP_EXTEND:
            gap_extend_explicit = true;
            gap_extend = parse<int>(optarg);
            if (gap_extend < 1 || gap_extend > 127) {
                cerr << "error [vg call]: --gap-extend must be between 1 and 127" << endl;
                return 1;
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
        case OPT_NESTED:
            nested_calling = true;
            nested_explicit = true;
            break;
        case OPT_NO_OFF_REF_NESTING:
            no_off_ref_nesting = true;
            break;
        case OPT_NO_NESTED:
            nested_calling = false;
            nested_explicit = true;

            break;
        case OPT_NO_PHASED:
            phased_output = false;
            phased_explicit = true;
            break;
        case OPT_LINKAGE_WEIGHT:
            linkage_weight = parse<double>(optarg);
            linkage_weight_explicit = true;
            break;
        case OPT_LINKAGE_SCALE:
            linkage_scale = parse<double>(optarg);
            break;
        case OPT_LINKAGE_FREQ_PRIOR:
            linkage_freq_prior = parse<double>(optarg);
            break;
        case OPT_HP_PRIOR:
            hp_prior = parse<double>(optarg);
            hp_prior_explicit = true;
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
            mismap_min_explicit = true;
            min_mismap_prob = parse<double>(optarg);
            break;
        case OPT_INSERTION_GAP_NATS:
            insertion_nats_explicit = true;
            insertion_gap_nats = parse<double>(optarg);
            break;
        case OPT_OPTIMAL_PAIRING:
            optimal_pairing_explicit = true;
            optimal_pairing = true;
            break;
        case OPT_ANCHORS_HOM_SPLIT:
            anchor_params.hom_split = true;
            break;
        case OPT_ANCHORS_PHASE_HETS:
            anchor_params.phase_hets = true;
            phase_hets_explicit = true;
            break;
        case OPT_NO_ANCHORS_PHASE_HETS:
            anchor_params.phase_hets = false;
            phase_hets_explicit = true;
            break;
        case OPT_ANCHORS_STRICT_HETS:
            anchor_params.strict_hets = true;
            phase_hets_explicit = true;
            break;
        case OPT_ANCHORS_PHASE_MIN:
            anchor_params.phase_min = parse<double>(optarg);
            if (anchor_params.phase_min < 0) {
                cerr << "error [vg call]: --split-min-q must be at least 0" << endl;
                return 1;
            }
            break;
        case OPT_ANCHORS_PHASE_MIN_SIDE:
            anchor_params.phase_min_side = parse<size_t>(optarg);
            if (anchor_params.phase_min_side < 1) {
                cerr << "error [vg call]: --split-min-side must be at least 1" << endl;
                return 1;
            }
            break;
        case OPT_NO_OPTIMAL_PAIRING:
            optimal_pairing_explicit = true;
            optimal_pairing = false;
            break;
        case OPT_READ_MIN_MAPQ:
            read_min_mapq_explicit = true;
            read_min_mapq = parse<int>(optarg);
            break;
        case OPT_GAM_INDEX:
            gam_index_filename = require_exists(logger, optarg);
            break;
        case OPT_GAF_INDEX:
            gaf_index_filename = require_exists(logger, optarg);
            break;
        case OPT_GAF_BASE:
            gaf_base_filename = require_exists(logger, optarg);
            break;
        case OPT_GBZ_BASE:
            gbz_base_filename = require_exists(logger, optarg);
            break;
        case OPT_GAF_BASE_BINARY:
            gaf_base_binary = optarg;
            break;
        case OPT_READ_WINDOW:
            read_window_size = parse<size_t>(optarg);
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

    if (hp_prior < 0.0 || hp_prior_run < 1) {
        cerr << "error [vg call]: --hp-prior takes a value >= 0, and --hp-prior-run a run of at least 1"
             << endl;
        return 1;
    }
    if (!preset.empty()) {
        // A preset only sets the options that were not given explicitly. The values for `ont`
        // were fitted to Oxford Nanopore reads.
        auto chosen = std::find_if(PRESETS.begin(), PRESETS.end(),
                                   [&](const pair<string, vector<PresetSetting>>& p) {
                                       return p.first == preset;
                                   });
        if (chosen == PRESETS.end()) {
            cerr << "error [vg call]: unknown --preset '" << preset << "'; known presets:";
            for (size_t i = 0; i < PRESETS.size(); ++i) {
                cerr << (i ? ", " : " ") << PRESETS[i].first;
            }
            cerr << endl;
            return 1;
        }
        for (const PresetSetting& setting : chosen->second) {
            const string option = setting.option;
            const string value = setting.value;
            if (option == "--read-phasing") {
                if (!read_phasing_explicit) {
                    read_phasing = true;
                }
            } else if (option == "--regenotype") {
                if (!regenotype_explicit) {
                    regenotype = true;
                }
            } else if (option == "--gap-open") {
                if (!gap_open_explicit) {
                    gap_open = parse<int>(value);
                }
            } else if (option == "--gap-extend") {
                if (!gap_extend_explicit) {
                    gap_extend = parse<int>(value);
                }
            } else if (option == "--mismap-min") {
                if (!mismap_min_explicit) {
                    min_mismap_prob = parse<double>(value);
                }
            } else if (option == "--read-min-mapq") {
                if (!read_min_mapq_explicit) {
                    read_min_mapq = parse<int>(value);
                }
            } else if (option == "--insertion-nats") {
                if (!insertion_nats_explicit) {
                    insertion_gap_nats = parse<double>(value);
                }
            } else if (option == "--hp-prior") {
                if (!hp_prior_explicit) {
                    hp_prior = parse<double>(value);
                }
            } else {
                // Every option in PRESETS needs a branch here.
                logger.error() << "--preset " << preset << " sets " << option
                               << ", which has no preset handling" << endl;
            }
        }
    }
    // --atomize-blocks is on by default, so where it cannot apply we turn it off, and only an
    // explicit --atomize-blocks is an error. These checks depend only on the options, so they run
    // before the graph is loaded.
    if (atomize_blocks && genotype_snarls) {
        // -a writes the same records, one per snarl, whatever the sample, but block emission
        // writes one record per difference in the called haplotypes.
        if (atomize_explicit) {
            logger.error() << "--atomize-blocks cannot be combined with -a/--genotype-snarls: "
                           << "a block list depends on the called haplotypes, so the record set "
                           << "would stop being sample-independent" << endl;
        }
        atomize_blocks = false;
    }
    if (atomize_blocks && (legacy || bottom_up || top_down)) {
        if (atomize_explicit) {
            logger.error() << "--atomize-blocks needs the default calling path "
                           << "(not --legacy, --bottom-up or --top-down)" << endl;
        }
        atomize_blocks = false;
    }
    // -I calls each piece of a chain as if it were a snarl. Only a piece holding a single snarl
    // matches a snarl in the snarl manager, and nested calling, on which block emission depends,
    // only works at those, so it would apply to only some sites.
    if (atomize_blocks && call_chains) {
        if (atomize_explicit) {
            logger.error() << "--atomize-blocks cannot be combined with -I/--chains, which calls "
                           << "pieces of chains; blocks would be written only at pieces that hold "
                           << "a single snarl" << endl;
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
    // -L merges alleles as VCFOutputCaller::emit_variant writes each record, and neither the VCF
    // genotyper (-v) nor GAF output writes records through it.
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

    // Reported ahead of the subsystem check below, since it is an error whatever else is given.
    if (!gam_index_filename.empty() && gam_filename.empty()) {
        logger.error() << "--gam-index requires --gam" << endl;
    }
    if (!gaf_index_filename.empty() && gaf_filename.empty()) {
        logger.error() << "--gaf-index requires --gaf-reads" << endl;
    }

    // Refuse options for a subsystem that is not turned on, since they would have no effect.
    {
        // Each check names the option that turns some subsystems on, and says whether they are
        // on. The read-likelihood genotyper contains all the other subsystems, so its check
        // covers them all.
        struct SubsystemCheck {
            const char* needs;
            const char* preposition;    // joins the error message to `needs`
            bool on;
            set<string> subsystems;
        };
        const SubsystemCheck checks[] = {
            {"--read-likelihood", "to", read_likelihood,
             {"read-likelihood", "anchors", "mosaic", "regenotype"}},
            {"--anchors-out", "with", !anchors_out.empty(), {"anchors"}},
            {"--mosaic-out", "with", !mosaic_out.empty(), {"mosaic"}},
            // A preset can turn --regenotype on and --no-read-phasing then turns it off again
            // (below), so we check the value that will take effect.
            {"--regenotype", "with", regenotype && read_phasing, {"regenotype"}},
        };
        for (const SubsystemCheck& check : checks) {
            if (check.on) {
                continue;
            }
            vector<string> offenders;
            for (int val : options_seen) {
                auto subsystem = option_subsystem.find(val);
                if (subsystem != option_subsystem.end() && check.subsystems.count(subsystem->second)) {
                    offenders.push_back("--" + option_name.at(val));
                }
            }
            if (offenders.empty()) {
                continue;
            }
            stringstream joined;
            for (size_t i = 0; i < offenders.size(); ++i) {
                joined << (i ? ", " : "") << offenders[i];
            }
            // logger.error() exits, so only the first failing check reports. We check
            // --read-likelihood first, because the user has to add it before any of the other
            // switches can take effect.
            logger.error() << joined.str()
                           << (offenders.size() == 1 ? " only applies " : " only apply ")
                           << check.preposition << " " << check.needs
                           << ", which was not given" << endl;
        }
    }

    // Read phasing, re-genotyping and hom splitting all start from the phase the linkage model
    // gives each site. Where there is no such phase, each is an error if given explicitly, and is
    // turned off if a preset turned it on. `why` ends the error's sentence, saying why there is
    // no phase.
    auto refuse_phase_dependents = [&](const string& why) {
        vector<string> offenders;
        if (read_phasing && read_phasing_explicit) {
            offenders.push_back("--read-phasing");
        }
        if (regenotype && regenotype_explicit) {
            offenders.push_back("--regenotype");
        }
        if (anchor_params.hom_split) {
            offenders.push_back("--anchors-hom-split");
        }
        if (!offenders.empty()) {
            stringstream joined;
            for (size_t i = 0; i < offenders.size(); ++i) {
                joined << (i ? ", " : "") << offenders[i];
            }
            logger.error() << joined.str() << (offenders.size() == 1 ? " starts" : " start")
                           << " from the phase the linkage model gives each site, " << why
                           << endl;
        }
        read_phasing = false;
        regenotype = false;
    };
    if (phased_explicit && !phased_output) {
        refuse_phase_dependents("which --no-phased turns off");
    }

    // --read-likelihood needs exactly one read source, and cannot be combined with the ratio or
    // legacy support callers.
    if (read_likelihood) {
        int read_source_count = (gam_filename.empty() ? 0 : 1) + (gaf_filename.empty() ? 0 : 1) +
                                (gaf_base_filename.empty() ? 0 : 1);
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
        if (legacy) {
            logger.error() << "--read-likelihood cannot be used with --legacy" << endl;
        }
    }

    // Under --read-likelihood on a GBZ, take candidate alleles from the GBZ's haplotypes by default,
    // as -z does, rather than from read support. This needs no pack file, but it can only offer
    // alleles that some haplotype carries; --enumerate-support turns it off.
    if (read_likelihood && !gbz_paths && !enumerate_support && gbz_graph &&
        gbwt_filename.empty() && vcf_filename.empty()) {
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
    if (enumerate_support && gbz_paths_explicit) {
        logger.error() << "--enumerate-support and -z/--gbz select different allele "
                       << "enumeration strategies: choose one or the other" << endl;
    }
    if (enumerate_support && !gbwt_filename.empty()) {
        logger.error() << "--enumerate-support and -g/--gbwt select different allele "
                       << "enumeration strategies: choose one or the other" << endl;
    }

    // The linkage model needs --read-likelihood, haplotypes to enumerate alleles from (-z or -g)
    // and a positive --linkage-weight. This is checked before the likelihood calculator is built,
    // so that it does not keep read-phasing evidence for nothing. A GBWT with too few haplotypes
    // is found only once it is loaded, below.
    const bool haplotypes_given = gbwt_index != nullptr || !gbwt_filename.empty();
    if (!(read_likelihood && haplotypes_given && linkage_weight > 0.0)) {
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
    if (optimal_pairing && !phase_min_q_explicit) {
        read_phasing_params.reliability = 8.5;
    }
    // Under read phasing, refuse a --phase-min-q above the heterozygous score ceiling,
    // phred(e / (e + (1 - e) / 2)) for e = --mismap-min: the confidence of a read that fits one of
    // two equal-length alleles perfectly and the other not at all. Sites whose alleles are of
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

    // --gbz-base only configures --gaf-base queries.
    if (!gbz_base_filename.empty() && gaf_base_filename.empty()) {
        logger.error() << "--gbz-base requires --gaf-base" << endl;
    }

    // Validation: -A, --top-down, and --bottom-up are mutually exclusive
    int nested_mode_count = (all_snarls ? 1 : 0) + (top_down ? 1 : 0) + (bottom_up ? 1 : 0);
    if (nested_mode_count > 1) {
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

    // Prefix given: find all paths matching it
    if (!ref_path_prefixes.empty()) {
        graph->for_each_path_of_sense({PathSense::REFERENCE, PathSense::GENERIC, PathSense::HAPLOTYPE}, [&](const path_handle_t& path_handle) {
            string path_name = graph->get_path_name(path_handle);
            // Never include alt paths in reference paths
            if (Paths::is_alt(path_name)) {
                return;
            }
            for (auto& prefix : ref_path_prefixes) {
                if (path_name.compare(0, prefix.size(), prefix) == 0) {
                    ref_paths.push_back(path_name);
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
                ref_paths.push_back(name);
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

    // Whether to genotype off-reference chains: nested chains that the reference paths do not
    // pass through. Their variants have no position on the reference, so they only get VCF records
    // when a gref fragment path (see gref.hpp) gives them a contig of their own. Selecting such a
    // path as a reference turns this on, as does --anchors-out (below). VG_CALL_NO_REF_NESTED
    // turns it on without either, for testing.
    bool off_ref_nesting = getenv("VG_CALL_NO_REF_NESTED") != nullptr;

    for (const string& ref_path : ref_paths) {
        // Only a fragment counts. The gref copy of a base contig is the reference under another
        // name, and gives no chain a contig of its own.
        if (GrefCover::is_gref_name(ref_path)) {
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
    if (!anchors_out.empty() && !no_off_ref_nesting) {
        off_ref_nesting = true;
        if (show_progress) {
            logger.info() << "anchors requested: descending into chains the reference does not "
                          << "cross, which have no VCF record but do have anchors" << endl;
        }
    }

    // For INFO/CH: how many levels of non-reference sequence each selected gref contig lies in,
    // keyed by the contig name the VCF uses. A base contig is level 0, a fragment attached to the
    // base reference is level 1, a fragment attached to a level-1 fragment is level 2, and so on.
    // The cover's paths share no nodes, so a fragment's parent is the gref path that owns a node
    // next to one of the fragment's ends.
    map<string, int> gref_levels;
    if (off_ref_nesting && graph != nullptr) {
        auto gref_owner = [&](handle_t h) {
            string owner;
            graph->for_each_step_on_handle(h, [&](const step_handle_t& step) {
                const string n = graph->get_path_name(graph->get_path_handle_of_step(step));
                if (GrefCover::is_gref_derived(n)) {
                    owner = n;
                    return false;
                }
                return true;
            });
            return owner;
        };
        map<string, int> by_path;
        unordered_set<string> in_progress;
        std::function<int(const string&)> level_of = [&](const string& path_name) -> int {
            if (!GrefCover::is_gref_name(path_name)) {
                return 0;   // a base contig, or the gref copy of one
            }
            auto seen = by_path.find(path_name);
            if (seen != by_path.end()) {
                return seen->second;
            }
            if (!in_progress.insert(path_name).second) {
                return 1;   // a cover is a tree, so this cannot happen; do not hang if it does
            }
            int level = 0;
            if (graph->has_path(path_name)) {
                const path_handle_t p = graph->get_path_handle(path_name);
                if (!graph->is_empty(p)) {
                    // Either end will do: a fragment is a maximal run of uncovered nodes, so both
                    // of its neighbours are boundary nodes of the same enclosing snarl.
                    for (bool left : {true, false}) {
                        const handle_t end = graph->get_handle_of_step(
                            left ? graph->path_begin(p) : graph->path_back(p));
                        graph->follow_edges(end, left, [&](const handle_t& next) {
                            const string owner = gref_owner(next);
                            if (!owner.empty() && owner != path_name) {
                                level = level_of(owner) + 1;
                                return false;
                            }
                            return true;
                        });
                        if (level > 0) {
                            break;
                        }
                    }
                }
            }
            if (level == 0) {
                // No neighbour is on another gref path, as when the neighbouring fragment was
                // shorter than the minimum fragment length and was not written. It is still a
                // fragment, so it is at least level 1.
                level = 1;
            }
            in_progress.erase(path_name);
            by_path[path_name] = level;
            return level;
        };
        for (const string& ref_path : ref_paths) {
            const int level = level_of(ref_path);
            if (level > 0) {
                const string locus = PathMetadata::parse_locus_name(ref_path);
                gref_levels[locus != PathMetadata::NO_LOCUS_NAME ? locus : ref_path] = level;
            }
        }
        if (show_progress && !gref_levels.empty()) {
            map<int, size_t> hist;
            for (const auto& kv : gref_levels) {
                ++hist[kv.second];
            }
            auto& l = logger.info() << "gref nesting levels:";
            for (const auto& kv : hist) {
                l << " " << kv.first << ":" << kv.second;
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
        // the callers only implement ploidy 1 and 2 (see the assert in
        // PoissonSupportSnarlCaller::genotype), and -d is checked for this below.  An unchecked -R
        // reaches the caller as an unsupported ploidy and aborts there instead.
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
            // Skip altpaths (they shouldn't influence snarl decomposition)
            if (GrefCover::is_gref_name(refpath_name)) {
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
    // The read-likelihood genotyper can run without a pack file if allele enumeration does not
    // need support either: GBWTTraversalFinder enumerates from haplotypes and needs none, while
    // FlowTraversalFinder works from node and edge support.
    bool gbwt_enumeration = !gbwt_filename.empty() || gbz_paths;
    bool support_free = read_likelihood && pack_filename.empty();

    if (!pack_filename.empty() || support_free) {
        if (support_free) {
            // Nothing downstream consults support, so a finder that reports none stands in
            // for a pack file.
            support_finder.reset(new NullTraversalSupportFinder(*graph, *snarl_manager));
        } else {
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
                
        // need to use average support when genotyping as small differences in between sample and graph
        // will lead to spots with 0-support, espeically in and around SVs. 
        support_finder->set_support_switch_threshold(avg_trav_threshold, avg_node_threshold);

        // upweight breakpoint edges even when taking average support otherwise
        support_finder->set_min_bp_edge_override(expect_bp_edges);

        // todo: toggle between min / average (or thresholds) via command line
        
        SupportBasedSnarlCaller* packed_caller = nullptr;

        if (read_likelihood) {
            // Read-likelihood genotyping. Support from a pack file, if given, is only used to
            // enumerate candidate alleles.
            if (show_progress) logger.info() << "Loading reads for read-level genotyping" << endl;

            SiteReadFilter read_filter;
            read_filter.min_mapq = read_min_mapq;

            if (!gaf_base_filename.empty()) {
                // GAF-Base: reads are fetched per window by running gbz-base, which resolves
                // node IDs against --gbz-base or else the input graph. A GBZ-Base is read
                // randomly, while a plain GBZ is loaded in full on every query.
                string query_graph = gbz_base_filename.empty() ? graph_filename
                                                               : gbz_base_filename;
                if (read_window_size == 0) {
                    read_window_size = DEFAULT_GAF_BASE_WINDOW;
                }
                // Four windows per thread, held in one cache that all threads share.
                auto gaf_base_source = new GafBaseSiteReadSource(*graph, gaf_base_filename,
                                                                 query_graph, read_filter,
                                                                 read_window_size,
                                                                 4 * (size_t)vg::get_thread_count(),
                                                                 gaf_base_binary);
                read_source.reset(gaf_base_source);
                // Check the setup now, so that a missing binary or unreadable database is
                // reported as a user error rather than as a crash.
                try {
                    gaf_base_source->check_setup();
                } catch (const std::exception& e) {
                    logger.error() << e.what() << endl;
                }
                if (show_progress) {
                    logger.info() << "Using GAF-Base " << gaf_base_filename
                                  << " queried against " << query_graph
                                  << " via " << gaf_base_binary << endl;
                    logger.info() << "GAF-Base: databases opened by gbz-base "
                                  << (gaf_base_source->immutable_databases()
                                      ? "as SQLite immutable URIs, without file locking"
                                      : "as plain paths, with SQLite file locking (VG_GAFBASE_LOCKING=1)")
                                  << endl;
                    if (gbz_base_filename.empty()) {
                        logger.info() << "Consider building a GBZ-Base ('gbz-base construct') and "
                                      << "passing --gbz-base: a plain GBZ is reloaded on every query"
                                      << endl;
                    }
                }
            } else if (!gaf_index_filename.empty()) {
                // Indexed GAF: reads are fetched as the sites need them, from the sorted GAF
                // through its tabix index, in vg's own threads.
                if (read_window_size == 0) {
                    read_window_size = DEFAULT_GAF_INDEX_WINDOW;
                }
                // Four windows per thread, held in one cache that all threads share.
                try {
                    read_source.reset(new TabixGafSiteReadSource(*graph, gaf_filename,
                                                                 gaf_index_filename, read_filter,
                                                                 read_window_size,
                                                                 4 * (size_t)vg::get_thread_count()));
                } catch (const std::exception& e) {
                    logger.error() << e.what() << endl;
                }
                if (show_progress) {
                    logger.info() << "Using indexed GAF " << gaf_filename
                                  << " with index " << gaf_index_filename << endl;
                }
            } else if (!gam_index_filename.empty()) {
                // Indexed: reads are fetched as the sites need them, so memory is bounded by
                // what the sites need rather than by the size of the read set.
                if (read_window_size == 0) {
                    read_window_size = DEFAULT_GAM_INDEX_WINDOW;
                }
                // Four windows per thread, held in one cache that all threads share.
                read_source.reset(new IndexedGamSiteReadSource(gam_filename, gam_index_filename,
                                                               read_filter, read_window_size,
                                                               4 * (size_t)vg::get_thread_count()));
                if (show_progress) {
                    logger.info() << "Using indexed GAM " << gam_filename
                                  << " with index " << gam_index_filename << endl;
                }
            } else {
                auto in_memory_source = new InMemorySiteReadSource();
                read_source.reset(in_memory_source);
                if (!gam_filename.empty()) {
                    in_memory_source->load_gam(gam_filename, read_filter);
                } else {
                    in_memory_source->load_gaf(*graph, gaf_filename, read_filter);
                }
                if (show_progress) {
                    logger.info() << "Loaded " << in_memory_source->get_read_count()
                                  << " reads (" << in_memory_source->get_filtered_count()
                                  << " filtered out)" << endl;
                }
            }

            // Two scorers: quality-adjusted for reads that have base qualities,
            // plain for reads that do not. Picking per read avoids either
            // fabricating qualities or mis-scoring.
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

            auto* graph_calculator = new GraphAlignedAlleleLikelihoodCalculator(
                *graph, *snarl_manager, *read_source, *qual_scorer, *plain_scorer,
                likelihood_params);
            likelihood_calculator.reset(graph_calculator);
            // Place the depth-rate windows on the reference paths, so that they do not depend
            // on the node numbering. A graph without positions keeps node-ID windows.
            if (auto* position_graph = dynamic_cast<PathPositionHandleGraph*>(graph)) {
                vector<path_handle_t> rate_paths;
                for (const string& ref_path : ref_paths) {
                    if (position_graph->has_path(ref_path)) {
                        rate_paths.push_back(position_graph->get_path_handle(ref_path));
                    }
                }
                graph_calculator->set_rate_reference(position_graph, rate_paths);
            }

            auto rl_caller = new ReadLikelihoodSnarlCaller(*graph, *snarl_manager, *support_finder,
                                                           *likelihood_calculator);

            if (!dump_likelihoods_filename.empty()) {
                likelihood_dump.reset(new ofstream(dump_likelihoods_filename));
                if (!(*likelihood_dump)) {
                    logger.error() << "could not open " << dump_likelihoods_filename
                                   << " for writing" << endl;
                }
                rl_caller->set_likelihood_dump(likelihood_dump.get());
            }

            // Without a pack file the support finder reports zero for everything, so
            // the caller must not prune alleles on support.
            rl_caller->set_support_available(!support_free);
            rl_caller->set_share_discount(!no_share_quality);
            rl_caller->set_depth_quality(depth_quality);
            rl_caller->set_min_confidence(min_confidence);

            packed_caller = rl_caller;
        } else if (ratio_caller == false) {
            // Make a depth index
            if (show_progress) logger.info() << "Computing coverage statistics" << endl;
            depth_index = vg::algorithms::binned_packed_depth_index(*packer, ref_paths, min_depth_bin_width, max_depth_bin_width,
                                                                depth_scale_fac, 0, true, true);
            if (show_progress) logger.info() << "Computed coverage statistics" << endl;
            // Make a new-stype probablistic caller
            auto poisson_caller = new PoissonSupportSnarlCaller(*graph, *snarl_manager, *support_finder, depth_index,
                                                                //todo: qualities need to be used
                                                                //better in conjunction with
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
        if (min_allele_support >= 0) {
            packed_caller->set_min_supports(min_allele_support, min_allele_support, min_site_support);
        }
        
        snarl_caller = unique_ptr<SnarlCaller>(packed_caller);
    }

    if (!snarl_caller) {
        logger.error() << "pack file (-k) is required" << endl;
    }

    // Guard the pack-free path: it is only sound where nothing consults support.
    if (support_free) {
        if (!gbwt_enumeration) {
            logger.error() << "--read-likelihood without -k/--pack requires haplotype-based allele "
                           << "enumeration (-g/--gbwt or -z/--gbz); support-based enumeration needs "
                           << "a pack file" << endl;
        }
        if (!vcf_filename.empty()) {
            // VCFTraversalFinder prunes alt paths on support before its brute-force
            // enumeration, so -v needs a pack file.
            logger.error() << "-v/--vcf with --read-likelihood requires -k/--pack" << endl;
        }
        if (bottom_up) {
            // NestedFlowCaller downcasts the support finder to a nested packed one.
            logger.error() << "--bottom-up with --read-likelihood requires -k/--pack" << endl;
        }
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

    unique_ptr<GraphCaller> graph_caller;
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
            
            // todo: do we ever want to toggle in min-support?
            function<double(handle_t)> node_support = [&] (handle_t h) {
                return support_finder->support_val(support_finder->get_avg_node_support(graph->get_id(h)));
            };
            
            function<double(edge_t)> edge_support = [&] (edge_t e) {
                return support_finder->support_val(support_finder->get_edge_support(e));
            };

            // create the flow traversal finder
            FlowTraversalFinder* flow_traversal_finder = new FlowTraversalFinder(*graph, *snarl_manager, max_yens_traversals,
                                                                                 node_support, edge_support,
                                                                                 max_allele_len);
            traversal_finder = unique_ptr<TraversalFinder>(flow_traversal_finder);
        }

        if (top_down) {
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
    if (!max_snarl_edges_explicit) {
        max_snarl_edges_opt = (read_likelihood && gbwt_enumeration) ? 0 : 10000;
    }
    // The FlowCaller constructors do not take the cap, so we set it on whichever one was built.
    if (FlowCaller* flow_caller = dynamic_cast<FlowCaller*>(graph_caller.get())) {
        flow_caller->set_max_snarl_edges(max_snarl_edges_opt);
    }

    // The caller as a VCFOutputCaller, or null if it does not write VCF.
    VCFOutputCaller* const vcf_out = dynamic_cast<VCFOutputCaller*>(graph_caller.get());

    // Per-region ploidy, if given.
    if (!ploidy_bed_filename.empty()) {
        VCFOutputCaller* ploidy_target = vcf_out;
        if (ploidy_target == nullptr) {
            cerr << "error [vg call]: --ploidy-bed needs a caller that emits VCF" << endl;
            return 1;
        }
        ploidy_target->set_ploidy_regions(ploidy_bed_filename);
    }

    // Nested calling: a called traversal that takes the reference's route through a snarl,
    // differing only inside nested snarls, is called as the reference allele, and the differences
    // are called at the nested snarls. It is on by default only for the read-likelihood genotyper;
    // other callers use it when --nested is given.
    if (nested_calling && !nested_explicit && !read_likelihood) {
        nested_calling = false;
    }
    // An explicit --regenotype without read phasing is an error. One set by a preset is turned
    // off, as with --preset ont --no-read-phasing.
    if (regenotype && !read_phasing) {
        if (regenotype_explicit) {
            cerr << "error [vg call]: --regenotype needs --read-phasing, which gives each read"
                 << " the strand log-odds it uses" << endl;
            return 1;
        }
        regenotype = false;
    }
    // Re-genotyping happens in FlowCaller::phase_and_regenotype(). --bottom-up uses
    // NestedFlowCaller, which does not have it. --top-down gives each child snarl candidate
    // traversals derived from its parent's called genotype, so changing that genotype afterwards
    // would leave the child genotyped against the wrong alleles. An explicit --regenotype is an
    // error with either; one set by a preset is turned off.
    if (regenotype && (top_down || bottom_up) && !regenotype_explicit) {
        regenotype = false;
    }
    if (regenotype && (top_down || bottom_up)) {
        cerr << "error [vg call]: --regenotype cannot be combined with "
             << (top_down ? "--top-down" : "--bottom-up") << "; "
             << (top_down
                 ? "--top-down derives each child's candidate traversals from its parent's called"
                   " genotype, so re-genotyping a parent would leave its children genotyped against"
                   " alleles it no longer carries"
                 : "--bottom-up uses a caller that cannot re-genotype")
             << endl;
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
    if ((anchor_params.phase_hets || anchor_params.strict_hets) && phase_hets_explicit
        && !read_phasing) {
        cerr << "error [vg call]: --anchors-phase-hets / --anchors-strict-hets need --read-phasing,"
             << " which gives each read the strand log-odds they use" << endl;
        return 1;
    }
    // -A already calls every nested snarl as a snarl of its own, so with nested calling each would
    // be called twice.
    if (nested_calling && all_snarls) {
        if (nested_explicit) {
            cerr << "error [vg call]: -A/--all-snarls calls every snarl independently while"
                 << " --nested calls children through their parent; they are alternatives, not a"
                 << " combination" << endl;
            return 1;
        }
        nested_calling = false;
    }
    {
        // The linkage model can change a record's GQN, so the output caller re-applies the lowconf
        // filter with the same threshold.
        VCFOutputCaller* confidence_target = vcf_out;
        if (confidence_target != nullptr) {
            confidence_target->set_linkage_min_confidence(min_confidence);
        }
    }
    if (nested_calling) {
        VCFOutputCaller* nested_target = vcf_out;
        if (nested_target == nullptr) {
            if (nested_explicit) {
                cerr << "error [vg call]: --nested needs a caller that emits VCF" << endl;
                return 1;
            }
            nested_calling = false;
        }
    }
    // Block emission splits up the records of nested calling, so it needs nested calling and a
    // caller that writes VCF. The checks that depend only on the options were made before
    // the graph was loaded.
    if (atomize_blocks) {
        if (!nested_calling) {
            if (atomize_explicit) {
                cerr << "error [vg call]: --atomize-blocks needs nested calling (--nested), which is"
                     << " on by default under --read-likelihood" << endl;
                return 1;
            }
            atomize_blocks = false;
        }
        VCFOutputCaller* atomize_target =
            atomize_blocks ? vcf_out : nullptr;
        if (atomize_blocks && atomize_target == nullptr) {
            if (atomize_explicit) {
                cerr << "error [vg call]: --atomize-blocks needs a caller that emits VCF" << endl;
                return 1;
            }
            atomize_blocks = false;
        }
        if (atomize_blocks) {
            atomize_target->set_atomize_blocks(true);
        }
    }

    // Nested calling would apply to only some sites under -I, for the reason given for
    // --atomize-blocks above. This check does not depend on VCF output, since -I is mostly used
    // with GAF output.
    if (nested_calling && call_chains) {
        if (nested_explicit) {
            cerr << "error [vg call]: --nested cannot be combined with -I/--chains, which calls"
                 << " pieces of chains; nested calling would apply only at pieces that hold a"
                 << " single snarl" << endl;
            return 1;
        }
        nested_calling = false;
        if (show_progress) {
            logger.info() << "Nested calling is off under -I/--chains, which calls pieces of chains"
                          << " rather than snarls" << endl;
        }
    }

    if (nested_calling) {
        VCFOutputCaller* nested_target = vcf_out;
        nested_target->set_symbolic_collapsing(snarl_manager.get());
    }

    // Owned here because write_variants(), at the very end of main, consumes the collector.
    unique_ptr<LinkageCollector> linkage_collector;
    vector<size_t> linkage_sequence_to_haplotype;

    string header;
    if (!gaf_output) {
        // Init The VCF       
        VCFOutputCaller* vcf_caller = vcf_out;
        assert(vcf_caller != nullptr);
        // Write the nesting INFO tags (such as LV and PS) with -A, --top-down or --bottom-up, and
        // when off-reference chains are genotyped: a record on a gref fragment contig lies inside
        // a non-reference allele, and only these tags say so.
        vcf_caller->set_nested(all_snarls || top_down || bottom_up
                               || off_ref_nesting);
        vcf_caller->set_off_reference_nesting(off_ref_nesting);
        vcf_caller->set_gref_levels(std::move(gref_levels));
        vcf_caller->set_translation(translation.get());

        // The linkage model compares genotypes with the haplotypes (-z or -g), and uses the
        // read-likelihood genotyper's likelihoods. Without both, the default weight is set to 0
        // and an explicit weight is an error. This is decided first because phasing, the mosaic
        // and nested calling depend on it.
        if (linkage_weight > 0.0 && !(gbwt_index != nullptr && read_likelihood)) {
            if (linkage_weight_explicit) {
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
        if (!mosaic_out.empty() && phased_explicit && !phased_output) {
            // The mosaic is written from the phasing, so it cannot be written without it.
            logger.error() << "--mosaic-out cannot be combined with --no-phased, because the mosaic "
                           << "is written from the phasing" << endl;
        }
        if (phased_output && linkage_weight <= 0.0) {
            // Phasing comes from the linkage model. Without it, an explicit --phased is an error,
            // and the default is turned off.
            if (phased_explicit) {
                logger.error() << "--phased needs the linkage model, which needs haplotype "
                               << "enumeration (-z or -g) and --read-likelihood" << endl;
            }
            phased_output = false;
        }
        if (nested_calling && linkage_weight > 0.0 && !phased_output) {
            // Where the linkage model runs, a nested site's ploidy and strand come from its
            // parent's phased genotype, so nested calling needs phasing. An explicit --nested is
            // an error; the default is turned off.
            if (nested_explicit) {
                cerr << "error [vg call]: --nested needs --phased where the linkage model runs, "
                     << "because a nested site's ploidy and strand come from its parent's phased "
                     << "genotype" << endl;
                return 1;
            }
            // Symbolic collapsing, which turns on nested calling in the caller, was set up above,
            // so it has to be turned off there too.
            nested_calling = false;
            VCFOutputCaller* nested_target = vcf_out;
            if (nested_target != nullptr) {
                nested_target->set_symbolic_collapsing(nullptr);
            }
            if (show_progress) {
                logger.info() << "Nested calling is off under --no-phased: a nested site's ploidy "
                              << "and strand come from its parent's phased genotype" << endl;
            }
        }
        if (linkage_weight > 0.0) {
            // Map each GBWT sequence to the index of its haplotype in the panel. A haplotype is a
            // (sample, phase) pair, and may be stored as several paths (one per contig or
            // fragment), each with a sequence per orientation. They all get the same index, so a
            // haplotype counts once.
            const gbwt::Metadata& meta = gbwt_index->metadata;
            // The samples that are not gref-derived; see below.
            unordered_set<string> base_samples;
            for (gbwt::size_type i = 0; i < meta.sample_names.size(); ++i) {
                const string s = meta.sample(i);
                if (!GrefCover::is_gref_derived(s)) {
                    base_samples.insert(s);
                }
            }
            map<pair<size_t, size_t>, size_t> hap_index;
            linkage_sequence_to_haplotype.assign(gbwt_index->sequences(), 0);
            for (gbwt::size_type path = 0; path < meta.paths(); ++path) {
                const gbwt::PathName& name = meta.path(path);
                // Leave gref cover paths out of the panel. The fragment paths all share one sample
                // and phase, so they would form one haplotype stitched together from many donors.
                // The gref copy of a base reference path is left out only if the base path's
                // sample is also present, so that the reference is in the panel once. Sequences
                // left out are mapped to WILDCARD; the vector's default of 0 would put them in
                // haplotype 0.
                const string path_sample = (size_t)name.sample < meta.sample_names.size()
                                               ? meta.sample(name.sample) : string();
                const string path_contig = (size_t)name.contig < meta.contig_names.size()
                                               ? meta.contig(name.contig) : string();
                const bool gref_fragment = GrefCover::is_gref_derived(path_sample)
                                           && GrefCover::is_gref_name(path_contig);
                const bool gref_shadowing_base =
                    GrefCover::is_gref_derived(path_sample) && !gref_fragment
                    && base_samples.count(path_sample.substr(GrefCover::gref_prefix.size())) > 0;
                if (gref_fragment || gref_shadowing_base) {
                    for (gbwt::size_type orientation = 0; orientation < 2; ++orientation) {
                        gbwt::size_type seq = gbwt::Path::encode(path, orientation);
                        if (seq < linkage_sequence_to_haplotype.size()) {
                            linkage_sequence_to_haplotype[seq] = LinkageModel::WILDCARD;
                        }
                    }
                    continue;
                }
                auto key = make_pair((size_t)name.sample, (size_t)name.phase);
                auto found = hap_index.find(key);
                size_t index;
                if (found == hap_index.end()) {
                    index = hap_index.size();
                    hap_index[key] = index;
                } else {
                    index = found->second;
                }
                for (gbwt::size_type orientation = 0; orientation < 2; ++orientation) {
                    gbwt::size_type seq = gbwt::Path::encode(path, orientation);
                    if (seq < linkage_sequence_to_haplotype.size()) {
                        linkage_sequence_to_haplotype[seq] = index;
                    }
                }
            }
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
                linkage_params.hp_prior = hp_prior;
                linkage_params.hp_prior_run = (size_t)hp_prior_run;
                linkage_collector.reset(new LinkageCollector(linkage_params, hap_index.size()));
                vcf_caller->set_linkage(linkage_collector.get(), gbwt_index,
                                        &linkage_sequence_to_haplotype);
                vcf_caller->set_emit_phasing(phased_output);
                // Name each panel haplotype "sample#phase", for the mosaic.
                vector<string> hap_names(hap_index.size());
                for (const auto& kv : hap_index) {
                    string sample = kv.first.first < meta.sample_names.size()
                                        ? meta.sample(kv.first.first)
                                        : string("sample") + std::to_string(kv.first.first);
                    hap_names[kv.second] = sample + "#" + std::to_string(kv.first.second);
                }
                // The mosaic's rows name only the contig, so it also records the full names of
                // the reference paths.
                vcf_caller->set_mosaic_out(mosaic_out, graph_filename, hap_names, ref_paths,
                                           mosaic_patch_gaps, mosaic_keep_nested,
                                           mosaic_connect_unexplained);
            }
            if (show_progress) {
                logger.info() << "Linkage: " << hap_index.size() << " panel haplotypes over "
                              << gbwt_index->sequences() << " GBWT sequences" << endl;
            }
        }
        if (!anchors_out.empty()) {
            // The anchor file names the read file, since its read offsets refer to those reads.
            string reads_source = !gaf_base_filename.empty()
                                      ? gaf_base_filename
                                      : (!gaf_filename.empty() ? gaf_filename : gam_filename);
            vcf_caller->set_anchors_out(anchors_out, anchor_params, graph_filename, reads_source,
                                        min_mismap_prob);
        }
        vcf_caller->set_read_phasing(read_phasing, read_phasing_params);
        vcf_caller->set_regenotype(regenotype, regenotype_params, regenotype_passes,
                                   regenotype_ledger);
        // one call covers FlowCaller (both ctors, so plain vg call gets it too), NestedFlowCaller
        // and LegacyCaller, since the merge lives on the shared VCFOutputCaller base
        vcf_caller->set_allele_merge(cluster_threshold, cluster_min_allele_len);
        // Make sure the basepath information we inferred above goes directy to the VCF header
        // (and that it does *not* try to read it from the graph paths)
        vector<string> header_ref_paths;
        vector<size_t> header_ref_lengths;
        bool need_overrides = dynamic_cast<VCFGenotyper*>(graph_caller.get()) == nullptr;
        for (const auto& path_len : basepath_length_map) {
            // Use the locus name (eg "chrI") to match what emit_variant()
            // writes to the CHROM column in data lines
            string contig_name = PathMetadata::parse_locus_name(path_len.first);
            header_ref_paths.push_back(contig_name != PathMetadata::NO_LOCUS_NAME ? contig_name : path_len.first);
            if (need_overrides) {
                header_ref_lengths.push_back(path_len.second);
            }
        }
        header = vcf_caller->vcf_header(*graph, header_ref_paths, header_ref_lengths);
    }

    graph_caller->set_show_progress(show_progress);
    
    // Call the graph
    // Determine recursion strategy based on mode:
    // - top_down: FlowCaller handles recursion internally, so RecurseNever
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

    // A read source that fetches reads by node-ID window works best when snarls are visited in
    // node-ID order, so that each fetched window serves many sites in a row.
    if (dynamic_cast<WindowedSiteReadSource*>(read_source.get()) != nullptr) {
        graph_caller->set_node_id_ordering(true, read_window_size);
        if (show_progress) {
            logger.info() << "Visiting snarls in node-ID order, window " << read_window_size << endl;
        }
    }

    // After the direct pass, choose the genotypes, parents before their nested children (the
    // linkage pass, FlowCaller::run_linkage_pass), and build each record from its settled
    // genotype. Nested calling needs this because a child's ploidy depends on its parent's
    // genotype, and the linkage model needs it because it can change genotypes after they are
    // first called.
    FlowCaller* deferring_caller = nullptr;
    if (nested_calling || linkage_collector != nullptr) {
        deferring_caller = dynamic_cast<FlowCaller*>(graph_caller.get());
        if (deferring_caller != nullptr) {
            deferring_caller->set_stage_records(true);
        }
    }

    if (!call_chains) {
        // Call each snarl
        if (show_progress) logger.info() << "Calling top-level snarls" << endl;
        graph_caller->call_top_level_snarls(*graph, recurse_type);
    } else {
        // Attempt to call chains instead of snarls so that the output traversals are longer
        // Todo: this could probably help in some cases when making VCFs too
        if (show_progress) logger.info() << "Calling top-level chains" << endl;
        graph_caller->call_top_level_chains(*graph, max_chain_edges, max_chain_trivial_travs, recurse_type);
    }

    // Calling is done, and the passes below work from what it kept, not from reads: free the
    // windows the read source still caches, before those passes reach the run's highest memory use.
    if (auto* windowed = dynamic_cast<WindowedSiteReadSource*>(read_source.get())) {
        size_t freed = windowed->drop_cached_windows();
        if (show_progress) {
            logger.info() << "Freed the read cache: " << freed << " reads" << endl;
        }
    }

    if (deferring_caller != nullptr) {
        // Round 1's linkage pass, then read phasing and any further rounds, then the render.
        deferring_caller->run_linkage_pass();
        deferring_caller->phase_and_regenotype();
        deferring_caller->render_retained_records();
    }

    // Anchors are collected while records are built, so they are written afterwards.
    {
        auto* anchor_caller = vcf_out;
        if (anchor_caller != nullptr) {
            anchor_caller->write_anchors();
        }
    }

    // Report how the windowed read source's fetches and cache performed.
    if (show_progress) {
        auto* windowed = dynamic_cast<WindowedSiteReadSource*>(read_source.get());
        if (windowed != nullptr) {
            size_t hits = windowed->get_cache_hits();
            size_t misses = windowed->get_cache_misses();
            size_t total = hits + misses;
            auto* gaf_base = dynamic_cast<GafBaseSiteReadSource*>(windowed);
            auto* tabix_gaf = dynamic_cast<TabixGafSiteReadSource*>(windowed);
            const char* label = gaf_base != nullptr ? "GAF-Base: "
                                : (tabix_gaf != nullptr ? "Indexed GAF: " : "Indexed GAM: ");
            logger.info() << label
                          << windowed->get_read_count() << " reads fetched, "
                          << hits << "/" << total << " site queries served from cache"
                          << (total > 0 ? " (" + std::to_string((int)(100.0 * hits / total)) + "%)" : "")
                          << endl;
            // Reads examined by site queries, against reads delivered to sites.
            size_t seen = windowed->get_scanned_count();
            size_t used = windowed->get_delivered_count();
            logger.info() << label
                          << seen << " read candidates examined, " << used << " delivered"
                          << (seen > 0 ? " (" + std::to_string((int)(100.0 * used / seen)) + "%)" : "")
                          << endl;
            // Sites too big for one window are fetched uncached, by their exact node ranges.
            logger.info() << label
                          << windowed->get_straddle_count() << " site queries too wide for a "
                          << "window, fetched uncached over " << windowed->get_straddle_wanted()
                          << " node IDs (spanning " << windowed->get_straddle_nodes() << ")"
                          << endl;
            logger.info() << label
                          << windowed->get_whole_fetches() << " windows fetched whole, holding "
                          << windowed->get_whole_fetch_reads() << " reads" << endl;
            if (gaf_base != nullptr) {
                // Each query runs a subprocess, so this count largely determines run time.
                logger.info() << "GAF-Base: " << gaf_base->get_query_count()
                              << " subprocess queries" << endl;
                // Alignments returned twice by one query, which would otherwise count twice.
                logger.info() << "GAF-Base: " << gaf_base->get_duplicate_count()
                              << " duplicate reads dropped" << endl;
                // Where the query time went. File locking, for one, shows up as kernel time.
                auto totals = gaf_base->get_query_totals();
                double cpu = totals.child_user_s + totals.child_system_s;
                logger.info() << "GAF-Base: threads waited " << totals.wait_s / 3600 << " h for "
                              << totals.queries << " gbz-base queries, which used "
                              << totals.child_user_s / 3600 << " h of CPU in their own code and "
                              << totals.child_system_s / 3600 << " h in the kernel ("
                              << (cpu > 0 ? (int)(100 * totals.child_system_s / cpu) : 0)
                              << "%); parsing their " << totals.gaf_bytes / 1e9 << " GB of GAF took "
                              << totals.parse_s / 3600 << " h" << endl;
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
    }

    if (show_progress) {
        auto* graph_calculator =
            dynamic_cast<GraphAlignedAlleleLikelihoodCalculator*>(likelihood_calculator.get());
        if (graph_calculator != nullptr && graph_calculator->rate_id_fallbacks() > 0) {
            logger.info() << graph_calculator->rate_id_fallbacks() << " sites had no reference "
                          << "position for their depth-rate window and used a node-ID window"
                          << endl;
        }
    }

    if (show_progress) logger.info() << "Calling complete" << endl;

    if (!gaf_output) {
        // Output VCF
        VCFOutputCaller* vcf_caller = vcf_out;
        assert(vcf_caller != nullptr);
        // Prune here rather than where the header string was built: the contig set is only
        // known once calling is done, and write_variants below drains the buffer it reads.
        cout << vcf_caller->prune_header_contigs(header, vcf_caller->get_output_contigs()) << flush;
        if (show_progress) logger.info() << "Writing VCF Variants" << endl;
        vcf_caller->write_variants(cout, snarl_manager.get());
        if (show_progress) logger.info() << "VCF complete" << endl;        
    }

    // Everything is written. Freeing what calling built -- the graph, the snarls, every site's
    // records and evidence -- takes minutes on a whole genome and changes nothing, so close what is
    // still open and exit without it. exit() still runs the handlers registered for exit, such as
    // the one that removes temporary files, and flushes C streams.
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

