/** \file call_main.cpp
 *
 * Defines the "vg call" subcommand, which calls variation from an augmented graph
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <regex>
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

/// Default --read-window per backend. A GAF-Base query costs a process spawn, so it
/// wants a window big enough to amortise one; an indexed GAM query is a seek into an
/// already-open file, and a large window there just over-fetches.
///
/// The GAF-Base default is set by long reads, because they are what a narrow window
/// punishes: a window has to be much wider than a read's node-ID span or every read
/// straddling a boundary is fetched twice. Measured on chr20, 44x ONT (33 kb mean) with
/// anchors: 4096 -> 293 s and 1.62 fetches per read, 8192 -> 265 s, 16384 -> 182 s and
/// 1.44, 32768 -> 195 s, 65536 -> 295 s at 10.7 GB. 30x Illumina over the same window
/// range moves from 156 s to 154 s, so the wider window costs short reads nothing.
///
/// It is not output-neutral, and the reason is worth stating: a window wider than
/// `max_query_nodes` is fetched as several concurrent queries whose results are
/// concatenated, which reorders the reads a site sees and so reorders a floating-point
/// sum. On chr20 short reads that moves `GL` in the sixth decimal for 3,402 of 115,410
/// records and changes **no genotype at all**; on ONT the output is byte-identical.
const size_t DEFAULT_GAM_INDEX_WINDOW = 256;
const size_t DEFAULT_GAF_BASE_WINDOW = 16384;

/// Count the distinct haplotypes a GBZ's GBWT can offer GBWTTraversalFinder.
///
/// Counts (sample, haplotype) pairs over HAPLOTYPE-sense paths only. Reference and
/// generic paths live in the same GBWT but are not panel members: a GBZ carrying only
/// GRCh38 and CHM13 holds two paths and no panel at all, and enumerating from it would
/// offer nothing but the reference allele at every site. A haplotype broken into
/// fragments owns several paths, so collapse onto the (sample, haplotype) key rather
/// than counting paths.
///
/// Deliberately not the same count as the linkage layer's, which collapses every path
/// including the reference ones: there a reference assembly is a legitimate panel
/// member to impute against. The question here is narrower, and is only whether there
/// is any alternative to the reference to enumerate at all.
static size_t count_panel_haplotypes(const PathHandleGraph& graph) {
    set<pair<string, size_t>> haplotypes;
    graph.for_each_path_of_sense(PathSense::HAPLOTYPE, [&](const path_handle_t& path) {
        haplotypes.emplace(graph.get_sample_name(path), graph.get_haplotype(path));
    });
    return haplotypes.size();
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
         << "read-likelihood calling options (all require --read-likelihood):" << endl
         << "  every term and parameter below: doc/read-likelihood-genotyping.md" << endl
         << "      --read-likelihood     genotype from an explicit P(reads|genotype) model" << endl
         << "                            instead of aggregate depth (needs one read source" << endl
         << "                            below)" << endl
         << "" << endl
         << "  reads in (one source required):" << endl
         << "      --gam FILE            read alignments for --read-likelihood" << endl
         << "      --gaf-reads FILE      read alignments for --read-likelihood, as GAF" << endl
         << "      --gam-index FILE      .gai index for --gam, so reads are fetched per site" << endl
         << "                            instead of all held in memory (from vg gamsort -i)" << endl
         << "      --gaf-base FILE       GAF-Base of read alignments, fetched per site by" << endl
         << "                            running gbz-base (needs it on the PATH)" << endl
         << "      --gbz-base FILE       graph to resolve --gaf-base queries against, as a" << endl
         << "                            GBZ-Base or GBZ [the input graph]" << endl
         << "      --gaf-base-binary P   gbz-base executable to run [gbz-base]" << endl
         << "      --read-window N       node-ID window for indexed read fetches" << endl
         << "                            [16384 for --gaf-base, 256 for --gam-index]" << endl
         << "      --read-min-mapq N     ignore reads with MAPQ below N [0]" << endl
         << "" << endl
         << "  allele enumeration:" << endl
         << "      --enumerate-support   enumerate candidate alleles from read support rather" << endl
         << "                            than from the GBZ haplotype panel. The panel is the" << endl
         << "                            default here, as -z gives elsewhere: it measured" << endl
         << "                            better on every small-variant class tested and needs" << endl
         << "                            no pack file. This flag opts out, and then -k/--pack" << endl
         << "                            is needed again. Worth it where the sample's" << endl
         << "                            variation is poorly represented in the panel, since" << endl
         << "                            panel enumeration can never spell an allele no" << endl
         << "                            haplotype carries. A panel of under 2 haplotypes" << endl
         << "                            falls back to support on its own" << endl
         << "" << endl
         << "  model terms:" << endl
         << "      --depth-term W        add W * ln P(N reads | genotype) to the likelihood," << endl
         << "                            so a genotype is also judged on whether it predicts" << endl
         << "                            the number of reads seen. The rate is measured over" << endl
         << "                            the read source's own fetch window, so it costs no" << endl
         << "                            extra I/O. 0 disables it; DR is emitted either way" << endl
         << "                            [0.1]" << endl
         << "      --depth-count-raw     count every read as one read of depth, instead of as" << endl
         << "                            1 - e_r, the probability it came from this locus" << endl
         << "      --no-mismap-term      disable the MAPQ-derived mismapping term" << endl
         << "      --mismap-max P        upper clamp on the MAPQ-derived mismapping" << endl
         << "                            probability. Governs how much a read's placement" << endl
         << "                            ambiguity counts, so it matters most on graphs with" << endl
         << "                            many similar haplotypes [0.7]" << endl
         << "      --mismap-min P        lower clamp: floor on how unreliable any read may" << endl
         << "                            be, capping one read's veto at ln(P). Covers local" << endl
         << "                            misalignment, which MAPQ does not measure. Mainly an" << endl
         << "                            indel knob; interacts with --mismap-max [0.02]" << endl
         << "      --preset NAME         a fitted parameter set for one read type." << endl
         << "                            `ont`: --gap-open 1 --gap-extend 1 --mismap-min" << endl
         << "                            0.05, which on 43x ONT chr20 moves indel GT F1 0.749" << endl
         << "                            -> 0.816 and ALL 0.926 -> 0.945 at no cost to SNVs," << endl
         << "                            reproducing at +0.061 on a held-out contig and" << endl
         << "                            +0.069 at matched 30x coverage. An explicit flag" << endl
         << "                            overrides the preset either side of it" << endl
         << "      --gap-open N          read scorer's gap-open penalty [6]" << endl
         << "      --read-phasing        phase from the reads that span consecutive hets," << endl
         << "                            not from the haplotype panel alone. On under" << endl
         << "                            `--preset ont`" << endl
         << "      --no-read-phasing     leave the phase to the haplotype panel" << endl
         << "      --phase-min-q N       a site below this per-read confidence may not" << endl
         << "                            carry a phase link [9.5]" << endl
         << "      --phase-break N       break the chain below this many log10 units [10]" << endl
         << "      --phase-relink N      reliable sites either side of a break [3]" << endl
         << "      --phase-hang N        neighbours to hang an unreliable site from [4]" << endl
         << "      --phase-prior N       weight of the panel when hanging a site [3]" << endl
         << "      --phase-cap N         clamp one pair's contribution, 0 to disable [0]" << endl
         << "      --regenotype          let the reads' phase decide the genotype, not only" << endl
         << "                            the order of an already-settled pair. Needs" << endl
         << "                            `--read-phasing`; on under `--preset ont`" << endl
         << "      --no-regenotype       leave the genotype to the per-site likelihoods" << endl
         << "      --regeno-temper N     how hard to believe a read's strand. 0 reproduces" << endl
         << "                            the uncorrected caller byte for byte; the default" << endl
         << "                            is fitted from the run's own data, with no truth" << endl
         << "      --regeno-passes N     cap on barrier passes [2], i.e. one correction" << endl
         << "                            round. The iteration stops early if the genotypes" << endl
         << "                            settle or start repeating a state; on chr20 they do" << endl
         << "                            neither, and more rounds score slightly worse. 1" << endl
         << "                            scores the correction without acting on it" << endl
         << "      --regeno-ledger FILE  write one line per site the correction would move" << endl
         << "      --regeno-shuffle      DEBUG. Randomise each read's strand sign, keeping" << endl
         << "                            the magnitude, so peakedness survives and only the" << endl
         << "                            phase information is destroyed" << endl
         << "      --gap-extend N        read scorer's gap-extension penalty [1]. Together" << endl
         << "                            these set how hard a read votes against an allele" << endl
         << "                            that differs from it by an indel, and they are the" << endl
         << "                            one scoring primitive base quality never softens. At" << endl
         << "                            6/1 a 1 bp difference is a saturated vote, which" << endl
         << "                            suits a read whose errors are substitutions and not" << endl
         << "                            a read whose modal error is a single-base indel" << endl
         << "      --flat-mixture        weight each haplotype of a genotype equally" << endl
         << "                            (1/ploidy) instead of by the reads it is expected to" << endl
         << "                            contribute. The flat weight is wrong wherever the" << endl
         << "                            alleles differ in length: it loses large" << endl
         << "                            heterozygous deletions and mis-genotypes large" << endl
         << "                            heterozygous insertions. Restores the pre-correction" << endl
         << "                            model exactly" << endl
         << "" << endl
         << "  linkage between sites (needs panel enumeration, so off" << endl
         << "  under --enumerate-support):" << endl
         << "      --linkage-weight W    re-decide genotypes with a Li-Stephens model over" << endl
         << "                            the GBWT haplotypes, so consecutive calls are judged" << endl
         << "                            against combinations the panel carries. Declines" << endl
         << "                            quietly where enumeration is absent, unless asked" << endl
         << "                            for explicitly. 0 is off and reproduces the per-site" << endl
         << "                            caller exactly. Tuned on a 34-haplotype panel;" << endl
         << "                            roughly neutral on 4 [2]" << endl
         << "      --linkage-scale N     distance scale of the linkage decay, in bp [10000]" << endl
         << "      --mosaic-break-unexplained" << endl
         << "                            break the path where the panel cannot explain a" << endl
         << "                            stretch, instead of carrying the flanking haplotype" << endl
         << "                            through it. Connecting is the default and is rare" << endl
         << "      --no-mosaic-nested    merge haplotype switches that happen inside a nested" << endl
         << "                            chain into the enclosing run: fewer switches, still" << endl
         << "                            one contiguous walk, but the parent's route through" << endl
         << "                            the child" << endl
         << "      --no-mosaic-patch-gaps" << endl
         << "                            leave a gap where no panel haplotype can be carried" << endl
         << "                            across it, instead of filling it with the reference." << endl
         << "                            Patching is on by default and marked `ref` in the" << endl
         << "                            file" << endl
         << "      --anchors-out FILE    write pangenome-guided assembly anchors to FILE. An" << endl
         << "                            anchor is a zero-length pin at a snarl boundary," << endl
         << "                            holding the reads that cross it partitioned by which" << endl
         << "                            called allele they fit. Implies --read-likelihood" << endl
         << "      --anchors-end-new N   emit the end pin only where it holds at least N" << endl
         << "                            reads the start pin does not. Both pins carry the" << endl
         << "                            SAME partition, so where their read sets agree the" << endl
         << "                            two are joined by all the same reads and the end pin" << endl
         << "                            can offer no linkage the start pin does not [0," << endl
         << "                            always emit it]" << endl
         << "      --anchors-het-only    only heterozygous sites. By default homozygous and" << endl
         << "                            haploid ones are emitted too: they carry no" << endl
         << "                            haplotype information, but an anchor graph is built" << endl
         << "                            out of contiguity as much as out of phasing" << endl
         << "      --anchors-leaf-only   only leaf snarls [every genotyped site]" << endl
         << "      --anchors-reads N     minimum reads per anchor [2]" << endl
         << "      --anchors-min-gqn F   minimum site GQN for its partition to be trusted [0]" << endl
         << "      --anchors-min-q F     minimum per-read assignment phred. 3 is the measured" << endl
         << "                            knee -- ~2 points of purity for ~9% of reads -- but" << endl
         << "                            0 by default, since a read with no MAPQ cannot clear" << endl
         << "                            any positive threshold. Bounded above" << endl
         << "                            by --mismap-min: at 0.02 a perfectly discriminating" << endl
         << "                            read on a balanced het scores 14, not 60 [0]" << endl
         << "      --anchors-keep-off-call" << endl
         << "                            keep reads whose best-fitting allele is not a called" << endl
         << "                            one. Off by default: for assembly anchors purity" << endl
         << "                            beats yield, and such a read fits neither called" << endl
         << "                            allele" << endl
         << "      --mosaic-out FILE     write the inferred genome as a mosaic of panel" << endl
         << "                            haplotypes: one line per maximal run of one strand" << endl
         << "                            on one haplotype, anchored on node IDs so it is read" << endl
         << "                            back against the graph rather than a reference. Only" << endl
         << "                            switch points are stored, so it is smaller than" << endl
         << "                            explicit paths by about two orders of magnitude." << endl
         << "                            Implies --phased" << endl
         << "      --no-phased           emit unphased genotypes (0/1) and no FORMAT/PS." << endl
         << "                            Phasing is on by default wherever the linkage model" << endl
         << "                            runs" << endl
         << "      --phased              emit phased genotypes (0|1) and a FORMAT/PS phase" << endl
         << "                            set, from the linkage layer's most probable path of" << endl
         << "                            haplotype pairs. Phase comes from the panel, not" << endl
         << "                            from reads spanning sites, so a phase set is a whole" << endl
         << "                            chain rather than a read-length block. Constrained" << endl
         << "                            to the genotypes actually emitted, so GT stays a" << endl
         << "                            permutation of the unphased call." << endl
         << "                            Needs --linkage-weight above 0" << endl
         << "      --linkage-prior F     exponent on the panel allele-frequency prior implied" << endl
         << "                            by the state space. Only acts with --linkage-weight." << endl
         << "                            0 removes it, 1 keeps it as the states present it," << endl
         << "                            and above 1 amplifies it; measured best near 5 on a" << endl
         << "                            34- haplotype panel, inverting past 8. Mostly an" << endl
         << "                            indel effect [5]" << endl
         << "" << endl
         << "  quality reporting (ranking only; these never change a genotype):" << endl
         << "      --no-share-quality    report GQ as the raw likelihood ratio, without" << endl
         << "                            scaling it by the fraction of reads the called" << endl
         << "                            genotype explains. GQI always carries the raw value" << endl
         << "      --depth-quality A     scale GQ by exp(-A * |ln DR|) at records whose" << endl
         << "                            called alleles change length by 50 bp or more, so a" << endl
         << "                            call whose read count is implausible for the" << endl
         << "                            sequence it claims ranks lower. Ranking only; no" << endl
         << "                            genotype changes. 0 is off [0]" << endl
         << "      --min-confidence X    mark records whose GQN is below X as FILTER=lowconf." << endl
         << "                            GQN is a fraction of what the site could achieve, so" << endl
         << "                            one threshold means the same thing at any depth and" << endl
         << "                            ploidy; a raw GQ threshold does not, and GQ >= 10" << endl
         << "                            costs a 5x diploid contig a third of its F1. 0.05" << endl
         << "                            raises precision on every arm measured, for 1-2% of" << endl
         << "                            recall. Marks, never drops. There is no good" << endl
         << "                            default: it helps haploid F1 and hurts diploid, so" << endl
         << "                            the choice is yours. 0 is off [0]" << endl
         << "" << endl
         << "  debugging:" << endl
         << "      --dump-likelihoods F  write the per-site read/allele matrix to F as TSV" << endl
         << "  -b, --het-bias M,N        homozygous alt/ref allele must have >= M/N times" << endl
         << "                            more support than the next best allele [6,6]" << endl
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
         << "  -A, --all-snarls          call all snarls including nested (each independent)" << endl
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
         << "                            (applies only if input graph is GBZ; already the" << endl
         << "                            default under --read-likelihood)" << endl
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
         << "      --no-nested           genotype each snarl against its own full traversals," << endl
         << "                            without collapsing or descent. Nested calling is on" << endl
         << "                            by default under --read-likelihood, where it is" << endl
         << "                            measured: it takes genome-wide SNV F1 from 0.9752 to" << endl
         << "                            0.9833 and SV F1 from 0.5134 to 0.5467" << endl
         << "      --nested              treat a called traversal that differs from the" << endl
         << "                            reference only inside a nested chain as the" << endl
         << "                            reference allele, so its differences are left to the" << endl
         << "                            nested sites that contain them rather than emitted" << endl
         << "                            as one long substitution. A chain is descended into" << endl
         << "                            once its parent's genotype is settled: after the" << endl
         << "                            linkage pass where linkage can still move it, and" << endl
         << "                            immediately where it cannot" << endl
         << "      --atomize-blocks      on by default under --read-likelihood. Aligns the" << endl
         << "                            reference and each called haplotype as symbolic" << endl
         << "                            alleles and emits one record per difference block," << endl
         << "                            so a snarl differing in two separated places reports" << endl
         << "                            two variants rather than one substitution spanning" << endl
         << "                            both. Worth +0.0043 SV F1 genome-wide on HG002" << endl
         << "                            (0.5577 -> 0.5620 over 23 contigs: 48 more true SVs," << endl
         << "                            266 fewer false), though the per-contig spread is" << endl
         << "                            wide -- +0.010 on chr20 against +0.002 on chr6." << endl
         << "                            Higher under truvari refine, which controls for the" << endl
         << "                            representation change splitting causes, so the" << endl
         << "                            record-matching metric understates it. Small-variant" << endl
         << "                            F1 is unchanged. Declines on any other calling path," << endl
         << "                            and with -a, whose record set has to stay" << endl
         << "                            sample-independent; asking for it there by name is" << endl
         << "                            an error rather than a decline. CAVEAT: every block" << endl
         << "                            of a snarl repeats the site's AD, GL, GQ, GQI, GP" << endl
         << "                            and QUAL, so evidence summed across the records of" << endl
         << "                            one snarl is counted more than once. INFO/SB" << endl
         << "                            identifies the set. Pass --no-atomize-blocks for one" << endl
         << "                            record per snarl" << endl
         << "      --no-atomize-blocks   one record per snarl, whatever the shape of the" << endl
         << "                            difference." << endl
         << "      --ploidy-bed FILE     BED of CHROM START END PLOIDY setting ploidy per" << endl
         << "                            region, overriding -d/-R where an interval covers a" << endl
         << "                            site. CHROM is the contig as the output VCF spells" << endl
         << "                            it (chrX, not CHM13#0#chrX); intervals are 0-based" << endl
         << "                            half-open and must not overlap. Lets one run call a" << endl
         << "                            male chrX haploid outside the pseudoautosomal" << endl
         << "                            regions and diploid inside them, which per-contig" << endl
         << "                            ploidy cannot express and which otherwise needs two" << endl
         << "                            runs spliced together. Linkage and the mosaic break" << endl
         << "                            at each ploidy boundary, as they did at the splice" << endl
         << "  -R, --ploidy-regex RULES  use this comma-separated list of colon-delimited" << endl
         << "                            REGEX:PLOIDY rules to assign ploidies to contigs" << endl
         << "                            not visited by the selected samples, or to all" << endl
         << "                            contigs simulated from if no samples are used." << endl
         << "                            Unmatched contigs get ploidy 2 (or that from -d)." << endl
         << "      --top-down            top-down nested calling with genotype propagation" << endl
         << "                            from parent to child snarls (writes LV/PS tags)" << endl
         << "      --bottom-up           bottom-up nested calling with snarl merging" << endl
         << "  -I, --chains              call chains instead of snarls (experimental)" << endl
         << "  -L, --cluster F           merge called alt alleles whose length-weighted" << endl
         << "                            similarity is >= F, so 1/2 of two effectively" << endl
         << "                            identical alleles becomes 1/1 [1.0; experimental]" << endl
         << "      --cluster-min-len N   only apply -L merging at sites whose core length" << endl
         << "                            -- the longest allele after stripping the prefix" << endl
         << "                            and suffix common to all alleles -- is >= N bp" << endl
         << "                            (0 = every site) [50]" << endl
         << "                            in a nested run (-A/--top-down) a merged parent" << endl
         << "                            disagrees with its own child records by design" << endl
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
    // On by default: phasing is the linkage layer's Viterbi path, and the layer is already on by
    // default wherever it can run. Declines with the layer rather than erroring, so a run without a
    // panel is quietly unphased instead of refusing to start.
    bool   phased_output = true;
    bool   phased_explicit = false;
    string mosaic_out;
    string anchors_out;
    AnchorParams anchor_params;
    // Fill a gap no panel haplotype can be carried across with the reference, so a strand stays one
    // walk. On by default: it is the only path contiguous across such a region -- on chr20 all 37
    // remaining boundaries -- and the fill is marked `ref` in the file rather than passed off as an
    // evidenced haplotype. Off leaves the gap, for a consumer that would rather see a discontinuity
    // than an assertion.
    bool mosaic_patch_gaps = true;
    // Keep a switch inside a nested chain as its own segment. Off merges across it: fewer
    // switches, still one contiguous walk, and the parent's route through the child snarl.
    bool mosaic_keep_nested = true;
    // Carry the flanking haplotype through a stretch the panel cannot explain rather than
    // breaking the path there. On by default; rare, and it keeps threads contiguous.
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
    // On by default, and only ever armed where symbolic nested calling is -- which in practice
    // means --read-likelihood. Every other calling path declines it rather than erroring, so
    // the support-based default's output is untouched.
    bool atomize_blocks = true;
    bool atomize_explicit = false;
    // On by default under --read-likelihood, where it was measured: genome-wide it takes SNV F1
    // from 0.9752 to 0.9833 and SV F1 from 0.5134 to 0.5467 at no runtime or memory cost. It
    // declines rather than errors where its preconditions are absent, the way --linkage-weight does.
    bool nested_calling = true;
    bool nested_explicit = false;
    bool call_chains = false;
    bool all_snarls = false;
    size_t min_allele_len = 0;
    size_t max_allele_len = numeric_limits<size_t>::max();
    bool show_progress = false;

    // Nested calling option (for use with --top-down)
    bool star_allele = false;

    // Post-genotyping alt-allele merging, mirroring vg deconstruct's -L/--cluster-min-len
    double cluster_threshold = 1.0;
    // see the note in deconstruct_main.cpp: at 0 the similarity is dominated by sequence every
    // allele shares, so -L merges small variants it should leave alone.  50 is the SV size cutoff.
    int64_t cluster_min_allele_len = 50;
    bool cluster_min_len_set = false;

    // Read-level genotyping options. Opt-in: without --read-likelihood none of
    // this code runs and the default caller is unchanged.
    bool read_likelihood = false;
    string gam_filename;
    string gaf_filename;
    string dump_likelihoods_filename;
    string gam_index_filename;
    string gaf_base_filename;
    string gbz_base_filename;
    string gaf_base_binary = "gbz-base";
    // 0 means "let the backend choose": the two on-demand backends want different
    // windows, because a GAF-Base query is a process spawn where a .gai group scan is
    // a seek. See DEFAULT_GAM_INDEX_WINDOW / DEFAULT_GAF_BASE_WINDOW below.
    size_t read_window_size = 0;
    bool no_mismap_term = false;
    bool no_share_quality = false;
    double depth_quality = 0.0;
    // The read scorer's gap penalties. Exposed because they are the only primitive the
    // quality-adjusted scorer does NOT override -- score_gap is absent from
    // QualAdjAlignmentScorer's override list, so a gap is scored quality-blind -- and
    // because at the shipped 6/1 a 1 bp indel difference is worth 7 units, which at
    // log base 1.3833 nats/unit is a likelihood ratio of 6.2e-5: past the -ln(0.02)/1.3833
    // = 2.83-unit window the mismap floor leaves visible, so every indel difference of a
    // single base is a saturated vote whatever the base qualities say. That is the right
    // answer for a 0.1%-substitution read and the wrong one for a read whose modal error
    // IS a single-base homopolymer indel.
    int gap_open = default_gap_open;
    int gap_extend = default_gap_extension;
    // Charge a gap's open against the read's own confidence at it. Off by default because it moves
    // genotypes, and because the functional form is still being chosen by measurement.
    // Read-backed phasing. Off by default: with it off the phase is the panel's, byte for byte.
    bool read_phasing = false;
    bool read_phasing_explicit = false;
    bool regenotype = false;
    bool regenotype_explicit = false;
    RegenotypeParams regenotype_params;
    // Two barrier passes, i.e. ONE correction round, and that is a measured default rather than a
    // cautious one. The iteration is implemented and runs to a fixed point when there is one;
    // on chr20 there is not -- it enters a period-3 limit cycle at round 7 -- and eleven rounds
    // score WORSE than one: ALL F1 0.95136 against 0.95150, SV 0.56234 against 0.56577. Raise it
    // to watch the iteration; leave it here to get the answer.
    size_t regenotype_passes = 2;
    string regenotype_ledger;
    ReadPhasingParams read_phasing_params;
    // Which of the preset's values the user set for themselves. The preset is applied after the
    // whole option loop, so `--preset ont --gap-open 3` and `--gap-open 3 --preset ont` mean the
    // same thing: an explicit flag always wins, whichever side of the preset it is written on.
    string preset;
    bool gap_open_explicit = false, gap_extend_explicit = false, mismap_min_explicit = false;
    double min_confidence = 0.0;
    double linkage_weight = 2.0;
    /// Whether a weight was asked for, as opposed to inherited from the default. The two must
    /// behave differently where linkage is impossible: an explicit request has to fail loudly, and
    /// the default has to decline quietly, or every run without -z would error.
    bool linkage_weight_explicit = false;
    double linkage_scale = 10000.0;
    double linkage_freq_prior = 5.0;
    bool flat_mixture = false;
    double depth_weight = 0.1;
    bool depth_count_raw = false;
    double max_mismap_prob = 0.7;
    double min_mismap_prob = 0.02;
    int read_min_mapq = 0;

    // constants
    const size_t avg_trav_threshold = 50;
    const size_t avg_node_threshold = 50;
    const size_t min_depth_bin_width = 50;
    const size_t max_depth_bin_width = 50000000;
    const double depth_scale_fac = 1.5;
    const size_t max_yens_traversals = traversals_only ? 100 : 50;
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
    constexpr int OPT_MISMAP_MAX = 1019;
    constexpr int OPT_MISMAP_MIN = 1020;
    constexpr int OPT_NO_SHARE_QUALITY = 1021;
    constexpr int OPT_FLAT_MIXTURE = 1023;
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
    constexpr int OPT_REGENOTYPE = 1082;
    constexpr int OPT_NO_REGENOTYPE = 1087;
    constexpr int OPT_REGENO_TEMPER = 1083;
    constexpr int OPT_REGENO_PASSES = 1084;
    constexpr int OPT_REGENO_SHUFFLE = 1085;
    constexpr int OPT_REGENO_LEDGER = 1086;
    constexpr int OPT_MIN_CONFIDENCE = 1042;
    constexpr int OPT_PLOIDY_BED = 1043;
    constexpr int OPT_NESTED = 1044;
    constexpr int OPT_NO_NESTED = 1045;
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
    // Hoisted out of the getopt loop so the requires---read-likelihood scan below can
    // resolve abbreviated long options exactly the way getopt_long does.
    static const struct option long_options[] = {
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
        {"no-phased", no_argument, 0, OPT_NO_PHASED},
        {"gaf", no_argument, 0, 'G'},
        {"traversals", no_argument, 0, 'T'},
        {"trav-padding", required_argument, 0, 'M'},
        {"legacy", no_argument, 0, OPT_LEGACY},
        {"top-down", no_argument, 0, OPT_TOP_DOWN},
        {"bottom-up", no_argument, 0, OPT_BOTTOM_UP},
        {"atomize-blocks", no_argument, 0, OPT_ATOMIZE_BLOCKS},
        {"no-atomize-blocks", no_argument, 0, OPT_NO_ATOMIZE_BLOCKS},
        {"read-likelihood", no_argument, 0, OPT_READ_LIKELIHOOD},
        {"gam", required_argument, 0, OPT_GAM},
        {"gaf-reads", required_argument, 0, OPT_GAF},
        {"dump-likelihoods", required_argument, 0, OPT_DUMP_LIKELIHOODS},
        {"no-mismap-term", no_argument, 0, OPT_NO_MISMAP_TERM},
        {"mismap-max", required_argument, 0, OPT_MISMAP_MAX},
        {"mismap-min", required_argument, 0, OPT_MISMAP_MIN},
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
        {"regenotype", no_argument, 0, OPT_REGENOTYPE},
        {"no-regenotype", no_argument, 0, OPT_NO_REGENOTYPE},
        {"regeno-temper", required_argument, 0, OPT_REGENO_TEMPER},
        {"regeno-passes", required_argument, 0, OPT_REGENO_PASSES},
        {"regeno-shuffle", no_argument, 0, OPT_REGENO_SHUFFLE},
        {"regeno-ledger", required_argument, 0, OPT_REGENO_LEDGER},
        {"min-confidence", required_argument, 0, OPT_MIN_CONFIDENCE},
        {"linkage-weight", required_argument, 0, OPT_LINKAGE_WEIGHT},
        {"linkage-scale", required_argument, 0, OPT_LINKAGE_SCALE},
        {"linkage-prior", required_argument, 0, OPT_LINKAGE_FREQ_PRIOR},
        {"enumerate-support", no_argument, 0, OPT_ENUMERATE_SUPPORT},
        {"phased", no_argument, 0, OPT_PHASED},
        {"mosaic-out", required_argument, 0, OPT_MOSAIC_OUT},
        {"anchors-out", required_argument, 0, OPT_ANCHORS_OUT},
        {"anchors-het-only", no_argument, 0, OPT_ANCHORS_HET_ONLY},
        {"anchors-leaf-only", no_argument, 0, OPT_ANCHORS_LEAF_ONLY},
        {"anchors-reads", required_argument, 0, OPT_ANCHORS_MIN_READS},
        {"anchors-min-gqn", required_argument, 0, OPT_ANCHORS_MIN_GQN},
        {"anchors-min-q", required_argument, 0, OPT_ANCHORS_MIN_SCORE},
        {"anchors-keep-off-call", no_argument, 0, OPT_ANCHORS_KEEP_OFF_CALL},
        {"anchors-end-new", required_argument, 0, OPT_ANCHORS_END_PIN_MIN_NEW},
        {"mosaic-patch-gaps", no_argument, 0, OPT_MOSAIC_PATCH},
        {"no-mosaic-patch-gaps", no_argument, 0, OPT_NO_MOSAIC_PATCH},
        {"no-mosaic-nested", no_argument, 0, OPT_NO_MOSAIC_NESTED},
        {"mosaic-break-unexplained", no_argument, 0, OPT_MOSAIC_BREAK_UNEXPL},
        {"read-min-mapq", required_argument, 0, OPT_READ_MIN_MAPQ},
        {"gam-index", required_argument, 0, OPT_GAM_INDEX},
        {"gaf-base", required_argument, 0, OPT_GAF_BASE},
        {"gbz-base", required_argument, 0, OPT_GBZ_BASE},
        {"gaf-base-binary", required_argument, 0, OPT_GAF_BASE_BINARY},
        {"read-window", required_argument, 0, OPT_READ_WINDOW},
        {"chains", no_argument, 0, 'I'},
        {"cluster", required_argument, 0, 'L'},
        {"cluster-min-len", required_argument, 0, OPT_CLUSTER_MIN_LEN},
        // deprecated: shipped through v1.76 as an accepted no-op.  Kept accepted (and absent
        // from the helptext, which check_options.py allows) so pipelines carrying it do not die
        // on an unrecognized option.  Remove after one release.
        {"cluster-post", no_argument, 0, OPT_CLUSTER_POST},
        {"star-allele", no_argument, 0, 'Y'},
        {"threads", required_argument, 0, 't'},
        {"progress", no_argument, 0, OPT_PROGRESS },
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    while (true) {


        int option_index = 0;

        c = getopt_long (argc, argv, "k:Be:b:m:v:aAc:C:f:i:s:r:g:zN:Op:P:S:o:l:d:R:GTM:IL:Yt:h?",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

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
                logger.error() << "--regeno-temper must be >= 0; it is a temper, not a switch."
                               << " 0 reproduces --no-regenotype exactly" << endl;
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
        case OPT_READ_MIN_MAPQ:
            read_min_mapq = parse<int>(optarg);
            break;
        case OPT_GAM_INDEX:
            gam_index_filename = require_exists(logger, optarg);
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

    // Block emission is on by default, so these checks must DECLINE rather than refuse when the
    // setting is implicit -- refusing would break `vg call -a` and `--legacy` for everyone, on a
    // flag they never passed. Asked for by name, they still refuse, which is the --nested pattern.
    //
    // Checked here rather than after the graph is loaded, because these are option-compatibility
    // facts and making a user wait for a 22 GB load to be told the combination is invalid is waste.
    if (!preset.empty()) {
        // Applied after the option loop, so an explicit flag wins wherever it is written.
        //
        // `ont` is fitted on chr20 of HG002 at 43x against the T2T-Q100 benchmark and validated on
        // chr6 of the same read set and at 30x; it is not a default because the same values cost
        // 150 bp reads 0.0037 indel F1. See the gap-penalty help above for why a long read wants a
        // different gap scale, and note the two values are near-additive rather than alternatives.
        if (preset == "ont") {
            if (!read_phasing_explicit) {
                // On for long reads because the reads carry the answer and the panel does not: at
                // 33 kb against a 349 bp median het spacing, 99.96% of adjacent het pairs share a
                // read. chr20 switch error 3.7942% -> 0.5180% and chr6, held out at these same
                // values, 3.1849% -> 0.3447%, with the genotypes provably unmoved. Not a global
                // default: at 151 bp most adjacent pairs are not spanned by one read, and that
                // case is unmeasured.
                read_phasing = true;
            }
            if (!gap_open_explicit) {
                gap_open = 1;
            }
            if (!gap_extend_explicit) {
                gap_extend = 1;
            }
            if (!mismap_min_explicit) {
                min_mismap_prob = 0.05;
            }
            if (!regenotype_explicit) {
                // The reads' phase decides the genotype, not only the order of a settled pair.
                // chr20 ALL F1 0.94477 -> 0.95151 and indel 0.81568 -> 0.83725, with precision
                // AND recall both up; chr6, held out and fitting its own temper, 0.95342 ->
                // 0.95960 and 0.83781 -> 0.86019. The 6,350 sites it moves on chr20 run a 40.6%
                // false-positive rate against a 6.7% background, so it aims at what is broken,
                // and the shuffled-sign control -- same |Lambda|, phase destroyed -- loses 0.107
                // F1, so the gain is the phase and not a change in peakedness.
                //
                // Costs +13.2% wall and about +110 MB of retention. Requires read phasing, which
                // this preset also turns on; see the decline below for the one combination that
                // leaves it off.
                regenotype = true;
            }
        } else {
            cerr << "error [vg call]: unknown --preset '" << preset << "'; known presets: ont"
                 << endl;
            return 1;
        }
    }
    if (atomize_blocks && genotype_snarls) {
        // -a's record set is meant to be sample-independent -- one line per snarl, and the harness
        // reads a `-a -A` run as a snarl inventory. A block list is a function of the called
        // haplotypes, so it cannot be sample-independent. Fatal only if asked for by name;
        // otherwise the default steps aside, because -a must keep working without a flag.
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
    // -I is a sharper incompatibility than the three above, which merely pick a different calling
    // path. It calls a *fabricated* snarl spanning a whole chain piece, and whether the SnarlManager
    // recognises that snarl depends on how `break_chain` packed the chain: a piece holding one snarl
    // is value-identical to the managed snarl and resolves, a piece holding several does not. So
    // `symbolic_site_resolvable` -- the gate governing symbolic collapsing, is_symbolically_reference,
    // nested descent and block emission alike -- passes on an arbitrary subset of sites.
    //
    // Measured on the 18_vg_call.t `x` fixture: 70 records with a 4 bp longest ALT by default, 1
    // record with a 997 bp ALT under -I. That is the swallowed-variant shape nested calling exists
    // to undo, and applying the cure to only some sites is not a defensible middle.
    if (atomize_blocks && call_chains) {
        if (atomize_explicit) {
            logger.error() << "--atomize-blocks cannot be combined with -I/--chains: a chain piece "
                           << "is called as a fabricated snarl, and only a piece holding a single "
                           << "snarl resolves, so decomposition would apply to an arbitrary subset "
                           << "of sites" << endl;
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

    // -L/--cluster validation, deliberately ahead of the graph load below: it is all command-line
    // state, and a typo'd "-L 5" should not cost a multi-gigabyte graph load before it is rejected.
    // vg deconstruct clamps out-of-range values; we reject them, since -L 5 is a plausible typo for
    // -L 0.5 and clamping would silently disable merging altogether.
    if (cluster_threshold < 0.0 || cluster_threshold > 1.0 || std::isnan(cluster_threshold)) {
        logger.error() << "-L/--cluster threshold must be in range [0.0, 1.0]" << endl;
    }
    // only for an explicit --cluster-min-len: the default is nonzero, so testing the value alone
    // would warn on every run that does not pass -L
    if (cluster_min_len_set && cluster_min_allele_len > 0 && cluster_threshold >= 1.0) {
        logger.warn() << "--cluster-min-len has no effect without -L (cluster threshold < 1.0)" << endl;
    }
    // -L merges in VCFOutputCaller::emit_variant, which these paths never reach: the VCF genotyper
    // builds its records by hand and the GAF path has no VCF at all.
    if (cluster_threshold < 1.0 && !vcf_filename.empty()) {
        logger.error() << "-L/--cluster cannot be used when genotyping a VCF (-v)" << endl;
    }
    if (cluster_threshold < 1.0 && (gaf_output || traversals_only)) {
        logger.error() << "-L/--cluster cannot be used with GAF output (-G/-T)" << endl;
    }
    // the ratio caller's QUAL and its XADL/lowxadl describe a het that no longer exists after a merge
    if (cluster_threshold < 1.0 && ratio_caller) {
        logger.error() << "-L/--cluster cannot be used with the ratio caller (-B)" << endl;
    }
    // NestedFlowCaller represents a child snarl as a Visit carrying a Snarl rather than a node, which
    // the clusterer has no handle for, so the merge would be inert at exactly the nested sites
    // --bottom-up exists for.  (Rejected on master too; the check has only moved ahead of the load.)
    if (cluster_threshold < 1.0 && bottom_up) {
        logger.error() << "-L/--cluster cannot be used with --bottom-up mode" << endl;
    }
    // -Y writes "*" in a child record to mean "an upstream deletion covers this site".  If the merge
    // absorbs the parent allele that deletion came from, the "*" refers to nothing in the file.  That
    // is a malformed record, not merely a lossy one -- unlike the nested case, which is allowed: in a
    // nested run a merged parent deliberately disagrees with its own child records, the parent giving
    // the collapsed view of a large variant and the children the precise one.  MAT records it.
    if (cluster_threshold < 1.0 && star_allele) {
        logger.error() << "-L/--cluster cannot be used with -Y/--star-allele" << endl;
    }
    // LegacyCaller has its own traversal finder and support model; -L has never been run through it
    if (cluster_threshold < 1.0 && legacy) {
        logger.error() << "-L/--cluster cannot be used with the legacy caller (--legacy)" << endl;
    }
    // the merge needs two distinct called ALT alleles, which a haploid genotype can never have.
    // -R/--ploidy-regex sets ploidy per contig and we cannot tell here which contigs will be
    // called, so only -d is worth warning about.
    if (cluster_threshold < 1.0 && ploidy == 1) {
        logger.warn() << "-L/--cluster has no effect at ploidy 1 (-d 1)" << endl;
    }
    // NestedFlowCaller's Snarl-carrying Visits also break the GAF emitters: --bottom-up -T aborts in
    // to_mapping, and --bottom-up -G emits a header and no records.
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
            // Under --read-likelihood this is decided below rather than suggested, so the
            // hint would either be wrong or be immediately contradicted.
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

    // An index without the thing it indexes is always a mistake, whatever else was
    // passed. Checked before the read-likelihood validation so the message is
    // deterministic rather than depending on which error is reached first.
    if (!gam_index_filename.empty() && gam_filename.empty()) {
        logger.error() << "--gam-index requires --gam" << endl;
    }

    // Validation: every option below only means something under --read-likelihood, so passing
    // one without it is a command line that does not do what it says. Silently ignoring them is
    // how a run gets analysed under the wrong assumptions -- and the inconsistency was real
    // before this check existed: an explicit --linkage-weight was refused while --depth-term,
    // --flat-mixture and the rest were accepted and dropped.
    //
    // Read from argv rather than from an `explicit` bool per option. Most of these have non-zero
    // defaults, so "was it set" is not recoverable from the parsed value, and eighteen more
    // tracking bools would be worse than one scan. Tokens are matched whole, or up to an '=',
    // so a filename containing one of these strings cannot trigger it.
    if (!read_likelihood) {
        static const vector<string> read_likelihood_only = {
            "--gam", "--gaf-reads", "--gam-index", "--gaf-base", "--gbz-base",
            "--gaf-base-binary", "--read-window", "--read-min-mapq", "--no-mismap-term",
            "--depth-term", "--depth-count-raw", "--linkage-weight", "--linkage-scale",
            "--linkage-prior", "--depth-quality", "--min-confidence", "--flat-mixture",
            "--gap-open", "--gap-extend", "--preset",
            "--read-phasing", "--no-read-phasing", "--phase-min-q", "--phase-break",
            "--regenotype", "--no-regenotype", "--regeno-temper", "--regeno-passes",
            "--regeno-shuffle", "--regeno-ledger",
            "--phase-relink",
            "--phase-hang", "--phase-prior", "--phase-cap",
            "--no-share-quality",
            "--mismap-max", "--mismap-min", "--dump-likelihoods", "--enumerate-support",
            "--phased", "--no-phased", "--mosaic-out", "--anchors-out", "--anchors-het-only",
            "--anchors-leaf-only", "--anchors-reads", "--anchors-min-gqn",
            "--anchors-min-q", "--anchors-keep-off-call", "--anchors-end-new"};
        vector<string> offenders;
        for (int i = 1; i < argc; ++i) {
            string arg(argv[i]);
            if (arg.size() < 3 || arg[0] != '-' || arg[1] != '-') {
                continue;
            }
            size_t eq = arg.find('=');
            string name = (eq == string::npos ? arg : arg.substr(0, eq)).substr(2);
            // Resolved the way getopt_long resolves it: an exact match, or an unambiguous prefix.
            // Matching whole tokens only let every abbreviated spelling through -- "--mosaic" set
            // --mosaic-out and sailed past this scan, silently reproducing the very ignored-flag
            // failure it exists to close.
            const char* resolved = nullptr;
            bool ambiguous = false;
            for (const struct option* o = long_options; o->name != nullptr; ++o) {
                string oname(o->name);
                if (oname == name) {
                    resolved = o->name;
                    ambiguous = false;
                    break;
                }
                if (oname.compare(0, name.size(), name) == 0) {
                    ambiguous = resolved != nullptr;
                    resolved = o->name;
                }
            }
            if (resolved == nullptr || ambiguous) {
                continue;   // unknown or ambiguous: getopt itself already rejected it
            }
            string full = string("--") + resolved;
            for (const string& flag : read_likelihood_only) {
                if (full == flag) {
                    offenders.push_back(flag);
                    break;
                }
            }
        }
        if (!offenders.empty()) {
            stringstream joined;
            for (size_t i = 0; i < offenders.size(); ++i) {
                joined << (i ? ", " : "") << offenders[i];
            }
            logger.error() << joined.str()
                           << (offenders.size() == 1 ? " only applies" : " only apply")
                           << " to --read-likelihood, which was not given" << endl;
        }
    }

    // Validation: --read-likelihood needs reads, and cannot be combined with the
    // support-based model selection flags. Failing here rather than later keeps a
    // read-free "read-level" genotyping run from silently happening.
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

    // Default allele enumeration for --read-likelihood on a GBZ that carries a haplotype
    // panel: enumerate the traversals the panel actually spells, rather than every
    // traversal the support flow permits.
    //
    // Measured on HG002 against HPRC graphs (chr20 and chr6, 4- and 34-haplotype), as
    // haplotype enumeration against support enumeration under this caller: better small
    // variant F1 in all four (0.9487 -> 0.9507, 0.9513 -> 0.9645, 0.9583 -> 0.9602,
    // 0.9588 -> 0.9689), and better SV F1 in three of four (+0.0352, +0.0144, +0.0269),
    // the exception being chr20 4-haplotype at -0.0018, which is under two events on a
    // 765 event benchmark and below what it resolves. It also drops the pack file
    // requirement, since nothing then consults support.
    //
    // Not made the default for the support caller, where the same comparison loses SV F1
    // on all four datasets (0.4954 -> 0.4930, 0.4535 -> 0.4391, 0.5490 -> 0.5478,
    // 0.4944 -> 0.4881); flipping it there would regress existing callers by as much as
    // 0.0144, so -z stays opt-in outside this mode.
    //
    // Caveat worth knowing before trusting the default: enumeration from a panel cannot
    // spell an allele no haplotype carries, so the ceiling is the panel's content. The
    // numbers above are one sample against a panel that excludes it but represents its
    // variation well, and a sample poorly represented in the panel would fare worse.
    // --enumerate-support is the way out.
    if (read_likelihood && !gbz_paths && !enumerate_support && gbz_graph &&
        gbwt_filename.empty() && vcf_filename.empty()) {
        // Only worth it if there is a panel to enumerate. A GBZ always has a GBWT, but it
        // may hold nothing but reference paths, and enumerating from that would offer the
        // reference allele and nothing else: silently near-zero alt recall, which is far
        // worse than the support enumeration it replaced. One haplotype is not an error
        // but is too thin to choose automatically; -z still forces it.
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
                // A pack is only ever consulted by support enumeration, so under this default
                // it is dead weight. Passing one is a fair signal that support enumeration was
                // what was wanted, and taking it without using it would leave the run quietly
                // doing something other than what the command line asks for.
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

    if (max_mismap_prob <= 0.0 || max_mismap_prob >= 1.0) {
        logger.error() << "--mismap-max must be in (0, 1)" << endl;
    }
    if (min_mismap_prob <= 0.0 || min_mismap_prob > max_mismap_prob) {
        logger.error() << "--mismap-min must be in (0, --mismap-max]" << endl;
    }

    // --gbz-base only says where to point the query; it means nothing without the read
    // database that is being queried.
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

    // A gref cover changes what "nested" can mean. Chains the linear reference does not cross are
    // skipped by default because their records have no REF and no POS -- but a gref fragment IS a
    // contig, with its own coordinates, covering exactly the inserted sequence those chains live in.
    // So selecting one is the signal to descend into them: it is the thing that makes them
    // reportable, and asking for the cover and then not using it has no other reading.
    //
    // Prefix-selected covers are the normal case (`-P gref_CHM13#0#chr20_`), so this tests the
    // resolved path list rather than the flags.
    for (const string& ref_path : ref_paths) {
        // A FRAGMENT, not merely a gref-derived name. The gref copy of a base contig is just the
        // reference under another name and makes nothing new reportable, so selecting it alone must
        // not turn the descent on -- on a graph that ships only the gref view of its reference,
        // which is the ordinary case, that would enable it for every run including the control.
        if (GrefCover::is_gref_name(ref_path)) {
            enable_off_reference_nesting();
            if (show_progress) {
                logger.info() << "gref reference selected: descending into chains the reference "
                              << "does not cross, and reporting them against their gref contig"
                              << endl;
            }
            break;
        }
    }

    // How deep in non-reference sequence each selected gref contig sits, so INFO/CH can say so
    // without depending on which of a record's ancestors happened to be emitted.
    //
    // The cover is node-disjoint, so the node a fragment's first node hangs off belongs to exactly
    // one gref interval, and that interval is its parent -- which is `assign_nesting()`'s rule
    // ("walk up until reaching a snarl whose boundary nodes are owned by some gref interval")
    // answered locally instead of over the snarl tree. A fragment off the base reference is 1, one
    // inside a level-1 fragment is 2, and the base contig itself is 0.
    //
    // Keyed by LOCUS, because that is what reaches the VCF's CHROM column.
    map<string, int> gref_levels;
    if (off_reference_nesting_enabled() && graph != nullptr) {
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
                return 0;   // a base contig, or the gref copy of one: the outermost system
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
                // Its neighbours are in no interval -- short runs the length filter dropped. It is
                // still a fragment, so it is still at least one layer in.
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
        // After the decomposition, not after the finder's constructor. The finder is cheap and the
        // decomposition is not: on chr20 it is 46 s of a 197 s run, and printing "Computed snarls"
        // in front of it left the log claiming the step was done while a single thread spent a
        // quarter of the wall clock on it.
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
    // Read-level genotyping can run without a pack file, but only when allele
    // enumeration does not need support either. GBWTTraversalFinder enumerates from
    // recorded haplotypes, so it needs none; FlowTraversalFinder is driven entirely
    // by node and edge weights, so it does.
    bool gbwt_enumeration = !gbwt_filename.empty() || gbz_paths;
    bool support_free = read_likelihood && pack_filename.empty();

    if (!pack_filename.empty() || support_free) {
        if (support_free) {
            // Nothing downstream will consult support: stand in a finder that
            // reports none rather than requiring a pack file for its own sake.
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
                // Make a packed traversal support finder (using cached version important for poisson caller)
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
            // Read-level genotyping. The pack file is still required, but only so
            // FlowTraversalFinder has node/edge weights to *enumerate* alleles
            // with; the genotyping itself uses no support at all.
            if (show_progress) logger.info() << "Loading reads for read-level genotyping" << endl;

            SiteReadFilter read_filter;
            read_filter.min_mapq = read_min_mapq;

            if (!gaf_base_filename.empty()) {
                // GAF-Base: reads are fetched per window by running gbz-base. A runtime
                // dependency on that binary, and no build dependency at all.
                //
                // The query needs a graph to resolve node IDs against. Default to the
                // graph vg call was given, which is usually the right one and is
                // certainly the right *graph*; but a GBZ-Base is random-access where a
                // plain GBZ is loaded in full on every query, so say so.
                string query_graph = gbz_base_filename.empty() ? graph_filename
                                                               : gbz_base_filename;
                if (read_window_size == 0) {
                    read_window_size = DEFAULT_GAF_BASE_WINDOW;
                }
                auto gaf_base_source = new GafBaseSiteReadSource(*graph, gaf_base_filename,
                                                                 query_graph, read_filter,
                                                                 read_window_size, 2,
                                                                 gaf_base_binary);
                read_source.reset(gaf_base_source);
                // Probe now, on the main thread. A missing binary or an unreadable
                // database is the user's setup, not a vg bug, so report it as an error
                // and exit rather than letting the exception out to the crash handler --
                // which would print a bug-report banner for "install gbz-base".
                try {
                    gaf_base_source->check_setup();
                } catch (const std::exception& e) {
                    logger.error() << e.what() << endl;
                }
                if (show_progress) {
                    logger.info() << "Using GAF-Base " << gaf_base_filename
                                  << " queried against " << query_graph
                                  << " via " << gaf_base_binary << endl;
                    if (gbz_base_filename.empty()) {
                        logger.info() << "Consider building a GBZ-Base ('gbz-base construct') and "
                                      << "passing --gbz-base: a plain GBZ is reloaded on every query"
                                      << endl;
                    }
                }
            } else if (!gam_index_filename.empty()) {
                // Indexed: reads are fetched per site, so memory is bounded by what a
                // site needs rather than by the size of the read set.
                if (read_window_size == 0) {
                    read_window_size = DEFAULT_GAM_INDEX_WINDOW;
                }
                read_source.reset(new IndexedGamSiteReadSource(gam_filename, gam_index_filename,
                                                               read_filter, read_window_size));
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
            likelihood_params.collect_anchors = anchor_params.enabled;
            likelihood_params.collect_read_phasing = read_phasing;

            likelihood_calculator.reset(new GraphAlignedAlleleLikelihoodCalculator(
                *graph, *snarl_manager, *read_source, *qual_scorer, *plain_scorer,
                likelihood_params));

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
                           << "enumeration (-g/--gbwt or -z/--gbz); otherwise a pack file is needed "
                           << "for the flow traversal finder's node and edge weights" << endl;
        }
        if (!vcf_filename.empty()) {
            // VCFTraversalFinder prunes alt paths on support before its brute-force
            // enumeration, so -v genuinely needs a pack file.
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

    // Per-region ploidy, if given. Applied to whichever caller was built: every one of them
    // derives from VCFOutputCaller, which is where the override map and its lookup live.
    if (!ploidy_bed_filename.empty()) {
        VCFOutputCaller* ploidy_target = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        if (ploidy_target == nullptr) {
            cerr << "error [vg call]: --ploidy-bed needs a caller that emits VCF" << endl;
            return 1;
        }
        ploidy_target->set_ploidy_regions(ploidy_bed_filename);
    }

    // Symbolic collapsing. A called traversal that takes the same route through a snarl as the
    // reference, differing only inside child chains, is the reference allele there; emitting it as
    // a long ALT is what buries the nested variants it contains.
    //
    // On by default under --read-likelihood, and declining rather than erroring elsewhere. It is
    // only measured there: every arm behind the numbers in doc/read-likelihood-genotyping.md runs
    // that model, and turning it on for the support-based caller would change the legacy default's
    // output on no evidence. An explicit --nested still works anywhere, as it did when it was opt-in.
    if (nested_calling && !nested_explicit && !read_likelihood) {
        nested_calling = false;
    }
    // -A and symbolic nested calling are two different ways to reach a child snarl, and running
    // both visits every nested snarl TWICE. `-A` sets RecurseAlways, so GraphCaller queues each
    // child as its own snarl and it lands at generation 0; the symbolic descent then retains the
    // same child as a generation-1 chain. Both hash the same print_snarl string, so both take the
    // same record_key -- and LinkageCollector::live_index returns the FIRST non-retracted entry
    // for a key, so which of the two describes the site depends on insertion order, which is
    // thread assignment. Measured on the nestblk fixture at -t 1: six VCF lines in three identical
    // pairs, the top-level record absent, and nested chains called `1|1` where --nested gives
    // `1|.`.
    //
    // --top-down resolves the same conflict the other way, taking RecurseNever with the note that
    // "FlowCaller handles recursion internally". -A cannot: its contract is that every snarl is an
    // INDEPENDENT record, which is precisely a statement that the nesting relationship is not
    // used, and the descent's pruning would make "all snarls" untrue.
    //
    // -A is from 2021 and predates both the descent and its default, so this restores what -A did
    // before nested calling became automatic under --read-likelihood.
    // Refused rather than declined: `--regenotype` without read phasing is not a weaker version
    // of the feature, it is the identity. `Lambda` comes from the phasing chain, and with no chain
    // every read's strand log-odds is zero, every tilted weight is the site's own weight, and the
    // correction is exactly zero at every site. Running it would burn a pass to change nothing
    // while the flag says otherwise.
    // Refused under --top-down, and this is the one place the nested story genuinely breaks.
    //
    // The read-likelihood descent hands a child `nullptr` for its traversal sets, so the child
    // enumerates its own candidates and a parent moving changes only its PLOIDY -- which the
    // barrier re-derives. `--top-down` is the other path: it builds `ChildTraversalSets` from
    // `trav_genotype`, the parent's CALLED genotype, and those sets decide the child's ploidy,
    // are merged into its candidate list, and supply its pseudo-reference. Move the parent there
    // and the child was never scored against the alleles the parent now carries -- which no
    // amount of re-resolving repairs, because the columns were never in `rel`.
    //
    // --bottom-up is refused for a duller reason: it builds NestedFlowCaller, which is not a
    // FlowCaller, so `render_retained_records` never runs and the flag would do nothing at all.
    if (regenotype && (top_down || bottom_up)) {
        cerr << "error [vg call]: --regenotype cannot be combined with "
             << (top_down ? "--top-down" : "--bottom-up") << "; "
             << (top_down
                 ? "--top-down derives each child's candidate traversals from its parent's called"
                   " genotype, so moving a parent leaves children scored against alleles it no"
                   " longer carries, and nothing downstream can repair that"
                 : "--bottom-up uses a caller that does not run the render pass this needs, so the"
                   " flag would be silently inert")
             << endl;
        return 1;
    }
    if (regenotype && !read_phasing) {
        if (regenotype_explicit) {
            cerr << "error [vg call]: --regenotype needs --read-phasing; the correction is"
                 << " computed from the phasing chain, so without one it is the identity" << endl;
            return 1;
        }
        // Armed by the preset, and the same preset's phasing has been switched off by hand --
        // `--preset ont --no-read-phasing`. Declining is right: the user turned off the thing this
        // is computed from, and erroring on a flag they never typed would make a documented
        // combination fail.
        regenotype = false;
    }
    if (nested_calling && all_snarls) {
        if (nested_explicit) {
            cerr << "error [vg call]: -A/--all-snarls calls every snarl independently while"
                 << " --nested calls children through their parent; they are alternatives, not a"
                 << " combination" << endl;
            return 1;
        }
        nested_calling = false;   // the default declines, as it does for the support caller
    }
    {
        // The confidence threshold also to the output layer: a record whose GQN the linkage layer
        // re-derives has to be re-labelled against the same number the per-site emission used,
        // rather than cleared to PASS regardless.
        VCFOutputCaller* confidence_target = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        if (confidence_target != nullptr) {
            confidence_target->set_linkage_min_confidence(min_confidence);
        }
    }
    if (nested_calling) {
        VCFOutputCaller* nested_target = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        if (nested_target == nullptr) {
            if (nested_explicit) {
                cerr << "error [vg call]: --nested needs a caller that emits VCF" << endl;
                return 1;
            }
            nested_calling = false;   // the default declines
        }
    }
    if (atomize_blocks) {
        // Refused, not silently declined. Each of these would produce a record set the flag's own
        // contract does not describe, and a silent decline is how a measurement ends up being of
        // something other than what it was labelled.
        //
        // -a/--genotype-snarls is the sharpest: its record set is meant to be sample-independent,
        // one line per snarl, and the harness reads a `-a -A` run as a snarl inventory. A block
        // list is a function of the called haplotypes, so it cannot be sample-independent. Note
        // that nested calling is ON by default under --read-likelihood and is NOT cleared by -a,
        // so this cannot be left to the nested-calling gate below.
        // The two purely-option refusals are checked at option-validation time above, so they
        // fire before a graph is loaded. What is left here depends on constructed state.
        if (!nested_calling) {
            if (atomize_explicit) {
                cerr << "error [vg call]: --atomize-blocks needs symbolic nested calling, which is"
                     << " on by default under --read-likelihood; pass --nested to enable it"
                     << " elsewhere" << endl;
                return 1;
            }
            atomize_blocks = false;   // the default declines, as --nested does
        }
        VCFOutputCaller* atomize_target =
            atomize_blocks ? dynamic_cast<VCFOutputCaller*>(graph_caller.get()) : nullptr;
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

    if (nested_calling && call_chains) {
        // Declined BEFORE the arming below rather than undone after it, and outside the
        // `if (!gaf_output)` block where the first version of this sat: -I exists for GAF output
        // (commit 8d34aceab, "call chains instead of snarls when writing GAF"), so a decline that
        // fires only for VCF misses the mode the flag is actually used in. Measured before the
        // move: `-I --nested` exited 1 while `-G -I --nested` exited 0 and descended.
        //
        // The reason is NOT "no site resolves". That was wrong, and measurably so. `break_chain`
        // packs a chain into pieces of up to 1,000 edges, and a piece holding exactly ONE snarl
        // yields a fabricated Snarl value-identical to the managed one, which `resolve_site`
        // accepts -- so collapsing and descent work there. On the `nestblk` fixture `-G -I` retains
        // 2 nested chains and emits 6 GAF lines against 4 with --no-nested.
        //
        // What is true is worse for being intermittent: whether a site resolves depends on how
        // `break_chain` happened to pack the chain, which is graph topology against a hard-coded
        // edge budget, not anything the data says. So nested calling silently works on some sites
        // and not others within one run, and no output distinguishes the two.
        if (nested_explicit) {
            cerr << "error [vg call]: --nested cannot be combined with -I/--chains, which calls a"
                 << " fabricated snarl spanning a chain piece; only a piece holding a single snarl"
                 << " resolves, so nested calling would apply to an arbitrary subset of sites"
                 << endl;
            return 1;
        }
        nested_calling = false;
        if (show_progress) {
            logger.info() << "Nested calling declines under -I/--chains: a chain piece is called as"
                          << " a fabricated snarl, and only a single-snarl piece resolves" << endl;
        }
    }

    if (nested_calling) {
        VCFOutputCaller* nested_target = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        nested_target->set_symbolic_collapsing(snarl_manager.get());
        // A nested site's ploidy comes from a parent genotype linkage can afterwards invalidate,
        // so descent asks the genotyper for both ploidies' answers at once, from the matrix it has
        // already built -- per call, via set_want_alt_ploidy, so nothing needs arming here. The
        // barrier then settles the ploidy and renders the record at it, which is why the second
        // answer is load-bearing rather than diagnostic and why the INFO/NGT2 field that used to
        // report it is gone: the caller acts on it instead of describing it.
    }

    // Owned here because write_variants(), at the very end of main, consumes the collector.
    unique_ptr<LinkageCollector> linkage_collector;
    vector<size_t> linkage_sequence_to_haplotype;

    string header;
    if (!gaf_output) {
        // Init The VCF       
        VCFOutputCaller* vcf_caller = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        assert(vcf_caller != nullptr);
        // Make sure we get the LV/PS tags with -A, --top-down, or --bottom-up -- and under a gref
        // cover, where a record's contig no longer says where it sits. On a linear reference every
        // record is on the contig and nesting is a detail; with a cover, a record on
        // gref_CHM13#0#chr20_4_alt is INSIDE an insertion, and LV/CH/PS are how a consumer says so.
        // The gref wiki's own filtering recipe is `bcftools view -i 'INFO/CH==1'`, and without this
        // it selected nothing.
        //
        // NOT on for every --read-likelihood run. Its descent emits nested records on the ordinary
        // path too and they are equally untagged, so there is a case for it -- but it would add INFO
        // to every nested record of every existing run, which is a change to make deliberately and
        // measure, not to fold into gref support.
        vcf_caller->set_nested(all_snarls || top_down || bottom_up
                               || off_reference_nesting_enabled());
        vcf_caller->set_gref_levels(std::move(gref_levels));
        vcf_caller->set_translation(translation.get());

        // Linkage pass. Needs both -z and --read-likelihood: the panel matrix comes from asking
        // which GBWT haplotypes traverse each allele, so without haplotype enumeration there is no
        // panel, and the emission is the read-likelihood model's own ln P(reads | G), so without
        // that there is nothing for the transition to reweight. Kept alive to the end of main
        // because write_variants() consumes it.
        //
        // The --read-likelihood half is asserted here rather than left to the dynamic_cast in
        // emit_variant. That cast already fails on the Poisson path, so the layer was inert there
        // by accident: calls came out byte-identical, but the setup ran and the diagnostic reported
        // "0 sites" on every -z run. Inert on purpose reads better than inert by luck, and it
        // survives someone giving the Poisson path a richer CallInfo later.
        // The layer declines first, so everything built on it can see whether it survived. This
        // used to run after the two checks below, which was harmless only while phasing was opt-in:
        // with phasing on by default, testing linkage_weight before the decline would refuse every
        // run without a panel instead of quietly emitting unphased calls.
        if (linkage_weight > 0.0 && !(gbwt_index != nullptr && read_likelihood)) {
            if (linkage_weight_explicit) {
                cerr << "error [vg call]: --linkage-weight needs haplotype enumeration (-z or -g) "
                     << "and --read-likelihood" << endl;
                return 1;
            }
            linkage_weight = 0.0;   // the default declines; only an explicit request is an error
        }
        if (!anchors_out.empty() && !read_likelihood) {
            // Always an explicit request -- a path was named -- so always an error rather than a
            // silent downgrade. The partition is the per-read/per-allele responsibility, which only
            // the read-likelihood caller computes.
            logger.error() << "--anchors-out needs --read-likelihood: the anchor partition is the "
                           << "read/allele likelihood matrix, which no other caller builds" << endl;
        }
        if (!mosaic_out.empty() && linkage_weight <= 0.0) {
            // Always an explicit request -- a path was named -- so always an error.
            logger.error() << "--mosaic-out needs the linkage model, which needs haplotype "
                           << "enumeration (-z or -g) and --read-likelihood" << endl;
        }
        if (phased_output && linkage_weight <= 0.0) {
            // Phasing is the linkage layer's Viterbi path, so without the layer there is no path to
            // emit. Asked for outright that is an error, because silently writing unphased output
            // would look like the flag had worked; on by default it declines the same way the layer
            // itself does.
            if (phased_explicit) {
                logger.error() << "--phased needs the linkage model, which needs haplotype "
                               << "enumeration (-z or -g) and --read-likelihood" << endl;
            }
            phased_output = false;
        }
        if (nested_calling && linkage_weight > 0.0 && !phased_output) {
            // Linkage on, phasing off: the one configuration nested descent cannot serve.
            //
            // A child's ploidy and strand come from its parent's *settled* genotype, and where the
            // linkage layer runs that genotype is only pinned down once the layer has phased it --
            // the settled allele pair lives in the phasing. Without it the choices are to descend
            // from the pre-linkage genotype, which is the incoherence this design exists to remove,
            // or to defer and then drop every child whose parent cannot be found.
            //
            // Declines or errors on the same rule the rest of this option layer uses: an explicit
            // --nested has to fail loudly, because silently dropping it would look like it worked,
            // while the default turns itself off and says so. Before this it printed a warning and
            // carried on, which left a run differing from its neighbour only by --no-phased quietly
            // emitting nested records at a ploidy their own parents contradicted.
            if (nested_explicit) {
                cerr << "error [vg call]: --nested needs --phased where the linkage model runs, "
                     << "because a nested site's ploidy and strand come from its parent's phased "
                     << "genotype" << endl;
                return 1;
            }
            // Undone on the caller, not just in this flag. The caller was configured for nested
            // calling further up -- symbolic collapsing on -- because whether phasing would run was
            // not settled then. Clearing the flag alone would leave symbolic collapsing armed,
            // which is what *enables* descent, so the run would still descend, inline, from
            // genotypes linkage then rewrote: precisely the configuration being refused.
            nested_calling = false;
            VCFOutputCaller* nested_target = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
            if (nested_target != nullptr) {
                // Disarming collapsing is what disarms descent, and descent is the only thing that
                // ever requests the alternate ploidy (set_want_alt_ploidy is per call), so no
                // second knob needs unwinding here.
                nested_target->set_symbolic_collapsing(nullptr);
            }
            if (show_progress) {
                logger.info() << "Nested calling declines under --no-phased: a nested site's ploidy "
                              << "and strand come from its parent's phased genotype" << endl;
            }
        }
        if (linkage_weight > 0.0) {
            // A haplotype is (sample, phase); a GBWT sequence is one orientation of one path, and
            // a haplotype in several fragments owns several paths. Collapse all of them onto one
            // index, or the panel would double-count a fragmented haplotype and treat its pieces
            // as independent evidence.
            const gbwt::Metadata& meta = gbwt_index->metadata;
            // Which base samples the graph carries in their own right, so a gref copy of one can be
            // told from a gref copy of one that is absent.
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
                // The cover's FRAGMENTS are not a haplotype, and they reach this loop because the
                // panel is built from every GBWT path regardless of sense. They all carry sample
                // "gref_<REF>" and phase 0, so they would collapse onto ONE panel index: an
                // individual stitched greedily from whichever donor had the longest uncovered run at
                // each site. Letting Li-Stephens copy from that chimera moved 8.82% of chr20's
                // genotypes and lost 298 sites, and it was the ONLY cause -- the 13,711
                // off-reference chains gref makes reportable left the contig byte-identical.
                //
                // The gref COPY OF A BASE CONTIG is a different thing: it is the reference, renamed,
                // one contig with one walk. It is kept unless the base path is also in the graph, so
                // the reference sits in the panel exactly ONCE either way. Both layouts occur --
                // `GrefCover::apply()` documents writing the reference twice, but a graph can ship
                // carrying only the gref view, and then dropping it would take a real haplotype out
                // of the panel. Measured on hprc-v2.1-mc-chm13-eval.gref.HG002.hap32.gbz, which has
                // no CHM13 path at all: keeping it is the difference between panel 34 and 33.
                //
                // WILDCARD rather than skipping: the vector defaults to 0, so a gref sequence left
                // unwritten would reach panel_alleles as haplotype 0 and write the cover's allele
                // into a real haplotype's slot. panel_alleles' `hap < out.size()` guard drops it.
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
            // A pair of haplotypes is the state, so two is the minimum that can express anything:
            // below that every state is the wildcard and the posterior collapses back to the
            // per-site likelihood. Harmless, but running an HMM over an empty panel and reporting
            // it as active is worse than declining and saying so.
            if (hap_index.size() < 2) {
                cerr << "warning [vg call]: linkage disabled -- the GBWT carries "
                     << hap_index.size() << " haplotype(s) and the model needs at least 2" << endl;
                linkage_weight = 0.0;
            } else {
                // Not an error, but worth saying: the default weight was tuned against panels of
                // tens of haplotypes and measures as roughly neutral on four, where a genotype is
                // spelled by too few haplotype pairs for the frequency prior to have anything to
                // act on. On a thin panel it costs runtime and buys close to nothing.
                if (hap_index.size() < 8) {
                    cerr << "warning [vg call]: linkage on a " << hap_index.size()
                         << "-haplotype panel; the default --linkage-weight 2 was tuned on 34 and"
                         << " measures as roughly neutral on 4 (consider --linkage-weight 0)"
                         << endl;
                }
                LinkageModel::Params linkage_params;
                linkage_params.weight = linkage_weight;
                linkage_params.scale = linkage_scale;
                linkage_params.freq_prior = linkage_freq_prior;
                linkage_collector.reset(new LinkageCollector(linkage_params, hap_index.size()));
                vcf_caller->set_linkage(linkage_collector.get(), gbwt_index,
                                        &linkage_sequence_to_haplotype);
                vcf_caller->set_emit_phasing(phased_output);
                // Panel index -> "sample#phase", so the mosaic names haplotypes rather than
                // only numbering them. Built by inverting hap_index, which is the same numbering
                // the model itself uses.
                vector<string> hap_names(hap_index.size());
                for (const auto& kv : hap_index) {
                    string sample = kv.first.first < meta.sample_names.size()
                                        ? meta.sample(kv.first.first)
                                        : string("sample") + std::to_string(kv.first.first);
                    hap_names[kv.second] = sample + "#" + std::to_string(kv.first.second);
                }
                // The reference paths this run called against, in full. The mosaic rows carry
                // only the locus part, and this graph has two reference samples (CHM13 and
                // GRCh38), so without these the coordinates name no particular assembly.
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
            // The read source, named in full so a consumer can tell which reads the offsets index
            // into -- the offsets are in the read AS SEQUENCED, so the file they came from is the
            // one an assembler must be given.
            string reads_source = !gaf_base_filename.empty()
                                      ? gaf_base_filename
                                      : (!gaf_filename.empty() ? gaf_filename : gam_filename);
            vcf_caller->set_anchors_out(anchors_out, anchor_params, graph_filename, reads_source,
                                        min_mismap_prob);
        }
        // one call covers FlowCaller (both ctors, so plain vg call gets it too), NestedFlowCaller
        // and LegacyCaller, since the merge lives on the shared VCFOutputCaller base
        vcf_caller->set_read_phasing(read_phasing, read_phasing_params);
        vcf_caller->set_regenotype(regenotype, regenotype_params, regenotype_passes,
                                   regenotype_ledger);
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

    // Ordered visits only help a read source that fetches by node-ID window, and they
    // change the traversal order of code the default caller shares -- so gate on such a
    // source actually being in use. With it off, the default path is bit-for-bit what it
    // was. GAF-Base needs this more than the GAM index does: a query there is a process
    // spawn, so an unordered visit pays milliseconds per site rather than a rescan.
    if (dynamic_cast<WindowedSiteReadSource*>(read_source.get()) != nullptr) {
        graph_caller->set_node_id_ordering(true, read_window_size);
        if (show_progress) {
            logger.info() << "Visiting snarls in node-ID order, window " << read_window_size << endl;
        }
    }

    // Descend into a nested chain once its parent's genotype is settled. Armed unconditionally
    // wherever nested calling runs, because it costs nothing where nothing defers: a snarl with no
    // linkage entry cannot be moved, so its per-site genotype is already final and its children are
    // visited inline. With no panel at all -- --enumerate-support, or a GBZ holding only reference
    // paths -- no site is ever recorded and every descent takes that inline path, which is the same
    // rule rather than a fallback to a different one.
    FlowCaller* deferring_caller = nullptr;
    // Armed wherever the linkage layer is, not only where nesting is. The layer can be armed with
    // nesting off -- `--no-nested`, and `--no-phased` too, which turns nesting off for its own
    // reasons above -- and in that configuration nothing was staged, the render pass was empty, and
    // patching a line after the fact was the only mechanism there was. Keeping that path alive for
    // one configuration would leave two ways to write a record, which is the thing this phase exists
    // to remove; arming here instead means a non-nested run resolves generation 0 over an empty
    // pending set and renders every record from the settled genotype, which is the same rule rather
    // than a second one.
    if (nested_calling || linkage_collector != nullptr) {
        deferring_caller = dynamic_cast<FlowCaller*>(graph_caller.get());
        if (deferring_caller != nullptr) {
            deferring_caller->set_defer_nested_descent(true);
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

    // Resolve, descend, repeat. Nothing below this needs to know which mode ran: the loop leaves
    // the linkage pass resolved, and write_variants' own resolve is idempotent.
    if (deferring_caller != nullptr) {
        // Barrier first, render second. The barrier settles every generation's genotypes, so the
        // render can build each record from the answer rather than from the reads' first guess -- and
        // a record built from the settled genotype needs no patch, which is what stage 11 removes.
        deferring_caller->run_deferred_descent();
        deferring_caller->render_retained_records();
    }

    // After both passes, so every settled record has had its chance to contribute. Written here
    // rather than from write_variants because it is not a VCF and shares none of its machinery.
    {
        auto* anchor_caller = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        if (anchor_caller != nullptr) {
            anchor_caller->write_anchors();
        }
    }

    // Report the indexed read source's cache behaviour, now that calling is done and
    // the counters mean something. The index over-fetches, so a low hit rate means the
    // parent-then-descendants locality the cache relies on is not materialising.
    if (show_progress) {
        auto* windowed = dynamic_cast<WindowedSiteReadSource*>(read_source.get());
        if (windowed != nullptr) {
            size_t hits = windowed->get_cache_hits();
            size_t misses = windowed->get_cache_misses();
            size_t total = hits + misses;
            auto* gaf_base = dynamic_cast<GafBaseSiteReadSource*>(windowed);
            logger.info() << (gaf_base != nullptr ? "GAF-Base: " : "Indexed GAM: ")
                          << windowed->get_read_count() << " reads fetched, "
                          << hits << "/" << total << " site queries served from cache"
                          << (total > 0 ? " (" + std::to_string((int)(100.0 * hits / total)) + "%)" : "")
                          << endl;
            // Candidates a site query looked at against reads it delivered. Inside a
            // window a candidate is one node-index entry, so a ratio near 1 means sites
            // name few nodes each read touches and a high one means long reads crossing
            // a site repeatedly; a straddling query examines whole reads instead.
            size_t seen = windowed->get_scanned_count();
            size_t used = windowed->get_delivered_count();
            logger.info() << (gaf_base != nullptr ? "GAF-Base: " : "Indexed GAM: ")
                          << seen << " read candidates examined, " << used << " delivered"
                          << (seen > 0 ? " (" + std::to_string((int)(100.0 * used / seen)) + "%)" : "")
                          << endl;
            // Sites too big for one window are fetched uncached, by their exact node
            // ranges. The two totals are the ranges asked for against the span they
            // sit in: if a change ever collapses one to the other, this is where it
            // shows, and on chr20 that difference was 133 k node IDs against 13.2 M.
            logger.info() << (gaf_base != nullptr ? "GAF-Base: " : "Indexed GAM: ")
                          << windowed->get_straddle_count() << " site queries too wide for a "
                          << "window, fetched uncached over " << windowed->get_straddle_wanted()
                          << " node IDs (spanning " << windowed->get_straddle_nodes() << ")"
                          << endl;
            if (gaf_base != nullptr) {
                // The count that governs run time: each one is a process spawn, so this
                // is the number to watch if a run is slow.
                logger.info() << "GAF-Base: " << gaf_base->get_query_count()
                              << " subprocess queries" << endl;
            }
        }
    }

    if (show_progress) logger.info() << "Calling complete" << endl;

    if (!gaf_output) {
        // Output VCF
        VCFOutputCaller* vcf_caller = dynamic_cast<VCFOutputCaller*>(graph_caller.get());
        assert(vcf_caller != nullptr);
        // Prune here rather than where the header string was built: the contig set is only
        // known once calling is done, and write_variants below drains the buffer it reads.
        cout << vcf_caller->prune_header_contigs(header, vcf_caller->get_output_contigs()) << flush;
        if (show_progress) logger.info() << "Writing VCF Variants" << endl;
        vcf_caller->write_variants(cout, snarl_manager.get());
        if (show_progress) logger.info() << "VCF complete" << endl;        
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_call("call", "call or genotype VCF variants", PIPELINE, 10, main_call);

