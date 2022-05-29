/**
 * \file giraffe_main.cpp: G(ir)AF (Graph Alignment Format) Fast Emitter: a fast short-read-to-haplotypes mapper
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <cassert>
#include <cstring>
#include <ctime>
#include <map>
#include <vector>
#include <unordered_set>
#include <chrono>
#include <mutex>

#include "subcommand.hpp"

#include "../snarl_seed_clusterer.hpp"
#include "../mapper.hpp"
#include "../annotation.hpp"
#include <vg/io/vpkg.hpp>
#include <vg/io/stream.hpp>
#include "../hts_alignment_emitter.hpp"
#include "../gapless_extender.hpp"
#include "../minimizer_mapper.hpp"
#include "../index_registry.hpp"
#include <bdsg/overlays/overlay_helper.hpp>

#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>

//#define USE_CALLGRIND

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#endif

#include <sys/ioctl.h>
#ifdef __linux__
#include <linux/perf_event.h>
#include <asm/unistd.h>
/// Bind perf_event_open for counting instructions.
/// See <https://stackoverflow.com/a/64863392>
static long perf_event_open(struct perf_event_attr* hw_event, pid_t pid, int cpu, int group_fd, unsigned long flags) {
    return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
}
#endif

using namespace std;
using namespace vg;
using namespace vg::subcommand;

namespace vg {

// Define a system for grid searching

/// This defines a range of values to test, from start to <=end, going up by step
template<typename Number>
struct Range {
    
    // Expose the thing we are a range of
    using type = Number;

    /// Represents the start of the range
    Number start = 0;
    /// Represents the inclusive end of the range
    Number end = 0;
    /// Represents the step to move by each tick
    Number step = 1;
    
    /// Represents the current value the range is at
    Number here = 0;
    /// Determines if we are running or not (i.e. is here valid)
    bool running = false;
    
    /// This will be called when we want to reset_chain what we are chained onto.
    function<void(void)> reset_chain_parent = []() {
    };
    /// This will be called when we need to tick_chain our parent
    function<bool(void)> tick_chain_parent = []() {
        return false;
    };
    
    /// Default constructor
    Range() {
        // Nothing to do!
    }
    
    /// Construct from a single value
    Range(const Number& val): start(val), end(val) {
        // Nothing to do!
    }
    
    /// Copy, preserving destination links
    Range(const Range& other): start(other.start), end(other.end), step(other.step) {
        // Nothing to do
    }
    
    /// Move, preserving destination links
    Range(Range&& other): start(other.start), end(other.end), step(other.step) {
        // Nothing to do
    }
    
    /// Copy assignment, preserving destination links
    Range& operator=(const Range& other) {
        start = other.start;
        end = other.end;
        step = other.step;
        return *this;
    }
    
    /// Move assignment, preserving destination links
    Range& operator=(Range&& other) {
        start = other.start;
        end = other.end;
        step = other.step;
        return *this;
    }
    
    /// Check the range for usefulness
    inline bool is_valid() {
        if (start != end && step == 0) {
            // We'll never make it
            cerr << "Invalid range (no movement): " << start << " to " << end << " step " << step << endl;
            return false;
        }
        
        if (start > end && step > 0) {
            // We're going the wrong way
            cerr << "Invalid range (need to go down): " << start << " to " << end << " step " << step << endl;
            return false;
        }
        
        if (start < end && step < 0) {
            // We're going the other wrong way
            cerr << "Invalid range (need to go up): " << start << " to " << end << " step " << step << endl;
            return false;
        }
        
        return true;
    }
    
    /// Convert to Number with the current value
    operator Number() const {
        if (running) {
            return here;
        } else {
            return start;
        }
    }
    
    /// Start at our start value
    void reset() {
        here = start;
        running = true;
    }
    
    /// Start us and all the things we are chained onto at their start values
    void reset_chain() {
        reset();
        reset_chain_parent();
    }
    
    /// Increment our value.
    /// Returns true if the new value needs processing, and false if we have left or would leave the range.
    bool tick() {
        if (here == end) {
            // We are at the end
            return false;
        }
        
        here += step;
        if ((step > 0 && here > end) || (step < 0 && here < end)) {
            // We have passed the end (for things like double)
            return false;
        }
        
        return true;
    }
    
    /// Increment our value.
    /// If it overflows, tock_chain whatever we are chained onto, and reset and succeed if that succeeds.
    bool tick_chain() {
        if (tick()) {
            // We could change
            return true;
        } else {
            // We couldn't change.
            if (tick_chain_parent()) {
                // We have a parent we could advance.
                reset();
                return true;
            } else {
                // Our parent couldn't advance either.
                return false;
            }
        }
    }
    
    /// Chain the given range onto this one.
    /// Return the passed-in range.
    /// Neither range may be moved away!
    template<typename Other>
    Range<Other>& chain(Range<Other>& next) {
        
        // Attach next to us
        next.reset_chain_parent = [&]() {
            this->reset_chain();
        };
        next.tick_chain_parent = [&]() {
            return this->tick_chain();
        };
        
        return next;
    }
    
    /// Get a function that runs another function for each combination of
    /// values for this Range and all Ranges it has been chained onto.
    function<void(const function<void(void)>&)> get_iterator() {
        return [&](const function<void(void)>& iteratee) {
            // Start
            reset_chain();
            
            do {
                // Run iteratee
                iteratee();
                // And tick the whole chain before running again
            } while(tick_chain());
        };
    }
};

// Define a way to test if a type is a Result<T>
// See https://stackoverflow.com/a/25803794

// In general, things aren't instantiations of things
template <typename Subject, template<typename...> class Predicate>
struct is_instantiation_of : std::false_type {
};

// Except things that are instantiations of things with some arguments
template <template<typename... > class Predicate, class... PredicateArgs>
struct is_instantiation_of<Predicate<PredicateArgs...>, Predicate> : std::true_type {
};

/// Parse a range as start[:end[:step]]
template<typename Result>
inline bool parse(const string& arg, typename enable_if<is_instantiation_of<Result, Range>::value, Result>::type& dest) {

    auto colon1 = arg.find(':');
    
    if (colon1 == string::npos) {
        // No colons here. Parse one number.
        if (!parse<typename Result::type>(arg, dest.start)) {
            return false;
        }
        dest.end = dest.start;
        dest.step = 0;
        return dest.is_valid();
    } else if (colon1 == arg.size()) {
        // Can't end in a colon
        return false;
    } else {
        // Look for another colon
        auto colon2 = arg.find(':', colon1 + 1);
        if (colon2 == string::npos) {
            // Just a range of two things
            if (!parse<typename Result::type>(arg.substr(0, colon1), dest.start)) {
                return false;
            }
            if (!parse<typename Result::type>(arg.substr(colon1 + 1), dest.end)) {
                return false;
            }
            dest.step = 1;
            return dest.is_valid();
        } else if (colon2 == arg.size()) {
            // Can't end in a colon
            return false;
        } else {
            // We have 3 numbers
            if (!parse<typename Result::type>(arg.substr(0, colon1), dest.start)) {
                return false;
            }
            if (!parse<typename Result::type>(arg.substr(colon1 + 1, colon2 - colon1 - 1), dest.end)) {
                return false;
            }
            if (!parse<typename Result::type>(arg.substr(colon2 + 1), dest.step)) {
                return false;
            }
            
            return dest.is_valid();
        }
    }
}
}

// Try stripping all suffixes in the vector, one at a time, and return on failure.
std::string strip_suffixes(std::string filename, const std::vector<std::string>& suffixes) {
    for (const std::string& suffix : suffixes) {
        if (filename.length() > suffix.length() && filename.substr(filename.length() - suffix.length()) == suffix) {
            filename = filename.substr(0, filename.length() - suffix.length());
        } else {
            break;
        }
    }
    return filename;
}

void help_giraffe(char** argv) {
    cerr
    << "usage: " << argv[0] << " giraffe [options] [ref.fa [variants.vcf.gz]] > output.gam" << endl
    << "Fast haplotype-aware short read mapper." << endl
    << endl
    << "basic options:" << endl
    << "  -Z, --gbz-name FILE           use this GBZ file (GBWT index + GBWTGraph)" << endl
    << "  -m, --minimizer-name FILE     use this minimizer index" << endl
    << "  -d, --dist-name FILE          cluster using this distance index" << endl
    << "  -p, --progress                show progress" << endl
    << "input options:" << endl
    << "  -G, --gam-in FILE             read and realign GAM-format reads from FILE" << endl
    << "  -f, --fastq-in FILE           read and align FASTQ-format reads from FILE (two are allowed, one for each mate)" << endl
    << "  -i, --interleaved             GAM/FASTQ input is interleaved pairs, for paired-end alignment" << endl
    << "alternate indexes:" << endl
    << "  -x, --xg-name FILE            use this xg index or graph" << endl
    << "  -g, --graph-name FILE         use this GBWTGraph" << endl
    << "  -H, --gbwt-name FILE          use this GBWT index" << endl
    << "output options:" << endl
    << "  -M, --max-multimaps INT       produce up to INT alignments for each read [1]" << endl
    << "  -N, --sample NAME             add this sample name" << endl
    << "  -R, --read-group NAME         add this read group" << endl
    << "  -o, --output-format NAME      output the alignments in NAME format (gam / gaf / json / tsv / SAM / BAM / CRAM) [gam]" << endl
    << "  --ref-paths FILE              ordered list of paths in the graph, one per line or HTSlib .dict, for HTSLib @SQ headers" << endl
    << "  --named-coordinates           produce GAM outputs in named-segment (GFA) space" << endl
    << "  -P, --prune-low-cplx          prune short and low complexity anchors during linear format realignment" << endl
    << "  -n, --discard                 discard all output alignments (for profiling)" << endl
    << "  --output-basename NAME        write output to a GAM file beginning with the given prefix for each setting combination" << endl
    << "  --report-name NAME            write a TSV of output file and mapping speed to the given file" << endl
    << "  --show-work                   log how the mapper comes to its conclusions about mapping locations" << endl
    << "algorithm presets:" << endl
    << "  -b, --parameter-preset NAME   set computational parameters (fast / default) [default]" << endl
    << "computational parameters:" << endl
    << "  -c, --hit-cap INT             use all minimizers with at most INT hits [10]" << endl
    << "  -C, --hard-hit-cap INT        ignore all minimizers with more than INT hits [500]" << endl
    << "  -F, --score-fraction FLOAT    select minimizers between hit caps until score is FLOAT of total [0.9]" << endl
    << "  -D, --distance-limit INT      cluster using this distance limit [200]" << endl
    << "  -e, --max-extensions INT      extend up to INT clusters [800]" << endl
    << "  -a, --max-alignments INT      align up to INT extensions [8]" << endl
    << "  -s, --cluster-score INT       only extend clusters if they are within INT of the best score [50]" << endl
    << "  -S, --pad-cluster-score INT   also extend clusters within INT of above threshold to get a second-best cluster [0]" << endl
    << "  -u, --cluster-coverage FLOAT  only extend clusters if they are within FLOAT of the best read coverage [0.3]" << endl
    << "  -v, --extension-score INT     only align extensions if their score is within INT of the best score [1]" << endl
    << "  -w, --extension-set INT       only align extension sets if their score is within INT of the best score [20]" << endl
    << "  -O, --no-dp                   disable all gapped alignment" << endl
    << "  -r, --rescue-attempts         attempt up to INT rescues per read in a pair [15]" << endl
    << "  -A, --rescue-algorithm NAME   use algorithm NAME for rescue (none / dozeu / gssw / haplotypes) [dozeu]" << endl
    << "  -L, --max-fragment-length INT assume that fragment lengths should be smaller than INT when estimating the fragment length distribution" << endl
    << "  --fragment-mean FLOAT         force the fragment length distribution to have this mean (requires --fragment-stdev)" << endl
    << "  --fragment-stdev FLOAT        force the fragment length distribution to have this standard deviation (requires --fragment-mean)" << endl
    << "  --paired-distance-limit FLOAT cluster pairs of read using a distance limit FLOAT standard deviations greater than the mean [2.0]" << endl
    << "  --rescue-subgraph-size FLOAT  search for rescued alignments FLOAT standard deviations greater than the mean [4.0]" << endl
    << "  --rescue-seed-limit INT       attempt rescue with at most INT seeds [100]" << endl
    << "  --track-provenance            track how internal intermediate alignment candidates were arrived at" << endl
    << "  --track-correctness           track if internal intermediate alignment candidates are correct (implies --track-provenance)" << endl
    << "  -t, --threads INT             number of mapping threads to use" << endl;
}

int main_giraffe(int argc, char** argv) {

    std::chrono::time_point<std::chrono::system_clock> launch = std::chrono::system_clock::now();

    if (argc == 2) {
        help_giraffe(argv);
        return 1;
    }

    #define OPT_OUTPUT_BASENAME 1001
    #define OPT_REPORT_NAME 1002
    #define OPT_TRACK_PROVENANCE 1003
    #define OPT_TRACK_CORRECTNESS 1004
    #define OPT_FRAGMENT_MEAN 1005
    #define OPT_FRAGMENT_STDEV 1006
    #define OPT_CLUSTER_STDEV 1007
    #define OPT_RESCUE_STDEV 1008
    #define OPT_RESCUE_SEED_LIMIT 1009
    #define OPT_REF_PATHS 1010
    #define OPT_SHOW_WORK 1011
    #define OPT_NAMED_COORDINATES 1012
    

    // initialize parameters with their default options
    
    // This holds and manages finding our indexes.
    IndexRegistry registry = VGIndexes::get_vg_index_registry();
    string output_basename;
    string report_name;
    // How close should two hits be to be in the same cluster?
    Range<size_t> distance_limit = 200;
    Range<size_t> hit_cap = 10, hard_hit_cap = 500;
    Range<double> minimizer_score_fraction = 0.9;
    bool show_progress = false;
    // Should we try chaining or just give up if we can't find a full length gapless alignment?
    bool do_dp = true;
    // What GAM should we realign?
    string gam_filename;
    // What FASTQs should we align.
    // Note: multiple FASTQs are not interpreted as paired.
    string fastq_filename_1;
    string fastq_filename_2;
    // Is the input interleaved/are we in paired-end mode?
    bool interleaved = false;
    // True if fastq_filename_2 or interleaved is set.
    bool paired = false;
    string param_preset = "default";
    // How many mappings per read can we emit?
    Range<size_t> max_multimaps = 1;
    // How many clusters should we extend?
    Range<size_t> max_extensions = 800;
    // How many extended clusters should we align, max?
    Range<size_t> max_alignments = 8;
    //Throw away cluster with scores that are this amount below the best
    Range<double> cluster_score = 50;
    //Unless they are the second best and within this amount beyond that
    Range<double> pad_cluster_score = 0;
    //Throw away clusters with coverage this amount below the best 
    Range<double> cluster_coverage = 0.3;
    //Throw away extension sets with scores that are this amount below the best
    Range<double> extension_set = 20;
    //Throw away extensions with scores that are this amount below the best
    Range<int> extension_score = 1;
    //Attempt up to this many rescues of reads with no pairs
    bool forced_rescue_attempts = false;
    Range<int> rescue_attempts = 15;
    //Don't let distances larger than this contribute to the fragment length distribution
    size_t fragment_length = 2000;
    // Which rescue algorithm do we use?
    MinimizerMapper::RescueAlgorithm rescue_algorithm = MinimizerMapper::rescue_dozeu;
    //Did we force the fragment length distribution?
    bool forced_mean = false;
    //And if so what is it?
    double fragment_mean = 0.0;
    bool forced_stdev = false;
    double fragment_stdev = 0.0;
    //How many sdevs to we look out when clustering pairs?
    double cluster_stdev = 2.0;
    //How many stdevs do we look out when rescuing? 
    double rescue_stdev = 4.0;
    // Attempt rescue with up to this many seeds.
    size_t rescue_seed_limit = 100;
    // How many pairs should we be willing to buffer before giving up on fragment length estimation?
    size_t MAX_BUFFERED_PAIRS = 100000;
    // What sample name if any should we apply?
    string sample_name;
    // What read group if any should we apply?
    string read_group;
    // Should we throw out our alignments instead of outputting them?
    bool discard_alignments = false;
    // Should we track candidate provenance?
    bool track_provenance = false;
    // Should we track candidate correctness?
    bool track_correctness = false;
    // Should we log our mapping decision making?
    bool show_work = false;

    // Chain all the ranges and get a function that loops over all combinations.
    auto for_each_combo = distance_limit
        .chain(hit_cap)
        .chain(hard_hit_cap)
        .chain(minimizer_score_fraction)
        .chain(rescue_attempts)
        .chain(max_multimaps)
        .chain(max_extensions)
        .chain(max_alignments)
        .chain(cluster_score)
        .chain(pad_cluster_score)
        .chain(cluster_coverage)
        .chain(extension_set)
        .chain(extension_score)
        .get_iterator();
    

    // Formats for alignment output.
    std::string output_format = "GAM";
    std::set<std::string> output_formats = { "GAM", "GAF", "JSON", "TSV", "SAM", "BAM", "CRAM" };

    // For HTSlib formats, where do we get sequence header info?
    std::string ref_paths_name;
    // And should we drop low complexity anchors when surjectng?
    bool prune_anchors = false;
    
    // For GAM format, should we report in named-segment space instead of node ID space?
    bool named_coordinates = false;

    // Map algorithm names to rescue algorithms
    std::map<std::string, MinimizerMapper::RescueAlgorithm> rescue_algorithms = {
        { "none", MinimizerMapper::rescue_none },
        { "dozeu", MinimizerMapper::rescue_dozeu },
        { "gssw", MinimizerMapper::rescue_gssw },
        { "haplotypes", MinimizerMapper::rescue_haplotypes }
    };
    std::map<MinimizerMapper::RescueAlgorithm, std::string> algorithm_names =  {
        { MinimizerMapper::rescue_none, "none" },
        { MinimizerMapper::rescue_dozeu, "dozeu" },
        { MinimizerMapper::rescue_gssw, "gssw" },
        { MinimizerMapper::rescue_haplotypes, "haplotypes" }
    };
    //TODO: Right now there can be two versions of the distance index. This ensures that the correct minimizer type gets built

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"gbz-name", required_argument, 0, 'Z'},
            {"xg-name", required_argument, 0, 'x'},
            {"graph-name", required_argument, 0, 'g'},
            {"gbwt-name", required_argument, 0, 'H'},
            {"minimizer-name", required_argument, 0, 'm'},
            {"dist-name", required_argument, 0, 'd'},
            {"progress", no_argument, 0, 'p'},
            {"gam-in", required_argument, 0, 'G'},
            {"fastq-in", required_argument, 0, 'f'},
            {"interleaved", no_argument, 0, 'i'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"output-format", required_argument, 0, 'o'},
            {"ref-paths", required_argument, 0, OPT_REF_PATHS},
            {"prune-low-cplx", no_argument, 0, 'P'},
            {"named-coordinates", no_argument, 0, OPT_NAMED_COORDINATES},
            {"discard", no_argument, 0, 'n'},
            {"output-basename", required_argument, 0, OPT_OUTPUT_BASENAME},
            {"report-name", required_argument, 0, OPT_REPORT_NAME},
            {"fast-mode", no_argument, 0, 'b'},
            {"hit-cap", required_argument, 0, 'c'},
            {"hard-hit-cap", required_argument, 0, 'C'},
            {"distance-limit", required_argument, 0, 'D'},
            {"max-extensions", required_argument, 0, 'e'},
            {"max-alignments", required_argument, 0, 'a'},
            {"cluster-score", required_argument, 0, 's'},
            {"pad-cluster-score", required_argument, 0, 'S'},
            {"cluster-coverage", required_argument, 0, 'u'},
            {"extension-score", required_argument, 0, 'v'},
            {"extension-set", required_argument, 0, 'w'},
            {"score-fraction", required_argument, 0, 'F'},
            {"no-dp", no_argument, 0, 'O'},
            {"rescue-attempts", required_argument, 0, 'r'},
            {"rescue-algorithm", required_argument, 0, 'A'},
            {"paired-distance-limit", required_argument, 0, OPT_CLUSTER_STDEV },
            {"rescue-subgraph-size", required_argument, 0, OPT_RESCUE_STDEV },
            {"rescue-seed-limit", required_argument, 0, OPT_RESCUE_SEED_LIMIT},
            {"max-fragment-length", required_argument, 0, 'L' },
            {"fragment-mean", required_argument, 0, OPT_FRAGMENT_MEAN },
            {"fragment-stdev", required_argument, 0, OPT_FRAGMENT_STDEV },
            {"track-provenance", no_argument, 0, OPT_TRACK_PROVENANCE},
            {"track-correctness", no_argument, 0, OPT_TRACK_CORRECTNESS},
            {"show-work", no_argument, 0, OPT_SHOW_WORK},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hZ:x:g:H:m:s:d:pG:f:iM:N:R:o:Pnb:c:C:D:F:e:a:S:u:v:w:Ot:r:A:L:",
                         long_options, &option_index);


        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
            case 'Z':
                if (!optarg || !*optarg) {
                    cerr << "error:[vg giraffe] Must provide GBZ file with -Z." << endl;
                    exit(1);
                }
                if (!std::ifstream(optarg).is_open()) {
                    cerr << "error:[vg giraffe] Couldn't open GBZ file " << optarg << endl;
                    exit(1);
                }
                registry.provide("Giraffe GBZ", optarg);

                // If we have a GBZ we probably want to use its name as the base name.
                // But see -g.
                registry.set_prefix(strip_suffixes(std::string(optarg), { ".gbz", ".giraffe" }));

                break;

            case 'x':
                if (!optarg || !*optarg) {
                    cerr << "error:[vg giraffe] Must provide graph file with -x." << endl;
                    exit(1);
                }
                if (!std::ifstream(optarg).is_open()) {
                    cerr << "error:[vg giraffe] Couldn't open graph file " << optarg << endl;
                    exit(1); 
                }
                registry.provide("XG", optarg);
                
                // If we have an xg we probably want to use its name as the base name.
                // But see -g.
                registry.set_prefix(split_ext(optarg).first);
                
                break;

            case 'g':
                if (!optarg || !*optarg) {
                    cerr << "error:[vg giraffe] Must provide GBWTGraph file with -g." << endl;
                    exit(1);
                }
                if (!std::ifstream(optarg).is_open()) {
                    cerr << "error:[vg giraffe] Couldn't open GBWTGraph file " << optarg << endl;
                    exit(1); 
                }
                registry.provide("GBWTGraph", optarg);
                
                // But if we have a GBWTGraph we probably want to use *its* name as the base name.
                // Whichever is specified last will win, unless we also have a FASTA input name.
                registry.set_prefix(split_ext(optarg).first);
                
                break;

            case 'H':
                if (!optarg || !*optarg) {
                    cerr << "error:[vg giraffe] Must provide GBWT file with -H." << endl;
                    exit(1);
                }
                if (!std::ifstream(optarg).is_open()) {
                    cerr << "error:[vg giraffe] Couldn't open GBWT file " << optarg << endl;
                    exit(1); 
                }
                registry.provide("Giraffe GBWT", optarg);
                break;
                
            case 'm':
                if (!optarg || !*optarg) {
                    cerr << "error:[vg giraffe] Must provide minimizer file with -m." << endl;
                    exit(1);
                }
                if (!std::ifstream(optarg).is_open()) {
                    cerr << "error:[vg giraffe] Couldn't open minimizer file " << optarg << endl;
                    exit(1); 
                }
                registry.provide("Minimizers", optarg);
                break;
                
            case 'd':
                if (!optarg || !*optarg) {
                    cerr << "error:[vg giraffe] Must provide distance index file with -d." << endl;
                    exit(1);
                }
                if (!std::ifstream(optarg).is_open()) {
                    cerr << "error:[vg giraffe] Couldn't open distance index file " << optarg << endl;
                    exit(1); 
                }
                registry.provide("Giraffe Distance Index", optarg);
                break;

            case 'p':
                show_progress = true;
                break;
                
            case 'G':
                gam_filename = optarg;
                if (gam_filename.empty()) {
                    cerr << "error:[vg giraffe] Must provide GAM file with -G." << endl;
                    exit(1);
                }
                break;
            
            case 'f':
                if (fastq_filename_1.empty()) {
                    fastq_filename_1 = optarg;
                    if (fastq_filename_1.empty()) {
                        cerr << "error:[vg giraffe] Must provide FASTQ file with -f." << endl;
                        exit(1);
                    }
                }
                else if (fastq_filename_2.empty()) {
                    fastq_filename_2 = optarg;
                    if (fastq_filename_2.empty()) {
                        cerr << "error:[vg giraffe] Must provide FASTQ file with -f." << endl;
                        exit(1);
                    }
                    paired = true;
                } else {
                    cerr << "error:[vg giraffe] Cannot specify more than two FASTQ files." << endl;
                    exit(1);
                }
                break;

            case 'i':
                interleaved = true;
                paired = true;
                break;
                
            case 'M':
                max_multimaps = parse<Range<size_t>>(optarg);
                break;
            
            case 'N':
                sample_name = optarg;
                break;
                
            case 'R':
                read_group = optarg;
                break;

            case 'o':
                {
                    output_format = optarg;
                    for (char& c : output_format) {
                        c = std::toupper(c);
                    }
                    if (output_formats.find(output_format) == output_formats.end()) {
                        std::cerr << "error: [vg giraffe] Invalid output format: " << optarg << std::endl;
                        std::exit(1);
                    }
                }
                break;
                
            case OPT_REF_PATHS:
                ref_paths_name = optarg;
                break;
                
            case 'P':
                prune_anchors = true;
                break;
                
            case OPT_NAMED_COORDINATES:
                named_coordinates = true;
                break;

            case 'n':
                discard_alignments = true;
                break;
                
            case OPT_OUTPUT_BASENAME:
                output_basename = optarg;
                break;
            
            case OPT_REPORT_NAME:
                report_name = optarg;
                break;
            case 'b':
                param_preset = optarg;

                if (param_preset == "fast" ) {

                    hit_cap = 10;
                    hard_hit_cap = 500;
                    minimizer_score_fraction = 0.5;
                    max_multimaps = 1;
                    max_extensions = 400;
                    max_alignments = 8;
                    cluster_score = 50;
                    pad_cluster_score = 0;
                    cluster_coverage = 0.2;
                    extension_set = 20;
                    extension_score = 1;
                } else if (param_preset != "default" ) {
                    std::cerr << "error: [vg giraffe] invalid parameter preset: " << optarg << std:: endl;
                }
                break;

            case 'c':
                {
                    auto cap = parse<Range<size_t>>(optarg);
                    if (cap <= 0) {
                        cerr << "error: [vg giraffe] Hit cap (" << cap << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    hit_cap = cap;
                }
                break;

            case 'C':
                {
                    auto cap = parse<Range<size_t>>(optarg);
                    if (cap <= 0) {
                        cerr << "error: [vg giraffe] Hard hit cap (" << cap << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    hard_hit_cap = cap;
                }
                break;
                
            case 'D':
                {
                    auto limit = parse<Range<size_t>>(optarg);
                    if (limit <= 0) {
                        cerr << "error: [vg giraffe] Distance limit (" << limit << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    distance_limit = limit;
                }
                break;

            case 'F':
                minimizer_score_fraction = parse<Range<double>>(optarg);
                break;

            case 'e':
                {
                    auto extensions = parse<Range<size_t>>(optarg);
                    if (extensions <= 0) {
                        cerr << "error: [vg giraffe] Number of extensions (" << extensions << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    max_extensions = extensions;
                }
                break;

            case 'a':
                {
                    auto alignments = parse<Range<size_t>>(optarg);
                    if (alignments <= 0) {
                        cerr << "error: [vg giraffe] Number of alignments (" << alignments << ") must be a positive integer" << endl;
                        exit(1);
                    }
                    max_alignments = alignments;
                }
                break;

            case 's':
                {
                    auto score = parse<Range<double>>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg giraffe] Cluster score threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    cluster_score = score;
                }
                break;
                
            case 'S':
                {
                    auto score = parse<Range<double>>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg giraffe] Second best cluster score threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    pad_cluster_score = score;
                }
                break;

            case 'u':
                {
                    auto score = parse<Range<double>>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg giraffe] Cluster coverage threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    cluster_coverage = score;
                }
                break;
            case 'v':
                {
                    auto score = parse<Range<int>>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg giraffe] Extension score threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    extension_score = score;
                }
                break;
            case 'w':
                {
                    auto score = parse<Range<double>>(optarg);
                    if (score < 0) {
                        cerr << "error: [vg giraffe] Extension set score threshold (" << score << ") must be positive" << endl;
                        exit(1);
                    }
                    extension_set = score;
                }
                break;
                
            case 'O':
                do_dp = false;
                break;
                
            case 'r':
                {
                    forced_rescue_attempts = true;
                    rescue_attempts = parse<Range<int>>( optarg);
                    if (rescue_attempts < 0) {
                        cerr << "error: [vg giraffe] Rescue attempts must be positive" << endl;
                        exit(1);
                    }
                }
                break;

            case 'A':
                {
                    std::string algo_name = optarg;
                    for (char& c : algo_name) {
                        c = std::tolower(c);
                    }
                    auto iter = rescue_algorithms.find(algo_name);
                    if (iter == rescue_algorithms.end()) {
                        std::cerr << "error: [vg giraffe] Invalid rescue algorithm: " << optarg << std::endl;
                        std::exit(1);
                    }
                    rescue_algorithm = iter->second;
                }
                break;

            case OPT_FRAGMENT_MEAN:
                forced_mean = true;
                fragment_mean = parse<double>(optarg);
                break;

            case OPT_FRAGMENT_STDEV:
                forced_stdev = true;
                fragment_stdev = parse<double>(optarg);
                break;
            case 'L':
                fragment_length = parse<size_t>(optarg);
                break;

            case OPT_CLUSTER_STDEV:
                cluster_stdev = parse<double>(optarg);
                break;

            case OPT_RESCUE_STDEV:
                rescue_stdev = parse<double>(optarg);
                break;

            case OPT_RESCUE_SEED_LIMIT:
                rescue_seed_limit = parse<size_t>(optarg);
                break;

            case OPT_TRACK_PROVENANCE:
                track_provenance = true;
                break;
            
            case OPT_TRACK_CORRECTNESS:
                track_provenance = true;
                track_correctness = true;
                break;
                
            case OPT_SHOW_WORK:
                show_work = true;
                break;
                
            case 't':
            {
                int num_threads = parse<int>(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[vg giraffe] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
            }
                break;
                
            case 'h':
            case '?':
            default:
                help_giraffe(argv);
                exit(1);
                break;
        }
    }

   
    // Get positional arguments before validating user intent
    if (have_input_file(optind, argc, argv)) {
        // Must be the FASTA, but check.
        
        string fasta_filename = get_input_file_name(optind, argc, argv);
        
        auto fasta_parts = split_ext(fasta_filename);
        if (fasta_parts.second == "gz") {
            fasta_parts = split_ext(fasta_parts.first);
        }
        if (fasta_parts.second != "fa" && fasta_parts.second != "fasta" && fasta_parts.second != "fna") {
            cerr << "error:[vg giraffe] FASTA file " << fasta_filename << " is not named like a FASTA" << endl;
            exit(1);
        }
        
        registry.provide("Reference FASTA", fasta_filename);
        // Everything else should be named like the FASTA by default
        registry.set_prefix(fasta_parts.first);
        
        if (have_input_file(optind, argc, argv)) {
            // Next one must be VCF, but check.
            // TODO: Unify with FASTA check?
            
            string vcf_filename = get_input_file_name(optind, argc, argv);
            
            auto vcf_parts = split_ext(vcf_filename);
            if (vcf_parts.second == "gz") {
                vcf_parts = split_ext(vcf_parts.first);
            }
            if (vcf_parts.second != "vcf") {
                cerr << "error:[vg giraffe] VCF file " << vcf_filename << " is not named like a VCF" << endl;
                exit(1);
            }
            
            // Determine if it is phased or not
            string file_type = IndexRegistry::vcf_is_phased(vcf_filename) ? "VCF w/ Phasing" : "VCF";
            
            // Feed it to the index registry to maybe use
            registry.provide(file_type, vcf_filename);
        }
    }

    // If we don't want rescue, let the user see we don't try it.
    if (rescue_attempts == 0 || rescue_algorithm == MinimizerMapper::rescue_none) {
        rescue_attempts = 0;
        rescue_algorithm = MinimizerMapper::rescue_none;
    }
    
    // Now all the arguments are parsed, so see if they make sense
    
    // Decide if we are outputting to an htslib format
    bool hts_output = (output_format == "SAM" || output_format == "BAM" || output_format == "CRAM");
    
    if (!ref_paths_name.empty() && !hts_output) {
        cerr << "warning:[vg giraffe] Reference path file (--ref-paths) is only used when output format (-o) is SAM, BAM, or CRAM." << endl;
        ref_paths_name = "";
    }
    
    if (output_format != "GAM" && !output_basename.empty()) {
        cerr << "error:[vg giraffe] Using an output basename (--output-basename) only makes sense for GAM format (-o)" << endl;
        exit(1);
    }
    
    if (interleaved && !fastq_filename_2.empty()) {
        cerr << "error:[vg giraffe] Cannot designate both interleaved paired ends (-i) and separate paired end file (-f)." << endl;
        exit(1);
    }

    if (!fastq_filename_1.empty() && !gam_filename.empty()) {
        cerr << "error:[vg giraffe] Cannot designate both FASTQ input (-f) and GAM input (-G) in same run." << endl;
        exit(1);
    }
    
    if (have_input_file(optind, argc, argv)) {
        // TODO: work out how to interpret additional files as reads.
        cerr << "error:[vg giraffe] Extraneous input file: " << get_input_file_name(optind, argc, argv) << endl;
        exit(1);
    }

    if ((forced_mean && ! forced_stdev) || (!forced_mean && forced_stdev)) {
        cerr << "warning:[vg giraffe] Both a mean and standard deviation must be specified for the fragment length distribution" << endl;
        cerr << "                   Detecting fragment length distribution automatically" << endl;
        forced_mean = false;
        forced_stdev = false;
        fragment_mean = 0.0;
        fragment_stdev = 0.0;
    }
    if ((forced_mean || forced_stdev || forced_rescue_attempts) && (!paired)) {
        cerr << "warning:[vg giraffe] Attempting to set paired-end parameters but running in single-end mode" << endl;
    }
    
    // The IndexRegistry doesn't try to infer index files based on the
    // basename, so do that here. We can have multiple extension options that
    // we try in order of priority.
    unordered_map<string, vector<string>> indexes_and_extensions = {
        {"Giraffe GBZ", {"giraffe.gbz", "gbz"}},
        {"XG", {"xg"}},
        {"Giraffe GBWT", {"gbwt"}},
        {"GBWTGraph", {"gg"}},
        {"Giraffe Distance Index", {"dist"}},
        {"Minimizers", {"min"}}
    };
    for (auto& completed : registry.completed_indexes()) {
        // Drop anything we already got from the list
        indexes_and_extensions.erase(completed);
    }
    for (auto& index_and_extensions : indexes_and_extensions) {
        // For each index type
        for (auto& extension : index_and_extensions.second) {
            // For each extension in priority order
            string inferred_filename = registry.get_prefix() + "." + extension;
            if (ifstream(inferred_filename).is_open()) {
                // A file with the appropriate name exists and we can read it
                registry.provide(index_and_extensions.first, inferred_filename);
                // Report it because this may not be desired behavior
                cerr << "Guessing that " << inferred_filename << " is " << index_and_extensions.first << endl;
                // Skip other extension options for the index
                break;
            }
        }
    }

    // create in-memory objects
    
    // Don't try and use all the memory.
    // TODO: add memory options like autoindex?
    registry.set_target_memory_usage(IndexRegistry::get_system_memory() / 2);
    
    auto index_targets = VGIndexes::get_default_giraffe_indexes();

#ifdef debug
    for (auto& needed : index_targets) {
        cerr << "Want index: " << needed << endl;
    }
#endif
    
    try {
        registry.make_indexes(index_targets);
    }
    catch (InsufficientInputException ex) {
        cerr << "error:[vg giraffe] Input is not sufficient to create indexes" << endl;
        cerr << ex.what();
        return 1;
    }
   
#ifdef debug
    for (auto& completed : registry.completed_indexes()) {
        cerr << "Have index: " << completed << endl;
        for (auto& filename : registry.require(completed)) {
            cerr << "\tAt: " << filename << endl;
        }
    }
#endif
    
    // Grab the minimizer index
    auto minimizer_index = vg::io::VPKG::load_one<gbwtgraph::DefaultMinimizerIndex>(registry.require("Minimizers").at(0));

    // Grab the GBZ
    auto gbz = vg::io::VPKG::load_one<gbwtgraph::GBZ>(registry.require("Giraffe GBZ").at(0));

    // Grab the distance index
    std::unique_ptr<MinimumDistanceIndex> old_distance_index = nullptr;
    SnarlDistanceIndex distance_index;
    SnarlDistanceIndex* distance_index_ptr = nullptr;
    string distance_index_file_name = registry.require("Giraffe Distance Index").at(0);
    ifstream infile_dist (distance_index_file_name);
    if (vg::io::MessageIterator::sniff_tag(infile_dist) == "DistanceIndex") {
        old_distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(registry.require("Giraffe Distance Index").at(0));
    } else {
        distance_index.deserialize(distance_index_file_name);
        distance_index_ptr = &distance_index;
    }
    
    // If we are tracking correctness, we will fill this in with a graph for
    // getting offsets along ref paths.
    PathPositionHandleGraph* path_position_graph = nullptr;
    // If we need an overlay for position lookup, we might be pointing into
    // this overlay
    bdsg::PathPositionOverlayHelper overlay_helper;
    // And we might load an XG
    unique_ptr<PathHandleGraph> xg_graph;
    if (track_correctness || hts_output) {
        // Usually we will get our paths from the GBZ
        PathHandleGraph* base_graph = &gbz->graph;
        // But if an XG is around, we should use that instead. Otherwise, it's not possible to provide paths when using an old GBWT/GBZ that doesn't have them.
        if (registry.available("XG")) {
            xg_graph = vg::io::VPKG::load_one<PathHandleGraph>(registry.require("XG").at(0));
            base_graph = xg_graph.get();
        }
    
        // Apply the overlay if needed.
        path_position_graph = overlay_helper.apply(base_graph);
    }

    // Set up the mapper
    if (show_progress) {
        cerr << "Initializing MinimizerMapper" << endl;
    }
    MinimizerMapper minimizer_mapper(gbz->graph, *minimizer_index, &(*old_distance_index), distance_index_ptr, path_position_graph);
    if (forced_mean && forced_stdev) {
        minimizer_mapper.force_fragment_length_distr(fragment_mean, fragment_stdev);
    }

    
    std::chrono::time_point<std::chrono::system_clock> init = std::chrono::system_clock::now();
    std::chrono::duration<double> init_seconds = init - launch;
    if (show_progress) {
        cerr << "Loading and initialization: " << init_seconds.count() << " seconds" << endl;
    }
    
    // Set up to write a report of mapping speed if requested, instead of just dumping to stderr.
    ofstream report;
    if (!report_name.empty()) {
        // Open the report
        report.open(report_name);
        if (!report) {
            // Make sure it worked
            cerr << "error[vg giraffe]: Could not open report file " << report_name << endl;
            exit(1);
        }
        
        // Add a header
        report << "#file\treads/second/thread" << endl;
    }

    // We need to loop over all the ranges...
    for_each_combo([&]() {
    
        // Work out where to send the output. Default to stdout.
        string output_filename = "-";
        if (!output_basename.empty()) {
            // Compose a name using all the parameters.
            stringstream s;
            
            s << output_basename;
            
            if (interleaved) {
                s << "-i";
            }
            s << "-D" << distance_limit;
            s << "-c" << hit_cap;
            s << "-C" << hard_hit_cap;
            s << "-F" << minimizer_score_fraction;
            s << "-M" << max_multimaps;
            s << "-e" << max_extensions;
            s << "-a" << max_alignments;
            s << "-s" << cluster_score;
            s << "-u" << cluster_coverage;
            s << "-w" << extension_set;
            s << "-v" << extension_score;
            
            s << ".gam";
            
            output_filename = s.str();
        }
    
        if (show_progress) {
            if (discard_alignments) {
                cerr << "Discarding output alignments" << endl;
            } else {
                cerr << "Mapping reads to \"" << output_filename << "\" (" << output_format << ")" << endl;
            }
        }

        if (show_progress && interleaved) {
            cerr << "--interleaved" << endl;
        }
        
        if (show_progress && prune_anchors) {
            cerr << "--prune-low-cplx" << endl;
        }

        if (show_progress) {
            cerr << "--hit-cap " << hit_cap << endl;
        }
        minimizer_mapper.hit_cap = hit_cap;

        if (show_progress) {
            cerr << "--hard-hit-cap " << hard_hit_cap << endl;
        }
        minimizer_mapper.hard_hit_cap = hard_hit_cap;

        if (show_progress) {
            cerr << "--score-fraction " << minimizer_score_fraction << endl;
        }
        minimizer_mapper.minimizer_score_fraction = minimizer_score_fraction;

        if (show_progress) {
            cerr << "--max-extensions " << max_extensions << endl;
        }
        minimizer_mapper.max_extensions = max_extensions;

        if (show_progress) {
            cerr << "--max-alignments " << max_alignments << endl;
        }
        minimizer_mapper.max_alignments = max_alignments;

        if (show_progress) {
            cerr << "--cluster-score " << cluster_score << endl;
        }
        minimizer_mapper.cluster_score_threshold = cluster_score;
        
        if (show_progress) {
            cerr << "--pad-cluster-score " << pad_cluster_score << endl;
        }
        minimizer_mapper.pad_cluster_score_threshold = pad_cluster_score;

        if (show_progress) {
            cerr << "--cluster-coverage " << cluster_coverage << endl;
        }
        minimizer_mapper.cluster_coverage_threshold = cluster_coverage;

        if (show_progress) {
            cerr << "--extension-score " << extension_score << endl;
        }
        minimizer_mapper.extension_score_threshold = extension_score;

        if (show_progress) {
            cerr << "--extension-set " << extension_set << endl;
        }
        minimizer_mapper.extension_set_score_threshold = extension_set;

        if (show_progress && !do_dp) {
            cerr << "--no-dp " << endl;
        }
        minimizer_mapper.do_dp = do_dp;

        if (show_progress) {
            cerr << "--max-multimaps " << max_multimaps << endl;
        }
        minimizer_mapper.max_multimaps = max_multimaps;

        if (show_progress) {
            cerr << "--distance-limit " << distance_limit << endl;
        }
        minimizer_mapper.distance_limit = distance_limit;
        
        if (show_progress && track_provenance) {
            cerr << "--track-provenance " << endl;
        }
        minimizer_mapper.track_provenance = track_provenance;
        
        if (show_progress && track_correctness) {
            cerr << "--track-correctness " << endl;
        }
        minimizer_mapper.track_correctness = track_correctness;
        
        if (show_progress && show_work) {
            cerr << "--show-work " << endl;
        }
        minimizer_mapper.show_work = show_work;

        if (show_progress && paired) {
            if (forced_mean && forced_stdev) {
                cerr << "--fragment-mean " << fragment_mean << endl; 
                cerr << "--fragment-stdev " << fragment_stdev << endl;
            }
            cerr << "--paired-distance-limit " << cluster_stdev << endl;
            cerr << "--rescue-subgraph-size " << rescue_stdev << endl;
            cerr << "--rescue-seed-limit " << rescue_seed_limit << endl;
            cerr << "--rescue-attempts " << rescue_attempts << endl;
            cerr << "--rescue-algorithm " << algorithm_names[rescue_algorithm] << endl;
        }
        minimizer_mapper.max_fragment_length = fragment_length;
        minimizer_mapper.paired_distance_stdevs = cluster_stdev;
        minimizer_mapper.rescue_subgraph_stdevs = rescue_stdev;
        minimizer_mapper.rescue_seed_limit = rescue_seed_limit;
        minimizer_mapper.max_rescue_attempts = rescue_attempts;
        minimizer_mapper.rescue_algorithm = rescue_algorithm;

        minimizer_mapper.sample_name = sample_name;
        minimizer_mapper.read_group = read_group;

        // Work out the number of threads we will have
        size_t thread_count = omp_get_max_threads();

        // Set up counters per-thread for total reads mapped
        vector<size_t> reads_mapped_by_thread(thread_count, 0);
        
        // For timing, we may run one thread first and then switch to all threads. So track both start times.
        std::chrono::time_point<std::chrono::system_clock> first_thread_start;
        std::chrono::time_point<std::chrono::system_clock> all_threads_start;
        
        // We also time in terms of CPU time
        clock_t cpu_time_before;
        
        // We may also have access to perf stats.
        vector<int> perf_fds;
        
#ifdef __linux__
        // Set up a counter for executed instructions.
        // See <https://stackoverflow.com/a/64863392/402891>
        struct perf_event_attr perf_config;
        memset(&perf_config, 0, sizeof(struct perf_event_attr));
        perf_config.type = PERF_TYPE_HARDWARE;
        perf_config.size = sizeof(struct perf_event_attr);
        perf_config.config = PERF_COUNT_HW_INSTRUCTIONS;
        perf_config.exclude_kernel = 1;
        perf_config.exclude_hv = 1;
        
        perf_fds.resize(thread_count);
        
        perf_fds[omp_get_thread_num()] = perf_event_open(&perf_config, 0, -1, -1, 0);
        if (show_progress && perf_fds[omp_get_thread_num()] == -1) {
            int problem = errno;
            cerr << "Not counting CPU instructions because perf events are unavailable: " << strerror(problem) << endl;
            perf_fds.clear();
        }
        
        // Each OMP thread will call this to make sure perf is on.
        auto ensure_perf_for_thread = [&]() {
            if (!perf_fds.empty() && perf_fds[omp_get_thread_num()] == 0) {
                perf_fds[omp_get_thread_num()] = perf_event_open(&perf_config, 0, -1, -1, 0);
            }
        };
        
        // Main thread will call this to turn it off
        auto stop_perf_for_thread = [&]() {
            if (!perf_fds.empty() && perf_fds[omp_get_thread_num()] != 0) {
                ioctl(perf_fds[omp_get_thread_num()], PERF_EVENT_IOC_DISABLE, 0);
            }
        };
        
        // Main thread will call this when mapping starts to reset the counter.
        auto reset_perf_for_thread = [&]() {
            if (!perf_fds.empty() && perf_fds[omp_get_thread_num()] != 0) {
                ioctl(perf_fds[omp_get_thread_num()], PERF_EVENT_IOC_RESET, 0);
            }
        };
        
        // TODO: we won't count the output thread, but it will appear in CPU time!
#endif

        {
        
            // Look up all the paths we might need to surject to.
            vector<tuple<path_handle_t, size_t, size_t>> paths;
            if (hts_output) {
                // For htslib we need a non-empty list of paths.
                assert(path_position_graph != nullptr);
                paths = get_sequence_dictionary(ref_paths_name, *path_position_graph);
            }
            
            // Set up output to an emitter that will handle serialization and surjection.
            // Unless we want to discard all the alignments in which case do that.
            unique_ptr<AlignmentEmitter> alignment_emitter;
            if (discard_alignments) {
                alignment_emitter = make_unique<NullAlignmentEmitter>();
            } else {
                // We actually want to emit alignments.
                // Encode flags describing what we want to happen.
                int flags = ALIGNMENT_EMITTER_FLAG_NONE;
                if (prune_anchors) {
                    // When surjecting, do anchor pruning.
                    flags |= ALIGNMENT_EMITTER_FLAG_HTS_PRUNE_SUSPICIOUS_ANCHORS;
                }
                if (named_coordinates) {
                    // When not surjecting, use named segments instead of node IDs.
                    flags |= ALIGNMENT_EMITTER_FLAG_VG_USE_SEGMENT_NAMES;
                }
                
                // We send along the positional graph when we have it, and otherwise we send the GBWTGraph which is sufficient for GAF output.
                // TODO: What if we need both a positional graph and a NamedNodeBackTranslation???
                const HandleGraph* emitter_graph = path_position_graph ? (const HandleGraph*)path_position_graph : (const HandleGraph*)&(gbz->graph);
                
                alignment_emitter = get_alignment_emitter(output_filename, output_format,
                                                          paths, thread_count,
                                                          emitter_graph, flags);
            }
            
#ifdef USE_CALLGRIND
            // We want to profile the alignment, not the loading.
            CALLGRIND_START_INSTRUMENTATION;
#endif

            // Start timing overall mapping time now that indexes are loaded.
            first_thread_start = std::chrono::system_clock::now();
            cpu_time_before = clock();
            
#ifdef __linux__
            reset_perf_for_thread();
#endif

            if (interleaved || !fastq_filename_2.empty()) {
                //Map paired end from either one gam or fastq file or two fastq files

                // a buffer to hold read pairs that can't be unambiguously mapped before the fragment length distribution
                // is estimated
                // note: sufficient to have only one buffer because multithreading code enforces single threaded mode
                // during distribution estimation
                vector<pair<Alignment, Alignment>> ambiguous_pair_buffer;
                
                // Track whether the distribution was ready, so we can detect when it becomes ready and capture the all-threads start time.
                bool distribution_was_ready = false;

                // Define how to know if the paired end distribution is ready
                auto distribution_is_ready = [&]() {
                    bool is_ready = minimizer_mapper.fragment_distr_is_finalized();
                    if (is_ready && !distribution_was_ready) {
                        // It has become ready now.
                        distribution_was_ready = true;
                        
                        if (show_progress) {
                            // Report that it is now ready
                            #pragma omp critical (cerr)
                            {
                                cerr << "Using fragment length estimate: " << minimizer_mapper.get_fragment_length_mean() << " +/- " << minimizer_mapper.get_fragment_length_stdev() << endl;
                            }
                        }
                        
                        // Remember when now is.
                        all_threads_start = std::chrono::system_clock::now();
                    }
                    return is_ready;
                };
                
                // Define a way to force the distribution ready
                auto require_distribution_finalized = [&]() {
                    if (!minimizer_mapper.fragment_distr_is_finalized()){
                        cerr << "warning[vg::giraffe]: Finalizing fragment length distribution before reaching maximum sample size" << endl;
                        cerr << "                      mapped " << minimizer_mapper.get_fragment_length_sample_size() 
                             << " reads single ended with " << ambiguous_pair_buffer.size() << " pairs of reads left unmapped" << endl;
                        cerr << "                      mean: " << minimizer_mapper.get_fragment_length_mean() << ", stdev: " 
                             << minimizer_mapper.get_fragment_length_stdev() << endl;
                        minimizer_mapper.finalize_fragment_length_distr();
                    }
                };
                
                // Define how to align and output a read pair, in a thread.
                auto map_read_pair = [&](Alignment& aln1, Alignment& aln2) {
#ifdef __linux__
                    ensure_perf_for_thread();
#endif
                    toUppercaseInPlace(*aln1.mutable_sequence());
                    toUppercaseInPlace(*aln2.mutable_sequence());

                    pair<vector<Alignment>, vector<Alignment>> mapped_pairs = minimizer_mapper.map_paired(aln1, aln2, ambiguous_pair_buffer);
                    if (!mapped_pairs.first.empty() && !mapped_pairs.second.empty()) {
                        //If we actually tried to map this paired end
                        
                        // Work out whether it could be properly paired or not, if that is relevant.
                        // If we're here, let the read be properly paired in
                        // HTSlib terms no matter how far away it is in linear
                        // space (on the same contig), because it went into
                        // pair distribution estimation.
                        // TODO: The semantics are weird here. 0 means
                        // "properly paired at any distance" and
                        // numeric_limits<int64_t>::max() doesn't.
                        int64_t tlen_limit = 0;
                        if (hts_output && minimizer_mapper.fragment_distr_is_finalized()) {
                             tlen_limit = minimizer_mapper.get_fragment_length_mean() + 6 * minimizer_mapper.get_fragment_length_stdev();
                        }
                        // Emit it
                        alignment_emitter->emit_mapped_pair(std::move(mapped_pairs.first), std::move(mapped_pairs.second), tlen_limit);
                        // Record that we mapped a read.
                        reads_mapped_by_thread.at(omp_get_thread_num()) += 2;
                    }
                    
                    if (!minimizer_mapper.fragment_distr_is_finalized() && ambiguous_pair_buffer.size() >= MAX_BUFFERED_PAIRS) {
                        // We risk running out of memory if we keep this up.
                        cerr << "warning[vg::giraffe]: Encountered " << ambiguous_pair_buffer.size() << " ambiguously-paired reads before finding enough" << endl
                             << "                      unambiguously-paired reads to learn fragment length distribution. Are you sure" << endl
                             << "                      your reads are paired and your graph is not a hairball?" << endl;
                        require_distribution_finalized();
                    }
                };

                if (!gam_filename.empty()) {
                    // GAM file to remap
                    get_input_file(gam_filename, [&](istream& in) {
                        // Map pairs of reads to the emitter
                        vg::io::for_each_interleaved_pair_parallel_after_wait<Alignment>(in, map_read_pair, distribution_is_ready);
                    });
                } else if (!fastq_filename_2.empty()) {
                    //A pair of FASTQ files to map
                    fastq_paired_two_files_for_each_parallel_after_wait(fastq_filename_1, fastq_filename_2, map_read_pair, distribution_is_ready);


                } else if ( !fastq_filename_1.empty()) {
                    // An interleaved FASTQ file to map, map all its pairs in parallel.
                    fastq_paired_interleaved_for_each_parallel_after_wait(fastq_filename_1, map_read_pair, distribution_is_ready);
                }

                // Now map all the ambiguous pairs
                // Make sure fragment length distribution is finalized first.
                require_distribution_finalized();
                for (pair<Alignment, Alignment>& alignment_pair : ambiguous_pair_buffer) {

                    auto mapped_pairs = minimizer_mapper.map_paired(alignment_pair.first, alignment_pair.second);
                    // Work out whether it could be properly paired or not, if that is relevant.
                    int64_t tlen_limit = 0;
                    if (hts_output && minimizer_mapper.fragment_distr_is_finalized()) {
                         tlen_limit = minimizer_mapper.get_fragment_length_mean() + 6 * minimizer_mapper.get_fragment_length_stdev();
                    }
                    // Emit the read
                    alignment_emitter->emit_mapped_pair(std::move(mapped_pairs.first), std::move(mapped_pairs.second), tlen_limit);
                    // Record that we mapped a read.
                    reads_mapped_by_thread.at(omp_get_thread_num()) += 2;
                }
            } else {
                // Map single-ended

                // All the threads start at once.
                all_threads_start = first_thread_start;
            
                // Define how to align and output a read, in a thread.
                auto map_read = [&](Alignment& aln) {
#ifdef __linux__
                    ensure_perf_for_thread();
#endif
                    toUppercaseInPlace(*aln.mutable_sequence());
                
                    // Map the read with the MinimizerMapper.
                    minimizer_mapper.map(aln, *alignment_emitter);
                    // Record that we mapped a read.
                    reads_mapped_by_thread.at(omp_get_thread_num())++;
                };
                    
                if (!gam_filename.empty()) {
                    // GAM file to remap
                    get_input_file(gam_filename, [&](istream& in) {
                        // Open it and map all the reads in parallel.
                        vg::io::for_each_parallel<Alignment>(in, map_read);
                    });
                }
                
                if (!fastq_filename_1.empty()) {
                    // FASTQ file to map, map all its reads in parallel.
                    fastq_unpaired_for_each_parallel(fastq_filename_1, map_read);
                }
            }
        
        } // Make sure alignment emitter is destroyed and all alignments are on disk.
        
        // Now mapping is done
        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        clock_t cpu_time_after = clock();
#ifdef __linux__
        stop_perf_for_thread();
#endif
        
        // Compute wall clock elapsed
        std::chrono::duration<double> all_threads_seconds = end - all_threads_start;
        std::chrono::duration<double> first_thread_additional_seconds = all_threads_start - first_thread_start;
        
        // Compute CPU time elapsed
        double cpu_seconds = (cpu_time_after - cpu_time_before) / (double)CLOCKS_PER_SEC;
        
        // Compute instructions used
        long long total_instructions = 0;
        for (auto& perf_fd : perf_fds) {
            if (perf_fd > 0) {
                long long thread_instructions;
                if (read(perf_fd, &thread_instructions, sizeof(long long)) != sizeof(long long)) {
                    // Read failed for some reason.
                    cerr << "warning:[vg giraffe] Could not count CPU instructions executed" << endl;
                    thread_instructions = 0;
                }
                if (close(perf_fd)) {
                    int problem = errno;
                    cerr << "warning:[vg giraffe] Error closing perf event instruction counter: " << strerror(problem) << endl;
                }
                total_instructions += thread_instructions;
            }
        }
        
        // How many reads did we map?
        size_t total_reads_mapped = 0;
        for (auto& reads_mapped : reads_mapped_by_thread) {
            total_reads_mapped += reads_mapped;
        }
        
        // Compute speed (as reads per thread-second)
        double reads_per_second_per_thread = total_reads_mapped / (all_threads_seconds.count() * thread_count + first_thread_additional_seconds.count());
        // And per CPU second (including any IO threads)
        double reads_per_cpu_second = total_reads_mapped / cpu_seconds;
        double mega_instructions_per_read = total_instructions / (double)total_reads_mapped / 1E6;
        double mega_instructions_per_second = total_instructions / cpu_seconds / 1E6;
        
        if (show_progress) {
            // Log to standard error
            cerr << "Mapped " << total_reads_mapped << " reads across "
                << thread_count << " threads in "
                << all_threads_seconds.count() << " seconds with " 
                << first_thread_additional_seconds.count() << " additional single-threaded seconds." << endl;
            cerr << "Mapping speed: " << reads_per_second_per_thread
                << " reads per second per thread" << endl;
            
            cerr << "Used " << cpu_seconds << " CPU-seconds (including output)." << endl;
            cerr << "Achieved " << reads_per_cpu_second
                << " reads per CPU-second (including output)" << endl;
            
            if (total_instructions != 0) {
                cerr << "Used " << total_instructions << " CPU instructions (not including output)." << endl;
                cerr << "Mapping slowness: " << mega_instructions_per_read
                    << " M instructions per read at " << mega_instructions_per_second
                    << " M mapping instructions per inclusive CPU-second" << endl;
            }

            cerr << "Memory footprint: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GB" << endl;
        }
        
        
        if (report) {
            // Log output filename and mapping speed in reads/second/thread to report TSV
            report << output_filename << "\t" << reads_per_second_per_thread << endl;
        }
        
    });
        
    return 0;
}

// Register subcommand
static Subcommand vg_giraffe("giraffe", "fast haplotype-aware short read alignment", PIPELINE, 6, main_giraffe);
