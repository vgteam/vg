/** \file filter_main.cpp
 *
 * Defines the "vg filter" subcommand, which filters alignments.
 */


#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <list>
#include <fstream>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include "../readfilter.hpp"
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_filter(char** argv) {
    cerr << "usage: " << argv[0] << " filter [options] <alignment.gam> > out.gam" << endl
         << "Filter alignments by properties." << endl
         << endl
         << "options:" << endl
         << "    -M, --input-mp-alns        input is multipath alignments (GAMP) rather than GAM" << endl
         << "    -n, --name-prefix NAME     keep only reads with this prefix in their names [default='']" << endl
         << "    -N, --name-prefixes FILE   keep reads with names with one of many prefixes, one per nonempty line" << endl
         << "    -c, --exact-name           match read names exactly instead of by prefix" << endl
         << "    -a, --subsequence NAME     keep reads that contain this subsequence" << endl
         << "    -A, --subsequences FILE    keep reads that contain one of these subsequences, one per nonempty line" << endl
         << "    -p, --proper-pairs         keep reads that are annotated as being properly paired" << endl
         << "    -P, --only-mapped          keep reads that are mapped" << endl
         << "    -X, --exclude-contig REGEX drop reads with refpos annotations on contigs matching the given regex (may repeat)" << endl
         << "    -F, --exclude-feature NAME drop reads with the given feature in the \"features\" annotation (may repeat)" << endl
         << "    -s, --min-secondary N      minimum score to keep secondary alignment" << endl
         << "    -r, --min-primary N        minimum score to keep primary alignment" << endl
         << "    -L, --max-length N         drop reads with length > N" << endl
         << "    -O, --rescore              re-score reads using default parameters and only alignment information" << endl
         << "    -f, --frac-score           normalize score based on length" << endl
         << "    -u, --substitutions        use substitution count instead of score" << endl
         << "    -o, --max-overhang N       drop reads whose alignments begin or end with an insert > N [default=99999]" << endl
         << "    -m, --min-end-matches N    drop reads that don't begin with at least N matches on each end" << endl
         << "    -S, --drop-split           remove split reads taking nonexistent edges" << endl
         << "    -x, --xg-name FILE         use this xg index or graph (required for -S and -D)" << endl
         << "    -v, --verbose              print out statistics on numbers of reads dropped by what." << endl
         << "    -V, --no-output            print out statistics (as above) but do not write out filtered GAM." << endl
         << "    -q, --min-mapq N           drop alignments with mapping quality < N" << endl
         << "    -E, --repeat-ends N        drop reads with tandem repeat (motif size <= 2N, spanning >= N bases) at either end" << endl
         << "    -D, --defray-ends N        clip back the ends of reads that are ambiguously aligned, up to N bases" << endl
         << "    -C, --defray-count N       stop defraying after N nodes visited (used to keep runtime in check) [default=99999]" << endl
         << "    -d, --downsample S.P       drop all but the given portion 0.P of the reads. S may be an integer seed as in SAMtools" << endl
         << "    -i, --interleaved          assume interleaved input. both ends will be dropped if either fails filter" << endl
         << "    -I, --interleaved-all      assume interleaved input. both ends will be dropped if *both* fail filters" << endl
         << "    -b, --min-base-quality Q:F drop reads with where fewer than fraction F bases have base quality >= PHRED score Q." << endl
         << "    -U, --complement           apply the complement of the filter implied by the other arguments." << endl
         << "    -t, --threads N            number of threads [1]" << endl;
}

int main_filter(int argc, char** argv) {

    if (argc <= 2) {
        help_filter(argv);
        return 1;
    }
    
    bool input_gam = true;
    vector<string> name_prefixes;
    bool exact_name = false;
    vector<regex> excluded_refpos_contigs;
    unordered_set<string> excluded_features;
    vector<string> subsequences;
    bool set_min_primary = false;
    double min_primary;
    bool set_min_secondary = false;
    double min_secondary;
    size_t max_length = std::numeric_limits<size_t>::max();
    bool rescore = false;
    bool frac_score = false;
    bool sub_score = false;
    bool set_max_overhang = false;
    int max_overhang;
    bool set_min_end_matches = false;
    int min_end_matches;
    bool drop_split = false;
    bool set_min_mapq = false;
    int min_mapq;
    bool verbose = false;
    bool write_output = true;
    bool set_repeat_size = false;
    int repeat_size;
    bool set_defray_length = false;
    int defray_length;
    bool set_defray_count = false;
    int defray_count;
    bool set_downsample = false;
    uint64_t seed;
    double downsample_probability;
    bool interleaved = false;
    bool filter_on_all = false;
    bool set_min_base_quality = false;
    int min_base_quality;
    double min_base_quality_fraction;
    bool complement_filter = false;
    bool only_proper_pairs = false;
    bool only_mapped = false;

    // What XG index, if any, should we load to support the other options?
    string xg_name;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"input-mp-alns", no_argument, 0, 'M'},
                {"name-prefix", required_argument, 0, 'n'},
                {"name-prefixes", required_argument, 0, 'N'},
                {"exact-name", no_argument, 0, 'c'},
                {"subsequence", required_argument, 0, 'a'},
                {"subsequences", required_argument, 0, 'A'},
                {"proper-pairs", no_argument, 0, 'p'},
                {"only-mapped", no_argument, 0, 'P'},
                {"exclude-contig", required_argument, 0, 'X'},
                {"exclude-feature", required_argument, 0, 'F'},
                {"min-secondary", required_argument, 0, 's'},
                {"min-primary", required_argument, 0, 'r'},
                {"max-length", required_argument, 0, 'L'},
                {"rescore", no_argument, 0, 'O'},
                {"frac-score", required_argument, 0, 'f'},
                {"substitutions", required_argument, 0, 'u'},
                {"max-overhang", required_argument, 0, 'o'},
                {"min-end-matches", required_argument, 0, 'm'},
                {"drop-split",  no_argument, 0, 'S'},
                {"xg-name", required_argument, 0, 'x'},
                {"verbose",  no_argument, 0, 'v'},
                {"min-mapq", required_argument, 0, 'q'},
                {"repeat-ends", required_argument, 0, 'E'},
                {"defray-ends", required_argument, 0, 'D'},
                {"defray-count", required_argument, 0, 'C'},
                {"downsample", required_argument, 0, 'd'},
                {"interleaved", no_argument, 0, 'i'},
                {"interleaved-all", no_argument, 0, 'I'},
                {"min-base-quality", required_argument, 0, 'b'},
                {"complement", no_argument, 0, 'U'},
                {"threads", required_argument, 0, 't'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "Mn:N:ca:A:pPX:F:s:r:L:Od:e:fauo:m:Sx:vVq:E:D:C:d:iIb:Ut:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'M':
            input_gam = false;
            break;
        case 'n':
            name_prefixes.push_back(optarg);
            break;
        case 'N':
            get_input_file(optarg, [&](istream& in) {
                // Parse the input file right here in the option parsing.
                for (string line; getline(in, line);) {
                    // For each line
                    if (line.empty()) {
                        // No empty lines
                        break;
                    }
                    name_prefixes.push_back(line);
                }
            });
            break;
        case 'c':
            exact_name = true;
            break;
        case 'a':
            subsequences.push_back(optarg);
            break;
        case 'A':
            get_input_file(optarg, [&](istream& in) {
                // Parse the input file right here in the option parsing.
                for (string line; getline(in, line);) {
                    // For each line
                    if (line.empty()) {
                        // No empty lines
                        break;
                    }
                    subsequences.push_back(line);
                }
            });
            break;
        case 'p':
            only_proper_pairs = true;
            break;
        case 'P':
            only_mapped = true;
            break;
        case 'X':
            excluded_refpos_contigs.push_back(parse<std::regex>(optarg));
            break;
        case 'F':
            excluded_features.insert(optarg);
            break;
        case 's':
            set_min_secondary = true;
            min_secondary = parse<double>(optarg);
            break;
        case 'r':
            set_min_primary = true;
            min_primary = parse<double>(optarg);
            break;
        case 'L':
            max_length = parse<size_t>(optarg);
            break;
        case 'O':
            rescore = true;
            break;
        case 'f':
            frac_score = true;
            break;
        case 'u':
            sub_score = true;
            break;
        case 'o':
            set_max_overhang = true;
            max_overhang = parse<int>(optarg);
            break;
        case 'm':
            set_min_end_matches = true;
            min_end_matches = parse<int>(optarg);
            break;            
        case 'S':
            drop_split = true;
        case 'x':
            xg_name = optarg;
            break;
        case 'q':
            set_min_mapq = true;
            min_mapq = parse<double>(optarg);
            break;
        case 'v':
            verbose = true;
            break;
        case 'V':
            verbose = true;
            write_output = false;
            break;
        case 'E':
            set_repeat_size = true;
            repeat_size = parse<int>(optarg);
            break;
        case 'D':
            set_defray_length = true;
            defray_length = parse<int>(optarg);
            break;
        case 'C':
            set_defray_count = true;
            defray_count = parse<int>(optarg);
            break;
        case 'd':
            {
                set_downsample = true;
                // We need to split out the seed and the probability in S.P
                string opt_string(optarg);
                
                if (opt_string != "1") {
                    // We need to parse
                    auto point = opt_string.find('.');
                    
                    if (point == -1) {
                        cerr << "error: no decimal point in seed/probability " << opt_string << endl;
                        exit(1);
                    }
                    
                    // Everything including and after the decimal point is the probability
                    auto prob_string = opt_string.substr(point);
                    downsample_probability = parse<double>(prob_string);
                    
                    // Everything before that is the seed
                    auto seed_string = opt_string.substr(0, point);
                    // Use the 0 seed by default even if no actual seed shows up
                    if (seed_string != "") {
                        // If there was a seed (even 0), parse it
                        seed = parse<int32_t>(seed_string);
                    }
                }
            }
            break;
        case 'i':
            interleaved = true;
            break;
        case 'I':
            interleaved = true;
            filter_on_all = true;
            break;
        case 'b':
            {
                set_min_base_quality = true;
                vector<string> parts = split_delims(string(optarg), ":");
                if (parts.size() != 2) {
                    cerr << "[vg filter] Error: -b expects value in form of <INT>:<FLOAT>" << endl;
                    return 1;
                }
                min_base_quality = parse<int>(parts[0]);
                min_base_quality_fraction = parse<double>(parts[1]);
                if (min_base_quality_fraction < 0 || min_base_quality_fraction > 1) {
                    cerr << "[vg filter] Error: second part of -b input must be between 0 and 1" << endl;
                    return 1;
                }
            }
            break;
        case 'U':
            complement_filter = true;
            break;
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_filter(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }
    
    if (optind >= argc) {
        help_filter(argv);
        return 1;
    }

    // What should our return code be?
    int error_code = 0;
    
    // Sort the prefixes for reads we will accept, for efficient search
    sort(name_prefixes.begin(), name_prefixes.end());
    
     // If the user gave us an XG index, we probably ought to load it up.
    PathPositionHandleGraph* xindex = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::ReferencePathOverlayHelper overlay_helper;
    if (!xg_name.empty()) {
        // read the xg index
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        xindex = overlay_helper.apply(path_handle_graph.get());
    }
    
    // template lambda to set parameters
    auto set_params = [&](auto& filter) {
        filter.name_prefixes = name_prefixes;
        filter.exact_name = exact_name;
        filter.subsequences = subsequences;
        filter.excluded_refpos_contigs = excluded_refpos_contigs;
        filter.excluded_features = excluded_features;
        if (set_min_secondary) {
            filter.min_secondary = min_secondary;
        }
        if (set_min_primary) {
            filter.min_primary = min_primary;
        }
        filter.max_length = max_length;
        filter.rescore = rescore;
        filter.frac_score = frac_score;
        filter.sub_score = sub_score;
        if (set_max_overhang){
            filter.max_overhang = max_overhang;
        }
        if (set_min_end_matches) {
            filter.min_end_matches = min_end_matches;
        }
        filter.drop_split = drop_split;
        if (set_min_mapq) {
            filter.min_mapq = min_mapq;
        }
        filter.verbose = verbose;
        filter.write_output = write_output;
        if (set_repeat_size) {
            filter.repeat_size = repeat_size;
        }
        if (set_defray_length){
            filter.defray_length = defray_length;
        }
        if (set_defray_count) {
            filter.defray_count = defray_count;
        }
        if (set_downsample) {
            filter.downsample_probability = downsample_probability;
            if (seed != 0) {
                // We want a nonempty mask.
                
                // Use the C rng like Samtools does to get a mask.
                // See https://github.com/samtools/samtools/blob/60138c42cf04c5c473dc151f3b9ca7530286fb1b/sam_view.c#L298-L299
                srand(seed);
                filter.downsample_seed_mask = rand();
            }
        }
        filter.only_proper_pairs = only_proper_pairs;
        filter.only_mapped = only_mapped;
        filter.interleaved = interleaved;
        filter.filter_on_all = filter_on_all;
        if (set_min_base_quality) {
            filter.min_base_quality = min_base_quality;
            filter.min_base_quality_fraction = min_base_quality_fraction;
        }
        filter.complement_filter = complement_filter;
        filter.threads = get_thread_count();
        filter.graph = xindex;
    };
    
    // Read in the alignments and filter them.
    get_input_file(optind, argc, argv, [&](istream& in) {
        // Open up the alignment stream
        
        // Read in the alignments and filter them.
        if (input_gam) {
            ReadFilter<Alignment> filter;
            set_params(filter);
            error_code = filter.filter(&in);
        }
        else {
            ReadFilter<MultipathAlignment> filter;
            set_params(filter);
            error_code = filter.filter(&in);
        }
    });

    return error_code;
}

// Register subcommand
static Subcommand vg_filter("filter", "filter reads", main_filter);

