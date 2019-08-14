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

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_filter(char** argv) {
    cerr << "usage: " << argv[0] << " filter [options] <alignment.gam> > out.gam" << endl
         << "Filter alignments by properties." << endl
         << endl
         << "options:" << endl
         << "    -n, --name-prefix NAME     keep only reads with this prefix in their names [default='']" << endl
         << "    -N, --name-prefixes FILE   keep reads with names with one of many prefixes, one per nonempty line" << endl
         << "    -X, --exclude-contig REGEX drop reads with refpos annotations on contigs matching the given regex (may repeat)" << endl
         << "    -F, --exclude-feature NAME drop reads with the given feature in the \"features\" annotation (may repeat)" << endl
         << "    -s, --min-secondary N      minimum score to keep secondary alignment [default=0]" << endl
         << "    -r, --min-primary N        minimum score to keep primary alignment [default=0]" << endl
         << "    -O, --rescore              re-score reads using default parameters and only alignment information" << endl
         << "    -f, --frac-score           normalize score based on length" << endl
         << "    -u, --substitutions        use substitution count instead of score" << endl
         << "    -o, --max-overhang N       filter reads whose alignments begin or end with an insert > N [default=99999]" << endl
         << "    -m, --min-end-matches N    filter reads that don't begin with at least N matches on each end" << endl
         << "    -S, --drop-split           remove split reads taking nonexistent edges" << endl
         << "    -x, --xg-name FILE         use this xg index (required for -S and -D)" << endl
         << "    -A, --append-regions       append to alignments created with -RB" << endl
         << "    -v, --verbose              print out statistics on numbers of reads filtered by what." << endl
         << "    -V, --no-output            print out statistics (as above) but do not write out filtered GAM." << endl
         << "    -q, --min-mapq N           filter alignments with mapping quality < N" << endl
         << "    -E, --repeat-ends N        filter reads with tandem repeat (motif size <= 2N, spanning >= N bases) at either end" << endl
         << "    -D, --defray-ends N        clip back the ends of reads that are ambiguously aligned, up to N bases" << endl
         << "    -C, --defray-count N       stop defraying after N nodes visited (used to keep runtime in check) [default=99999]" << endl
         << "    -d, --downsample S.P       filter out all but the given portion 0.P of the reads. S may be an integer seed as in SAMtools" << endl
         << "    -i, --interleaved          assume interleaved input. both ends will be filtered out if either fails filter" << endl
         << "    -I, --interleaved-all      assume interleaved input. both ends will be filtered out if *both* fail filters" << endl
         << "    -b, --min-base-quality Q:F filter reads with where fewer than fraction F bases have base quality >= PHRED score Q." << endl
         << "    -U, --complement           apply the complement of the filter implied by the other arguments." << endl
         << "    -t, --threads N            number of threads [1]" << endl;
}

int main_filter(int argc, char** argv) {

    if (argc <= 2) {
        help_filter(argv);
        return 1;
    }

    // This is the better design for a subcommand: we have a class that
    // implements it and encapsulates all the default parameters, and then we
    // just feed in overrides in the option parsing code. Thsi way we don't have
    // multiple defaults all over the place.
    ReadFilter filter;

    // What XG index, if any, should we load to support the other options?
    string xg_name;

    int c;
    optind = 2; // force optind past command positional arguments
    while (true) {
        static struct option long_options[] =
            {
                {"name-prefix", required_argument, 0, 'n'},
                {"name-prefixes", required_argument, 0, 'N'},
                {"exclude-contig", required_argument, 0, 'X'},
                {"exclude-feature", required_argument, 0, 'F'},
                {"min-secondary", required_argument, 0, 's'},
                {"min-primary", required_argument, 0, 'r'},
                {"rescore", no_argument, 0, 'O'},
                {"frac-score", required_argument, 0, 'f'},
                {"substitutions", required_argument, 0, 'u'},
                {"max-overhang", required_argument, 0, 'o'},
                {"min-end-matches", required_argument, 0, 'm'},
                {"drop-split",  no_argument, 0, 'S'},
                {"xg-name", required_argument, 0, 'x'},
                {"append-regions", no_argument, 0, 'A'},
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
        c = getopt_long (argc, argv, "n:N:X:F:s:r:Od:e:fauo:m:Sx:AvVq:E:D:C:d:iIb:Ut:",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'n':
            filter.name_prefixes.push_back(optarg);
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
                    filter.name_prefixes.push_back(line);
                }
            });
            break;
        case 'X':
            filter.excluded_refpos_contigs.push_back(parse<std::regex>(optarg));
            break;
        case 'F':
            filter.excluded_features.insert(optarg);
            break;
        case 's':
            filter.min_secondary = parse<double>(optarg);
            break;
        case 'r':
            filter.min_primary = parse<double>(optarg);
            break;
        case 'O':
            filter.rescore = true;
            break;
        case 'f':
            filter.frac_score = true;
            break;
        case 'u':
            filter.sub_score = true;
            break;
        case 'o':
            filter.max_overhang = parse<int>(optarg);
            break;
        case 'm':
            filter.min_end_matches = parse<int>(optarg);
            break;            
        case 'S':
            filter.drop_split = true;
        case 'x':
            xg_name = optarg;
            break;
        case 'A':
            filter.append_regions = true;
            break;
        case 'q':
            filter.min_mapq = parse<double>(optarg);
            break;
        case 'v':
            filter.verbose = true;
            break;
        case 'V':
            filter.verbose = true;
            filter.write_output = false;
            break;
        case 'E':
            filter.repeat_size = parse<int>(optarg);
            break;
        case 'D':
            filter.defray_length = parse<int>(optarg);
            break;
        case 'C':
            filter.defray_count = parse<int>(optarg);
            break;
        case 'd':
            {
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
                    filter.downsample_probability = parse<double>(prob_string);
                    
                    // Everything before that is the seed
                    auto seed_string = opt_string.substr(0, point);
                    // Use the 0 seed by default even if no actual seed shows up
                    int32_t seed = 0;
                    if (seed_string != "") {
                        // If there was a seed (even 0), parse it
                        seed = parse<int32_t>(seed_string);
                    }
                    
                    if (seed != 0) {
                        // We want a nonempty mask.
                        
                        // Use the C rng like Samtools does to get a mask.
                        // See https://github.com/samtools/samtools/blob/60138c42cf04c5c473dc151f3b9ca7530286fb1b/sam_view.c#L298-L299
                        srand(seed);
                        filter.downsample_seed_mask = rand();
                    }
                }
            }
            break;
        case 'i':
            filter.interleaved = true;
            break;
        case 'I':
            filter.interleaved = true;
            filter.filter_on_all = true;
            break;
        case 'b':
            {
                vector<string> parts = split_delims(string(optarg), ":");
                if (parts.size() != 2) {
                    cerr << "[vg filter] Error: -b expects value in form of <INT>:<FLOAT>" << endl;
                    return 1;
                }
                filter.min_base_quality = parse<int>(parts[0]);
                filter.min_base_quality_fraction = parse<double>(parts[1]);
                if (filter.min_base_quality_fraction < 0 || filter.min_base_quality_fraction > 1) {
                    cerr << "[vg filter] Error: second part of -b input must be between 0 and 1" << endl;
                    return 1;
                }
            }
            break;
        case 'U':
            filter.complement_filter = true;
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

    filter.threads = get_thread_count();
    
    if (optind >= argc) {
        help_filter(argv);
        return 1;
    }

    // What should our return code be?
    int error_code = 0;
    
    // Sort the prefixes for reads we will accept, for efficient search
    sort(filter.name_prefixes.begin(), filter.name_prefixes.end());
    
     // If the user gave us an XG index, we probably ought to load it up.
    unique_ptr<PathPositionHandleGraph> xindex;
    if (!xg_name.empty()) {
        // read the xg index
        xindex = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_name);
    }
    filter.graph = xindex.get();
    
    // Read in the alignments and filter them.
    get_input_file(optind, argc, argv, [&](istream& in) {
        // Open up the alignment stream
        
        // Read in the alignments and filter them.
        error_code = filter.filter(&in);
    });

    return error_code;
}

// Register subcommand
static Subcommand vg_filter("filter", "filter reads", main_filter);

