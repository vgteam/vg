// score_main.cpp: defines a GAM to GAM annotation function

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include <subcommand.hpp>

#include "../alignment.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_score(char** argv) {
    cerr << "usage: " << argv[0] << " score aln.gam " << endl;
}

int main_score(int argc, char** argv) {

    if (argc == 2) {
        help_score(argv);
        exit(1);
    }

    int threads = 1;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'h':
        case '?':
            help_score(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    // We need to read the second argument first, so we can't use get_input_file with its free error checking.
    string read_file_name = get_input_file_name(optind, argc, argv);

    vector<size_t> mapq_counts ( 61, 0);
    vector<size_t> correct_counts (61, 0);;

    size_t total_read_count = 0;
    function<void(Alignment&)> count_mapqs = [&](Alignment& aln) {
        auto mapq = aln.mapping_quality();
        if (mapq) {
            total_read_count ++;
            mapq_counts[mapq] ++;
            if (aln.correctly_mapped()) {
                correct_counts[mapq]++;
            }
        }

    };


    ifstream truth_file_in(read_file_name);
    if (!truth_file_in) {
        cerr << "error[vg score]: Unable to read " << read_file_name << " when looking for true reads" << endl;
        exit(1);
    }
    //Get the mapqs
    vg::io::for_each_parallel(truth_file_in, count_mapqs);
   
    size_t accumulated_count = 0;
    size_t accumulated_correct_count= 0;
    float mapping_goodness_score = 0.0;
    for (int i = 60 ; i >= 0 ; i-- ) {
        accumulated_count += mapq_counts[i];
        accumulated_correct_count += correct_counts[i];
        double fraction_incorrect = accumulated_count == 0 ? 0.0 :
            (float) (accumulated_count - accumulated_correct_count) / (float) accumulated_count;
        fraction_incorrect = fraction_incorrect == 0.0 ? 1.0/(float)total_read_count : fraction_incorrect;

        mapping_goodness_score -= log10(fraction_incorrect) * mapq_counts[i];
    }


    
    cerr << mapping_goodness_score / total_read_count << endl;
    
    return 0;
}

// Register subcommand
static Subcommand vg_score("score", "score compared alignment", main_score);
