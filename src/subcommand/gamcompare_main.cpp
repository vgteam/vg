// gamcompare_main.cpp: defines a GAM to GAM annotation function

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include <subcommand.hpp>

#include "../alignment.hpp"
#include "../vg.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_gamcompare(char** argv) {
    cerr << "usage: " << argv[0] << " gamcompare aln.gam truth.gam >output.gam" << endl
         << endl
         << "options:" << endl
         << "    -r, --range N            distance within which to consider reads correct" << endl
         << "    -t, --threads N          number of threads to use" << endl;
}

int main_gamcompare(int argc, char** argv) {

    if (argc == 2) {
        help_gamcompare(argv);
    }

    int threads = 1;
    int64_t range = -1;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"range", required_argument, 0, 'r'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "ht:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'r':
            range = atoi(optarg);
            break;

        case 't':
            threads = atoi(optarg);
            omp_set_num_threads(threads);
            break;

        case 'h':
        case '?':
            help_gamcompare(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    string test_file_name = get_input_file_name(optind, argc, argv);
    string truth_file_name = get_input_file_name(optind, argc, argv);

    // We will collect all the truth positions
    string_hash_map<string, map<string ,vector<pair<size_t, bool> > > > true_positions;
    function<void(Alignment&)> record_truth = [&true_positions](Alignment& aln) {
        auto val = alignment_refpos_to_path_offsets(aln);
#pragma omp critical (truth_table)
        true_positions[aln.name()] = val;
    };
    if (truth_file_name == "-") {
        assert(test_file_name != "-");
        stream::for_each_parallel(std::cin, record_truth);
    } else {
        ifstream truth_file_in(truth_file_name);
        stream::for_each_parallel(truth_file_in, record_truth);
    }

    vector<Alignment> buf;
    function<void(Alignment&)> annotate_test = [&buf,&true_positions,&range](Alignment& aln) {
        auto f = true_positions.find(aln.name());
        if (f != true_positions.end()) {
            auto& true_position = f->second;
            alignment_set_distance_to_correct(aln, true_position);
            
            if (range != -1) {
                // We are flagging reads correct/incorrect.
                // It is correct if there is a path for its minimum distance and it is in range on that path.
                aln.set_correctly_mapped(aln.to_correct().name() != "" && aln.to_correct().offset() <= range);
            }
            
        }
#pragma omp critical (buf)
        {
            buf.push_back(aln);
            if (buf.size() > 1000) {
                write_alignments(std::cout, buf);
                buf.clear();
            }
        }
    };

    if (test_file_name == "-") {
        assert(truth_file_name != "-");
        stream::for_each_parallel(std::cin, annotate_test);
    } else {
        ifstream test_file_in(test_file_name);
        stream::for_each_parallel(test_file_in, annotate_test);
    }

    write_alignments(std::cout, buf);
    buf.clear();
    cout.flush();

    return 0;
}

// Register subcommand
static Subcommand vg_gamcompare("gamcompare", "compare alignment positions", main_gamcompare);
