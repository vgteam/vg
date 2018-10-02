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
#include "../stream.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_gamcompare(char** argv) {
    cerr << "usage: " << argv[0] << " gamcompare aln.gam truth.gam >output.gam" << endl
         << endl
         << "options:" << endl
         << "    -r, --range N            distance within which to consider reads correct" << endl
         << "    -T, --tsv                output TSV (correct, mq, aligner, read) comaptible with plot-qq.R instead of GAM" << endl
         << "    -a, --aligner            aligner name for TSV output [\"vg\"]" << endl
         << "    -t, --threads N          number of threads to use" << endl;
}

int main_gamcompare(int argc, char** argv) {

    if (argc == 2) {
        help_gamcompare(argv);
        exit(1);
    }

    int threads = 1;
    int64_t range = -1;
    bool output_tsv = false;
    string aligner_name = "vg";

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"range", required_argument, 0, 'r'},
            {"tsv", no_argument, 0, 'T'},
            {"aligner", required_argument, 0, 'a'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hr:Ta:t:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'r':
            range = parse<int>(optarg);
            break;
            
        case 'T':
            output_tsv = true;
            break;
            
        case 'a':
            aligner_name = optarg;
            break;

        case 't':
            threads = parse<int>(optarg);
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

    // We have a buffered emitter for annotated alignments, if we're not outputting text
    std::unique_ptr<stream::ProtobufEmitter<Alignment>> emitter;
    if (!output_tsv) {
        emitter = std::unique_ptr<stream::ProtobufEmitter<Alignment>>(new stream::ProtobufEmitter<Alignment>(cout));
    }
    
    // We have an ordinary buffer we use for text output
    vector<Alignment> text_buffer;
    
    // We have an output function to dump all the reads in the text buffer in TSV
    auto flush_text_buffer = [&text_buffer,&output_tsv,&aligner_name]() {
        // We print exactly one header line.
        static bool header_printed = false;
        // Output TSV to standard out in the format plot-qq.R needs.
        if (!header_printed) {
            // It needs a header
            cout << "correct\tmq\taligner\tread" << endl;
            header_printed = true;
        }
        
        for (auto& aln : text_buffer) {
            // Dump each alignment
            cout << (aln.correctly_mapped() ? "1" : "0") << "\t";
            cout << aln.mapping_quality() << "\t";
            cout << aligner_name << "\t";
            cout << aln.name() << endl;
        }
        text_buffer.clear();
    };
    
    // This function annotates every read with distance and correctness, and batch-outputs them.
    function<void(Alignment&)> annotate_test = [&](Alignment& aln) {
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
#pragma omp critical
        {
            if (output_tsv) {
                text_buffer.emplace_back(std::move(aln));
                if (text_buffer.size() > 1000) {
                    flush_text_buffer();
                }
            } else {
                emitter->write(std::move(aln));
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

    if (output_tsv) {
        // Save whatever's in the buffer at the end.
        flush_text_buffer();
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_gamcompare("gamcompare", "compare alignment positions", main_gamcompare);
