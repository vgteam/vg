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
#include <vg/io/stream.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_gamcompare(char** argv) {
    cerr << "usage: " << argv[0] << " gamcompare aln.gam truth.gam >output.gam" << endl
         << endl
         << "options:" << endl
         << "    -r, --range N            distance within which to consider reads correct" << endl
         << "    -T, --tsv                output TSV (correct, mq, aligner, read) compatible with plot-qq.R instead of GAM" << endl
         << "    -a, --aligner            aligner name for TSV output [\"vg\"]" << endl
         << "    -s, --score-alignment    get a correctness score of the alignment (higher is better)" << endl
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
    bool score_alignment = false;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"range", required_argument, 0, 'r'},
            {"tsv", no_argument, 0, 'T'},
            {"aligner", required_argument, 0, 'a'},
            {"score-alignment", no_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hr:Ta:st:",
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

        case 's':
            score_alignment = true;
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

    // We need to read the second argument first, so we can't use get_input_file with its free error checking.
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
        // Read truth fropm standard input, if it looks good.
        if (test_file_name == "-") {
            cerr << "error[vg gamcompare]: Standard input can only be used for truth or test file, not both" << endl;
            exit(1);
        }
        if (!std::cin) {
            cerr << "error[vg gamcompare]: Unable to read standard input when looking for true reads" << endl;
            exit(1);
        }
        vg::io::for_each_parallel(std::cin, record_truth);
    } else {
        // Read truth from this file, if it looks good.
        ifstream truth_file_in(truth_file_name);
        if (!truth_file_in) {
            cerr << "error[vg gamcompare]: Unable to read " << truth_file_name << " when looking for true reads" << endl;
            exit(1);
        }
        vg::io::for_each_parallel(truth_file_in, record_truth);
    }
    if (score_alignment && range == -1) {
        cerr << "error[vg gamcompare]: Score-alignment requires range" << endl;
        exit(1);
    }

    // We have a buffered emitter for annotated alignments, if we're not outputting text
    std::unique_ptr<vg::io::ProtobufEmitter<Alignment>> emitter;
    if (!output_tsv) {
        emitter = std::unique_ptr<vg::io::ProtobufEmitter<Alignment>>(new vg::io::ProtobufEmitter<Alignment>(cout));
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
   
    // We want to count correct reads
    vector<size_t> correct_counts(get_thread_count(), 0);

    //Get stats for calculating the score
    vector<size_t> read_count_by_thread (get_thread_count(), 0);
    vector<vector<size_t>> mapq_count_by_thread (get_thread_count());
    vector<vector<size_t>> correct_count_by_mapq_by_thread(get_thread_count());
    for (size_t i = 0 ; i < get_thread_count() ; i++) {
        mapq_count_by_thread[i].resize(61, 0);
        correct_count_by_mapq_by_thread[i].resize(61,0);
    }
   
    // This function annotates every read with distance and correctness, and batch-outputs them.
    function<void(Alignment&)> annotate_test = [&](Alignment& aln) {
        auto f = true_positions.find(aln.name());
        if (f != true_positions.end()) {
            auto& true_position = f->second;
            alignment_set_distance_to_correct(aln, true_position);
            
            if (range != -1) {
                // We are flagging reads correct/incorrect.
                // It is correct if there is a path for its minimum distance and it is in range on that path.
                bool correctly_mapped = (aln.to_correct().name() != "" && aln.to_correct().offset() <= range);
                
                // Annotate it as such
                aln.set_correctly_mapped(correctly_mapped);
                
                if (correctly_mapped) {
                    correct_counts.at(omp_get_thread_num()) += 1;
                }
                auto mapq = aln.mapping_quality();
                if (mapq) {
                    if (mapq >= mapq_count_by_thread.at(omp_get_thread_num()).size()) {
                        mapq_count_by_thread.at(omp_get_thread_num()).resize(mapq+1);
                        correct_count_by_mapq_by_thread.at(omp_get_thread_num()).resize(mapq+1);
                    }

                    read_count_by_thread.at(omp_get_thread_num()) += 1;
                    mapq_count_by_thread.at(omp_get_thread_num()).at(mapq) += 1;
                    if (correctly_mapped) {
                        correct_count_by_mapq_by_thread.at(omp_get_thread_num()).at(mapq) += 1;
                    }
                }
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
        if (!std::cin) {
            cerr << "error[vg gamcompare]: Unable to read standard input when looking for reads under test" << endl;
            exit(1);
        }
        vg::io::for_each_parallel(std::cin, annotate_test);
    } else {
        ifstream test_file_in(test_file_name);
        if (!test_file_in) {
            cerr << "error[vg gamcompare]: Unable to read " << test_file_name << " when looking for reads under test" << endl;
            exit(1);
        }
        vg::io::for_each_parallel(test_file_in, annotate_test);
    }

    if (output_tsv) {
        // Save whatever's in the buffer at the end.
        flush_text_buffer();
    }

    
    if (range != -1) {
        // We are flagging reads correct/incorrect. So report the total correct.
        size_t total_correct = 0;
        for (auto& count : correct_counts) {
            total_correct += count;
        }
        
        cerr << total_correct << " reads correct" << endl;
    }

    if (score_alignment) {
        //Get a goodness score of the alignment that takes into account correctness and mapq calibration
        size_t total_reads = 0;
        vector<size_t> mapq_count (61, 0);
        vector<size_t> correct_count_by_mapq (61, 0);
        for (size_t i = 0 ; i < get_thread_count() ; i++) {
            total_reads += read_count_by_thread.at(i);
            for (size_t mq = 0 ; mq < mapq_count_by_thread.at(i).size() ; mq++) {
                if (mq >= mapq_count.size()) {
                    mapq_count.resize(mq+1);
                    correct_count_by_mapq.resize(mq+1);
                }

                mapq_count.at(mq) += mapq_count_by_thread.at(i).at(mq);
                correct_count_by_mapq.at(mq) += correct_count_by_mapq_by_thread.at(i).at(mq);
            }
        }
        size_t accumulated_count = 0;
        size_t accumulated_correct_count = 0;
        float mapping_goodness_score = 0.0;
        for (int i = mapq_count.size() ; i >= 0 ; i--) {
            accumulated_count += mapq_count[i];
            accumulated_correct_count += correct_count_by_mapq[i];
            double fraction_incorrect = accumulated_count == 0 ? 0.0 :
                (float) (accumulated_count - accumulated_correct_count) / (float) accumulated_count;
            fraction_incorrect = fraction_incorrect == 0.0 ? 1.0/ (float) total_reads : fraction_incorrect;
            mapping_goodness_score -= log10(fraction_incorrect) * mapq_count[i];
        }
        cerr << "mapping goodness score: " << mapping_goodness_score / total_reads << endl;

    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_gamcompare("gamcompare", "compare alignment positions", main_gamcompare);
