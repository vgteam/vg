// gamcompare_main.cpp: defines a GAM to GAM annotation function

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include "subcommand.hpp"

#include "../alignment.hpp"
#include "../min_distance.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_gamcompare(char** argv) {
    cerr << "usage: " << argv[0] << " gamcompare aln.gam truth.gam >output.gam" << endl
         << endl
         << "options:" << endl
         << "    -d, --distance-index FILE  use distances from this distance index instead of path position annotations" << endl
         << "    -r, --range N              distance within which to consider reads correct" << endl
         << "    -T, --tsv                  output TSV (correct, mq, aligner, read) compatible with plot-qq.R instead of GAM" << endl
         << "    -a, --aligner              aligner name for TSV output [\"vg\"]" << endl
         << "    -s, --score-alignment      get a correctness score of the alignment (higher is better)" << endl
         << "    -t, --threads N            number of threads to use" << endl;
}

// A gapless alignment between a read and a single node.
struct MappingRun {
    pos_t  start; // Starting position in the graph.
    size_t read_offset; // Starting position in the read.
    size_t length; // Length of the alignment.

    size_t limit() const {
        return this->read_offset + this->length;
    }

    // Get the graph position at read offset `offset >= this->read_offset`.
    pos_t pos_at(size_t offset) const {
        pos_t result = this->start;
        get_offset(result) += offset - this->read_offset;
        return result;
    }
};

// Returns the maximal MappingRuns for the alignment.
std::vector<MappingRun> base_mappings(const Alignment& aln) {
    std:vector<MappingRun> result;
    size_t read_offset = 0;
    const Path& path = aln.path();
    for (size_t i = 0; i < path.mapping_size(); i++) {
        const Mapping& mapping = path.mapping(i);
        pos_t start = make_pos_t(mapping.position());
        size_t length = 0; // Number of consecutive matches/mismatches.
        for (size_t j = 0; j < mapping.edit_size(); j++) {
            const Edit& edit = mapping.edit(j);
            if (edit.from_length() == edit.to_length()) {
                length += edit.to_length();
            } else {
                if (length > 0) {
                    result.push_back({ start, read_offset, length });
                    get_offset(start) += length;
                    read_offset += length;
                    length = 0;
                }
                get_offset(start) += edit.from_length();
                read_offset += edit.to_length();
            }
        }
        if (length > 0) {
            result.push_back({ start, read_offset, length });
            read_offset += length;
        }
    }
    return result;
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
    string distance_name;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"distance-index", required_argument, 0, 'd'},
            {"range", required_argument, 0, 'r'},
            {"tsv", no_argument, 0, 'T'},
            {"aligner", required_argument, 0, 'a'},
            {"score-alignment", no_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hd:r:Ta:st:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'r':
            range = parse<int>(optarg);
            break;

        case 'd':
            distance_name = optarg;
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

    // True path positions. For each alignment name, store a mapping from reference path names
    // to sets of (sequence offset, is_reverse). There is usually either one position per
    // alignment or one position per node.
    string_hash_map<string, map<string, vector<pair<size_t, bool> > > > true_path_positions;
    function<void(Alignment&)> record_path_positions = [&true_path_positions](Alignment& aln) {
        auto val = alignment_refpos_to_path_offsets(aln);
#pragma omp critical (truth_table)
        true_path_positions[aln.name()] = val;
    };

    // True graph positions. For each alignment name, we find the maximal read intervals that correspond
    // to a gapless alignment between the read and a single node.
    string_hash_map<string, std::vector<MappingRun>> true_graph_positions;
    function<void(Alignment&)> record_graph_positions = [&true_graph_positions](Alignment& aln) {
        if (aln.path().mapping_size() > 0) {
#pragma omp critical (truth_table)
            true_graph_positions[aln.name()] = base_mappings(aln);
        }
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
        if (distance_name.empty()) {
            vg::io::for_each_parallel(std::cin, record_path_positions);
        } else {
            vg::io::for_each_parallel(std::cin, record_graph_positions);
        }
    } else {
        // Read truth from this file, if it looks good.
        ifstream truth_file_in(truth_file_name);
        if (!truth_file_in) {
            cerr << "error[vg gamcompare]: Unable to read " << truth_file_name << " when looking for true reads" << endl;
            exit(1);
        }
        if (distance_name.empty()) {
            vg::io::for_each_parallel(truth_file_in, record_path_positions);
        } else {
            vg::io::for_each_parallel(truth_file_in, record_graph_positions);
        }
    }
    if (score_alignment && range == -1) {
        cerr << "error[vg gamcompare]: Score-alignment requires range" << endl;
        exit(1);
    }

    // Load the distance index.
    std::unique_ptr<MinimumDistanceIndex> distance_index;
    if (!distance_name.empty()) {
        distance_index = vg::io::VPKG::load_one<MinimumDistanceIndex>(distance_name);
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
        bool found = false;
        if (distance_index == nullptr) {
            auto iter = true_path_positions.find(aln.name());
            if (iter != true_path_positions.end()) {
                alignment_set_distance_to_correct(aln, iter->second);
                found = true;
            }
        } else {
            auto iter = true_graph_positions.find(aln.name());
            if (iter != true_graph_positions.end() && aln.path().mapping_size() > 0) {
                std::vector<MappingRun> read_mappings = base_mappings(aln);
                int64_t distance = std::numeric_limits<int64_t>::max();
                auto read_iter = read_mappings.begin();
                auto truth_iter = iter->second.begin();
                // Break the read into maximal intervals such that each interval corresponds
                // to a gapless alignment between the read and a single node both in the true
                // alignment and the candidate alignment. Compute the distance for each
                // interval and use the minimum distance over all intervals.
                while (read_iter != read_mappings.end() && truth_iter != iter->second.end()) {
                    size_t start = std::max(read_iter->read_offset, truth_iter->read_offset);
                    size_t limit = std::min(read_iter->limit(), truth_iter->limit());
                    if (start < limit) {
                        pos_t read_pos = read_iter->pos_at(start);
                        pos_t truth_pos = truth_iter->pos_at(start);
                        int64_t forward = distance_index->min_distance(read_pos, truth_pos);
                        if (forward != -1) {
                            distance = std::min(forward, distance);
                        }
                        int64_t reverse = distance_index->min_distance(truth_pos, read_pos);
                        if (reverse != -1) {
                            distance = std::min(reverse, distance);
                        }
                    }
                    if (read_iter->limit() <= limit) {
                        ++read_iter;
                    }
                    if (truth_iter->limit() <= limit) {
                        ++truth_iter;
                    }
                }
                Position result;
                result.set_name("graph");
                result.set_offset(distance);
                *aln.mutable_to_correct() = result;
                found = true;
            }
        }
        if (found && range != -1) {
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
                    mapq_count_by_thread.at(omp_get_thread_num()).resize(mapq+1, 0);
                    correct_count_by_mapq_by_thread.at(omp_get_thread_num()).resize(mapq+1, 0);
                }

                read_count_by_thread.at(omp_get_thread_num()) += 1;
                mapq_count_by_thread.at(omp_get_thread_num()).at(mapq) += 1;
                if (correctly_mapped) {
                    correct_count_by_mapq_by_thread.at(omp_get_thread_num()).at(mapq) += 1;
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
                    mapq_count.resize(mq+1, 0);
                    correct_count_by_mapq.resize(mq+1, 0);
                }

                mapq_count.at(mq) += mapq_count_by_thread.at(i).at(mq);
                correct_count_by_mapq.at(mq) += correct_count_by_mapq_by_thread.at(i).at(mq);
            }
        }
        size_t accumulated_count = 0;
        size_t accumulated_correct_count = 0;
        float mapping_goodness_score = 0.0;
        for (int i = mapq_count.size()-1 ; i >= 0 ; i--) {
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
