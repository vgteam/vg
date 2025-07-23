// gamcompare_main.cpp: defines a GAM to GAM annotation function

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>
#include <iomanip>

#include "subcommand.hpp"

#include "../alignment.hpp"
#include "../annotation.hpp"
#include "../snarl_distance_index.hpp"
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
         << "  -d, --distance-index FILE  use distances from this distance index" << endl
         << "                             instead of path position annotations" << endl
         << "  -r, --range N              distance within which to consider reads correct" << endl
         << "  -n, --rename Q=T           treat query contig Q as truth contig T (may repeat)" << endl
         << "  -I, --ignore T             ignore the given truth contig name (may repeat)" << endl
         << "  -o, --output-gam FILE      output GAM to FILE instead of standard output" << endl
         << "  -T, --tsv                  output TSV (correct, mq, aligner, read)" << endl
         << "                             compatible with plot-qq.R to standard output" << endl
         << "  -a, --aligner STR          aligner name for TSV output [\"vg\"]" << endl
         << "  -s, --score-alignment      get alignment correctness score (higher is better)" << endl
         << "  -t, --threads N            number of threads to use" << endl
         << "  -h, --help                 print this help message to stderr and exit" << endl;
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
    string output_gam;
    bool output_tsv = false;
    string aligner_name = "vg";
    bool score_alignment = false;
    string distance_name;
    // Map from query contigs to corresponding truth contigs
    std::unordered_map<string, string> renames;
    // Keep a set of ignored truth contigs
    std::unordered_set<std::string> ignores;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"distance-index", required_argument, 0, 'd'},
            {"range", required_argument, 0, 'r'},
            {"rename", required_argument, 0, 'n'},
            {"ignore", required_argument, 0, 'I'},
            {"output-gam", required_argument, 0, 'o'},
            {"tsv", no_argument, 0, 'T'},
            {"aligner", required_argument, 0, 'a'},
            {"score-alignment", no_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "h?d:r:I:n:o:Ta:st:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'r':
            range = parse<int>(optarg);
            break;
            
        case 'n':
            {
                // Parse the rename old=new
                string key_value(optarg);
                auto found = key_value.find('=');
                if (found == string::npos || found == 0 || found + 1 == key_value.size()) {
                    cerr << "error:[vg gamcompare] could not parse rename " << key_value << endl;
                    exit(1);
                }
                // Parse out the two parts
                string query_contig = key_value.substr(0, found);
                string truth_contig = key_value.substr(found + 1);
                // Add the name mapping
                renames.emplace(query_contig, truth_contig);
            }
            break;

        case 'I':
            ignores.insert(optarg);
            break;

        case 'd':
            distance_name = optarg;
            break;

        case 'o':
            output_gam = optarg;
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
    vg::string_hash_map<string, map<string, vector<pair<size_t, bool>>>> true_path_positions;
    function<void(Alignment&)> record_path_positions = [&true_path_positions,&ignores](Alignment& aln) {
        if (aln.refpos_size() > 0) {
            std::map<std::string, std::vector<std::pair<size_t, bool>>> val = alignment_refpos_to_path_offsets(aln);

            // TODO: Is it faster to poll all the contigs against the ignores
            // list and drop them as we go, or look up and remove each ignored
            // contig?
            auto it = val.begin();
            while(it != val.end()) {
                // See if each contig we have a position on is ignored.
                if (ignores.count(it->first)) {
                    // Drop this contig
                    it = val.erase(it);
                } else {
                    // Keep this contig
                    ++it;
                }
            }

            if (!val.empty()) {
                #pragma omp critical (truth_table)
                true_path_positions[aln.name()] = val;
            }
        }
    };

    // True graph positions. For each alignment name, we find the maximal read intervals that correspond
    // to a gapless alignment between the read and a single node.
    vg::string_hash_map<string, std::vector<MappingRun>> true_graph_positions;
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

    // Count eligible reads that actually have positions that could be got.
    size_t eligible_reads = distance_name.empty() ? true_path_positions.size() : true_graph_positions.size();

    // Load the distance index.
    unique_ptr<SnarlDistanceIndex> distance_index;
    if (!distance_name.empty()) {
        distance_index = vg::io::VPKG::load_one<SnarlDistanceIndex>(distance_name);
    }

    // We have a buffered emitter for annotated alignments, if we're not outputting text.
    // Start out with this empty so we output nowhere.
    std::unique_ptr<vg::io::ProtobufEmitter<Alignment>> emitter;
    std::ofstream output_gam_stream;
    if (!output_gam.empty()) {
        // Output to specified location
        output_gam_stream.open(output_gam, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
        if (output_gam_stream.fail() || !output_gam_stream.is_open()) {
            cerr << "error[vg gamcompare]: Cannot output to " << output_gam << endl;
            exit(1);
        }
        emitter = std::unique_ptr<vg::io::ProtobufEmitter<Alignment>>(new vg::io::ProtobufEmitter<Alignment>(output_gam_stream));
    } else if (!output_tsv) {
        // Output to standard output.
        emitter = std::unique_ptr<vg::io::ProtobufEmitter<Alignment>>(new vg::io::ProtobufEmitter<Alignment>(cout));
    }
    
    // We have an ordinary buffer we use for text output
    vector<Alignment> text_buffer;
    
    // We have an output function to dump all the reads in the text buffer in TSV
    auto flush_text_buffer = [&text_buffer,&aligner_name]() {
        // We print exactly one header line.
        static bool header_printed = false;
        // Output TSV to standard out in the format plot-qq.R needs.
        if (!header_printed) {
            // It needs a header
            cout << "correct\tmq\taligner\tread\teligible" << endl;
            header_printed = true;
        }
        
        for (auto& aln : text_buffer) {
            // Dump each alignment
            cout << (aln.correctly_mapped() ? "1" : "0") << "\t";
            cout << aln.mapping_quality() << "\t";
            cout << aligner_name << "\t";
            cout << aln.name() << "\t";
            cout << (has_annotation(aln, "no_truth") ? "0" : "1") << endl;
        }
        text_buffer.clear();
    };
   
    // We want to count correct reads
    vector<size_t> correct_counts(vg::get_thread_count(), 0);

    //Get stats for calculating the score
    vector<size_t> read_count_by_thread (vg::get_thread_count(), 0);
    vector<vector<size_t>> mapq_count_by_thread (vg::get_thread_count());
    vector<vector<size_t>> correct_count_by_mapq_by_thread(vg::get_thread_count());
    for (size_t i = 0 ; i < vg::get_thread_count() ; i++) {
        mapq_count_by_thread[i].resize(61, 0);
        correct_count_by_mapq_by_thread[i].resize(61,0);
    }
   
    // This function annotates every read with distance and correctness, and batch-outputs them.
    function<void(Alignment&)> annotate_test = [&](Alignment& aln) {
        bool found = false;
        if (distance_name.empty()) {
            //If the distance index isn't used
            auto iter = true_path_positions.find(aln.name());
            if (iter != true_path_positions.end()) {
                alignment_set_distance_to_correct(aln, iter->second, &renames);
                found = true;
            }
        } else {
            //If the distance index gets used
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
                        size_t forward = minimum_distance(*distance_index, read_pos, truth_pos);
                        if (forward != std::numeric_limits<size_t>::max()) {
                            distance = std::min((int64_t)forward, distance);
                        }
                        size_t reverse = minimum_distance(*distance_index, truth_pos, read_pos);
                        if (reverse != std::numeric_limits<size_t>::max()) {
                            distance = std::min((int64_t)reverse, distance);
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
            // And make sure we say it was possible to get
            clear_annotation(aln, "no_truth");
            
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
        } else if (range != -1) {
            // We are flagging reads correct/incorrect, but this read has no truth position.
            // Remember that it was impossible to get.
            set_annotation(aln, "no_truth", true);
        }
#pragma omp critical
        {
            if (output_tsv) {
                if (emitter) {
                    // Copy the alignment since we need it twice
                    text_buffer.emplace_back(aln);
                } else {
                    text_buffer.emplace_back(std::move(aln));
                }
                if (text_buffer.size() > 1000) {
                    flush_text_buffer();
                }
            }
            if (emitter) {
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
        
        cerr << total_correct << " reads correct, " << eligible_reads << " reads eligible";
        if (eligible_reads > 0 && eligible_reads >= total_correct) {
            std::ios state(nullptr);
            state.copyfmt(cerr);
            cerr << ", " << std::fixed << std::setprecision(2) << (double)total_correct / eligible_reads * 100 << "% accuracy";
            cerr.copyfmt(state);
        }
        cerr << endl;
    }

    if (score_alignment) {
        //Get a goodness score of the alignment that takes into account correctness and mapq calibration
        size_t total_reads = 0;
        vector<size_t> mapq_count (61, 0);
        vector<size_t> correct_count_by_mapq (61, 0);
        for (size_t i = 0 ; i < vg::get_thread_count() ; i++) {
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

    if (emitter) {
        // Make sure to get rid of the emitter before the file it might write to
        emitter.reset();
    }
    if (output_gam_stream.is_open()) {
        output_gam_stream.close();
    }
    
    return 0;
}

// Register subcommand
static Subcommand vg_gamcompare("gamcompare", "compare alignment positions", main_gamcompare);
