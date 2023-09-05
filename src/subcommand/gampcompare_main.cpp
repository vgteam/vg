// gampcompare_main.cpp: defines a GAMP to GAM annotation function

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include <subcommand.hpp>

#include "../algorithms/alignment_path_offsets.hpp"
#include "../multipath_alignment.hpp"
#include "../alignment.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_gampcompare(char** argv) {
    cerr << "usage: " << argv[0] << " gampcompare [options] alngraph.xg aln.gamp truth.gam > output.tsv" << endl
         << endl
         << "options:" << endl
         << "    -G, --gam                alignments are in GAM format rather than GAMP" << endl
         << "    -r, --range N            distance within which to consider reads correct [100]" << endl
         << "    -a, --aligner STR        aligner name for TSV output [\"vg\"]" << endl
         << "    -d, --distance           report minimum distance along a path rather than correctness" << endl
         << "    -t, --threads N          number of threads to use [1]" << endl;
}

int main_gampcompare(int argc, char** argv) {

    if (argc == 2) {
        help_gampcompare(argv);
        exit(1);
    }

    int threads = 1;
    int64_t range = 100;
    string aligner_name = "vg";
    int buffer_size = 10000;
    bool gam_input = false;
    bool report_distance = false;

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"range", required_argument, 0, 'r'},
            {"gam", no_argument, 0, 'G'},
            {"aligner", required_argument, 0, 'a'},
            {"distance", required_argument, 0, 'd'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hr:a:t:Gd",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1) break;

        switch (c)
        {

        case 'r':
            range = parse<int>(optarg);
            break;
            
        case 'a':
            aligner_name = optarg;
            break;
                
        case 'd':
            report_distance = true;
            break;

        case 't':
            threads = parse<int>(optarg);
            omp_set_num_threads(threads);
            break;
                
        case 'G':
            gam_input = true;
            break;

        case 'h':
        case '?':
            help_gampcompare(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    // We need to read the second argument first, so we can't use get_input_file with its free error checking.
    string graph_file_name = get_input_file_name(optind, argc, argv);
    string test_file_name = get_input_file_name(optind, argc, argv);
    string truth_file_name = get_input_file_name(optind, argc, argv);

    if ((truth_file_name == "-") + (test_file_name == "-") + (graph_file_name == "-") > 1) {
        cerr << "error[vg gampcompare]: Standard input can only be used for one input file" << endl;
        exit(1);
    }
    
    // Load the graph we mapped to
    unique_ptr<PathHandleGraph> path_handle_graph;
    if (graph_file_name == "-") {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(std::cin);
    }
    else {
        ifstream graph_stream(graph_file_name);
        if (!graph_stream) {
            cerr << "error:[vg mpmap] Cannot open graph file " << graph_file_name << endl;
            exit(1);
        }
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(graph_stream);
    }
    
    bdsg::PathPositionOverlayHelper overlay_helper;
    PathPositionHandleGraph* path_position_handle_graph = overlay_helper.apply(path_handle_graph.get());
    
    // We will collect all the truth positions
    string_hash_map<string, map<string ,vector<pair<size_t, bool> > > > true_positions;
    function<void(Alignment&)> record_truth = [&true_positions](Alignment& aln) {
        auto val = alignment_refpos_to_path_offsets(aln);
#pragma omp critical (truth_table)
        true_positions[move(*aln.mutable_name())] = move(val);
    };
    
    if (truth_file_name == "-") {
        // Read truth fropm standard input, if it looks good.
        if (!std::cin) {
            cerr << "error[vg gampcompare]: Unable to read standard input when looking for true reads" << endl;
            exit(1);
        }
        vg::io::for_each_parallel(std::cin, record_truth);
    }
    else {
        // Read truth from this file, if it looks good.
        ifstream truth_file_in(truth_file_name);
        if (!truth_file_in) {
            cerr << "error[vg gampcompare]: Unable to read " << truth_file_name << " when looking for true reads" << endl;
            exit(1);
        }
        vg::io::for_each_parallel(truth_file_in, record_truth);
    }
    
    // A buffer we use for the TSV output
    vector<vector<tuple<int64_t, bool, int64_t, int64_t, string>>> buffers(get_thread_count());
    
    // We have an output function to dump all the reads in the text buffer in TSV
    auto flush_buffer = [&](vector<tuple<int64_t, bool, int64_t, int64_t, string>>& buffer) {
        // We print exactly one header line.
        static bool header_printed = false;
        // Output TSV to standard out in the format plot-qq.R needs.
        if (!header_printed) {
            // It needs a header
            if (report_distance) {
                cout << "distance";
            }
            else {
                cout << "correct";
            }
            cout << "\tmapped\tmq\tgroupmq\taligner\tread" << endl;
            header_printed = true;
        }
        for (auto& result : buffer) {
            // Dump each alignment
            if (report_distance) {
                cout << get<0>(result);
            }
            else {
                cout << (get<0>(result) <= range);
            }
            cout << '\t' << get<1>(result) << '\t' << get<2>(result) << '\t' << get<3>(result) << '\t' << aligner_name << '\t' << get<4>(result) << endl;
        }
        buffer.clear();
    };
   
    // We want to count correct reads
    vector<size_t> correct_counts(get_thread_count(), 0);
   
    // This function annotates every read with distance and correctness, and batch-outputs them.
    function<void(MultipathAlignment&)> evaluate_correctness = [&](MultipathAlignment& proto_mp_aln) {
        
        // check the multipath mapping for correctness
        int64_t abs_dist = numeric_limits<int64_t>::max();
        auto f = true_positions.find(proto_mp_aln.name());
        if (f != true_positions.end()) {
            
            multipath_alignment_t mp_aln;
            from_proto_multipath_alignment(proto_mp_aln, mp_aln);
            
            auto& true_positions = f->second;
            auto mapped_positions = algorithms::multipath_alignment_path_offsets(*path_position_handle_graph,
                                                                                 mp_aln);
            for (auto it = true_positions.begin(); it != true_positions.end(); ++it) {
                // TODO: it really should be possible to do this with only path handles instead of names
                auto path_handle = path_position_handle_graph->get_path_handle(it->first);
                if (mapped_positions.count(path_handle)) {
                    // the true and mapped positions share this path
                    auto& path_true_positions = it->second;
                    auto& path_mapped_positions = mapped_positions[path_handle];
                    // check all pairs of positions
                    for (size_t i = 0; i < path_true_positions.size(); ++i) {
                        for (size_t j = 0; j < path_mapped_positions.size(); ++j) {
                            if (path_true_positions[i].second == path_mapped_positions[j].second) {
                                // there is a pair of positions on the same strand of the same path
                                abs_dist = min<int64_t>(abs_dist,
                                                        abs<int64_t>(path_true_positions[i].first - path_mapped_positions[j].first));
                            }
                        }
                    }
                }
            }
        }
        
        if (abs_dist <= range) {
            correct_counts[omp_get_thread_num()]++;
        }
        
        // group mapq defaults to regular mapq
        int64_t group_mapq = proto_mp_aln.mapping_quality();
        if (has_annotation(proto_mp_aln, "group_mapq")) {
            group_mapq = get_annotation<double>(proto_mp_aln, "group_mapq");
        }
        
        // put the result on the IO buffer
        auto& buffer = buffers[omp_get_thread_num()];
        buffer.emplace_back(abs_dist, proto_mp_aln.subpath_size() > 0, proto_mp_aln.mapping_quality(), group_mapq, move(*proto_mp_aln.mutable_name()));
        if (buffer.size() > buffer_size) {
#pragma omp critical
            flush_buffer(buffer);
        }
    };
    
    function<void(Alignment&)> evaluate_gam_correctness = [&](Alignment& aln) {
        // TODO: kinda ugly, but whatever
        multipath_alignment_t mp_aln;
        to_multipath_alignment(aln, mp_aln);
        MultipathAlignment proto_mp_aln;
        to_proto_multipath_alignment(mp_aln, proto_mp_aln);
        // we need the name to survive transit in multipath_alignment_t
        proto_mp_aln.set_name(aln.name());
        // now use the same evaluation function
        evaluate_correctness(proto_mp_aln);
    };

    if (test_file_name == "-") {
        if (!std::cin) {
            cerr << "error[vg gampcompare]: Unable to read standard input when looking for mapped reads" << endl;
            exit(1);
        }
        if (gam_input) {
            vg::io::for_each_parallel(std::cin, evaluate_gam_correctness);
        }
        else {
            vg::io::for_each_parallel(std::cin, evaluate_correctness);
        }
    } else {
        ifstream test_file_in(test_file_name);
        if (!test_file_in) {
            cerr << "error[vg gampcompare]: Unable to read " << test_file_name << " when looking for mapped reads" << endl;
            exit(1);
        }
        if (gam_input) {
            vg::io::for_each_parallel(test_file_in, evaluate_gam_correctness);
        }
        else {
            vg::io::for_each_parallel(test_file_in, evaluate_correctness);
        }
    }

    // Empty out whatever's in the buffers at the end.
    for (auto& buffer : buffers) {
        flush_buffer(buffer);
    }
    
    // We are flagging reads correct/incorrect. So report the total correct.
    size_t total_correct = 0;
    for (auto& count : correct_counts) {
        total_correct += count;
    }
    
    cerr << total_correct << " reads correct" << endl;
    
    return 0;
}

// Register subcommand
static Subcommand vg_gampcompare("gampcompare", "compare multipath alignment positions", main_gampcompare);
