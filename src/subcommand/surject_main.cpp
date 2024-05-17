// surject_main.cpp: define the "vg surject" subcommand, which forces alignments into linear space

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include "subcommand.hpp"

#include <bdsg/hash_graph.hpp>
#include <bdsg/overlays/path_position_overlays.hpp>
#include <bdsg/overlays/overlay_helper.hpp>

#include "../vg.hpp"
#include "../xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../utility.hpp"
#include "../surjector.hpp"
#include "../hts_alignment_emitter.hpp"
#include "../multipath_alignment_emitter.hpp"
#include "../crash.hpp"
#include "../watchdog.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_surject(char** argv) {
    cerr << "usage: " << argv[0] << " surject [options] <aln.gam> >[proj.cram]" << endl
         << "Transforms alignments to be relative to particular paths." << endl
         << endl
         << "options:" << endl
         << "  -x, --xg-name FILE       use this graph or xg index (required)" << endl
         << "  -t, --threads N          number of threads to use" << endl
         << "  -p, --into-path NAME     surject into this path or its subpaths (many allowed, default: reference, then non-alt generic)" << endl
         << "  -F, --into-paths FILE    surject into path names listed in HTSlib sequence dictionary or path list FILE" << endl
         << "  -n, --into-sample NAME   surject into pahts from this sample (multiple allowed)" << endl
         << "  -i, --interleaved        GAM is interleaved paired-ended, so when outputting HTS formats, pair reads" << endl
         << "  -M, --multimap           include secondary alignments to all overlapping paths instead of just primary" << endl
         << "  -G, --gaf-input          input file is GAF instead of GAM" << endl
         << "  -m, --gamp-input         input file is GAMP instead of GAM" << endl
         << "  -c, --cram-output        write CRAM to stdout" << endl
         << "  -b, --bam-output         write BAM to stdout" << endl
         << "  -s, --sam-output         write SAM to stdout" << endl
         << "  -l, --subpath-local      let the multipath mapping surjection produce local (rather than global) alignments" << endl
         << "  -P, --prune-low-cplx     prune short and low complexity anchors during realignment" << endl
         << "  -a, --max-anchors N      use no more than N anchors per target path (default: 200)" << endl
         << "  -S, --spliced            interpret long deletions against paths as spliced alignments" << endl
         << "  -A, --qual-adj           adjust scoring for base qualities, if they are available" << endl
         << "  -N, --sample NAME        set this sample name for all reads" << endl
         << "  -R, --read-group NAME    set this read group for all reads" << endl
         << "  -f, --max-frag-len N     reads with fragment lengths greater than N will not be marked properly paired in SAM/BAM/CRAM" << endl
         << "  -L, --list-all-paths     annotate SAM records with a list of all attempted re-alignments to paths in SS tag" << endl
         << "  -C, --compression N      level for compression [0-9]" << endl
         << "  -V, --no-validate        skip checking whether alignments plausibly are against the provided graph" << endl
         << "  -w, --watchdog-timeout N warn when reads take more than the given number of seconds to surject" << endl;
}

/// If the given alignment doesn't make sense against the given graph (i.e.
/// doesn't agree with the nodes in the graph), print a message and stop the
/// program. Is thread-safe.
static void ensure_alignment_is_for_graph(const Alignment& aln, const HandleGraph& graph) {
    AlignmentValidity validity = alignment_is_valid(aln, &graph);
    if (!validity) {
        #pragma omp critical (cerr)
        {
            std::cerr << "error:[vg surject] Alignment " << aln.name() << " cannot be interpreted against this graph: " << validity.message << std::endl;
            std::cerr << "Make sure that you are using the same graph that the reads were mapped to!" << std::endl;
        }
        exit(1);
    }
}

int main_surject(int argc, char** argv) {
    
    if (argc == 2) {
        help_surject(argv);
        return 1;
    }
    
    string xg_name;
    vector<string> path_names;
    vector<string> path_sample_names;
    string path_file;
    string output_format = "GAM";
    string input_format = "GAM";
    bool spliced = false;
    bool interleaved = false;
    string sample_name;
    string read_group;
    int32_t max_frag_len = 0;
    int compress_level = 9;
    int min_splice_length = 20;
    size_t watchdog_timeout = 10;
    bool subpath_global = true; // force full length alignments in mpmap resolution
    bool qual_adj = false;
    bool prune_anchors = false;
    size_t max_anchors = 200;
    bool annotate_with_all_path_scores = false;
    bool multimap = false;
    bool validate = true;

    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"xg-name", required_argument, 0, 'x'},
            {"threads", required_argument, 0, 't'},
            {"into-path", required_argument, 0, 'p'},
            {"into-paths", required_argument, 0, 'F'},
            {"into-sample", required_argument, 0, 'n'},
            {"ref-paths", required_argument, 0, 'F'}, // Now an alias for --into-paths
            {"subpath-local", required_argument, 0, 'l'},
            {"interleaved", no_argument, 0, 'i'},
            {"multimap", no_argument, 0, 'M'},
            {"gaf-input", no_argument, 0, 'G'},
            {"gamp-input", no_argument, 0, 'm'},
            {"cram-output", no_argument, 0, 'c'},
            {"bam-output", no_argument, 0, 'b'},
            {"sam-output", no_argument, 0, 's'},
            {"spliced", no_argument, 0, 'S'},
            {"prune-low-cplx", no_argument, 0, 'P'},
            {"max-anchors", required_argument, 0, 'a'},
            {"qual-adj", no_argument, 0, 'A'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"max-frag-len", required_argument, 0, 'f'},
            {"list-all-paths", no_argument, 0, 'L'},
            {"compress", required_argument, 0, 'C'},
            {"no-validate", required_argument, 0, 'V'},
            {"watchdog-timeout", required_argument, 0, 'w'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:p:F:n:liGmcbsN:R:f:C:t:SPa:ALMVw:",
                long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'x':
            xg_name = optarg;
            break;

        case 'p':
            path_names.push_back(optarg);
            break;

        case 'F':
            path_file = optarg;
            break;

        case 'n':
            path_sample_names.push_back(optarg);
            break;            
        
        case 'l':
            subpath_global = false;
            break;

        case 'i':
            interleaved = true;
            break;
                
        case 'M':
            multimap = true;
            break;

        case 'G':
            input_format = "GAF";
            break;
                
        case 'm':
            input_format = "GAMP";
            break;
            
        case 'c':
            output_format = "CRAM";
            break;

        case 'b':
            output_format = "BAM";
            break;

        case 's':
            compress_level = -1;
            output_format = "SAM";
            break;
                
        case 'S':
            spliced = true;
            break;
                
        case 'P':
            prune_anchors = true;
            break;
            
        case 'a':
            max_anchors = parse<size_t>(optarg);
            break;
            
        case 'A':
            qual_adj = true;
            break;

        case 'N':
            sample_name = optarg;
            break;
            
        case 'R':
            read_group = optarg;
            break;
            
        case 'f':
            max_frag_len = parse<int32_t>(optarg);
            break;

        case 'C':
            compress_level = parse<int>(optarg);
            break;
            
        case 'V':
            validate = false;
            break;
            
        case 'w':
            watchdog_timeout = parse<size_t>(optarg);
            break;
            
        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;
                
        case 'L':
            annotate_with_all_path_scores = true;
            break;

        case 'h':
        case '?':
            help_surject(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    // Create a preprocessor to apply read group and sample name overrides in place
    auto set_metadata = [&](Alignment& update) {
        if (!sample_name.empty()) {
            update.set_sample_name(sample_name);
        }
        if (!read_group.empty()) {
            update.set_read_group(read_group);
        }
    };

    string file_name = get_input_file_name(optind, argc, argv);
    
    PathPositionHandleGraph* xgidx = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    // If we add an overlay for path position queries, use one optimized for
    // use with reference paths.
    bdsg::ReferencePathOverlayHelper overlay_helper;
    if (!xg_name.empty()) {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        xgidx = overlay_helper.apply(path_handle_graph.get());
    } else {
        // We need an XG index for the rest of the algorithm
        cerr << "error[vg surject] XG index (-x) is required for surjection" << endl;
        exit(1);
    }
    
    // Get the paths to surject into and their length information, either from
    // the given file, or from the provided list, or from sniffing the graph.
    for (const string& sample_name : path_sample_names) {
        path_handle_graph->for_each_path_of_sample(sample_name, [&](path_handle_t path_handle) {
            path_names.push_back(path_handle_graph->get_path_name(path_handle));
        });
    }
    vector<tuple<path_handle_t, size_t, size_t>> sequence_dictionary = get_sequence_dictionary(path_file, path_names, *xgidx);
    // Clear out path_names so we don't accidentally use it
    path_names.clear();

    // Convert to a set for membership testing 
    unordered_set<path_handle_t> paths;
    paths.reserve(sequence_dictionary.size());
    for (auto& entry : sequence_dictionary) {
        paths.insert(get<0>(entry));
    }
    
    // Make a single thread-safe Surjector.
    Surjector surjector(xgidx);
    surjector.adjust_alignments_for_base_quality = qual_adj;
    surjector.prune_suspicious_anchors = prune_anchors;
    surjector.max_anchors = max_anchors;
    if (spliced) {
        surjector.min_splice_length = min_splice_length;
        // we have to bump this up to be sure to align most splice junctions
        surjector.max_subgraph_bases = 16 * 1024 * 1024;
    }
    else {
        surjector.min_splice_length = numeric_limits<int64_t>::max();
    }
    surjector.annotate_with_all_path_scores = annotate_with_all_path_scores;
    
    // Count our threads
    int thread_count = vg::get_thread_count();
    
    // Prepare the watchdog
    unique_ptr<Watchdog> watchdog(new Watchdog(thread_count, chrono::seconds(watchdog_timeout)));
    
    if (input_format == "GAM" || input_format == "GAF") {
        
        // Give helpful warning if someone tries to surject an un-surjectable GAF
        auto check_gaf_aln = [&](const Alignment& src) {
            if (src.has_path() && src.sequence().empty()) {
#pragma omp critical
                {
                    cerr << "error:[surject] Read " << src.name() << " is aligned but does not have a sequence and therefore cannot be surjected. Was it derived from a GAF without a base-level alignment? Or a GAF with a CIGAR string in the 'cg' tag (which does not provide enough information to reconstruct the sequence)?" << endl;
                    exit(1);
                }
            }
        };
        
        // Set up output to an emitter that will handle serialization.
        // It should process output raw, without any surjection, and it should
        // respect our parameter for whether to think with splicing.
        unique_ptr<AlignmentEmitter> alignment_emitter = get_alignment_emitter("-", 
            output_format, sequence_dictionary, thread_count, xgidx,
            ALIGNMENT_EMITTER_FLAG_HTS_RAW | (spliced * ALIGNMENT_EMITTER_FLAG_HTS_SPLICED));

        if (interleaved) {
            // GAM input is paired, and for HTS output reads need to know their pair partners' mapping locations.
            // TODO: We don't preserve order relationships (like primary/secondary) beyond the interleaving.
            function<void(Alignment&, Alignment&)> lambda = [&](Alignment& src1, Alignment& src2) {
                try {
                    set_crash_context(src1.name() + ", " + src2.name());
                    size_t thread_num = omp_get_thread_num();
                    if (watchdog) {
                        watchdog->check_in(thread_num, src1.name() + ", " + src2.name());
                    }
                    // Make sure that the alignments are actually paired with each other
                    // (proper fragment_prev/fragment_next). We want to catch people giving us
                    // un-interleaved GAMs as interleaved.
                    // TODO: Integrate into for_each_interleaved_pair_parallel when running on Alignments.
                    if (src1.has_fragment_next()) {
                        // Alignment 1 comes first in fragment
                        if (src1.fragment_next().name() != src2.name() ||
                            !src2.has_fragment_prev() ||
                            src2.fragment_prev().name() != src1.name()) {
                            
#pragma omp critical (cerr)
                            cerr << "[vg surject] error: alignments " << src1.name()
                            << " and " << src2.name() << " are adjacent but not paired" << endl;
                            
                            exit(1);
                            
                        }
                    } else if (src2.has_fragment_next()) {
                        // Alignment 2 comes first in fragment
                        if (src2.fragment_next().name() != src1.name() ||
                            !src1.has_fragment_prev() ||
                            src1.fragment_prev().name() != src2.name()) {
                            
#pragma omp critical (cerr)
                            cerr << "[vg surject] error: alignments " << src1.name()
                            << " and " << src2.name() << " are adjacent but not paired" << endl;
                            
                            exit(1);
                            
                        }
                    } else {
                        // Alignments aren't paired up at all
#pragma omp critical (cerr)
                        cerr << "[vg surject] error: alignments " << src1.name()
                        << " and " << src2.name() << " are adjacent but not paired" << endl;
                        
                        exit(1);
                    }
                    
                    if (validate) {
                        ensure_alignment_is_for_graph(src1, *xgidx);
                        ensure_alignment_is_for_graph(src2, *xgidx);
                    }
                    
                    // Preprocess read to set metadata before surjection
                    set_metadata(src1);
                    set_metadata(src2);
                    
                    // Surject and emit.
                    if (multimap) {
                        
                        auto surjected1 = surjector.multi_surject(src1, paths, subpath_global, spliced);
                        auto surjected2 = surjector.multi_surject(src2, paths, subpath_global, spliced);
                        
                        // we have to pair these up manually
                        unordered_map<pair<string, bool>, size_t> strand_idx1, strand_idx2;
                        for (size_t i = 0; i < surjected1.size(); ++i) {
                            const auto& pos = surjected1[i].refpos(0);
                            strand_idx1[make_pair(pos.name(), pos.is_reverse())] = i;
                        }
                        for (size_t i = 0; i < surjected2.size(); ++i) {
                            const auto& pos = surjected2[i].refpos(0);
                            strand_idx2[make_pair(pos.name(), pos.is_reverse())] = i;
                        }
                        
                        for (size_t i = 0; i < surjected1.size(); ++i) {
                            const auto& pos = surjected1[i].refpos(0);
                            auto it = strand_idx2.find(make_pair(pos.name(), !pos.is_reverse()));
                            if (it != strand_idx2.end()) {
                                // the alignments are paired on this strand
                                alignment_emitter->emit_pair(move(surjected1[i]), move(surjected2[it->second]), max_frag_len);
                            }
                            else {
                                // this strand's surjection is unpaired
                                alignment_emitter->emit_single(move(surjected1[i]));
                            }
                        }
                        for (size_t i = 0; i < surjected2.size(); ++i) {
                            const auto& pos = surjected2[i].refpos(0);
                            if (!strand_idx1.count(make_pair(pos.name(), !pos.is_reverse()))) {
                                // this strand's surjection is unpaired
                                alignment_emitter->emit_single(move(surjected2[i]));
                            }
                        }
                    }
                    else {
                        // FIXME: these aren't forced to be on the same path, which could be fucky
                        alignment_emitter->emit_pair(surjector.surject(src1, paths, subpath_global, spliced),
                                                     surjector.surject(src2, paths, subpath_global, spliced),
                                                     max_frag_len);
                    }
                    if (watchdog) {
                        watchdog->check_out(thread_num);
                    }
                    clear_crash_context();
                } catch (const std::exception& ex) {
                    report_exception(ex);
                }
            };
            if (input_format == "GAM") {
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each_interleaved_pair_parallel<Alignment>(in, lambda);
                });
            } else {
                auto gaf_checking_lambda = [&](Alignment& src1, Alignment& src2) {
                    check_gaf_aln(src1);
                    check_gaf_aln(src2);
                    return lambda(src1, src2);
                };
                vg::io::gaf_paired_interleaved_for_each_parallel(*xgidx, file_name, gaf_checking_lambda);
            }
        } else {
            // We can just surject each Alignment by itself.
            // TODO: We don't preserve order relationships (like primary/secondary).
            function<void(Alignment&)> lambda = [&](Alignment& src) {
                try {
                    set_crash_context(src.name());
                    size_t thread_num = omp_get_thread_num();
                    if (watchdog) {
                        watchdog->check_in(thread_num, src.name());
                    }
                    if (validate) {
                        ensure_alignment_is_for_graph(src, *xgidx);
                    }
                    
                    // Preprocess read to set metadata before surjection
                    set_metadata(src);
                    
                    // Surject and emit the single read.
                    if (multimap) {
                        alignment_emitter->emit_singles(surjector.multi_surject(src, paths, subpath_global, spliced));
                    }
                    else {
                        alignment_emitter->emit_single(surjector.surject(src, paths, subpath_global, spliced));
                    }
                    if (watchdog) {
                        watchdog->check_out(thread_num);
                    }
                    clear_crash_context();
                } catch (const std::exception& ex) {
                    report_exception(ex);
                }
            };
            if (input_format == "GAM") {
                get_input_file(file_name, [&](istream& in) {
                    vg::io::for_each_parallel<Alignment>(in,lambda);
                });
            } else {
                auto gaf_checking_lambda = [&](Alignment& src) {
                    check_gaf_aln(src);
                    return lambda(src);
                };
                vg::io::gaf_unpaired_for_each_parallel(*xgidx, file_name, gaf_checking_lambda);
            }
        }
    } else if (input_format == "GAMP") {
        // Working on multipath alignments. We need to set the emitter up ourselves.
        auto path_order_and_length = extract_path_metadata(sequence_dictionary, *xgidx).first;
        MultipathAlignmentEmitter mp_alignment_emitter("-", thread_count, output_format, xgidx, &path_order_and_length);
        mp_alignment_emitter.set_read_group(read_group);
        mp_alignment_emitter.set_sample_name(sample_name);
        mp_alignment_emitter.set_min_splice_length(spliced ? min_splice_length : numeric_limits<int64_t>::max());
        
        // TODO: largely repetitive with GAM
        get_input_file(file_name, [&](istream& in) {
            if (interleaved) {
                
                // GAMP input is paired, and for HTS output reads need to know their pair partners' mapping locations.
                // TODO: We don't preserve order relationships (like primary/secondary) beyond the interleaving.
                vg::io::for_each_interleaved_pair_parallel<MultipathAlignment>(in, [&](MultipathAlignment& src1, MultipathAlignment& src2) {
                    try {
                        set_crash_context(src1.name() + ", " + src2.name());
                        size_t thread_num = omp_get_thread_num();
                        if (watchdog) {
                            watchdog->check_in(thread_num, src1.name() + ", " + src2.name());
                        }
                    
                        // Make sure that the alignments are actually paired with each other
                        // (proper fragment_prev/fragment_next). We want to catch people giving us
                        // un-interleaved GAMs as interleaved.
                        // TODO: Integrate into for_each_interleaved_pair_parallel when running on Alignments.
                        if (src1.paired_read_name() != src2.name() || src2.paired_read_name() != src1.name()) {
                            
#pragma omp critical (cerr)
                            cerr << "[vg surject] error: alignments " << src1.name()
                            << " and " << src2.name() << " are adjacent but not paired" << endl;
                            
                            exit(1);
                            
                        }
                        else if (src1.paired_read_name().empty() || src2.paired_read_name().empty()) {
                            // Alignments aren't paired up at all
#pragma omp critical (cerr)
                            cerr << "[vg surject] error: alignments " << src1.name()
                            << " and " << src2.name() << " are adjacent but not paired" << endl;
                            
                            exit(1);
                        }
                        
                        // convert out of protobuf
                        multipath_alignment_t mp_src1, mp_src2;
                        from_proto_multipath_alignment(src1, mp_src1);
                        from_proto_multipath_alignment(src2, mp_src2);
                        
                        
                        vector<pair<tuple<string, bool, int64_t>, tuple<string, bool, int64_t>>> positions;
                        vector<pair<multipath_alignment_t, multipath_alignment_t>> surjected;
                        
                        vector<tuple<string, bool, int64_t>> positions_unpaired1, positions_unpaired2;
                        vector<multipath_alignment_t> surjected_unpaired1, surjected_unpaired2;
                        
                        // surject and record path positions
                        if (multimap) {
                            
                            // TODO: highly repetitive with the version above for Alignments
                            
                            vector<tuple<string, int64_t, bool>> positions1, positions2;
                            auto surjected1 = surjector.multi_surject(mp_src1, paths, positions1, subpath_global, spliced);
                            auto surjected2 = surjector.multi_surject(mp_src2, paths, positions2, subpath_global, spliced);
                            
                            // we have to pair these up manually
                            unordered_map<pair<string, bool>, size_t> strand_idx1, strand_idx2;
                            for (size_t i = 0; i < surjected1.size(); ++i) {
                                strand_idx1[make_pair(get<0>(positions1[i]), get<2>(positions1[i]))] = i;
                            }
                            for (size_t i = 0; i < surjected2.size(); ++i) {
                                strand_idx2[make_pair(get<0>(positions2[i]), get<2>(positions2[i]))] = i;
                            }
                                                    
                            for (size_t i = 0; i < surjected1.size(); ++i) {
                                auto it = strand_idx2.find(make_pair(get<0>(positions1[i]), !get<2>(positions1[i])));
                                if (it != strand_idx2.end()) {
                                    // the alignments are paired on this strand
                                    size_t j = it->second;
                                    surjected.emplace_back(move(surjected1[i]), move(surjected2[j]));
                                    
                                    // reorder the positions to deal with the mismatch in the interfaces
                                    positions.emplace_back();
                                    get<0>(positions.back().first) = get<0>(positions1[i]);
                                    get<1>(positions.back().first) = get<2>(positions1[i]);
                                    get<2>(positions.back().first) = get<1>(positions1[i]);
                                    get<0>(positions.back().second) = get<0>(positions2[j]);
                                    get<1>(positions.back().second) = get<2>(positions2[j]);
                                    get<2>(positions.back().second) = get<1>(positions2[j]);
                                }
                                else {
                                    // this strand's surjection is unpaired
                                    surjected_unpaired1.emplace_back(move(surjected1[i]));
                                    
                                    // reorder the position to deal with the mismatch in the interfaces
                                    positions_unpaired1.emplace_back();
                                    get<0>(positions_unpaired1.back()) = move(get<0>(positions1[i]));
                                    get<1>(positions_unpaired1.back()) = get<2>(positions1[i]);
                                    get<2>(positions_unpaired1.back()) = get<1>(positions1[i]);
                                }
                            }
                            for (size_t i = 0; i < surjected2.size(); ++i) {
                                if (!strand_idx1.count(make_pair(get<0>(positions2[i]), !get<2>(positions2[i])))) {
                                    // this strand's surjection is unpaired
                                    surjected_unpaired2.emplace_back(move(surjected2[i]));
                                    
                                    // reorder the position to deal with the mismatch in the interfaces
                                    positions_unpaired2.emplace_back();
                                    get<0>(positions_unpaired2.back()) = move(get<0>(positions2[i]));
                                    get<1>(positions_unpaired2.back()) = get<2>(positions2[i]);
                                    get<2>(positions_unpaired2.back()) = get<1>(positions2[i]);
                                }
                            }
                        }
                        else {
                            
                            // FIXME: these aren't required to be on the same path...
                            positions.emplace_back();
                            surjected.emplace_back(surjector.surject(mp_src1, paths, get<0>(positions.front().first),
                                                                     get<2>(positions.front().first), get<1>(positions.front().first),
                                                                     subpath_global, spliced),
                                                   surjector.surject(mp_src2, paths, get<0>(positions.front().second),
                                                                     get<2>(positions.front().second), get<1>(positions.front().second),
                                                                     subpath_global, spliced));
                        }
                                            
                        // write to output
                        vector<int64_t> tlen_limits(surjected.size(), max_frag_len);
                        mp_alignment_emitter.emit_pairs(src1.name(), src2.name(), move(surjected), &positions, &tlen_limits);
                        mp_alignment_emitter.emit_singles(src1.name(), move(surjected_unpaired1), &positions_unpaired1);
                        mp_alignment_emitter.emit_singles(src2.name(), move(surjected_unpaired2), &positions_unpaired2);
                        
                        if (watchdog) {
                            watchdog->check_out(thread_num);
                        }
                        clear_crash_context();
                    } catch (const std::exception& ex) {
                        report_exception(ex);
                    }
                });
            } else {
                // TODO: We don't preserve order relationships (like primary/secondary).
                vg::io::for_each_parallel<MultipathAlignment>(in, [&](MultipathAlignment& src) {
                    try {
                        set_crash_context(src.name());
                        size_t thread_num = omp_get_thread_num();
                        if (watchdog) {
                            watchdog->check_in(thread_num, src.name());
                        }

                        multipath_alignment_t mp_src;
                        from_proto_multipath_alignment(src, mp_src);
                        
                        // surject and record path positions
                        vector<tuple<string, bool, int64_t>> positions;
                        vector<multipath_alignment_t> surjected;
                        
                        if (multimap) {
                            
                            vector<tuple<string, int64_t, bool>> multi_positions;
                            surjected = surjector.multi_surject(mp_src, paths, multi_positions, subpath_global, spliced);
                            
                            // positions are in different orders in these two interfaces
                            for (auto& position : multi_positions) {
                                positions.emplace_back(move(get<0>(position)), get<2>(position), get<1>(position));
                            }
                        }
                        else {
                            positions.emplace_back();
                            surjected.emplace_back(surjector.surject(mp_src, paths, get<0>(positions.front()),
                                                                     get<2>(positions.front()), get<1>(positions.front()),
                                                                     subpath_global, spliced));
                        }
                        
                        // write to output
                        mp_alignment_emitter.emit_singles(src.name(), move(surjected), &positions);
                    
                        if (watchdog) {
                            watchdog->check_out(thread_num);
                        }
                        clear_crash_context();
                    } catch (const std::exception& ex) {
                        report_exception(ex);
                    }
                });
            }
        });
    } else {
        cerr << "[vg surject] Unimplemented input format " << input_format << endl;
        exit(1);
    }
    
    cout.flush();
    
    return 0;
}

// Register subcommand
static Subcommand vg_surject("surject", "map alignments onto specific paths", main_surject);
