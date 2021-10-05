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


using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_surject(char** argv) {
    cerr << "usage: " << argv[0] << " surject [options] <aln.gam> >[proj.cram]" << endl
         << "Transforms alignments to be relative to particular paths." << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE      use this graph or xg index (required)" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << "    -p, --into-path NAME    surject into this path (many allowed, default: all in xg)" << endl
         << "    -F, --into-paths FILE   surject into nonoverlapping path names listed in FILE (one per line)" << endl
         << "    --ref-paths FILE        ordered list of paths in the graph, one per line or HTSlib .dict, for HTSLib @SQ headers" << endl    
         << "    -l, --subpath-local     let the multipath mapping surjection produce local (rather than global) alignments" << endl
         << "    -i, --interleaved       GAM is interleaved paired-ended, so when outputting HTS formats, pair reads" << endl
         << "    -G, --gaf-input         input file is GAF instead of GAM" << endl
         << "    -m, --gamp-input        input file is GAMP instead of GAM" << endl
         << "    -c, --cram-output       write CRAM to stdout" << endl
         << "    -b, --bam-output        write BAM to stdout" << endl
         << "    -s, --sam-output        write SAM to stdout" << endl
         << "    -P, --prune-low-cplx    prune short and low complexity anchors during realignment" << endl
         << "    -S, --spliced           interpret long deletions against paths as spliced alignments" << endl
         << "    -A, --qual-adj          adjust scoring for base qualities, if they are available" << endl
         << "    -N, --sample NAME       set this sample name for all reads" << endl
         << "    -R, --read-group NAME   set this read group for all reads" << endl
         << "    -f, --max-frag-len N    reads with fragment lengths greater than N will not be marked properly paired in SAM/BAM/CRAM" << endl
         << "    -C, --compression N     level for compression [0-9]" << endl;
}

int main_surject(int argc, char** argv) {
    
    if (argc == 2) {
        help_surject(argv);
        return 1;
    }
    
    #define OPT_REF_PATHS 1001

    string xg_name;
    set<string> path_names;
    string path_file;
    string ref_paths_name;
    string output_format = "GAM";
    string input_format = "GAM";
    bool spliced = false;
    bool interleaved = false;
    string sample_name;
    string read_group;
    int32_t max_frag_len = 0;
    int compress_level = 9;
    int min_splice_length = 20;
    bool subpath_global = true; // force full length alignments in mpmap resolution
    bool qual_adj = false;
    bool prune_anchors = false;

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
            {"ref-paths", required_argument, 0, OPT_REF_PATHS},
            {"subpath-local", required_argument, 0, 'l'},
            {"interleaved", no_argument, 0, 'i'},
            {"gaf-input", no_argument, 0, 'G'},
            {"gamp-input", no_argument, 0, 'm'},
            {"cram-output", no_argument, 0, 'c'},
            {"bam-output", no_argument, 0, 'b'},
            {"sam-output", no_argument, 0, 's'},
            {"spliced", no_argument, 0, 'S'},
            {"prune-low-cplx", no_argument, 0, 'P'},
            {"qual-adj", no_argument, 0, 'A'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"max-frag-len", required_argument, 0, 'f'},
            {"compress", required_argument, 0, 'C'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:p:F:liGmcbsN:R:f:C:t:SPA",
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
            path_names.insert(optarg);
            break;

        case 'F':
            path_file = optarg;
            break;
        
        case OPT_REF_PATHS:
            ref_paths_name = optarg;
            break;

        case 'l':
            subpath_global = false;
            break;

        case 'i':
            interleaved = true;
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
            
        case 't':
            omp_set_num_threads(parse<int>(optarg));
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

    if (!path_file.empty()){
        // open the file
        ifstream in(path_file);
        string line;
        while (std::getline(in,line)) {
            // Load up all the paths to surject into from the file
            path_names.insert(line);
        }
    }

    PathPositionHandleGraph* xgidx = nullptr;
    unique_ptr<PathHandleGraph> path_handle_graph;
    bdsg::PathPositionOverlayHelper overlay_helper;
    if (!xg_name.empty()) {
        path_handle_graph = vg::io::VPKG::load_one<PathHandleGraph>(xg_name);
        xgidx = overlay_helper.apply(path_handle_graph.get());
    } else {
        // We need an XG index for the rest of the algorithm
        cerr << "error[vg surject] XG index (-x) is required for surjection" << endl;
        exit(1);
    }

    // add the target paths.  if there are no path names, add them all.  otherwise, take into account possible subpath names.
    unordered_set<path_handle_t> paths;
    xgidx->for_each_path_handle([&](path_handle_t path_handle) {
        if (path_names.empty() || path_names.count(Paths::get_base_name(xgidx->get_path_name(path_handle)))) {
            paths.insert(path_handle);
        }
        });
    // don't use this anymore: it's no longer updated to include all paths when none given. 
    path_names.clear();    

    // Make a single thread-safe Surjector.
    Surjector surjector(xgidx);
    surjector.adjust_alignments_for_base_quality = qual_adj;
    surjector.prune_suspicious_anchors = prune_anchors;
    if (spliced) {
        surjector.min_splice_length = min_splice_length;
        // we have to bump this up to be sure to align most splice junctions
        surjector.max_subgraph_bases = 1024 * 1024;
    }
    else {
        surjector.min_splice_length = numeric_limits<int64_t>::max();
    }
    
    // Get the paths to use in the HTSLib header sequence dictionary
    vector<tuple<path_handle_t, size_t, size_t>> sequence_dictionary = get_sequence_dictionary(ref_paths_name, *xgidx); 
   
    // Count our threads
    int thread_count = get_thread_count();
    
    if (input_format == "GAM" || input_format == "GAF") {
        
        // Give helpful warning if someone tries to surject an un-surjectable GAF
        auto check_gaf_aln = [&](const Alignment& src) {
            if (src.has_path() && src.sequence().empty()) {
#pragma omp critical
                {
                    cerr << "error:[surject] Read " << src.name() << " is aligned but does not have a sequence and therefore cannot be surjected. Was it derived from a GAF without a base-level alignment? or a GAF with a CIGAR string in the 'cg' tag (which does not provide enough information to reconstruct the sequence)?" << endl;
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
                
                
                // Preprocess read to set metadata before surjection
                set_metadata(src1);
                set_metadata(src2);
                
                // Surject and emit.
                alignment_emitter->emit_pair(surjector.surject(src1, paths, subpath_global, spliced),
                                             surjector.surject(src2, paths, subpath_global, spliced),
                                             max_frag_len);
                
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
                
                // Preprocess read to set metadata before surjection
                set_metadata(src);
                
                // Surject and emit the single read.
                alignment_emitter->emit_single(surjector.surject(src, paths, subpath_global, spliced));
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
                // we can reuse this
                vector<int64_t> tlen_limits(1, max_frag_len);
                
                // GAMP input is paired, and for HTS output reads need to know their pair partners' mapping locations.
                // TODO: We don't preserve order relationships (like primary/secondary) beyond the interleaving.
                vg::io::for_each_interleaved_pair_parallel<MultipathAlignment>(in, [&](MultipathAlignment& src1, MultipathAlignment& src2) {
                    
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
                    
                    // surject and record path positions
                    vector<pair<tuple<string, bool, int64_t>, tuple<string, bool, int64_t>>> positions(1);
                    vector<pair<multipath_alignment_t, multipath_alignment_t>> surjected;
                    surjected.emplace_back(surjector.surject(mp_src1, paths, get<0>(positions.front().first),
                                                             get<2>(positions.front().first), get<1>(positions.front().first),
                                                             subpath_global, spliced),
                                           surjector.surject(mp_src2, paths, get<0>(positions.front().second),
                                                             get<2>(positions.front().second), get<1>(positions.front().second),
                                                             subpath_global, spliced));
                    
                    // write to output
                    mp_alignment_emitter.emit_pairs(src1.name(), src2.name(), move(surjected), &positions, &tlen_limits);
                });
            } else {
                // TODO: We don't preserve order relationships (like primary/secondary).
                vg::io::for_each_parallel<MultipathAlignment>(in, [&](MultipathAlignment& src) {

                    multipath_alignment_t mp_src;
                    from_proto_multipath_alignment(src, mp_src);
                    
                    // surject and record path positions
                    vector<tuple<string, bool, int64_t>> positions(1);
                    vector<multipath_alignment_t> surjected;
                    surjected.emplace_back(surjector.surject(mp_src, paths, get<0>(positions.front()),
                                                             get<2>(positions.front()), get<1>(positions.front()),
                                                             subpath_global, spliced));
                    
                    // write to output
                    mp_alignment_emitter.emit_singles(src.name(), move(surjected), &positions);
                    
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
