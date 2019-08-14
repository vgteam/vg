// surject_main.cpp: define the "vg surject" subcommand, which forces alignments into linear space

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../xg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include "../utility.hpp"
#include "../surjector.hpp"
#include "../alignment_emitter.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_surject(char** argv) {
    cerr << "usage: " << argv[0] << " surject [options] <aln.gam> >[proj.cram]" << endl
         << "Transforms alignments to be relative to particular paths." << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE      use the graph in this xg index (required)" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << "    -p, --into-path NAME    surject into this path (many allowed, default: all in xg)" << endl
         << "    -F, --into-paths FILE   surject into nonoverlapping path names listed in FILE (one per line)" << endl
         << "    -l, --subpath-local     let the multipath mapping surjection produce local (rather than global) alignments" << endl
         << "    -i, --interleaved       GAM is interleaved paired-ended, so when outputting HTS formats, pair reads" << endl
         << "    -c, --cram-output       write CRAM to stdout" << endl
         << "    -b, --bam-output        write BAM to stdout" << endl
         << "    -s, --sam-output        write SAM to stdout" << endl
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

    string xg_name;
    set<string> path_names;
    string path_file;
    string output_format = "GAM";
    string input_format = "GAM";
    bool interleaved = false;
    string sample_name;
    string read_group;
    int32_t max_frag_len = 0;
    int compress_level = 9;
    bool subpath_global = true; // force full length alignments in mpmap resolution

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
            {"subpath-local", required_argument, 0, 'l'},
            {"interleaved", no_argument, 0, 'i'},
            {"cram-output", no_argument, 0, 'c'},
            {"bam-output", no_argument, 0, 'b'},
            {"sam-output", no_argument, 0, 's'},
            {"sample", required_argument, 0, 'N'},
            {"read-group", required_argument, 0, 'R'},
            {"max-frag-len", required_argument, 0, 'f'},
            {"compress", required_argument, 0, 'C'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:p:F:licbsN:R:f:C:t:",
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

        case 'l':
            subpath_global = false;
            break;

        case 'i':
            interleaved = true;
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

    unique_ptr<PathPositionHandleGraph> xgidx;
    if (!xg_name.empty()) {
        xgidx = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_name);
    } else {
        // We need an XG index for the rest of the algorithm
        cerr << "error[vg surject] XG index (-x) is required for surjection" << endl;
        exit(1);
    }

    // if no paths were given take all of those in the index
    if (path_names.empty()) {
        xgidx->for_each_path_handle([&](path_handle_t path_handle) {
                path_names.insert(xgidx->get_path_name(path_handle));
            });
    }

    // Make a single therad-safe Surjector.
    Surjector surjector(xgidx.get());
    
    // Get the lengths of all the paths in the XG to populate the HTS headers
    map<string, int64_t> path_length;
    xgidx->for_each_path_handle([&](path_handle_t path_handle) {
            path_length[xgidx->get_path_name(path_handle)] = xgidx->get_path_length(path_handle);
        });
   
    // Count our threads
    int thread_count = get_thread_count();
   
    // Set up output to an emitter that will handle serialization
    unique_ptr<AlignmentEmitter> alignment_emitter = get_alignment_emitter("-", output_format, path_length, thread_count);

    if (input_format == "GAM") {
        get_input_file(file_name, [&](istream& in) {
            if (interleaved) {
                // GAM input is paired, and for HTS output reads need to know their pair partners' mapping locations.
                // TODO: We don't preserve order relationships (like primary/secondary) beyond the interleaving.
                vg::io::for_each_interleaved_pair_parallel<Alignment>(in, [&](Alignment& src1, Alignment& src2) {
               
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
                    alignment_emitter->emit_pair(surjector.surject(src1, path_names, subpath_global),
                                                 surjector.surject(src2, path_names, subpath_global));
                
                });
            } else {
                // We can just surject each Alignment by itself.
                // TODO: We don't preserve order relationships (like primary/secondary).
                vg::io::for_each_parallel<Alignment>(in, [&](Alignment& src) {
                
                    // Preprocess read to set metadata before surjection
                    set_metadata(src);
                    
                    // Surject and emit the single read.
                    alignment_emitter->emit_single(surjector.surject(src, path_names, subpath_global));
                        
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
