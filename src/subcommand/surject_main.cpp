// surject_main.cpp: define the "vg surject" subcommand, which forces alignments into linear space

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>

#include "subcommand.hpp"

#include "../vg.hpp"
#include "../stream.hpp"
#include "../utility.hpp"
#include "../surjector.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

void help_surject(char** argv) {
    cerr << "usage: " << argv[0] << " surject [options] <aln.gam> >[proj.cram]" << endl
         << "Transforms alignments to be relative to particular paths." << endl
         << endl
         << "options:" << endl
         << "    -x, --xg-name FILE      use the graph in this xg index" << endl
         << "    -t, --threads N         number of threads to use" << endl
         << "    -p, --into-path NAME    surject into this path (many allowed, default: all in xg)" << endl
         << "    -F, --into-paths FILE   surject into nonoverlapping path names listed in FILE (one per line)" << endl
         << "    -f, --full-length       use the full length of the graph alignment, even if it makes a negative score" << endl
         << "    -i, --interleaved       GAM is interleaved paired-ended, so when outputting HTS formats, pair reads" << endl
         << "    -c, --cram-output       write CRAM to stdout" << endl
         << "    -b, --bam-output        write BAM to stdout" << endl
         << "    -s, --sam-output        write SAM to stdout" << endl
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
    string output_type = "gam";
    string input_type = "gam";
    bool interleaved = false;
    int compress_level = 9;
    bool full_length = false;

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
            {"full-length", required_argument, 0, 'f'},
            {"interleaved", no_argument, 0, 'i'},
            {"cram-output", no_argument, 0, 'c'},
            {"bam-output", no_argument, 0, 'b'},
            {"sam-output", no_argument, 0, 's'},
            {"compress", required_argument, 0, 'C'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long (argc, argv, "hx:p:F:ficbsC:t:",
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

        case 'f':
            full_length = true;
            break;

        case 'i':
            interleaved = true;
            break;
            
        case 'c':
            output_type = "cram";
            break;

        case 'b':
            output_type = "bam";
            break;

        case 's':
            compress_level = -1;
            output_type = "sam";
            break;

        case 't':
            omp_set_num_threads(parse<int>(optarg));
            break;

        case 'C':
            compress_level = parse<int>(optarg);
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

    string file_name = get_input_file_name(optind, argc, argv);

    if (!path_file.empty()){
        // open the file
        ifstream in(path_file);
        string line;
        while (std::getline(in,line)) {
            path_names.insert(line);
        }
    }

    xg::XG* xgidx = nullptr;
    ifstream xg_stream(xg_name);
    if(xg_stream) {
        xgidx = new xg::XG(xg_stream);
    }
    if (!xg_stream || xgidx == nullptr) {
        cerr << "[vg surject] error: could not open xg index" << endl;
        return 1;
    }

    // if no paths were given take all of those in the index
    if (path_names.empty()) {
        for (size_t i = 1; i <= xgidx->path_count; ++i) {
            path_names.insert(xgidx->path_name(i));
        }
    }

    map<string, int64_t> path_length;
    int num_paths = xgidx->max_path_rank();
    for (int i = 1; i <= num_paths; ++i) {
        auto name = xgidx->path_name(i);
        path_length[name] = xgidx->path_length(name);
    }
    
    int thread_count = get_thread_count();
    vector<Surjector*> surjectors(thread_count);
    for (int i = 0; i < surjectors.size(); i++) {
        surjectors[i] = new Surjector(xgidx);
    }

    if (input_type == "gam") {
        if (output_type == "gam") {
            vector<vector<Alignment> > buffer;
            buffer.resize(thread_count);
            function<void(Alignment&)> lambda = [&](Alignment& src) {
                int tid = omp_get_thread_num();
                // Since we're outputting full GAM, we ignore all this info
                // about where on the path the alignment falls. But we need to
                // provide the space to the surject call anyway.
                string path_name;
                int64_t path_pos;
                bool path_reverse;
                buffer[tid].push_back(surjectors[omp_get_thread_num()]->surject(src,
                                                                                path_names,
                                                                                path_name,
                                                                                path_pos,
                                                                                path_reverse,
                                                                                full_length));
                stream::write_buffered(cout, buffer[tid], 100);
            };
            get_input_file(file_name, [&](istream& in) {
                stream::for_each_parallel(in, lambda);
            });
            for (int i = 0; i < thread_count; ++i) {
                stream::write_buffered(cout, buffer[i], 0); // flush
            }
        } else {
            char out_mode[5];
            string out_format = "";
            strcpy(out_mode, "w");
            if (output_type == "bam") { out_format = "b"; }
            else if (output_type == "cram") { out_format = "c"; }
            else { out_format = ""; }
            strcat(out_mode, out_format.c_str());
            if (compress_level >= 0) {
                char tmp[2];
                tmp[0] = compress_level + '0'; tmp[1] = '\0';
                strcat(out_mode, tmp);
            }
            
            int thread_count = get_thread_count();
            
            // bam/sam/cram output
            
            // Define a string to hold the SAM header, to be generated later.
            string header;
            // To generate the header, we need to know the read group for each sample name.
            map<string, string> rg_sample;
            
            samFile* out = nullptr;
            int buffer_limit = 100;

            bam_hdr_t* hdr = nullptr;
            int64_t count = 0;
            // TODO: What good is this lock if we continue without getting it if the buffer is overfull???
            omp_lock_t output_lock;
            omp_init_lock(&output_lock);
            
            // We define a type to represent a surjected alignment, ready for
            // HTSlib output. It consists of surjected path name (or ""),
            // surjected position (or -1), surjected orientation, and the
            // actual Alignment.
            using surjected_t = tuple<string, int64_t, bool, Alignment>;
            // You make one with make_tuple()
            
            // We define a basic surject function, which also fills in the read group info we need to make the header
            auto surject_alignment = [&](const Alignment& src) {
                
                // Set out some variables to populate with the linear position
                // info we need for SAM/BAM/CRAM
                string path_name;
                // Make sure to initialize pos to -1 since it may not be set by
                // surject_alignment if the read is unmapped, and unmapped
                // reads need to come out with a 0 1-based position.
                int64_t path_pos = -1; 
                bool path_reverse = false;
                auto surj = surjectors[omp_get_thread_num()]->surject(src,
                                                                      path_names,
                                                                      path_name,
                                                                      path_pos,
                                                                      path_reverse,
                                                                      full_length);
                // Always use the surjected alignment, even if it surjects to unmapped.
                
                if (!hdr && !surj.read_group().empty() && !surj.sample_name().empty()) {
                    // There's no header yet (although we race its
                    // construction) and we have a sample and a read group.
                    
                    // Record the read group for the sample that this read
                    // represents, so that when we build the header we list it.
#pragma omp critical (hts_header)
                    rg_sample[surj.read_group()] = surj.sample_name();
                }
                
                return make_tuple(path_name, path_pos, path_reverse, surj);
            };
            
            // We also define a function to emit the header if it hasn't been made already.
            // Note that the header will only list the samples and read groups in reads we have encountered so far!
            auto ensure_header = [&]() {
#pragma omp critical (hts_header)
                {
                    if (!hdr) {
                        hdr = hts_string_header(header, path_length, rg_sample);
                        if ((out = sam_open("-", out_mode)) == 0) {
#pragma omp critical (cerr)
                            cerr << "[vg surject] error: failed to open stdout for writing HTS output" << endl;
                            exit(1);
                        } else {
                            // write the header
                            if (sam_hdr_write(out, hdr) != 0) {
#pragma omp critical (cerr)
                                cerr << "[vg surject] error: failed to write the SAM header" << endl;
                            }
                        }
                    }
                }
            };
            
            // Finally, we have a little widget function to write a BAM record and check for errors.
            // Consumes the passed record.
            auto write_bam_record = [&](bam1_t* b) {
                assert(out != nullptr);
                int r = 0;
#pragma omp critical (cout)
                r = sam_write1(out, hdr, b);
                if (r == 0) {
#pragma omp critical (cerr)
                    cerr << "[vg surject] error: writing to stdout failed" << endl;
                    exit(1);
                }
                bam_destroy1(b);
            };
            
            if (interleaved) {
                // GAM input is paired, and for HTS output reads need to know their pair partners' mapping locations              
                
                // We keep a buffer, one per thread, of pairs of surjected alignments
                vector<vector<pair<surjected_t, surjected_t>>> buffer;
                buffer.resize(thread_count);
                
                // Define a function that handles buffers, possibly opening the
                // output file if we're on the first record
                auto handle_buffer = [&](vector<pair<surjected_t, surjected_t>>& buf) {
                     if (buf.size() >= buffer_limit) {
                            // We have enough data to start the file.
                            
                            // Make sure we have emitted the header
                            ensure_header();

                            // try to get a lock, and force things if we've built up a huge buffer waiting
                            // TODO: Is continuing without the lock safe? And if so why do we have the lock in the first place?
                            if (omp_test_lock(&output_lock) || buf.size() > 10*buffer_limit) {
                                for (auto& surjected_pair : buf) {
                                    // For each pair of surjected reads
                                
                                    // Unpack the first read
                                    auto& name1 = get<0>(surjected_pair.first);
                                    auto& pos1 = get<1>(surjected_pair.first);
                                    auto& reverse1 = get<2>(surjected_pair.first);
                                    auto& surj1 = get<3>(surjected_pair.first);
                                    
                                    // Unpack the second read
                                    auto& name2 = get<0>(surjected_pair.second);
                                    auto& pos2 = get<1>(surjected_pair.second);
                                    auto& reverse2 = get<2>(surjected_pair.second);
                                    auto& surj2 = get<3>(surjected_pair.second);
                                    
                                    // Compute CIGAR strings if actually surjected
                                    string cigar1 = "", cigar2 = "";
                                    if (name1 != "") {
                                        size_t path_len1 = xgidx->path_length(name1);
                                        cigar1 = cigar_against_path(surj1, reverse1, pos1, path_len1, 0);
                                    }
                                    if (name2 != "") {
                                        size_t path_len2 = xgidx->path_length(name2);
                                        cigar2 = cigar_against_path(surj2, reverse2, pos2, path_len2, 0);
                                    }
                                    
                                    // TODO: compute template length based on
                                    // pair distance and alignment content.
                                    int template_length = 0;
                                    
                                    // Create and write paired BAM records referencing each other
                                    write_bam_record(alignment_to_bam(header, surj1, name1, pos1, reverse1, cigar1,
                                        name2, pos2, template_length));
                                    write_bam_record(alignment_to_bam(header, surj2, name2, pos2, reverse2, cigar2,
                                        name1, pos1, template_length));
                                
                                }
                                
                                omp_unset_lock(&output_lock);
                                buf.clear();
                            }
                        }
                };
                
                // Define a function to surject the pair and fill in the
                // HTSlib-required crossreferences
                function<void(Alignment&,Alignment&)> lambda = [&](Alignment& aln1, Alignment& aln2) {
                    // Make sure that the alignments being surjected are
                    // actually paired with each other (proper
                    // fragment_prev/fragment_next). We want to catch people
                    // giving us un-interleaved GAMs as interleaved.
                    if (aln1.has_fragment_next()) {
                        // Alignment 1 comes first in fragment
                        if (aln1.fragment_next().name() != aln2.name() ||
                            !aln2.has_fragment_prev() ||
                            aln2.fragment_prev().name() != aln1.name()) {
                        
#pragma omp critical (cerr)
                            cerr << "[vg surject] error: alignments " << aln1.name()
                                 << " and " << aln2.name() << " are adjacent but not paired" << endl;
                                 
                            exit(1);
                        
                        }
                    } else if (aln2.has_fragment_next()) {
                        // Alignment 2 comes first in fragment
                        if (aln2.fragment_next().name() != aln1.name() ||
                            !aln1.has_fragment_prev() ||
                            aln1.fragment_prev().name() != aln2.name()) {
                        
#pragma omp critical (cerr)
                            cerr << "[vg surject] error: alignments " << aln1.name()
                                 << " and " << aln2.name() << " are adjacent but not paired" << endl;
                                 
                            exit(1);
                        
                        }
                    } else {
                        // Alignments aren't paired up at all
#pragma omp critical (cerr)
                        cerr << "[vg surject] error: alignments " << aln1.name()
                             << " and " << aln2.name() << " are adjacent but not paired" << endl;
                             
                        exit(1);
                    }
                    
                    // Find our buffer
                    auto& thread_buffer = buffer[omp_get_thread_num()];
                    // Surject each of the pair and buffer the surjected pair
                    thread_buffer.emplace_back(surject_alignment(aln1), surject_alignment(aln2));
                    // Spit out the buffer if (over)full
                    handle_buffer(thread_buffer);
                };
                
                // now apply the alignment processor to the stream
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_interleaved_pair_parallel(in, lambda);
                });
                
                // Spit out any remaining data
                buffer_limit = 0;
                for (auto& buf : buffer) {
                    handle_buffer(buf);
                }
                
            } else {
                // GAM input is single-ended, so each read can be surjected
                // independently
                
                // We keep a buffer, one per thread, of surjected alignments.
                vector<vector<surjected_t>> buffer;
                buffer.resize(thread_count);

                // Define a function that handles buffers, possibly opening the
                // output file if we're on the first record
                auto handle_buffer = [&](vector<surjected_t>& buf) {
                    if (buf.size() >= buffer_limit) {
                        // We have enough data to start the file.
                        
                        // Make sure we have emitted the header
                        ensure_header();

                        // try to get a lock, and force things if we've built up a huge buffer waiting
                        // TODO: Is continuing without the lock safe? And if so why do we have the lock in the first place?
                        if (omp_test_lock(&output_lock) || buf.size() > 10*buffer_limit) {
                            for (auto& s : buf) {
                                // For each alignment in the buffer
                                
                                // Unpack it
                                auto& name = get<0>(s);
                                auto& pos = get<1>(s);
                                auto& reverse = get<2>(s);
                                auto& surj = get<3>(s);
                                
                                // Generate a CIGAR string for it
                                string cigar = "";
                                if (name != "") {
                                    size_t path_len = xgidx->path_length(name);
                                    cigar = cigar_against_path(surj, reverse, pos, path_len, 0);
                                }
                                
                                // Create and write a single unpaired BAM record
                                write_bam_record(alignment_to_bam(header, surj, name, pos, reverse, cigar));
                                
                            }
                            
                            omp_unset_lock(&output_lock);
                            buf.clear();
                        }
                    }
                };

                function<void(Alignment&)> lambda = [&](Alignment& src) {
                    auto& thread_buffer = buffer[omp_get_thread_num()];
                    thread_buffer.push_back(surject_alignment(src));
                    handle_buffer(thread_buffer);
                };


                // now apply the alignment processor to the stream
                get_input_file(file_name, [&](istream& in) {
                    stream::for_each_parallel(in, lambda);
                });
                buffer_limit = 0;
                for (auto& buf : buffer) {
                    handle_buffer(buf);
                }
                
            }
            
            
            if (hdr != nullptr) {
                bam_hdr_destroy(hdr);
            }
            assert(out != nullptr);
            sam_close(out);
            omp_destroy_lock(&output_lock);
        }
    }
    cout.flush();
    
    for (Surjector* surjector : surjectors) {
        delete surjector;
    }

    return 0;
}

// Register subcommand
static Subcommand vg_surject("surject", "map alignments onto specific paths", main_surject);
