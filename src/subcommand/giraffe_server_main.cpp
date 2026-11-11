/**
 * \file giraffe_server_main.cpp
 * vg giraffe-server: load the Giraffe indexes once, then map reads streamed on
 * stdin and write GAF to stdout.
 */

#include <getopt.h>
#include <unistd.h>

#include <condition_variable>
#include <deque>
#include <functional>
#include <iostream>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "subcommand.hpp"
#include "giraffe_options.hpp"

#include "../giraffe_engine.hpp"
#include <vg/io/alignment_io.hpp>
#include "../utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

namespace {

// Codes for the long options with no short form. They start above 255, so
// they cannot clash with a short option character.
constexpr int OPT_EMIT_HEADER   = 1000;
constexpr int OPT_FRAMED_OUTPUT = 1001;
constexpr int OPT_SURJECT_THREADS = 1002;
constexpr int OPT_PARAMETER_PRESET = 1003;
constexpr int OPT_SET = 1004;
constexpr int OPT_NO_REC_MODE = 1005;

/// A fixed set of threads that run SURJECT_WITH_ANCHORS requests, so several
/// surjections run at once while the main thread keeps reading stdin.
class SurjectPool {
public:
    explicit SurjectPool(size_t n) {
        for (size_t i = 0; i < n; ++i) {
            workers_.emplace_back([this]() { run(); });
        }
    }
    ~SurjectPool() { drain(); }

    void submit(function<void()> task) {
        {
            lock_guard<mutex> lk(mu_);
            tasks_.push_back(move(task));
        }
        cv_.notify_one();
    }

    /// Finish every queued task, then stop the threads.
    void drain() {
        {
            lock_guard<mutex> lk(mu_);
            if (stopping_) return;
            stopping_ = true;
        }
        cv_.notify_all();
        for (auto& w : workers_) w.join();
    }

private:
    void run() {
        while (true) {
            function<void()> task;
            {
                unique_lock<mutex> lk(mu_);
                cv_.wait(lk, [this]() { return stopping_ || !tasks_.empty(); });
                if (tasks_.empty()) return;          // stopping, and nothing is left
                task = move(tasks_.front());
                tasks_.pop_front();
            }
            task();
        }
    }

    vector<thread> workers_;
    deque<function<void()>> tasks_;
    mutex mu_;
    condition_variable cv_;
    bool stopping_ = false;
};

void help_giraffe_server(char** argv) {
    cerr << "usage: " << argv[0] << " giraffe-server [options]" << endl
         << "Loads Giraffe indexes once and maps sequence batches from stdin." << endl
         << endl
         << "required options:" << endl
         << "  -Z, --gbz FILE               GBZ index path" << endl
         << "  -m, --minimizer FILE         minimizer index path" << endl
         << "  -d, --distance FILE          distance index path" << endl
         << "  -z, --zipcodes FILE          zipcode index path" << endl
         << endl
         << "optional:" << endl
         << "  -t, --threads N              mapping threads [1]" << endl
         << "  -M, --max-multimaps N        max mappings per read (reporting is capped by" << endl
         << "                               max-alignments) [the preset's value]" << endl
         << "  -b, --batch-size N           input reads per mapping batch [256]" << endl
         << "      --emit-header            emit GAF header lines before output" << endl
         << "      --framed-output          emit read-grouped framed output for middleware" << endl
         << "      --surject-threads N      run SURJECT_WITH_ANCHORS requests on N worker" << endl
         << "                               threads (framed output only; 0 = inline, one at" << endl
         << "                               a time) [--threads]" << endl
         << "      --parameter-preset NAME  vg giraffe parameter preset [hifi]" << endl
         << "      --set NAME[=VALUE]       set any vg giraffe option by its long name," << endl
         << "                               e.g. mapq-score-window=150; a flag takes no" << endl
         << "                               value (may repeat)" << endl
         << "  -E, --rec-mode               map recombination-aware, as vg giraffe does" << endl
         << "      --no-rec-mode            don't (default: on for long-read presets when" << endl
         << "                               the minimizer index stores paths)" << endl
         << "  -S, --surject-target NAME    pre-index this haplotype path for surjection" << endl
         << "                               (may be repeated; required for per-read targets)" << endl
         << "  -h, --help                   show help" << endl
         << endl
         << "stdin protocol (one read per line):" << endl
         << "  SEQUENCE" << endl
         << "  NAME<TAB>SEQUENCE" << endl
         << "  NAME<TAB>SEQUENCE<TAB>QUALITY" << endl
         << "  NAME<TAB>SEQUENCE<TAB>QUALITY<TAB>SURJ_TARGET" << endl
         << "                               (QUALITY may be empty when SURJ_TARGET is set)" << endl
         << "  QUALITY is FASTQ phred+33 ASCII; if non-empty, must match SEQUENCE length." << endl
         << "  Without QUALITY, MAPQ is not capped by base quality (as for FASTA input)." << endl
         << "  PROCESS_BATCH   (map any buffered reads now; alias: FLUSH_NOW)" << endl
         << endl
         << "  SURJECT_WITH_ANCHORS<TAB>NAME<TAB>TARGET<TAB>PATH_LEN<TAB>N_ANCHORS" << endl
         << "    followed by N_ANCHORS tab-separated anchor lines (10 fields each):" << endl
         << "      step_begin_node, step_begin_offset, step_end_node, step_end_offset," << endl
         << "      path_off_begin, path_off_end, read_begin, read_end," << endl
         << "      src_mapping_begin, src_mapping_end (all tab-separated)" << endl
         << "    followed by ONE graph-alignment GAF line." << endl
         << "    Response is framed: READ<TAB>NAME<TAB>1 then one GAF line" << endl
         << "    carrying the surjection tags (sj:Z:..., sn:Z:..., etc.)." << endl
         << endl
         << "stdout:" << endl
         << "  GAF header lines first, with --emit-header. Then, per read in input order," << endl
         << "  one GAF line per reported alignment (an unmapped read gets one unmapped line)." << endl
         << "  Default: the GAF lines only; reads with errors are skipped (see stderr)." << endl
         << "  --framed-output: READ<TAB>NAME<TAB>N, then the read's N GAF lines. N is 0" << endl
         << "    for a read that could not be mapped because of an error (see stderr)." << endl
         << "  Each line has extra tags:" << endl
         << "    hp:Z haplotypes carrying the path (comma-separated)  hn:i their count" << endl
         << "    hb:Z representative haplotype  hf:i 1 if it is a reference" << endl
         << "    hc:i carrying GBWT sequences  ht:i 1 if hp was capped" << endl
         << "    hl:i, hr:i enclosing snarl's boundary nodes (only when found)" << endl
         << "    hq:i 1 if one haplotype carries the whole path  hv:i % of the path it covers" << endl
         << "    hs:i haplotype-consistent segments  hm:Z bp:count:name per segment" << endl
         << "  With a SURJ_TARGET: sj:Z status (ok, not_indexed, incompatible, failed, ...)," << endl
         << "    and if ok sn:Z path, sp:i position, sr:i reverse, ss:i score, sm:i MAPQ," << endl
         << "    sc:Z CIGAR, sM/sX/sI/sD:i matches, mismatches, insertions, deletions." << endl
         << endl;
}

/// Characters allowed in a phred+33 quality string.
const string PHRED33_CHARS = [] {
    string chars;
    for (char c = '!'; c <= '~'; ++c) chars.push_back(c);
    return chars;
}();

/// Return true if seq is non-empty and holds only ATCGN
bool is_nucleotide_sequence(const string& seq) {
    static const string allowed = "ATCGNatcgn";
    return !seq.empty() && seq.find_first_not_of(allowed) == string::npos;
}

vector<string> split_tab_fields(const string& line) {
    vector<string> fields;
    size_t start = 0;
    while (true) {
        size_t tab = line.find('\t', start);
        if (tab == string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    return fields;
}

/// Return the number in a request field, or nullopt if the field is not a
/// decimal number that fits in a size_t. vg's parse<T>(string) exits the
/// program on bad input, and a bad request must not end the server.
optional<size_t> parse_request_number(const string& field) {
    size_t value;
    if (field.empty() || field.find_first_not_of("0123456789") != string::npos) {
        return nullopt;
    }
    try {
        if (parse<size_t>(field, value)) {
            return value;
        }
    } catch (const exception&) {
        // std::stoull throws on a number too big for its type.
    }
    return nullopt;
}

int main_giraffe_server(int argc, char** argv) {
    GiraffeEngineConfig config;
    size_t batch_size = 256;
    bool emit_header = false;
    bool framed_output = false;
    long surject_threads = -1;   // -1 means the --threads value

    GiraffeParameters params;
    params.preset = "hifi";

    int c;
    optind = 2;
    while (true) {
        static struct option long_options[] = {
            {"gbz", required_argument, 0, 'Z'},
            {"minimizer", required_argument, 0, 'm'},
            {"distance", required_argument, 0, 'd'},
            {"zipcodes", required_argument, 0, 'z'},
            {"threads", required_argument, 0, 't'},
            {"max-multimaps", required_argument, 0, 'M'},
            {"batch-size", required_argument, 0, 'b'},
            {"emit-header", no_argument, 0, OPT_EMIT_HEADER},
            {"framed-output", no_argument, 0, OPT_FRAMED_OUTPUT},
            {"surject-threads", required_argument, 0, OPT_SURJECT_THREADS},
            {"parameter-preset", required_argument, 0, OPT_PARAMETER_PRESET},
            {"set", required_argument, 0, OPT_SET},
            {"rec-mode", no_argument, 0, 'E'},
            {"no-rec-mode", no_argument, 0, OPT_NO_REC_MODE},
            {"surject-target", required_argument, 0, 'S'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "Z:m:d:z:t:M:b:S:Eh?", long_options, &option_index);
        if (c == -1) {
            break;
        }

        switch (c) {
            case 'Z':
                config.gbz_path = optarg;
                break;
            case 'm':
                config.minimizer_path = optarg;
                break;
            case 'd':
                config.distance_path = optarg;
                break;
            case 'z':
                config.zipcode_path = optarg;
                break;
            case 't':
                config.threads = parse<size_t>(optarg);
                break;
            case 'M':
                params.options.emplace_back("max-multimaps", optarg);
                break;
            case 'b':
                batch_size = parse<size_t>(optarg);
                break;
            case OPT_EMIT_HEADER:
                emit_header = true;
                break;
            case OPT_FRAMED_OUTPUT:
                framed_output = true;
                break;
            case OPT_SURJECT_THREADS:
                surject_threads = parse<long>(optarg);
                break;
            case OPT_PARAMETER_PRESET:
                params.preset = optarg;
                break;
            case OPT_SET: {
                // NAME=VALUE sets an option, and NAME alone turns on a flag.
                std::string kv = optarg;
                size_t eq = kv.find('=');
                if (kv.empty() || eq == 0) {
                    cerr << "error [vg giraffe-server]: --set expects NAME[=VALUE], got " << kv << endl;
                    return 1;
                }
                if (eq == std::string::npos) {
                    params.options.emplace_back(kv, "");
                } else {
                    params.options.emplace_back(kv.substr(0, eq), kv.substr(eq + 1));
                }
                break;
            }
            case 'S':
                config.surjection_target_paths.emplace_back(optarg);
                break;
            case 'E':
                params.rec_mode = true;
                break;
            case OPT_NO_REC_MODE:
                params.rec_mode = false;
                break;
            case 'h':
            case '?':
            default:
                help_giraffe_server(argv);
                return 1;
        }
    }

    if (argc > optind || config.gbz_path.empty() || config.minimizer_path.empty()
        || config.distance_path.empty() || config.zipcode_path.empty() || batch_size == 0) {
        help_giraffe_server(argv);
        return 1;
    }

    // Fill vg giraffe's option parser now, so bad parameters fail before any
    // index is loaded.
    std::unique_ptr<GroupedOptionGroup> parser = make_giraffe_parser(params);
    cerr << "[vg giraffe-server] parameter preset: " << params.preset;
    for (const auto& kv : params.options) {
        cerr << " --" << kv.first;
        if (!kv.second.empty()) {
            cerr << " " << kv.second;
        }
    }
    cerr << endl;

    try {
        GiraffeEngine engine;
        engine.load(config);

        // Set the mapper's parameters. Recombination-aware mode depends on the
        // loaded minimizer index.
        bool use_rec_mode = giraffe_rec_mode(engine.minimizer_index(),
            giraffe_preset_maps_long_reads(params.preset), params.rec_mode);
        set_up_giraffe_mapper(*parser, use_rec_mode, engine.mapper());
        cerr << "[vg giraffe-server] recombination-aware mode: " << (use_rec_mode ? "on" : "off") << endl;

        if (emit_header) {
            for (const auto& header : engine.gaf_header_lines()) {
                cout << header << '\n';
            }
            cout.flush();
        }

        // Every write to stdout holds out_mu, so a READ line and the read's GAF
        // lines never mix with another thread's output.
        mutex out_mu;

        // Surjections run on threads only with framed output. Each frame names
        // its read, so a client can match up answers that come out of order.
        const size_t n_surject = (surject_threads < 0)
            ? static_cast<size_t>(config.threads) : static_cast<size_t>(surject_threads);
        unique_ptr<SurjectPool> surject_pool;
        if (framed_output && n_surject > 0) {
            surject_pool = make_unique<SurjectPool>(n_surject);
        }

        vector<GiraffeServerRequest> batch;
        batch.reserve(batch_size);

        string line;
        size_t auto_id = 0;
        auto flush_batch = [&]() {
            if (batch.empty()) {
                return;
            }
            auto mapped = engine.map_reads(batch);
            lock_guard<mutex> lk(out_mu);
            for (size_t i = 0; i < mapped.size(); ++i) {
                const auto& read_mappings = mapped[i];
                if (framed_output) {
                    cout << "READ\t" << batch[i].read.name() << '\t' << read_mappings.size() << '\n';
                }
                for (const auto& gaf_line : read_mappings) {
                    cout << gaf_line << '\n';
                }
            }
            cout.flush();
            batch.clear();
        };

        // Report a bad read on stderr. With framed output, also send an empty
        // frame, so a client waiting for the read gets an answer.
        auto emit_read_error = [&](const string& name, const string& reason) {
            cerr << "warning [vg giraffe-server]: skipping read"
                 << (name.empty() ? "" : " '" + name + "'") << ": " << reason << endl;
            if (framed_output) {
                lock_guard<mutex> lk(out_mu);
                cout << "READ\t" << name << "\t0\n";
                cout.flush();
            }
        };

        auto emit_empty_frame = [&](const string& name) {
            if (framed_output) {
                lock_guard<mutex> lk(out_mu);
                cout << "READ\t" << name << "\t0\n";
                cout.flush();
            }
        };

        size_t line_no = 0;
        while (getline(cin, line)) {
            ++line_no;
            if (line.empty()) {
                continue;
            }
            // PROCESS_BATCH maps the reads buffered so far. FLUSH_NOW is its old name.
            if (line == "PROCESS_BATCH" || line == "FLUSH_NOW") {
                flush_batch();
                continue;
            }
            // A SURJECT_WITH_ANCHORS request is a header line, N anchor lines and
            // one GAF line. Map the buffered reads first, so their frames come out
            // before this one.
            if (line.compare(0, sizeof("SURJECT_WITH_ANCHORS\t") - 1,
                             "SURJECT_WITH_ANCHORS\t") == 0) {
                flush_batch();
                auto header_fields = split_tab_fields(line);
                // The header holds the command, read name, target path, path
                // length and anchor count.
                if (header_fields.size() != 5) {
                    cerr << "error [vg giraffe-server]: malformed SURJECT_WITH_ANCHORS header "
                            "(need 5 tab fields, got " << header_fields.size() << ")" << endl;
                    if (header_fields.size() >= 2) emit_empty_frame(header_fields[1]);
                    continue;
                }
                const string& sa_name = header_fields[1];
                const string& sa_target = header_fields[2];
                optional<size_t> sa_path_len = parse_request_number(header_fields[3]);
                optional<size_t> sa_n_anchors = parse_request_number(header_fields[4]);
                if (!sa_n_anchors) {
                    // Without the anchor count we cannot skip the request's
                    // other lines. They are not reads, so the read parser
                    // ignores them with a warning.
                    cerr << "error [vg giraffe-server]: malformed SURJECT_WITH_ANCHORS header "
                            "(anchor count must be a number)" << endl;
                    emit_empty_frame(sa_name);
                    continue;
                }

                // Read all of the request's lines even after an error, so the
                // next request starts on its own line.
                vector<WireAnchor> sa_anchors;
                string sa_error;
                if (!sa_path_len) {
                    sa_error = "malformed SURJECT_WITH_ANCHORS header (path length must be a number)";
                }
                bool sa_eof = false;
                for (size_t ai = 0; ai < *sa_n_anchors; ++ai) {
                    string anchor_line;
                    if (!getline(cin, anchor_line)) {
                        cerr << "error [vg giraffe-server]: SURJECT_WITH_ANCHORS hit EOF "
                             << "while reading anchor " << ai << "/" << *sa_n_anchors << endl;
                        sa_eof = true;
                        break;
                    }
                    ++line_no;
                    if (!sa_error.empty()) {
                        continue;
                    }
                    auto af = split_tab_fields(anchor_line);
                    vector<size_t> nums;
                    for (const auto& field : af) {
                        optional<size_t> num = parse_request_number(field);
                        if (!num) {
                            break;
                        }
                        nums.push_back(*num);
                    }
                    if (af.size() != 10 || nums.size() != 10) {
                        sa_error = "malformed anchor line " + to_string(line_no)
                                   + " (need 10 tab-separated numbers)";
                        continue;
                    }
                    WireAnchor wa;
                    wa.step_begin_node        = nums[0];
                    wa.step_begin_offset      = nums[1];
                    wa.step_end_node          = nums[2];
                    wa.step_end_offset        = nums[3];
                    wa.path_offset_step_begin = nums[4];
                    wa.path_offset_step_end   = nums[5];
                    wa.read_begin_offset      = nums[6];
                    wa.read_end_offset        = nums[7];
                    wa.source_mapping_begin   = nums[8];
                    wa.source_mapping_end     = nums[9];
                    sa_anchors.push_back(wa);
                }

                string sa_gaf_line;
                if (sa_eof || !getline(cin, sa_gaf_line)) {
                    if (!sa_eof) {
                        cerr << "error [vg giraffe-server]: SURJECT_WITH_ANCHORS hit EOF "
                             << "while reading the GAF line" << endl;
                    }
                    emit_empty_frame(sa_name);
                    continue;
                }
                ++line_no;
                if (!sa_error.empty()) {
                    cerr << "error [vg giraffe-server]: " << sa_error << endl;
                    emit_empty_frame(sa_name);
                    continue;
                }

                // Parsing stays on this thread, which reads stdin. The surjection
                // needs nothing from other requests, so it can run on the pool.
                auto task = [&engine, &out_mu, framed_output,
                             name = sa_name, target = sa_target,
                             gaf = move(sa_gaf_line), anchors = move(sa_anchors),
                             path_len = *sa_path_len]() {
                    vector<string> out_lines;
                    try {
                        out_lines = engine.surject_with_anchors(
                            name, gaf, anchors, target, path_len);
                    } catch (const exception& e) {
                        lock_guard<mutex> lk(out_mu);
                        cerr << "error [vg giraffe-server]: surject_with_anchors failed for "
                             << "read '" << name << "': " << e.what() << endl;
                        if (framed_output) {
                            cout << "READ\t" << name << "\t0\n";
                            cout.flush();
                        }
                        return;
                    } catch (...) {
                        // An exception that escapes a pool thread terminates the
                        // server.
                        lock_guard<mutex> lk(out_mu);
                        cerr << "error [vg giraffe-server]: surject_with_anchors failed for "
                             << "read '" << name << "' (non-standard exception)" << endl;
                        if (framed_output) {
                            cout << "READ\t" << name << "\t0\n";
                            cout.flush();
                        }
                        return;
                    }
                    lock_guard<mutex> lk(out_mu);
                    if (framed_output) {
                        cout << "READ\t" << name << "\t" << out_lines.size() << "\n";
                    }
                    for (const auto& gl : out_lines) {
                        cout << gl << "\n";
                    }
                    cout.flush();
                };
                if (surject_pool) {
                    surject_pool->submit(move(task));
                } else {
                    task();
                }
                continue;
            }
            auto fields = split_tab_fields(line);
            // A read line has 1 to 4 tab-separated fields:
            //   1: SEQ
            //   2: NAME, SEQ
            //   3: NAME, SEQ, QUAL
            //   4: NAME, SEQ, QUAL, SURJ_TARGET (QUAL may be empty)
            if (fields.size() < 1 || fields.size() > 4) {
                cerr << "warning [vg giraffe-server]: line " << line_no
                     << ": ignored (expected 1-4 tab-separated fields or a known command, got "
                     << fields.size() << ")" << endl;
                continue;
            }

            // A lone field must be a sequence. Anything else is an unknown command.
            if (fields.size() == 1 && !is_nucleotide_sequence(fields[0])) {
                cerr << "warning [vg giraffe-server]: line " << line_no
                     << ": ignored (not a known command or a nucleotide sequence)" << endl;
                continue;
            }

            const string name = fields.size() == 1 ? "read_" + to_string(auto_id++) : fields[0];
            const string& sequence = fields.size() == 1 ? fields[0] : fields[1];
            const string quality = fields.size() >= 3 ? fields[2] : "";

            if (sequence.empty()) {
                emit_read_error(name, "empty SEQUENCE");
                continue;
            }
            if (!is_nucleotide_sequence(sequence)) {
                emit_read_error(name, "SEQUENCE has characters that are not IUPAC nucleotides");
                continue;
            }
            if (quality.find_first_not_of(PHRED33_CHARS) != string::npos) {
                emit_read_error(name, "QUALITY has characters outside phred+33 ('!' to '~')");
                continue;
            }
            if (!quality.empty() && quality.size() != sequence.size()) {
                emit_read_error(name,
                    "QUALITY length (" + to_string(quality.size())
                    + ") does not match SEQUENCE length ("
                    + to_string(sequence.size()) + ")");
                continue;
            }

            GiraffeServerRequest read;
            read.read.set_name(name);
            read.read.set_sequence(sequence);
            read.read.set_quality(vg::io::string_quality_char_to_short(quality));
            if (fields.size() == 4) {
                read.surjection_target = fields[3];
            }
            batch.emplace_back(move(read));
            if (batch.size() >= batch_size) {
                flush_batch();
            }
        }
        flush_batch();
        // Answer every surjection still in flight before exiting.
        if (surject_pool) surject_pool->drain();
    } catch (const exception& e) {
        cerr << "error [vg giraffe-server]: " << e.what() << endl;
        return 1;
    }

    return 0;
}

} // namespace

static Subcommand vg_giraffe_server("giraffe-server", "server-style in-process giraffe mapper", TOOLKIT, main_giraffe_server);
