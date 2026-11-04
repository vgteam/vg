/**
 * \file giraffe_server_main.cpp
 * Server-style in-process Giraffe mapping:
 * - load indexes once
 * - map batches in parallel
 * - stream GAF output to stdout
 */

#include <getopt.h>
#include <unistd.h>

#include <condition_variable>
#include <cstdlib>
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

// Named option codes for long-only options (>255 to avoid clashing with short-option chars).
// Keep these in sync with the long_options table and the switch in main_giraffe_server.
constexpr int OPT_EMIT_HEADER   = 1000;
constexpr int OPT_FRAMED_OUTPUT = 1001;
constexpr int OPT_SURJECT_THREADS = 1002;
constexpr int OPT_PARAMETER_PRESET = 1003;
constexpr int OPT_SET = 1004;
constexpr int OPT_NO_REC_MODE = 1005;

/// Fixed-size worker pool for SURJECT_WITH_ANCHORS requests, so a batch of
/// targets surjects concurrently instead of one at a time on the stdin thread.
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

    /// Finish every queued task, then stop the workers.
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
                if (tasks_.empty()) return;          // stopping and drained
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
         << "  -M, --max-multimaps N        max mappings per read (0 = the preset's value;" << endl
         << "                               reporting is capped by max-alignments)" << endl
         << "  -b, --batch-size N           input reads per mapping batch [256]" << endl
         << "      --emit-header            emit GAF header lines before output" << endl
         << "      --framed-output          emit read-grouped framed output for middleware" << endl
         << "      --surject-threads N      run SURJECT_WITH_ANCHORS requests on N worker" << endl
         << "                               threads (framed output only; 0 = inline, one at" << endl
         << "                               a time) [--threads]" << endl
         << "      --parameter-preset NAME  vg giraffe parameter preset [hifi]" << endl
         << "      --set NAME=VALUE         set any vg giraffe option by its long name," << endl
         << "                               e.g. mapq-score-window=150 (may repeat)" << endl
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

/// True if seq is non-empty and only IUPAC nucleotide codes, in either case.
bool is_nucleotide_sequence(const string& seq) {
    static const string iupac = "ACGTUNRYKMSWBDHVacgtunrykmswbdhv";
    return !seq.empty() && seq.find_first_not_of(iupac) == string::npos;
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

int main_giraffe_server(int argc, char** argv) {
    GiraffeEngineConfig config;
    size_t batch_size = 256;
    bool emit_header = false;
    bool framed_output = false;
    long surject_threads = -1;   // -1 = same as --threads
    std::string preset;          // empty = $PANGENOME_GIRAFFE_PRESET, else "hifi"
    std::vector<std::pair<std::string, std::string>> overrides;
    std::optional<bool> rec_mode;   // unset = vg giraffe's default

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
                config.max_multimaps = parse<size_t>(optarg);
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
                preset = optarg;
                break;
            case OPT_SET: {
                std::string kv = optarg;
                size_t eq = kv.find('=');
                if (eq == std::string::npos || eq == 0) {
                    cerr << "error [vg giraffe-server]: --set expects NAME=VALUE, got " << kv << endl;
                    return 1;
                }
                overrides.emplace_back(kv.substr(0, eq), kv.substr(eq + 1));
                break;
            }
            case 'S':
                config.surjection_target_paths.emplace_back(optarg);
                break;
            case 'E':
                rec_mode = true;
                break;
            case OPT_NO_REC_MODE:
                rec_mode = false;
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

    // Mapping parameters come from vg giraffe's own presets and options. The
    // PANGENOME_* variables still work for existing deployments; --set wins.
    if (preset.empty()) {
        const char* e = std::getenv("PANGENOME_GIRAFFE_PRESET");
        preset = e ? e : "hifi";
    }
    std::vector<std::pair<std::string, std::string>> env_overrides;
    if (const char* e = std::getenv("PANGENOME_MAPQ_SCORE_SCALE")) {
        env_overrides.emplace_back("mapq-score-scale", e);
    }
    if (const char* e = std::getenv("PANGENOME_MAPQ_SCORE_WINDOW")) {
        env_overrides.emplace_back("mapq-score-window", e);
    }
    overrides.insert(overrides.begin(), env_overrides.begin(), env_overrides.end());
    const bool map_long_reads = giraffe_preset_maps_long_reads(preset);
    try {
        check_giraffe_parameters(preset, overrides);   // fail before loading indexes
        if (rec_mode.value_or(false) && !map_long_reads) {
            throw std::invalid_argument("Recombination-aware mode cannot be used with short reads yet.");
        }
    } catch (const std::invalid_argument& e) {
        cerr << "error [vg giraffe-server]: " << e.what() << endl;
        return 1;
    }
    config.configure_mapper = [preset, overrides, map_long_reads, rec_mode](
            MinimizerMapper& mapper, const gbwtgraph::DefaultMinimizerIndex& minimizer_index) {
        configure_giraffe_mapper(mapper, preset, overrides);
        bool use_rec_mode = giraffe_rec_mode(minimizer_index, map_long_reads, rec_mode);
        mapper.use_payload_paths &= use_rec_mode;
        cerr << "[vg giraffe-server] recombination-aware mode: " << (use_rec_mode ? "on" : "off") << endl;
    };
    cerr << "[vg giraffe-server] parameter preset: " << preset;
    for (const auto& kv : overrides) {
        cerr << " --" << kv.first << " " << kv.second;
    }
    cerr << endl;

    try {
        GiraffeEngine engine;
        engine.load(config);

        if (emit_header) {
            for (const auto& header : engine.gaf_header_lines()) {
                cout << header << '\n';
            }
            cout.flush();
        }

        // Every write to stdout goes through out_mu, so a frame (READ header plus
        // its lines) is never interleaved with another worker's.
        mutex out_mu;

        // Requests are demultiplexed by NAME, so surjections may finish out of
        // order -- but only framed output names its records, so concurrency is
        // limited to that mode.
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

        // Helper: report a per-read error so framed-output clients still get a record
        // for the read and don't hang waiting on it.
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
            // Batch-processing command (PROCESS_BATCH; FLUSH_NOW kept as deprecated alias).
            if (line == "PROCESS_BATCH" || line == "FLUSH_NOW") {
                flush_batch();
                continue;
            }
            // Anchored-surjection command. Multi-line: header + N anchor lines + 1 GAF line.
            // Flushes any pending batch first so framed output stays in order.
            if (line.compare(0, sizeof("SURJECT_WITH_ANCHORS\t") - 1,
                             "SURJECT_WITH_ANCHORS\t") == 0) {
                flush_batch();
                auto header_fields = split_tab_fields(line);
                // Expected: ["SURJECT_WITH_ANCHORS", name, target, path_len, n_anchors]
                if (header_fields.size() != 5) {
                    cerr << "error [vg giraffe-server]: malformed SURJECT_WITH_ANCHORS header "
                            "(need 5 tab fields, got " << header_fields.size() << ")" << endl;
                    if (header_fields.size() >= 2) emit_empty_frame(header_fields[1]);
                    continue;
                }
                const string& sa_name = header_fields[1];
                const string& sa_target = header_fields[2];
                size_t sa_path_len = parse<size_t>(header_fields[3]);
                size_t sa_n_anchors = parse<size_t>(header_fields[4]);

                vector<WireAnchor> sa_anchors;
                sa_anchors.reserve(sa_n_anchors);
                bool sa_read_ok = true;
                for (size_t ai = 0; ai < sa_n_anchors; ++ai) {
                    string anchor_line;
                    if (!getline(cin, anchor_line)) {
                        cerr << "error [vg giraffe-server]: SURJECT_WITH_ANCHORS hit EOF "
                             << "while reading anchor " << ai << "/" << sa_n_anchors << endl;
                        sa_read_ok = false;
                        break;
                    }
                    ++line_no;
                    auto af = split_tab_fields(anchor_line);
                    if (af.size() != 10) {
                        cerr << "error [vg giraffe-server]: malformed anchor line " << line_no
                             << " (need 10 tab fields, got " << af.size() << ")" << endl;
                        sa_read_ok = false;
                        break;
                    }
                    WireAnchor wa;
                    wa.step_begin_node       = parse<uint64_t>(af[0]);
                    wa.step_begin_offset     = parse<uint64_t>(af[1]);
                    wa.step_end_node         = parse<uint64_t>(af[2]);
                    wa.step_end_offset       = parse<uint64_t>(af[3]);
                    wa.path_offset_step_begin = parse<size_t>(af[4]);
                    wa.path_offset_step_end   = parse<size_t>(af[5]);
                    wa.read_begin_offset      = parse<size_t>(af[6]);
                    wa.read_end_offset        = parse<size_t>(af[7]);
                    wa.source_mapping_begin   = parse<size_t>(af[8]);
                    wa.source_mapping_end     = parse<size_t>(af[9]);
                    sa_anchors.push_back(wa);
                }
                if (!sa_read_ok) {
                    emit_empty_frame(sa_name);
                    continue;
                }

                string sa_gaf_line;
                if (!getline(cin, sa_gaf_line)) {
                    cerr << "error [vg giraffe-server]: SURJECT_WITH_ANCHORS hit EOF "
                         << "while reading the GAF line" << endl;
                    emit_empty_frame(sa_name);
                    continue;
                }
                ++line_no;

                // Parsing above stays on this thread (it reads stdin); the
                // surjection itself is independent per request.
                auto task = [&engine, &out_mu, framed_output,
                             name = sa_name, target = sa_target,
                             gaf = move(sa_gaf_line), anchors = move(sa_anchors),
                             path_len = sa_path_len]() {
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
                        // A worker thread must never let an exception escape:
                        // that would terminate the whole server.
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
            // Field layouts (tab-separated):
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

            // A lone field must be a sequence; anything else is a bad command.
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
