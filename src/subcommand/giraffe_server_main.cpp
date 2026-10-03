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
#include <deque>
#include <functional>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "subcommand.hpp"

#include "../giraffe_engine.hpp"
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
         << "  -M, --max-multimaps N        max mappings per read (0 = BLAT default 100)" << endl
         << "  -b, --batch-size N           input reads per mapping batch [256]" << endl
         << "      --emit-header            emit GAF header lines before output" << endl
         << "      --framed-output          emit read-grouped framed output for middleware" << endl
         << "      --surject-threads N      run SURJECT_WITH_ANCHORS requests on N worker" << endl
         << "                               threads (framed output only; 0 = inline, one at" << endl
         << "                               a time) [--threads]" << endl
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
         << "  QUALITY is FASTQ phred+33 ASCII; if non-empty, must match SEQUENCE length" << endl
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
         << endl;
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
    GiraffeEnginePaths paths;
    GiraffeEngineConfig config;
    size_t batch_size = 256;
    bool emit_header = false;
    bool framed_output = false;
    long surject_threads = -1;   // -1 = same as --threads

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
            {"surject-target", required_argument, 0, 'S'},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "Z:m:d:z:t:M:b:S:h?", long_options, &option_index);
        if (c == -1) {
            break;
        }

        switch (c) {
            case 'Z':
                paths.gbz_path = optarg;
                break;
            case 'm':
                paths.minimizer_path = optarg;
                break;
            case 'd':
                paths.distance_path = optarg;
                break;
            case 'z':
                paths.zipcode_path = optarg;
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
            case 'S':
                config.surjection_target_paths.emplace_back(optarg);
                break;
            case 'h':
            case '?':
            default:
                help_giraffe_server(argv);
                return 1;
        }
    }

    if (argc > optind || paths.gbz_path.empty() || paths.minimizer_path.empty()
        || paths.distance_path.empty() || paths.zipcode_path.empty() || batch_size == 0) {
        help_giraffe_server(argv);
        return 1;
    }

    try {
        GiraffeEngine engine;
        engine.load(paths, config);

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

        vector<GiraffeFastqRead> batch;
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
                    cout << "READ\t" << batch[i].name << '\t' << read_mappings.size() << '\n';
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

            GiraffeFastqRead read;
            if (fields.size() == 1) {
                read.name = "read_" + to_string(auto_id++);
                read.sequence = fields[0];
            } else if (fields.size() == 2) {
                read.name = fields[0];
                read.sequence = fields[1];
            } else if (fields.size() == 3) {
                read.name = fields[0];
                read.sequence = fields[1];
                read.quality = fields[2];
            } else {  // fields.size() == 4
                read.name = fields[0];
                read.sequence = fields[1];
                read.quality = fields[2];
                read.surjection_target = fields[3];
            }

            // Validate: empty sequence is malformed; respond rather than dropping silently.
            if (read.sequence.empty()) {
                emit_read_error(read.name, "empty SEQUENCE");
                continue;
            }
            // Validate: if QUALITY is supplied, it must match SEQUENCE length (phred+33 ASCII).
            if (!read.quality.empty() && read.quality.size() != read.sequence.size()) {
                emit_read_error(read.name,
                    "QUALITY length (" + to_string(read.quality.size())
                    + ") does not match SEQUENCE length ("
                    + to_string(read.sequence.size()) + ")");
                continue;
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
