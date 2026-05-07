/**
 * \file giraffe_server_main.cpp
 * Server-style in-process Giraffe mapping:
 * - load indexes once
 * - map batches in parallel
 * - stream GAF output to stdout
 */

#include <getopt.h>
#include <unistd.h>

#include <iostream>
#include <stdexcept>
#include <string>
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
         << "  -M, --max-multimaps N        max mappings per read [1]" << endl
         << "  -b, --batch-size N           input reads per mapping batch [256]" << endl
         << "      --emit-header            emit GAF header lines before output" << endl
         << "      --framed-output          emit read-grouped framed output for middleware" << endl
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

        vector<GiraffeFastqRead> batch;
        batch.reserve(batch_size);

        string line;
        size_t auto_id = 0;
        auto flush_batch = [&]() {
            if (batch.empty()) {
                return;
            }
            auto mapped = engine.map_reads(batch);
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
    } catch (const exception& e) {
        cerr << "error [vg giraffe-server]: " << e.what() << endl;
        return 1;
    }

    return 0;
}

} // namespace

static Subcommand vg_giraffe_server("giraffe-server", "server-style in-process giraffe mapper", TOOLKIT, main_giraffe_server);
