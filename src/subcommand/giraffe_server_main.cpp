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
         << "      --emit-graph-path        after each GAF line, emit @GRAPH with node/rev/len triples" << endl
         << "  -h, --help                   show help" << endl
         << endl
         << "stdin protocol (one read per line):" << endl
         << "  SEQUENCE" << endl
         << "  NAME<TAB>SEQUENCE" << endl
         << "  NAME<TAB>SEQUENCE<TAB>QUALITY" << endl
         << "  @@FLUSH@@   (force immediate processing of buffered reads)" << endl
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
    bool emit_graph_path = false;

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
            {"emit-header", no_argument, 0, 1000},
            {"framed-output", no_argument, 0, 1001},
            {"emit-graph-path", no_argument, 0, 1002},
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "Z:m:d:z:t:M:b:h?", long_options, &option_index);
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
            case 1000:
                emit_header = true;
                break;
            case 1001:
                framed_output = true;
                break;
            case 1002:
                emit_graph_path = true;
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
            if (emit_graph_path) {
                auto detailed = engine.map_reads_detailed(batch);
                for (size_t i = 0; i < detailed.size(); ++i) {
                    const auto& read_mappings = detailed[i].alignments;
                    if (framed_output) {
                        cout << "@READ\t" << batch[i].name << '\t' << read_mappings.size() << '\n';
                    }
                    for (const auto& rec : read_mappings) {
                        cout << rec.gaf_line << '\n';
                        cout << "@GRAPH";
                        for (const auto& step : rec.graph_path) {
                            cout << '\t' << step.node_id << '\t' << (step.is_reverse ? 1 : 0) << '\t' << step.from_length;
                        }
                        cout << '\n';
                    }
                }
            } else {
                auto mapped = engine.map_reads(batch);
                for (size_t i = 0; i < mapped.size(); ++i) {
                    const auto& read_mappings = mapped[i];
                    if (framed_output) {
                        cout << "@READ\t" << batch[i].name << '\t' << read_mappings.size() << '\n';
                    }
                    for (const auto& gaf_line : read_mappings) {
                        cout << gaf_line << '\n';
                    }
                }
            }
            cout.flush();
            batch.clear();
        };

        while (getline(cin, line)) {
            if (line.empty()) {
                continue;
            }
            if (line == "@@FLUSH@@") {
                flush_batch();
                continue;
            }
            auto fields = split_tab_fields(line);
            GiraffeFastqRead read;
            if (fields.size() == 1) {
                read.name = "read_" + to_string(auto_id++);
                read.sequence = fields[0];
            } else if (fields.size() == 2) {
                read.name = fields[0];
                read.sequence = fields[1];
            } else {
                read.name = fields[0];
                read.sequence = fields[1];
                read.quality = fields[2];
            }

            if (read.sequence.empty()) {
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
