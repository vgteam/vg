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
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "subcommand.hpp"
#include "giraffe_options.hpp"

#include "../giraffe_engine.hpp"
#include <vg/io/alignment_io.hpp>
#include <vg/io/gafkluge.hpp>
#include "../utility.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

namespace {

// Codes for the long options with no short form. They start above 255, so
// they cannot clash with a short option character.
constexpr int OPT_EMIT_HEADER   = 1000;
constexpr int OPT_SURJECT_THREADS = 1001;
constexpr int OPT_PARAMETER_PRESET = 1002;
constexpr int OPT_SET = 1003;
constexpr int OPT_NO_REC_MODE = 1004;

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

/// Writes the server's answers as frames. A frame is written whole, so frames
/// from different threads never mix.
class FrameWriter {
public:
    explicit FrameWriter(ostream& out) : out(out) {}

    /// Write READ<TAB>NAME<TAB>N and then the read's N GAF lines.
    void write_result(const string& name, const vector<string>& gaf_lines) {
        lock_guard<mutex> lock(mu);
        out << "READ\t" << name << '\t' << gaf_lines.size() << '\n';
        for (const string& gaf_line : gaf_lines) {
            out << gaf_line << '\n';
        }
        out.flush();
    }

    /// Write ERROR<TAB>NAME<TAB>MESSAGE. A frame is one line, so tabs and
    /// newlines in the message become spaces.
    void write_error(const string& name, string message) {
        for (char& c : message) {
            if (c == '\t' || c == '\n' || c == '\r') {
                c = ' ';
            }
        }
        lock_guard<mutex> lock(mu);
        out << "ERROR\t" << name << '\t' << message << '\n';
        out.flush();
    }

private:
    ostream& out;
    mutex mu;
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
         << "      --surject-threads N      run SURJECT_WITH_ANCHORS requests on N worker" << endl
         << "                               threads (0 = inline, one at a time) [--threads]" << endl
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
         << "    Response: READ<TAB>NAME<TAB>1 then one GAF line" << endl
         << "    carrying the surjection tags (sj:Z:..., sn:Z:..., etc.)." << endl
         << endl
         << "stdout:" << endl
         << "  GAF header lines first, with --emit-header. Then one frame per request:" << endl
         << "    READ<TAB>NAME<TAB>N, then N GAF lines, one per reported alignment (an" << endl
         << "      unmapped read gets one unmapped line)" << endl
         << "    ERROR<TAB>NAME<TAB>MESSAGE for a request that failed. NAME is * for a" << endl
         << "      line that is not a request." << endl
         << "  Reads are answered in input order. The answer to a SURJECT_WITH_ANCHORS" << endl
         << "  request can come after the answers to later requests." << endl
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

/// Append the GAF tags for a haplotype assignment.
void append_haplotype_tags(stringstream& ss, const HaplotypeAssignmentResult& hap) {
    // hp:Z lists the distinct haplotype names and hn:i counts them. hb:Z is
    // the representative, and hf:i is 1 if it is a reference. hc:i counts the
    // candidate sequences before the cap.
    ss << "\thp:Z:" << join_delim(hap.haplotype_names, ',')
       << "\thn:i:" << hap.haplotype_names.size()
       << "\thb:Z:" << (hap.representative_name.empty()
                        ? string("*") : hap.representative_name)
       << "\thf:i:" << (hap.representative_is_reference ? 1 : 0)
       << "\thc:i:" << hap.candidate_count;
    // hl:i and hr:i are the enclosing snarl's boundary nodes, if one was
    // found.
    if (hap.left_boundary_node != 0) {
        ss << "\thl:i:" << hap.left_boundary_node;
    }
    if (hap.right_boundary_node != 0) {
        ss << "\thr:i:" << hap.right_boundary_node;
    }
    // ht:i is 1 if the cap cut the list. hq:i is 1 if one haplotype carries
    // the whole read. hv:i is the percentage of the read's path that the
    // longest segment covers, and hs:i counts the segments.
    ss << "\tht:i:" << (hap.seq_ids.size() < hap.candidate_count ? 1 : 0)
       << "\thq:i:" << (hap.fully_covered ? 1 : 0)
       << "\thv:i:" << (hap.path_length_bp
                        ? (100 * hap.best_covered_bp / hap.path_length_bp)
                        : 0)
       << "\ths:i:" << hap.segments.size();
    // hm:Z lists the segments in path order, if there are several. Each
    // entry is covered_bp:haplotype_count:representative, and ';' separates
    // them.
    if (hap.segments.size() > 1) {
        ss << "\thm:Z:";
        for (size_t s = 0; s < hap.segments.size(); ++s) {
            const auto& seg = hap.segments[s];
            ss << (s ? ";" : "")
               << seg.covered_bp << ":"
               << seg.haplotype_names.size() << ":"
               << (seg.representative_name.empty()
                   ? string("*") : seg.representative_name);
        }
    }
}

/// Return the token that stands for a status in the sj:Z tag.
const char* surjection_status_token(HaplotypeSurjectionResult::Status status) {
    using S = HaplotypeSurjectionResult::Status;
    switch (status) {
        case S::OK:                return "ok";
        case S::EMPTY_INPUT:       return "empty";
        case S::UNKNOWN_PATH:      return "unknown_path";
        case S::PATH_NOT_INDEXED:  return "not_indexed";
        case S::INCOMPATIBLE:      return "incompatible";
        case S::SURJECTION_FAILED: return "failed";
    }
    return "unknown";
}

/// Append the GAF tags for a surjection result.
void append_surjection_tags(stringstream& ss, const HaplotypeSurjectionResult& r) {
    // sj:Z always holds the status. Only a successful result has the others.
    ss << "\tsj:Z:" << surjection_status_token(r.status);
    if (r.status != HaplotypeSurjectionResult::Status::OK) {
        return;
    }
    // sn:Z is the path name, sp:i the start, sr:i 1 if on the reverse strand,
    // ss:i the score, sm:i the MAPQ and sc:Z the CIGAR. sM:i, sX:i, sI:i and
    // sD:i count the matched, mismatched, inserted and deleted bases.
    ss << "\tsn:Z:" << r.path_name
       << "\tsp:i:" << r.path_position
       << "\tsr:i:" << (r.path_reverse ? 1 : 0)
       << "\tss:i:" << r.score
       << "\tsm:i:" << r.mapping_quality
       << "\tsc:Z:" << (r.cigar.empty() ? string("*") : r.cigar)
       << "\tsM:i:" << r.matches
       << "\tsX:i:" << r.mismatches
       << "\tsI:i:" << r.inserted
       << "\tsD:i:" << r.deleted;
}

/// Return the GAF line of a mapped alignment, with its haplotype tags and
/// surjection tags.
string mapped_alignment_to_gaf(const HandleGraph& graph, const MappedAlignment& mapped) {
    stringstream ss;
    ss << io::alignment_to_gaf(graph, mapped.alignment);
    if (mapped.haplotypes) {
        append_haplotype_tags(ss, *mapped.haplotypes);
    }
    if (mapped.surjection) {
        append_surjection_tags(ss, *mapped.surjection);
    }
    return ss.str();
}

/// Surject the alignment of a GAF line with anchors. Return the line with the
/// surjection tags. Throws if the line is empty or is not GAF, or if the
/// surjector throws.
string surject_gaf_line(GiraffeEngine& engine, const string& read_name, const string& gaf_line,
                        const vector<PathAnchor>& anchors, const string& target,
                        size_t target_path_length) {
    if (gaf_line.empty()) {
        throw runtime_error("surject_with_anchors: empty GAF input");
    }
    gafkluge::GafRecord gaf_record;
    gafkluge::parse_gaf_record(gaf_line, gaf_record);
    Alignment aln;
    io::gaf_to_alignment(engine.gbz().graph, gaf_record, aln);
    if (!read_name.empty()) {
        aln.set_name(read_name);
    }

    using S = HaplotypeSurjectionResult::Status;
    stringstream ss;
    ss << gaf_line;
    HaplotypeSurjectionResult result =
        engine.surject_with_anchors(aln, anchors, target, target_path_length);
    append_surjection_tags(ss, result);
    if (result.status == S::OK) {
        // an:i counts the anchors and ap:Z is the target path's name.
        ss << "\tan:i:" << anchors.size()
           << "\tap:Z:" << result.path_name;
    }
    return ss.str();
}

int main_giraffe_server(int argc, char** argv) {
    GiraffeEngineConfig config;
    size_t batch_size = 256;
    bool emit_header = false;
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

    // Fill vg giraffe's option parser now, so bad parameters fail before any index is loaded.
    std::unique_ptr<GroupedOptionGroup> parser = make_giraffe_parser(params);
    // this variable belongs to surject, so we handle it a bit differently
    config.prune_suspicious_anchors = parser->get_option_value<bool>("prune-low-cplx");
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
            gbwtgraph::GraphName graph_name = engine.gbz().graph_name();
            for (const auto& header : graph_name.gaf_header_lines()) {
                cout << header << '\n';
            }
            cout.flush();
        }

        // Everything after the header goes to stdout through the frame writer.
        FrameWriter frames(cout);

        // Surjections run on threads. Each frame names its read, so a client
        // can match up answers that come out of order.
        const size_t n_surject = (surject_threads < 0)
            ? static_cast<size_t>(config.threads) : static_cast<size_t>(surject_threads);
        unique_ptr<SurjectPool> surject_pool;
        if (n_surject > 0) {
            surject_pool = make_unique<SurjectPool>(n_surject);
        }

        vector<GiraffeServerRequest> batch;
        batch.reserve(batch_size);

        string line;
        size_t auto_id = 0;
        // Report a bad read on stderr, and answer it with an error frame.
        auto emit_read_error = [&](const string& name, const string& reason) {
            cerr << "warning [vg giraffe-server]: skipping read"
                 << (name.empty() ? "" : " '" + name + "'") << ": " << reason << endl;
            frames.write_error(name, reason);
        };

        auto flush_batch = [&]() {
            if (batch.empty()) {
                return;
            }
            vector<MappedRead> mapped_reads = engine.map_reads(batch);
            for (size_t i = 0; i < mapped_reads.size(); ++i) {
                const string& name = batch[i].read.name();
                if (mapped_reads[i].error) {
                    emit_read_error(name, *mapped_reads[i].error);
                    continue;
                }
                vector<string> gaf_lines;
                for (const MappedAlignment& mapped_aln : mapped_reads[i].alignments) {
                    gaf_lines.push_back(mapped_alignment_to_gaf(engine.gbz().graph, mapped_aln));
                }
                frames.write_result(name, gaf_lines);
            }
            batch.clear();
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
                auto header_fields = split_delims_keep_empty(line, "\t");
                // The header holds the command, read name, target path, path
                // length and anchor count.
                if (header_fields.size() != 5) {
                    string error = "malformed SURJECT_WITH_ANCHORS header (need 5 tab fields, got "
                                   + to_string(header_fields.size()) + ")";
                    cerr << "error [vg giraffe-server]: " << error << endl;
                    frames.write_error(header_fields.size() >= 2 ? header_fields[1] : "*", error);
                    continue;
                }
                const string& sa_name = header_fields[1];
                const string& sa_target = header_fields[2];
                optional<size_t> sa_path_len = parse_request_number(header_fields[3]);
                optional<size_t> sa_n_anchors = parse_request_number(header_fields[4]);
                if (!sa_n_anchors) {
                    // Without the anchor count we cannot skip the request's
                    // other lines. They are not reads, so each of them gets
                    // an error frame of its own.
                    string error = "malformed SURJECT_WITH_ANCHORS header (anchor count must be a number)";
                    cerr << "error [vg giraffe-server]: " << error << endl;
                    frames.write_error(sa_name, error);
                    continue;
                }

                // Read all of the request's lines even after an error, so the
                // next request starts on its own line.
                vector<PathAnchor> sa_anchors;
                string sa_error;
                if (!sa_path_len) {
                    sa_error = "malformed SURJECT_WITH_ANCHORS header (path length must be a number)";
                }
                bool sa_eof = false;
                for (size_t ai = 0; ai < *sa_n_anchors; ++ai) {
                    string anchor_line;
                    if (!getline(cin, anchor_line)) {
                        sa_error = "SURJECT_WITH_ANCHORS hit EOF while reading anchor "
                                   + to_string(ai) + "/" + to_string(*sa_n_anchors);
                        sa_eof = true;
                        break;
                    }
                    ++line_no;
                    if (!sa_error.empty()) {
                        continue;
                    }
                    auto af = split_delims_keep_empty(anchor_line, "\t");
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
                    // The last four numbers are the anchor's read and source
                    // mapping ranges, which the surjection does not use.
                    PathAnchor anchor;
                    anchor.step_begin_node        = nums[0];
                    anchor.step_begin_offset      = nums[1];
                    anchor.step_end_node          = nums[2];
                    anchor.step_end_offset        = nums[3];
                    anchor.path_offset_step_begin = nums[4];
                    anchor.path_offset_step_end   = nums[5];
                    sa_anchors.push_back(anchor);
                }

                string sa_gaf_line;
                if (!sa_eof) {
                    if (getline(cin, sa_gaf_line)) {
                        ++line_no;
                    } else {
                        sa_error = "SURJECT_WITH_ANCHORS hit EOF while reading the GAF line";
                    }
                }
                if (!sa_error.empty()) {
                    cerr << "error [vg giraffe-server]: " << sa_error << endl;
                    frames.write_error(sa_name, sa_error);
                    continue;
                }

                // Parsing stays on this thread, which reads stdin. The surjection
                // needs nothing from other requests, so it can run on the pool.
                auto task = [&engine, &frames,
                             name = sa_name, target = sa_target,
                             gaf = move(sa_gaf_line), anchors = move(sa_anchors),
                             path_len = *sa_path_len]() {
                    // Each stderr message is one string, so messages from two
                    // threads do not mix.
                    string gaf_line;
                    try {
                        gaf_line = surject_gaf_line(engine, name, gaf, anchors, target, path_len);
                    } catch (const exception& e) {
                        cerr << "error [vg giraffe-server]: surject_with_anchors failed for read '"
                                + name + "': " + e.what() + "\n";
                        frames.write_error(name, e.what());
                        return;
                    } catch (...) {
                        // An exception that escapes a pool thread terminates the
                        // server.
                        cerr << "error [vg giraffe-server]: surject_with_anchors failed for read '"
                                + name + "' (non-standard exception)\n";
                        frames.write_error(name, "non-standard exception");
                        return;
                    }
                    frames.write_result(name, {gaf_line});
                };
                if (surject_pool) {
                    surject_pool->submit(move(task));
                } else {
                    task();
                }
                continue;
            }
            auto fields = split_delims_keep_empty(line, "\t");
            // A read line has 1 to 4 tab-separated fields:
            //   1: SEQ
            //   2: NAME, SEQ
            //   3: NAME, SEQ, QUAL
            //   4: NAME, SEQ, QUAL, SURJ_TARGET (QUAL may be empty)
            if (fields.size() < 1 || fields.size() > 4) {
                string reason = "expected 1-4 tab-separated fields or a known command, got "
                                + to_string(fields.size());
                cerr << "warning [vg giraffe-server]: line " << line_no
                     << ": ignored (" << reason << ")" << endl;
                frames.write_error("*", "line " + to_string(line_no) + ": " + reason);
                continue;
            }

            // A lone field must be a sequence. Anything else is an unknown command.
            if (fields.size() == 1 && !is_nucleotide_sequence(fields[0])) {
                string reason = "not a known command or a nucleotide sequence";
                cerr << "warning [vg giraffe-server]: line " << line_no
                     << ": ignored (" << reason << ")" << endl;
                frames.write_error("*", "line " + to_string(line_no) + ": " + reason);
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
