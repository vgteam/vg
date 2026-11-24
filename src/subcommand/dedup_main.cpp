/** \file dedup_main.cpp
 *
 * Defines the "vg dedup" subcommand, which uses Content-Defined Chunking for
 * finding good chopping points in long nodes to enable sequence deduplication
 * in GBZ.
 */

#include "subcommand.hpp"

#include "../dedup.hpp"
#include "../gfa.hpp"
#include "../handle.hpp"
#include "../log.hpp"
#include "../utility.hpp"
#include "../vg.hpp"
#include "../io/save_handle_graph.hpp"

#include <vg/io/vpkg.hpp>
#include <gbwt/utils.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <string>

#include <getopt.h>
#include <omp.h>

using namespace vg;
using namespace vg::subcommand;

//----------------------------------------------------------------------------

namespace {

const std::string DEFAULT_NAME_PREFIX = "superstring_";

void help_dedup(char** argv) {
    ContentChopParameters defaults;
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " [options] graph" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Content-aware node chopping to enable sequence deduplication in GBZ." << std::endl;
    std::cerr << "For best results, unchop the graph first." << std::endl;
    std::cerr << std::endl;

    std::cerr << "Output options:" << std::endl;
    std::cerr << "  -o, --output FILE       write the chopped graph to FILE in the input format" << std::endl;
    std::cerr << "  -f, --fasta FILE        write the stored strings to FILE in FASTA format" << std::endl;
    std::cerr << "  -G, --gfa FILE          write the chopped graph to FILE in GFA format, with" << std::endl;
    std::cerr << "                          the sequences of chopped nodes stored externally" << std::endl;
    std::cerr << "  -p, --prefix STR        name stored string i STR followed by i [" << DEFAULT_NAME_PREFIX << "]" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Chopping options:" << std::endl;
    std::cerr << "  -m, --max-length N      chop nodes over N bp into at most N bp pieces [" << defaults.max_node_length << "]" << std::endl;
    std::cerr << "  -w, --piece-weight X    cost of each piece in bp of stored sequence [" << defaults.piece_weight << "]" << std::endl;
    std::cerr << "  -k, --kmer-length N     k-mers length for picking candidate cuts [" << defaults.kmer_length << "]" << std::endl;
    std::cerr << "  -r, --anchor-radius N   a candidate cut is the minimum k-mer within N bp [" << defaults.anchor_radius << "]" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Other options:" << std::endl;
    std::cerr << "  -t, --threads N         use N threads [" << omp_get_max_threads() << "]" << std::endl;
    std::cerr << "  -v, --verbose           print statistics, running time, and memory usage" << std::endl;
    std::cerr << "  -h, --help              print this help message to stderr and exit" << std::endl;
    std::cerr << std::endl;
}

struct DedupConfig {
    std::string graph_file;

    // Output files.
    std::string graph_output;
    std::string fasta_output;
    std::string gfa_output;

    /// Name prefix for stored strings.
    std::string name_prefix = DEFAULT_NAME_PREFIX;

    ContentChopParameters parameters;
    bool verbose = false;

    DedupConfig(int argc, char** argv, Logger& logger);
};

DedupConfig::DedupConfig(int argc, char** argv, Logger& logger) {
    if (argc < 3) {
        help_dedup(argv);
        std::exit(EXIT_FAILURE);
    }

    static struct option long_options[] = {
        { "output", required_argument, 0, 'o' },
        { "fasta", required_argument, 0, 'f' },
        { "gfa", required_argument, 0, 'G' },
        { "prefix", required_argument, 0, 'p' },
        { "max-length", required_argument, 0, 'm' },
        { "piece-weight", required_argument, 0, 'w' },
        { "kmer-length", required_argument, 0, 'k' },
        { "anchor-radius", required_argument, 0, 'r' },
        { "threads", required_argument, 0, 't' },
        { "verbose", no_argument, 0, 'v' },
        { "help", no_argument, 0, 'h' },
        { 0, 0, 0, 0 }
    };
    optind = 2; // force optind past command positional argument
    while (true) {
        int option_index = 0;
        int c = getopt_long(argc, argv, "o:f:G:p:m:w:k:r:t:vh?", long_options, &option_index);
        if (c == -1) { break; } // End of options.
        switch (c)
        {
        case 'o':
            this->graph_output = ensure_writable(logger, optarg);
            break;
        case 'f':
            this->fasta_output = ensure_writable(logger, optarg);
            break;
        case 'G':
            this->gfa_output = ensure_writable(logger, optarg);
            break;
        case 'p':
            this->name_prefix = optarg;
            break;
        case 'm':
            this->parameters.max_node_length = parse<size_t>(optarg);
            break;
        case 'w':
            this->parameters.piece_weight = parse<double>(optarg);
            break;
        case 'k':
            this->parameters.kmer_length = parse<size_t>(optarg);
            break;
        case 'r':
            this->parameters.anchor_radius = parse<size_t>(optarg);
            break;
        case 't':
            set_thread_count(logger, optarg);
            break;
        case 'v':
            this->verbose = true;
            break;
        case 'h':
        case '?':
            help_dedup(argv);
            std::exit(EXIT_FAILURE);
        default:
            std::abort();
        }
    }

    if (optind + 1 != argc) {
        help_dedup(argv);
        std::exit(EXIT_FAILURE);
    }
    this->graph_file = require_exists(logger, argv[optind]);

    if (this->graph_output.empty() && this->fasta_output.empty() && this->gfa_output.empty() && !this->verbose) {
        logger.error() << "no output requested; use --output, --fasta, --gfa, and/or --verbose" << std::endl;
    }

    try {
        this->parameters.validate();
    } catch (const std::invalid_argument& e) {
        logger.error() << e.what() << std::endl;
    }
}

//----------------------------------------------------------------------------

struct ChoppingSummary {
    // Original graph.
    size_t input_nodes = 0;
    size_t input_length = 0;
    size_t max_node_length = 0;

    // Nodes that will be chopped in the original graph.
    size_t long_nodes = 0;
    size_t long_length = 0;

    // Chopped nodes.
    size_t chopped_length = 0;
    size_t pieces = 0;
    size_t max_piece_length = 0;

    // Stored superstrings.
    size_t stored_strings = 0;
    size_t stored_length = 0;
    size_t max_string_length = 0;

    ChoppingSummary(const HandleGraph& graph, const ChoppingPlan& plan, const ContentChopParameters& parameters);

    void report(Logger& logger) const;
};

ChoppingSummary::ChoppingSummary(const HandleGraph& graph, const ChoppingPlan& plan, const ContentChopParameters& parameters) {
    graph.for_each_handle([&](const handle_t& handle) {
        size_t length = graph.get_length(handle);
        this->input_nodes++;
        this->input_length += length;
        this->max_node_length = std::max(this->max_node_length, length);
        if (length > parameters.max_node_length) {
            this->long_nodes++;
            this->long_length += length;
        }
    });

    this->pieces = plan.pieces.size();
    for (const PieceLocation& piece : plan.pieces) {
        this->chopped_length += piece.length;
        this->max_piece_length = std::max(this->max_piece_length, piece.length);
    }

    this->stored_strings = plan.strings.size();
    for (const std::string& s : plan.strings) {
        this->stored_length += s.size();
        this->max_string_length = std::max(this->max_string_length, s.size());
    }
}

double ratio(size_t numerator, size_t denominator) {
    return (denominator > 0 ? static_cast<double>(numerator) / denominator : 0.0);
}

void ChoppingSummary::report(Logger& logger) const {
    logger.info() << "Input graph: " << this->input_nodes << " nodes, " << this->input_length << " bp" << std::endl;
    logger.info() << "Long nodes: " << this->long_nodes << " nodes, " << this->long_length << " bp" << std::endl;
    logger.info() << "Chopped nodes: " << this->pieces << " pieces (mean " << ratio(this->chopped_length, this->pieces)
                  << " bp, max " << this->max_piece_length << " bp)" << std::endl;
    logger.info() << "Stored strings: " << this->stored_strings << " strings, " << this->stored_length << " bp ("
                  << ratio(this->stored_length, this->chopped_length) << "x; max "
                  << this->max_string_length << " bp)" << std::endl;

    // Unchopped nodes keep their sequences, and the pieces share the stored strings.
    size_t output_nodes = this->input_nodes - this->long_nodes + this->pieces;
    size_t output_length = this->input_length - this->chopped_length + this->stored_length;
    logger.info() << "Output graph: " << output_nodes << " nodes (" << ratio(output_nodes, this->input_nodes)
                  << "x), " << output_length << " bp ("
                  << ratio(output_length, this->input_length) << "x)" << std::endl;
}

//----------------------------------------------------------------------------

std::string stored_string_name(const std::string& prefix, size_t string_id) {
    return prefix + std::to_string(string_id);
}

void write_stored_strings_fasta(const ChoppingPlan& plan, const std::string& name_prefix, std::ostream& out) {
    constexpr size_t LINE_WIDTH = 80;
    for (size_t string_id = 0; string_id < plan.strings.size(); string_id++) {
        out << ">" << stored_string_name(name_prefix, string_id) << "\n";
        const std::string& sequence = plan.strings[string_id];
        for (size_t offset = 0; offset < sequence.size(); offset += LINE_WIDTH) {
            size_t length = std::min(LINE_WIDTH, sequence.size() - offset);
            out.write(sequence.data() + offset, length);
            out << "\n";
        }
    }
}

// Custom version of write_gfa_nodes().
void write_gfa_piece_nodes(
    const handlegraph::HandleGraph& graph, const PieceDirectory& directory,
    const std::string& name_prefix, std::ostream& out
) {
    graph.for_each_handle([&](const handle_t& handle) {
        nid_t node_id = graph.get_id(handle);
        out << "S\t" << node_id << "\t";
        const PieceLocation* piece = directory.find(node_id);
        if (piece == nullptr) {
            out << graph.get_sequence(handle);
        } else {
            out << "*"
                << "\tsn:Z:" << stored_string_name(name_prefix, piece->string_id)
                << "\tso:i:" << piece->string_offset
                << "\tLN:i:" << piece->length;
        }
        out << "\n";
    });
}

// Custom version of graph_to_gfa() that uses write_gfa_piece_nodes()
// instead of write_gfa_nodes().
void write_chopped_gfa(
    const handlegraph::PathHandleGraph& graph, const PieceDirectory& directory,
    const std::string& name_prefix, std::ostream& out
) {
    write_gfa_header(&graph, out);
    write_gfa_piece_nodes(graph, directory, name_prefix, out);
    write_gfa_edges(&graph, out);
    write_gfa_paths(&graph, out);
}

void write_to_file(const std::string& filename, const std::function<void(std::ostream&)>& writer, Logger& logger) {
    std::ofstream out(filename, std::ios_base::binary);
    if (!out) {
        logger.error() << "cannot open " << filename << " for writing" << std::endl;
    }
    writer(out);
    out.close();
    if (!out) {
        logger.error() << "cannot write to " << filename << std::endl;
    }
}

} // anonymous namespace

//----------------------------------------------------------------------------

int main_dedup(int argc, char** argv) {
    Logger logger("vg dedup");
    DedupConfig config(argc, argv, logger);
    double start = gbwt::readTimer();

    double checkpoint = start;
    if (config.verbose) {
        logger.info() << "Loading input graph " << config.graph_file << std::endl;
    }
    std::unique_ptr<MutablePathDeletableHandleGraph> graph =
        vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(config.graph_file);
    if (config.verbose) {
        logger.info() << "Loaded the graph in " << (gbwt::readTimer() - checkpoint) << " seconds" << std::endl;
    }

    checkpoint = gbwt::readTimer();
    if (config.verbose) {
        logger.info() << "Making a chopping plan with max node length " << config.parameters.max_node_length << std::endl;
    }
    ChoppingPlan plan = plan_content_chopping(*graph, config.parameters);
    if (config.verbose) {
        logger.info() << "Planned the chopping in " << (gbwt::readTimer() - checkpoint) << " seconds with "
                      << omp_get_max_threads() << " threads" << std::endl;
        ChoppingSummary summary(*graph, plan, config.parameters);
        summary.report(logger);
    }

    if (!config.fasta_output.empty()) {
        checkpoint = gbwt::readTimer();
        if (config.verbose) {
            logger.info() << "Writing stored strings to " << config.fasta_output << std::endl;
        }
        write_to_file(config.fasta_output, [&](std::ostream& out) {
            write_stored_strings_fasta(plan, config.name_prefix, out);
        }, logger);
        if (config.verbose) {
            logger.info() << "Wrote the stored strings in " << (gbwt::readTimer() - checkpoint) << " seconds" << std::endl;
        }
    }

    checkpoint = gbwt::readTimer();
    if (config.verbose) {
        // TODO: Memory bottleneck appears to be plan_content_chopping().
        logger.info() << "Memory usage so far: " << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
        logger.info() << "Chopping the graph" << std::endl;
    }
    PieceDirectory directory;
    try {
        directory = chop_by_plan(*graph, plan);
    } catch (const std::runtime_error& e) {
        logger.error() << e.what() << std::endl;
    }
    if (VG* vg_graph = dynamic_cast<VG*>(graph.get())) {
        // Chopping leaves gaps in the ranks of path steps in a VG graph.
        vg_graph->paths.compact_ranks();
    }
    if (config.verbose) {
        logger.info() << "Chopped the graph in " << (gbwt::readTimer() - checkpoint) << " seconds" << std::endl;
    }

    if (!config.graph_output.empty()) {
        checkpoint = gbwt::readTimer();
        if (config.verbose) {
            logger.info() << "Writing the chopped graph to " << config.graph_output << std::endl;
        }
        vg::io::save_handle_graph(graph.get(), config.graph_output);
        if (config.verbose) {
            logger.info() << "Wrote the chopped graph in " << (gbwt::readTimer() - checkpoint) << " seconds" << std::endl;
        }
    }

    if (!config.gfa_output.empty()) {
        checkpoint = gbwt::readTimer();
        if (config.verbose) {
            logger.info() << "Writing the GFA to " << config.gfa_output << std::endl;
        }
        write_to_file(config.gfa_output, [&](std::ostream& out) {
            write_chopped_gfa(*graph, directory, config.name_prefix, out);
        }, logger);
        if (config.verbose) {
            logger.info() << "Wrote the GFA in " << (gbwt::readTimer() - checkpoint) << " seconds" << std::endl;
        }
    }

    if (config.verbose) {
        logger.info() << "Used " << (gbwt::readTimer() - start) << " seconds, "
                      << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }
    return 0;
}

//----------------------------------------------------------------------------

// Register subcommand
static Subcommand vg_dedup("dedup", "content-aware node chopping to enable sequence deduplication", WIDGET, main_dedup);
