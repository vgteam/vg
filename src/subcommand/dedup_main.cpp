/** \file dedup_main.cpp
 *
 * Defines the "vg dedup" subcommand, which chops long nodes at places chosen
 * by their content, so that a GBZ graph can store repeated sequence once.
 */

#include "subcommand.hpp"

#include "../../artifacts/piece_directory.hpp"

#include "../dedup.hpp"
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

const Logger logger("vg dedup");

const std::string DEFAULT_NAME_PREFIX = "superstring_";

void help_dedup(char** argv) {
    ContentChopParameters defaults;
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " [options] graph" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Chop long nodes at places chosen by their sequence, so that repeated sequence" << std::endl;
    std::cerr << "can be stored once. Pieces of chopped nodes take their sequences from stored" << std::endl;
    std::cerr << "strings. For best results, unchop the graph first." << std::endl;
    std::cerr << std::endl;

    std::cerr << "Output options:" << std::endl;
    std::cerr << "  -o, --output FILE       write the chopped graph to FILE in the input format" << std::endl;
    std::cerr << "  -f, --fasta FILE        write the stored strings to FILE in FASTA format" << std::endl;
    std::cerr << "  -G, --gfa FILE          write the chopped graph to FILE in GFA format, with" << std::endl;
    std::cerr << "                          the sequences of pieces as locations in stored strings" << std::endl;
    std::cerr << "  -p, --prefix STR        name stored string i STR followed by i [" << DEFAULT_NAME_PREFIX << "]" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Chopping options:" << std::endl;
    std::cerr << "  -m, --max-length N      chop nodes over N bp into pieces of at most N bp" << std::endl;
    std::cerr << "                          [" << defaults.max_node_length << "]" << std::endl;
    std::cerr << "  -w, --piece-weight X    cost of each piece in bp of stored sequence [" << defaults.piece_weight << "]" << std::endl;
    std::cerr << "  -k, --kmer-length N     length of the k-mers that pick candidate cuts [" << defaults.kmer_length << "]" << std::endl;
    std::cerr << "  -r, --anchor-radius N   a candidate cut is the minimum k-mer within N bp [" << defaults.anchor_radius << "]" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Other options:" << std::endl;
    std::cerr << "  -t, --threads N         use N threads [" << omp_get_max_threads() << "]" << std::endl;
    std::cerr << "  -v, --verbose           print statistics, running time, and memory usage" << std::endl;
    std::cerr << "  -h, --help              print this help message to stderr and exit" << std::endl;
    std::cerr << std::endl;
}

/// The settings of one run of the subcommand, read from the command line.
struct DedupConfig {
    std::string graph_file;

    /// Names of the output files. We skip each output whose name is empty.
    std::string graph_output;
    std::string fasta_output;
    std::string gfa_output;

    /// Stored string `i` is named `name_prefix` followed by `i`.
    std::string name_prefix = DEFAULT_NAME_PREFIX;

    ContentChopParameters parameters;
    bool verbose = false;

    /// Read the settings from the command line. Exit with an error if they
    /// cannot be used.
    DedupConfig(int argc, char** argv);

    /// Return true if some output needs the chopped graph.
    bool needs_chopped_graph() const {
        return !this->graph_output.empty() || !this->gfa_output.empty();
    }
};

DedupConfig::DedupConfig(int argc, char** argv) {
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
        logger.error() << "no output requested; use --output, --fasta, --gfa, or --verbose" << std::endl;
    }

    try {
        this->parameters.validate();
    } catch (const std::invalid_argument& e) {
        logger.error() << e.what() << std::endl;
    }
}

//----------------------------------------------------------------------------

/// A record of what a chopping plan does to the nodes of a graph and to the
/// sequence that the graph stores.
struct ChoppingSummary {
    size_t input_nodes = 0;
    size_t input_length = 0;
    size_t max_node_length = 0;

    /// Nodes longer than the maximum length in the chopping parameters.
    size_t long_nodes = 0;
    size_t long_length = 0;

    size_t chopped_nodes = 0;
    size_t chopped_length = 0;
    size_t pieces = 0;
    size_t max_piece_length = 0;

    size_t stored_strings = 0;
    size_t stored_length = 0;
    size_t max_string_length = 0;

    /// Summarize the plan that was made for the graph with these parameters.
    ChoppingSummary(const HandleGraph& graph, const ChoppingPlan& plan, const ContentChopParameters& parameters);

    /// Log the summary as info messages.
    void report() const;
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

    this->chopped_nodes = plan.chop_points.size();
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

/// Return the ratio, or 0 if the denominator is 0.
double ratio(size_t numerator, size_t denominator) {
    return (denominator > 0 ? static_cast<double>(numerator) / denominator : 0.0);
}

void ChoppingSummary::report() const {
    logger.info() << "Input graph: " << this->input_nodes << " nodes, " << this->input_length << " bp"
                  << " (max node length " << this->max_node_length << " bp)" << std::endl;
    logger.info() << "Long nodes: " << this->long_nodes << " nodes, " << this->long_length << " bp" << std::endl;
    logger.info() << "Chopped nodes: " << this->chopped_nodes << " nodes, " << this->chopped_length << " bp, into "
                  << this->pieces << " pieces (mean " << ratio(this->chopped_length, this->pieces)
                  << " bp, max " << this->max_piece_length << " bp)" << std::endl;
    logger.info() << "Stored strings: " << this->stored_strings << " strings, " << this->stored_length << " bp ("
                  << ratio(this->stored_length, this->chopped_length) << "x of chopped sequence, max "
                  << this->max_string_length << " bp)" << std::endl;

    // Unchopped nodes keep their sequences, and the pieces share the stored strings.
    size_t output_nodes = this->input_nodes - this->chopped_nodes + this->pieces;
    size_t output_length = this->input_length - this->chopped_length + this->stored_length;
    logger.info() << "Output graph: " << output_nodes << " nodes (" << ratio(output_nodes, this->input_nodes)
                  << "x), " << output_length << " bp of sequence to store ("
                  << ratio(output_length, this->input_length) << "x)" << std::endl;
}

//----------------------------------------------------------------------------

/// Open the file and pass it to the writer. Exit with an error if the file
/// cannot be opened or written.
void write_to_file(const std::string& filename, const std::function<void(std::ostream&)>& writer) {
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
    DedupConfig config(argc, argv);
    double start = gbwt::readTimer();
    double checkpoint = start;
    std::unique_ptr<MutablePathDeletableHandleGraph> graph =
        vg::io::VPKG::load_one<MutablePathDeletableHandleGraph>(config.graph_file);
    if (config.verbose) {
        logger.info() << "Loaded the graph in " << (gbwt::readTimer() - checkpoint) << " s" << std::endl;
    }

    checkpoint = gbwt::readTimer();
    ChoppingPlan plan = plan_content_chopping(*graph, config.parameters);
    if (config.verbose) {
        logger.info() << "Planned the chopping in " << (gbwt::readTimer() - checkpoint) << " s with "
                      << omp_get_max_threads() << " threads" << std::endl;
        ChoppingSummary(*graph, plan, config.parameters).report();
    }

    if (!config.fasta_output.empty()) {
        checkpoint = gbwt::readTimer();
        write_to_file(config.fasta_output, [&](std::ostream& out) {
            write_stored_strings_fasta(plan, config.name_prefix, out);
        });
        if (config.verbose) {
            logger.info() << "Wrote the stored strings in " << (gbwt::readTimer() - checkpoint) << " s" << std::endl;
        }
    }

    if (config.needs_chopped_graph()) {
        checkpoint = gbwt::readTimer();
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
            logger.info() << "Chopped the graph in " << (gbwt::readTimer() - checkpoint) << " s" << std::endl;
        }

        if (!config.graph_output.empty()) {
            checkpoint = gbwt::readTimer();
            vg::io::save_handle_graph(graph.get(), config.graph_output);
            if (config.verbose) {
                logger.info() << "Wrote the chopped graph in " << (gbwt::readTimer() - checkpoint) << " s" << std::endl;
            }
        }

        if (!config.gfa_output.empty()) {
            checkpoint = gbwt::readTimer();
            write_to_file(config.gfa_output, [&](std::ostream& out) {
                write_chopped_gfa(*graph, directory, config.name_prefix, out);
            });
            if (config.verbose) {
                logger.info() << "Wrote the GFA in " << (gbwt::readTimer() - checkpoint) << " s" << std::endl;
            }
        }
    }

    if (config.verbose) {
        logger.info() << "Used " << (gbwt::readTimer() - start) << " s and "
                      << gbwt::inGigabytes(gbwt::memoryUsage()) << " GiB" << std::endl;
    }
    return 0;
}

//----------------------------------------------------------------------------

// Register subcommand
static Subcommand vg_dedup("dedup", "chop long nodes so that repeated sequence can be stored once", TOOLKIT, main_dedup);
