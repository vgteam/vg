/** \file chains_main.cpp
 *
 * Defines the "vg chains" subcommand, which extracts the handles in top-level chains from a distance index or a snarl file.
 *
 * NOTE: If a component (~chromosome) does not start with a shared node, the top-level chain is not unique.
 * For example, if the component starts with edges (u, w) and (v, w), the chain may start with either u or v.
 *
 * TODO: Option to get chain (~contig, chromosome) names from a graph.
 *
 * TODO: Full GFA format with segments and jumps.
 */

#include "subcommand.hpp"

#include "../gbwt_helper.hpp"
#include "../snarl_distance_index.hpp"
#include "../snarls.hpp"

#include <sdsl/int_vector.hpp>
#include <vg/io/vpkg.hpp>

#include <fstream>
#include <iostream>
#include <tuple>
#include <vector>

#include <getopt.h>

using namespace vg;

//----------------------------------------------------------------------------

struct ChainsConfig {
    std::string graph_file;
    std::string input_file, output_file;
    enum Format {
        BINARY,
        GFA
    } output_format = BINARY;
    bool progress = false;

    ChainsConfig(int argc, char** argv, const Logger& logger);
};

void write_binary(const std::vector<sdsl::int_vector<>>& chains, std::ostream& out);

void write_gfa_paths(const std::vector<sdsl::int_vector<>>& chains, std::ostream& out);

gbwt::vector_type extract_chain(const SnarlDistanceIndex& distance_index, const HandleGraph& graph, net_handle_t chain, size_t chain_id);

gbwt::vector_type extract_chain(const SnarlManager& snarls, const HandleGraph& graph, const Chain& chain, size_t chain_id);

sdsl::int_vector<> normalize_chain(gbwt::vector_type& chain);

//----------------------------------------------------------------------------

int main_chains(int argc, char** argv) {
    Logger logger("vg chains");
    ChainsConfig config(argc, argv, logger);

    if (config.progress) {
        logger.info() << "Loading graph from " << config.graph_file << std::endl;
    }
    std::unique_ptr<HandleGraph> graph = io::VPKG::load_one<HandleGraph>(config.graph_file);

    if (config.progress) {
        logger.info() << "Loading distance index or snarl file from " << config.input_file << std::endl;
    }
    std::unique_ptr<SnarlDistanceIndex> distance_index;
    std::unique_ptr<SnarlManager> snarls;
    std::tie(distance_index, snarls) = io::VPKG::try_load_first<SnarlDistanceIndex, SnarlManager>(config.input_file);
    if (distance_index != nullptr) {
        if (config.progress) {
            logger.info() << "Found a distance index" << std::endl;
        }
    } else if (snarls != nullptr) {
        if (config.progress) {
            logger.info() << "Found a snarl file" << std::endl;
        }
    } else {
        logger.error() << "unable to load distance index or snarl file from " 
                       << config.input_file << std::endl;
    }

    if (config.progress) {
        logger.info() << "Extracting chains" << std::endl;
    }
    std::vector<sdsl::int_vector<>> chains;
    if (distance_index != nullptr) {
        size_t chain_id = 0;
        distance_index->for_each_child(distance_index->get_root(), [&](const handlegraph::net_handle_t& chain) {
            gbwt::vector_type buffer = extract_chain(*distance_index, *graph, chain, chain_id);
            if (!buffer.empty()) {
                chains.push_back(normalize_chain(buffer));
            }
            chain_id++;
        });
    } else {
        size_t chain_id = 0;
        snarls->for_each_top_level_chain([&](const Chain* chain) {
            gbwt::vector_type buffer = extract_chain(*snarls, *graph, *chain, chain_id);
            if (!buffer.empty()) {
                chains.push_back(normalize_chain(buffer));
            }
            chain_id++;
        });
    }
    if (config.progress) {
        size_t total_handles = 0;
        for (const auto& chain : chains) {
            total_handles += chain.size();
        }
        logger.info() << "Extracted " << chains.size() << " chains with " << total_handles << " handles" << std::endl;
    }
    std::sort(chains.begin(), chains.end());

    std::ostream* out_stream;
    std::ofstream out_file;
    if (config.output_file.empty()) {
        if (config.progress) {
            logger.info() << "Writing output to stdout" << std::endl;
        }
        out_stream = &std::cout;
    } else {
        if (config.progress) {
            logger.info() << "Writing output to " << config.output_file << std::endl;
        }
        out_file.open(config.output_file, std::ios_base::binary);
        out_stream = &out_file;
    }
    if (config.output_format == ChainsConfig::BINARY) {
        if (config.progress) {
            logger.info() << "Writing binary format" << std::endl;
        }
        write_binary(chains, *out_stream);
    } else if (config.output_format == ChainsConfig::GFA) {
        if (config.progress) {
            logger.info() << "Writing GFA paths" << std::endl;
        }
        write_gfa_paths(chains, *out_stream);
    } else {
        logger.error() << "unknown output format" << std::endl;
        return 1;
    }
    out_file.close();

    return 0;
}

static vg::subcommand::Subcommand vg_chains("chains", "extract handles in top-level chains", 
                                            vg::subcommand::WIDGET, main_chains);

//----------------------------------------------------------------------------

void help_chains(char** argv) {
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " [options] graph input > output" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Extracts handles in top-level chains from a distance index or a snarl file." << std::endl;
    std::cerr << std::endl;

    std::cerr << "Output options:" << std::endl;
    std::cerr << "  -o, --output FILE  write the output to FILE" << std::endl;
    std::cerr << "  -b, --binary       output binary format (default)" << std::endl;
    std::cerr << "  -g, --gfa          output GFA paths using jumps" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Other options:" << std::endl;
    std::cerr << "  -h, --help         print this help message to stderr and exit" << std::endl;
    std::cerr << "  -p, --progress     report progress to stderr" << std::endl;
    std::cerr << std::endl;
}

ChainsConfig::ChainsConfig(int argc, char** argv, const Logger& logger) {
    static struct option long_options[] = {
        { "output", required_argument, nullptr, 'o' },
        { "binary", no_argument, nullptr, 'b' },
        { "gfa", no_argument, nullptr, 'g' },
        { "help", no_argument, nullptr, 'h' },
        { "progress", no_argument, nullptr, 'p' },
        { 0, 0, 0, 0 }
    };

    optind = 2; // force optind past command positional argument
    while (true) {
        int option_index = 0;
        int c = getopt_long(argc, argv, "o:bgh?p", long_options, &option_index);
        if (c == -1) { break; } // End of options.

        switch (c)
        {
        case 'o':
            this->output_file = ensure_writable(logger, optarg);
            break;
        case 'b':
            this->output_format = BINARY;
            break;
        case 'g':
            this->output_format = GFA;
            break;

        case 'h':
        case '?':
            help_chains(argv);
            std::exit(EXIT_FAILURE);
        case 'p':
            this->progress = true;
            break;

        default:
            std::abort();
        }
    }

    if (optind + 2 != argc) {
        help_chains(argv);
        std::exit(EXIT_FAILURE);
    }
    this->graph_file = argv[optind]; optind++;
    this->input_file = argv[optind]; optind++;
}

//----------------------------------------------------------------------------

void write_binary(const std::vector<sdsl::int_vector<>>& chains, std::ostream& out) {
    std::uint64_t num_chains = chains.size();
    sdsl::simple_sds::serialize_value(num_chains, out);
    for (const auto& chain : chains) {
        chain.simple_sds_serialize(out);
    }
}

void write_gfa_paths(const std::vector<sdsl::int_vector<>>& chains, std::ostream& out) {
    for (size_t i = 0; i < chains.size(); i++) {
        const auto& chain = chains[i];
        out << "P\t" << i << "\t";
        for (size_t j = 0; j < chain.size(); j++) {
            nid_t id = gbwt::Node::id(chain[j]);
            bool is_reverse = gbwt::Node::is_reverse(chain[j]);
            if (j > 0) {
                // These are typically jumps, not edges.
                out << ";";
            }
            out << id << (is_reverse ? "-" : "+");
        }
        out << "\t*" << std::endl;
    }
}

net_handle_t follow_chain(const SnarlDistanceIndex& distance_index, const HandleGraph& graph, 
                          size_t chain_id, net_handle_t curr) {
    net_handle_t next = curr;
    size_t successors = 0;
    distance_index.follow_net_edges(next, &graph, false, [&](const net_handle_t& child) {
        successors++;
        next = child;
    });
    if (successors != 1) {
        logging::error("follow_chain()") << "chain " << chain_id 
                                         << " has " << successors << " successors for a child";
    }
    return next;
}

void try_append(gbwt::vector_type& chain, gbwt::node_type start, gbwt::node_type end) {
    if (start == gbwt::ENDMARKER || end == gbwt::ENDMARKER) {
        return;
    }
    if (chain.empty() || chain.back() != start) {
        chain.push_back(start);
    }
    if (chain.empty() || chain.back() != end) {
        chain.push_back(end);
    }
}

gbwt::vector_type extract_chain(const SnarlDistanceIndex& distance_index, const HandleGraph& graph, 
                                net_handle_t chain, size_t chain_id) {
    gbwt::vector_type result;

    // Closed interval of net handles.
    net_handle_t curr = distance_index.get_bound(chain, false, true);
    net_handle_t chain_end = distance_index.get_bound(chain, true, false);
    gbwt::node_type prev = gbwt::ENDMARKER; // Possible start of a snarl.
    bool was_snarl = false;
    while (true) {
        if (distance_index.is_node(curr)) {
            handle_t handle = distance_index.get_handle(curr, &graph);
            gbwt::node_type node = handle_to_gbwt(graph, handle);
            if (was_snarl) {
                try_append(result, prev, node);
            }
            prev = node;
            was_snarl = false;
        } else if (distance_index.is_snarl(curr)) {
            // Trivial snarls do not show up in the traversal.
            was_snarl = true;
        } else {
            was_snarl = false;
        }
        if (curr == chain_end) {
            break;
        }
        curr = follow_chain(distance_index, graph, chain_id, curr);
    }

    for (auto handle : result) {
        nid_t id = gbwt::Node::id(handle);
        if (!graph.has_node(id)) {
            logging::error("extract_chain()") << "chain " << chain_id
                                              << " has a handle for missing node " << id << std::endl;
        }
    }

    return result;
}

gbwt::vector_type extract_chain(const SnarlManager& snarls, const HandleGraph& graph, 
                                const Chain& chain, size_t chain_id) {
    gbwt::vector_type result;

    for (auto iter = chain_begin(chain); iter != chain_end(chain); ++iter) {
        const Snarl* snarl = iter->first;
        if (snarls.is_trivial(snarl, graph)) {
            // Skip trivial snarls.
            continue;
        }
        gbwt::node_type start = gbwt::Node::encode(snarl->start().node_id(), snarl->start().backward());
        gbwt::node_type end = gbwt::Node::encode(snarl->end().node_id(), snarl->end().backward());
        if (iter->second) {
            // Reverse snarl; we must swap and flip the endpoints.
            std::swap(start, end);
            start = gbwt::Node::reverse(start);
            end = gbwt::Node::reverse(end);
        }
        try_append(result, start, end);
    }

    return result;
}

sdsl::int_vector<> normalize_chain(gbwt::vector_type& chain) {
    // Normalize the orientation.
    // TODO: What if there is the same number of forward and reverse nodes?
    size_t reverse_nodes = 0;
    for (const auto& node : chain) {
        if (gbwt::Node::is_reverse(node)) {
            reverse_nodes++;
        }
    }
    if (reverse_nodes > chain.size() / 2) {
        gbwt::reversePath(chain);
    }

    // Bit compress the vector.
    size_t width = 1;
    auto iter = std::max_element(chain.begin(), chain.end());
    if (iter != chain.end()) {
        width = sdsl::bits::length(*iter);
    }
    sdsl::int_vector<> result(chain.size(), 0, width);
    for (size_t i = 0; i < chain.size(); i++) {
        result[i] = chain[i];
    }

    return result;
}

//----------------------------------------------------------------------------
