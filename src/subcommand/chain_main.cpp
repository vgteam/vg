/** \file chain_main.cpp
 *
 * Defines the "vg chain" subcommand, which runs a serialized hit chaining problem.
 */

#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <iostream>

#include "subcommand.hpp"

#include "../algorithms/chain_items.hpp"
#include "../integrated_snarl_finder.hpp"


using namespace std;
using namespace vg;
using namespace vg::subcommand;

/// Define a single-graph-node item we can use for chaining
struct LiteralItem {
    size_t start;
    size_t length;
    pos_t pos;
    
    inline bool operator==(const LiteralItem& other) const {
        return (this->start == other.start && this->length == other.length && this->pos == other.pos);
    }
    
    inline bool operator!=(const LiteralItem& other) const {
        return !(*this == other);
    }
};

std::ostream& operator<<(std::ostream& out, const LiteralItem& item) {
    return out << "{read " << item.start << "-" << (item.start + item.length) << " = graph " << item.pos << "}"; 
}

namespace vg {

namespace algorithms {

template<>
struct ChainingSpace<LiteralItem, void> : public BaseChainingSpace<LiteralItem> {
    using Item = LiteralItem;
    
    ChainingSpace(const ChainingScorer& scorer,
                  const SnarlDistanceIndex* distance_index,
                  const HandleGraph* graph) :
        BaseChainingSpace<Item>(scorer, distance_index, graph) {
        
        // Nothing to do!
    }
    
    virtual size_t read_start(const Item& item) const {
        return item.start;
    }
    
    virtual size_t read_end(const Item& item) const {
        return item.start + item.length;
    }
    
    virtual size_t read_length(const Item& item) const {
        return item.length;
    }
    
    pos_t graph_start(const Item& item) const {
        return item.pos;
    }
    
    pos_t graph_end(const Item& item) const {
        pos_t end = item.pos;
        get_offset(end) += this->graph_length(item);
        return end;
    }
    
    size_t graph_path_size(const Item& item) const {
        return 1;
    }
    
    handle_t graph_path_at(const Item& item, size_t index) const {
        return this->graph->get_handle(id(item.pos), is_rev(item.pos));
    }
    
    size_t graph_path_offset(const Item& item) const {
        return offset(item.pos);
    }
};

}

}

void help_chain(char** argv) {
    cerr << "usage: " << argv[0] << " chain [options] input.json" << endl
         << "options:" << endl
         << "    -p, --progress         show progress" << endl;
}

int main_chain(int argc, char** argv) {

    bool show_progress = false;
    
    int c;
    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"progress",  no_argument, 0, 'p'},
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "ph?",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {

        case 'p':
            show_progress = true;
            break;
            
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help_chain(argv);
            exit(1);
            break;

        default:
            abort ();

        }
    }
    
    if (optind == argc) {
        // No positional arguments
        help_chain(argv);
        exit(1);
    }
    
    string problem_filename = get_input_file_name(optind, argc, argv);
    
    if (optind != argc) {
        // Too many positional arguments
        help_chain(argv);
        exit(1);
    }
    
    // Load up the JSON file we are supposed to have.
    json_error_t json_error;
    json_t* problem_json = json_load_file(problem_filename.c_str(), 0, &json_error);
    if (!problem_json) {
        throw std::runtime_error(json_error.text);
    }
    assert(json_is_object(problem_json));
    
    if (show_progress) {
        std::cerr << "Loaded problem from " << problem_filename << std::endl;
    }
    
    // Populate the graph.
    // TODO: use more libvgio stuff somehow.
    HashGraph graph;
    json_t* graph_json = json_object_get(problem_json, "subgraph");
    if (graph_json && json_is_object(graph_json)) {
        json_t* nodes_json = json_object_get(graph_json, "node");
        if (nodes_json && json_is_array(nodes_json)) {
            // Go through the node array
            for (size_t i = 0; i < json_array_size(nodes_json); i++) {
                json_t* node_json = json_array_get(nodes_json, i);
                if (node_json && json_is_object(node_json)) {
                    // Make each node
                    const char* node_id = nullptr;
                    const char* sequence = nullptr;
                    
                    if (json_unpack_ex(node_json, &json_error, 0, "{s:s, s:s}", "id", &node_id, "sequence", &sequence) == 0) {
                        // This record is well-formed, so make the node
                        assert(node_id != nullptr);
                        assert(sequence != nullptr);
                        graph.create_handle(sequence, vg::parse<nid_t>(node_id));
                    } else {
                        std::cerr << "warning:[vg chain] Unreadable node object at index " << i << ": " << json_error.text << std::endl;
                    }
                } else {
                    std::cerr << "warning:[vg chain] No node object at index " << i << std::endl;
                }
            }
        } else {
            std::cerr << "warning:[vg chain] No nodes" << std::endl;
        }
        json_t* edges_json = json_object_get(graph_json, "edge");
        if (edges_json && json_is_array(edges_json)) {
            // Go through the edge array
            for (size_t i = 0; i < json_array_size(edges_json); i++) {
                json_t* edge_json = json_array_get(edges_json, i);
                if (edge_json && json_is_object(edge_json)) {
                    // Make each edge
                    const char* from_id = nullptr;
                    bool from_start = false;
                    const char* to_id = nullptr;
                    bool to_end = false;
                    
                    if (json_unpack_ex(edge_json, &json_error, 0, "{s:s, s?b, s:s, s?b}", "from", &from_id, "from_start", &from_start, "to", &to_id, "to_end", &to_end) == 0) {
                        // This record is well-formed, so make the edge
                        assert(from_id != nullptr);
                        assert(to_id != nullptr);
                        handle_t from_handle = graph.get_handle(vg::parse<nid_t>(from_id), from_start);
                        handle_t to_handle = graph.get_handle(vg::parse<nid_t>(to_id), to_end);
                        graph.create_edge(from_handle, to_handle);
                    } else {
                        std::cerr << "warning:[vg chain] Unreadable edge object at index " << i << ": " << json_error.text << std::endl;
                    }
                } else {
                    std::cerr << "warning:[vg chain] No edge object at index " << i << std::endl;
                }
            }
        } else {
            std::cerr << "warning:[vg chain] No edges" << std::endl;
        }
    } else {
        std::cerr << "warning:[vg chain] No graph" << std::endl;
    }
    if (show_progress) {
        std::cerr << "Reconstructed " << graph.get_node_count() << " nodes and " << graph.get_edge_count() << " edges" << std::endl;
    }
    
    if (graph.get_node_count() == 0) {
        std::cerr << "error:[vg chain] Cannot build indexes for an empty graph" << std::endl;
        exit(1);
    }
    
    // Create the chaining space based on it
    IntegratedSnarlFinder snarl_finder(graph);
    SnarlDistanceIndex distance_index;
    fill_in_distance_index(&distance_index, &graph, &snarl_finder);
    if (show_progress) {
        std::cerr << "Built distance index" << std::endl;
    }
    
    vg::algorithms::IndelOnlyChainingScorer scorer;
    vg::algorithms::ChainingSpace<LiteralItem> space(scorer, &distance_index, &graph);
    if (show_progress) {
        std::cerr << "Built chaining space" << std::endl;
    }
    
    // Create all the items to chain
    std::vector<LiteralItem> items;
    json_t* items_json = json_object_get(problem_json, "items");
    if (items_json && json_is_array(items_json)) {
        items.reserve(json_array_size(items_json));
        for (size_t i = 0; i < json_array_size(items_json); i++) {
            json_t* item_json = json_array_get(items_json, i);
            if (item_json && json_is_object(item_json)) {
                // For each chainable item we got
                const char* read_start = nullptr;
                const char* read_end = nullptr;
                json_t* graph_start = nullptr;
                const char* graph_start_id = nullptr;
                const char* graph_start_offset = "0"; 
                bool graph_start_is_reverse = false;
                json_t* graph_end = nullptr;
                const char* graph_end_id = nullptr;
                const char* graph_end_offset = "0";
                bool graph_end_is_reverse = false;
                if (json_unpack_ex(item_json, &json_error, 0, "{s:s, s:s, s:o, s:o}",
                                   "read_start", &read_start, 
                                   "read_end", &read_end,
                                   "graph_start", &graph_start,
                                   "graph_end", &graph_end) == 0 &&
                    json_unpack_ex(graph_start, &json_error, 0, "{s:s, s?s, s?b}",
                                   "node_id", &graph_start_id, "offset", &graph_start_offset, "is_reverse", &graph_start_is_reverse) == 0 &&
                    json_unpack_ex(graph_end, &json_error, 0, "{s:s, s?s, s?b}",
                                   "node_id", &graph_end_id, "offset", &graph_end_offset, "is_reverse", &graph_end_is_reverse) == 0) {
                    
                    // We have an item record.
                    assert(read_start != nullptr);
                    assert(read_end != nullptr);
                    assert(graph_start_id != nullptr);
                    assert(graph_end_id != nullptr);
                   
                    // We can only handle items where they occupy space on just one node.
                    assert(strcmp(graph_start_id, graph_end_id) == 0);
                    
                    // And we assume the lengths match.
                    size_t start = vg::parse<size_t>(read_start);
                    size_t length = vg::parse<size_t>(read_end) - start;
                    
                    // Pack up into an item
                    items.emplace_back(LiteralItem {start, length, make_pos_t(vg::parse<nid_t>(graph_start_id), graph_start_is_reverse, vg::parse<size_t>(graph_start_offset))});
                } else {
                    std::cerr << "warning:[vg chain] Unreadable item object at index " << i << ": " << json_error.text << std::endl;
                }
            } else {
                std::cerr << "warning:[vg chain] No item object at index " << i << std::endl;
            }
        }
    } else {
        std::cerr << "warning:[vg chain] No items" << std::endl;
    }
    if (show_progress) {
        std::cerr << "Reconstructed " << items.size() << " chainable items" << std::endl;
    }
    
    // Now we have parsed the JSON, so throw it out.
    json_decref(problem_json);
    problem_json = nullptr;
    
    // Do the chaining. We assume items is already sorted right.
    std::pair<int, std::vector<size_t>> score_and_chain = vg::algorithms::find_best_chain(items, space);
    
    std::cout << "Best chain gets score " << score_and_chain.first << std::endl;
    
    return 0;
}

// Register subcommand
static Subcommand vg_chain("chain", "run a serialized chaining problem", DEVELOPMENT, main_chain);

