#include "convert_handle.hpp"
namespace vg {
using namespace std;
    void convert_handle_graph(const HandleGraph* converting, MutableHandleGraph* converted) {
        assert(converted->node_size() == 0);
        
        // If xg is a null pointer, throw a runtime error
        if (converting == nullptr) {
            throw runtime_error("There is no graph to convert from");
        }
        if (converted == nullptr) {
            throw runtime_error("There is no graph to convert to");
        }
        assert(converted->node_size() == 0);

        // Iterate through each handle in xg and create the same handle in mutable graph
        converting->for_each_handle([&](const handle_t& here) {
            // Get the id of the graph handle
            id_t converting_id = converting->get_id(here);
            // Get the sequence of the graph handle
            string converting_seq = converting->get_sequence(here);
            // Create a handle in mutable graph using the graph id and sequence
            converted->create_handle(converting_seq, converting_id);
        });
        // add any edges that are not yet present
        converting->for_each_edge([&](const edge_t& edge) {
            edge_t converted_edge(converted->get_handle(converting->get_id(edge.first), converting->get_is_reverse(edge.first)),
                             converted->get_handle(converting->get_id(edge.second), converting->get_is_reverse(edge.second)));
            if (!converted->has_edge(converted_edge)) {
                converted->create_edge(converted_edge);
            }
            // always keep going
            return true;
        });
    }

    void convert_path_handle_graph(const PathHandleGraph* converting, MutablePathDeletableHandleGraph* converted) {
        assert(converted->get_path_count() == 0);
        
        // Must convert nodes and edges before converting paths.
        convert_handle_graph(converting, converted);
        
        // Iterate through paths and add paths from converting to converted
        converting->for_each_path_handle([&](const path_handle_t& path) {
            // create path
            string path_name = converting->get_path_name(path);
            path_handle_t converted_path = converted->create_path_handle(path_name);
            // find first occurence in path
            occurrence_handle_t converting_occurence = converting->get_first_occurrence(path);
            handle_t converting_handle = converting->get_occurrence(converting_occurence);
            handle_t converted_handle = converted->get_handle(converting->get_id(converting_handle), converting->get_is_reverse(converting_handle));
            occurrence_handle_t converted_occurance_handle = converted->append_occurrence(converted_path, converted_handle);
            // iterate through path and add each handle to converted graph
            while (converting->has_next_occurrence(converting_occurence)){
                converting_occurence = converting->get_next_occurrence(converting_occurence);
                converting_handle = converting->get_occurrence(converting_occurence);
                converted_handle = converted->get_handle(converting->get_id(converting_handle), converting->get_is_reverse(converting_handle));
                occurrence_handle_t converted_occurance_handle = converted->append_occurrence(converted_path, converted_handle);
            }
        });
    }
}
