#include "random_graph.hpp"

#include <random>
#include <time.h>
#include <iostream>

#include "randomness.hpp"

namespace vg {
namespace unittest {

using namespace std;

void random_graph(int64_t seq_size, int64_t variant_len, int64_t variant_count,
                  MutablePathMutableHandleGraph* graph) {
    vector<int64_t> seq_sizes = {seq_size};
    return random_graph(seq_sizes, variant_len, variant_count, graph);
}
    
 
void random_graph(vector<int64_t> seq_sizes, int64_t variant_len, int64_t total_variant_count,
                  MutablePathMutableHandleGraph* graph) {           
#ifdef debug
    cerr << "Make random graphs with " << total_variant_count << " ~" << variant_len << "-bp large variants" << endl;
#endif

    //How many variants for each component?
    int64_t total_sequence_length = 0;
    for (int64_t len : seq_sizes) {
        total_sequence_length += len;
    }
    vector<int64_t> variant_counts;
    for (int64_t len : seq_sizes) {
        variant_counts.emplace_back(round((len* total_variant_count) / total_sequence_length) );
    }

    for (size_t i = 0 ; i < seq_sizes.size() ; i++) {

        int64_t seq_size = seq_sizes[i];
        int64_t variant_count = variant_counts[i];

        //Create a random graph for a sequence of length seq_size
        //variant_len is the mean length of a larger variation and variationCount
        //is the number of variations in the graph

        map<size_t, id_t> index_to_id;
                      //Index of original sequence to node that starts at that index
        
        //Get random number generator for lengths and variation types
        default_random_engine generator(test_seed_source());
        
        uniform_int_distribution<int> base_distribution(0, 3);
        poisson_distribution<size_t> length_distribution(variant_len);
        uniform_int_distribution<int> variant_distribution(0, 4);
        uniform_int_distribution<int> index_distribution(0, seq_size-1);
        
        // Init the string with placeholder bases
        string seq(seq_size, 'A');
        // Set the string to random bases
        string alphabet = "ACGT";
        for (size_t i = 0; i < seq.size(); i++) {
            seq[i] = alphabet[base_distribution(generator)];
        }
        
        handle_t h = graph->create_handle(seq);
        
        path_handle_t p = graph->create_path_handle("path");
        graph->append_step(p, h);
        
        index_to_id[0] = graph->get_id(h);

        auto do_split = [&] (size_t index) -> pair<handle_t, handle_t> {
            //Split graph at index and update the path index
            
#ifdef debug
            cerr << "Split original sequence at base " << index << "/" << seq_size << endl;
#endif
            
            // Make sure we aren't trying to split out of bounds
            assert(index < seq_size);
            
            auto n = --index_to_id.upper_bound(index); //orig node containing pos
            size_t first_index = n->first;               //Index of first node
            handle_t first_handle = graph->get_handle(n->second);    //handle of first node
            
#ifdef debug
            cerr << "\tFalls in node " << graph->get_id(first_handle) << " length " << graph->get_length(first_handle) << " at " << first_index << endl;
#endif
            
            pair<handle_t, handle_t> return_val;
            
            if (index > first_index) {
                // it's in the middle of a node, do a split
                
                size_t split_length = index - first_index;
                
#ifdef debug
                cerr << "\t\tSplit at " << split_length << " along node" << endl;
#endif
                
                return_val = graph->divide_handle(first_handle, split_length);
                
#ifdef debug
                cerr << "\t\t\tCreate " << graph->get_id(return_val.first) << " at " << first_index << endl;
                cerr << "\t\t\tCreate " << graph->get_id(return_val.second) << " at " << index << endl;
#endif
                
                index_to_id[first_index] = graph->get_id(return_val.first);
                index_to_id[index] = graph->get_id(return_val.second);

            }
            else {
                return_val.second = first_handle;
                --n;
                return_val.first = graph->get_handle(n->second);
                
            }
            
            return return_val;
        };

        enum VariationType {SNP = 0, POINT_INDEL = 1, STRUCTURAL_INDEL = 2,
                            CNV = 3, INVERSION = 4};

        for (int j = 0; j < variant_count; j++) {
            //add variants
            int start_index = index_distribution(generator);
            VariationType variation_type = (VariationType) variant_distribution(generator);

            switch (variation_type) {
                case SNP:
                {
                    
                    if (start_index == 0) {

                        handle_t new_node = graph->create_handle(string(1, alphabet[base_distribution(generator)]));
                        
                        pair<handle_t, handle_t> end_nodes = do_split(start_index+1);
                        graph->create_edge(new_node, end_nodes.second);
                        
                    }
                    else if (start_index < seq_size - 2) {
                        
                        handle_t new_node = graph->create_handle(string(1, alphabet[base_distribution(generator)]));
                        
                        pair<handle_t, handle_t> start_nodes = do_split(start_index);
                        pair<handle_t, handle_t> end_nodes = do_split(start_index+1);
                        
                        graph->create_edge(start_nodes.first, new_node);
                        graph->create_edge(new_node, end_nodes.second);
                        
                    }
                    else if (start_index == seq_size - 2) {
                        
                        handle_t new_node = graph->create_handle(string(1, alphabet[base_distribution(generator)]));
                        
                        pair<handle_t, handle_t> start_nodes = do_split(start_index);
                        
                        graph->create_edge(start_nodes.first, new_node);
                        
                    }
                    break;
                }
                    
                case POINT_INDEL:
                {
                    if (start_index > 0 && start_index < seq_size-1) {
                        pair<handle_t, handle_t> start_nodes = do_split(start_index);
                        pair<handle_t, handle_t> end_nodes = do_split(start_index+1);
                        
                        graph->create_edge(start_nodes.first, end_nodes.second);
                        
                    }
                    break;
                }
                    
                case STRUCTURAL_INDEL:
                {
                    //long indel
                    size_t length = length_distribution(generator);
                    
                    if (length > 0 && start_index > 0 && length + start_index < seq_size-1) {
                        
                        pair<handle_t, handle_t> start_nodes = do_split(start_index);
                        pair<handle_t, handle_t> end_nodes = do_split(start_index+length);
                        
                        graph->create_edge(start_nodes.first, end_nodes.second);
                        
                    }
                    break;
                }
                    
                case CNV:
                {
                    //Copy number variation
                    size_t length = length_distribution(generator);
                    
                    if (length > 0) {
                        if (start_index == 0 && length + start_index < seq_size - 1) {
                            // Very first base is involved, but very last base isn't (we don't handle that case)
                            pair<handle_t, handle_t> end_nodes = do_split(start_index+length);
                            
                            auto node_pair = index_to_id.begin();//first node
                            handle_t first_handle = graph->get_handle(node_pair->second);
                            
                            graph->create_edge(end_nodes.first, first_handle);
                            
                        }
                        else if (length + start_index < seq_size - 1) {
                            // Only interior bases are involved
                            pair<handle_t, handle_t> start_nodes = do_split(start_index);
                            pair<handle_t, handle_t> end_nodes = do_split(start_index+length);
                            // hack: this won't redo the split since it's already there, but it
                            // will get the second handle again, which may have been invalidated
                            // if we split the second node
                            start_nodes = do_split(start_index);
                            
                            graph->create_edge(end_nodes.first, start_nodes.second);
                            
                        }
                        else if (length + start_index == seq_size - 1) {
                            // Very last base is involved (but not very first)
                            pair<handle_t, handle_t> start_nodes = do_split(start_index);
                            
                            //last node
                            auto node_pair = --index_to_id.end();
                            handle_t last_handle = graph->get_handle(node_pair->second);
                            
                            graph->create_edge(last_handle, start_nodes.second);
                        }
                    }
                    break;
                }
                    
                case INVERSION:
                {
                    //Inversion
                    size_t length = length_distribution(generator);
                    if (length > 0) {
                        if (start_index == 0 && length + start_index < seq_size - 1) {
                            // Very first base is involved, but very last base isn't (we don't handle that case)
                            pair<handle_t, handle_t> end_nodes = do_split(start_index+length);
                            
                            graph->create_edge(graph->flip(end_nodes.first), end_nodes.second);
                            
                            
                        } else if (length + start_index < seq_size-1) {
                            // Only interior bases are involved
                            pair<handle_t, handle_t> start_nodes = do_split(start_index);
                            pair<handle_t, handle_t> end_nodes = do_split(start_index+length);
                            // hack: this won't redo the split since it's already there, but it
                            // will get the second handle again, which may have been invalidated
                            // if we split the second node
                            start_nodes = do_split(start_index);
                            
                            graph->create_edge(start_nodes.first, graph->flip(end_nodes.first));
                            graph->create_edge(graph->flip(start_nodes.second), end_nodes.second);
                            
                        } else if (length + start_index == seq_size - 1) {
                            // Very last base is involved (but not very first)
                            pair<handle_t, handle_t> start_nodes = do_split(start_index);
                            
                            graph->create_edge(start_nodes.first, graph->flip(start_nodes.second));
                        }
                    }
                    break;
                }
                    
                default:
                    break;
            }
        }
    }
};

vector<vector<size_t>> random_adjacency_list(size_t node_count, size_t edge_count) {

    // Work out how to randomly pick a node
    default_random_engine generator(test_seed_source());
    uniform_int_distribution<size_t> node_distribution(0, node_count - 1);
    
    vector<vector<size_t>> to_return(node_count);
    
    if (node_count == 0) {
        // Can't have edges without nodes
        return to_return;
    }
    
    for (size_t i = 0; i < edge_count; i++) {
        // Pick two nodes
        size_t a = node_distribution(generator);
        size_t b = node_distribution(generator);
    
        // Record the edge in both directions
        to_return[a].push_back(b);
        if (a != b) {
            // Unless it's a self loop
            to_return[b].push_back(a);
        }
    }
    
    return to_return;
    
}
    

}
}
