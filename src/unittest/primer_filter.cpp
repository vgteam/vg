//
//  primer_filter.cpp
//
//  Unit tests for primer filter 
//

#include <stdio.h>
#include <iostream>
#include <regex>
#include <fstream>
#include <vector>
#include <sstream>
#include <set>
#include "vg/io/json2pb.h"
#include <vg/vg.pb.h>
#include "catch.hpp"
#include "random_graph.hpp"
#include "randomness.hpp"
#include "../snarl_distance_index.hpp"
#include "../integrated_snarl_finder.hpp"
#include "../genotypekit.hpp"
#include "../traversal_finder.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>
#include "xg.hpp"

namespace vg {
    namespace unittest {

        TEST_CASE( "filter simple primers",
                  "[primers]" ) {
            
            struct Primer {
                string sequence;
                bool left = true;
                size_t position;
                size_t length;
                size_t offset;
                vector<size_t> mapped_nodes_ids;
            };

            struct Primer_pair {
                Primer left_primer;
                Primer right_primer;
                size_t linear_product_size;
                size_t min_product_size;
                size_t max_product_size;
            };
            
            class Primer_finder {
            private:
                vector<Primer_pair> primer_pairs;
                // HandleGraph* graph;
                PathPositionHandleGraph* graph;
                SnarlDistanceIndex* distance_index;
                path_handle_t reference_path_handle; 
                vector<Primer_pair> selected_primer_pairs;
                
            public:
                Primer_finder() = default;
                Primer_finder(
                //unique_ptr<handlegraph::HandleGraph>& graph_param,
                unique_ptr<handlegraph::PathPositionHandleGraph>& graph_param,
                string reference_path_name,
                SnarlDistanceIndex* distance_index_param) {
                    graph = graph_param.get();
                    reference_path_handle = graph->get_path_handle("y");
                    distance_index = distance_index_param;
                }
                ~Primer_finder() = default;
                

                void load_primers(string path_to_primers) {
                    regex left_seq_pattern("PRIMER_LEFT_\\d+_SEQUENCE=(\\w+)");
                    regex right_seq_pattern("PRIMER_RIGHT_\\d+_SEQUENCE=(\\w+)");
                    regex left_pos_pattern("PRIMER_LEFT_\\d+=(\\d+,\\d+)");
                    regex right_pos_pattern("PRIMER_RIGHT_\\d+=(\\d+,\\d+)");
                    
                    Primer left_primer {""};
                    Primer right_primer {"", false};

                    ifstream file_handle(path_to_primers);
                    if (file_handle.is_open()) {
                        string line;
                        while (getline(file_handle, line)) {
                            line = rstrip(line);
                            smatch match;
                            if (regex_search(line, match, left_seq_pattern)) {
                                if (right_primer.sequence != "") {
                                    map_to_nodes(left_primer);
                                    map_to_nodes(right_primer);
                                    Primer_pair primer_pair{left_primer, right_primer,
                                        right_primer.position - left_primer.position + 1};
                                    primer_pairs.push_back(primer_pair);
                                    left_primer = {""};
                                    right_primer = {"", false};
                                }
                                left_primer.sequence = match[1];
                            } else if (regex_search(line, match, right_seq_pattern)) {
                                right_primer.sequence = match[1];
                            } else if (regex_search(line, match, left_pos_pattern)) {
                                vector<string> pos_and_len = split(match[1], ",");
                                left_primer.position = stoi(pos_and_len[0]);
                                left_primer.length = stoi(pos_and_len[1]);
                            } else if (regex_search(line, match, right_pos_pattern)) {
                                vector<string> pos_and_len = split(match[1], ",");
                                right_primer.length = stoi(pos_and_len[1]);
                                right_primer.position = stoi(pos_and_len[0]) - right_primer.length;
                                
                            }
                        }
                        map_to_nodes(left_primer);
                        map_to_nodes(right_primer);
                        Primer_pair primer_pair{left_primer, right_primer, 
                            right_primer.position - left_primer.position + 1};
                        primer_pairs.push_back(primer_pair);
                    }
                }


                void map_to_nodes(Primer& primer) {
                    string primer_seq;
                    if (primer.left) {
                        primer_seq = primer.sequence;
                    } else {
                        primer_seq = revcomp(primer.sequence);
                    }
                    step_handle_t  cur_node_step_handle = graph->get_step_at_position(reference_path_handle, primer.position);
                    handle_t cur_node_handle = graph->get_handle_of_step(cur_node_step_handle);
                    primer.mapped_nodes_ids.push_back(graph->get_id(cur_node_handle));
                    string cur_node_sequence = graph->get_sequence(cur_node_handle);
                    size_t primer_matched_index =  longest_match_len(primer, cur_node_sequence, primer_seq) - 1;
                    while (primer_matched_index < primer_seq.size()-1) {
                        cur_node_step_handle = graph->get_next_step(cur_node_step_handle);
                        cur_node_handle = graph->get_handle_of_step(cur_node_step_handle);
                        primer.mapped_nodes_ids.push_back(graph->get_id(cur_node_handle));
                        cur_node_sequence = graph->get_sequence(cur_node_handle);
                        string primer_substr = primer_seq.substr(primer_matched_index + 1, primer.length - primer_matched_index - 1);
                        primer_matched_index += longest_match_len(primer, primer_substr, cur_node_sequence);
                    }
                }


                void run_test() {
                        cout << "testing HandleGraph..." << endl;
                        nid_t min_node_id = graph->min_node_id();
                        nid_t max_node_id = graph->max_node_id();
                        handle_t min_node = graph->get_handle(min_node_id);
                        handle_t max_node = graph->get_handle(max_node_id);
                        cout << "min node id: " << min_node_id << endl;
                        cout << "sequence: " << graph->get_sequence(min_node) << endl;
                        cout << "max node id: " << max_node_id << endl;
                        cout << "sequence: " << graph->get_sequence(max_node) << endl;
                        cout << "HandleGraph works! :)" << endl;

                        cout << "-------------------------------------" << endl;

                        cout << "testing PathHandleGraph..." << endl;
                        cout << graph->get_path_count() << endl;
                        cout << "path with name y exists: " << graph->has_path("y") << endl;
                        cout << "reference path has " << graph->get_step_count(reference_path_handle) << 
                            " node steps" << endl;
                        cout << "PathHandleGraph works! :)" << endl;

                        cout << "-------------------------------------" << endl;

                        cout << "testing PathPositionHandleGraph..." << endl;
                        cout << "referecne path length: " << graph->get_path_length(reference_path_handle) << endl;
                        step_handle_t  step_handle_lprimer = graph->get_step_at_position(reference_path_handle, 362); // 362 is the position of the left primer of the first primer pair
                        handle_t handle_lprimer = graph->get_handle_of_step(step_handle_lprimer);
                        cout << "left primer sequence begin node seq: " << graph->get_sequence(handle_lprimer) << endl;
                        step_handle_t step_handle_lprimer_next = graph->get_next_step(step_handle_lprimer);
                        handle_t handle_lprimer_next = graph->get_handle_of_step(step_handle_lprimer_next);
                        cout << "left priemr sequence second node seq: " << graph->get_sequence(handle_lprimer_next) << endl;

                        // Get node id for position 0
                        step_handle_t position_0_step_handle = graph->get_step_at_position(reference_path_handle, 0);
                        handle_t position_0_handle = graph->get_handle_of_step(position_0_step_handle);
                        cout << "node id for position 0 " << graph->get_id(position_0_handle) << endl;
                        // Get node id for position 1
                        step_handle_t position_1_step_handle = graph->get_step_at_position(reference_path_handle, 1);
                        handle_t position_1_handle = graph->get_handle_of_step(position_1_step_handle);
                        cout << "node id for position 1 " << graph->get_id(position_1_handle) << endl;

                        // Get node id for position 31
                        step_handle_t position_31_step_handle = graph->get_step_at_position(reference_path_handle, 31);
                        handle_t position_31_handle = graph->get_handle_of_step(position_31_step_handle);
                        cout << "node id for position 31 " << graph->get_id(position_31_handle) << endl;
                        
                        // Get node id for position 32
                        step_handle_t position_32_step_handle = graph->get_step_at_position(reference_path_handle, 32);
                        handle_t position_32_handle = graph->get_handle_of_step(position_32_step_handle);
                        cout << "node id for position 32 " << graph->get_id(position_32_handle) << endl;

                        cout << "PathPositionHandleGraph works! :)" << endl;

                        cout << "-------------------------------------" << endl;
                        
                        cout << "testing SnarlDistanceIndex..." << endl; 
                        net_handle_t root_node = distance_index->get_root();
                        cout << "is root a root? " << distance_index->is_root(root_node) << endl;
                        cout << "is root a node? " << distance_index->is_node(root_node) << endl;
                        cout << "is root a snarl? " << distance_index->is_snarl(root_node) << endl;
                        cout << "is root a chain? " << distance_index->is_chain(root_node) << endl;
                        cout << "depth of root is: " << distance_index->get_depth(root_node) << endl;

                        net_handle_t min_node_net_handle = distance_index->get_net(min_node, graph);
                        cout << "depth of min node is: " << distance_index->get_depth(min_node_net_handle) << endl;
                        cout << "make sure that min node net handle is a node: " << distance_index->is_node(min_node_net_handle) << endl;
                        
                        size_t min_dist_12_17 = distance_index->minimum_distance(12, false, 1, 17, false, 2);
                        size_t max_dist_12_17 = distance_index->maximum_distance(12, false, 3, 17, false, 3);

                        cout << "min dist between node 12 and node 17: " << min_dist_12_17 << endl;
                        cout << "max dist between node 12 and node 17: " << max_dist_12_17 << endl;

                        cout << "SnarlDistanceIndex works! :)" << endl;
                        cout << "-------------------------------------" << endl;

                        cout << "testing load_primers..." << endl;
                        for (vector<Primer_pair>::iterator it = primer_pairs.begin(); it != primer_pairs.end(); ++it) {
                            Primer left_primer = it->left_primer;
                            Primer right_primer = it->right_primer;
                            cout << "product size: " << it->linear_product_size << endl;
                            cout << left_primer.left << " " << left_primer.position << " " << 
                                left_primer.length << " " << left_primer.sequence << endl;
                            for (int i = 0; i < left_primer.mapped_nodes_ids.size(); i++) {
                                size_t cur_node_id = left_primer.mapped_nodes_ids[i];
                                handle_t cur_node_handle = graph->get_handle(cur_node_id);
                                cout << graph->get_sequence(cur_node_handle) << " ";
                            }
                            cout << endl;
                            cout << right_primer.left << " " << right_primer.position << " " << 
                                right_primer.length << " " << right_primer.sequence << " " << revcomp(right_primer.sequence) << endl;
                            for (int i = 0; i < right_primer.mapped_nodes_ids.size(); i++) {
                                size_t cur_node_id = right_primer.mapped_nodes_ids[i];
                                handle_t cur_node_handle = graph->get_handle(cur_node_id);
                                cout << graph->get_sequence(cur_node_handle) << " ";
                            }
                            cout << endl;
                            cout << endl;
                        }
                        cout << "load_primers works! :)" << endl;
                        cout << "-------------------------------------" << endl;

                }
                // void filter_primer() {

                //}
            
            private:
                // Functions only used in load_primers().. Not sure where to put them for now
                string rstrip(string const s) {
                    const string WHITESPACE = " \n\r\t\f\v";
                    size_t end = s.find_last_not_of(WHITESPACE);
                    if (end == string::npos) {
                            return "";
                    }
                    return s.substr(0, end+1);
                }

                size_t longest_match_len(Primer& primer, string const left_seq, string const right_seq) {
                    size_t llen = left_seq.size(), rlen = right_seq.size();
                    size_t length = min(llen, rlen);
                    size_t longest_match = 0;

                    // Change .. can be done in one for loop
                    if (llen >= rlen) {
                        for (size_t i = 0; i <= llen - rlen; i++) {
                            if (left_seq.substr(i, rlen) == right_seq) {
                                longest_match = rlen;
                                return longest_match;
                            }
                        }
                    }

                    for (size_t i = 1; i <= length; i++) {
                        if (left_seq.substr(llen - i, i) == right_seq.substr(0, i)) {
                            longest_match = i;
                        }
                    }

                    return longest_match;
                }

                char complement(char nt) {
                    switch(nt) {
                        case 'A': return 'T';
                        case 'C': return 'G';
                        case 'G': return 'C';
                        case 'T': return 'A';
                    }
                    return 'N';
                }

                string revcomp(string const seq) {
                    string revcomp_seq;
                    for (int i = seq.size()-1; i >= 0; i--) {
                        revcomp_seq += complement(seq[i]);
                    }
                    return revcomp_seq;
                }

                vector<string> split(string str, string const delim) {
                    size_t cur_pos = 0;
                    string word;
                    vector<string> word_list;
                    while ((cur_pos = str.find(delim)) != string::npos) {
                        word = str.substr(0, cur_pos);
                        word_list.push_back(word);
                        str.erase(0, cur_pos + delim.length());
                    }
                    word = str;
                    word_list.push_back(word);
                    return word_list;
                }
            };
            
            SnarlDistanceIndex distance_index;
            //unique_ptr<handlegraph::HandleGraph> graph;
            unique_ptr<handlegraph::PathPositionHandleGraph> graph;
            string snarl_index_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.dist";
            string xg_graph_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.xg";
            string primers_path = "/home/azhang/rotations/rotation_1/vg/alan/small/y.primer3.out";
            distance_index.deserialize(snarl_index_path);
            //graph = vg::io::VPKG::load_one<HandleGraph>(xg_graph_path);
            graph = vg::io::VPKG::load_one<PathPositionHandleGraph>(xg_graph_path);

            //Primer_finder primer_finder;
            Primer_finder primer_finder(graph, "y", &distance_index);
            primer_finder.load_primers(primers_path);
            primer_finder.run_test();
            



            // SnarlDistanceIndex distance_index;
            
            // string snarl_index_path = "/home/azhang/rotations/rotation_1/vg/alan/tiny/tiny.dist";
            // string xg_graph_path = "/home/azhang/rotations/rotation_1/vg/alan/tiny/tiny.xg";
            
            // distance_index.deserialize(snarl_index_path);
            // unique_ptr<HandleGraph> graph = vg::io::VPKG::load_one<HandleGraph>(xg_graph_path);

            // net_handle_t node = distance_index.get_node_net_handle(1);
            // cout << distance_index.net_handle_as_string(node) << endl;

            // net_handle_t root_node = distance_index.get_root();
            // cout << "root: " << distance_index.net_handle_as_string(root_node) << endl;
            
            // cout << distance_index.is_root(node) << endl;
            // cout << distance_index.is_snarl(node) << endl;
            // cout << distance_index.is_chain(node) << endl;
            // cout << distance_index.is_node(node) << endl;
            // cout << "--------------------------------------------" << endl;
            // cout << graph->get_node_count() << endl;
            // cout << graph->get_edge_count() << endl;
            // nid_t min_node_id = graph->min_node_id();
            // cout << graph->has_node(min_node_id) << endl;
            // nid_t max_node_id = graph->max_node_id();
            // cout << min_node_id << endl;
            // cout << max_node_id << endl;
            // cout << "--------------------------------------------" << endl;
            // handle_t min_node = graph->get_handle(min_node_id);
            // cout << graph->get_sequence(min_node) << endl;
            // handle_t max_node = graph->get_handle(max_node_id);
            // cout << graph->get_sequence(max_node) << endl;


        }

    }

}