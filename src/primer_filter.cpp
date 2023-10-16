#include "primer_filter.hpp"

namespace vg {

using namespace std;

Primer_finder::Primer_finder( unique_ptr<handlegraph::PathPositionHandleGraph>& graph_param,
                string reference_path_name, SnarlDistanceIndex* distance_index_param) {
    graph = graph_param.get();
    reference_path_handle = graph->get_path_handle("y");
    distance_index = distance_index_param;
}

Primer_finder::~Primer_finder() {
    // nothing to do
}

vector<Primer_pair> Primer_finder::get_primer_pairs() {
    return primer_pairs;
}

vector<Primer_pair> Primer_finder::get_selected_primer_pairs() {
    return selected_primer_pairs;
}

void Primer_finder::add_primer_pair(size_t left_primer_starting_node_id,
                    size_t left_primer_offset, size_t left_primer_length,
                    size_t right_primer_starting_node_id,
                    size_t right_primer_offset, size_t right_primer_length) {
    Primer left_primer;
    Primer right_primer;
    make_primer(left_primer, left_primer_starting_node_id, left_primer_offset, left_primer_length, true);
    make_primer(right_primer, right_primer_starting_node_id, right_primer_offset, right_primer_length, false);
    Primer_pair new_primer_pair {left_primer, right_primer,
        right_primer.position - left_primer.position + right_primer.length};
    update_min_max_product_size(new_primer_pair);
    primer_pairs.push_back(new_primer_pair);
    if (no_variation(new_primer_pair)) {
        selected_primer_pairs.push_back(new_primer_pair);
    }
}

void Primer_finder::load_primers(string path_to_primers) {
    regex left_seq_pattern("PRIMER_LEFT_\\d+_SEQUENCE=(\\w+)");
    regex right_seq_pattern("PRIMER_RIGHT_\\d+_SEQUENCE=(\\w+)");
    regex left_pos_pattern("PRIMER_LEFT_\\d+=(\\d+,\\d+)");
    regex right_pos_pattern("PRIMER_RIGHT_\\d+=(\\d+,\\d+)");
    
    Primer left_primer {""};
    Primer right_primer {"", false};

    ifstream file_handle(path_to_primers);
    assert(file_handle.is_open());
    
    string line;
    while (getline(file_handle, line)) {
        line = rstrip(line);
        smatch match;
        if (regex_search(line, match, left_seq_pattern)) {
            if (right_primer.sequence != "") {
                map_to_nodes(left_primer);
                map_to_nodes(right_primer);
                Primer_pair primer_pair {left_primer, right_primer,
                    right_primer.position - left_primer.position + right_primer.length};
                update_min_max_product_size(primer_pair);
                primer_pairs.push_back(primer_pair);
                if (no_variation(primer_pair)) {
                    selected_primer_pairs.push_back(primer_pair);
                }
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
            right_primer.position = stoi(pos_and_len[0]) - right_primer.length + 1;
            
        }
    }

    // Process and store the last pair of primers
    map_to_nodes(left_primer);
    map_to_nodes(right_primer);
    Primer_pair primer_pair {left_primer, right_primer, 
        right_primer.position - left_primer.position + right_primer.length};
    update_min_max_product_size(primer_pair);
    if (no_variation(primer_pair)) {
        selected_primer_pairs.push_back(primer_pair);
    }
    primer_pairs.push_back(primer_pair);
}

void Primer_finder::make_primer(Primer& primer, size_t starting_node_id, size_t offset, size_t length, bool is_left) {
    if (is_left) {
        primer.left = true;
    } else {
        primer.left = false;
    }
    primer.length = length;
    string sequence = "";
    handle_t cur_handle = graph->get_handle(starting_node_id);
    step_handle_t cur_step_handle = graph->steps_of_handle(cur_handle)[0];
    primer.position = graph->get_position_of_step(cur_step_handle) + offset;
    if (graph->get_length(cur_handle) - offset > length) {
        sequence += graph->get_sequence(cur_handle).substr(offset, length);
    } else {
        sequence += graph->get_sequence(cur_handle).substr(offset, graph->get_length(cur_handle) - offset);
        while (sequence.size() < length) {
            cur_step_handle = graph->get_next_step(cur_step_handle);
            cur_handle = graph->get_handle_of_step(cur_step_handle);
            sequence += graph->get_sequence(cur_handle).substr(0, min(graph->get_length(cur_handle), length-sequence.size()));
        }
    }

    if (is_left) {
        primer.sequence = sequence;
    } else {
        primer.sequence = revcomp(sequence);
    }
    map_to_nodes(primer);
}

void Primer_finder::update_min_max_product_size(Primer_pair& primer_pair) {
    Primer left_primer = primer_pair.left_primer;
    Primer right_primer = primer_pair.right_primer;
    
    primer_pair.min_product_size = distance_index->minimum_distance(left_primer.mapped_nodes_ids[0],
        false, left_primer.offset, right_primer.mapped_nodes_ids[right_primer.mapped_nodes_ids.size()-1], 
        false, right_primer.offset);

    primer_pair.max_product_size = distance_index->maximum_distance(left_primer.mapped_nodes_ids[0],
        false, left_primer.offset, right_primer.mapped_nodes_ids[right_primer.mapped_nodes_ids.size()-1], 
        false, right_primer.offset);
}

void Primer_finder::map_to_nodes(Primer& primer) {
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
    size_t primer_matched_index =  longest_match_len(primer, cur_node_sequence, primer_seq, true) - 1;
    while (primer_matched_index < primer_seq.size()-1) {
        cur_node_step_handle = graph->get_next_step(cur_node_step_handle);
        cur_node_handle = graph->get_handle_of_step(cur_node_step_handle);
        primer.mapped_nodes_ids.push_back(graph->get_id(cur_node_handle));
        cur_node_sequence = graph->get_sequence(cur_node_handle);
        string primer_substr = primer_seq.substr(primer_matched_index + 1, primer.length - primer_matched_index - 1);
        primer_matched_index += longest_match_len(primer, primer_substr, cur_node_sequence, false);
    }
}

size_t Primer_finder::longest_match_len(Primer& primer, string const left_seq, string const right_seq,
                        bool const first_node) {
    size_t llen = left_seq.size(), rlen = right_seq.size();
    size_t length = min(llen, rlen);
    size_t longest_match = 0;

    if (first_node) {
        if (llen >= rlen) {
            for (size_t i = 0; i <= llen - rlen; i++) {
                if (left_seq.substr(i, rlen) == right_seq) {
                    longest_match = rlen;
                    primer.offset = (primer.left) ? i : i + primer.sequence.size();
                    return longest_match;
                }
            }
        }
        for (size_t i = 1; i <= length; i++) {
            if (left_seq.substr(llen - i, i) == right_seq.substr(0, i)) {
                longest_match = i;
                primer.offset = (primer.left && first_node) ? llen - i : i;
            }
        }
    } else {
        for (size_t i = 1; i <= length; i++) {
            if (left_seq.substr(0, i) == right_seq.substr(0, i)) {
                longest_match = i;
                primer.offset = (!primer.left) ? i : primer.offset;
            }
        }
    }

    return longest_match;
}
        
string Primer_finder::rstrip(string const s) {
    const string WHITESPACE = " \n\r\t\f\v";
    size_t end = s.find_last_not_of(WHITESPACE);
    if (end == string::npos) {
            return "";
    }
    return s.substr(0, end+1);
}


bool Primer_finder::no_variation(const Primer_pair& primer_pair) {
    Primer left_primer = primer_pair.left_primer;
    Primer right_primer = primer_pair.right_primer; 
    for (vector<size_t>::iterator node_id = left_primer.mapped_nodes_ids.begin(); node_id != left_primer.mapped_nodes_ids.end(); ++node_id) {
        handle_t cur_handle = graph->get_handle(*node_id);
        net_handle_t cur_net_handle = distance_index->get_net(cur_handle, graph);
        size_t depth = distance_index->get_depth(cur_net_handle);
        if (depth != 1) {
            return false;
        }
    }
    return true;
}


char Primer_finder::complement(char nt) {
    switch(nt) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
    }
    return 'N';
}


string Primer_finder::revcomp(string const seq) {
    string revcomp_seq;
    for (int i = seq.size()-1; i >= 0; i--) {
        revcomp_seq += complement(seq[i]);
    }
    return revcomp_seq;
}


vector<string> Primer_finder::split(string str, string const delim) {
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

}