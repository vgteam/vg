#include "primer_filter.hpp"

namespace vg {

using namespace std;

// Constructor
PrimerFinder::PrimerFinder(const unique_ptr<handlegraph::PathPositionHandleGraph>& graph_param,
                const string& reference_path_name, const SnarlDistanceIndex* distance_index_param) {
    graph = graph_param.get();
    reference_path_handle = graph->get_path_handle(reference_path_name);
    distance_index = distance_index_param;
}

// Destructor
PrimerFinder::~PrimerFinder() {
    // nothing to do
}

const vector<PrimerPair>& PrimerFinder::get_primer_pairs() const {
    return primer_pairs;
}

const vector<PrimerPair>& PrimerFinder::get_selected_primer_pairs() const {
    return selected_primer_pairs;
}


// Make a new pair of primers with given attributes. Primers are processed and 
// added to primer_pairs and selected_primer_pairs.
void PrimerFinder::add_primer_pair(const size_t& left_primer_starting_node_id,
        const size_t& left_primer_offset, const size_t& left_primer_length,
        const size_t& right_primer_starting_node_id,
        const size_t& right_primer_offset, const size_t& right_primer_length) {
    
    primer_pairs.emplace_back();
    PrimerPair& primer_pair = primer_pairs.back();
    primer_pair.right_primer.left = false;

    make_primer(primer_pair.left_primer, left_primer_starting_node_id, left_primer_offset, left_primer_length, true);
    make_primer(primer_pair.right_primer, right_primer_starting_node_id, right_primer_offset, right_primer_length, false);
    primer_pair.linear_product_size = primer_pair.right_primer.position - primer_pair.left_primer.position + primer_pair.right_primer.length;
    update_min_max_product_size(primer_pair);
    if (no_variation(primer_pair)) {
        primer_pair.no_variation = true;
        selected_primer_pairs.push_back(primer_pairs.back());
    }

}

void PrimerFinder::load_primers(const string& path_to_primers) {

    // regular expression patterns to look for primers' sequences, positions on
    // the reference genome, and lengths
    regex left_seq_pattern("PRIMER_LEFT_\\d+_SEQUENCE=(\\w+)"); // e.g. PRIMER_LEFT_0_SEQUENCE=ACCGT
    regex right_seq_pattern("PRIMER_RIGHT_\\d+_SEQUENCE=(\\w+)");
    regex left_pos_pattern("PRIMER_LEFT_\\d+=(\\d+,\\d+)"); // e.g. PRIMER_LEFT_0_=125,20
    regex right_pos_pattern("PRIMER_RIGHT_\\d+=(\\d+,\\d+)");
    
    // iniate first primer pair
    primer_pairs.emplace_back();
    primer_pairs.back().right_primer.left = false;
    ifstream file_handle(path_to_primers);
    assert(file_handle.is_open());
    
    string line;
    while (getline(file_handle, line)) {
        line = rstrip(line);
        smatch match;
        if (regex_search(line, match, left_seq_pattern)) {
            if (primer_pairs.back().right_primer.sequence != "") {
                // primers' attributes are processed and stored into primer_pairs here
                map_to_nodes(primer_pairs.back().left_primer);
                map_to_nodes(primer_pairs.back().right_primer);
                primer_pairs.back().linear_product_size = primer_pairs.back().right_primer.position
                    - primer_pairs.back().left_primer.position + primer_pairs.back().right_primer.length;
                update_min_max_product_size(primer_pairs.back());
                if (no_variation(primer_pairs.back())) {
                    primer_pairs.back().no_variation = true;
                    selected_primer_pairs.push_back(primer_pairs.back());
                    PrimerPair& pp = primer_pairs.back();
                }
                primer_pairs.emplace_back();
                primer_pairs.back().right_primer.left = false;
            }
            primer_pairs.back().left_primer.sequence = match[1];
        } else if (regex_search(line, match, right_seq_pattern)) {
            primer_pairs.back().right_primer.sequence = match[1];
        } else if (regex_search(line, match, left_pos_pattern)) {
            const vector<string> pos_and_len = split(match[1], ",");
            primer_pairs.back().left_primer.position = stoi(pos_and_len[0]);
            primer_pairs.back().left_primer.length = stoi(pos_and_len[1]);
        } else if (regex_search(line, match, right_pos_pattern)) {
            const vector<string> pos_and_len = split(match[1], ",");
            primer_pairs.back().right_primer.length = stoi(pos_and_len[1]);
            primer_pairs.back().right_primer.position = stoi(pos_and_len[0]) - stoi(pos_and_len[1]) + 1;
        }
    }

    // Process and store the last pair of primers
    map_to_nodes(primer_pairs.back().left_primer);
    map_to_nodes(primer_pairs.back().right_primer);
    primer_pairs.back().linear_product_size = primer_pairs.back().right_primer.position
        - primer_pairs.back().left_primer.position + primer_pairs.back().right_primer.length;
    update_min_max_product_size(primer_pairs.back());
    if (no_variation(primer_pairs.back())) {
        primer_pairs.back().no_variation = true;
        selected_primer_pairs.push_back(primer_pairs.back());
    }
}

void PrimerFinder::make_primer(Primer& primer, const size_t& starting_node_id,
        const size_t& offset, const size_t& length, const bool& is_left) {
    if (is_left) {
        primer.left = true;
    } else {
        primer.left = false;
    }
    primer.length = length;
    string sequence = "";
    handle_t cur_handle = graph->get_handle(starting_node_id); // get the starting node handle
    step_handle_t cur_step_handle = graph->steps_of_handle(cur_handle)[0];
    primer.position = graph->get_position_of_step(cur_step_handle) + offset;
    // Walk down the path and get the sequence of primer
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
        primer.sequence = reverse_complement(sequence); // Take the reverse complement for right primer
    }
    map_to_nodes(primer); // Search and store corresponding nodes ids 
}

void PrimerFinder::update_min_max_product_size(PrimerPair& primer_pair) {
    const Primer& left_primer = primer_pair.left_primer;
    const Primer& right_primer = primer_pair.right_primer;
    
    primer_pair.min_product_size = distance_index->minimum_distance(left_primer.mapped_nodes_ids[0],
        false, left_primer.offset, right_primer.mapped_nodes_ids[right_primer.mapped_nodes_ids.size()-1], 
        false, right_primer.offset);

    primer_pair.max_product_size = distance_index->maximum_distance(left_primer.mapped_nodes_ids[0],
        false, left_primer.offset, right_primer.mapped_nodes_ids[right_primer.mapped_nodes_ids.size()-1], 
        false, right_primer.offset);
}

void PrimerFinder::map_to_nodes(Primer& primer) {
    string primer_seq;
    if (primer.left) {
        primer_seq = primer.sequence;
    } else {
        primer_seq = reverse_complement(primer.sequence);
    }

    step_handle_t  cur_node_step_handle = graph->get_step_at_position(reference_path_handle, primer.position);
    handle_t cur_node_handle = graph->get_handle_of_step(cur_node_step_handle);
    primer.mapped_nodes_ids.push_back(graph->get_id(cur_node_handle));
    string cur_node_sequence = graph->get_sequence(cur_node_handle);
    // Get the index at which primer.sequence[0:index] maps to the first node.
    // Stop here if the first node contains the entire primer sequence 
    size_t primer_matched_index =  longest_match_len(primer, cur_node_sequence, primer_seq, true) - 1;
    // If the first node containly a prefix of primer sequence, walk down the path and keep adding
    // node until the entire primer sequence is covered
    while (primer_matched_index < primer_seq.size()-1) {
        cur_node_step_handle = graph->get_next_step(cur_node_step_handle);
        cur_node_handle = graph->get_handle_of_step(cur_node_step_handle);
        primer.mapped_nodes_ids.push_back(graph->get_id(cur_node_handle));
        cur_node_sequence = graph->get_sequence(cur_node_handle);
        string primer_substr = primer_seq.substr(primer_matched_index + 1, primer.length - primer_matched_index - 1);
        primer_matched_index += longest_match_len(primer, primer_substr, cur_node_sequence, false);
    }
}

size_t PrimerFinder::longest_match_len(Primer& primer, const string& left_seq,
        const string& right_seq, const bool& first_node) {
    size_t llen = left_seq.size(), rlen = right_seq.size();
    size_t length = min(llen, rlen);
    size_t longest_match = 0;

    if (first_node) {
        if (llen >= rlen) {
            // Check if the first node contains the entire sequence of the priemr
            for (size_t i = 0; i <= llen - rlen; i++) {
                if (left_seq.substr(i, rlen) == right_seq) {
                    longest_match = rlen;
                    primer.offset = (primer.left) ? i : i + primer.sequence.size();
                    return longest_match;
                }
            }
        }
        for (size_t i = 1; i <= length; i++) {
            // Find the length of match between first node sequence's suffix and primer sequnece's prefix
            if (left_seq.substr(llen - i, i) == right_seq.substr(0, i)) {
                longest_match = i;
                primer.offset = (primer.left && first_node) ? llen - i : i;
            }
        }
    } else {
        for (size_t i = 1; i <= length; i++) {
            // Find the length of match between downstream nodes seqeunces and primer sequence
            if (left_seq.substr(0, i) == right_seq.substr(0, i)) {
                longest_match = i;
                primer.offset = (!primer.left) ? i : primer.offset;
            }
        }
    }

    return longest_match;
}
        
const string PrimerFinder::rstrip(const string& s) const {
    const string WHITESPACE = " \n\r\t\f\v";
    size_t end = s.find_last_not_of(WHITESPACE);
    if (end == string::npos) {
            return "";
    }
    return s.substr(0, end+1);
}


const bool PrimerFinder::no_variation(const PrimerPair& primer_pair) const {
    Primer left_primer = primer_pair.left_primer;
    Primer right_primer = primer_pair.right_primer; 
    for (vector<size_t>::iterator node_id = left_primer.mapped_nodes_ids.begin(); node_id != left_primer.mapped_nodes_ids.end(); ++node_id) {
        // Check if any node has depth more than 1 (i.e. inside a bubble)
        handle_t cur_handle = graph->get_handle(*node_id);
        net_handle_t cur_net_handle = distance_index->get_net(cur_handle, graph);
        size_t depth = distance_index->get_depth(cur_net_handle);
        if (depth != 1) {
            return false;
        }
    }
    return true;
}

const vector<string> PrimerFinder::split(string str, const string& delim) const {
    // Works like python split() function
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