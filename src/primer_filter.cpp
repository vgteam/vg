#include "primer_filter.hpp"
#include <unordered_set>

namespace vg {

using namespace std;

// Constructor
PrimerFinder::PrimerFinder(const unique_ptr<handlegraph::PathPositionHandleGraph>& graph_param,
    const SnarlDistanceIndex* distance_index_param, ifstream& primers_file_handle,
    const gbwtgraph::GBWTGraph& gbwt_graph_param, const gbwt::GBWT& gbwt_index_param,
    const gbwt::FastLocate& r_index_param)
    : graph(graph_param.get()),
      distance_index(distance_index_param),
      gbwt_graph(gbwt_graph_param),
      gbwt_index(gbwt_index_param),
      r_index(r_index_param) {

    load_primers(primers_file_handle);
}

// Destructor
PrimerFinder::~PrimerFinder() {
    // nothing to do
}

const vector<PrimerPair>& PrimerFinder::get_primer_pairs_of_chrom(const string& chrom_name) const {
    return chroms.at(chrom_name);
}

// Make a new pair of primers with given attributes. Primers are processed and 
// added to primer_pairs and selected_primer_pairs.
void PrimerFinder::add_primer_pair(const string& path_name,
    const size_t& left_primer_starting_node_id, const size_t& left_primer_offset,
    const size_t& left_primer_length, const size_t& right_primer_starting_node_id,
    const size_t& right_primer_offset, const size_t& right_primer_length) {

    chroms.at(path_name).emplace_back();
    PrimerPair& primer_pair = chroms.at(path_name).back();
    primer_pair.chromosome_name   = path_name;
    primer_pair.template_position = 0;
    primer_pair.right_primer.left = false;

    make_primer(primer_pair.left_primer, path_name, left_primer_starting_node_id,
        left_primer_offset, left_primer_length, true);
    make_primer(primer_pair.right_primer, path_name, right_primer_starting_node_id,
        right_primer_offset, right_primer_length, false);
    primer_pair.linear_product_size = primer_pair.right_primer.position_template
        - primer_pair.left_primer.position_template + primer_pair.right_primer.length;
    update_variation(primer_pair, path_name);
    update_min_max_product_size(primer_pair);
}

void PrimerFinder::load_primers(ifstream& file_handle) {

    //ifstream file_handle(path_to_primers);
    assert(file_handle.is_open());
    
    vector<string> cur_fields;
    size_t cur_template_offset;
    string cur_template_info;
    string cur_template_feature;
    string cur_path;
    string line;
    while (getline(file_handle, line)) {
        line = strip(line);
        size_t left_primer_line_start  = line.find("LEFT PRIMER");
        size_t right_primer_line_start = line.find("RIGHT PRIMER");
        
        if (startswith(line, "PRIMER PICKING RESULTS FOR")) {
            if (chroms.size() != 0) {
                assert(chroms[cur_path].back().right_primer.sequence.empty());
                chroms[cur_path].pop_back();
            }
            cur_fields           = move(split(line));
            cur_template_info    = cur_fields[cur_fields.size()-1];
            cur_fields           = move(split(cur_template_info,'|'));
            cur_template_feature = cur_fields[1] + "|" + cur_fields[2];
            cur_template_offset  = stoi(cur_fields[3]);
            cur_path             = cur_fields[0];
            chroms[cur_path].emplace_back();
            chroms[cur_path].back().chromosome_name   = cur_path;
            chroms[cur_path].back().template_position = cur_template_offset;
            chroms[cur_path].back().template_feature  = cur_template_feature;
            chroms[cur_path].back().right_primer.left = false;
        } else if (left_primer_line_start != string::npos) {
            cur_fields = move(split(line.substr(left_primer_line_start, line.size())));
            PrimerPair& primer_pair = chroms[cur_path].back();
            primer_pair.left_primer.position_chromosome = stoi(cur_fields[2]) + cur_template_offset;
            primer_pair.left_primer.position_template   = stoi(cur_fields[2]);
            primer_pair.left_primer.sequence            = cur_fields[9];
            primer_pair.left_primer.length              = stoi(cur_fields[3]);
        } else if (startswith(line, "RIGHT PRIMER")) {
            cur_fields = move(split(line.substr(right_primer_line_start, line.size())));
            PrimerPair& primer_pair = chroms[cur_path].back();
            primer_pair.right_primer.position_chromosome = stoi(cur_fields[2]) - stoi(cur_fields[3]) + 1 + cur_template_offset;
            primer_pair.right_primer.position_template   = stoi(cur_fields[2]) - stoi(cur_fields[3]) + 1;
            primer_pair.right_primer.sequence            = cur_fields[9];
            primer_pair.right_primer.length              = stoi(cur_fields[3]);

            assert(!primer_pair.left_primer.sequence.empty());
            map_to_nodes(primer_pair.left_primer, cur_path);
            map_to_nodes(primer_pair.right_primer, cur_path);
            primer_pair.linear_product_size = primer_pair.right_primer.position_template
                - primer_pair.left_primer.position_template + primer_pair.right_primer.length;
            update_variation(primer_pair, cur_path);
            update_min_max_product_size(primer_pair);
            chroms[cur_path].emplace_back();
            chroms[cur_path].back().chromosome_name   = cur_path;
            chroms[cur_path].back().template_position = cur_template_offset;
            chroms[cur_path].back().template_feature  = cur_template_feature;
            chroms[cur_path].back().right_primer.left = false;
        }
    }
    assert(chroms[cur_path].back().right_primer.sequence.empty());
    chroms[cur_path].pop_back();
}

const size_t PrimerFinder::total_reference_paths() const {
    return chroms.size();
}

vector<string> PrimerFinder::get_reference_paths() {
    vector<string> reference_paths;
    for (const auto& chrom : chroms) {
        reference_paths.push_back(chrom.first);
    }
    return reference_paths;
}

void PrimerFinder::make_primer(Primer& primer, const string& path_name,
    const size_t& starting_node_id, const size_t& offset, const size_t& length,
    const bool& is_left) {
    
    if (is_left) {
        primer.left = true;
    } else {
        primer.left = false;
    }
    primer.length = length;
    string sequence = "";
    handle_t cur_handle = graph->get_handle(starting_node_id); // get the starting node handle
    step_handle_t cur_step_handle = graph->steps_of_handle(cur_handle)[0];
    primer.position_template   = graph->get_position_of_step(cur_step_handle) + offset;
    primer.position_chromosome = primer.position_template;
    // Walk down the path and get the sequence of primer
    if (graph->get_length(cur_handle) - offset > length) {
        sequence += graph->get_sequence(cur_handle).substr(offset, length);
    } else {
        sequence += graph->get_sequence(cur_handle).substr(offset, graph->get_length(cur_handle) - offset);
        while (sequence.size() < length) {
            cur_step_handle = graph->get_next_step(cur_step_handle);
            cur_handle      = graph->get_handle_of_step(cur_step_handle);
            sequence       += graph->get_sequence(cur_handle).substr(0, min(graph->get_length(cur_handle), length-sequence.size()));
        }
    }

    if (is_left) {
        primer.sequence = sequence;
    } else {
        primer.sequence = reverse_complement(sequence); // Take the reverse complement for right primer
    }
    map_to_nodes(primer, path_name); // Search and store corresponding nodes ids 
}

static string get_haplotype_sequence(gbwt::size_type sequence_visit_offset, handle_t start_handle,
    handle_t end_handle, const gbwtgraph::GBWTGraph& gbwt_graph, size_t start_max, size_t end_max) {

    string haplotype;
    gbwt::edge_type pos = gbwt::edge_type(gbwtgraph::GBWTGraph::handle_to_node(start_handle), sequence_visit_offset);
    
    if (pos == gbwt::invalid_edge() || pos.first == gbwt::ENDMARKER) {
        return haplotype;
    }
    
    handle_t curr = gbwt_graph.node_to_handle(pos.first);
    if (curr == end_handle) {
        return haplotype;
    }
    gbwtgraph::view_type view = gbwt_graph.get_sequence_view(curr);
    size_t offset = (view.second > start_max ? view.second - start_max : 0);
    haplotype.append(view.first + offset, view.second - offset);

    while (true) {
        pos = gbwt_graph.index->LF(pos);
        if (pos.first == gbwt::ENDMARKER) {
                break;
        }
        curr = gbwtgraph::GBWTGraph::node_to_handle(pos.first);
        view = gbwt_graph.get_sequence_view(curr);
        if (curr == end_handle) {
            haplotype.append(view.first, std::min(view.second, end_max));
            break;
        } else {
            haplotype.append(view.first, view.second);
        }
    }
    return haplotype;
}

void PrimerFinder::update_min_max_product_size(PrimerPair& primer_pair) {
    const auto& sequence_visits = primer_pair.sequence_visits;

    handle_t start_handle = gbwt_graph.get_handle(primer_pair.left_primer.mapped_nodes_ids.front());
    handle_t end_handle   = gbwt_graph.get_handle(primer_pair.right_primer.mapped_nodes_ids.back());
    if (start_handle == end_handle) {
        primer_pair.min_product_size = primer_pair.linear_product_size;
        primer_pair.max_product_size = primer_pair.linear_product_size;
        return;
    }
    
    size_t start_max = gbwt_graph.get_length(start_handle) - primer_pair.left_primer.offset;
    size_t end_max   = primer_pair.right_primer.offset;
    size_t minimum_distance = numeric_limits<size_t>::max();
    size_t maximum_distance = 0;
    for (const auto& visit : sequence_visits) {
        string haplotype = get_haplotype_sequence(visit.second, start_handle, end_handle, gbwt_graph, start_max, end_max);
        if (haplotype.size() < minimum_distance) {
            minimum_distance = haplotype.size();
        }
        if (haplotype.size() > maximum_distance) {
            maximum_distance = haplotype.size();
        }
    }
    primer_pair.min_product_size = minimum_distance;
    primer_pair.max_product_size = maximum_distance;
}

void PrimerFinder::map_to_nodes(Primer& primer, const string& path_name) {
    path_handle_t reference_path_handle = graph->get_path_handle(path_name);
    string primer_seq;
    if (primer.left) {
        primer_seq = primer.sequence;
    } else {
        primer_seq = reverse_complement(primer.sequence);
    }

    step_handle_t  cur_node_step_handle = graph->get_step_at_position(reference_path_handle, primer.position_chromosome);
    handle_t cur_node_handle = graph->get_handle_of_step(cur_node_step_handle);
    size_t cur_node_length   = graph->get_length(cur_node_handle);
    size_t cur_node_position = graph->get_position_of_step(cur_node_step_handle);
    size_t cur_offset = primer.position_chromosome - cur_node_position;
    primer.mapped_nodes_ids.push_back(graph->get_id(cur_node_handle));
    if (primer.left) {
        primer.offset = cur_offset;
    }
    size_t matched_length = 0;
    while (cur_node_length - cur_offset < primer.length - matched_length) {
        assert(graph->get_sequence(cur_node_handle).substr(cur_offset, cur_node_length - cur_offset)
            == primer_seq.substr(matched_length, cur_node_length - cur_offset));
        matched_length += cur_node_length - cur_offset;
        cur_offset = 0;
        cur_node_step_handle = graph->get_next_step(cur_node_step_handle);
        cur_node_handle = graph->get_handle_of_step(cur_node_step_handle);
        cur_node_length   = graph->get_length(cur_node_handle);
        primer.mapped_nodes_ids.push_back(graph->get_id(cur_node_handle));
    }
    assert(graph->get_sequence(cur_node_handle).substr(cur_offset, primer.length - matched_length)
        == primer_seq.substr(matched_length, primer.length - matched_length));
    if (!primer.left) {
        primer.offset = cur_offset + primer.length - matched_length;
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
        
const string PrimerFinder::strip(const string& s) const {
    const string WHITESPACE = " \n\r\t\f\v";
    size_t end   = s.find_last_not_of(WHITESPACE);
    size_t start = s.find_first_not_of(WHITESPACE); 
    if (end == string::npos) {
            return "";
    }
    return s.substr(start, end+1);
}

vector<HaplotypePartitioner::sequence_type> get_sequence_visits(const handle_t& handle,
    const gbwt::FastLocate& r_index, const gbwtgraph::GBWTGraph& gbwt_graph) {

    vector<gbwt::size_type> sa = r_index.decompressSA(gbwt_graph.handle_to_node(handle));
    vector<HaplotypePartitioner::sequence_type> result;
    result.reserve(sa.size());
    for (size_t i = 0; i < sa.size(); i++) {
        result.push_back({ sa[i], i });
    }
    std::sort(result.begin(), result.end(), [&](HaplotypePartitioner::sequence_type a, HaplotypePartitioner::sequence_type b) -> bool {
        gbwt::size_type a_id = r_index.seqId(a.first);
        gbwt::size_type a_offset = r_index.seqOffset(a.first);
        gbwt::size_type b_id = r_index.seqId(b.first);
        gbwt::size_type b_offset = r_index.seqOffset(b.first);
        return ((a_id < b_id) || ((a_id == b_id) && (a_offset > b_offset)));
    });
    return result;
}

static void sa_to_da(std::vector<HaplotypePartitioner::sequence_type>& sequences, const gbwt::FastLocate& r_index) {
    for (auto& sequence : sequences) {
        sequence.first = r_index.seqId(sequence.first);
    }
}

void PrimerFinder::update_variation(PrimerPair& primer_pair, const string& path_name) {
    const vector<size_t>& left_primer_node_ids  = primer_pair.left_primer.mapped_nodes_ids;
    const vector<size_t>& right_primer_node_ids = primer_pair.right_primer.mapped_nodes_ids;
    vector<size_t> nodes_id;
    merge(left_primer_node_ids.begin(),  left_primer_node_ids.end(),
          right_primer_node_ids.begin(), right_primer_node_ids.end(), back_inserter(nodes_id));
    handle_t cur_handle = gbwt_graph.get_handle(nodes_id[0]);
    auto sequence_visits = get_sequence_visits(cur_handle, r_index, gbwt_graph);
    sa_to_da(sequence_visits, r_index);

    for (size_t i = 1; i < nodes_id.size(); i++) {
        cur_handle = gbwt_graph.get_handle(nodes_id[i]);
        auto cur_sequence_visits = get_sequence_visits(cur_handle, r_index, gbwt_graph);
        sa_to_da(cur_sequence_visits, r_index);
        unordered_set<gbwt::size_type> seq_ids;
        for (const auto& seq_visit : cur_sequence_visits) {
            seq_ids.insert(seq_visit.first);
        }

        sequence_visits.erase(
            remove_if(
                sequence_visits.begin(), 
                sequence_visits.end(),
                [&seq_ids](const HaplotypePartitioner::sequence_type& visit) {
                    return seq_ids.find(visit.first) == seq_ids.end();
                }
            ),
            sequence_visits.end()
        );
    }

    auto unique_haplotypes = sequence_visits;
    auto it = unique(unique_haplotypes.begin(), unique_haplotypes.end(), [this](const auto& a, const auto& b) {
        const gbwt::PathName& path_name_a = this->gbwt_graph.index->metadata.path(gbwt::Path::id(a.first));
        const gbwt::PathName& path_name_b = this->gbwt_graph.index->metadata.path(gbwt::Path::id(b.first));
        return (path_name_a.sample == path_name_b.sample) && (path_name_a.phase == path_name_b.phase);
    });
    unique_haplotypes.erase(it, unique_haplotypes.end());

    primer_pair.sequence_visits = sequence_visits;
    primer_pair.variation_level = static_cast<double>(unique_haplotypes.size()) / static_cast<double>(gbwt_graph.index->metadata.haplotypes());

}

vector<string> PrimerFinder::split(const string& str) {
    istringstream iss(str);
    string field;
    vector<string> fields;

    while (iss >> field) {
        fields.push_back(field);
    }
    return fields;
}

vector<string> PrimerFinder::split(const string& str, const char& delim) {
    istringstream iss(str);
    string field;
    vector<string> fields;

    while (getline(iss, field, delim)) {
        fields.push_back(field);
    }

    return fields;
}


bool PrimerFinder::startswith(const string& str, const string& prefix) {
    return str.compare(0, prefix.length(), prefix) == 0;
}

}