#include "primer_filter.hpp"
#include <unordered_set>
#include "surjector.hpp"
#include "hts_alignment_emitter.hpp"

namespace vg {

using namespace std;

//#define DEBUG_PRIMER_FILTER

// Constructor
PrimerFinder::PrimerFinder(const handlegraph::PathPositionHandleGraph* graph_param,
    const SnarlDistanceIndex* distance_index_param, istream& primers_file_handle,
    const gbwtgraph::GBWTGraph& gbwt_graph_param, const gbwt::GBWT& gbwt_index_param,
    const gbwt::FastLocate& r_index_param, MinimizerMapper* giraffe_mapper_param)
    : graph(graph_param),
      distance_index(distance_index_param),
      gbwt_graph(gbwt_graph_param),
      gbwt_index(gbwt_index_param),
      r_index(r_index_param),
      giraffe_mapper(giraffe_mapper_param) {

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

void PrimerFinder::load_primers(istream& file_handle) {

    //ifstream file_handle(path_to_primers);
    assert(file_handle);
    

    // Regular expressions for matching fields with numbers
    std::regex left_seq ("^PRIMER_LEFT_[0-9]*_SEQUENCE.*");
    std::regex right_seq ("^PRIMER_RIGHT_[0-9]*_SEQUENCE.*");
    std::regex left_primer_position ("^PRIMER_LEFT_[0-9]*=.*");
    std::regex right_primer_position ("^PRIMER_RIGHT_[0-9]*=.*");

    vector<PrimerPair>::iterator curr_primer_iterator;

    string chromosome_name = "";
    string template_feature = "";
    size_t template_position = std::numeric_limits<size_t>::max();
    bool has_path = false;

    string line;
    while (getline(file_handle, line)) {
        line = strip(line);

        if (line == "=") {
            //End of the record for one primer pair
            chromosome_name = "";
            template_feature = "";
            template_position = std::numeric_limits<size_t>::max();
            has_path = false;
        } else if (startswith(line, "SEQUENCE_ID")) {
            //Get the path, path offset, and features from the sequence_id of the primer pair
            //This will be the same for all primer pairs up to the next "="
            vector<string> cur_fields = move(split(split(line,'=')[1], '|'));

            if (cur_fields.size() == 4) {
                //If the sequence id is correctly formatted
                chromosome_name = cur_fields[0];
                template_feature = cur_fields[1] + "|" + cur_fields[2];
                template_position = stoi(cur_fields[3]);
                has_path = graph->has_path(chromosome_name);
                if (!has_path) {
                    cerr << "warning: primer finder can't find a path named " << chromosome_name << " in the graph" << endl << "\tfalling back on mapping the template sequence" << endl;
                }
            } else {
                template_feature = line;
                has_path = false;
                cerr << "warning: primer finder " << line << " is not formatted with a path and offset" << endl << "\tfalling back on mapping the template sequence" << endl;
            }
#ifdef DEBUG_PRIMER_FILTER
            cerr << "FIND PRIMERS FOR INPUT " << line << ": " << chromosome_name << ", " << template_feature << ", " << template_position << endl;
#endif

        } else if (startswith(line, "SEQUENCE_TEMPLATE") && !has_path) {
            //If the path from the sequence id isn't in the graph, then get the path and path offset by mapping the sequence
            string seq = split(line,'=')[1];
            if (giraffe_mapper == nullptr) {
                throw std::runtime_error("error: primer filter doesn't have a minimizer file to map the template");
            }
            std::tie(chromosome_name, template_position) = get_graph_coordinates_from_sequence(seq);

        } else if (startswith(line, "PRIMER_PAIR_NUM_RETURNED")) {
            //How many primer pairs for this sequence template?

            size_t primer_pair_count = stoi(split(line,'=')[1]);
            size_t new_vector_start = chroms[chromosome_name].size();

            //Add all new primer pairs for this template
            chroms.reserve(new_vector_start + primer_pair_count);
            for (size_t i = 0 ; i < primer_pair_count ; i++) {
                chroms[chromosome_name].emplace_back();
                chroms[chromosome_name].back().chromosome_name   = chromosome_name;
                chroms[chromosome_name].back().template_position = template_position;
                chroms[chromosome_name].back().template_feature  = template_feature;
                chroms[chromosome_name].back().right_primer.left = false;
            }

            //Set the current primer pair iterator to the first new pair
            curr_primer_iterator = chroms[chromosome_name].begin() + new_vector_start; 
        } else if (std::regex_match(line, left_seq)) {
            curr_primer_iterator->left_primer.sequence = split(line, '=')[1];
#ifdef DEBUG_PRIMER_FILTER
            cerr << "\tGet left sequence " << line << ": " << curr_primer_iterator->left_primer.sequence << endl;
#endif
        } else if (std::regex_match(line, right_seq)) {
            curr_primer_iterator->right_primer.sequence = split(line, '=')[1];
#ifdef DEBUG_PRIMER_FILTER
            cerr << "\tGet right sequence " << line << ": " << curr_primer_iterator->left_primer.sequence << endl;
#endif
        } else if (std::regex_match(line, left_primer_position)) {
            //Start position and length of the left primer
            curr_primer_iterator->left_primer.position_template   = stoi(split(split(line, '=')[1], ',')[0]);
            curr_primer_iterator->left_primer.length              = stoi(split(split(line, '=')[1], ',')[1]);
            curr_primer_iterator->left_primer.position_chromosome = curr_primer_iterator->left_primer.position_template + template_position;
#ifdef DEBUG_PRIMER_FILTER
            cerr << "old template position " << template_position << endl;
            cerr << "\tGet left primer position" << line << ": " << curr_primer_iterator->left_primer.position_template << ", " 
                                                                 << curr_primer_iterator->left_primer.length << ", " 
                                                                 << curr_primer_iterator->left_primer.position_chromosome << endl;
#endif
        } else if (std::regex_match(line, right_primer_position)) {
#ifdef DEBUG_PRIMER_FILTER
            cerr << "\tGet right primer position" << line << ": " <<  curr_primer_iterator->left_primer.position_chromosome << endl;
#endif
            //Start position and length of the right primer
            size_t right_primer_offset = stoi(split(split(line, '=')[1], ',')[0]);
            curr_primer_iterator->right_primer.length              = stoi(split(split(line, '=')[1], ',')[1]);
            curr_primer_iterator->right_primer.position_chromosome = right_primer_offset - curr_primer_iterator->right_primer.length + 1 + template_position;
            curr_primer_iterator->right_primer.position_template   = right_primer_offset - curr_primer_iterator->right_primer.length + 1;

            //This is the last thing for this primer pair, so update the primer pair
            map_to_nodes(curr_primer_iterator->left_primer, chromosome_name);
            map_to_nodes(curr_primer_iterator->right_primer, chromosome_name);

            curr_primer_iterator->linear_product_size = curr_primer_iterator->right_primer.position_template
                - curr_primer_iterator->left_primer.position_template + curr_primer_iterator->right_primer.length;
            update_variation(*curr_primer_iterator, chromosome_name);
            update_min_max_product_size(*curr_primer_iterator);

            //Iterator to the new primer pair
            curr_primer_iterator++;
        }

    }
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
        std::cerr << "Only use rest of " << graph->get_id(cur_handle) << " to get " << sequence << std::endl;
    } else {
        sequence += graph->get_sequence(cur_handle).substr(offset, graph->get_length(cur_handle) - offset);
        std::cerr << "Start with rest of " << graph->get_id(cur_handle) << " to get " << sequence << std::endl;
        while (sequence.size() < length) {
            cur_step_handle = graph->get_next_step(cur_step_handle);
            cur_handle      = graph->get_handle_of_step(cur_step_handle);
            sequence       += graph->get_sequence(cur_handle).substr(0, min(graph->get_length(cur_handle), length-sequence.size()));
            std::cerr << "Add at least part of " << graph->get_id(cur_handle) << " to get " << sequence << std::endl;
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

std::pair<string, size_t> PrimerFinder::get_graph_coordinates_from_sequence(const string& seq) {
    string ref_name;
    int64_t ref_offset;
    bool ref_rev;

    //Make an alignment from the sequence
    Alignment aln;
    aln.set_sequence(seq);
    aln.set_name("primer_template");

    //Map the alignment
    vector<Alignment> mapped = giraffe_mapper->map(aln);

    //If there wasn't an alignment, error
    if (mapped.empty() || mapped.front().mapping_quality() < 30) {
        cerr << "error: Primer filter could not map template sequence " << seq << endl;
        return std::make_pair(ref_name, std::numeric_limits<size_t>::max());
    }


    //Get the reference paths we want to align to
    //This is done automatically
    //TODO: These are empty but they could be command line arguments
    string path_file;
    vector<string> path_names;
    vector<tuple<path_handle_t, size_t, size_t>> sequence_dictionary = get_sequence_dictionary(path_file, path_names, *graph);
    unordered_set<path_handle_t> reference_paths;
    reference_paths.reserve(sequence_dictionary.size());
    for (auto& entry : sequence_dictionary) {
        reference_paths.insert(get<0>(entry));
    }

    //Surject the alignment onto the reference paths
    Surjector surjector(graph);
    surjector.surject(mapped.front(), reference_paths, ref_name, ref_offset, ref_rev);

    //TODO: Double check that this is correct. idk why ref_offset is an int and not a size_t
    if (ref_rev) {
        ref_offset -= seq.size();
    }
    assert (graph->has_path(ref_name));
#ifdef DEBUG_PRIMER_FILTER
    cerr << "\tmapped sequence to " << ref_name << " at offset " << ref_offset << endl;
#endif

    return std::make_pair(ref_name, (size_t)ref_offset);
}


void PrimerFinder::update_min_max_product_size(PrimerPair& primer_pair) {
    if (primer_pair.chromosome_name.empty()) {
        return;
    }
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
#ifdef DEBUG_PRIMER_FILTER
    cerr << "Map to nodes for primer " << primer.sequence << endl;
#endif
    if (path_name.empty()) {
        return;
    }
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
        std::string graph_sequence = graph->get_sequence(cur_node_handle).substr(cur_offset, cur_node_length - cur_offset);
        std::string primer_sequence = primer_seq.substr(matched_length, cur_node_length - cur_offset);
        if (graph_sequence != primer_sequence) {
            std::stringstream ss;
            ss << "Graph sequence " << graph_sequence
                << " at node " << graph->get_id(cur_node_handle) << (graph->get_is_reverse(cur_node_handle) ? "-" : "+")
                << "(" << graph->get_sequence(cur_node_handle) << ")"
                << " range " << cur_offset << "-" << cur_node_length
                << " did not match primer " << primer.sequence << " sequence " + primer_sequence
                << " at range " << matched_length << "-" << (matched_length + cur_node_length - cur_offset);
            throw std::runtime_error(ss.str());
        }
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
#ifdef DEBUG_PRIMER_FILTER
    cerr << "Update variation" << endl;
#endif
    if (path_name.empty()) {
        return;
    }
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
