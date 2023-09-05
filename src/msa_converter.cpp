/**
 * \file
 * msa_converter.cpp: contains a class that can construct VGs from clustal MSAs
 */


#include "vg.hpp"
#include "msa_converter.hpp"


//#define debug_msa_converter

namespace vg {

using namespace std;
    
    MSAConverter::MSAConverter() {
        // nothing to do
    }
    
    MSAConverter::~MSAConverter() {
        // nothing to do
    }

    void MSAConverter::load_alignments(istream& in, string format){
        
        create_progress("loading MSA into memory", 1);
        
        auto tokenize = [](string str) {
            string buf;
            vector<string> tokens;
            stringstream strm(str);
            while (strm >> buf) {
                tokens.push_back(buf);
            }
            return tokens;
        };
        
        if (format == "clustal") {
            
            unordered_set<char> conservation_chars{'.', ':', '*'};
            
            auto is_conservation_line = [&](string& line) {
                bool conservation_line = false;
                for (char c : line) {
                    if (!isspace(c)) {
                        if (conservation_chars.count(c)) {
                            conservation_line = true;
                        }
                        else {
                            conservation_line = false;
                        }
                        break;
                    }
                }
                return conservation_line;
            };
            
            auto is_blank = [](const string& str) {
                return all_of(str.begin(), str.end(), [](char c){return isspace(c);});
            };
            
            auto get_next_sequence_line = [&](istream& in) {
                string next;
                
                bool got_data = getline(in, next).good();
                bool conservation_line = is_conservation_line(next);
                bool blank = is_blank(next);
                
                while (got_data && (next.empty() || conservation_line || blank)) {
                    
                    got_data = getline(in, next).good();
                    conservation_line = is_conservation_line(next);
                    blank = is_blank(next);
                    
                }
                return next;
            };
            
            // make an alignment block
            alignments.emplace_back();
            auto& alignment = alignments.back();
            
            // skip the header line
            get_next_sequence_line(in);
            
            string line = get_next_sequence_line(in);
            while (!line.empty()) {
                vector<string> tokens = tokenize(line);
                
                if (tokens.size() != 2) {
                    continue;
                }
                
                auto iter = alignment.find(tokens[0]);
                if (iter != alignment.end()) {
                    iter->second.append(tokens[1]);
                }
                else {
                    alignment[tokens[0]] = tokens[1];
                }
                
                line = get_next_sequence_line(in);
            }
        }
        else if (format == "fasta") {
            
            // make an alignment block
            alignments.emplace_back();
            auto& alignment = alignments.back();
            
            string curr_seq_name;
            string line;
            
            bool got_data = getline(in, line).good();
            while (got_data && !line.empty()) {
                if (line[0] == '>') {
                    curr_seq_name = tokenize(line)[0].substr(1, line.size() - 1);
                    if (alignment.count(curr_seq_name)) {
                        cerr << "error:[MSAConverter] repeated sequence name '" << curr_seq_name << "' within an alignment, sequence names must be unique" << endl;
                        exit(1);
                    }
                }
                else {
                    alignment[curr_seq_name].append(line);
                }
                got_data = getline(in, line).good();
            }
        }
        else {
            cerr << "error:[MSAConverter] unsupported MSA format '" << format << "'" << endl;
            exit(1);
        }
        
        
#ifdef debug_msa_converter
        cerr << "alignments:" << endl;
        for (const auto& aln : alignments) {
            for (const auto& seq : aln) {
                cerr << seq.first << "\t" << seq.second << endl;
            }
            cerr << endl;
        }
#endif
        
        for (const auto& aln : alignments) {
            size_t aln_len = aln.begin()->second.size();
            for (const auto& seq : aln) {
                if (seq.second.size() != aln_len) {
                    cerr << "error:[MSAConverter] aligned sequences must be the same length, any unaligned sequences must be fully specified with '-' characters" << endl;
                    exit(1);
                }
            }
        }
        
        increment_progress();
        destroy_progress();
    }
    
    VG MSAConverter::make_graph(bool keep_paths, size_t max_node_length) {
        
        
        unordered_set<char> alphabet{'A', 'C', 'T', 'G', 'N', '-'};
        
        VG graph;
        
        // detect sequences with duplicate names and determine the size of the conversion
        bool contains_duplicate_seq_names = false;
        size_t total_size = 0;
        unordered_map<string, pair<size_t, size_t>> seq_name_count;
        for (const unordered_map<string, string>& alignment : alignments) {
            for (const pair<string, string>& sequence : alignment) {
                total_size += sequence.second.size();
                seq_name_count[sequence.first].first++;
                contains_duplicate_seq_names |= (seq_name_count[sequence.first].first > 1);
            }
        }
        
        create_progress("converting to graph", total_size);
        
        // append a number to each duplicate name
        if (contains_duplicate_seq_names) {
            for (unordered_map<string, string>& alignment : alignments) {
                
                // collect the names in this block that have duplicates
                vector<string> duplicate_names;
                for (const pair<string, string>& sequence : alignment) {
                    if (seq_name_count[sequence.first].first > 1) {
                        duplicate_names.push_back(sequence.first);
                    }
                }
                
                // TODO: it's technically possible that the new key might collide with one already
                // in the map, although this would have required some weird naming
                // swap the key in the map for the modified one
                for (string& seq_name : duplicate_names) {
                    size_t name_num = ++seq_name_count[seq_name].second;
                    stringstream sstrm;
                    sstrm << seq_name << "." << name_num;
                    alignment[sstrm.str()] = move(alignment[seq_name]);
                    alignment.erase(seq_name);
                }
            }
        }
        
        for (const unordered_map<string, string>& alignment : alignments) {
            
            // the node that each input sequence is extending
            unordered_map<string, Node*> current_node;
            
            // the path we're building for each aligned sequence
            unordered_map<string, Path*> aln_path;
            
            // start all of the alignments on a dummy node
            Node* dummy_node = graph.create_node("N");
            for (const auto& seq : alignment) {
                current_node[seq.first] = dummy_node;
                
                if (keep_paths) {
                    Path* path = graph.graph.add_path();
                    aln_path[seq.first] = path;
                    path->set_name(seq.first);
                }
            }
            
            // nodes that we don't want to extend any more
            // (we never want to extend the dummy node)
            unordered_set<Node*> completed_nodes{dummy_node};
            
            size_t aln_len = alignment.begin()->second.size();
            for (size_t i = 0; i < aln_len; i++) {
#ifdef debug_msa_converter
                cerr << "## beginning column " << i << endl;
#endif
                unordered_map<Node*, char> forward_transitions;
                unordered_map<char, pair<unordered_set<Node*>, vector<string>>> transitions;
                for (const auto& seq : alignment) {
                    char aln_char = toupper(seq.second[i]);
                    
                    if (!alphabet.count(aln_char)) {
                        cerr << "error:[MSAConverter] MSA contains non-nucleotide characters" << endl;
                        exit(1);
                        
                    }
                    
                    Node* node_here = current_node[seq.first];
                    
                    if (aln_char != '-') {
                        // this alignment is transitioning to a new aligned character
                        transitions[aln_char].first.insert(node_here);
                        transitions[aln_char].second.push_back(seq.first);
                        
                        auto iter = forward_transitions.find(node_here);
                        if (iter != forward_transitions.end()) {
                            if (iter->second != aln_char) {
                                // this node splits in the current column, so don't extend it anymore
                                completed_nodes.insert(node_here);
                            }
                        }
                        else {
                            forward_transitions[node_here] = aln_char;
                        }
                    }
                    else {
                        // this alignment isn't transitioning anywhere yet
                        
                        // we don't want to extend nodes where we'll need to attach a gap edge later
                        completed_nodes.insert(node_here);
                    }
                }
                
                for (const auto& transition : transitions) {
#ifdef debug_msa_converter
                    cerr << "transition to " << transition.first << endl;
                    cerr << "from nodes:" << endl;
                    for (const Node* n : transition.second.first) {
                        cerr << "\t" << n->id() << ": " << n->sequence() << endl;
                    }
                    cerr << "on sequences:" << endl;
                    for (const string& s : transition.second.second) {
                        cerr << "\t" << s << endl;
                    }
#endif
                    Node* at_node;
                    
                    if (transition.second.first.size() > 1) {
                        Node* new_node = graph.create_node(string(1, transition.first));
                        
                        for (Node* attaching_node : transition.second.first) {
                            graph.create_edge(attaching_node, new_node);
                            // we don't want to extend nodes that already have edges out of their ends
                            completed_nodes.insert(attaching_node);
                        }
                        
                        // keep track of the fact that now we are on the new node
                        at_node = new_node;
                        
                    }
                    else {
                        // there's only one node that wants to transition to this
                        // character, so we might be able to just extend the node
                        at_node = *transition.second.first.begin();
                        
                        if (at_node->sequence().size() >= max_node_length ||
                            completed_nodes.count(at_node)) {
                            // we either want to split this node just because of length or because
                            // we've already marked it as unextendable
                            
                            Node* new_node = graph.create_node(string(1, transition.first));
                            
                            graph.create_edge(at_node, new_node);
                            completed_nodes.insert(at_node);
                            
                            // keep track of the fact that now we are on the new node
                            at_node = new_node;
                        }
                        else {
                            at_node->mutable_sequence()->append(1, transition.first);
                        }
                    }
                    
                    // update which node the paths are currently extending
                    for (const string& name : transition.second.second) {
                        current_node[name] = at_node;
                        
                        if (keep_paths) {
                            Path* path = aln_path[name];
                            if (path->mapping_size() == 0 ? true :
                                at_node->id() != path->mapping(path->mapping_size() - 1).position().node_id()) {
                                Mapping* mapping = path->add_mapping();
                                mapping->mutable_position()->set_node_id(at_node->id());
                                mapping->set_rank(path->mapping_size());
                            }
                        }
                    }
                }
                
                update_progress(alignment.size());
                
#ifdef debug_msa_converter
                cerr << "graph: " << pb2json(graph.graph) << endl;
                cerr << "node locations of sequences:" << endl;
                for (const auto& curr : current_node) {
                    cerr << "\t" << curr.first << " " << curr.second->id() << endl;
                }
#endif
            }
            
            for (pair<const string, Path*> path_record : aln_path) {
                Path* path = path_record.second;
                for (size_t i = 0; i < path->mapping_size(); i++) {
                    Mapping* mapping = path->mutable_mapping(i);
                    assert(mapping->edit_size() == 0);
                    Edit* edit = mapping->add_edit();
                    
                    size_t node_length = graph.get_node(mapping->position().node_id())->sequence().size();
                    edit->set_from_length(node_length);
                    edit->set_to_length(node_length);
                }
                graph.paths.extend(*path);
            }
            
            graph.destroy_node(dummy_node);
            graph.sync_paths();
        }
        
        destroy_progress();
        
        return graph;
    }

}


