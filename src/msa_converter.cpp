/**
 * \file
 * msa_converter.cpp: contains a class that can construct VGs from clustal MSAs
 */


#include "vg.hpp"
#include "msa_converter.hpp"


//#define debug_msa_converter

namespace vg {

using namespace std;

    MSAConverter::MSAConverter(istream& in, string format, size_t max_node_length) : max_node_length(max_node_length) {
        
        
        auto tokenize = [](string str) {
            string buf;
            vector<string> tokens;
            stringstream strm(str);
            while (strm >> buf) {
                tokens.push_back(buf);
            }
            return tokens;
        };
        
        if (format == "maf") {
            auto get_next_sequence_line = [](istream& in) {
                string next;
                
                bool got_data = getline(in, next).good();
                while (got_data && (next.empty() ? true : next[0] != 's')) {
                    got_data = getline(in, next).good();
                }
                return next;
            };
            
            string line = get_next_sequence_line(in);
            
            while (!line.empty()) {
                vector<string> tokens = tokenize(line);
                
                assert(tokens.size() >= 7);
                
                auto iter = alignments.find(tokens[1]);
                if (iter != alignments.end()) {
                    iter->second.append(tokens[6]);
                }
                else {
                    alignments[tokens[1]] = tokens[6];
                }
                
                line = get_next_sequence_line(in);
            }
        }
        else if (format == "clustal") {
            
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
            
            auto get_next_sequence_line = [&](istream& in) {
                string next;
                
                bool got_data = getline(in, next).good();
                bool conservation_line = is_conservation_line(next);
                
                while (got_data && (next.empty() || conservation_line)) {
                    
                    got_data = getline(in, next).good();
                    conservation_line = is_conservation_line(next);
                    
                }
                if (conservation_line) {
                    // hack for edge case that the final line is a conservation line
                    next.clear();
                }
                return next;
            };
            
            // skip the header line
            get_next_sequence_line(in);
            
            string line = get_next_sequence_line(in);
            while (!line.empty()) {
                vector<string> tokens = tokenize(line);
                
                if (tokens.size() != 2) {
                    continue;
                }
                
                auto iter = alignments.find(tokens[0]);
                if (iter != alignments.end()) {
                    iter->second.append(tokens[1]);
                }
                else {
                    alignments[tokens[0]] = tokens[1];
                }
                
                
                line = get_next_sequence_line(in);
            }
        }
        else {
            cerr << "error:[MSAConverter] unsupported MSA format" << endl;
            exit(1);
        }
        
        
#ifdef debug_msa_converter
        cerr << "alignments:" << endl;
        for (const auto& aln : alignments) {
            cerr << aln.first << "\t" << aln.second << endl;
        }
#endif
        
        size_t aln_len = alignments.begin()->second.size();
        for (const auto& aln : alignments) {
            assert(aln.second.size() == aln_len);
        }
    }
    
    MSAConverter::~MSAConverter() {
        // nothing to do
    }
    
    VG MSAConverter::make_graph(bool keep_paths) {
        
        unordered_set<char> alphabet{'A', 'C', 'T', 'G', 'N', '-'};
        
        VG graph;
        
        // the node that each input sequence is extending
        unordered_map<string, Node*> current_node;
        
        // the path we're building for each aligned sequence
        unordered_map<string, Path*> aln_path;
        
        // start all of the alignments on a dummy node
        Node* dummy_node = graph.create_node("N");
        for (const auto& aln : alignments) {
            current_node[aln.first] = dummy_node;
            
            if (keep_paths) {
                Path* path = graph.graph.add_path();
                aln_path[aln.first] = path;
                path->set_name(aln.first);
            }
        }
        
        // nodes that we don't want to extend any more
        // (we never want to extend the dummy node)
        unordered_set<Node*> completed_nodes{dummy_node};
        
        size_t aln_len = alignments.begin()->second.size();
        for (size_t i = 0; i < aln_len; i++) {
#ifdef debug_msa_converter
            cerr << "## beginning column " << i << endl;
#endif
            unordered_map<Node*, char> forward_transitions;
            unordered_map<char, pair<unordered_set<Node*>, vector<string>>> transitions;
            for (const auto& aln : alignments) {
                char aln_char = toupper(aln.second[i]);
                
                if (!alphabet.count(aln_char)) {
                    cerr << "error:[MSAConverter] MSA contains non-nucleotide characters" << endl;
                    exit(1);
                    
                }
                
                Node* node_here = current_node[aln.first];
                
                if (aln_char != '-') {
                    // this alignment is transitioning to a new aligned character
                    transitions[aln_char].first.insert(node_here);
                    transitions[aln_char].second.push_back(aln.first);
                    
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
                        }
                    }
                }
            }
            
#ifdef debug_msa_converter
            cerr << "graph: " << pb2json(graph.graph) << endl;
            cerr << "node locations of sequences:" << endl;
            for (const auto& curr : current_node) {
                cerr << "\t" << curr.first << " " << curr.second->id() << endl;
            }
#endif
        }
        
        graph.destroy_node(dummy_node);
        
        return graph;
    }

}


