//
//  suffix_tree.hpp
//  
// Suffix tree implementation including Ukkonen's linear time construction algorithm
//


#ifndef suffix_tree_hpp
#define suffix_tree_hpp

#include <stdio.h>
#include <unordered_map>
#include <list>
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;

namespace vg {
    
    /**
     * An implementation of a suffix tree with linear time and space complexity for construction.
     *
     */
    class SuffixTree {
    
    public:
        /// Linear time constructor.
        ///
        /// Note: string cannot have null characters, but this is not checked.
        SuffixTree(string::const_iterator begin, string::const_iterator end);
        ~SuffixTree() = default;
        
        /// Returns the length of the longest prefix of str that exactly matches
        /// a suffix of the string used to construct the suffix tree.
        ///
        /// Note: string cannot have null characters, but this is not checked.
        size_t longest_overlap(const string& str);
        size_t longest_overlap(string::const_iterator begin, string::const_iterator end);
        
        /// Retuns a vector of all of the indices where a string occurs as a substring
        /// of the string used to construct the suffix tree. Indices are ordered arbitrarily.
        ///
        /// Note: string cannot have null characters, but this is not checked.
        vector<size_t> substring_locations(const string& str);
        vector<size_t> substring_locations(string::const_iterator begin, string::const_iterator end);
        
        /// Beginning of string used to make tree
        const string::const_iterator begin;
        /// End of string used to make tree
        const string::const_iterator end;
        
    private:
        struct STNode;
        
        /// All nodes in the tree (in a list to avoid difficulties with pointers and reallocations)
        list<STNode> nodes;
        
        /// The edges from the root node
        unordered_map<char, STNode*> root;
        
        /// Returns a char of the string or null char at past-the-last index
        inline char get_char(size_t i);
        
        // debugging functions for constructor
        string to_string();
        string partial_tree_to_string(int64_t phase);
        string label_string(unordered_map<char, STNode*>* branch_point);
        string node_string(STNode* node, int64_t phase);
        string active_point_string(unordered_map<char, STNode*>* active_branch_point,
                                   STNode* active_node,
                                   int64_t active_length,
                                   int64_t phase);
        string suffix_links_string(unordered_map<unordered_map<char, STNode*>*,
                                                 unordered_map<char, STNode*>*>& suffix_links);
    };
    
    /**
     * A node of a suffix tree corresponding a substring of the string
     *
     */
    struct SuffixTree::STNode {
        /// Constructor
        STNode(int64_t first, int64_t last);
        ~STNode() = default;
        
        /// Edges down the tree
        unordered_map<char, STNode*> children;
        
        /// First index of string on this node
        int64_t first;
        /// Last index of string on this node, inclusive (-1 indicates end sentinel during consruction)
        int64_t last;
        
        /// The length of the the node during a phase of construction
        inline int64_t length(int64_t phase) {
            return last >= 0 ? last - first + 1 : phase - first + 1;
        }
        
        /// The last index contained on the this node during a phase of construction
        inline int64_t final_index(int64_t phase) {
            return last >= 0 ? last - first : phase - first;
        }
    };
}


#endif /* suffix_tree_hpp */
