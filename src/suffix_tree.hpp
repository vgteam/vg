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
    
}


#endif /* suffix_tree_hpp */
