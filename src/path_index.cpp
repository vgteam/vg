#include "path_index.hpp"

namespace vg {

PathIndex::PathIndex(const Path& path) {
    // Just trace the path, which we assume has mapping lengths filled in.
    
    // What base are we at in the path?
    size_t path_base = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t last_rank = -1;
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        // For every mapping
        auto& mapping = path.mapping(i);
    
        if (!by_id.count(mapping.position().node_id())) {
            // This is the first time we have visited this node in the path.
            
            // Add in a mapping.
            by_id[mapping.position().node_id()] = 
                std::make_pair(path_base, mapping.position().is_reverse());
        }
        
        // Say that this node appears here along the reference in this
        // orientation.
        by_start[path_base] = NodeSide(mapping.position().node_id(), mapping.position().is_reverse());
        
        // Just advance and don't grab sequence.
        path_base += mapping_from_length(mapping);
    }
    
    // Record the length of the last mapping, since there's no next mapping to work it out from
    last_node_length = path.mapping_size() > 0 ? mapping_from_length(path.mapping(path.mapping_size() - 1)) : 0;

#ifdef debug    
    // Announce progress.
    #pragma omp critical (cerr)
    std::cerr << "Traced " << path_base << " bp path of " << path.mapping_size() << " mappings." << std::endl;
#endif
    
}

PathIndex::PathIndex(const Path& path, VG& vg) {
    // Trace the given path in the given VG graph, collecting sequence
    
    // We're going to build the sequence string
    std::stringstream seq_stream;
    
    // What base are we at in the path?
    size_t path_base = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t last_rank = -1;
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        auto& mapping = path.mapping(i);
    
        if (!by_id.count(mapping.position().node_id())) {
            // This is the first time we have visited this node in the path.
            
            // Add in a mapping.
            by_id[mapping.position().node_id()] = 
                std::make_pair(path_base, mapping.position().is_reverse());
#ifdef debug
            #pragma omp critical (cerr)
            std::cerr << "Node " << mapping.position().node_id() << " rank " << mapping.rank()
                << " starts at base " << path_base << " with "
                << vg.get_node(mapping.position().node_id())->sequence() << std::endl;
#endif
            
            // Make sure ranks are monotonically increasing along the path.
            assert(mapping.rank() > last_rank);
            last_rank = mapping.rank();
        }
        
        // Say that this node appears here along the reference in this
        // orientation.
        by_start[path_base] = NodeSide(mapping.position().node_id(), mapping.position().is_reverse());
    
    
        // Find the node's sequence
        std::string node_sequence = vg.get_node(mapping.position().node_id())->sequence();
    
        while(path_base == 0 && node_sequence.size() > 0 &&
            (node_sequence[0] != 'A' && node_sequence[0] != 'T' && node_sequence[0] != 'C' &&
            node_sequence[0] != 'G' && node_sequence[0] != 'N')) {
            
            // If the path leads with invalid characters (like "X"), throw them
            // out when computing path positions.
            
            // TODO: this is a hack to deal with the debruijn-brca1-k63 graph,
            // which leads with an X.
            #pragma omp critical (cerr)
            std::cerr << "Warning: dropping invalid leading character "
                << node_sequence[0] << " from node " << mapping.position().node_id()
                << std::endl;
                
            node_sequence.erase(node_sequence.begin());
        }
        
        if (mapping.position().is_reverse()) {
            // Put the reverse sequence in the path
            seq_stream << reverse_complement(node_sequence);
        } else {
            // Put the forward sequence in the path
            seq_stream << node_sequence;
        }
        
        // Whether we found the right place for this node in the reference or
        // not, we still need to advance along the reference path. We assume the
        // whole node (except any leading bogus characters) is included in the
        // path (since it sort of has to be, syntactically, unless it's the
        // first or last node).
        path_base += node_sequence.size();
        
        // TODO: handle leading bogus characters in calls on the first node.
    }
    
    // Record the length of the last mapping's node, since there's no next mapping to work it out from
    last_node_length = path.mapping_size() > 0 ? 
        vg.get_node(path.mapping(path.mapping_size() - 1).position().node_id())->sequence().size() :
        0;
    
    // Create the actual reference sequence we will use
    sequence = seq_stream.str();
    
#ifdef debug
    // Announce progress.
    #pragma omp critical (cerr)
    std::cerr << "Traced " << path_base << " bp path." << std::endl;
    
    if (sequence.size() < 100) {
        #pragma omp critical (cerr)
        std::cerr << "Sequence: " << sequence << std::endl;
    }
#endif
    
}

PathIndex::PathIndex(const Path& path, const xg::XG& index) {
    // Trace the given path in the given XG graph, collecting sequence
    
    // We're going to build the sequence string
    std::stringstream seq_stream;
    
    // What base are we at in the path?
    size_t path_base = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t last_rank = -1;
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        auto& mapping = path.mapping(i);
    
        if (!by_id.count(mapping.position().node_id())) {
            // This is the first time we have visited this node in the path.
            
            // Add in a mapping.
            by_id[mapping.position().node_id()] = 
                std::make_pair(path_base, mapping.position().is_reverse());
#ifdef debug
            #pragma omp critical (cerr)
            std::cerr << "Node " << mapping.position().node_id() << " rank " << mapping.rank()
                << " starts at base " << path_base << " with "
                << index.node_sequence(mapping.position().node_id()) << std::endl;
#endif
            
            // Make sure ranks are monotonically increasing along the path.
            assert(mapping.rank() > last_rank);
            last_rank = mapping.rank();
        }
        
        // Say that this node appears here along the reference in this
        // orientation.
        by_start[path_base] = NodeSide(mapping.position().node_id(), mapping.position().is_reverse());
    
    
        // Find the node's sequence
        std::string node_sequence = index.node_sequence(mapping.position().node_id());
    
        while(path_base == 0 && node_sequence.size() > 0 &&
            (node_sequence[0] != 'A' && node_sequence[0] != 'T' && node_sequence[0] != 'C' &&
            node_sequence[0] != 'G' && node_sequence[0] != 'N')) {
            
            // If the path leads with invalid characters (like "X"), throw them
            // out when computing path positions.
            
            // TODO: this is a hack to deal with the debruijn-brca1-k63 graph,
            // which leads with an X.
            #pragma omp critical (cerr)
            std::cerr << "Warning: dropping invalid leading character "
                << node_sequence[0] << " from node " << mapping.position().node_id()
                << std::endl;
                
            node_sequence.erase(node_sequence.begin());
        }
        
        if (mapping.position().is_reverse()) {
            // Put the reverse sequence in the path
            seq_stream << reverse_complement(node_sequence);
        } else {
            // Put the forward sequence in the path
            seq_stream << node_sequence;
        }
        
        // Whether we found the right place for this node in the reference or
        // not, we still need to advance along the reference path. We assume the
        // whole node (except any leading bogus characters) is included in the
        // path (since it sort of has to be, syntactically, unless it's the
        // first or last node).
        path_base += node_sequence.size();
        
        // TODO: handle leading bogus characters in calls on the first node.
    }
    
    // Record the length of the last mapping's node, since there's no next mapping to work it out from
    last_node_length = path.mapping_size() > 0 ?
        index.node_length(path.mapping(path.mapping_size() - 1).position().node_id()) :
        0;
    
    // Create the actual reference sequence we will use
    sequence = seq_stream.str();
    
#ifdef debug
    // Announce progress.
    #pragma omp critical (cerr)
    std::cerr << "Traced " << path_base << " bp path." << std::endl;
    
    if (sequence.size() < 100) {
        #pragma omp critical (cerr)
        std::cerr << "Sequence: " << sequence << std::endl;
    }
#endif

}

PathIndex::PathIndex(VG& vg, const string& path_name, bool extract_sequence) {
    // Make sure the path is present
    assert(vg.paths.has_path(path_name));
    
    if (extract_sequence) {
        // Constructor dispatch hack
        *this = PathIndex(vg.paths.path(path_name), vg);
    } else {
        *this = PathIndex(vg.paths.path(path_name));
    }
}

PathIndex::PathIndex(const xg::XG& index, const string& path_name, bool extract_sequence) {
    // Make sure the path is present
    assert(index.path_rank(path_name) != 0);
    
    if (extract_sequence) {
        // Constructor dispatch hack
        *this = PathIndex(index.path(path_name), index);
    } else {
        *this = PathIndex(index.path(path_name));
    }
}

NodeSide PathIndex::at_position(size_t position) const {
    return find_position(position)->second;
}

PathIndex::iterator PathIndex::begin() const {
    return by_start.begin();
}

PathIndex::iterator PathIndex::end() const {
    return by_start.end();
}

PathIndex::iterator PathIndex::find_position(size_t position) const {
    assert(!by_start.empty());
    
    // Look up the iterator to whatever starts after here
    auto starts_next = by_start.upper_bound(position);
    
    // This can't work if we try to look before the first node.
    assert(starts_next != by_start.begin());
    
    // Walk one back, to the node that has to own the position we asked about.
    starts_next--;
    
#ifdef debug
    cerr << "At " << position << " we have " << starts_next->second << endl;
#endif
    
    // Return that
    return starts_next;
}

size_t PathIndex::node_length(const iterator& here) const {
    assert(here != by_start.end());
    
    // Look at the next node visit
    iterator next = here;
    next++;
    
    if (next == by_start.end()) {
        // We want the length of the last visit
        return last_node_length;        
    } else {
        // We have a next node visit, so from our start to its start is our
        // length
        return next->first - here->first;
    }
}

}


















