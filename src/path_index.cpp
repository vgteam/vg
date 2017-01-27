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

PathIndex::PathIndex(const list<Mapping>& mappings, VG& vg) {
    // Trace the given path in the given VG graph, collecting sequence
    
    // We're going to build the sequence string
    std::stringstream seq_stream;
    
    // What base are we at in the path?
    size_t path_base = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t last_rank = -1;
    
    for (auto& mapping : mappings) {
    
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
    
        // Say this Mapping happens at this base along the path
        mapping_positions[&mapping] = path_base;
    
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
    last_node_length = mappings.empty() ?
        0 : 
        vg.get_node(mappings.back().position().node_id())->sequence().size();
    
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

    // Follow the path (again) and place all its Mappings
    
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
        *this = PathIndex(vg.paths.get_path(path_name), vg);
    } else {
        *this = PathIndex(vg.paths.path(path_name));
        update_mapping_positions(vg, path_name);
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

void PathIndex::update_mapping_positions(VG& vg, const string& path_name) {
    // Brute force recalculate mapping positions.
    // Ignores bad characters...
    
    // TODO: Don't make this brute force. Integrate into Paths or keep our own
    // original-rank-based index or something.
    
    mapping_positions.clear();
    
    // Where are we in the path?
    size_t path_base = 0;
    
    for (auto& mapping : vg.paths.get_path(path_name)) {
        // For each mapping currently in the path, remember its start position
        // along the path.
        mapping_positions[&mapping] = path_base;
        
        // Go right by its length
        path_base += mapping_to_length(mapping);
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

map<id_t, vector<Mapping>> PathIndex::parse_translation(const Translation& translation) {

    // We take as a precondition that the translation is replacing a set of old
    // nodes each with a nonempty set of new nodes. So we won't have to combine
    // nodes or parts of nodes.
    
    // We'll populate this with the mappings that partition each old node.
    map<id_t, vector<Mapping>> old_node_to_new_nodes;
    
    // We know the new Mappings are conceptually nested in the old Mappings, so
    // we can use nested loops.

    // How many bases in the old and new paths are accounted for?
    size_t old_bases = 0;
    size_t new_bases = 0;

    // This represents our index in the new path
    size_t j = 0;    
    
    for(size_t i = 0; i < translation.from().mapping_size(); i++) {
        // For every old mapping
        auto& from_mapping = translation.from().mapping(i);
        
        // Count up its bases
        old_bases += mapping_from_length(from_mapping);
        
        // Grab a reference to the list of replacement mappings
        auto& replacements = old_node_to_new_nodes[from_mapping.position().node_id()];
        
        // We know the old mapping must have at least one new mapping in it
        do {
            // For each mapping in the new path, copy it
            auto to_mapping = translation.to().mapping(j);
            
            if (from_mapping.position().is_reverse()) {
                // Flip its strand if the mapping we're partitioning is backward
                to_mapping.mutable_position()->set_is_reverse(!to_mapping.position().is_reverse());
            }
            
            // Account for its bases
            new_bases += mapping_from_length(to_mapping);
            
            // Copy it into the list for just this from node
            replacements.push_back(to_mapping);
            
            // Look at the next to mapping
            j++;
        } while (j < translation.to().mapping_size() && new_bases < old_bases);
        
        if (from_mapping.position().is_reverse()) {
            // Flip the order of the replacement mappings around
            reverse(replacements.begin(), replacements.end());
        }
        
#ifdef debug
        cerr << "Old node " << from_mapping.position().node_id() << " "
            << from_mapping.position().is_reverse() << " becomes: " << endl;
        for(auto& m : old_node_to_new_nodes[from_mapping.position().node_id()]) {
            cerr << "\t" << pb2json(m) << endl;
        }
#endif
    }
    
    return old_node_to_new_nodes;

}

void PathIndex::apply_translation(const Translation& translation) {
    
    // Parse the translation, to get a map form old node ID to vector of
    // replacement mappings.
    auto old_node_to_new_nodes = parse_translation(translation);
    
    
    // TODO: we would like to update mapping_positions efficiently, but we
    // can't, because it's full of potentially invalidated pointers.
    mapping_positions.clear();
    
    // We really need to update by_start (and maybe last_node_length)
    
    // TODO: We need another index from node ID to all its start positions to pull it off efficiently.
    
    // For now just do a dumb loop.
    
    for(iterator here = by_start.begin(); here != by_start.end(); ++here) {
        // For every by_start entry
        
        // Grab the node ID
        auto node_id = here->second.node;
        
        if (!old_node_to_new_nodes.count(node_id)) {
            // No replacement to be done. Try the next node occurrence.
            continue;
        }
        
        // Grab the Mappings describing the replacements
        auto replacements = old_node_to_new_nodes.at(node_id);
        
        // Replace this occurrence of this mapping with these partitioning
        // replacements.
        replace_occurrence(here, replacements);
        
    }
}

void PathIndex::replace_occurrence(iterator to_replace, const vector<Mapping>& replacements) {
    // Grab the node ID
    auto node_id = to_replace->second.node;
    
    // Grab it's start
    auto start = to_replace->first;
    
    // Determine if we want to insert replacement nodes forward or backward
    bool reverse = to_replace->second.is_end;
    
    // Peek ahead at whatever's next (so we know if we're replacing the last
    // node in the path).
    auto comes_next = to_replace;
    comes_next++;
    
    if (by_id.count(node_id) && by_id.at(node_id).first == start) {
        // We're removing the first occurrence of this node, so we need to
        // clean up by_id. But since we're removing all occurreences of this
        // node we don't have to point it at the second occurrence.
        by_id.erase(node_id);
    }
    
    for (int64_t i = reverse ? (replacements.size() - 1) : 0;
        i != (reverse ? -1 : replacements.size());
        i += (reverse ? -1 : 1)) {
        
        // For each replacement mapping in the appropriate order
        auto& mapping = replacements.at(i);
        
        // What ID do we put?
        auto new_id = mapping.position().node_id();
        
        // What orientation doies it go in?
        auto new_orientation = mapping.position().is_reverse() != reverse;
        
        // Stick the replacements in the map
        by_start[start] = NodeSide(new_id, new_orientation);
        
        if (!by_id.count(new_id) || by_id.at(new_id).first > start) {
            // We've created a new first mapping to this new node.
            // Record it.
            by_id[new_id] = make_pair(start, new_orientation);
        }
        
        // Budge start up so the next mapping gets inserted after this one.
        start += mapping_from_length(mapping);
        
        if (comes_next == by_start.end() && i == (reverse ? 0 : replacements.size() - 1)) {
            // We just added the last mapping replacing what the old last
            // mapping was. So update the length of the last node to reflect
            // this new last node.
            last_node_length = mapping_from_length(mapping);
        }
    }
}

}


















