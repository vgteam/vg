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
        
        // Remember that occurrence by node ID.
        node_occurrences[mapping.position().node_id()].push_back(by_start.find(path_base));
        
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

PathIndex::PathIndex(const list<mapping_t>& mappings, VG& vg) {
    // Trace the given path in the given VG graph, collecting sequence
    
    // We're going to build the sequence string
    std::stringstream seq_stream;
    
    // What base are we at in the path?
    size_t path_base = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t last_rank = -1;
    
    for (auto& mapping : mappings) {
    
        if (!by_id.count(mapping.node_id())) {
            // This is the first time we have visited this node in the path.
            
            // Add in a mapping.
            by_id[mapping.node_id()] = 
                std::make_pair(path_base, mapping.is_reverse());
#ifdef debug
            #pragma omp critical (cerr)
            std::cerr << "Node " << mapping.node_id() << " rank " << mapping.rank()
                << " starts at base " << path_base << " with "
                << vg.get_node(mapping.node_id())->sequence() << std::endl;
#endif
            
            // Make sure ranks are monotonically increasing along the path, or
            // unset.
            assert(mapping.rank > last_rank || (mapping.rank == 0 && last_rank == 0));
            last_rank = mapping.rank;
        }
        
        // Say that this node appears here along the reference in this
        // orientation.
        by_start[path_base] = NodeSide(mapping.node_id(), mapping.is_reverse());
    
        // Remember that occurrence by node ID.
        node_occurrences[mapping.node_id()].push_back(by_start.find(path_base));
    
        // Say this Mapping happens at this base along the path
        mapping_positions[&mapping] = path_base;
    
        // Find the node's sequence
        std::string node_sequence = vg.get_node(mapping.node_id())->sequence();
    
        while(path_base == 0 && node_sequence.size() > 0 &&
            (node_sequence[0] != 'A' && node_sequence[0] != 'T' && node_sequence[0] != 'C' &&
            node_sequence[0] != 'G' && node_sequence[0] != 'N')) {
            
            // If the path leads with invalid characters (like "X"), throw them
            // out when computing path positions.
            
            // TODO: this is a hack to deal with the debruijn-brca1-k63 graph,
            // which leads with an X.
            #pragma omp critical (cerr)
            std::cerr << "Warning: dropping invalid leading character "
                << node_sequence[0] << " from node " << mapping.node_id()
                << std::endl;
                
            node_sequence.erase(node_sequence.begin());
        }
        
        if (mapping.is_reverse()) {
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
        vg.get_node(mappings.back().node_id())->sequence().size();
    
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

PathIndex::PathIndex(const Path& path, const HandleGraph& graph) {
    // Trace the given path in the given graph, collecting sequence
    
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
                << graph.get_sequence(graph.get_handle(mapping.position().node_id())) << std::endl;
#endif
            
            // Make sure ranks are monotonically increasing along the path, or
            // unset.
            assert(mapping.rank() > last_rank || (mapping.rank() == 0 && last_rank == 0));
            last_rank = mapping.rank();
        }
        
        // Say that this node appears here along the reference in this
        // orientation.
        by_start[path_base] = NodeSide(mapping.position().node_id(), mapping.position().is_reverse());
    
        // Remember that occurrence by node ID.
        node_occurrences[mapping.position().node_id()].push_back(by_start.find(path_base));
    
        // Find the node's sequence
        std::string node_sequence = graph.get_sequence(graph.get_handle(mapping.position().node_id()));
    
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
        graph.get_length(graph.get_handle(path.mapping(path.mapping_size() - 1).position().node_id())) :
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

PathIndex::PathIndex(const PathHandleGraph& graph, const string& path_name, bool extract_sequence) {
    // Make sure the path is present
    assert(graph.has_path(path_name));
    
    // Make a Protobuf path object
    auto path = path_from_path_handle(graph, graph.get_path_handle(path_name));

    if (extract_sequence) {
        // Constructor dispatch hack
        *this = PathIndex(path, graph);
    } else {
        *this = PathIndex(path);
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
        path_base += mapping.length;
    }
}

bool PathIndex::path_contains_node(int64_t node_id) const {
    if (by_id.find(node_id) != by_id.end()){
        return true;
    }
    return false;
}

bool PathIndex::path_contains_node_in_orientation(int64_t node_id, bool is_reverse) const {
    return find_in_orientation(node_id, is_reverse) != end();
}

PathIndex::iterator PathIndex::find_in_orientation(int64_t node_id, bool is_reverse) const {
    auto found_occurrences = node_occurrences.find(node_id);
    
    if (found_occurrences == node_occurrences.end()) {
        // There are no occurrences
        return end();
    }
    
    for (auto& occ : found_occurrences->second) {
        // Do a linear scan for the correct orientation.
        // TODO: index by orientation
        if (occ->second.is_end == is_reverse) {
            // We found an occurrence in the requested orientation.
            return occ;
        }
    }
    
    return end();
}

pair<bool, bool> PathIndex::get_contained_orientations(int64_t node_id) const {
    // TODO: Do scans manually to be twice as fast!
    return make_pair(path_contains_node_in_orientation(node_id, false), path_contains_node_in_orientation(node_id, true));
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

    // Make sure we didn't fall off the ends
    assert(position - starts_next->first < node_length(starts_next));
    
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

pair<size_t, size_t> PathIndex::round_outward(size_t start, size_t past_end) const {
    // Find the node occurrence the start position is on
    auto start_occurrence = find_position(start);
    // Seek to the start of that occurrence
    size_t start_rounded = start_occurrence->first;
    
    // Now try and round the end
    size_t past_end_rounded;
    if (past_end == 0) {
        // Range must have been empty anyway, so keep it ending before the first
        // node.
        past_end_rounded = 0;
    } else {
        // Look for the node holding the last included base
        auto end_occurrence = find_position(past_end - 1);
        // Then go out past its end.
        past_end_rounded = end_occurrence->first + node_length(end_occurrence);
    }
    
    return make_pair(start_rounded, past_end_rounded);
}

map<id_t, vector<Mapping>> PathIndex::parse_translation(const Translation& translation) {

    // We take as a precondition that the translation is replacing a set of old
    // nodes each with a nonempty set of new nodes. So we won't have to combine
    // nodes or parts of nodes.
    
#ifdef debug
    cerr << "Partitioning translation: " << pb2json(translation) << endl;
#endif
    
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
    
    for (auto kv : old_node_to_new_nodes) {
        // For every node that needs to be replaced with other nodes
        auto& old_id = kv.first;
        auto& replacements = kv.second;
        
        // Pull out all the occurrences of the old node.
        vector<iterator> occurrences{std::move(node_occurrences[old_id])};
        node_occurrences.erase(old_id);
        
        for (auto& occurrence : occurrences) {
            // Each time the old node appeared, replace it (and log occurrences
            // of the new nodes that partition it, one of which may re-use the
            // ID)
            replace_occurrence(occurrence, replacements);
        }
        
#ifdef debug
        cerr << "by_start is now: " << endl;
        for (auto kv2 : by_start) {
            cerr << "\t" << kv2.first << ": " << kv2.second << endl;
        }
#endif
    }
}

void PathIndex::apply_translations(const vector<Translation>& translations) {
    // Convert from normal to partitioning translations
    
    // For each original node ID, we keep a vector of pairs of from mapping and
    // to mapping. We only keep pairs where the from mapping isn't empty.
    map<id_t, vector<pair<Mapping, Mapping>>> collated;
    
    for (auto& t : translations) {
        if (t.from().mapping_size() < 1 || t.to().mapping_size() != 1) {
            // Ensure the translations are the format we expect. They always have
            // at least one from mapping (but maybe an insert too) and exactly 1
            // to mapping.
            cerr << "error:[vg::PathIndex] Bad translation: " << pb2json(t) << endl;
            throw runtime_error("Translation not in VG::edit() format");
        }
        
        if (mapping_from_length(t.from().mapping(0)) == 0) {
            // This is a novel node and can't be on our path
            continue;
        }
        
        if (t.from().mapping(0).position().is_reverse()) {
            // Wait for the forward-orientation version
            continue;
        }
        
        // Stick the from and to mappings in the list for the from node
        collated[t.from().mapping(0).position().node_id()].push_back(make_pair(t.from().mapping(0), t.to().mapping(0)));
    }
    
    for (auto& kv : collated) {
        // For every original node and its replacement nodes
        
        // Sort the replacement mappings
        std::sort(kv.second.begin(), kv.second.end(), [](const pair<Mapping, Mapping>& a, const pair<Mapping, Mapping>& b) {
            // Return true if the a pair belongs before the b pair along the path through the original node
            return a.first.position().offset() <= b.first.position().offset();
        });
        
        // Make a new translation to cover the original node
        Translation covering;
        
        for (auto mapping_pair : kv.second) {
            // Split across these parts of new nodes
            *(covering.mutable_to()->add_mapping()) = mapping_pair.second;
        }
        
        // Just assume we take up the whole original node
        auto* from_mapping = covering.mutable_from()->add_mapping();
        from_mapping->mutable_position()->set_node_id(kv.first);
        // Give it a full length perfect match
        auto* from_edit = from_mapping->add_edit();
        from_edit->set_from_length(path_from_length(covering.to()));
        from_edit->set_to_length(from_edit->from_length());
        
        // Apply this (single node) translation.
        // TODO: batch up a bit?
        apply_translation(covering);
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
        
        // Remember that occurrence by node ID.
        node_occurrences[new_id].push_back(by_start.find(start));
        
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


















