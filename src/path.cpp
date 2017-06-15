#include "path.hpp"
#include "stream.hpp"


namespace vg {

Paths::Paths(void) {
    // noop
}

void Paths::load(istream& in) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    function<void(Path&)> lambda = [this](Path& p) {
        this->extend(p);
    };
    stream::for_each(in, lambda, handle_count);
}

void Paths::write(ostream& out) {
    vector<string> path_names;
    for (auto& p : _paths) {
        const string& name = p.first;
        path_names.push_back(name);
    }
    function<Path(uint64_t)> lambda =
        [this, &path_names](uint64_t i) -> Path {
        list<Mapping>& mappings = _paths[path_names.at(i)];
        Path path;
        for (auto& m : mappings) {
            Mapping* nm = path.add_mapping();
            *nm = m;
        }
        path.set_name(path_names.at(i));
        if (circular.count(path_names.at(i))) {
            path.set_is_circular(true);
        }
        return path;
    };
    stream::write(out, _paths.size(), lambda);
}

void Paths::to_graph(Graph& g) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<Mapping>& mappings = p.second;
        Path* path = g.add_path();
        path->set_name(name);
        if (circular.count(name)) {
            path->set_is_circular(true);
        }
        for (auto& m : mappings) {
            Mapping* nm = path->add_mapping();
            *nm = m;
        }
    }
}

Path Paths::path(const string& name) {
    Path path;
    auto p = _paths.find(name);
    if (p == _paths.end()) {
        return path;
    }
    list<Mapping>& mappings = p->second;
    path.set_name(name);
    for (auto& m : mappings) {
        Mapping* nm = path.add_mapping();
        *nm = m;
    }
    if (circular.count(name)) {
        path.set_is_circular(true);
    }
    return path;
}

void Paths::for_each(const function<void(const Path&)>& lambda) {
    for (auto& p : _paths) {
        const string& name = p.first;
        lambda(path(name));
    }
}

void Paths::for_each_name(const function<void(const string&)>& lambda) {
    for (auto& p : _paths) {
        const string& name = p.first;
        lambda(name);
    }
}

void Paths::for_each_mapping(const function<void(Mapping*)>& lambda) {
    for (auto& p : _paths) {
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            lambda(&m);
        }
    }
}

void Paths::for_each_stream(istream& in, const function<void(Path&)>& lambda) {
    uint64_t count = 0;
    function<void(uint64_t)> handle_count = [this, &count](uint64_t c) { count = c; };
    stream::for_each(in, lambda, handle_count);
}

void Paths::make_circular(const string& name) {
    circular.insert(name);
}

void Paths::make_linear(const string& name) {
    circular.erase(name);
}

void Paths::extend(const Path& p) {
    const string& name = p.name();
    list<Mapping>& path = get_create_path(name);
    for (int i = 0; i < p.mapping_size(); ++i) {
        const Mapping& m = p.mapping(i);
        append_mapping(name, m);
    }
    if (p.is_circular()) {
        make_circular(name);
    }
    // re-sort?
    sort_by_mapping_rank();
    rebuild_mapping_aux();
}

// one of these should go away
void Paths::extend(Paths& p) {
    for (auto& l : p._paths) {
        const string& name = l.first;
        list<Mapping>& path = l.second;
        // Make sure we preserve empty paths
        get_create_path(name);
        for (auto& m : path) {
            append_mapping(name, m);
        }
        if (p.circular.count(name)) {
            make_circular(name);
        }
    }
    sort_by_mapping_rank();
    rebuild_mapping_aux();
}

void Paths::append(Paths& paths) {
    for (auto& p : paths._paths) {
        const string& name = p.first;
        const list<Mapping>& path = p.second;
        // Make sure we preserve empty paths
        get_create_path(name);
        for (auto& m : path) {
            append_mapping(name, m);
        }
        if (paths.circular.count(name)) {
            make_circular(name);
        }
    }
    sort_by_mapping_rank();
    rebuild_mapping_aux();
}

void Paths::append(Graph& g) {
    for (int i = 0; i < g.path_size(); ++i) {
        const Path& p = g.path(i);
        // Make sure we preserve empty paths
        get_create_path(p.name());
        for (int j = 0; j < p.mapping_size(); ++j) {
            const Mapping& m = p.mapping(j);
            append_mapping(p.name(), m);
            if (p.is_circular()) {
                make_circular(p.name());
            }
        }
    }
}

Path& append_path(Path& a, const Path& b) {
    a.mutable_mapping()->MergeFrom(b.mapping());
    return a;
}


bool Paths::has_mapping(const string& name, size_t rank) {
    return mappings_by_rank.count(name) && mappings_by_rank[name].count(rank);
}

void Paths::append_mapping(const string& name, const Mapping& m) {
    // get or create the path with this name
    list<Mapping>& pt = get_create_path(name);
    // now if we haven't already supplied a mapping
    // add it
    
    if (!m.rank() || !has_mapping(name, m.rank())) {
        // If we don't have a rank set or we don't have a mapping in this path
        // with that rank, we need to add the mapping.
        
        // First figure out what the previous rank on the path is.
        size_t last_rank = 0;
        if(!pt.empty()) {
            last_rank = (*pt.rbegin()).rank();
        }

        // Add this mapping at the end of the path. Note that this may not
        // actually be where it belongs according to its rank; if that is the
        // case, sort_by_mapping_rank() and rebuild_mapping_aux() need to be
        // called, after all the mappings are loaded, in order to put them in
        // order by rank.
        pt.push_back(m);
        Mapping* mp = &pt.back();
        
        if(mp->rank() == 0) {
            // 0 rank defaults to being ranked at the end of what's already
            // there. After all, that's where we just added it.
            if(&pt.front() != mp) {
                // There are other things on the path
                if(last_rank && !mappings_by_rank[name].count(last_rank + 1)) {
                    // We can go right after the thing that we're physically
                    // after.
                    mp->set_rank(last_rank + 1);
                } else {
                    // We're putting it after something which itself has no
                    // rank, or the next rank is already occupied. Leave rank
                    // unset; the caller will have to globally rebuild the rank
                    // info later with rebuild_mapping_aux()
                    
                    // Really this is undefined ordering according to the
                    // precondition, but this makes sense.
                }
            } else {
                // This is the first mapping, so give it minimal rank. Ordering
                // is undefined with respect to any pre-ranked mappings that may
                // be added to the path later.
                mp->set_rank(1);
            }
        }
        
        // add it to the node mappings
        auto& ms = get_node_mapping(mp->position().node_id());
        ms[name].insert(mp);
        // and record its position in this list
        list<Mapping>::iterator mi = pt.end(); --mi;
        mapping_itr[mp] = mi;
        mapping_path[mp] = name;
        if(mp->rank()) {
            // Only if we actually end up with a rank (i.e. all the existing
            // ranks weren't cleared) do we really index by rank.
            mappings_by_rank[name][mp->rank()] = mp;
        }
    }
}

void Paths::append_mapping(const string& name, id_t id, size_t rank, bool is_reverse) {
    Mapping m;
    m.mutable_position()->set_node_id(id);
    m.mutable_position()->set_is_reverse(is_reverse);
    // If the rank passed in is 0, it will get filled in by the other version of
    // append_mapping.
    m.set_rank(rank);
    append_mapping(name, m);
}

void Paths::prepend_mapping(const string& name, const Mapping& m) {
    // get or create the path with this name
    list<Mapping>& pt = get_create_path(name);
    
    // We can't prepend a mapping that doesn't have a rank set. We would like to
    // generate ranks, but we can't keep decrementing the first rank
    // indefinitely, and that might not be correct. Also, what rank would we use
    // for the only mapping in a path?
    assert(m.rank());
    
    // now if we haven't already supplied a mapping
    // add it
    if (!has_mapping(name, m.rank())) {
        // If we don't have a rank set or we don't have a mapping in this path
        // with that rank, we need to add the mapping.
        
        pt.push_front(m);
        Mapping* mp = &pt.front();
        // add it to the node mappings
        auto& ms = get_node_mapping(m.position().node_id());
        ms[name].insert(mp);
        // and record its position in this list
        list<Mapping>::iterator mi = pt.begin();
        mapping_itr[mp] = mi;
        mapping_path[mp] = name;
        mappings_by_rank[name][mp->rank()] = mp;
    }
}

void Paths::prepend_mapping(const string& name, id_t id, size_t rank, bool is_reverse) {
    Mapping m;
    m.mutable_position()->set_node_id(id);
    m.mutable_position()->set_is_reverse(is_reverse);
    if (rank) {
        m.set_rank(rank);
    } else {
        m.set_rank(get_path(name).size()+1); // rank is 1-based
    }
    prepend_mapping(name, m);
}

size_t Paths::get_next_rank(const string& name) {
    auto& p = get_path(name);
    //cerr << "next rank be " << p.size()+1 << " or " << (size_t) p.rend()->rank()+1 << endl;
    return max(p.size()+1, (size_t) (p.size() ? p.rend()->rank()+1 : 0));
}

// these will split a mapping into two
// NB: each submapping ends up with the same rank as the parent
// however, they will be ordered correctly in the path
// we will need to normalize path ranks to make this right
pair<Mapping*, Mapping*> Paths::divide_mapping(Mapping* m, const Position& pos) {
    // this is needed to split mappinsg during e.g. normalization
    // but still ensure that the mappings are out there
    // what do we do?
    // first we take the mapping and divide it as we do
    auto n = cut_mapping(*m, pos);
    return replace_mapping(m, n);
}

pair<Mapping*, Mapping*> Paths::divide_mapping(Mapping* m, size_t offset) {
    auto n = cut_mapping(*m, offset);
    return replace_mapping(m, n);
}

pair<Mapping*, Mapping*> Paths::replace_mapping(Mapping* m, pair<Mapping, Mapping> n) {
    // then we remove it from the node it's pointing to
    // and replace it with the other two mappings
    // we'll give them the same rank, but record them in the right order
    // this leaves an invalid graph
    // there are a few ways to fix this--- they involve changing the way we record ranks
    // but for now it's going to be simplest if the calling context manages this
    auto& path_name = mapping_path_name(m);
    n.first.set_rank(m->rank());
    n.second.set_rank(m->rank());
    // remove the mapping, getting an iterator pointing to the element that was after it
    if (!m->position().is_reverse()) {
        auto i = remove_mapping(m);
        // the insertion will happen in reverse order
        // because insert puts the position before the iterator it's given
        // so we first insert the second element
        auto j = insert_mapping(i, path_name, n.second);
        // and then the first
        auto k = insert_mapping(j, path_name, n.first);
        // and we return them in proper order
        return make_pair(&*j, &*k);
    } else {
        // things get flipped around for reversed mappings
        auto i = remove_mapping(m);
        auto j = insert_mapping(i, path_name, n.first);
        auto k = insert_mapping(j, path_name, n.second);
        // and we return them in proper order
        return make_pair(&*k, &*j);
    }
}

bool Paths::has_path(const string& name) {
    return _paths.find(name) != _paths.end();
}

void Paths::increment_node_ids(id_t inc) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            m.mutable_position()->set_node_id(m.position().node_id()+inc);
        }
    }
    rebuild_node_mapping();
}

void Paths::swap_node_ids(hash_map<id_t, id_t>& id_mapping) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            // Look up the replacement ID
            auto replacement = id_mapping.find(m.position().node_id());
            if(replacement != id_mapping.end()) {
                // If there is a replacement, use it.
                m.mutable_position()->set_node_id((*replacement).second);
            }
        }
    }
    rebuild_node_mapping();
}

void Paths::reassign_node(id_t new_id, Mapping* m) {
    // erase the old node id
    node_mapping[m->position().node_id()][mapping_path_name(m)].erase(m);
    // set the new node id
    m->mutable_position()->set_node_id(new_id);
    // and record it in the new node record
    node_mapping[m->position().node_id()][mapping_path_name(m)].insert(m);
}

void Paths::rebuild_node_mapping(void) {
    // starts with paths and rebuilds the index
    node_mapping.clear();
    for (auto& p : _paths) {
        const string& path_name = p.first;
        list<Mapping>& path = p.second;
        for (auto& m : path) {
            get_node_mapping(m.position().node_id())[path_name].insert(&m);
        }
    }
}

// attempt to sort the paths based on the recorded ranks of the mappings
void Paths::sort_by_mapping_rank(void) {
    for (auto p = _paths.begin(); p != _paths.end(); ++p) {
        list<Mapping>& path = p->second;
        path.sort([](const Mapping& m1, const Mapping& m2) {
                return m1.rank() < m2.rank();
            });
    }
}

// compact the ranks preserving the relative rank order
void Paths::compact_ranks(void) {
    // first ensure the storage order of the mappings is correct
    sort_by_mapping_rank();
    // clear the ranks
    clear_mapping_ranks();
    // and rebuild them and other aux data structures
    rebuild_node_mapping();
    rebuild_mapping_aux();
}

void Paths::rebuild_mapping_aux(void) {
    mapping_itr.clear();
    mapping_path.clear();
    mappings_by_rank.clear();
    for (auto& p : _paths) {
        const string& path_name = p.first;
        list<Mapping>& path = p.second;
        size_t order_in_path = 0;
        for (list<Mapping>::iterator i = path.begin(); i != path.end(); ++i) {
            mapping_itr[&*i] = i;
            mapping_path[&*i] = path_name;
            
            if(i->rank() > order_in_path + 1) {
                // Make sure that if we have to assign a rank to a node after
                // this one, it is greater than this node's rank. TODO: should
                // we just uniformly re-rank all the nodes starting at 0? Or
                // will we ever want to cut and paste things back together using
                // the old preserved ranks?
                order_in_path = i->rank() - 1;
            }
            
            if (i->rank() == 0 || i->rank() < order_in_path + 1) {
                // If we don't already have a rank, or if we see a rank that
                // can't be correct given the ranks we have already seen, we set
                // the rank based on what we've built
                i->set_rank(order_in_path+1);
            }
            
            // Save the mapping as being at the given rank in its path.
            mappings_by_rank[path_name][i->rank()] = &*i;
            
            ++order_in_path;
        }
    }
}

void Paths::remove_node(id_t id) {
    node_mapping.erase(id);
}

list<Mapping>::iterator Paths::find_mapping(Mapping* m) {
    return mapping_itr[m];
}

list<Mapping>::iterator Paths::remove_mapping(Mapping* m) {
    // The mapping has to exist
    assert(mapping_path.find(m) != mapping_path.end());
    const string& path_name = mapping_path[m];
    id_t id = m->position().node_id();
    auto& x = _paths[path_name];
    
    // This gets tricky because we're going to deallocate the storage pointed to
    // by m. We need to remove it from other things first.
    if(m->rank() && mappings_by_rank[path_name].count(m->rank()) && 
        mappings_by_rank[path_name][m->rank()] == m) {
        // If we have this node stored for its path and rank, kick it out.
        mappings_by_rank[path_name].erase(m->rank());
    }
    
    // Actually deallocate the mapping
    list<Mapping>::iterator p = _paths[path_name].erase(mapping_itr[m]);
    if (has_node_mapping(id)) {
        auto& node_path_mapping = get_node_mapping(id);
        node_path_mapping[path_name].erase(m);
        if (node_path_mapping.empty()) node_mapping.erase(id);
    }
    mapping_path.erase(m);
    mapping_itr.erase(m);
    
    return p;
}

list<Mapping>::iterator Paths::insert_mapping(list<Mapping>::iterator w, const string& path_name, const Mapping& m) {
    auto px = _paths.find(path_name);
    assert(px != _paths.end());
    list<Mapping>& path = px->second;
    list<Mapping>::iterator p;
    if (path.empty()) {
        path.push_front(m);
        p = path.begin();
    } else {
        p = path.insert(w, m);
    }
    get_node_mapping(m.position().node_id())[path_name].insert(&*p);
    mapping_path[&*p] = path_name;
    mapping_itr[&*p] = p;
    return p;
}

void Paths::to_json(ostream& out) {
    function<void(const Path&)> lambda = [this, &out](const Path& p) {
        out << pb2json(p) <<endl;
    };
    for_each(lambda);
}

size_t Paths::size(void) const {
    return _paths.size();
}

bool Paths::empty(void) const {
    return _paths.size() == 0;
}

void Paths::clear(void) {
    _paths.clear();
    node_mapping.clear();
    mapping_path.clear();
    mappings_by_rank.clear();
}

void Paths::clear_mapping_ranks(void) {
    for (auto p = _paths.begin(); p != _paths.end(); ++p) {
        list<Mapping>& path = p->second;
        for (auto m = path.begin(); m != path.end(); ++m) {
            Mapping& mapping = *m;
            mapping.set_rank(0);
        }
    }
    mappings_by_rank.clear();
}

list<Mapping>& Paths::get_path(const string& name) {
    return _paths[name];
}

void Paths::remove_paths(const set<string>& names) {
    for (auto& name : names) {
        _paths.erase(name);
    }
    rebuild_node_mapping();
    rebuild_mapping_aux();
}

void Paths::remove_path(const string& name) {
    auto& path = _paths[name];
    
    for(auto& mapping : path) {
        // Unindex all the mappings
        mapping_itr.erase(&mapping);
        mapping_path.erase(&mapping);
        if(node_mapping.count(mapping.position().node_id())) {
            // Throw out all the mappings for this path on this node
            node_mapping[mapping.position().node_id()].erase(name);
        }
    }

    // Delete the actual mappings
    _paths.erase(name);
    
    // Delete the indexes that are by path name
    mappings_by_rank.erase(name);

}

void Paths::keep_paths(const set<string>& names) {
    set<string> to_remove;
    for (auto& p : _paths) {
        if (!names.count(p.first)) {
            to_remove.insert(p.first);
        }
    }
    remove_paths(to_remove);
}

list<Mapping>& Paths::create_path(const string& name) {
    return _paths[name];
}

list<Mapping>& Paths::get_create_path(const string& name) {
    if (!has_path(name)) {
        return create_path(name);
    } else {
        return get_path(name);
    }
}

bool Paths::has_node_mapping(id_t id) {
    return node_mapping.find(id) != node_mapping.end();
}

bool Paths::has_node_mapping(Node* n) {
    return node_mapping.find(n->id()) != node_mapping.end();
}

map<string, set<Mapping*>>& Paths::get_node_mapping(id_t id) {
    return node_mapping[id];
}

map<string, set<Mapping*>>& Paths::get_node_mapping(Node* n) {
    return node_mapping[n->id()];
}

map<string, map<int, Mapping*>> Paths::get_node_mappings_by_rank(id_t id) {
    map<string, map<int, Mapping*>> by_ranks;
    for (auto& p : get_node_mapping(id)) {
        auto& name = p.first;
        auto& mp = p.second;
        for (auto* m : mp) by_ranks[name][m->rank()] = m;
    }
    return by_ranks;
}

map<string, map<int, Mapping>> Paths::get_node_mapping_copies_by_rank(id_t id) {
    map<string, map<int, Mapping>> by_ranks;
    for (auto& p : get_node_mapping(id)) {
        auto& name = p.first;
        auto& mp = p.second;
        for (auto* m : mp) by_ranks[name][m->rank()] = *m;
    }
    return by_ranks;
}

Mapping* Paths::traverse_left(Mapping* mapping) {
    // Get the iterator for this Mapping*
    list<Mapping>::iterator place = mapping_itr.at(mapping);

    // Get the path name for this Mapping*
    string path_name = mapping_path_name(mapping);

    // Get the list that the iterator is in
    list<Mapping>& path_list = _paths.at(path_name);

    // If we're already the beginning, return null.
    if(place == path_list.begin()) {
        return nullptr;
    }

    // Else walk left and return the address of the stored Mapping. std::list
    // iterators are bidirectional, so we will be able to do it.
    place--;
    return &(*place);
}

Mapping* Paths::traverse_right(Mapping* mapping) {
    // Get the iterator for this Mapping*
    list<Mapping>::iterator place = mapping_itr.at(mapping);

    // Get the path name for this Mapping*
    string path_name = mapping_path_name(mapping);

    // Get the list that the iterator is in
    list<Mapping>& path_list = _paths.at(path_name);

    // Advance the iterator right.
    place++;

    // If we're at the end, return null
    if(place == path_list.end()) {
        return nullptr;
    }

    // Else return the address of the stored Mapping.
    return &(*place);
}

vector<string> Paths::all_path_names(void) {
    vector<string> names;
    for (auto& p : mappings_by_rank) {
        names.push_back(p.first);
    }
    return names;
}

bool Paths::is_head_or_tail_node(id_t id) {
    return head_tail_nodes.count(id);
}

const string Paths::mapping_path_name(Mapping* m) {
    auto n = mapping_path.find(m);
    if (n == mapping_path.end()) {
        return "";
    } else {
        return n->second;
    }
}

map<string, int> Paths::node_path_traversal_counts(id_t id, bool rev) {
    map<string, int> path_counts;
    if (has_node_mapping(id)) {
        for (auto& p : get_node_mapping(id)) {
            path_counts[p.first]++;
        }
    } else if (is_head_or_tail_node(id)) {
        for (auto& n : all_path_names()) {
            path_counts[n]++;
        }
    }
    return path_counts;
}

vector<string> Paths::node_path_traversals(id_t id, bool rev) {
    vector<string> names;
    if (has_node_mapping(id)) {
        for (auto& p : get_node_mapping(id)) {
            names.push_back(p.first);
        }
    } else if (is_head_or_tail_node(id)) {
        names = all_path_names();
    }
    return names;
}

set<string> Paths::of_node(id_t id) {
    set<string> names;
    if (has_node_mapping(id)) {
        for (auto& p : get_node_mapping(id)) {
            names.insert(p.first);
        }
    }
    return names;
}

bool Paths::are_consecutive_nodes_in_path(id_t id1, id_t id2, const string& path_name) {
    if (of_node(id1).count(path_name) && of_node(id2).count(path_name)) {
        auto& p1 = get_node_mapping(id1);
        auto& p2 = get_node_mapping(id2);
        // is p1 directly before p2?
        vector<list<Mapping>::iterator> i1s, i2s;
        // note that this will get the first mapping in each path, not an arbitrary one
        // (we can have looping paths, so there could be several mappings per path)
        for (auto& mp : p1[path_name]) {
            i1s.push_back(mapping_itr[mp]);
        }
        for (auto& mp : p2[path_name]) {
            i2s.push_back(mapping_itr[mp]);
        }
        for (auto i1 : i1s) {
            ++i1; // increment the first node's mapping iterator
            for (auto i2 : i2s) {
                if (i1 == i2) return true;
            }
        }
    }
    return false;
}

vector<string> Paths::over_edge(id_t id1, bool rev1, id_t id2, bool rev2,
                                vector<string> following) {
    // try both ways
    auto forward = over_directed_edge(id1, rev1, id2, rev2, following);
    auto reverse = over_directed_edge(id2, !rev2, id1, !rev1, following);
    // take the union
    std::sort(forward.begin(), forward.end());
    std::sort(reverse.begin(), reverse.end());
    vector<string> continued;
    std::set_union(forward.begin(), forward.end(),
                   reverse.begin(), reverse.end(),
                   std::back_inserter(continued));
    return continued;
}

// among the set of followed paths which paths connect these two node strands?
vector<string> Paths::over_directed_edge(id_t id1, bool rev1, id_t id2, bool rev2,
                                         vector<string> following) {
    vector<string> consecutive;
    /* for future debugging
    cerr << "looking for edges " << id1 << " -> " << id2 << endl;
    cerr << "following ";
    std::for_each(following.begin(), following.end(), [](const string& s) { cerr << s << ", "; });
    cerr << endl;
    */
    
    // handle the head/tail node case
    // we treat these like catch-alls; every path that reaches them is assumed to continue on
    if (head_tail_nodes.count(id1)
        || head_tail_nodes.count(id2)) {
        return following;
    }

    // we want the mappings of the first node
    // which are on the same strand
    
    // if either of the nodes has no path, the result is the empty set
    if (!has_node_mapping(id1) || !has_node_mapping(id2)) return {};
    // otherwise, we can safely get reference to the mappings
    auto m1 = get_node_mappings_by_rank(id1);
    auto m2 = get_node_mappings_by_rank(id2);
    // iterate over the paths
    // and see how many pairs of consecutive ranks we have
    // in the expected direction between the nodes

    for (auto& p1 : m1) {
        auto& name = p1.first;
        auto& r1 = p1.second;
        // and the corresponding one
        auto p2 = m2.find(name);
        // matching path possible
        if (p2 == m2.end()) continue;
        // now iterate over the mappings in the first
        // use the traversal directions to determine the expected order
        // and check if there is a mapping at 
        //auto& n2 = p2.first;
        auto& r2 = p2->second;
        for (auto i1 : r1) {
            auto& m1 = i1.second;
            // only consider mappings that touch these traversals
            if (m1->position().is_reverse() != rev1) continue;
            auto rank1 = i1.first;
            // do we have something at the successive mapping in r2
            auto i2 = r2.find(rank1+1);
            if (i2 == r2.end()) continue;
            if (i2->second->position().is_reverse() != rev2) continue;
            consecutive.push_back(name);
        }
    }
    // if we aren't following anything, assume everything
    if (following.empty()) return consecutive;
    // otherwise find the paths which we continue following
    std::sort(following.begin(), following.end());
    std::sort(consecutive.begin(), consecutive.end());
    vector<string> continued;
    std::set_intersection(following.begin(), following.end(),
                          consecutive.begin(), consecutive.end(),
                          std::back_inserter(continued));
    return continued;
}

void parse_region(const string& target, string& name, id_t& start, id_t& end) {
    start = -1;
    end = -1;
    size_t foundFirstColon = target.find(":");
    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        name = target;
    } else {
        name = target.substr(0, foundFirstColon);
	    size_t foundRangeDash = target.find("-", foundFirstColon);
        if (foundRangeDash == string::npos) {
            start = atoi(target.substr(foundFirstColon + 1).c_str());
            end = start;
        } else {
            start = atoi(target.substr(foundFirstColon + 1, foundRangeDash - foundRangeDash - 1).c_str());
            end = atoi(target.substr(foundRangeDash + 1).c_str());
        }
    }

}

int path_to_length(const Path& path) {
    int l = 0;
    for (const auto& m : path.mapping()) {
        l += mapping_to_length(m);
    }
    return l;
}

int path_from_length(const Path& path) {
    int l = 0;
    for (const auto& m : path.mapping()) {
        l += mapping_from_length(m);
    }
    return l;
}

int mapping_to_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.to_length();
    }
    return l;
}

int mapping_from_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;

}

int softclip_start(const Mapping& mapping) {
    int from_length = 0;
    int to_length = 0;
    int i = 0;
    while (i < mapping.edit_size() && from_length == 0) {
        from_length += mapping.edit(i).from_length();
        if (from_length > 0) break;
        to_length += mapping.edit(i).to_length();
        ++i;
    }
    return to_length;
}

int softclip_end(const Mapping& mapping) {
    int from_length = 0;
    int to_length = 0;
    int i = mapping.edit_size()-1;
    while (i >= 0 && from_length == 0) {
        from_length += mapping.edit(i).from_length();
        if (from_length > 0) break;
        to_length += mapping.edit(i).to_length();
        --i;
    }
    return to_length;
}

// returns the first non-softclip position in the path
Position first_path_position(const Path& path) {
    // step through soft clips
    int i = 0;
    while (i < path.mapping_size()) {
        if (from_length(path.mapping(i))) break;
        ++i;
    }
    if (i == path.mapping_size()) { cerr << "[vg::Path] end of path without a from_length" << endl; exit(1); }
    const Mapping& mapping = path.mapping(i);
    // find the soft clip length here
    Position pos = mapping.position();
    pos.set_offset(pos.offset()+softclip_start(mapping));
    return pos;
}

Position last_path_position(const Path& path) {
    int i = path.mapping_size()-1;
    while (i >= 0) {
        if (from_length(path.mapping(i))) break;
        --i;
    }
    if (i < 0) { cerr << "[vg::Path] start of path without a from_length" << endl; exit(1); }
    const Mapping& mapping = path.mapping(i);
    // find the soft clip length here
    Position pos = mapping.position();
    pos.set_offset(pos.offset() + from_length(mapping));
    return pos;
}

int to_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.to_length();
    }
    return l;
}

int from_length(const Mapping& m) {
    int l = 0;
    for (int i = 0; i < m.edit_size(); ++i) {
        const Edit& e = m.edit(i);
        l += e.from_length();
    }
    return l;
}

Path& extend_path(Path& path1, const Path& path2) {

    auto& path1_back = path1.mapping(path1.mapping_size()-1);
    auto& path2_front = path2.mapping(0);

    // move insertions from the front of the second path onto the back of the first
    // this does not change positions or anything else

    // We can merge two mappings to the same node only if they're in the same
    // direction and at compatible positions on the node. TODO: determine that
    // and merge mappings.

    // check if we have to splice the last mapping together
    if (!path2_front.has_position() || !path1_back.has_position()) {
        // we build up the last mapping
        auto mapping = path1_back;
        // adapt unmapped paths (which look like insertions here)
        if (!path2_front.has_position() && path1_back.has_position()) {
            // If one mapping has no position, it can't reference reference sequence
            assert(mapping_from_length(path2_front) == 0);
            *mapping.mutable_position() = path1_back.position();
            // Copy the reverse flag
            mapping.mutable_position()->set_is_reverse(path1_back.position().is_reverse());
        } else if (!path1_back.has_position() && path2_front.has_position()) {
            // If one mapping has no position, it can't reference reference sequence
            assert(mapping_from_length(path1_back) == 0);
            *mapping.mutable_position() = path2_front.position();
            // Copy the reverse flag
            mapping.mutable_position()->set_is_reverse(path2_front.position().is_reverse());
        }
        // merge the edits from the second onto the last mapping
        for (size_t i = 0; i < path2_front.edit_size(); ++i) {
            *mapping.add_edit() = path2_front.edit(i);
        }
        // replace the last mapping with the merged one
        *path1.mutable_mapping(path1.mapping_size()-1) = mapping;
    } else {
        // just tack it on, it's on the next node
        *path1.add_mapping() = path2_front;
    }
    // add the rest of the mappings
    for (size_t i = 1; i < path2.mapping_size(); ++i) {
        *path1.add_mapping() = path2.mapping(i);
    }
    /*
    if (path_from_length(path1) != path_from_length(start) + path_from_length(path2)
        || path_to_length(path1) != path_to_length(start) + path_to_length(path2)) {
        cerr << "error:[vg::path.cpp] extend fails to produce a path with from_length and to_length "
             << "equal to the sum of those of its inputs" << endl
             << "path1  " << pb2json(start) << endl
             << "path2  " << pb2json(path2) << endl
             << "return " << pb2json(path1) << endl;
        exit(1);
    }
    */
    // and return a reference to the first path where we added the mappings of the second
    return path1;
}

// concatenates paths
Path concat_paths(const Path& path1, const Path& path2) {
    Path res = path1;
    //cerr << "-------------------- concat thing ------------------" << endl;
    //cerr << pb2json(path1) << endl << pb2json(path2) << endl;
    // tack on the edits from the last
    auto& path1_back = path1.mapping(path1.mapping_size()-1);
    auto& path2_front = path2.mapping(0);
    // move insertions from the front of the second path onto the back of the first
    // this does not change positions or anything else

    // We can merge two mappings to the same node only if they're in the same
    // direction and at compatible positions on the node. TODO: determine that
    // and merge mappings.

    // check if we have to splice the last mapping together
    if (!path2_front.has_position() || !path1_back.has_position()) {
        auto* mapping = res.mutable_mapping(res.mapping_size()-1);
        // adapt unmapped paths (which look like insertions here)
        if (!path2_front.has_position() && path1_back.has_position()) {
            // If one mapping has no position, it can't reference reference sequence
            assert(mapping_from_length(path2_front) == 0);
            *mapping->mutable_position() = path1_back.position();
            // Copy the reverse flag
            mapping->mutable_position()->set_is_reverse(path1_back.position().is_reverse());
        } else if (!path1_back.has_position() && path2_front.has_position()) {
            // If one mapping has no position, it can't reference reference sequence
            assert(mapping_from_length(path1_back) == 0);
            *mapping->mutable_position() = path2_front.position();
            // Copy the reverse flag
            mapping->mutable_position()->set_is_reverse(path2_front.position().is_reverse());
        }
        // merge the edits from the second onto the last mapping
        for (size_t i = 0; i < path2_front.edit_size(); ++i) {
            *mapping->add_edit() = path2_front.edit(i);
        }
    } else {
        // just tack it on, it's on the next node
        *res.add_mapping() = path2_front;
    }
    // add the rest of the mappings
    for (size_t i = 1; i < path2.mapping_size(); ++i) {
        *res.add_mapping() = path2.mapping(i);
    }
    if (path_from_length(res) != path_from_length(path1) + path_from_length(path2)
        || path_to_length(res) != path_to_length(path1) + path_to_length(path2)) {
        cerr << "error:[vg::path.cpp] concatenate fails to produce a path with from_length and to_length "
             << "equal to the sum of those of its inputs" << endl
             << "path1  " << pb2json(path1) << endl
             << "path2  " << pb2json(path1) << endl
             << "return " << pb2json(res) << endl;
        exit(1);
    }
    //cerr << ">>>>" << endl;
    //cerr << pb2json(res) << endl;
    //cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<< end " << endl;
    return res;
}

Path simplify(const Path& p) {
    Path s;
    s.set_name(p.name());
    //cerr << "simplifying " << pb2json(p) << endl;
    // loop over the mappings in the path, doing a few things
    // exclude mappings that are total deletions
    // when possible, merge a mapping with the previous mapping
    // push inserted sequences to the left
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        auto m = simplify(p.mapping(i));
        //cerr << "simplified " << pb2json(m) << endl;
        // remove wholly-deleted or empty mappings as these are redundant
        if ((m.edit_size() == 1 && edit_is_deletion(m.edit(0)))
            || m.edit_size() == 0) continue;
        if (s.mapping_size()
            && m.position().is_reverse() == s.mapping(s.mapping_size()-1).position().is_reverse()) {
            // if this isn't the first mapping
            // refer to the last mapping
            Mapping* l = s.mutable_mapping(s.mapping_size()-1);
            // split off any insertions from the start
            // and push them to the last mapping
            size_t ins_at_start = 0;
            for (size_t j = 0; j < m.edit_size(); ++j) {
                auto& e = m.edit(j);
                if (!edit_is_insertion(e)) break;
                ins_at_start += e.to_length();
            }
            // if there are insertions at the start, move them left
            if (ins_at_start) {
                auto p = cut_mapping(m, ins_at_start);
                auto& ins = p.first;
                // cerr << "insertion " << pb2json(ins) << endl;
                // take the position from the original mapping
                m = p.second;
                *m.mutable_position() = ins.position();
                // cerr << "before and after " << pb2json(ins) << " and " << pb2json(m) << endl;
                for (size_t j = 0; j < ins.edit_size(); ++j) {
                    auto& e = ins.edit(j);
                    *l->add_edit() = e;
                }
            }
            // if our last mapping has no position, but we do, merge
            if (!l->has_position() && m.has_position()) {
                *l->mutable_position() = m.position();
            // otherwise, if we end at exactly the start position of the next mapping, we can merge
            } else if (l->has_position() && m.has_position()
                       && l->position().node_id() == m.position().node_id()
                       && l->position().offset() + mapping_from_length(*l) == m.position().offset()) {
                // we can merge the current mapping onto the old one
                *l = concat_mappings(*l, m);
            } else {
                *s.add_mapping() = m;
            }
        } else {
            *s.add_mapping() = m;
        }
    }
    // we will return this path after a final cleanup
    Path r;
    r.set_name(p.name());
    // remove any edit-less mappings that may have resulted from left-shifting indels
    for (size_t i = 0; i < s.mapping_size(); ++i) {
        auto& m = s.mapping(i);
        if (!m.edit_size()) continue; // skips empty mappings
        *r.add_mapping() = m;
    }
    // now set ranks
    for (size_t i = 0; i < r.mapping_size(); ++i) {
        auto* m = r.mutable_mapping(i);
        m->set_rank(i+1);
    }
    //cerr << "simplified " << pb2json(s) << endl;
    return r;
}

// simple merge
Mapping concat_mappings(const Mapping& m, const Mapping& n) {
    Mapping c = m;
    // add the edits on
    for (size_t i = 0; i < n.edit_size(); ++i) {
        *c.add_edit() = n.edit(i);
    }
    // merge anything that's identical
    return simplify(c);
}

Mapping simplify(const Mapping& m) {
    Mapping n;
    if (m.rank()) n.set_rank(m.rank());
    //cerr << "pre simplify " << pb2json(m) << endl;
    // take the old position (which may be empty)
    *n.mutable_position() = m.position();

    size_t j = 0;
    // to simplify, we skip deletions at the very start of the node
    // these are implied by jumps in the path from other nodes
    if (m.position().offset() == 0) {
        for ( ; j < m.edit_size(); ++j) {
            if (!edit_is_deletion(m.edit(j))) {
                break;
            } else {
                if (n.position().node_id() == 0) {
                    // Complain if a Mapping has no position *and* has edit-initial
                    // deletions, which we need to remove but can't.
                    throw runtime_error(
                        "Cannot simplify Mapping with no position: need to update position when removing leading deletion");
                }
                
                // Adjust the offset by the size of the deletion.
                n.mutable_position()->set_offset(n.position().offset()
                                                 + m.edit(j).from_length());
            }
        }
    }

    // now go through the rest of the edits and see if we can merge them
    if (j < m.edit_size()) {
        Edit e = m.edit(j++);
        for ( ; j < m.edit_size(); ++j) {
            auto& f = m.edit(j);
            // if the edit types are the same, merge them
            if ((edit_is_match(e) && edit_is_match(f))
                || (edit_is_sub(e) && edit_is_sub(f))
                || (edit_is_deletion(e) && edit_is_deletion(f))
                || (edit_is_insertion(e) && edit_is_insertion(f))) {
                // will be 0 for insertions, and + for the rest
                e.set_from_length(e.from_length()+f.from_length());
                // will be 0 for deletions, and + for the rest
                e.set_to_length(e.to_length()+f.to_length());
                // will be empty for both or have sequence for both
                e.set_sequence(e.sequence() + f.sequence());
            } else {
                // mismatched types are just put on
                *n.add_edit() = e;
                e = f;
            }
        }
        // and keep the last edit
        // if it isn't a deletion
        if (!edit_is_deletion(e)) *n.add_edit() = e;
    }
    //cerr << "post simplify " << pb2json(n) << endl;
    return n;
}

Mapping merge_adjacent_edits(const Mapping& m) {

    Mapping n;
    if (m.rank()) n.set_rank(m.rank());
    // get the position
    if (m.has_position()) {
        // take the old position
        *n.mutable_position() = m.position();
    }

    // This tracks what edit we're merging into
    size_t j = 0;
    // Go through the edits and see if we can merge them
    if (j < m.edit_size()) {
        Edit e = m.edit(j++);
        for ( ; j < m.edit_size(); ++j) {
            auto& f = m.edit(j);
            // if the edit types are the same, merge them
            if ((edit_is_match(e) && edit_is_match(f))
                || (edit_is_sub(e) && edit_is_sub(f))
                || (edit_is_deletion(e) && edit_is_deletion(f))
                || (edit_is_insertion(e) && edit_is_insertion(f))) {
                // will be 0 for insertions, and + for the rest
                e.set_from_length(e.from_length()+f.from_length());
                // will be 0 for deletions, and + for the rest
                e.set_to_length(e.to_length()+f.to_length());
                // will be empty for both or have sequence for both
                e.set_sequence(e.sequence() + f.sequence());
            } else {
                // mismatched types are just put on
                *n.add_edit() = e;
                e = f;
            }
        }
        // Keep the last edit
        *n.add_edit() = e;
    }
    return n;

}

Path trim_hanging_ends(const Path& p) {

    if (p.is_circular()) {
        return p;
    }
    
    Path trimmed_path;
    int first_m = -1;
    int first_e = -1;
    int last_m = -2;
    int last_e = -2;

    // walk left-to-right until we find a match
    for (int mi = 0; mi < p.mapping_size() && first_m < 0; ++mi) {
        const Mapping& mapping = p.mapping(mi);
        for (int ei = 0; ei < mapping.edit_size() && first_e < 0; ++ei) {
            const Edit& edit = mapping.edit(ei);
            if (edit_is_match(edit)) {
                first_m = mi;
                first_e = ei;
            }
        }
    }

    // walk right-to-left until we find a match
    if (first_m >= 0) {
        for (int mi = p.mapping_size() - 1; mi >=0 && last_m < 0; --mi) {
            const Mapping& mapping = p.mapping(mi);
            for (int ei = mapping.edit_size() - 1; ei >= 0 && last_e < 0; --ei) {
                const Edit& edit = mapping.edit(ei);
                if (edit_is_match(edit)) {
                    last_m = mi;
                    last_e = ei;
                }
            }
        }
        // no reason to cross over since we search for same thing in both dirs
        assert(first_m < last_m || (first_m == last_m && first_e <= last_e));
    }

    // Make a new path beginning at first_m, first_e and ending at last_m, last_e
    Path r;
    r.set_name(p.name());
    // todo: handle length field?

    int rank = 1;
    for (int mi = first_m; mi <= last_m; ++mi) {
        const Mapping& mapping = p.mapping(mi);
        Mapping* new_mapping = r.add_mapping();
        *new_mapping->mutable_position() = mapping.position();
        new_mapping->set_rank(rank);
        int ei = mi == first_m ? first_e : 0;
        int lei = mi == last_m ? last_e : mapping.edit_size() - 1;
        for (; ei <= lei; ++ei) {
            const Edit& edit = mapping.edit(ei);
            Edit* new_edit = new_mapping->add_edit();
            *new_edit = edit;
        }
    }

    return r;
}

bool mapping_ends_in_deletion(const Mapping& m){
    return m.edit_size() >= 1 && edit_is_deletion(m.edit(m.edit_size()-1));
}

bool mapping_starts_in_deletion(const Mapping& m) {
    return m.edit_size() >= 1 && edit_is_deletion(m.edit(0));
}

bool mapping_is_total_deletion(const Mapping& m) {
    return m.edit_size() == 1 && edit_is_deletion(m.edit(0));
}

bool mapping_is_simple_match(const Mapping& m) {
    return m.edit_size() == 1 && edit_is_match(m.edit(0));
}

bool path_is_simple_match(const Path& p) {
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        if (!mapping_is_simple_match(p.mapping(i))) return false;
    }
    return true;
}

const string mapping_sequence(const Mapping& mp, const string& node_seq) {
    string seq;
    // todo reverse the mapping
    function<int64_t(id_t)> lambda = [&node_seq](id_t i){return node_seq.size();};
    Mapping m = (mp.position().is_reverse()
                 ? reverse_complement_mapping(mp, lambda) : mp);
    // then edit in the forward direction (easier)
    // and, if the mapping is reversed, finally reverse-complement the result
    size_t t = 0;
    size_t f = m.position().offset();
    for (size_t i = 0; i < m.edit_size(); ++i) {
        auto& e = m.edit(i);
        if (edit_is_match(e)) {
            seq.append(node_seq.substr(f, e.from_length()));
        } else if (edit_is_sub(e)) {
            seq.append(e.sequence());
        } else if (edit_is_insertion(e)) {
            seq.append(e.sequence());
        } else if (edit_is_deletion(e)) {
            // no-op
        }
        t += e.to_length();
        f += e.from_length();
    }
    // TODO: we must resolve these semantics
    // probably we shouldn't have perfect matches be represented this way
    // it is better to be explicit
    // perfect match
    if (m.edit_size() == 0) {
        seq = node_seq;
    }
    return (mp.position().is_reverse() ? reverse_complement(seq) : seq);
}

const string mapping_sequence(const Mapping& mp, const Node& n) {
    if (!mp.has_position()) {
        assert(mp.edit_size()==1);
        return mp.edit(0).sequence();
    }
    assert(mp.position().node_id() == n.id());
    auto& node_seq = n.sequence();
    return mapping_sequence(mp, node_seq);
}

Mapping reverse_complement_mapping(const Mapping& m,
                                   const function<int64_t(id_t)>& node_length) {
    // Make a new reversed mapping
    Mapping reversed = m;

    // switching around to the reverse strand requires us to change offsets
    // that are nonzero to count the unused bases on the other side of the block
    // of used bases.
    if(m.has_position() && m.position().node_id() != 0) {
        Position* p = reversed.mutable_position();
        
        if(m.edit_size() == 0) {
            // This is an old-style perfect match mapping with no edits.
            
            // The offset becomes unused bases after the mapping
            size_t unused_bases_after = p->offset();
            size_t used_bases = node_length(p->node_id()) - unused_bases_after;
            
            // There are now no unused bases before
            p->set_offset(0);
            // And we are on the other strand
            p->set_is_reverse(!p->is_reverse());
            
            // But we have to have a mapping in there in order to actually
            // specify that we're ending at the not-at-the-end position where we
            // used to start.
            Edit* edit = reversed.add_edit();
            edit->set_from_length(used_bases);
            edit->set_to_length(used_bases);
        } else {
            // We have edits; we may not be a full-length perfect match.
        
            // How many node bases are used by the mapping?
            size_t used_bases = mapping_from_length(m);
            // How many are taken up by the offset on the other strand?
            size_t unused_bases_after = p->offset();
            // The remainder ought to be taken up by the offset on this strand.
            size_t unused_bases_before = node_length(p->node_id()) - used_bases - unused_bases_after;
            
    #ifdef debug
            cerr << "Node " << p->node_id() << " breakdown: " << unused_bases_before << ", "
                << used_bases << ", " << unused_bases_after << endl;
    #endif
            
            // Adopt the new offset
            p->set_offset(unused_bases_before);
            // Toggle the reversed-ness flag
            p->set_is_reverse(!p->is_reverse());
            
        }
    }

    // Clear out all the edits. TODO: we wasted time copying them
    reversed.clear_edit();

    for(size_t i = m.edit_size() - 1; i != (size_t) -1; i--) {
        // For each edit in reverse order, put it in reverse complemented
        *reversed.add_edit() = reverse_complement_edit(m.edit(i));
    }

    return reversed;
}

Path reverse_complement_path(const Path& path,
                             const function<int64_t(id_t)>& node_length) {
    // Make a new reversed path
    Path reversed = path;

    // Clear out all the mappings. TODO: we wasted time copying them
    reversed.clear_mapping();

    for(size_t i = path.mapping_size() - 1; i != (size_t) -1; i--) {
        // For each mapping in reverse order, put it in reverse complemented and
        // measured from the other end of the node.
        *reversed.add_mapping() = reverse_complement_mapping(path.mapping(i), node_length);
    }
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        reversed.mutable_mapping(i)->set_rank(i+1);
    }

    return reversed;
}

// ref-relative
pair<Mapping, Mapping> cut_mapping(const Mapping& m, const Position& pos) {
    Mapping left, right;
    //cerr << "cutting mapping " << pb2json(m) << " at " << pb2json(pos) << endl;
    assert(m.has_position() && m.position().node_id());
    // the two mappings get the same rank
    left.set_rank(m.rank());
    right.set_rank(m.rank());

    // TODO: support reverse mappings
    assert(!m.position().is_reverse());

    //cerr << "cutting mapping " << pb2json(m) << " at pos " << pb2json(pos) << endl;
    // left always has the position of the input mapping
    *left.mutable_position() = m.position();
    if (m.position().node_id() != pos.node_id()) {
        left = m;
    } else { // we're on the node
        if (pos.offset() == m.position().offset()) {
            // we will get a 0-length left
            //cerr << "should get a 0-length left" << endl;
            right = m;
        } else if (pos.offset() >= mapping_from_length(m)) {
            //cerr << "should get a 0-length right" << endl;
            // or a 0-length right
            left = m;
        } else {
            //cerr << "we need to cut the mapping" << endl;
            // we need to cut the mapping
            // find the cut point and build the two mappings
            size_t seen = 0;
            size_t j = 0;
            // loop over those before our position
            for ( ; j < m.edit_size() && seen < pos.offset(); ++j) {
                auto& e = m.edit(j);
                if (seen + e.from_length() == pos.offset()) {
                    // this would be the last edit before the target position
                    // so we just drop it onto the last mapping of p1
                    //cerr << "at last edit before the target position" << endl;
                    *left.add_edit() = e;
                } else if (seen + e.from_length() > pos.offset()) {
                    //cerr << "we need to divide the edit" << endl;
                    // we need to divide this edit
                    auto edits = cut_edit_at_from(e, seen + e.from_length() - pos.offset());
                    //cerr << "got edits " << pb2json(edits.first) << " and " << pb2json(edits.second) << endl;
                    *left.add_edit() = edits.first;
                    *right.add_edit() = edits.second;
                }
                seen += e.from_length();
            }
            // now we add to the second path
            assert(seen >= pos.offset());
            for ( ; j < m.edit_size(); ++j) {
                *right.add_edit() = m.edit(j);
            }
        }
    }
    *right.mutable_position() = pos;
    assert(!m.has_position()
           || (left.has_position()
               && left.position().node_id()
               && right.has_position()
               && right.position().node_id()));
    return make_pair(left, right);
}

// mapping-relative
pair<Mapping, Mapping> cut_mapping(const Mapping& m, size_t offset) {
    Mapping left, right;
    //cerr << ".cutting mapping " << pb2json(m) << " at " << offset << endl;
    // both result mappings will be in the same orientation as the original
    left.mutable_position()->set_is_reverse(m.position().is_reverse());
    right.mutable_position()->set_is_reverse(m.position().is_reverse());
    left.set_rank(m.rank());
    right.set_rank(m.rank());

    //assert(m.has_position() && m.position().node_id());
    // left always has the position of the input mapping
    if (m.has_position()) *left.mutable_position() = m.position();
    // nothing to cut
    if (offset == 0) {
        // we will get a 0-length left
        right = m;
    } else if (offset >= mapping_to_length(m)) {
        // or a 0-length right
        left = m;
    } else {
        // we need to cut the mapping

        // find the cut point and build the two mappings
        size_t seen = 0;
        size_t j = 0;
        // loop over those before our position
        for ( ; j < m.edit_size() && seen < offset; ++j) {
            auto& e = m.edit(j);
            //cerr << "at edit " << pb2json(e) << endl;
            if (seen + e.to_length() > offset) {
                // we need to divide this edit
                auto edits = cut_edit_at_to(e, seen + e.to_length() - offset);
                *left.add_edit() = edits.first;
                *right.add_edit() = edits.second;
            } else {
                // this would be the last edit before the target position
                // so we just drop it onto the last mapping of p1
                *left.add_edit() = e;
            }
            seen += e.to_length();
        }
        // now we add to the second path
        assert(seen >= offset);
        for ( ; j < m.edit_size(); ++j) {
            *right.add_edit() = m.edit(j);
        }
    }
    if (m.has_position()) {
        // The right mapping has a position on this same node
        right.mutable_position()->set_node_id(m.position().node_id());
        right.mutable_position()->set_offset(left.position().offset()
                                             + mapping_from_length(left));
    }
    if (left.has_position()
        && !left.position().node_id()) {
        left.clear_position();
    }
    if (right.has_position()
        && !right.position().node_id()) {
        right.clear_position();
    }
    /*
    if (!(!m.has_position()
           || (left.has_position()
               && left.position().node_id()
               && right.has_position()
               && right.position().node_id()))) {
        cerr << "problem with cut alignment" << endl
             << "------left " << pb2json(left) << endl << "------right " << pb2json(right) << endl;
        assert(false);
    }
    */
    //cerr << "cut mappings " << endl
    //<< "------left " << pb2json(left) << endl << "------right " << pb2json(right) << endl;
    return make_pair(left, right);
}

// divide path at reference-relative position
// TODO make this work on the reverse strand / inverting paths
pair<Path, Path> cut_path(const Path& path, const Position& pos) {
    Path p1, p2;
    size_t i = 0;
    // seek forward to the cut point
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        //cerr << "seeking cut pos " << pb2json(pos) << " at mapping " << pb2json(m) << endl;
        // the position is in this node, so make the cut
        if (m.position().node_id() == pos.node_id()) {
            //cerr << "making cuts" << endl;
            auto mappings = cut_mapping(m, pos);
            //cerr << "left cut " << pb2json(mappings.first) << " and right " << pb2json(mappings.second) << endl;
            // and save the cuts
            *p1.add_mapping() = mappings.first;
            *p2.add_mapping() = mappings.second;
            ++i; // we don't increment our mapping index when we break here
            break;
        } else {
            // otherwise keep adding the mappings onto our first path
            *p1.add_mapping() = m;
        }
    }
    // add in the rest of the edits
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        *p2.add_mapping() = m;
    }
    assert(p1.mapping(0).has_position()
           && p1.mapping(0).position().node_id()
           && p2.mapping(0).has_position()
           && p2.mapping(0).position().node_id());
    //cerr << "---cut_path left " << pb2json(p1) << endl << "---and right " << pb2json(p2) << endl;
    return make_pair(p1, p2);
}

// divide the path at a path-relative offset as measured in to_length from start
pair<Path, Path> cut_path(const Path& path, size_t offset) {
    Path p1, p2;
    size_t seen = 0;
    size_t i = 0;
    if (!path.mapping_size()) {
        cerr << "[vg::Path] cannot cut an empty path" << endl;
        assert(false);
    }
#ifdef debug
    cerr << "cutting path at " << offset << endl << pb2json(m) << endl;
#endif
    // seek forward to the cut point
    for ( ; i < path.mapping_size() && seen < offset; ++i) {
        auto& m = path.mapping(i);
#ifdef debug
        cerr << "seeking cut offset " << offset << " at mapping " << pb2json(m) << endl;
#endif
        // the position is in this node, so make the cut
        if (seen + mapping_to_length(m) == offset) {
            *p1.add_mapping() = m;
        } else if (seen + mapping_to_length(m) > offset) {
#ifdef debug
            cerr << "making cuts" << endl;
#endif
            auto mappings = cut_mapping(m, offset - seen);
            // and save the cuts
            *p1.add_mapping() = mappings.first;
            *p2.add_mapping() = mappings.second;
#ifdef debug
            cerr << "left cut " << pb2json(mappings.first) << " and right " << pb2json(mappings.second) << endl;
#endif
            ++i; // we don't increment our mapping index when we break here
            seen += mapping_to_length(m); // same problem
            break;
        } else {
            // otherwise keep adding the mappings onto our first path
            *p1.add_mapping() = m;
        }
        seen += mapping_to_length(m);
    }
#ifdef debug
    cerr << "seen " << seen << " offset " << offset << endl;
#endif
    assert(seen >= offset);
    // add in the rest of the edits
    for ( ; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        *p2.add_mapping() = m;
    }
#ifdef debug
    cerr << "---cut_path left " << pb2json(p1) << endl << "---and right " << pb2json(p2) << endl;
#endif
    return make_pair(p1, p2);
}

bool maps_to_node(const Path& p, id_t id) {
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        if (p.mapping(i).position().node_id() == id) return true;
    }
    return false;
}

// returns the start position, or an empty position if the path has no mappings with positions
Position path_start(const Path& path) {
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        if (mapping.has_position()) return mapping.position();
    }
    Position p;
    return p; // empty
}

string path_to_string(Path p){
    stringstream ostr;
    int i;
    for (i = 0; i < p.mapping_size(); i++){
        Mapping m = p.mapping(i);
        Position pos = m.position();
        ostr << pos.node_id();
        ostr << (pos.is_reverse() ? "<-" : "->");
    }
    return ostr.str();
}

// determine the path end
Position path_end(const Path& path) {
    Position pos;
    if (!path.mapping_size()) return pos;
    auto& last = path.mapping(path.mapping_size()-1);
    pos = last.position();
    pos.set_offset(pos.offset()+mapping_from_length(last));
    return pos;
}

bool adjacent_mappings(const Mapping& m1, const Mapping& m2) {
    return abs(m1.rank() - m2.rank()) == 1;
}

bool mapping_is_match(const Mapping& m) {
    for(size_t i = 0; i < m.edit_size(); i++) {
        if(!edit_is_match(m.edit(i))) {
            // We have a non-match edit
            return false;
        }
    }
    // No non-match edits found
    return true;
}

// assumes complete description of mapped node by the edits
double divergence(const Mapping& m) {
    double from_length = 0;
    double to_length = 0;
    double matches = 0;
    double mismatches = 0;
    double insertions = 0;
    double deletions = 0;
    for (size_t i = 0; i < m.edit_size(); ++i) {
        auto& edit = m.edit(i);
        // what is the length
        from_length += edit.from_length();
        to_length += edit.to_length();
        if (edit_is_match(edit)) {
            matches += edit.to_length();
        } else if (edit_is_sub(edit)) {
            mismatches += edit.to_length();
        } else if (edit_is_insertion(edit)) {
            insertions += edit.to_length();
        } else if (edit_is_deletion(edit)) {
            deletions += edit.from_length();
        }
    }
    return 1 - (matches*2.0 / (from_length + to_length));
}

double identity(const Path& path) {
    double ident = 0;
    size_t total_length = path_to_length(path);
    size_t matched_length = 0;
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        for (size_t j = 0; j < mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
            if (edit_is_match(edit)) {
                matched_length += edit.from_length();
            }
        }
    }
    return total_length == 0 ? 0.0 : (double) matched_length / (double) total_length;
}


void
decompose(const Path& path,
          map<pos_t, int>& ref_positions,
          map<pos_t, Edit>& edits) {
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& mapping = path.mapping(i);
        pos_t ref_pos = make_pos_t(mapping.position());
        for (int j = 0; j < mapping.edit_size(); ++j) {
            auto& edit = mapping.edit(j);
            if (edit_is_match(edit)) {
                int p = offset(ref_pos);
                for (int q = p; q < edit.from_length()+p; ++q) {
                    ref_positions[ref_pos]++;
                    get_offset(ref_pos)++;
                }
            } else {
                edits[ref_pos] = edit;
                get_offset(ref_pos) += edit.from_length();
            }
        }
    }
}

double overlap(const Path& p1, const Path& p2) {
    map<pos_t, int> ref1, ref2;
    map<pos_t, Edit> edit1, edit2;
    decompose(p1, ref1, edit1);
    decompose(p2, ref2, edit2);
    // compare the two position sets
    int matching = 0;
    //int total = path_to_length(p1) + path_to_length(p2);
    // match positions from 1 in 2
    for (auto& p : ref1) {
        auto f = ref2.find(p.first);
        if (f != ref2.end()) {
            matching += min(p.second, f->second);
        }
    }
    // handle edits
    for (auto& e : edit1) {
        auto& edit = e.second;
        auto f = edit2.find(e.first);
        if (f != edit2.end()) {
            if (edit == f->second) {
                matching += edit.to_length();
            }
            edit2.erase(f);
        }
    }
    return (double) matching / (double) path_to_length(p1);
}
    
void translate_node_ids(Path& path, const unordered_map<id_t, id_t>& translator) {
    for (size_t i = 0; i < path.mapping_size(); i++) {
        Position* position = path.mutable_mapping(i)->mutable_position();
        position->set_node_id(translator.at(position->node_id()));
    }
}

void translate_oriented_node_ids(Path& path, const unordered_map<id_t, pair<id_t, bool>>& translator) {
    for (size_t i = 0; i < path.mapping_size(); i++) {
        Position* position = path.mutable_mapping(i)->mutable_position();
        const pair<id_t, bool>& translation = translator.at(position->node_id());
        position->set_node_id(translation.first);
        position->set_is_reverse(translation.second != position->is_reverse());
    }
}
    
pos_t initial_position(const Path& path) {
    if (!path.mapping_size()) {
        return pos_t();
    }
    return path.mapping_size() ? make_pos_t(path.mapping(0).position()) : pos_t();
}

pos_t final_position(const Path& path) {
    if (!path.mapping_size()) {
        return pos_t();
    }
    const Mapping& mapping = path.mapping(path.mapping_size() - 1);
    return make_pos_t(mapping.position().node_id(),
                      mapping.position().is_reverse(),
                      mapping.position().offset() + mapping_from_length(mapping) - 1);
}

Path path_from_node_traversals(const list<NodeTraversal>& traversals) {
    // We'll fill in this path
    Path toReturn;
    
    // What rank should the next node be?
    size_t rank = 1;
    
    for(const auto& traversal : traversals) {
        // Add a mapping for each NodeTraversal
        Mapping* mapping = toReturn.add_mapping();
        
        // Set up the position
        mapping->mutable_position()->set_node_id(traversal.node->id());
        mapping->mutable_position()->set_is_reverse(traversal.backward);
        
        // Set the rank
        mapping->set_rank(rank++);
        
        // Add an edit
        Edit* edit = mapping->add_edit();
        // Make it cover the full node
        edit->set_from_length(traversal.node->sequence().size());
        edit->set_to_length(traversal.node->sequence().size());
    }
    
    // We're done making the path
    return toReturn;
}

}
