#include "path.hpp"
#include <vg/io/stream.hpp>
#include "region.hpp"
#include <sstream>

using namespace vg::io;

namespace vg {

const std::function<bool(const string&)> Paths::is_alt = [](const string& path_name) {
    // Really we want things that match the regex "_alt_.+_[0-9]+"
    // But std::regex was taking loads and loads of time (probably matching .+) so we're replacing it with special-purpose code.
    
    string prefix("_alt_");
    
    if (path_name.length() < prefix.length() || !std::equal(prefix.begin(), prefix.end(), path_name.begin())) {
        // We lack the prefix
        return false;
    }
    
    // Otherwise it's almost certainly an alt, but make sure it ends with numbers after '_' to be sure.
    
    size_t found_digits = 0;
    for (auto it = path_name.rbegin(); it != path_name.rend() && *it != '_'; ++it) {
        // Scan in reverse until '_' (which we know exists)
        if (*it < '0' || *it > '9') {
            // Out of range character
            return false;
        }
        found_digits++;
    }
    
    // If there were any digits, and ony digits, it matches.
    return (found_digits > 0);
    
};

string Paths::strip_subrange(const string& path_name, subrange_t* out_subrange) {
    subrange_t subrange = PathMetadata::parse_subrange(path_name);
    string base_name;
    if (subrange == PathMetadata::NO_SUBRANGE) {
        base_name = path_name;
    } else {
        PathSense sense;
        string sample;
        string locus;
        size_t haplotype;
        size_t phase_block;
        PathMetadata::parse_path_name(path_name, sense, sample, locus, haplotype, phase_block, subrange);
        base_name = PathMetadata::create_path_name(sense, sample, locus, haplotype, phase_block, PathMetadata::NO_SUBRANGE);        
    }
    if (out_subrange) {
        *out_subrange = subrange;
    }
    return base_name;
}

mapping_t::mapping_t(void) : traversal(0), length(0), rank(1) { }

mapping_t::mapping_t(const Mapping& m) {
    traversal = m.position().node_id();
    set_is_reverse(m.position().is_reverse());
    length = mapping_from_length(m);
    // only import fully embedded mappings
    assert(length == mapping_to_length(m));
    // Note that in the case of no edits (deprecated shorthand for a full-length perfect match), length will be 0. 
    rank = m.rank();
}

Mapping mapping_t::to_mapping(void) const {
    Mapping m;
    Position* p = m.mutable_position();
    p->set_node_id(node_id());
    p->set_is_reverse(is_reverse());
    if (length != 0) {
        // We're storing an Edit that knows it length, not a deprecated shorthand full-length match.
        Edit* e = m.add_edit();
        e->set_to_length(length);
        e->set_from_length(length);
    }
    m.set_rank(rank);
    return m;
}

id_t mapping_t::node_id(void) const {
    return abs(traversal);
}

void mapping_t::set_node_id(id_t id) {
    bool is_rev = is_reverse();
    traversal = id * (is_rev ? -1 : 1);
}

bool mapping_t::is_reverse(void) const {
    return traversal < 0;
}

void mapping_t::set_is_reverse(bool is_rev) {
    traversal = abs(traversal) * (is_rev ? -1 : 1);
}

ostream& operator<<(ostream& out, mapping_t mapping) {
    return out << mapping.node_id() << " " << (mapping.is_reverse() ? "rev" : "fwd");
}

Paths::Paths(void) {
    max_path_id = 0;
    // noop
}

void Paths::load(istream& in) {
    function<void(Path&)> lambda = [this](Path& p) {
        this->extend(p);
    };
    vg::io::for_each(in, lambda);
}

void Paths::write(ostream& out) {
    vector<string> path_names;
    for (auto& p : _paths) {
        const string& name = p.first;
        path_names.push_back(name);
    }
    function<Path(size_t)> lambda =
        [this, &path_names](size_t i) -> Path {
        auto& mappings = _paths[path_names.at(i)];
        Path path;
        for (auto& m : mappings) {
            *path.add_mapping() = m.to_mapping();
        }
        path.set_name(path_names.at(i));
        if (circular.count(path_names.at(i))) {
            path.set_is_circular(true);
        }
        return path;
    };
    vg::io::write(out, _paths.size(), lambda);
    vg::io::finish(out);
}

void Paths::to_graph(Graph& g) {
    for (auto& p : _paths) {
        const string& name = p.first;
        auto& mappings = p.second;
        Path* path = g.add_path();
        path->set_name(name);
        if (circular.count(name)) {
            path->set_is_circular(true);
        }
        for (auto& m : mappings) {
            *path->add_mapping() = m.to_mapping();
        }
    }
}

Path Paths::path(const string& name) {
    Path path;
    auto p = _paths.find(name);
    if (p == _paths.end()) {
        return path;
    }
    auto& mappings = p->second;
    path.set_name(name);
    for (auto& m : mappings) {
        *path.add_mapping() = m.to_mapping();
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

void Paths::for_each_name(const function<void(const string&)>& lambda) const {
    for (auto& p : _paths) {
        const string& name = p.first;
        lambda(name);
    }
}

bool Paths::for_each_name_stoppable(const function<bool(const string&)>& lambda) const {
    for (auto& p : _paths) {
        const string& name = p.first;
        if (!lambda(name)) {
            return false;
        }
    }
    return true;
}

void Paths::for_each_mapping(const function<void(mapping_t&)>& lambda) {
    for (auto& p : _paths) {
        auto& path = p.second;
        for (auto& m : path) {
            lambda(m);
        }
    }
}

void Paths::for_each_stream(istream& in, const function<void(Path&)>& lambda) {
    vg::io::for_each(in, lambda);
}

void Paths::make_circular(const string& name) {
    circular.insert(name);
}

void Paths::make_linear(const string& name) {
    circular.erase(name);
}

void Paths::extend(const Path& p, bool warn_on_duplicates, bool rebuild_indexes) {
    const string& name = p.name();
    // Make sure we preserve empty paths
    get_create_path(name);
    for (int i = 0; i < p.mapping_size(); ++i) {
        const Mapping& m = p.mapping(i);
        append_mapping(name, m, warn_on_duplicates);
    }
    if (p.is_circular()) {
        make_circular(name);
    }
    if (rebuild_indexes) {
        // re-sort?
        sort_by_mapping_rank();
        rebuild_mapping_aux();
    }
}

// one of these should go away
void Paths::extend(const Paths& paths, bool warn_on_duplicates, bool rebuild_indexes) {
    for (auto& p : paths._paths) {
        const string& name = p.first;
        auto& path = p.second;
        // Make sure we preserve empty paths
        get_create_path(name);
        for (auto& m : path) {
            append_mapping(name, m.to_mapping(), warn_on_duplicates);
        }
        if (paths.circular.count(name)) {
            make_circular(name);
        }
    }
    if (rebuild_indexes) {
        sort_by_mapping_rank();
        rebuild_mapping_aux();
    }
}

void Paths::extend(const vector<Path> & paths, bool warn_on_duplicates, bool rebuild_indexes) {
    for (auto& p : paths) {
        extend(p, warn_on_duplicates, false);
    }
    if (rebuild_indexes) {
        sort_by_mapping_rank();
        rebuild_mapping_aux();
    }
}

void Paths::append(const Paths& paths, bool warn_on_duplicates, bool rebuild_indexes) {
    extend(paths, warn_on_duplicates, rebuild_indexes);
}

void Paths::append(const Graph& g, bool warn_on_duplicates, bool rebuild_indexes) {
    for (int i = 0; i < g.path_size(); ++i) {
        // Make sure we preserve empty paths
        extend(g.path(i), warn_on_duplicates, false);
    }
    if (rebuild_indexes) {
        sort_by_mapping_rank();
        rebuild_mapping_aux();
    }
}

Path& append_path(Path& a, const Path& b) {
    a.mutable_mapping()->MergeFrom(b.mapping());
    return a;
}

bool Paths::has_mapping(const string& name, int32_t rank) {
    auto iter = mappings_by_rank.find(name);
    if (iter != mappings_by_rank.end()) {
        return iter->second.count(rank);
    }
    return false;
}

void Paths::append_mapping(const string& name, const mapping_t& m, bool warn_on_duplicates) {
    // get or create the path with this name
    list<mapping_t>& pt = get_create_path(name);
    // now if we haven't already supplied a mapping
    // add it
    
    if (!m.rank || !has_mapping(name, m.rank)) {
        // If we don't have a rank set or we don't have a mapping in this path
        // with that rank, we need to add the mapping.
        
        // First figure out what the previous rank on the path is.
        size_t last_rank = 0;
        if(!pt.empty()) {
            last_rank = (*pt.rbegin()).rank;
        }

        // Add this mapping at the end of the path. Note that this may not
        // actually be where it belongs according to its rank; if that is the
        // case, sort_by_mapping_rank() and rebuild_mapping_aux() need to be
        // called, after all the mappings are loaded, in order to put them in
        // order by rank.
        pt.push_back(m);
        mapping_t* mp = &pt.back();
        
        if(mp->rank == 0) {
            // 0 rank defaults to being ranked at the end of what's already
            // there. After all, that's where we just added it.
            if(&pt.front() != mp) {
                // There are other things on the path
                if(last_rank && !mappings_by_rank[name].count(last_rank + 1)) {
                    // We can go right after the thing that we're physically
                    // after.
                    mp->rank = last_rank + 1;
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
                mp->rank = 1;
            }
        }
        
        // add it to the node mappings
        auto& ms = get_node_mapping(mp->node_id());
        ms[get_path_id(name)].insert(mp);
        // and record its position in this list
        list<mapping_t>::iterator mi = pt.end(); --mi;
        auto& mitr = mapping_itr[mp];
        mitr.first = mi;
        mitr.second = get_path_id(name);
        if(mp->rank) {
            // Only if we actually end up with a rank (i.e. all the existing
            // ranks weren't cleared) do we really index by rank.
            mappings_by_rank[name][mp->rank] = mp;
        }
    } else if (warn_on_duplicates) {
        // This mapping duplicates the rank of an existing mapping.
        // We're not going to keep it, so we should complain.
        cerr << "[vg] warning: path " << name << " rank " << m.rank << " appears multiple times. Skipping." << endl;
    }
}

int64_t Paths::get_path_id(const string& name) const {
    int64_t path_id;
#pragma omp critical (path_id_map)
    {
        // in order to keep the critical section inside above if (so it's only touched when initializing)
        // we need the second check here
        if (!name_to_id.count(name)) {
            // Assign an ID.
            // These members are mutable.
            ++max_path_id;
            id_to_name[max_path_id] = name;
            name_to_id[name] = max_path_id;
        }
        path_id = name_to_id[name];
    }
    return path_id;
}

const string& Paths::get_path_name(int64_t id) const {
    const string* name;
#pragma omp critical (path_id_map)
    {
        name = &id_to_name[id];
    }
    return *name;
}

void Paths::append_mapping(const string& name, id_t id, bool is_reverse, size_t length, size_t rank, bool warn_on_duplicates) {
    mapping_t m;
    m.set_node_id(id);
    m.set_is_reverse(is_reverse);
    m.length = length;
    
    // If the rank passed in is 0, it will get filled in by the other version of
    // append_mapping.
    m.rank = rank;
    
    append_mapping(name, m, warn_on_duplicates);
}

void Paths::prepend_mapping(const string& name, const Mapping& m, bool warn_on_duplicates) {
    // get or create the path with this name
    list<mapping_t>& pt = get_create_path(name);
   
    // TODO: I'm not sure if this is the best way for handling ranks, but the ranks
    // are really a chunked serialization thing, not an in-memory construct. Moreover,
    // we're ideally going to move away from using the VG graph in the future, so I don't
    // expect this will even come up. Mostly just trying to meet the HandleGraph interface
    // in the interim.
    int32_t rank = m.rank();
    if (rank == 0) {
        // no given rank, decrement the first rank, skipping over 0 to preserve it as
        // a sentinel
        if (pt.empty()) {
            rank = 1;
        }
        else if (pt.front().rank != 1) {
            rank = pt.front().rank - 1;
        }
        else {
            rank = -1;
        }
    }
    
    // now if we haven't already supplied a mapping
    // add it
    if (!has_mapping(name, rank)) {
        // If we don't have a rank set or we don't have a mapping in this path
        // with that rank, we need to add the mapping.
        
        pt.push_front(mapping_t(m));
        mapping_t* mp = &pt.front();
        // add it to the node mappings
        auto& ms = get_node_mapping(m.position().node_id());
        ms[get_path_id(name)].insert(mp);
        // and record its position in this list
        list<mapping_t>::iterator mi = pt.begin();
        auto& mitr = mapping_itr[mp];
        mitr.first = mi;
        mitr.second = get_path_id(name);
        mappings_by_rank[name][mp->rank] = mp;
    } else if (warn_on_duplicates) {
        // This mapping duplicates the rank of an existing mapping.
        // We're not going to keep it, so we should complain.
        cerr << "[vg] warning: path " << name << " rank " << rank << " appears multiple times. Skipping." << endl;
    }
}

void Paths::prepend_mapping(const string& name, id_t id, bool is_reverse, size_t length, size_t rank, bool warn_on_duplicates) {
    mapping_t m;
    m.set_node_id(id);
    m.set_is_reverse(is_reverse);
    m.length = length;
    m.rank = rank;
    
    prepend_mapping(name, m.to_mapping(), warn_on_duplicates);
}

size_t Paths::get_next_rank(const string& name) {
    auto& p = get_path(name);
    //cerr << "next rank be " << p.size()+1 << " or " << (size_t) p.rend()->rank()+1 << endl;
    return max(p.size()+1, (size_t) (p.size() ? p.rend()->rank+1 : 0));
}

// these will split a mapping into two
// NB: each submapping ends up with the same rank as the parent
// however, they will be ordered correctly in the path
// we will need to normalize path ranks to make this right
pair<mapping_t*, mapping_t*> Paths::divide_mapping(mapping_t* m, const Position& pos) {
    // this is needed to split mappinsg during e.g. normalization
    // but still ensure that the mappings are out there
    // what do we do?
    // first we take the mapping and divide it as we do
    auto n = cut_mapping(*m, pos);
    return replace_mapping(m, n);
}

pair<mapping_t*, mapping_t*> Paths::divide_mapping(mapping_t* m, size_t offset) {
    auto n = cut_mapping(*m, offset);
    return replace_mapping(m, n);
}

pair<mapping_t*, mapping_t*> Paths::replace_mapping(mapping_t* m, pair<mapping_t, mapping_t> n) {
    // then we remove it from the node it's pointing to
    // and replace it with the other two mappings
    // we'll give them the same rank, but record them in the right order
    // this leaves an invalid graph
    // there are a few ways to fix this--- they involve changing the way we record ranks
    // but for now it's going to be simplest if the calling context manages this
    auto& path_name = mapping_path_name(m);
    n.first.rank = m->rank;
    n.second.rank = m->rank;
    // remove the mapping, getting an iterator pointing to the element that was after it
    if (!m->is_reverse()) {
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

bool Paths::has_path(const string& name) const {
    return _paths.find(name) != _paths.end();
}

void Paths::increment_node_ids(id_t inc) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<mapping_t>& path = p.second;
        for (auto& m : path) {
            m.set_node_id(m.node_id()+inc);
        }
    }
    rebuild_node_mapping();
}

void Paths::swap_node_ids(const std::function<nid_t(const nid_t&)>& get_new_id) {
    for (auto& p : _paths) {
        const string& name = p.first;
        list<mapping_t>& path = p.second;
        for (auto& m : path) {
            // Look up the replacement ID
            auto replacement = get_new_id(m.node_id());
            if(replacement != 0) {
                // If there is a nonzero replacement, use it.
                m.set_node_id(replacement);
            }
        }
    }
    rebuild_node_mapping();
}

void Paths::swap_node_ids(hash_map<id_t, id_t>& id_mapping) {
    swap_node_ids([&](const nid_t& id) -> nid_t {
        auto it = id_mapping.find(id);
        if (it == id_mapping.end()) {
            // Not found
            return 0;
        } else {
            // Use the result
            return it->second;
        }
    });
}

void Paths::reassign_node(id_t new_id, mapping_t* m) {
    // erase the old node id
    node_mapping[m->node_id()][mapping_path_id(m)].erase(m);
    // set the new node id
    m->set_node_id(new_id);
    // and record it in the new node record
    node_mapping[m->node_id()][mapping_path_id(m)].insert(m);
}

void Paths::rebuild_node_mapping(void) {
    // starts with paths and rebuilds the index
    node_mapping.clear();
    for (auto& p : _paths) {
        const string& path_name = p.first;
        list<mapping_t>& path = p.second;
        for (auto& m : path) {
            get_node_mapping(m.node_id())[get_path_id(path_name)].insert(&m);
        }
    }
}

// attempt to sort the paths based on the recorded ranks of the mappings
void Paths::sort_by_mapping_rank(void) {
    for (auto p = _paths.begin(); p != _paths.end(); ++p) {
        list<mapping_t>& path = p->second;
        path.sort([](const mapping_t& m1, const mapping_t& m2) {
                return m1.rank < m2.rank;
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
    mappings_by_rank.clear();
    for (auto& p : _paths) {
        const string& path_name = p.first;
        list<mapping_t>& path = p.second;
        size_t order_in_path = 0;
        for (list<mapping_t>::iterator i = path.begin(); i != path.end(); ++i) {
            auto& mitr = mapping_itr[&*i];
            mitr.first = i;
            mitr.second = get_path_id(path_name);
            
            if(i->rank > order_in_path + 1) {
                // Make sure that if we have to assign a rank to a node after
                // this one, it is greater than this node's rank. TODO: should
                // we just uniformly re-rank all the nodes starting at 0? Or
                // will we ever want to cut and paste things back together using
                // the old preserved ranks?
                order_in_path = i->rank - 1;
            }
            
            if (i->rank == 0 || i->rank < order_in_path + 1) {
                // If we don't already have a rank, or if we see a rank that
                // can't be correct given the ranks we have already seen, we set
                // the rank based on what we've built
                i->rank = order_in_path+1;
            }
            
            // Save the mapping as being at the given rank in its path.
            mappings_by_rank[path_name][i->rank] = &*i;
            
            ++order_in_path;
        }
    }
}

void Paths::remove_node(id_t id) {
    node_mapping.erase(id);
}

list<mapping_t>::iterator Paths::find_mapping(mapping_t* m) {
    return mapping_itr[m].first;
}

list<mapping_t>::iterator Paths::remove_mapping(mapping_t* m) {
    // The mapping has to exist
    assert(mapping_itr.find(m) != mapping_itr.end());
    auto& mitr = mapping_itr[m];
    const string& path_name = get_path_name(mitr.second);
    id_t id = m->node_id();
    auto& x = _paths[path_name];
    
    // This gets tricky because we're going to deallocate the storage pointed to
    // by m. We need to remove it from other things first.
    if(m->rank && mappings_by_rank[path_name].count(m->rank) && 
        mappings_by_rank[path_name][m->rank] == m) {
        // If we have this node stored for its path and rank, kick it out.
        mappings_by_rank[path_name].erase(m->rank);
    }
    
    // Actually deallocate the mapping
    list<mapping_t>::iterator p = _paths[path_name].erase(mitr.first);
    if (has_node_mapping(id)) {
        auto& node_path_mapping = get_node_mapping(id);
        node_path_mapping[mitr.second].erase(m);
        if (node_path_mapping.empty()) node_mapping.erase(id);
    }
    mapping_itr.erase(m);
    
    return p;
}

list<mapping_t>::iterator Paths::insert_mapping(list<mapping_t>::iterator w, const string& path_name, const mapping_t& m) {
    auto px = _paths.find(path_name);
    assert(px != _paths.end());
    list<mapping_t>& path = px->second;
    list<mapping_t>::iterator p;
    if (path.empty()) {
        path.push_front(m);
        p = path.begin();
    } else {
        p = path.insert(w, m);
    }
    get_node_mapping(m.node_id())[get_path_id(path_name)].insert(&*p);
    auto& mitr = mapping_itr[&*p];
    mitr.first = p;
    mitr.second = get_path_id(path_name);
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
    mappings_by_rank.clear();
}

void Paths::clear_mapping_ranks(void) {
    for (auto p = _paths.begin(); p != _paths.end(); ++p) {
        list<mapping_t>& path = p->second;
        for (auto m = path.begin(); m != path.end(); ++m) {
            mapping_t& mapping = *m;
            mapping.rank = 0;
        }
    }
    mappings_by_rank.clear();
}

list<mapping_t>& Paths::get_path(const string& name) {
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
        if(node_mapping.count(mapping.node_id())) {
            // Throw out all the mappings for this path on this node
            node_mapping[mapping.node_id()].erase(get_path_id(name));
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

list<mapping_t>& Paths::create_path(const string& name) {
    return _paths[name];
}

list<mapping_t>& Paths::get_create_path(const string& name) {
    if (!has_path(name)) {
        return create_path(name);
    } else {
        return get_path(name);
    }
}

bool Paths::has_node_mapping(id_t id) const {
    return node_mapping.find(id) != node_mapping.end();
}

bool Paths::has_node_mapping(Node* n) const {
    return node_mapping.find(n->id()) != node_mapping.end();
}

map<int64_t, set<mapping_t*>>& Paths::get_node_mapping(id_t id) {
    return node_mapping[id];
}
    
const map<int64_t, set<mapping_t*>>& Paths::get_node_mapping(id_t id) const {
    return node_mapping.at(id);
}

map<int64_t, set<mapping_t*>>& Paths::get_node_mapping(Node* n) {
    return node_mapping[n->id()];
}

map<string, set<mapping_t*>> Paths::get_node_mapping_by_path_name(id_t id) {
    auto& nm = node_mapping[id];
    map<string, set<mapping_t*>> mm;
    for (auto& p : nm) {
        mm[get_path_name(p.first)] = p.second;
    }
    return mm;
}

map<string, set<mapping_t*>> Paths::get_node_mapping_by_path_name(Node* n) {
    return get_node_mapping_by_path_name(n->id());
}

map<string, map<int, mapping_t*>> Paths::get_node_mappings_by_rank(id_t id) {
    map<string, map<int, mapping_t*>> by_ranks;
    for (auto& p : get_node_mapping(id)) {
        auto& name = get_path_name(p.first);
        auto& mp = p.second;
        for (auto* m : mp) by_ranks[name][m->rank] = m;
    }
    return by_ranks;
}

map<string, map<int, mapping_t>> Paths::get_node_mapping_copies_by_rank(id_t id) {
    map<string, map<int, mapping_t>> by_ranks;
    for (auto& p : get_node_mapping(id)) {
        auto& name = get_path_name(p.first);
        auto& mp = p.second;
        for (auto* m : mp) by_ranks[name][m->rank] = *m;
    }
    return by_ranks;
}

mapping_t* Paths::traverse_left(mapping_t* mapping) {
    // Get the iterator for this Mapping*
    auto& mitr = mapping_itr.at(mapping);
    list<mapping_t>::iterator place = mitr.first;

    // Get the path name for this Mapping*
    const string& path_name = get_path_name(mitr.second);

    // Get the list that the iterator is in
    list<mapping_t>& path_list = _paths.at(path_name);

    // If we're already the beginning, return null.
    if(place == path_list.begin()) {
        return nullptr;
    }

    // Else walk left and return the address of the stored Mapping. std::list
    // iterators are bidirectional, so we will be able to do it.
    place--;
    return &(*place);
}

mapping_t* Paths::traverse_right(mapping_t* mapping) {
    // Get the iterator for this Mapping*
    auto& mitr = mapping_itr.at(mapping);
    list<mapping_t>::iterator place = mitr.first;

    // Get the path name for this Mapping*
    const string& path_name = get_path_name(mitr.second);

    // Get the list that the iterator is in
    list<mapping_t>& path_list = _paths.at(path_name);

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

const string Paths::mapping_path_name(mapping_t* m) {
    int64_t id = mapping_path_id(m);
    if (id == 0) {
        return "";
    } else {
        return get_path_name(id);
    }
}

int64_t Paths::mapping_path_id(mapping_t* m) {
    auto n = mapping_itr.find(m);
    if (n == mapping_itr.end()) {
        return 0;
    } else {
        return n->second.second;
    }
}

map<string, int> Paths::node_path_traversal_counts(id_t id, bool rev) {
    map<string, int> path_counts;
    if (has_node_mapping(id)) {
        for (auto& p : get_node_mapping(id)) {
            path_counts[get_path_name(p.first)]++;
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
            names.push_back(get_path_name(p.first));
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
            names.insert(get_path_name(p.first));
        }
    }
    return names;
}

bool Paths::are_consecutive_nodes_in_path(id_t id1, id_t id2, const string& path_name) {
    if (of_node(id1).count(path_name) && of_node(id2).count(path_name)) {
        auto& p1 = get_node_mapping(id1);
        auto& p2 = get_node_mapping(id2);
        // is p1 directly before p2?
        vector<list<mapping_t>::iterator> i1s, i2s;
        // note that this will get the first mapping in each path, not an arbitrary one
        // (we can have looping paths, so there could be several mappings per path)
        for (auto& mp : p1[get_path_id(path_name)]) {
            i1s.push_back(mapping_itr[mp].first);
        }
        for (auto& mp : p2[get_path_id(path_name)]) {
            i2s.push_back(mapping_itr[mp].first);
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
            if (m1->is_reverse() != rev1) continue;
            auto rank1 = i1.first;
            // do we have something at the successive mapping in r2
            auto i2 = r2.find(rank1+1);
            if (i2 == r2.end()) continue;
            if (i2->second->is_reverse() != rev2) continue;
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
    // set path ranks
    for (size_t i = 0; i < path1.mapping_size(); ++i) {
        auto* m = path1.mutable_mapping(i);
        m->set_rank(i+1);
    }
    // and return a reference to the first path where we added the mappings of the second
    return path1;
}

// concatenates paths
Path concat_paths(const Path& path1, const Path& path2) {
    
    if (path1.mapping_size() == 0) {
        return path2;
    } else if (path2.mapping_size() == 0) {
        return path1;
    }
    
    // Otherwise there are mappings in both and we have real work to do

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
    return simplify(res);
}

Path simplify(const Path& p, bool trim_internal_deletions) {
    Path s;
    s.set_name(p.name());
    //cerr << "simplifying " << pb2json(p) << endl;
    // loop over the mappings in the path, doing a few things
    // exclude mappings that are total deletions
    // when possible, merge a mapping with the previous mapping
    // push inserted sequences to the left
    for (size_t i = 0; i < p.mapping_size(); ++i) {
        auto m = simplify(p.mapping(i), trim_internal_deletions);
        // remove empty mappings as these are redundant
        if (trim_internal_deletions) {
            // remove wholly-deleted or empty mappings as these are redundant
            if ((m.edit_size() == 1 && edit_is_deletion(m.edit(0)))
                || m.edit_size() == 0) continue;
        } else {
            // remove empty mappings as these are redundant
            if (m.edit_size() == 0) continue;
        }
        if (s.mapping_size()) {
            //&& m.position().is_reverse() == s.mapping(s.mapping_size()-1).position().is_reverse()) {
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
            if ((!l->has_position() || l->position().node_id() == 0)
                && (m.has_position() && m.position().node_id() != 0)) {
                *l->mutable_position() = m.position();
                // if our last mapping has a position, and we don't, merge
            } else if ((!m.has_position() || m.position().node_id() == 0)
                       && (l->has_position() && l->position().node_id() != 0)) {
                *m.mutable_position() = *l->mutable_position();
                m.mutable_position()->set_offset(from_length(*l));
            }
            // if we end at exactly the start position of the next mapping, we can merge
            if ((!l->has_position() && !m.has_position())
                ||
                (l->has_position() && m.has_position()
                 && l->position().is_reverse() == m.position().is_reverse()
                 && l->position().node_id() == m.position().node_id()
                 && l->position().offset() + mapping_from_length(*l) == m.position().offset())) {
                // we can merge the current mapping onto the old one
                *l = concat_mappings(*l, m, trim_internal_deletions);
            } else {
                if (from_length(m) || to_length(m)) {
                    *s.add_mapping() = m;
                }
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
        // remove position if it's empty
        auto& l = *r.mutable_mapping(r.mapping_size()-1);
        if (l.has_position() && l.position().node_id()==0) {
            l.clear_position();
        }
    }
    Path q;
    // remove leading and trailing deletions (these might result from global alignment)
    int total_to_length = path_to_length(r);
    int seen_to_length = 0;
    for (size_t i = 0; i < r.mapping_size(); ++i) {
        auto& m = r.mapping(i);
        int curr_to_length = mapping_to_length(m);
        // skip bits at the beginning and end
        if ((!seen_to_length && !curr_to_length)
            || seen_to_length == total_to_length) continue;
        Mapping n;
        *n.mutable_position() = m.position();
        if (seen_to_length) {
            if (seen_to_length + curr_to_length == total_to_length) {
                // this is the last mapping before we should trim
                // so trim any dels from the end of the mapping
                size_t j = m.edit_size()-1;
                for ( ; j >= 0; --j) {
                    if (!edit_is_deletion(m.edit(j))) {
                        ++j; // pointer to end
                        break;
                    }
                }
                for (size_t h = 0; h < j; ++h) {
                    *n.add_edit() = m.edit(h);
                }
            } else {
                // internal segment
                n = m;
            }
        } else {
            // our first matching mapping
            if (mapping_to_length(m)) {
                size_t j = 0;
                size_t seen = 0;
                for ( ; j < m.edit_size(); ++j) {
                    if (!edit_is_deletion(m.edit(j))) {
                        break;
                    } else {
                        seen += m.edit(j).from_length();
                    }
                }
                // adjust position
                n.mutable_position()->set_offset(n.position().offset()+seen);
                for ( ; j < m.edit_size(); ++j) {
                    *n.add_edit() = m.edit(j);
                }
            }
        }
        *q.add_mapping() = n;
        seen_to_length += mapping_to_length(n);
    }
    q.set_name(r.name());
    assert(path_to_length(q) == path_to_length(r));

    // now set ranks and clear empty positions and edits
    for (size_t i = 0; i < q.mapping_size(); ++i) {
        auto& m = q.mapping(i);
        Mapping n;
        //auto* m = q.mutable_mapping(i);
        n.set_rank(i+1);
        if (m.position().node_id() != 0) {
            // this is an empty position, so let's remove it
            *n.mutable_position() = m.position();
        }
        for (auto& e : m.edit()) {
            if (!edit_is_empty(e)) {
                *n.add_edit() = e;
            }
        }
        *q.mutable_mapping(i) = n;
    }

    return q;
}

// simple merge
Mapping concat_mappings(const Mapping& m, const Mapping& n, bool trim_internal_deletions) {
    Mapping c = m;
    // add the edits on
    for (size_t i = 0; i < n.edit_size(); ++i) {
        *c.add_edit() = n.edit(i);
    }
    // merge anything that's identical
    return simplify(c, trim_internal_deletions);
}

Mapping simplify(const Mapping& m, bool trim_internal_deletions) {
    Mapping n;
    if (m.rank()) n.set_rank(m.rank());
    // take the old position (which may be empty)
    if (m.has_position()) {
        *n.mutable_position() = m.position();
    }

    size_t j = 0;
    if (trim_internal_deletions) {
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
    }
    
    // go through the rest of the edits and see if we can merge them
    if (j < m.edit_size()) {
        Edit e = m.edit(j++);
        for ( ; j < m.edit_size(); ++j) {
            auto& f = m.edit(j);
            // if the edit types are the same, merge them
            if (edit_is_empty(f)) {
                continue;
            } else if (edits_are_compatible(e, f)) {
                merge_edits_in_place(e, f);
            } else {
                // mismatched types are just put on
                *n.add_edit() = e;
                e = f;
            }
        }
        if (trim_internal_deletions) {
            if (!edit_is_deletion(e)) *n.add_edit() = e;
        } else {
            *n.add_edit() = e;
        }
    }
    return n;
}

bool edits_are_compatible(const Edit& e, const Edit& f) {
    return (edit_is_match(e) && edit_is_match(f))
            || (edit_is_sub(e) && edit_is_sub(f))
            || (edit_is_deletion(e) && edit_is_deletion(f))
            || (edit_is_insertion(e) && edit_is_insertion(f));
}

void merge_edits_in_place(Edit& e, const Edit& f) {
    // will be 0 for insertions, and + for the rest
    e.set_from_length(e.from_length() + f.from_length());
    // will be 0 for deletions, and + for the rest
    e.set_to_length(e.to_length() + f.to_length());
    // will be empty for both or have sequence for both
    e.set_sequence(e.sequence() + f.sequence());
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
        new_mapping->set_rank(rank++);
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

bool mappings_equivalent(const Mapping& m1, const Mapping& m2) {
    bool equivalent = (m1.position().node_id() == m2.position().node_id()
                       && m1.position().is_reverse() == m2.position().is_reverse()
                       && m1.position().offset() == m2.position().offset()
                       && m1.edit_size() == m2.edit_size());
    for (size_t i = 0; i < m1.edit_size() && equivalent; ++i) {
        const auto& e1 = m1.edit(i);
        const auto& e2 = m2.edit(i);
        equivalent = (e1.from_length() == e2.from_length()
                      && e1.to_length() == e2.to_length()
                      && e1.sequence() == e2.sequence());
    }
    return equivalent;
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

bool mapping_is_total_insertion(const Mapping& m) {
    return m.edit_size() == 1 && edit_is_insertion(m.edit(0));
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
    if (!mp.has_position() || !mp.position().node_id()) {
        // With no grap position we must be a pure insert.
        // n is undefined.
        // But we might have multiple edits.
        std::stringstream s;
        for (auto& e : mp.edit()) {
            // We can't have any from bases if we have no graph position.
            assert(e.from_length() == 0);
            s << e.sequence();
        }
        return s.str();
    }
    assert(mp.position().node_id() == n.id());
    auto& node_seq = n.sequence();
    return mapping_sequence(mp, node_seq);
}

// convert the path to a sequence
string path_sequence(const HandleGraph& graph, const Path& path) {
    string seq;
    for (int i = 0; i < path.mapping_size(); ++i) {
        auto& m = path.mapping(i);
        handle_t h = graph.get_handle(m.position().node_id(), m.position().is_reverse());
        seq.append(graph.get_sequence(h));
    }
    return seq;
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

    // Clear out all the edits. TODO: we wasted time copying them
    reversed.clear_edit();

    for (int64_t i = m.edit_size() - 1; i >= 0; i--) {
        // For each edit in reverse order, put it in reverse complemented
        *reversed.add_edit() = reverse_complement_edit(m.edit(i));
    }

    return reversed;
}
    
void reverse_complement_mapping_in_place(Mapping* m,
                                         const function<int64_t(id_t)>& node_length) {
        
    Position* pos = m->mutable_position();
    pos->set_is_reverse(!pos->is_reverse());
    pos->set_offset(node_length(pos->node_id()) - pos->offset() - mapping_from_length(*m));
    
    size_t swap_size = m->edit_size() / 2;
    for (size_t i = 0, j = m->edit_size() - 1; i < swap_size; i++, j--) {
        Edit* e1 = m->mutable_edit(i);
        Edit* e2 = m->mutable_edit(j);
        
        int64_t from_length_tmp = e1->from_length();
        int64_t to_length_tmp = e1->to_length();
        string sequence_tmp = e1->sequence();
        
        e1->set_from_length(e2->from_length());
        e1->set_to_length(e2->to_length());
        e1->set_sequence(reverse_complement(e2->sequence()));
        
        e2->set_from_length(from_length_tmp);
        e2->set_to_length(to_length_tmp);
        e2->set_sequence(reverse_complement(sequence_tmp));
    }
    
    
    if (m->edit_size() % 2) {
        Edit* e = m->mutable_edit(swap_size);
        reverse_complement_in_place(*e->mutable_sequence());
    }
}

Path reverse_complement_path(const Path& path,
                             const function<int64_t(id_t)>& node_length) {

    // Make a new reversed path
    Path reversed = path;

    // Clear out all the mappings. TODO: we wasted time copying them
    reversed.clear_mapping();

    for(int64_t i = path.mapping_size() - 1; i >= 0; i--) {
        // For each mapping in reverse order, put it in reverse complemented and
        // measured from the other end of the node.
        *reversed.add_mapping() = reverse_complement_mapping(path.mapping(i), node_length);
    }
    for (size_t i = 0; i < path.mapping_size(); ++i) {
        reversed.mutable_mapping(i)->set_rank(i+1);
    }

    return reversed;
}

void reverse_complement_path_in_place(Path* path,
                                      const function<int64_t(id_t)>& node_length) {
    
    size_t swap_size = path->mapping_size() / 2;
    for (size_t i = 0, j = path->mapping_size() - 1; i < swap_size; i++, j--) {
        Mapping* m1 = path->mutable_mapping(i);
        Mapping* m2 = path->mutable_mapping(j);
        
        reverse_complement_mapping_in_place(m1, node_length);
        reverse_complement_mapping_in_place(m2, node_length);
        
        int64_t rank_tmp = m1->rank();
        m1->set_rank(m2->rank());
        m2->set_rank(rank_tmp);
        
        std::swap(*m1, *m2);
    }
    
    if (path->mapping_size() % 2) {
        reverse_complement_mapping_in_place(path->mutable_mapping(swap_size), node_length);
    }
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

pair<mapping_t, mapping_t> cut_mapping(const mapping_t& m, const Position& pos) {
    auto r = cut_mapping(m.to_mapping(), pos);
    return make_pair(mapping_t(r.first), mapping_t(r.second));
}

// ref-relative with offset
pair<Mapping, Mapping> cut_mapping_offset(const Mapping& m, size_t offset) {
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
    } else if (offset >= mapping_from_length(m)) {
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
            if (seen + e.from_length() > offset) {
                // we need to divide this edit
                auto edits = cut_edit_at_from(e, seen + e.from_length() - offset);
                *left.add_edit() = edits.first;
                *right.add_edit() = edits.second;
            } else {
                // this would be the last edit before the target position
                // so we just drop it onto the last mapping of p1
                *left.add_edit() = e;
            }
            seen += e.from_length();
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

pair<mapping_t, mapping_t> cut_mapping_offset(const mapping_t& m, size_t offset) {
    auto r = cut_mapping_offset(m.to_mapping(), offset);
    return make_pair(mapping_t(r.first), mapping_t(r.second));
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

pair<mapping_t, mapping_t> cut_mapping(const mapping_t& m, size_t offset) {
    auto r = cut_mapping(m.to_mapping(), offset);
    return make_pair(mapping_t(r.first), mapping_t(r.second));
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
        return make_pair(path, path);
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
Position path_start_position(const Path& path) {
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
Position path_end_position(const Path& path) {
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
    if (p1.mapping_size() == 0 || p2.mapping_size() == 0) return 0;
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

void translate_node_ids(Path& path, const unordered_map<id_t, id_t>& translator, id_t cut_node, size_t bases_removed, bool from_right) {
    // First just translate the IDs
    translate_node_ids(path, translator);
    
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        // Scan the whole path again. We can't count on the cut node only being in the first and last mappings.
        Position* position = path.mutable_mapping(i)->mutable_position();
        if (position->node_id() == cut_node) {
            // Then adjust offsets to account for the cut on the original node
            
            // If the position in the path is counting from the same end of the
            // node that we didn't keep after the cut, we have to bump up its
            // offset.
            if ((!position->is_reverse() && !from_right) || // We cut off the left of the node, and we're counting from the left
                (position->is_reverse() && from_right)) { // We cut off the right of the node, and we're counting from the right
                // Update the offset to reflect the removed bases
                position->set_offset(position->offset() + bases_removed);
            }
        }
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

void translate_oriented_node_ids(Path& path, const function<pair<id_t, bool>(id_t)>& translator) {
    for (size_t i = 0; i < path.mapping_size(); i++) {
        Position* position = path.mutable_mapping(i)->mutable_position();
        const pair<id_t, bool>& translation = translator(position->node_id());
        position->set_node_id(translation.first);
        position->set_is_reverse(translation.second != position->is_reverse());
    }
}


void translate_node_ids(path_t& path, const unordered_map<id_t, id_t>& translator) {
    for (size_t i = 0; i < path.mapping_size(); i++) {
        position_t* position = path.mutable_mapping(i)->mutable_position();
        position->set_node_id(translator.at(position->node_id()));
    }
}
void translate_oriented_node_ids(path_t& path, const unordered_map<id_t, pair<id_t, bool>>& translator) {
    for (size_t i = 0; i < path.mapping_size(); i++) {
        position_t* position = path.mutable_mapping(i)->mutable_position();
        const pair<id_t, bool>& translation = translator.at(position->node_id());
        position->set_node_id(translation.first);
        position->set_is_reverse(translation.second != position->is_reverse());
    }
}

void translate_oriented_node_ids(path_t& path, const function<pair<id_t, bool>(id_t)>& translator) {
    for (size_t i = 0; i < path.mapping_size(); i++) {
        position_t* position = path.mutable_mapping(i)->mutable_position();
        const pair<id_t, bool>& translation = translator(position->node_id());
        position->set_node_id(translation.first);
        position->set_is_reverse(translation.second != position->is_reverse());
    }
}
    
pos_t initial_position(const Path& path) {
    pos_t pos;
    if (path.mapping_size()) {
        const Position& position = path.mapping(0).position();
        get_id(pos) = position.node_id();
        get_is_rev(pos) = position.is_reverse();
        get_offset(pos) = position.offset();
    }
    return pos;
}

pos_t final_position(const Path& path) {
    pos_t pos;
    if (path.mapping_size()) {
        const Mapping& mapping = path.mapping(path.mapping_size() - 1);
        const Position& position = mapping.position();
        get_id(pos) = position.node_id();
        get_is_rev(pos) = position.is_reverse();
        get_offset(pos) = position.offset() + mapping_from_length(mapping);
    }
    return pos;
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

void remove_paths(Graph& graph, const function<bool(const string&)>& paths_to_take, std::list<Path>* matching) {

    std::list<Path> non_matching;
    for (size_t i = 0; i < graph.path_size(); i++) {
        if (paths_to_take(graph.path(i).name())) {
            if (matching != nullptr) {
                matching->push_back(graph.path(i));
            }
        } else {
            non_matching.push_back(graph.path(i));
        }
    }
    graph.clear_path();

    for (Path& path : non_matching) {
        *(graph.add_path()) = path;
    }
}

Path path_from_path_handle(const PathHandleGraph& graph, path_handle_t path_handle) {
    Path path;
    path.set_name(graph.get_path_name(path_handle));
    size_t rank = 1;
    for (handle_t handle : graph.scan_path(path_handle)) {
        Mapping* mapping = path.add_mapping();
        mapping->mutable_position()->set_node_id(graph.get_id(handle));
        mapping->mutable_position()->set_is_reverse(graph.get_is_reverse(handle));
        mapping->set_rank(rank++);
        Edit* edit = mapping->add_edit();
        edit->set_from_length(graph.get_length(handle));
        edit->set_to_length(graph.get_length(handle));
    }
    return path;
}

// Wrap a Path in an Alignment
Alignment alignment_from_path(const HandleGraph& graph, const Path& path) {
    Alignment aln;
    *aln.mutable_path() = path;
    aln.set_name(aln.path().name());
    aln.set_sequence(path_sequence(graph, path));
    return aln;
}

void from_proto_edit(const Edit& proto_edit, edit_t& edit) {
    edit.set_from_length(proto_edit.from_length());
    edit.set_to_length(proto_edit.to_length());
    edit.set_sequence(proto_edit.sequence());
}

void to_proto_edit(const edit_t& edit, Edit& proto_edit) {
    proto_edit.set_from_length(edit.from_length());
    proto_edit.set_to_length(edit.to_length());
    proto_edit.set_sequence(edit.sequence());
}

void from_proto_mapping(const Mapping& proto_mapping, path_mapping_t& mapping) {
    const auto& position = proto_mapping.position();
    auto position_copy = mapping.mutable_position();
    position_copy->set_node_id(position.node_id());
    position_copy->set_offset(position.offset());
    position_copy->set_is_reverse(position.is_reverse());
    for (const auto& edit : proto_mapping.edit()) {
        from_proto_edit(edit, *mapping.add_edit());
    }
}

void to_proto_mapping(const path_mapping_t& mapping, Mapping& proto_mapping) {
    const auto& position = mapping.position();
    auto position_copy = proto_mapping.mutable_position();
    position_copy->set_node_id(position.node_id());
    position_copy->set_offset(position.offset());
    position_copy->set_is_reverse(position.is_reverse());
    for (const auto& edit : mapping.edit()) {
        to_proto_edit(edit, *proto_mapping.add_edit());
    }
}

void from_proto_path(const Path& proto_path, path_t& path) {
    for (const auto& mapping : proto_path.mapping()) {
        from_proto_mapping(mapping, *path.add_mapping());
    }
}
void to_proto_path(const path_t& path, Path& proto_path) {
    for (const auto& mapping : path.mapping()) {
        auto mapping_copy = proto_path.add_mapping();
        to_proto_mapping(mapping, *mapping_copy);
        mapping_copy->set_rank(proto_path.mapping_size());
    }
}

int mapping_from_length(const path_mapping_t& mapping) {
    int length = 0;
    for (const auto& edit : mapping.edit()) {
        length += edit.from_length();
    }
    return length;
}

int path_from_length(const path_t& path) {
    int length = 0;
    for (const auto& mapping : path.mapping()) {
        length += mapping_from_length(mapping);
    }
    return length;
}

int mapping_to_length(const path_mapping_t& mapping) {
    int length = 0;
    for (const auto& edit : mapping.edit()) {
        length += edit.to_length();
    }
    return length;
}

int path_to_length(const path_t& path) {
    int length = 0;
    for (const auto& mapping : path.mapping()) {
        length += mapping_to_length(mapping);
    }
    return length;
}


void reverse_complement_mapping_in_place(path_mapping_t* m,
                                         const function<int64_t(id_t)>& node_length) {
    
    position_t* pos = m->mutable_position();
    pos->set_is_reverse(!pos->is_reverse());
    pos->set_offset(node_length(pos->node_id()) - pos->offset() - mapping_from_length(*m));
    
    size_t swap_size = m->edit_size() / 2;
    for (size_t i = 0, j = m->edit_size() - 1; i < swap_size; i++, j--) {
        edit_t* e1 = m->mutable_edit(i);
        edit_t* e2 = m->mutable_edit(j);
        
        int64_t from_length_tmp = e1->from_length();
        int64_t to_length_tmp = e1->to_length();
        string sequence_tmp = e1->sequence();
        
        e1->set_from_length(e2->from_length());
        e1->set_to_length(e2->to_length());
        e1->set_sequence(reverse_complement(e2->sequence()));
        
        e2->set_from_length(from_length_tmp);
        e2->set_to_length(to_length_tmp);
        e2->set_sequence(reverse_complement(sequence_tmp));
    }
    
    
    if (m->edit_size() % 2) {
        edit_t* e = m->mutable_edit(swap_size);
        reverse_complement_in_place(*e->mutable_sequence());
    }
}

path_mapping_t reverse_complement_mapping(const path_mapping_t& m,
                                          const function<int64_t(id_t)>& node_length) {
    
    path_mapping_t reversed;
    position_t* rev_pos = reversed.mutable_position();
    rev_pos->set_node_id(m.position().node_id());
    rev_pos->set_is_reverse(!m.position().is_reverse());
    rev_pos->set_offset(node_length(m.position().node_id()) - m.position().offset() - mapping_from_length(m));
    
    for (int64_t i = m.edit_size() - 1; i >= 0; i--) {
        const edit_t& e = m.edit(i);
        edit_t* rev_edit = reversed.add_edit();
        rev_edit->set_from_length(e.from_length());
        rev_edit->set_to_length(e.to_length());
        rev_edit->set_sequence(reverse_complement(e.sequence()));
    }
    
    return reversed;
}

path_t reverse_complement_path(const path_t& path,
                               const function<int64_t(id_t)>& node_length) {
    
    // Make a new reversed path
    path_t reversed;
    
    for (int64_t i = path.mapping_size() - 1; i >= 0; i--) {
        // For each mapping in reverse order, put it in reverse complemented and
        // measured from the other end of the node.
        *reversed.add_mapping() = reverse_complement_mapping(path.mapping(i), node_length);
    }
    
    return reversed;
}

void reverse_complement_path_in_place(path_t* path,
                                      const function<int64_t(id_t)>& node_length) {
    
    size_t swap_size = path->mapping_size() / 2;
    for (size_t i = 0, j = path->mapping_size() - 1; i < swap_size; i++, j--) {
        path_mapping_t* m1 = path->mutable_mapping(i);
        path_mapping_t* m2 = path->mutable_mapping(j);
        
        reverse_complement_mapping_in_place(m1, node_length);
        reverse_complement_mapping_in_place(m2, node_length);
        
        std::swap(*m1, *m2);
    }
    
    if (path->mapping_size() % 2) {
        reverse_complement_mapping_in_place(path->mutable_mapping(swap_size), node_length);
    }
}

pos_t initial_position(const path_t& path) {
    pos_t pos;
    if (path.mapping_size()) {
        const position_t& position = path.mapping(0).position();
        get_id(pos) = position.node_id();
        get_is_rev(pos) = position.is_reverse();
        get_offset(pos) = position.offset();
    }
    return pos;
}

pos_t final_position(const path_t& path) {
    pos_t pos;
    if (path.mapping_size()) {
        const path_mapping_t& mapping = path.mapping(path.mapping_size() - 1);
        const position_t& position = mapping.position();
        get_id(pos) = position.node_id();
        get_is_rev(pos) = position.is_reverse();
        get_offset(pos) = position.offset() + mapping_from_length(mapping);
    }
    return pos;
}

string debug_string(const path_t& path) {
    string to_return = "{";
    if (!path.mapping().empty()) {
        to_return += "mapping: [";
        for (size_t i = 0; i < path.mapping_size(); ++i) {
            if (i > 0) {
                to_return += ", ";
            }
            to_return += debug_string(path.mapping(i));
        }
        to_return += "]";
    }
    to_return += "}";
    return to_return;
}

string debug_string(const path_mapping_t& mapping) {
    string to_return = "{pos: " + debug_string(mapping.position());
    if (!mapping.edit().empty()) {
        to_return += ", edit: [";
        for (size_t i = 0; i < mapping.edit_size(); ++i) {
            if (i > 0) {
                to_return += ", ";
            }
            to_return += debug_string(mapping.edit(i));
        }
        to_return += "]";
    }
    to_return += "}";
    return to_return;
}

string debug_string(const edit_t& edit) {
    string to_return = "{fl: " + to_string(edit.from_length()) + ", tl: " + to_string(edit.to_length());
    if (!edit.sequence().empty()) {
        to_return += ", seq: " + edit.sequence();
    }
    to_return += "}";
    return to_return;
}

int corresponding_length_internal(const path_t& path, int given_length, bool is_from_length, bool from_end) {
    int from_length = 0;
    if (path.mapping().empty()) {
        return from_length;
    }
    int incr, i_begin;
    if (from_end) {
        i_begin = path.mapping_size() - 1;
        incr = -1;
    }
    else {
        incr = 1;
        i_begin = 0;
    }
    int remaining = given_length;
    int other_length_total = 0;
    for (int i = i_begin; i >= 0 && i < path.mapping_size() && remaining != 0; i += incr) {
        const auto& mapping = path.mapping(i);
        int j_begin = from_end ? mapping.edit_size() - 1 : 0;
        for (int j = j_begin; j >= 0 && j < mapping.edit_size() && remaining != 0; j += incr) {
            const edit_t& edit = mapping.edit(j);
            int walking_length, other_length;
            if (is_from_length) {
                walking_length = edit.from_length();
                other_length = edit.to_length();
            }
            else {
                walking_length = edit.to_length();
                other_length = edit.from_length();
            }
            if (remaining >= walking_length) {
                remaining -= walking_length;
                other_length_total += other_length;
            }
            else {
                other_length_total += (remaining * other_length) / walking_length;
                remaining = 0;
            }
        }
    }
    return other_length_total;
}
int corresponding_to_length(const path_t& path, int from_length, bool from_end) {
    return corresponding_length_internal(path, from_length, true, from_end);
}

int corresponding_from_length(const path_t& path, int to_length, bool from_end) {
    return corresponding_length_internal(path, to_length, false, from_end);
}

}
