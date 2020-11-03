#include "genotypekit.hpp"
#include "cactus.hpp"
#include "traversal_finder.hpp"

//#define debug

namespace vg {

using namespace std;

SnarlTraversal get_traversal_of_snarl(VG& graph, const Snarl* snarl, const SnarlManager& manager, const Path& path) {

    // We'll fill this in
    SnarlTraversal to_return;

    auto contents = manager.deep_contents(snarl, graph, true);

    for(size_t i = 0; i < path.mapping_size(); i++) {
        const Mapping& mapping = path.mapping(i);

        if(contents.first.count(mapping.position().node_id())) {
            // We're inside the bubble. This is super simple when we have the contents!
            *to_return.add_visit() = to_visit(mapping, true);
        }
    }

    return to_return;
}

string traversal_to_string(VG& graph, const SnarlTraversal& path) {
    string seq;
    for (const auto& visit : path.visit()) {
        // For every visit
        if (visit.node_id() != 0) {
            // If it's to a node, but the node's sequenece
            const Node* node = graph.get_node(visit.node_id());
            seq += visit.backward() ? reverse_complement(node->sequence()) : node->sequence();
        } else {
            // Put a description of the child snarl
            stringstream s;
            s << "(" << visit << ")";
            seq += s.str();
        }
    }
    return seq;
}


pair<const Edge*, bool> AugmentedGraph::base_edge(const Edge* edge) {
    assert(base_graph != NULL);

    bool is_trivial = false;
    const Edge* found_edge = NULL;

    // check if from node is even in the base graph
    Position from_pos;
    from_pos.set_node_id(edge->from());
    if (!translator.has_translation(from_pos)) {
        return pair<const Edge*, bool>(NULL, false);
    }

    // check if the to node is even in the base graph
    Position to_pos;
    to_pos.set_node_id(edge->to());
    if (!translator.has_translation(to_pos)) {
        return pair<const Edge*, bool>(NULL, false);        
    }

    // work in forward strand since translator doesn't seem strand-aware
    from_pos.set_is_reverse(false);
    from_pos.set_offset(edge->from_start() ? 0 :
                        graph.get_node(edge->from())->sequence().length() - 1);

    to_pos.set_is_reverse(false);
    to_pos.set_offset(edge->to_end() == false ? 0 :
                      graph.get_node(edge->to())->sequence().length() - 1);

    // Map to the base graph using our translation table
    Position base_from_pos = translator.translate(from_pos);
    Position base_to_pos = translator.translate(to_pos);

    // Make sure we don't have an unexpected reversal
    assert(base_from_pos.is_reverse() == false);
    assert(base_to_pos.is_reverse() == false);

    // Test if we're a trivial edge in the base graph (just to consecutive positions)
    if (base_from_pos.node_id() == base_to_pos.node_id() &&
        base_from_pos.is_reverse() == base_to_pos.is_reverse() &&
        abs(base_from_pos.offset() - base_to_pos.offset()) == 1) {
        is_trivial = true;
    }
    // 2) or an existing edge in the base graph (translations are all forward strand,
    //    so we need to look up sequence lengths to check if we're coming out end)
    else if (edge->from_start() && !edge->to_end() &&
             base_from_pos.offset() == 0 &&
             base_to_pos.offset() == 0) {
        found_edge = base_graph->get_edge(NodeSide(base_from_pos.node_id(), false),
                                          NodeSide(base_to_pos.node_id(), false));
    }
    else if (edge->from_start() && edge->to_end() &&
             base_from_pos.offset() == 0 &&
             base_to_pos.offset() == base_graph->get_node(base_to_pos.node_id())->sequence().length() - 1) {
        found_edge = base_graph->get_edge(NodeSide(base_from_pos.node_id(), false),
                                          NodeSide(base_to_pos.node_id(), true));
    }
    else if (!edge->from_start() && !edge->to_end() &&
             base_from_pos.offset() == base_graph->get_node(base_from_pos.node_id())->sequence().length() - 1 &&
             base_to_pos.offset() == 0) {
        found_edge = base_graph->get_edge(NodeSide(base_from_pos.node_id(), true),
                                          NodeSide(base_to_pos.node_id(), false));
    }
    else if (!edge->from_start() && edge->to_end() &&
             base_from_pos.offset() == base_graph->get_node(base_from_pos.node_id())->sequence().length() - 1 &&
             base_to_pos.offset() == base_graph->get_node(base_to_pos.node_id())->sequence().length() - 1) {
        found_edge = base_graph->get_edge(NodeSide(base_from_pos.node_id(), true),
                                          NodeSide(base_to_pos.node_id(), true));
    }

    return pair<const Edge*, bool>(found_edge, is_trivial);
}

vector<const Alignment*> AugmentedGraph::get_alignments(id_t node_id) const {
    if (alignments_by_node.count(node_id)) {
        auto& found = alignments_by_node.at(node_id);
        return vector<const Alignment*>{found.begin(), found.end()};
    } else {
        return {};
    }
}

vector<const Alignment*> AugmentedGraph::get_alignments(pair<NodeSide, NodeSide> edge) const {
    auto it = alignments_by_edge.find(minmax(edge.first, edge.second));
    if (it != alignments_by_edge.end()) {
        return vector<const Alignment*>{it->second.begin(), it->second.end()};
    } else {
        return {};
    }
}

Support AugmentedGraph::get_support(id_t node) {
    Support support;
    support.set_forward(get_alignments(node).size());
    return support;
}

Support AugmentedGraph::get_support(edge_t edge) {
    Support support;
    NodeSide from(graph.get_id(edge.first), !graph.get_is_reverse(edge.first));
    NodeSide to(graph.get_id(edge.second), graph.get_is_reverse(edge.second));
    support.set_forward(get_alignments(make_pair(from, to)).size());
    return support;
}

bool AugmentedGraph::has_supports() const {
    return !alignments_by_edge.empty() || !alignments_by_node.empty();
}

vector<const Alignment*> AugmentedGraph::get_alignments() const {
    vector<const Alignment*> to_return;
    to_return.reserve(embedded_alignments.size());
    for (auto& alignment : embedded_alignments) {
        to_return.push_back(&alignment);
    }
    return to_return;
}

void AugmentedGraph::clear() {
    // Reset to default state
    *this = AugmentedGraph();
}

void AugmentedGraph::augment_from_alignment_edits(vector<Alignment>& alignments,
                                                  bool unique_names, bool leave_edits) {

    // This should only be called once.
    assert(embedded_alignments.empty());
    
    if (unique_names) { 
        // Make sure they have unique names.
        set<string> names_seen;
        // We warn about duplicate names, but only once.
        bool duplicate_names_warned = false;
        for(size_t i = 0; i < alignments.size(); i++) {
            if(alignments[i].name().empty()) {
                // Generate a name
                alignments[i].set_name("_unnamed_alignment_" + to_string(i));
            }
            if(names_seen.count(alignments[i].name())) {
                // This name is duplicated
                if(!duplicate_names_warned) {
                    // Warn, but only once
                    cerr << "Warning: duplicate alignment names present! Example: " << alignments[i].name() << endl;
                    duplicate_names_warned = true;
                }

                // Generate a new name
                // TODO: we assume this is unique
                alignments[i].set_name("_renamed_alignment_" + to_string(i));
                assert(!names_seen.count(alignments[i].name()));
            }
            names_seen.insert(alignments[i].name());
        }
        names_seen.clear();
    }
    
    for(auto& alignment : alignments) {
        // Trim the softclips off of every read
        // Work out were to cut
        int cut_start = softclip_start(alignment);
        int cut_end = softclip_end(alignment);
        // Cut the sequence and quality
        alignment.set_sequence(alignment.sequence().substr(cut_start, alignment.sequence().size() - cut_start - cut_end));
        if (alignment.quality().size() != 0) {
            alignment.set_quality(alignment.quality().substr(cut_start, alignment.quality().size() - cut_start - cut_end));
        }
        // Trim the path
        *alignment.mutable_path() = trim_hanging_ends(alignment.path());
    }
    
    if (!leave_edits) {
        // We want to actually modify the graph to encompass these reads.
        
        // To make the edits and copy them to/from the Alignments, we need a vector of just Paths.
        // TODO: improve this interface!
        vector<Path> paths;
        paths.reserve(alignments.size());
        for (auto& alignment : alignments) {
            paths.push_back(alignment.path());
        }
    
        // Run them through vg::edit() to modify the graph, but don't embed them
        // as paths. Update the paths in place, and save the translations.
        vector<Translation> augmentation_translations;
        graph.edit(paths, &augmentation_translations, false, true, false);
        
        for (size_t i = 0; i < paths.size(); i++) {
            // Copy all the modified paths back.
            *alignments[i].mutable_path() = paths[i];
        }
        paths.clear();
        
        // Send out the translation
        translator.load(augmentation_translations);
    } else {
        // No need to add edits from the reads to the graph. We may have done it already with vg augment.
#ifdef debug
        cerr << "Don't add edits!" << endl;
#endif
    }
    
    // Steal the alignments for ourselves.
    std::swap(alignments, embedded_alignments);
    
    // Prepare the index from node ID to alignments that touch the node
    for (auto& alignment : embedded_alignments) {
        // For every alignment, grab its path
        const auto& path = alignment.path();
        
        // Keep track of the nodes we already saw it visit so we don't add it twice for a node.
        unordered_set<id_t> seen;
        // and edges
        unordered_set<pair<NodeSide, NodeSide>> seen_edges;

        NodeSide prev_side;
        for (const auto& mapping : path.mapping()) {
            // For each mapping
            if (!seen.count(mapping.position().node_id())) {
                // If it is to a new node, remember that the alignment touches the node
                seen.insert(mapping.position().node_id());
                alignments_by_node[mapping.position().node_id()].push_back(&alignment);
            }
            NodeSide cur_side = NodeSide(mapping.position().node_id(), mapping.position().is_reverse());
            if (prev_side.node > 0) {
                pair<NodeSide, NodeSide> edge(minmax(prev_side, cur_side));
                if (!seen_edges.count(edge)) {
                    seen_edges.insert(edge);
                    alignments_by_edge[edge].push_back(&alignment);
                }
            }
            // our path traverses in one end and out the other
            prev_side = cur_side.flip();
        }
    }
}

void AugmentedGraph::load_translations(istream& in_file) {
    translator.translations.clear();
    function<void(Translation&)> lambda = [&](Translation& translation) {
        translator.translations.push_back(translation);
    };
    vg::io::for_each(in_file, lambda);
    translator.build_position_table();
}

void AugmentedGraph::write_translations(ostream& out_file) {
    vg::io::write_buffered(out_file, translator.translations, 0);
}

void SupportAugmentedGraph::clear() {
    // Reset to default state
    *this = SupportAugmentedGraph();
}

bool SupportAugmentedGraph::has_supports() const {
    return !node_supports.empty() || !edge_supports.empty();
}

Support SupportAugmentedGraph::get_support(id_t node) {
    return node_supports.count(node) ? node_supports.at(node) : Support();
}

Support SupportAugmentedGraph::get_support(edge_t edge) {
    return edge_supports.count(edge) ? edge_supports.at(edge) : Support();
}

void SupportAugmentedGraph::load_supports(istream& in_file) {
    // This loads LocationSupport objects. We use them instead of pileups.
    // TODO: We need a way to view them with vg view
    node_supports.clear();
    edge_supports.clear();
    function<void(LocationSupport&)> lambda = [&](LocationSupport& location_support) {
#ifdef debug
        cerr << pb2json(location_support) << endl;
#endif
        if (location_support.oneof_location_case() == LocationSupport::kNodeId) {
            node_supports[location_support.node_id()] = location_support.support();
        } else {
            const Edge& edge = location_support.edge();
            edge_t edge_handle = graph.edge_handle(graph.get_handle(edge.from(), edge.from_start()),
                                                   graph.get_handle(edge.to(), edge.to_end()));
            edge_supports[edge_handle] = location_support.support();
        }
    };
    vg::io::for_each(in_file, lambda);    
}

void SupportAugmentedGraph::load_pack_as_supports(const string& pack_file_name, const HandleGraph* vectorizable_graph) {
    Packer packer(vectorizable_graph);
    packer.load_from_file(pack_file_name);
    vectorizable_graph->for_each_handle([&](const handle_t& handle) {
            Position pos;
            pos.set_node_id(vectorizable_graph->get_id(handle));
            size_t sequence_offset = packer.position_in_basis(pos);
            size_t total_coverage = 0;
            size_t node_length = vectorizable_graph->get_length(handle);
            for (size_t i = 0; i < node_length; ++i) {
                total_coverage += packer.coverage_at_position(sequence_offset + i);
            }
            double avg_coverage = node_length > 0 ? (double)total_coverage / node_length : 0.;
            Support support;
            // we just get one value and put it in "forward".  can't fill out the rest of the Support object. 
            support.set_forward(avg_coverage);
            node_supports[vectorizable_graph->get_id(handle)] = support;
        });
    vectorizable_graph->for_each_edge([&](const edge_t& handle_edge) {
            Edge edge;
            edge.set_from(vectorizable_graph->get_id(handle_edge.first));
            edge.set_from_start(vectorizable_graph->get_is_reverse(handle_edge.first));
            edge.set_to(vectorizable_graph->get_id(handle_edge.second));
            edge.set_to_end(vectorizable_graph->get_is_reverse(handle_edge.second));
            Support support;
            support.set_forward(packer.edge_coverage(edge));
            edge_supports[graph.edge_handle(graph.get_handle(edge.from(), edge.from_start()),
                                            graph.get_handle(edge.to(), edge.to_end()))] = support;
            return true;
        });
}

void SupportAugmentedGraph::write_supports(ostream& out_file) {
    vector<LocationSupport> buffer;
    for (auto& node_support : node_supports) {
        LocationSupport location_support;
        *location_support.mutable_support() = node_support.second;
        location_support.set_node_id(node_support.first);
        buffer.push_back(location_support);
        vg::io::write_buffered(out_file, buffer, 500);
    }
    for (auto& edge_support : edge_supports) {
        LocationSupport location_support;
        *location_support.mutable_support() = edge_support.second;
        Edge edge;
        edge.set_from(graph.get_id(edge_support.first.first));
        edge.set_from_start(graph.get_is_reverse(edge_support.first.first));
        edge.set_to(graph.get_id(edge_support.first.second));
        edge.set_to_end(graph.get_is_reverse(edge_support.first.second));
        *location_support.mutable_edge() = edge;
        buffer.push_back(location_support);
        vg::io::write_buffered(out_file, buffer, 500);
    }
    vg::io::write_buffered(out_file, buffer, 0);
}


SimpleConsistencyCalculator::~SimpleConsistencyCalculator(){

}

vector<bool> SimpleConsistencyCalculator::calculate_consistency(const Snarl& site,
        const vector<SnarlTraversal>& traversals, const Alignment& read) const {

            std::function<set<int64_t>(Alignment a, SnarlTraversal s)> shared_sites = [&](Alignment a, SnarlTraversal s){
                set<int64_t> aln_ids;
                set<int64_t> trav_ids;
                for (int i = 0; i < a.path().mapping_size(); i++){
                    Mapping m = a.path().mapping(i);
                    Position pos = m.position();
                    if (pos.node_id() != 0){
                        aln_ids.insert(pos.node_id());
                    }
                }

                for (int i = 0; i < s.visit_size(); ++i){
                    Visit v = s.visit(i);
                    if (v.node_id() != 0){
                        trav_ids.insert(v.node_id());
                    }

                }
                vector<int64_t> ret;
                std::set_intersection(aln_ids.begin(), aln_ids.end(),
                        trav_ids.begin(), trav_ids.end(),
                        std::back_inserter(ret));
                
                return set<int64_t>(ret.begin(), ret.end());
            };


            //create a consistency bool for each traversal (i.e. possible path / theoretical allele)
            vector<bool> consistencies(traversals.size());

            // For each traversal
            for (int i = 0; i < traversals.size(); ++i){
                SnarlTraversal trav = traversals[i];
                // Our snarltraversals run forward (i.e. along increading node_ids)
                // Our Alignment path may run forward OR backward.
                // We can check if an alignment is on the reverse strand and
                // flip its path around to match our snarltraversal direction.

                // Cases of consistency:
                // 1. A read maps to either end of the snarl but no internal nodes
                // 2. A read maps to both ends of the snarl but no internal nodes
                // 3. A read maps to one end of the snarl and some internal nodes.
                // 4. A read maps to both ends of the snarl and some internal nodes.
                // 5. A read maps to internal nodes, but not the snarl ends
                // A read may map to a node multiple times, or it may skip a node
                // and put an insert there.
                set<int64_t> common_ids = shared_sites(read, trav);
                bool maps_to_front = common_ids.count(trav.visit(0).node_id());
                bool maps_to_end = common_ids.count(trav.visit(trav.visit_size() - 1).node_id());
                bool maps_internally = false;

                Path read_path;
                std::function<Path(Path)> reverse_path = [&](Path p){
                    Path ret;
                    for (int i = p.mapping_size() - 1; i >= 0; i--){
                        Mapping* new_mapping = ret.add_mapping();
                        Position* pos = new_mapping->mutable_position();
                        pos->set_node_id(p.mapping(i).position().node_id());
                        int offset = p.mapping(i).position().offset();
                        pos->set_offset(offset);
                        pos->set_is_reverse(!p.mapping(i).position().is_reverse());
                        for (int j = 0; j < p.mapping(i).edit_size(); j++){
                            Edit* new_edit = new_mapping->add_edit();
                        }
                    }
                    return ret;
                };
                if (false){
                    read_path = reverse_path(read.path());
                }
                else{
                    read_path = read.path();
                }

                bool is_forward = true;
                bool is_right = true;

                if ((common_ids.size() > 1 && (maps_to_front | maps_to_end)) ||
                        common_ids.size() > 2){
                        maps_internally = true;
                }

                if (maps_to_front && maps_to_end && maps_internally){
                    // The read is anchored on both ends of the Snarl. Check
                    // the internal nodes of the Path for matches against the SnarlTraversal..

                    consistencies[i] = true;
                }
                else if (maps_to_front && maps_to_end){
                    // the read may represent a deletion,
                    // which may be in our list of traversals
                    // Either way, it's consistent with a valid traversal

                    if (true){
                        consistencies[i] = true;
                    }
                    else{
                        consistencies[i] = false;
                    }

                }
                else if (maps_to_front && maps_internally){
                    // The read maps to either the first or last node of the Snarl
                    // Check its internal nodes for matches between path and snarl


                    consistencies[i] = true;
                }
                else if (maps_to_end && maps_internally){


                    consistencies[i] = true;
                }
                else if (maps_to_front | maps_to_end){
                    // maps to the front or end, but no internal nodes.
                    // The read cannot be informative for any SnarlTraversal in this case.
                    consistencies[i] = false;
                    continue;
                }
                else{
                    // maps to neither front nor end
                    // The read could map internally, or not at all.
                    // Unless we know that the internal sequence is unique we can't guaratee that
                    // the mapping is consistent.
                    consistencies[i] = false;

                }

            }
            
        return consistencies;
}

SimpleTraversalSupportCalculator::~SimpleTraversalSupportCalculator(){

}

vector<Support> SimpleTraversalSupportCalculator::calculate_supports(const Snarl& site,
        const vector<SnarlTraversal>& traversals, const vector<Alignment*>& reads,
        const vector<vector<bool>>& consistencies) const{
        // Calculate the number of reads that support
        // the Traversal, and how they support it.    
        vector<Support> site_supports(traversals.size());

        for (int i = 0; i < reads.size(); i++){
            vector<bool> cons = consistencies[i];
            for (int t = 0; t < traversals.size(); t++){
                Support s;
                if (cons[t] == true && !(reads[i]->read_on_reverse_strand())){
                   s.set_forward(s.forward() + 1);
                }
                else if (cons[t] == true && reads[i]->read_on_reverse_strand()){
                    s.set_reverse(s.reverse() + 1);
                }
                else{
                    continue;
                }
            }
            
        }

        return site_supports;
}

    
double FixedGenotypePriorCalculator::calculate_log_prior(const Genotype& genotype) {
    // Are all the alleles the same?
    bool all_same = true;
    // What is the common allele number (-1 for unset)
    int allele_value = -1;
    for(size_t i = 0; i < genotype.allele_size(); i++) {
        // For each allele in the genotype
        if(allele_value == -1) {
            // On the first one, grab it. Everyone else will have to match.
            allele_value = genotype.allele(i);
        }
        
        if(allele_value != genotype.allele(i)) {
            // There are two distinct allele values in here
            all_same = false;
            break;
        }
    }
    
    // Return the appropriate prior depending on whether the alleles are all the
    // same (homozygous) or not (heterozygous).
    return all_same ? homozygous_prior_ln : heterozygous_prior_ln;
}


Support make_support(double forward, double reverse, double quality) {
    Support to_return;
    to_return.set_forward(forward);
    to_return.set_reverse(reverse);
    to_return.set_quality(quality);
    return to_return;
}

double total(const Support& support) {
    return support.forward() + support.reverse();
}

Support support_min(const Support& a, const Support& b) {
    Support to_return;
    to_return.set_forward(min(a.forward(), b.forward()));
    to_return.set_reverse(min(a.reverse(), b.reverse()));
    to_return.set_quality(min(a.quality(), b.quality()));
    return to_return;
}

Support support_max(const Support& a, const Support& b) {
    Support to_return;
    to_return.set_forward(max(a.forward(), b.forward()));
    to_return.set_reverse(max(a.reverse(), b.reverse()));
    to_return.set_quality(max(a.quality(), b.quality()));    
    return to_return;
}

Support flip(const Support& to_flip) {
    Support flipped = to_flip;
    flipped.set_forward(to_flip.reverse());
    flipped.set_reverse(to_flip.forward());
    return flipped;
}

Support operator+(const Support& one, const Support& other) {
    Support sum;
    sum.set_forward(one.forward() + other.forward());
    sum.set_reverse(one.reverse() + other.reverse());
    sum.set_left(one.left() + other.left());
    sum.set_right(one.right() + other.right());
    
    // log-scaled quality can just be added
    sum.set_quality(one.quality() + other.quality());
    
    return sum;
}

Support& operator+=(Support& one, const Support& other) {
    one.set_forward(one.forward() + other.forward());
    one.set_reverse(one.reverse() + other.reverse());
    one.set_left(one.left() + other.left());
    one.set_right(one.right() + other.right());
    
    // log-scaled quality can just be added
    one.set_quality(one.quality() + other.quality());
    
    return one;
}

bool operator< (const Support& a, const Support& b) {
    return total(a) < total(b);
}

bool operator> (const Support& a, const Support& b) {
    return total(a) > total(b);
}

ostream& operator<<(ostream& stream, const Support& support) {
    return stream << support.forward() << "," << support.reverse();
}

string to_vcf_genotype(const Genotype& gt) {
    // Emit parts into this stream
    stringstream stream;
    
    for (size_t i = 0; i < gt.allele_size(); i++) {
        // For each allele called as present in the genotype
        
        // Put it in the string
        stream << gt.allele(i);
        
        if (i + 1 != gt.allele_size()) {
            // Write a separator after all but the last one
            stream << (gt.is_phased() ? '|' : '/');
        }
    }
    
    return stream.str();
}


}
