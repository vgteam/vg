#include "aligner.hpp"

#include "hash_map.hpp"
#include "alignment.hpp"
#include "path.hpp"
#include "utility.hpp"
#include "statistics.hpp"
#include "banded_global_aligner.hpp"
#include "reverse_graph.hpp"
#include "null_masking_graph.hpp"
#include "dozeu_pinning_overlay.hpp"
#include "crash.hpp"
#include "algorithms/distance_to_tail.hpp"

//#define debug_print_score_matrices

namespace vg {

using namespace std;
using namespace vg::io;

GSSWAligner::~GSSWAligner() = default;

GSSWAligner::GSSWAligner(std::unique_ptr<MatrixAlignmentScorer> owned_scorer)
    : scorer(std::move(owned_scorer)),
      mapq_calc(std::make_unique<MappingQualityCalculator>(*scorer)),
      deletion_aligner(scorer->gap_open, scorer->gap_extension) {
}

gssw_graph* GSSWAligner::create_gssw_graph(const HandleGraph& g) const {
    
    vector<handle_t> topological_order = handlealgs::lazier_topological_order(&g);
    
    gssw_graph* graph = gssw_graph_create(g.get_node_count());
    unordered_map<int64_t, gssw_node*> nodes;
    
    // compute the topological order
    for (const handle_t& handle : topological_order) {
        auto cleaned_seq = nonATGCNtoN(g.get_sequence(handle));
        gssw_node* node = gssw_node_create(nullptr,       // TODO: the ID should be enough, don't need Node* too
                                           g.get_id(handle),
                                           cleaned_seq.c_str(),
                                           scorer->nt_table,
                                           scorer->score_matrix); // TODO: this arg isn't used, could edit
                                                          // in gssw
        nodes[node->id] = node;
        gssw_graph_add_node(graph, node);
    }
    
    g.for_each_edge([&](const edge_t& edge) {
        if(!g.get_is_reverse(edge.first) && !g.get_is_reverse(edge.second)) {
            // This is a normal end to start edge.
            gssw_nodes_add_edge(nodes[g.get_id(edge.first)], nodes[g.get_id(edge.second)]);
        }
        else if (g.get_is_reverse(edge.first) && g.get_is_reverse(edge.second)) {
            // This is a start to end edge, but isn't reversing and can be converted to a normal end to start edge.
            
            // Flip the start and end
            gssw_nodes_add_edge(nodes[g.get_id(edge.second)], nodes[g.get_id(edge.first)]);
        }
        else {
            // TODO: It's a reversing edge, which gssw doesn't support yet. What
            // we should really do is do a topological sort to break cycles, and
            // then flip everything at the lower-rank end of this edge around,
            // so we don't have to deal with its reversing-ness. But for now we
            // just die so we don't get nonsense into gssw.
#pragma omp critical
            {
                // We need the critical section so we don't throw uncaught
                // exceptions in multiple threads at once, leading to C++ trying
                // to run termiante in parallel. This doesn't make it safe, just
                // slightly safer.
                cerr << "Can't gssw over reversing edge " << g.get_id(edge.first) << (g.get_is_reverse(edge.first) ? "-" : "+") << " -> " << g.get_id(edge.second) << (g.get_is_reverse(edge.second) ? "-" : "+") << endl;
                // TODO: there's no safe way to kill the program without a way
                // to signal the master to do it, via a shared variable in the
                // clause that made us parallel.
            }
            exit(1);
        }
        return true;
    });
    
    return graph;
    
}

unordered_set<vg::id_t> GSSWAligner::identify_pinning_points(const HandleGraph& graph) const {
    
    unordered_set<vg::id_t> return_val;
    
    // start at the sink nodes
    vector<handle_t> sinks = handlealgs::tail_nodes(&graph);
    
    // walk backwards to find non-empty nodes if necessary
    for (const handle_t& handle : sinks) {
        vector<handle_t> stack(1, handle);
        while (!stack.empty()) {
            handle_t here =  stack.back();
            stack.pop_back();
            
            if (graph.get_length(here) > 0) {
                return_val.insert(graph.get_id(here));
            }
            else {
                graph.follow_edges(here, true, [&](const handle_t& prev) {
                    // TODO: technically this won't filter out all redundant walks, but it should
                    // handle all cases we're practically interested in and it doesn't require a
                    // second set object
                    if (!return_val.count(graph.get_id(prev))) {
                        stack.push_back(prev);
                    }
                });
            }
        }
    }
    
    return return_val;
}

void GSSWAligner::gssw_mapping_to_alignment(gssw_graph* graph,
                                            gssw_graph_mapping* gm,
                                            Alignment& alignment,
                                            bool pinned,
                                            bool pin_left) const {
    alignment.clear_path();
    alignment.set_score(gm->score);
    alignment.set_query_position(0);
    Path* path = alignment.mutable_path();
    //alignment.set_cigar(graph_cigar(gm));
    
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* ncs = gc->elements;
    //cerr << "gm->position " << gm->position << endl;
    string& to_seq = *alignment.mutable_sequence();
    //cerr << "-------------" << endl;
    
#ifdef debug_print_score_matrices
    gssw_graph_print_score_matrices(graph, to_seq.c_str(), to_seq.size(), stderr);
#endif
    
    int to_pos = 0;
    int from_pos = gm->position;
    
    for (int i = 0; i < gc->length; ++i) {
        // check that the current alignment has a non-zero length
        gssw_cigar* c = ncs[i].cigar;
        int l = c->length;
        if (l == 0) continue;
        gssw_cigar_element* e = c->elements;
        
        gssw_node* node = ncs[i].node;
        Mapping* mapping = path->add_mapping();
        
        if (i > 0) {
            // reset for each node after the first
            from_pos = 0;
        }
        
        mapping->mutable_position()->set_node_id(node->id);
        mapping->mutable_position()->set_offset(from_pos);
        mapping->set_rank(path->mapping_size());
        
        //cerr << node->id << ":" << endl;
        
        for (int j=0; j < l; ++j, ++e) {
            int32_t length = e->length;
            //cerr << e->length << e->type << endl;
            
            Edit* edit;
            switch (e->type) {
                case 'M':
                case 'X':
                case 'N': {
                    //cerr << "j = " << j << ", type = " << e->type << endl;
                    // do the sequences match?
                    // emit a stream of "SNPs" and matches
                    int h = from_pos;
                    int last_start = from_pos;
                    int k = to_pos;
                    for ( ; h < from_pos + length; ++h, ++k) {
                        //cerr << h << ":" << k << " " << node->seq[h] << " " << to_seq[k] << endl;
                        if (node->seq[h] != to_seq[k]) {
                            // emit the last "match" region
                            if (h - last_start > 0) {
                                edit = mapping->add_edit();
                                edit->set_from_length(h-last_start);
                                edit->set_to_length(h-last_start);
                            }
                            // set up the SNP
                            edit = mapping->add_edit();
                            edit->set_from_length(1);
                            edit->set_to_length(1);
                            edit->set_sequence(to_seq.substr(k,1));
                            last_start = h+1;
                        }
                    }
                    // handles the match at the end or the case of no SNP
                    if (h - last_start > 0) {
                        edit = mapping->add_edit();
                        edit->set_from_length(h-last_start);
                        edit->set_to_length(h-last_start);
                    }
                    to_pos += length;
                    from_pos += length;
                } break;
                case 'D':
                    edit = mapping->add_edit();
                    edit->set_from_length(length);
                    edit->set_to_length(0);
                    from_pos += length;
                    break;
                case 'I':
                    edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(length);
                    edit->set_sequence(to_seq.substr(to_pos, length));
                    to_pos += length;
                    break;
                case 'S':
                    // note that soft clips and insertions are semantically equivalent
                    // and can only be differentiated by their position in the read
                    // with soft clips coming at the start or end
                    edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(length);
                    edit->set_sequence(to_seq.substr(to_pos, length));
                    to_pos += length;
                    break;
                default:
                    cerr << "error:[Aligner::gssw_mapping_to_alignment] "
                    << "unsupported cigar op type " << e->type << endl;
                    exit(1);
                    break;
                    
            }
        }
    }
    
    // compute and set identity
    alignment.set_identity(identity(alignment.path()));
}

void GSSWAligner::unreverse_graph(gssw_graph* graph) const {
    // this is only for getting correct reference-relative edits, so we can get away with only
    // reversing the sequences and not paying attention to the edges
    
    for (size_t i = 0; i < graph->size; i++) {
        gssw_node* node = graph->nodes[i];
        for (int j = 0, stop = node->len / 2; j < stop; j++) {
            std::swap(node->seq[j], node->seq[node->len - j - 1]);
        }
    }
}

void GSSWAligner::unreverse_graph_mapping(gssw_graph_mapping* gm) const {
    
    gssw_graph_cigar* graph_cigar = &(gm->cigar);
    gssw_node_cigar* node_cigars = graph_cigar->elements;
    
    // reverse the order of the node cigars
    int32_t num_switching_nodes = graph_cigar->length / 2;
    int32_t last_idx = graph_cigar->length - 1;
    for (int32_t i = 0; i < num_switching_nodes; i++) {
        std::swap(node_cigars[i], node_cigars[last_idx - i]);
    }
    
    // reverse the actual cigar string for each node cigar
    for (int32_t i = 0; i < graph_cigar->length; i++) {
        gssw_cigar* node_cigar = node_cigars[i].cigar;
        gssw_cigar_element* elements = node_cigar->elements;
        
        int32_t num_switching_elements = node_cigar->length / 2;
        last_idx = node_cigar->length - 1;
        for (int32_t j = 0; j < num_switching_elements; j++) {
            std::swap(elements[j], elements[last_idx - j]);
        }
    }
    
    // compute the position in the first node
    if (graph_cigar->length > 0) {
        gssw_cigar_element* first_node_elements = node_cigars[0].cigar->elements;
        int32_t num_first_node_elements = node_cigars[0].cigar->length;
        uint32_t num_ref_aligned = 0; // the number of characters on the node sequence that are aligned
        for (int32_t i = 0; i < num_first_node_elements; i++) {
            switch (first_node_elements[i].type) {
                case 'M':
                case 'X':
                case 'N':
                case 'D':
                    num_ref_aligned += first_node_elements[i].length;
                    break;
                    
            }
        }
        gm->position = node_cigars[0].node->len - num_ref_aligned - (graph_cigar->length == 1 ? gm->position : 0);
    }
    else {
        gm->position = 0;
    }
}

string GSSWAligner::graph_cigar(gssw_graph_mapping* gm) const {

    stringstream s;
    gssw_graph_cigar* gc = &gm->cigar;
    gssw_node_cigar* nc = gc->elements;
    int to_pos = 0;
    int from_pos = gm->position;
    //string& to_seq = *alignment.mutable_sequence();
    s << from_pos << '@';
    for (int i = 0; i < gc->length; ++i, ++nc) {
        if (i > 0) from_pos = 0; // reset for each node after the first
        Node* from_node = (Node*) nc->node->data;
        s << from_node->id() << ':';
        gssw_cigar* c = nc->cigar;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (int j=0; j < l; ++j, ++e) {
            s << e->length << e->type;
        }
        if (i + 1 < gc->length) {
            s << ",";
        }
    }
    return s.str();
}


Aligner::Aligner(const int8_t* score_matrix,
                 int8_t gap_open,
                 int8_t gap_extension,
                 int8_t full_length_bonus,
                 double gc_content)
    : GSSWAligner(std::make_unique<MatrixAlignmentScorer>(score_matrix, gap_open, gap_extension, full_length_bonus, gc_content))
{
    // make an XdropAligner for each thread
    int num_threads = get_thread_count();
    xdrops.reserve(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        xdrops.emplace_back(score_matrix, gap_open, gap_extension);
    }
}

void Aligner::align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, const HandleGraph& g,
                             bool pinned, bool pin_left,int32_t max_alt_alns, bool traceback_aln) const {
    // bench_start(bench);
    // check input integrity
    if (pin_left && !pinned) {
        cerr << "error:[Aligner] cannot choose pinned end in non-pinned alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (multi_alignments && !pinned) {
        cerr << "error:[Aligner] multiple traceback is not implemented in local alignment, only pinned and global" << endl;
        exit(EXIT_FAILURE);
    }
    if (!multi_alignments && max_alt_alns != 1) {
        cerr << "error:[Aligner] cannot specify maximum number of alignments in single alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (max_alt_alns <= 0) {
        cerr << "error:[Aligner] cannot do less than 1 alignment" << endl;
        exit(EXIT_FAILURE);
    }

    // alignment pinning algorithm is based on pinning in bottom right corner, if pinning in top
    // left we need to reverse all the sequences first and translate the alignment back later
    
    // make a place to reverse the graph and sequence if necessary
    ReverseGraph reversed_graph(&g, false);
    string reversed_sequence;

    // choose forward or reversed objects
    const HandleGraph* oriented_graph = &g;
    const string* align_sequence = &alignment.sequence();
    if (pin_left) {
        // choose the reversed graph
        oriented_graph = &reversed_graph;
        
        // make and assign the reversed sequence
        reversed_sequence.resize(align_sequence->size());
        reverse_copy(align_sequence->begin(), align_sequence->end(), reversed_sequence.begin());
        align_sequence = &reversed_sequence;
    }
    
    // to save compute, we won't make these unless we're doing pinning
    unordered_set<vg::id_t> pinning_ids;
    NullMaskingGraph* null_masked_graph = nullptr;
    const HandleGraph* align_graph = oriented_graph;
    if (pinned) {
        pinning_ids = identify_pinning_points(*oriented_graph);
        null_masked_graph = new NullMaskingGraph(oriented_graph);
        align_graph = null_masked_graph;
    }
    
    // convert into gssw graph
    gssw_graph* graph = create_gssw_graph(*align_graph);
    
    // perform dynamic programming
    gssw_graph_fill_pinned(graph, align_sequence->c_str(),
                           scorer->nt_table, scorer->score_matrix,
                           scorer->gap_open, scorer->gap_extension, scorer->full_length_bonus,
                           pinned ? 0 : scorer->full_length_bonus, 15, 2, traceback_aln);

    // traceback either from pinned position or optimal local alignment
    if (traceback_aln) {
        if (pinned) {
            // we can only run gssw's DP on non-empty graphs, but we may have masked the entire graph
            // if it consists of only empty nodes, so don't both with the DP in that case
            gssw_graph_mapping** gms = nullptr;
            if (align_graph->get_node_count() > 0) {
                gssw_node** pinning_nodes = (gssw_node**) malloc(pinning_ids.size() * sizeof(gssw_node*));
                crash_unless(pinning_nodes != nullptr);
                size_t j = 0;
                for (size_t i = 0; i < graph->size; i++) {
                    gssw_node* node = graph->nodes[i];
                    if (pinning_ids.count(node->id)) {
                        pinning_nodes[j] = node;
                        j++;
                    }
                }
                
                // trace back pinned alignment
                gms = gssw_graph_trace_back_pinned_multi (graph,
                                                          max_alt_alns,
                                                          true,
                                                          align_sequence->c_str(),
                                                          align_sequence->size(),
                                                          pinning_nodes,
                                                          pinning_ids.size(),
                                                          scorer->nt_table,
                                                          scorer->score_matrix,
                                                          scorer->gap_open,
                                                          scorer->gap_extension,
                                                          scorer->full_length_bonus,
                                                          0);
                
                free(pinning_nodes);
            }
            
            // did we both 1) do DP (i.e. the graph is non-empty), and 2) find a traceback with positive score?
            if (gms ? gms[0]->score > 0 : false) {
                
                if (pin_left) {
                    // translate nodes and mappings into original sequence so that the cigars come out right
                    unreverse_graph(graph);
                    for (int32_t i = 0; i < max_alt_alns; i++) {
                        unreverse_graph_mapping(gms[i]);
                    }
                }
                
                // have a mapping, can just convert normally
                gssw_mapping_to_alignment(graph, gms[0], alignment, pinned, pin_left);
                
                if (multi_alignments) {
                    // determine how many non-null alignments were returned
                    int32_t num_non_null = max_alt_alns;
                    for (int32_t i = 1; i < max_alt_alns; i++) {
                        if (gms[i]->score <= 0) {
                            num_non_null = i;
                            break;
                        }
                    }
                    
                    // reserve to avoid illegal access errors that occur when the vector reallocates
                    multi_alignments->reserve(num_non_null);
                    
                    // copy the primary alignment
                    multi_alignments->emplace_back(alignment);
                    
                    // convert the alternate alignments and store them at the back of the vector (this will not
                    // execute if we are doing single alignment)
                    for (int32_t i = 1; i < num_non_null; i++) {
                        // make new alignment object
                        multi_alignments->emplace_back();
                        Alignment& next_alignment = multi_alignments->back();
                        
                        // copy over sequence information from the primary alignment
                        next_alignment.set_sequence(alignment.sequence());
                        next_alignment.set_quality(alignment.quality());
                        
                        // get path of the alternate alignment
                        gssw_mapping_to_alignment(graph, gms[i], next_alignment, pinned, pin_left);
                    }
                }
            }
            else if (g.get_node_count() > 0) {
                // we didn't get any alignments either because the graph was empty and we couldn't run
                // gssw DP or because they had score 0 and gssw didn't want to do traceback. however,
                // we can infer the location of softclips based on the pinning nodes, so we'll just make
                // those manually
                
                // find the sink nodes of the oriented graph, which may be empty
                auto pinning_points = handlealgs::tail_nodes(oriented_graph);
                // impose a consistent ordering for machine independent behavior
                sort(pinning_points.begin(), pinning_points.end(), [&](const handle_t& h1, const handle_t& h2) {
                    return oriented_graph->get_id(h1) < oriented_graph->get_id(h2);
                });
                
                for (size_t i = 0; i < max_alt_alns && i < pinning_points.size(); i++) {
                    // make a record in the multi alignments if we're using them
                    if (multi_alignments) {
                        multi_alignments->emplace_back();
                    }
                    // choose an alignment object to construct the path in
                    Alignment& softclip_alignment = i == 0 ? alignment : multi_alignments->back();
                    
                    handle_t& pinning_point = pinning_points[i];
                    
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    mapping->set_rank(1);
                    
                    // locate at the beginning or end of the node
                    Position* position = mapping->mutable_position();
                    position->set_node_id(oriented_graph->get_id(pinning_point));
                    position->set_offset(pin_left ? 0 : oriented_graph->get_length(pinning_point));
                    
                    // soft clip
                    Edit* edit = mapping->add_edit();
                    edit->set_to_length(alignment.sequence().length());
                    edit->set_sequence(alignment.sequence());
                    
                    // we want to also have the first alignment in the multi-alignment vector
                    if (i == 0 && multi_alignments) {
                        multi_alignments->back() = alignment;
                    }
                }
            }
            if (gms) {
                for (int32_t i = 0; i < max_alt_alns; i++) {
                    gssw_graph_mapping_destroy(gms[i]);
                }
                free(gms);
            }
        }
        else {
            // trace back local alignment
            gssw_graph_mapping* gm = gssw_graph_trace_back (graph,
                                                            align_sequence->c_str(),
                                                            align_sequence->size(),
                                                            scorer->nt_table,
                                                            scorer->score_matrix,
                                                            scorer->gap_open,
                                                            scorer->gap_extension,
                                                            scorer->full_length_bonus,
                                                            scorer->full_length_bonus);
        
            gssw_mapping_to_alignment(graph, gm, alignment, pinned, pin_left);
            gssw_graph_mapping_destroy(gm);
        }
    } else {
        // get the alignment position and score
        alignment.set_score(graph->max_node->alignment->score1);
        Mapping* m = alignment.mutable_path()->add_mapping();
        Position* p = m->mutable_position();
        p->set_node_id(graph->max_node->id);
        p->set_offset(graph->max_node->alignment->ref_end1); // mark end position; for de-duplication
    }
        
    // this might be null if we're not doing pinned alignment, but delete doesn't care
    delete null_masked_graph;
    
    gssw_graph_destroy(graph);
    // bench_end(bench);
}

void Aligner::align(Alignment& alignment, const HandleGraph& g, bool traceback_aln) const {
    
    align_internal(alignment, nullptr, g, false, false, 1, traceback_aln);
}

void Aligner::align(Alignment& alignment, const HandleGraph& g,
                    const std::vector<handle_t>& topological_order) const {

    // Create a gssw_graph and a mapping from handles to nodes.
    gssw_graph* graph = gssw_graph_create(topological_order.size());
    hash_map<handle_t, gssw_node*> nodes;
    nodes.reserve(topological_order.size());

    // Create the nodes. Use offsets in the topological order as node ids.
    for (size_t i = 0; i < topological_order.size(); i++) {
        handle_t handle = topological_order[i];
        auto cleaned_seq = nonATGCNtoN(g.get_sequence(handle));
        gssw_node* node = gssw_node_create(nullptr,
                                           i,
                                           cleaned_seq.c_str(),
                                           scorer->nt_table,
                                           scorer->score_matrix);
        nodes[handle] = node;
        gssw_graph_add_node(graph, node);
    }

    // Create the edges.
    for (const handle_t& from : topological_order) {
        gssw_node* from_node = nodes[from];
        g.follow_edges(from, false, [&](const handle_t& to) {
            auto iter = nodes.find(to);
            if (iter != nodes.end()) {
                gssw_nodes_add_edge(from_node, iter->second);
            }
        });
    }

    // Align the read to the subgraph.
    gssw_graph_fill_pinned(graph, alignment.sequence().c_str(),
                           scorer->nt_table, scorer->score_matrix,
                           scorer->gap_open, scorer->gap_extension, scorer->full_length_bonus, scorer->full_length_bonus,
                           15, 2, true);
    gssw_graph_mapping* gm = gssw_graph_trace_back(graph,
                                                   alignment.sequence().c_str(), alignment.sequence().length(),
                                                   scorer->nt_table, scorer->score_matrix,
                                                   scorer->gap_open, scorer->gap_extension, scorer->full_length_bonus, scorer->full_length_bonus);

    // Convert the mapping to Alignment.
    this->gssw_mapping_to_alignment(graph, gm, alignment, false, false);
    Path& path = *(alignment.mutable_path());
    for (size_t i = 0; i < path.mapping_size(); i++) {
        Position& pos = *(path.mutable_mapping(i)->mutable_position());
        handle_t handle = topological_order[pos.node_id()];
        pos.set_node_id(g.get_id(handle));
        pos.set_is_reverse(g.get_is_reverse(handle));
    }

    // Destroy the temporary objects.
    gssw_graph_mapping_destroy(gm);
    gssw_graph_destroy(graph);
}

void Aligner::align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left, bool xdrop,
                           uint16_t xdrop_max_gap_length) const {
    
    if (xdrop) {
        // XdropAligner manages its own stack, so it can never be threadsafe without be recreated
        // for every alignment, which meshes poorly with its stack implementation. We achieve
        // thread-safety by having one per thread, which makes this method const-ish.
        XdropAligner& xdrop = const_cast<XdropAligner&>(xdrops[omp_get_thread_num()]);
        
        // dozeu declines to produce an alignment when the gap is set to 0
        xdrop_max_gap_length = max<uint16_t>(xdrop_max_gap_length, 1);
        
        // wrap the graph so that empty pinning points are handled correctly
        DozeuPinningOverlay overlay(&g, !pin_left);
        
        if (overlay.get_node_count() == 0 && g.get_node_count() != 0) {
            // the only nodes in the graph are empty nodes for pinning, which got masked.
            // we can still infer a pinned alignment based purely on the pinning point but
            // dozeu won't handle this correctly
            g.for_each_handle([&](const handle_t& handle) {
                bool can_pin = g.follow_edges(handle, pin_left, [&](const handle_t& next) {return false;});
                if (can_pin) {
                    // manually make the softclip
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    Position* pos = mapping->mutable_position();
                    pos->set_node_id(g.get_id(handle));
                    pos->set_is_reverse(false);
                    pos->set_offset(pin_left ? 0 : g.get_length(handle));
                    
                    mapping->set_rank(1);
                    
                    Edit* edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(alignment.sequence().size());
                    edit->set_sequence(alignment.sequence());
                    alignment.set_score(0);
                    return false;
                }
                return true;
            });
        }
        else {
            // do the alignment
            xdrop.align_pinned(alignment, overlay, pin_left, scorer->full_length_bonus, xdrop_max_gap_length);
            
            if (overlay.performed_duplications()) {
                // the overlay is not a strict subset of the underlying graph, so we may
                // need to translate some node IDs
                translate_oriented_node_ids(*alignment.mutable_path(), [&](id_t node_id) {
                    handle_t under = overlay.get_underlying_handle(overlay.get_handle(node_id));
                    return make_pair(g.get_id(under), g.get_is_reverse(under));
                });
            }
        }
    }
    else {
        align_internal(alignment, nullptr, g, true, pin_left, 1, true);
    }
}

void Aligner::align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                 bool pin_left, int32_t max_alt_alns) const {
    
    if (alt_alignments.size() != 0) {
        cerr << "error:[Aligner::align_pinned_multi] output vector must be empty for pinned multi-aligning" << endl;
        exit(EXIT_FAILURE);
    }
    
    align_internal(alignment, &alt_alignments, g, true, pin_left, max_alt_alns, true);
}

void Aligner::align_global_banded(Alignment& alignment, const HandleGraph& g,
                                  int32_t band_padding, bool permissive_banding,
                                  uint64_t max_cells) const {
    
    if (alignment.sequence().empty()) {
        // we can save time by using a specialized deletion aligner for empty strings
        deletion_aligner.align(alignment, g);
        return;
    }
    
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * scorer->match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(scorer->mismatch, scorer->gap_open), scorer->gap_extension);
    
    // TODO: put this all into another template somehow?
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               band_padding,
                                               permissive_banding,
                                               false,
                                               max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                false,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    }
}

void Aligner::align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                        int32_t max_alt_alns, int32_t band_padding, bool permissive_banding,
                                        uint64_t max_cells) const {
                              
    if (alignment.sequence().empty()) {
        // we can save time by using a specialized deletion aligner for empty strings
        deletion_aligner.align_multi(alignment, alt_alignments, g, max_alt_alns);
        return;
    }
    
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * scorer->match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(scorer->mismatch, scorer->gap_open), scorer->gap_extension);
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               alt_alignments,
                                               max_alt_alns,
                                               band_padding,
                                               permissive_banding,
                                               false,
                                               max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                false,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    }
}

void Aligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<MaximalExactMatch>& mems,
                          bool reverse_complemented, uint16_t max_gap_length) const
{
    align_xdrop(alignment, g, handlealgs::lazier_topological_order(&g), mems, reverse_complemented,
                max_gap_length);
}

void Aligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order,
                          const vector<MaximalExactMatch>& mems, bool reverse_complemented, uint16_t max_gap_length) const
{
    // XdropAligner manages its own stack, so it can never be threadsafe without be recreated
    // for every alignment, which meshes poorly with its stack implementation. We achieve
    // thread-safety by having one per thread, which makes this method const-ish.
    XdropAligner& xdrop = const_cast<XdropAligner&>(xdrops[omp_get_thread_num()]);
    xdrop.align(alignment, g, order, mems, reverse_complemented, scorer->full_length_bonus, max_gap_length);
    if (!alignment.has_path() && mems.empty()) {
        // dozeu couldn't find an alignment, probably because it's seeding heuristic failed
        // we'll just fall back on GSSW
        // TODO: This is a bit inconsistent. GSSW gives a full-length bonus at both ends, while
        // dozeu only gives it once.
        align(alignment, g, order);
    }
}



QualAdjAligner::QualAdjAligner(const int8_t* score_matrix,
                               int8_t gap_open,
                               int8_t gap_extension,
                               int8_t full_length_bonus,
                               double gc_content)
    : GSSWAligner(std::make_unique<QualAdjAlignmentScorer>(score_matrix, gap_open, gap_extension,
                                                           full_length_bonus, gc_content))
{
    // make a QualAdjXdropAligner for each thread
    int num_threads = get_thread_count();
    xdrops.reserve(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        xdrops.emplace_back(score_matrix, scorer->score_matrix, gap_open, gap_extension);
    }
}


void QualAdjAligner::align_internal(Alignment& alignment, vector<Alignment>* multi_alignments, const HandleGraph& g,
                                    bool pinned, bool pin_left, int32_t max_alt_alns, bool traceback_aln) const {
    
    // check input integrity
    if (pin_left && !pinned) {
        cerr << "error:[Aligner] cannot choose pinned end in non-pinned alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (multi_alignments && !pinned) {
        cerr << "error:[Aligner] multiple traceback is not implemented in local alignment, only pinned and global" << endl;
        exit(EXIT_FAILURE);
    }
    if (!multi_alignments && max_alt_alns != 1) {
        cerr << "error:[Aligner] cannot specify maximum number of alignments in single alignment" << endl;
        exit(EXIT_FAILURE);
    }
    if (max_alt_alns <= 0) {
        cerr << "error:[Aligner] cannot do less than 1 alignment" << endl;
        exit(EXIT_FAILURE);
    }
    
    // alignment pinning algorithm is based on pinning in bottom right corner, if pinning in top
    // left we need to reverse all the sequences first and translate the alignment back later
    
    // make a place to reverse the graph and sequence if necessary
    ReverseGraph reversed_graph(&g, false);
    string reversed_sequence;
    string reversed_quality;
    
    // choose forward or reversed objects
    const HandleGraph* oriented_graph = &g;
    const string* align_sequence = &alignment.sequence();
    const string* align_quality = &alignment.quality();
    if (pin_left) {
        // choose the reversed graph
        oriented_graph = &reversed_graph;
        
        // make and assign the reversed sequence
        reversed_sequence.resize(align_sequence->size());
        reverse_copy(align_sequence->begin(), align_sequence->end(), reversed_sequence.begin());
        align_sequence = &reversed_sequence;
        
        // make and assign the reversed quality
        reversed_quality.resize(align_quality->size());
        reverse_copy(align_quality->begin(), align_quality->end(), reversed_quality.begin());
        align_quality = &reversed_quality;
    }
    
    if (align_quality->size() != align_sequence->size()) {
        cerr << "error:[QualAdjAligner] Read " << alignment.name() << " has sequence and quality strings with different lengths. Cannot perform base quality adjusted alignment. Consider toggling off base quality adjusted alignment at the command line." << endl;
        exit(EXIT_FAILURE);
    }
    
    // to save compute, we won't make these unless we're doing pinning
    unordered_set<vg::id_t> pinning_ids;
    NullMaskingGraph* null_masked_graph = nullptr;
    const HandleGraph* align_graph = oriented_graph;
    if (pinned) {
        pinning_ids = identify_pinning_points(*oriented_graph);
        null_masked_graph = new NullMaskingGraph(oriented_graph);
        align_graph = null_masked_graph;
    }
    
    // convert into gssw graph
    gssw_graph* graph = create_gssw_graph(*align_graph);
    
    int8_t front_full_length_bonus = qa_scorer()->qual_adj_full_length_bonuses[align_quality->front()];
    int8_t back_full_length_bonus = qa_scorer()->qual_adj_full_length_bonuses[align_quality->back()];
    
    // perform dynamic programming
    // offer a full length bonus on each end, or only on the left if the right end is pinned.
    gssw_graph_fill_pinned_qual_adj(graph, align_sequence->c_str(), align_quality->c_str(),
                                    scorer->nt_table, scorer->score_matrix,
                                    scorer->gap_open, scorer->gap_extension,
                                    front_full_length_bonus,
                                    pinned ? 0 : back_full_length_bonus,
                                    15, 2, traceback_aln);
    
    // traceback either from pinned position or optimal local alignment
    if (traceback_aln) {
        if (pinned) {
            gssw_graph_mapping** gms = nullptr;
            if (align_graph->get_node_count() > 0) {
                
                gssw_node** pinning_nodes = (gssw_node**) malloc(pinning_ids.size() * sizeof(gssw_node*));
                crash_unless(pinning_nodes != nullptr);
                size_t j = 0;
                for (size_t i = 0; i < graph->size; i++) {
                    gssw_node* node = graph->nodes[i];
                    if (pinning_ids.count(node->id)) {
                        pinning_nodes[j] = node;
                        j++;
                    }
                }
                
                // trace back pinned alignment
                gms = gssw_graph_trace_back_pinned_qual_adj_multi (graph,
                                                                   max_alt_alns,
                                                                   true,
                                                                   align_sequence->c_str(),
                                                                   align_quality->c_str(),
                                                                   align_sequence->size(),
                                                                   pinning_nodes,
                                                                   pinning_ids.size(),
                                                                   scorer->nt_table,
                                                                   scorer->score_matrix,
                                                                   scorer->gap_open,
                                                                   scorer->gap_extension,
                                                                   front_full_length_bonus,
                                                                   0);
                
                free(pinning_nodes);
            }
            
            // did we both 1) do DP (i.e. the graph is non-empty), and 2) find a traceback with positive score?
            if (gms && gms[0]->score > 0) {
                
                if (pin_left) {
                    // translate graph and mappings into original node space
                    unreverse_graph(graph);
                    for (int32_t i = 0; i < max_alt_alns; i++) {
                        unreverse_graph_mapping(gms[i]);
                    }
                }
                
                // have a mapping, can just convert normally
                gssw_mapping_to_alignment(graph, gms[0], alignment, pinned, pin_left);
                
                if (multi_alignments) {
                    // determine how many non-null alignments were returned
                    int32_t num_non_null = max_alt_alns;
                    for (int32_t i = 1; i < max_alt_alns; i++) {
                        if (gms[i]->score <= 0) {
                            num_non_null = i;
                            break;
                        }
                    }
                    
                    // reserve to avoid illegal access errors that occur when the vector reallocates
                    multi_alignments->reserve(num_non_null);
                    
                    // copy the primary alignment
                    multi_alignments->emplace_back(alignment);
                    
                    // convert the alternate alignments and store them at the back of the vector (this will not
                    // execute if we are doing single alignment)
                    for (int32_t i = 1; i < num_non_null; i++) {
                        // make new alignment object
                        multi_alignments->emplace_back();
                        Alignment& next_alignment = multi_alignments->back();
                        
                        // copy over sequence information from the primary alignment
                        next_alignment.set_sequence(alignment.sequence());
                        next_alignment.set_quality(alignment.quality());
                        
                        // get path of the alternate alignment
                        gssw_mapping_to_alignment(graph, gms[i], next_alignment, pinned, pin_left);
                    }
                }
            }
            else if (g.get_node_count() > 0) {
                /// we didn't get any alignments either because the graph was empty and we couldn't run
                // gssw DP or because they had score 0 and gssw didn't want to do traceback. however,
                // we can infer the location of softclips based on the pinning nodes, so we'll just make
                // those manually
                
                // find the sink nodes of the oriented graph, which may be empty
                auto pinning_points = handlealgs::tail_nodes(oriented_graph);
                // impose a consistent ordering for machine independent behavior
                sort(pinning_points.begin(), pinning_points.end(), [&](const handle_t& h1, const handle_t& h2) {
                    return oriented_graph->get_id(h1) < oriented_graph->get_id(h2);
                });
                
                for (size_t i = 0; i < max_alt_alns && i < pinning_points.size(); i++) {
                    // make a record in the multi alignments if we're using them
                    if (multi_alignments) {
                        multi_alignments->emplace_back();
                    }
                    // choose an alignment object to construct the path in
                    Alignment& softclip_alignment = i == 0 ? alignment : multi_alignments->back();
                    
                    handle_t& pinning_point = pinning_points[i];
                    
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    mapping->set_rank(1);
                    
                    // locate at the beginning or end of the node
                    Position* position = mapping->mutable_position();
                    position->set_node_id(oriented_graph->get_id(pinning_point));
                    position->set_offset(pin_left ? 0 : oriented_graph->get_length(pinning_point));
                    
                    // soft clip
                    Edit* edit = mapping->add_edit();
                    edit->set_to_length(alignment.sequence().length());
                    edit->set_sequence(alignment.sequence());
                    
                    // we want to also have the first alignment in the multi-alignment vector
                    if (i == 0 && multi_alignments) {
                        multi_alignments->back() = alignment;
                    }
                }
            }
            
            if (gms) {
                for (int32_t i = 0; i < max_alt_alns; i++) {
                    gssw_graph_mapping_destroy(gms[i]);
                }
                free(gms);
            }
        }
        else {
            // trace back local alignment
            gssw_graph_mapping* gm = gssw_graph_trace_back_qual_adj (graph,
                                                                     align_sequence->c_str(),
                                                                     align_quality->c_str(),
                                                                     align_sequence->size(),
                                                                     scorer->nt_table,
                                                                     scorer->score_matrix,
                                                                     scorer->gap_open,
                                                                     scorer->gap_extension,
                                                                     front_full_length_bonus,
                                                                     back_full_length_bonus);
        
            gssw_mapping_to_alignment(graph, gm, alignment, pinned, pin_left);
            gssw_graph_mapping_destroy(gm);
        }
    } else {
        // get the alignment position and score
        alignment.set_score(graph->max_node->alignment->score1);
        Mapping* m = alignment.mutable_path()->add_mapping();
        Position* p = m->mutable_position();
        p->set_node_id(graph->max_node->id);
        p->set_offset(graph->max_node->alignment->ref_end1); // mark end position; for de-duplication
    }
        
    // this might be null if we're not doing pinned alignment, but delete doesn't care
    delete null_masked_graph;
    
    gssw_graph_destroy(graph);
    
}

void QualAdjAligner::align(Alignment& alignment, const HandleGraph& g, bool traceback_aln) const {
    
    align_internal(alignment, nullptr, g, false, false, 1, traceback_aln);
}

void QualAdjAligner::align_pinned(Alignment& alignment, const HandleGraph& g, bool pin_left, bool xdrop,
                                  uint16_t xdrop_max_gap_length) const {
    if (xdrop) {
        // QualAdjXdropAligner manages its own stack, so it can never be threadsafe without be recreated
        // for every alignment, which meshes poorly with its stack implementation. We achieve
        // thread-safety by having one per thread, which makes this method const-ish.
        QualAdjXdropAligner& xdrop = const_cast<QualAdjXdropAligner&>(xdrops[omp_get_thread_num()]);
        
        // wrap the graph so that empty pinning points are handled correctly
        DozeuPinningOverlay overlay(&g, !pin_left);
        if (overlay.get_node_count() == 0 && g.get_node_count() != 0) {
            // the only nodes in the graph are empty nodes for pinning, which got masked.
            // we can still infer a pinned alignment based purely on the pinning point but
            // dozeu won't handle this correctly
            g.for_each_handle([&](const handle_t& handle) {
                bool can_pin = g.follow_edges(handle, pin_left, [&](const handle_t& next) {return false;});
                if (can_pin) {
                    // manually make the softclip
                    Mapping* mapping = alignment.mutable_path()->add_mapping();
                    Position* pos = mapping->mutable_position();
                    pos->set_node_id(g.get_id(handle));
                    pos->set_is_reverse(false);
                    pos->set_offset(pin_left ? 0 : g.get_length(handle));
                    
                    mapping->set_rank(1);
                    
                    Edit* edit = mapping->add_edit();
                    edit->set_from_length(0);
                    edit->set_to_length(alignment.sequence().size());
                    edit->set_sequence(alignment.sequence());
                    alignment.set_score(0);
                    return false;
                }
                return true;
            });
        }
        else {
            
            // dozeu declines to produce an alignment when the gap is set to 0
            xdrop_max_gap_length = max<uint16_t>(xdrop_max_gap_length, 1);
            
            // get the quality adjusted bonus
            int8_t bonus = qa_scorer()->qual_adj_full_length_bonuses[pin_left ? alignment.quality().back() : alignment.quality().front()];
            
            xdrop.align_pinned(alignment, overlay, pin_left, bonus, xdrop_max_gap_length);
            
            if (overlay.performed_duplications()) {
                // the overlay is not a strict subset of the underlying graph, so we may
                // need to translate some node IDs
                translate_oriented_node_ids(*alignment.mutable_path(), [&](id_t node_id) {
                    handle_t under = overlay.get_underlying_handle(overlay.get_handle(node_id));
                    return make_pair(g.get_id(under), g.get_is_reverse(under));
                });
            }
        }
    }
    else {
        align_internal(alignment, nullptr, g, true, pin_left, 1, true);
    }
}

void QualAdjAligner::align_pinned_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                        bool pin_left, int32_t max_alt_alns) const {
    align_internal(alignment, &alt_alignments, g, true, pin_left, max_alt_alns, true);
}

void QualAdjAligner::align_global_banded(Alignment& alignment, const HandleGraph& g,
                                         int32_t band_padding, bool permissive_banding,
                                         uint64_t max_cells) const {
    
    if (alignment.sequence().empty()) {
        // we can save time by using a specialized deletion aligner for empty strings
        deletion_aligner.align(alignment, g);
        return;
    }
    
    int64_t best_score = alignment.sequence().size() * scorer->match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(scorer->mismatch, scorer->gap_open), scorer->gap_extension);
    
    // TODO: put this all into another template somehow?
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               band_padding,
                                               permissive_banding,
                                               true,
                                               max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                true,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                true,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                band_padding,
                                                permissive_banding,
                                                true,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    }
}

void QualAdjAligner::align_global_banded_multi(Alignment& alignment, vector<Alignment>& alt_alignments, const HandleGraph& g,
                                               int32_t max_alt_alns, int32_t band_padding, bool permissive_banding,
                                               uint64_t max_cells) const {
    
    if (alignment.sequence().empty()) {
        // we can save time by using a specialized deletion aligner for empty strings
        deletion_aligner.align_multi(alignment, alt_alignments, g, max_alt_alns);
        return;
    }
    
    // We need to figure out what size ints we need to use.
    // Get upper and lower bounds on the scores. TODO: if these overflow int64 we're out of luck
    int64_t best_score = alignment.sequence().size() * scorer->match;
    size_t total_bases = 0;
    g.for_each_handle([&](const handle_t& handle) {
        total_bases += g.get_length(handle);
    });
    int64_t worst_score = (alignment.sequence().size() + total_bases) * -max(max(scorer->mismatch, scorer->gap_open), scorer->gap_extension);
    
    if (best_score <= numeric_limits<int8_t>::max() && worst_score >= numeric_limits<int8_t>::min()) {
        // We'll fit in int8
        BandedGlobalAligner<int8_t> band_graph(alignment,
                                               g,
                                               alt_alignments,
                                               max_alt_alns,
                                               band_padding,
                                               permissive_banding,
                                               true,
                                               max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int16_t>::max() && worst_score >= numeric_limits<int16_t>::min()) {
        // We'll fit in int16
        BandedGlobalAligner<int16_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                true,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else if (best_score <= numeric_limits<int32_t>::max() && worst_score >= numeric_limits<int32_t>::min()) {
        // We'll fit in int32
        BandedGlobalAligner<int32_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                true,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    } else {
        // Fall back to int64
        BandedGlobalAligner<int64_t> band_graph(alignment,
                                                g,
                                                alt_alignments,
                                                max_alt_alns,
                                                band_padding,
                                                permissive_banding,
                                                true,
                                                max_cells);
        
        band_graph.align(scorer->score_matrix, scorer->nt_table, scorer->gap_open, scorer->gap_extension);
    }
}

void QualAdjAligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<MaximalExactMatch>& mems,
                                 bool reverse_complemented, uint16_t max_gap_length) const
{
    align_xdrop(alignment, g, handlealgs::lazier_topological_order(&g), mems, reverse_complemented, max_gap_length);
}

void QualAdjAligner::align_xdrop(Alignment& alignment, const HandleGraph& g, const vector<handle_t>& order,
                                 const vector<MaximalExactMatch>& mems, bool reverse_complemented,
                                 uint16_t max_gap_length) const
{
    // QualAdjXdropAligner manages its own stack, so it can never be threadsafe without being recreated
    // for every alignment, which meshes poorly with its stack implementation. We achieve
    // thread-safety by having one per thread, which makes this method const-ish.
    QualAdjXdropAligner& xdrop = const_cast<QualAdjXdropAligner&>(xdrops[omp_get_thread_num()]);
    
    // get the quality adjusted bonus
    int8_t bonus = qa_scorer()->qual_adj_full_length_bonuses[reverse_complemented ? alignment.quality().front() : alignment.quality().back()];
    
    xdrop.align(alignment, g, order, mems, reverse_complemented, bonus, max_gap_length);
    if (!alignment.has_path() && mems.empty()) {
        // dozeu couldn't find an alignment, probably because it's seeding heuristic failed
        // we'll just fall back on GSSW
        // TODO: This is a bit inconsistent. GSSW gives a full-length bonus at both ends, while
        // dozeu only gives it once.
        align(alignment, g, true);
    }
}


AlignerClient::AlignerClient(double gc_content_estimate) : gc_content_estimate(gc_content_estimate) {
   
    // Adopt the default scoring parameters and make the aligners
    AlignerClient::set_alignment_scores(default_score_matrix,
                                        default_gap_open, default_gap_extension,
                                        default_full_length_bonus);
    // Note that we can't make a virtual call here; if you override the method
    // you also need to do extra work in your constructor!
}

const GSSWAligner* AlignerClient::get_aligner(bool have_qualities) const {
    return (have_qualities && adjust_alignments_for_base_quality) ?
        (GSSWAligner*) get_qual_adj_aligner() :
        (GSSWAligner*) get_regular_aligner();
}

const QualAdjAligner* AlignerClient::get_qual_adj_aligner() const {
    assert(qual_adj_aligner.get() != nullptr);
    return qual_adj_aligner.get();
}

const Aligner* AlignerClient::get_regular_aligner() const {
    assert(regular_aligner.get() != nullptr);
    return regular_aligner.get();
}

int8_t* AlignerClient::parse_matrix(istream& matrix_stream) {
    int8_t* matrix = (int8_t*) malloc(16 * sizeof(int8_t));
    crash_unless(matrix != nullptr);
    for (size_t i = 0; i < 16; i++) {
        if (!matrix_stream.good()) {
            std::cerr << "error: vg Aligner::parse_matrix requires a 4x4 whitespace separated integer matrix\n";
            throw "";
        }
        int score;
        matrix_stream >> score;
        if (score > 127 || score < -127) {
            std::cerr << "error: vg Aligner::parse_matrix requires values in the range [-127,127]\n";
            throw "";
        }
        matrix[i] = score;
    }
    return matrix;
}

void AlignerClient::set_alignment_scores(int8_t match, int8_t mismatch, int8_t gap_open, int8_t gap_extend, 
                                         int8_t full_length_bonus) {
    
    int8_t* matrix = (int8_t*) malloc(sizeof(int8_t) * 16);
    crash_unless(matrix != nullptr);
    for (size_t i = 0; i < 16; ++i) {
        if (i % 5 == 0) {
            // Fill the diagonal with matches
            matrix[i] = match;
        }
        else {
            // And everywhere else with mismatches
            matrix[i] = -mismatch;
        }
    }

    this->set_alignment_scores(matrix, gap_open, gap_extend, full_length_bonus);

    free(matrix);
}


void AlignerClient::set_alignment_scores(const int8_t* score_matrix, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus) {
    
    qual_adj_aligner = unique_ptr<QualAdjAligner>(new QualAdjAligner(score_matrix, gap_open, gap_extend,
                                                                     full_length_bonus, gc_content_estimate));
    regular_aligner = unique_ptr<Aligner>(new Aligner(score_matrix, gap_open, gap_extend,
                                                      full_length_bonus, gc_content_estimate));
    
}

void AlignerClient::set_alignment_scores(std::istream& matrix_stream, int8_t gap_open, int8_t gap_extend, int8_t full_length_bonus) {
    int8_t* score_matrix = parse_matrix(matrix_stream);
    this->set_alignment_scores(score_matrix, gap_open, gap_extend, full_length_bonus);
    free(score_matrix);
}

}
