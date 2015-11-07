#include <cstdlib>
#include <stdexcept>
#include "json2pb.h"
#include "caller.hpp"
#include "stream.hpp"

using namespace std;

namespace vg {

const double Caller::Log_zero = (double)-1e100;

// these values pretty arbitrary at this point
const double Caller::Default_het_prior = 0.001; // from MAQ
const int Caller::Default_min_depth = 20;
const int Caller::Default_max_depth = 5000;
const int Caller::Default_min_support = 20;
const double Caller::Default_min_frac = 0.75;
const double Caller::Default_min_likelihood = 1e-50;
const char Caller::Default_default_quality = 30;

Caller::Caller(VG* graph,
               double het_prior,
               int min_depth,
               int max_depth,
               int min_support,
               double min_frac,
               double min_likelihood, 
               bool leave_uncalled,
               int default_quality):
    _graph(graph),
    _het_log_prior(safe_log(het_prior)),
    _hom_log_prior(safe_log(.5 * (1. - het_prior))),
    _min_depth(min_depth),
    _max_depth(max_depth),
    _min_support(min_support),
    _min_frac(min_frac),
    _min_log_likelihood(safe_log(min_likelihood)),
    _leave_uncalled(leave_uncalled),
    _default_quality(default_quality) {
    _max_id = _graph->max_node_id();
}

// delete contents of table
Caller::~Caller() {
    clear();
}

void Caller::clear() {
    _node_calls.clear();
    _call_graph = VG();
    _start_node_map.clear();
    _end_node_map.clear();
    _visited_nodes.clear();
}

void Caller::write_call_graph(ostream& out, bool json) {
    if (json) {
        out << pb2json(_call_graph.graph);
    } else {
        _call_graph.serialize_to_ostream(out);
    }
}

void Caller::call_node_pileup(const NodePileup& pileup) {

    _node = _graph->get_node(pileup.node_id());
    assert(_node != NULL);
    
    _node_calls.clear();
    char def_char = _leave_uncalled ? '.' : '-';
    _node_calls.assign(_node->sequence().length(), Genotype(def_char, def_char));

    // todo: parallelize this loop
    // process each base in pileup individually
    for (int i = 0; i < pileup.base_pileup_size(); ++i) {
        if (pileup.base_pileup(i).num_bases() >= _min_depth &&
            pileup.base_pileup(i).num_bases() <= _max_depth) {
            call_base_pileup(pileup, i);
        }
    }

    // add nodes and edges created when making calls to the output graph
    // (_side_map gets updated)
    create_node_calls(pileup);

    _visited_nodes.insert(_node->id());
}

void Caller::update_call_graph() {
    // if we're leaving uncalled nodes, add'em:
    if (_leave_uncalled) {
        function<void(Node*)> add_node = [&](Node* node) {
            if (_visited_nodes.find(node->id()) == _visited_nodes.end()) {
                _call_graph.create_node(node->sequence(), node->id());
                _start_node_map[node->id()] = NodePair(node->id(), -1);
                _end_node_map[node->id()] = NodePair(node->id(), -1);
            }
        };
        _graph->for_each_node(add_node);
    }

    // map every edge in the original graph to equivalent sides
    // in the call graph. if both sides exist, make an edge in the call graph
    function<void(Edge*)> map_edge = [&](Edge* edge) {
        const NodeMap& from_map = edge->from_start() ? _start_node_map : _end_node_map;
        auto mappedNode1 = from_map.find(edge->from());
        if (mappedNode1 != from_map.end()) {
            const NodeMap& to_map = edge->to_end() ? _end_node_map : _start_node_map;
            auto mappedNode2 = to_map.find(edge->to());
            if (mappedNode2 != to_map.end()) {
                // up to 2 mappings for each node, so 4 possible edges:
                if (mappedNode1->second.first != -1 && mappedNode2->second.first != -1) {
                    _call_graph.create_edge(mappedNode1->second.first, mappedNode2->second.first,
                                            edge->from_start(), edge->to_end());
                }
                if (mappedNode1->second.first != -1 && mappedNode2->second.second != -1) {
                    _call_graph.create_edge(mappedNode1->second.first, mappedNode2->second.second,
                                            edge->from_start(), edge->to_end());
                }
                if (mappedNode1->second.second != -1 && mappedNode2->second.first != -1) {
                    _call_graph.create_edge(mappedNode1->second.second, mappedNode2->second.first,
                                            edge->from_start(), edge->to_end());
                }
                if (mappedNode1->second.second != -1 && mappedNode2->second.second != -1) {
                    _call_graph.create_edge(mappedNode1->second.second, mappedNode2->second.second,
                                            edge->from_start(), edge->to_end());
                }
            }
        }
    };
    _graph->for_each_edge(map_edge);

    // make sure paths are saved
    _call_graph.paths.rebuild_node_mapping();
    _call_graph.paths.rebuild_mapping_aux();
    _call_graph.paths.to_graph(_call_graph.graph);

    // fix ids
    //_call_graph.sort();
    //_call_graph.compact_ids();
}

void Caller::call_base_pileup(const NodePileup& np, int64_t offset) {
    const BasePileup& bp = np.base_pileup(offset);
    
    // parse the pilueup structure
    vector<pair<int, int> > base_offsets;
    vector<pair<int, int> > indel_offsets;
    Pileups::parse_base_offsets(bp, base_offsets, indel_offsets, indel_offsets);

    // compute top two most frequent bases and their counts
    char top_base;
    int top_count;
    char second_base;
    int second_count;
    compute_top_frequencies(bp, base_offsets, top_base, top_count, second_base, second_count);

    // note first and second base will be upper case too
    char ref_base = ::toupper(bp.ref_base());

    // test against thresholding heuristics
    if ((double)(top_count + second_count) / (double)base_offsets.size() >= _min_frac) {

        // compute max likelihood snp genotype.  it will be one of the three combinations
        // of the top two bases (we don't care about case here)
        pair<char, char> g = mp_snp_genotype(bp, base_offsets, top_base, second_base);


        // update the node calls
        if (top_count >= _min_support) {
            if (g.first != ref_base) {
                _node_calls[offset].first = g.first;
            } else {
                _node_calls[offset].first = '.';
            }
        }
        if (second_count >= _min_support) {
            if (g.second != ref_base && g.second != g.first) {
                _node_calls[offset].second = g.second;
            } else {
                _node_calls[offset].second = '.';
            }
        }
    }
}

void Caller::compute_top_frequencies(const BasePileup& bp,
                                     const vector<pair<int, int> >& base_offsets,
                                     char& top_base, int& top_count,
                                     char& second_base, int& second_count) {

    const string& bases = bp.bases();
    // frequency of each base (nidx / idxn converts to from / int)
    int hist[5] = {0};
    for (auto i : base_offsets) {
        char base = Pileups::extract_match(bp, i.first);
        ++hist[nidx(base)];
    }
    
    int first = max_element(hist, hist + 4) - hist;
    int second = first == 0 ? 1 : 0;
    for (int i = 0; i < 4; ++i) {
        if (i != first && hist[i] > hist[second]) {
            second = i;
        }
    }

    // break ties with reference
    int refidx = nidx(bp.ref_base());
    if (hist[first] == hist[refidx]) {
        first = refidx;
    } else if (hist[second] == hist[refidx]) {
        second = refidx;
    }

    top_base = idxn(first);
    top_count = hist[first];
    second_base = idxn(second);
    second_count = hist[second];
}

// Estimate the most probable snp genotype
pair<char, char> Caller::mp_snp_genotype(const BasePileup& bp,
                                         const vector<pair<int, int> >& base_offsets,
                                         char top_base, char second_base) {
    char ref_base = ::toupper(bp.ref_base());

    // gotta do better than this:
    pair<char, char> mp_genotype(ref_base, ref_base);
    double mp = _min_log_likelihood + _hom_log_prior;

    // genotype with 0 top_bases
    double gl = genotype_log_likelihood(bp, base_offsets, 0, top_base, second_base);
    double p = _hom_log_prior + gl;
    if (p > mp) {
        mp = p;
        mp_genotype = make_pair(second_base, second_base);
    }

    // genotype with 1 top_base
    gl = genotype_log_likelihood(bp, base_offsets, 1, top_base, second_base);
    p = _het_log_prior + gl;
    if (p > mp) {
        mp = p;
        mp_genotype = make_pair(top_base, second_base);
    }

    // genotype with 2 top_bases
    gl = genotype_log_likelihood(bp, base_offsets, 2, top_base, second_base);
    p = _hom_log_prior + gl;
    if (p > mp) {
        mp = p;
        mp_genotype = make_pair(top_base, top_base);
    }
    
    // note, we're throwing away the probabilty value (mp) here.
    // should figure out where to stick it in the output. 
    return mp_genotype;
}

// This is Equation 2 (tranformed to log) from
// A statistical framework for SNP calling ... , Heng Li, Bioinformatics, 2011
// http://bioinformatics.oxfordjournals.org/content/27/21/2987.full
double Caller::genotype_log_likelihood(const BasePileup& bp,
                                       const vector<pair<int, int> >& base_offsets,
                                       double g, char first, char second) {
    double m = 2.; // always assume two alleles

    double log_likelihood = log(0.25); // 1 / m^2, where m = ploidy = 2;

    const string& bases = bp.bases();
    const string& quals = bp.qualities();
    double perr;

    for (int i = 0; i < base_offsets.size(); ++i) {
        char base = Pileups::extract_match(bp, base_offsets[i].first);
        char qual = base_offsets[i].second >= 0 ? quals[base_offsets[i].second] : _default_quality;
        perr = phred2prob(qual);
        if (base == first) {
            log_likelihood += safe_log((m - g) * perr + g * (1. - perr));
        } else if (base == second) {
            log_likelihood += safe_log((m - g) * (1. - perr) + g * perr);
        } else {
            log_likelihood += safe_log(perr * perr);
        }
    }

    return log_likelihood;
}


void Caller::create_node_calls(const NodePileup& np) {
    
    int n = _node->sequence().length();
    const string& seq = _node->sequence();
    int cur = 0;
    int cat = call_cat(_node_calls[cur]);
    NodePair prev_nodes(-1, -1);

    // scan contiguous chunks of a node with same call
    // (note: snps will always be 1-base -- never merged)
    for (int next = 1; next <= n; ++next) {
        int next_cat = next == n ? -1 : call_cat(_node_calls[next]);
        if (cat == 2 || cat != next_cat) {
            NodePair new_nodes(-1, -1);
            bool secondary_snp = false;

            // process first genotype if it's not missing
            if (_node_calls[cur].first == '.') {
                // add single node for stretch of reference node
                string new_seq = seq.substr(cur, next - cur);
                new_nodes.first = ++_max_id;
                _call_graph.create_node(new_seq, new_nodes.first);
            } else if (_node_calls[cur].first != '-') {
                // add snp node
                assert(next - cur == 1);
                string new_seq(1, _node_calls[cur].first);
                new_nodes.first = ++_max_id;
                _call_graph.create_node(new_seq, new_nodes.first);
                create_snp_path(new_nodes.first, secondary_snp);
                secondary_snp = true;
            }

            // process second genotype if difference from first
            if (_node_calls[cur].second != _node_calls[cur].first) {
                if (_node_calls[cur].second == '.') {
                    // add single node for stretch of reference node
                    string new_seq = seq.substr(cur, next - cur);
                    new_nodes.second = ++_max_id;
                    _call_graph.create_node(new_seq, new_nodes.second);
                } else if (_node_calls[cur].second != '-') {
                    // add snp node
                    assert(next - cur == 1);
                    string new_seq(1, _node_calls[cur].second);
                    new_nodes.second = ++_max_id;
                    _call_graph.create_node(new_seq, new_nodes.second);
                    create_snp_path(new_nodes.second, secondary_snp);
                }
            }
            
            // update maps if new node abuts end of original node
            // so that edges can be updated later on:
            if (new_nodes.first != -1 || new_nodes.second != -1) {
                if (cur == 0) {
                    _start_node_map[_node->id()] = new_nodes;
                }
                if (next == n) {
                    _end_node_map[_node->id()] = new_nodes;
                }
            }

            // add edges
            if (prev_nodes.first != -1 && new_nodes.first != -1) {
                _call_graph.create_edge(prev_nodes.first, new_nodes.first);
            }
            if (prev_nodes.first != -1 && new_nodes.second != -1) {
                _call_graph.create_edge(prev_nodes.first, new_nodes.second);
            }
            if (prev_nodes.second != -1 && new_nodes.first != -1) {
                _call_graph.create_edge(prev_nodes.second, new_nodes.first);
            }
            if (prev_nodes.second != -1 && new_nodes.second != -1) {
                _call_graph.create_edge(prev_nodes.second, new_nodes.second);
            }

            // shift right
            cur = next;
            cat = next_cat;
            prev_nodes = new_nodes;
        }
    }
}

void Caller::create_snp_path(int64_t snp_node, bool secondary_snp) {

    // for now we don't write secdonary snp, so we have 1 path per *site*
    // and counting paths will give us somethign comparable to snp count
    // from bcftools
    if (!secondary_snp) {
        stringstream name;
        name << "SNP_" << snp_node;

        Mapping mapping;
        Position* pos = mapping.mutable_position();
        // make path that covers node forward with no edits.  not super
        // useful but will use to count snps... 
        pos->set_node_id(snp_node);
        pos->set_offset(0);
        mapping.set_is_reverse(false);
        
        // note: create_path doesn't seem to work.. too rushed to look into
        //list<Mapping>& mappings = _call_graph.paths.create_path(name.str());

        list<Mapping> mappings;
        mappings.push_back(mapping);
        _call_graph.paths._paths.insert(make_pair(name.str(), mappings));
    }
}

}
