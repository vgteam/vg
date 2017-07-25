// copied from https://github.com/adamnovak/glenn2vcf/blob/master/main.cpp
// as this logic really belongs in vg call

#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <set>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <getopt.h>

#include "vg.hpp"
#include "index.hpp"
#include "Variant.h"
#include "genotypekit.hpp"
#include "snarls.hpp"
#include "path_index.hpp"
#include "caller.hpp"
#include "stream.hpp"
#include "nested_traversal_finder.hpp"

//#define debug

namespace vg {

// How many bases may we put in an allele in VCF if we expect GATK to be able to
// parse it?
// 0 means no maximum is enforced.
const static int MAX_ALLELE_LENGTH = 0;

// Minimum log likelihood
const static double LOG_ZERO = (double)-1e100;

// convert to string using stringstream (to replace to_string when we want sci. notation)
template <typename T>
string to_string_ss(T val) {
    stringstream ss;
    ss << val;
    return ss.str();
}

/**
 * We need to suppress overlapping variants, but interval trees are hard to
 * write. This accomplishes the collision check with a massive bit vector.
 */
struct IntervalBitfield {
    // Mark every position that's used in a variant
    vector<bool> used;
    
    /**
     * Make a new IntervalBitfield covering a region of the specified length.
     */
    inline IntervalBitfield(size_t length) : used(length) {
        // Nothing to do
    }
    
    /**
     * Scan for a collision (O(n) in interval length)
     */
    inline bool collides(size_t start, size_t pastEnd) {
        for(size_t i = start; i < pastEnd; i++) {
            if(used[i]) {
                return true;
            }
        }
        return(false);
    }
    
    /**
     * Take up an interval.
     */
    inline void add(size_t start, size_t pastEnd) {
        for(size_t i = start; i < pastEnd; i++) {
            used[i] = true;
        }
    }
};

/**
 * Get the strand bias of a Support.
 */
double strand_bias(const Support& support) {
    return max(support.forward(), support.reverse()) / (support.forward() + support.reverse());
}

/**
 * Make a letter into a full string because apparently that's too fancy for the
 * standard library.
 */
string char_to_string(const char& letter) {
    string toReturn;
    toReturn.push_back(letter);
    return toReturn;
}

/**
 * Write a minimal VCF header for a file with the given samples, and the given
 * contigs with the given lengths.
 */
void write_vcf_header(ostream& stream, const vector<string>& sample_names,
    const vector<string>& contig_names, const vector<size_t>& contig_sizes,
    int min_mad_for_filter, int max_dp_for_filter, double max_dp_multiple_for_filter,
    double max_local_dp_multiple_for_filter, double min_ad_log_likelihood_for_filter) {
    
    stream << "##fileformat=VCFv4.2" << endl;
    stream << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" << endl;
    stream << "##INFO=<ID=XREF,Number=0,Type=Flag,Description=\"Present in original graph\">" << endl;
    stream << "##INFO=<ID=XSEE,Number=.,Type=String,Description=\"Original graph node:offset cross-references\">" << endl;
    stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << endl;
    stream << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << endl;
    stream << "##FILTER=<ID=lowad,Description=\"Variant does not meet minimum allele read support threshold of " << min_mad_for_filter << "\">" <<endl;
    stream << "##FILTER=<ID=highabsdp,Description=\"Variant has total depth greater than " << max_dp_for_filter << "\">" <<endl;
    stream << "##FILTER=<ID=highreldp,Description=\"Variant has total depth greater than "
        << max_dp_multiple_for_filter << " times global baseline\">" <<endl;
    stream << "##FILTER=<ID=highlocaldp,Description=\"Variant has total depth greater than "
        << max_local_dp_multiple_for_filter << " times local baseline\">" <<endl;
    stream << "##FILTER=<ID=lowxadl,Description=\"Variant has AD log likelihood less than "
        << min_ad_log_likelihood_for_filter << "\">" <<endl;
    stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
    stream << "##FORMAT=<ID=XDP,Number=2,Type=Integer,Description=\"Expected Local and Global Depth\">" << endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    stream << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << endl;
    stream << "##FORMAT=<ID=XADL,Number=1,Type=Float,Description=\"Likelihood of allelic depths for called alleles\">" << endl;
    stream << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Forward and reverse support for ref and alt alleles.\">" << endl;
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    stream << "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">" << endl;
    stream << "##FORMAT=<ID=AL,Number=.,Type=Float,Description=\"Allelic likelihoods for the ref and alt alleles in the order listed\">" << endl;
    
    for(size_t i = 0; i < contig_names.size(); i++) {
        // Announce the contigs as well.
        stream << "##contig=<ID=" << contig_names.at(i) << ",length=" << contig_sizes.at(i) << ">" << endl;
    }
    
    // Now the column header line
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (auto& sample_name : sample_names) {
        // Append columns for all the samples
        stream << "\t" << sample_name;
    }
    // End the column header line
    stream << endl;
}

/**
 * Return true if a variant may be output, or false if this variant is valid but
 * the GATK might choke on it.
 *
 * Mostly used to throw out variants with very long alleles, because GATK has an
 * allele length limit. How alleles that really *are* 1 megabase deletions are
 * to be specified to GATK is left as an exercise to the reader.
 */
bool can_write_alleles(vcflib::Variant& variant) {

    for(auto& allele : variant.alleles) {
        if(MAX_ALLELE_LENGTH > 0 && allele.size() > MAX_ALLELE_LENGTH) {
            return false;
        }
    }
    return true;
}

/**
 * Return true if a mapping is a perfect match, and false if it isn't.
 */
bool mapping_is_perfect_match(const Mapping& mapping) {
    for (auto edit : mapping.edit()) {
        if (edit.from_length() != edit.to_length() || !edit.sequence().empty()) {
            // This edit isn't a perfect match
            return false;
        }
    }
    
    // If we get here, all the edits are perfect matches.
    // Note that Mappings with no edits at all are full-length perfect matches.
    return true;
}

/**
 * Given a collection of pileups by original node ID, and a set of original node
 * id:offset cross-references in both ref and alt categories, produce a VCF
 * comment line giving the pileup for each of those positions on those nodes.
 * Includes a trailing newline if nonempty.
 *
 * TODO: VCF comments aren't really a thing.
 */
string get_pileup_line(const map<int64_t, NodePileup>& node_pileups,
    const set<pair<int64_t, size_t>>& refCrossreferences,
    const set<pair<int64_t, size_t>>& altCrossreferences) {
    // We'll make a stringstream to write to.
    stringstream out;
    
    out << "#";
    
    for(const auto& xref : refCrossreferences) {
        // For every cross-reference
        if(node_pileups.count(xref.first) && node_pileups.at(xref.first).base_pileup_size() > xref.second) {
            // If we have that base pileup, grab it
            auto basePileup = node_pileups.at(xref.first).base_pileup(xref.second);
            
            out << xref.first << ":" << xref.second << " (ref) " << basePileup.bases() << "\t";
        }
        // Nodes with no pileups (either no pileups were provided or they didn't
        // appear/weren't visited by reads) will not be mentioned on this line
    }
    
    for(const auto& xref : altCrossreferences) {
        // For every cross-reference
        if(node_pileups.count(xref.first) && node_pileups.at(xref.first).base_pileup_size() > xref.second) {
            // If we have that base pileup, grab it
            auto basePileup = node_pileups.at(xref.first).base_pileup(xref.second);
            
            out << xref.first << ":" << xref.second << " (alt) " << basePileup.bases() << "\t";
        }
        // Nodes with no pileups (either no pileups were provided or they didn't
        // appear/weren't visited by reads) will not be mentioned on this line
    }
    // TODO: make these nearly-identical loops a loop or a lambda or something.
    
    if(out.str().size() > 1) {
        // We actually found something. Send it out with a trailing newline
        out << endl;
        return out.str();
    } else {
        // Give an empty string.
        return "";
    }
}

Call2Vcf::PrimaryPath::PrimaryPath(AugmentedGraph& augmented, const string& ref_path_name, size_t ref_bin_size):
    ref_bin_size(ref_bin_size), index(augmented.graph, ref_path_name, true), name(ref_path_name)  {

    // Follow the reference path and extract indexes we need: index by node ID,
    // index by node start, and the reconstructed path sequence.
    PathIndex index(augmented.graph, ref_path_name, true);

    if (index.sequence.size() == 0) {
        // No empty reference paths allowed
        throw runtime_error("Reference path cannot be empty");
    }

    // Store support binned along reference path;
    // Last bin extended to include remainder
    ref_bin_size = min(ref_bin_size, index.sequence.size());
    if (ref_bin_size <= 0) {
        // No zero-sized bins allowed
        throw runtime_error("Reference bin size must be 1 or larger");
    }
    // Start out all the bins empty.
    binned_support = vector<Support>(max(1, int(index.sequence.size() / ref_bin_size)), Support());
    
    // Crunch the numbers on the reference and its read support. How much read
    // support in total (node length * aligned reads) does the primary path get?
    total_support = Support();
    for(auto& pointerAndSupport : augmented.node_supports) {
        if(index.by_id.count(pointerAndSupport.first->id())) {
            // This is a primary path node. Add in the total read bases supporting it
            total_support += pointerAndSupport.first->sequence().size() * pointerAndSupport.second;
            
            // We also update the total for the appropriate bin
            size_t bin = index.by_id[pointerAndSupport.first->id()].first / ref_bin_size;
            if (bin == binned_support.size()) {
                --bin;
            }
            binned_support[bin] = binned_support[bin] + 
                pointerAndSupport.first->sequence().size() * pointerAndSupport.second;
        }
    }
    
    // Average out the support bins too (in place)
    min_bin = 0;
    max_bin = 0;
    for (int i = 0; i < binned_support.size(); ++i) {
        // Compute the average over the bin's actual size
        binned_support[i] = binned_support[i] / (
            i < binned_support.size() - 1 ? (double)ref_bin_size :
            (double)(ref_bin_size + index.sequence.size() % ref_bin_size));
            
        // See if it's a min or max
        if (binned_support[i] < binned_support[min_bin]) {
            min_bin = i;
        }
        if (binned_support[i] > binned_support[max_bin]) {
            max_bin = i;
        }
    }

}

const Support& Call2Vcf::PrimaryPath::get_support_at(size_t primary_path_offset) const {
    return get_bin(get_bin_index(primary_path_offset));
}
        
size_t Call2Vcf::PrimaryPath::get_bin_index(size_t primary_path_offset) const {
    // Find which coordinate bin the position is in
    size_t bin = primary_path_offset / ref_bin_size;
    if (bin == get_total_bins()) {
        --bin;
    }
    return bin;
}
    
size_t Call2Vcf::PrimaryPath::get_min_bin() const {
    return min_bin;
}
    
size_t Call2Vcf::PrimaryPath::get_max_bin() const {
    return max_bin;
}
    
const Support& Call2Vcf::PrimaryPath::get_bin(size_t bin) const {
    return binned_support[bin];
}
        
size_t Call2Vcf::PrimaryPath::get_total_bins() const {
    return binned_support.size();
}
        
Support Call2Vcf::PrimaryPath::get_average_support() const {
    return get_total_support() / get_index().sequence.size();
}

Support Call2Vcf::PrimaryPath::get_average_support(const map<string, PrimaryPath>& paths) {
    // Track the total support overall
    Support total;
    // And the total number of bases
    size_t bases;
    
    for (auto& kv : paths) {
        // Sum over all paths
        total += kv.second.get_total_support();
        bases += kv.second.get_index().sequence.size();
    }
    
    // Then divide
    return total / bases;
}
        
Support Call2Vcf::PrimaryPath::get_total_support() const {
    return total_support;
}
  
PathIndex& Call2Vcf::PrimaryPath::get_index() {
    return index;
}
    
const PathIndex& Call2Vcf::PrimaryPath::get_index() const {
    return index;
}

const string& Call2Vcf::PrimaryPath::get_name() const {
    return name;
}

map<string, Call2Vcf::PrimaryPath>::iterator Call2Vcf::find_path(const Snarl& site, map<string, PrimaryPath>& primary_paths) {
    for(auto i = primary_paths.begin(); i != primary_paths.end(); ++i) {
        // Scan the whole map with an iterator
        
        if (i->second.get_index().by_id.count(site.start().node_id()) &&
            i->second.get_index().by_id.count(site.end().node_id())) {
            // This path threads through this site
            return i;
        }
    }
    // Otherwise we hit the end and found no path that this site can be strung
    // on.
    return primary_paths.end();
}

/**
 * Trace out the given traversal, handling nodes, child snarls, and edges
 * associated with particular visit numbers.
 */
void trace_traversal(const SnarlTraversal& traversal, const Snarl& site, function<void(size_t,id_t)> handle_node,
    function<void(size_t,NodeSide,NodeSide)> handle_edge, function<void(size_t,Snarl)> handle_child) {

    // Must at least have start and end
    assert(traversal.visits_size() >= 2);
    
    // Look at the edge leading from the start (also handles deletion traversals)
    handle_edge(0, to_right_side(traversal.visits(0)), to_left_side(traversal.visits(1)));
    
    for(int64_t i = 1; i < traversal.visits_size() - 1; i++) {
        // For all the (internal) visits...
        auto& visit = traversal.visits(i);
        
        if (visit.node_id() != 0) {
            // This is a visit to a node
            
            // Find the node
            handle_node(i - 1, visit.node_id());
        } else {
            // This is a snarl
            handle_child(i - 1, visit.snarl());
        }
        
        auto& next_visit = traversal.visits(i + 1);
        
        if (visit.node_id() == 0 && next_visit.node_id() == 0 &&
            to_right_side(visit).flip() == to_left_side(next_visit)) {
            
            // These are two back-to-back child snarl visits, which
            // share a node and have no connecting edge.
#ifdef debug
            cerr << "No edge needed for back-to-back child snarls" << endl;
#endif
            
        }
        else {
            // Do the edge to it
            handle_edge(i - 1, to_right_side(visit), to_left_side(next_visit));
        }
    }

}

/**
 * Get the min support, total support, bp size (to divide total by for average
 * support), and min likelihood for a traversal, optionally special-casing the
 * material used by another traversal. Material used by another traversal only
 * makes half its coverage available to this traversal.
 */
tuple<Support, Support, size_t> get_traversal_support(AugmentedGraph& augmented,
    SnarlManager& snarl_manager, const Snarl& site, const SnarlTraversal& traversal,
    const SnarlTraversal* already_used = nullptr) {

#ifdef debug
    cerr << "Evaluate traversal: " << endl;
    for (size_t i = 0; i < traversal.visits_size(); i++) {
        cerr << "\t" << pb2json(traversal.visits(i)) << endl;
    }
    if (already_used != nullptr) {
        cerr << "Need to share: " << endl;
        for (size_t i = 0; i < already_used->visits_size(); i++) {
            cerr << "\t" << pb2json(already_used->visits(i)) << endl;
        }
    }
#endif

    // First work out the stuff we need to share
    set<id_t> shared_nodes;
    set<Snarl> shared_children;
    set<Edge*> shared_edges;
    if (already_used != nullptr) {
        // Mark all the nodes and edges that the other traverasl uses.
        trace_traversal(*already_used, site, [&](size_t i, id_t node) {
            shared_nodes.insert(node);
        }, [&](size_t i, NodeSide end1, NodeSide end2) {
            shared_edges.insert(augmented.graph.get_edge(end1, end2));
        }, [&](size_t i, Snarl child) {
            shared_children.insert(child);
        });
    }
    
    
    // Compute min and total supports, and bp sizes, for all the visits by
    // number.
    size_t record_count = max(1, traversal.visits_size() - 2);
    // What's the min support observed at every visit (inclusing edges)?
    vector<Support> min_supports(record_count, make_support(INFINITY, INFINITY, INFINITY));
    // And the total support (ignoring edges)?
    vector<Support> total_supports(record_count, Support());
    // And the bp size of each visit
    vector<size_t> visit_sizes(record_count, 0);
    
    // Don't count nodes shared between child snarls more than once.
    set<Node*> coverage_counted;
    
    trace_traversal(traversal, site, [&](size_t i, id_t node_id) {
        // Find the node
        Node* node = augmented.graph.get_node(node_id);
    
        // Grab this node's total support along its length
        // Make sure to only use half the support if the node is shared
        total_supports[i] += augmented.get_support(node) * node->sequence().size() * (shared_nodes.count(node_id) ? 0.5 : 1.0);
        // And its size
        visit_sizes[i] += node->sequence().size();
        
        // And update its min support
        min_supports[i] = support_min(min_supports[i], augmented.get_support(node) * (shared_nodes.count(node_id) ? 0.5 : 1.0));
        
    }, [&](size_t i, NodeSide end1, NodeSide end2) {
        // This is an edge
        Edge* edge = augmented.graph.get_edge(end1, end2);
        assert(edge != nullptr);
        
        // Count as 1 base worth for the total/average support
        // Make sure to only use half the support if the edge is shared
        total_supports[i] += augmented.get_support(edge) * (shared_edges.count(edge) ? 0.5 : 1.0);
        visit_sizes[i] += 1;
        
        // Min in its support
        min_supports[i] = support_min(min_supports[i], augmented.get_support(edge) * (shared_edges.count(edge) ? 0.5 : 1.0));
    }, [&](size_t i, Snarl child) {
        // This is a child snarl, so get its max support.
        
        Support child_max;
        size_t child_size = 0;
        for (Node* node : snarl_manager.deep_contents(snarl_manager.manage(child),
            augmented.graph, true).first) {
            // For every node in the child
            
            if (coverage_counted.count(node)) {
                // Already used by another child snarl on this traversal
                continue;
            }
            // Claim this node for this child.
            coverage_counted.insert(node);
            
            // How many distinct reads must use the child, given the distinct reads on this node?
            child_max = support_max(child_max, augmented.get_support(node));
            
            // Add in the node's size to the child
            child_size += node->sequence().size();
            
#ifdef debug
            cerr << "From child snarl node " << node->id() << " get "
                << augmented.get_support(node) << " for distinct " << child_max << endl;
#endif
        }
        
        if (shared_children.count(child)) {
            // Make sure to halve the support if the child is shared
            child_max *= 0.5;
        }
        
        // Smoosh support over the whole child
        total_supports[i] += child_max * child_size;
        visit_sizes[i] += child_size;
        
        if (child_size != 0) {
            // We actually have some nodes to our name.
            min_supports[i] = support_min(min_supports[i], child_max);
        }
        
    });
    
    // Now aggregate across visits and their edges

    // What's the total support for this traversal?
    Support total_support;
    for (auto& support : total_supports) {
        total_support += support;
    }
    
    // And the length over which we have it (for averaging)
    size_t total_size = 0;
    for (auto& size : visit_sizes) {
        total_size += size;
    }
    
    // And the min support?
    Support min_support = make_support(INFINITY, INFINITY, INFINITY);
    for (auto& support : min_supports) {
        min_support = support_min(min_support, support);
    }
        
    if (min_support.forward() == INFINITY || min_support.reverse() == INFINITY) {
        // If we have actually no material, say we have actually no support
        min_support = Support();
    }
        
    // Spit out the supports, the size in bases observed.
    return tie(min_support, total_support, total_size);
        
}

/** Get the support for each traversal in a list, using average_support_switch_threshold
    to decide if we use the minimum or average */
tuple<vector<Support>, vector<size_t> > Call2Vcf::get_traversal_supports_and_sizes(
    AugmentedGraph& augmented, SnarlManager& snarl_manager, const Snarl& site,
    const vector<SnarlTraversal>& traversals, const SnarlTraversal* minus_traversal) {

    // How long is the longest traversal?
    // Sort of approximate because of the way nested site sizes are estimated.
    size_t longest_traversal_length = 0;
    
    // And the shortest one?
    size_t shortest_traversal_length = numeric_limits<size_t>::max();

    // Calculate average and min support for all the traversals of this snarl.
    vector<Support> min_supports;
    vector<Support> average_supports;
    vector<size_t> sizes;
    for(auto& traversal : traversals) {
        // Go through all the SnarlTraversals for this Snarl
        
        // What's the total support for this traversal?
        Support total_support;
        // And the length over which we have it (for averaging)
        size_t total_size;
        // And the min support?
        Support min_support;
        // Trace the traversal and get its support
        tie(min_support, total_support, total_size) = get_traversal_support(
            augmented, snarl_manager, site, traversal, minus_traversal);
            
        // Add average and min supports to vectors. Note that average support
        // ignores edges.
        min_supports.push_back(min_support);
        average_supports.push_back(total_size != 0 ? total_support / total_size : Support());
        
#ifdef debug
        cerr << "Min: " << min_support << " Total: " << total_support << " Average: " << average_supports.back() << endl;
#endif

        // Remember a new longest traversal length
        longest_traversal_length = max(longest_traversal_length, total_size);  
        // And a new shortest one
        shortest_traversal_length = min(shortest_traversal_length, total_size);
        // and the current size
        sizes.push_back(total_size);
    }
    
#ifdef debug
    cerr << "Min vs. average" << endl;
#endif
    for (size_t i = 0; i < average_supports.size(); i++) {
#ifdef debug
        cerr << "\t" << min_supports.at(i) << " vs. " << average_supports.at(i) << endl;
#endif
        // We should always have a higher average support than minumum support
        assert(support_val(average_supports.at(i)) >= support_val(min_supports.at(i)));
    }

    return (longest_traversal_length > average_support_switch_threshold || use_average_support) ?
        tie(average_supports, sizes) :
        tie(min_supports, sizes);
}

vector<SnarlTraversal> Call2Vcf::find_best_traversals(AugmentedGraph& augmented,
        SnarlManager& snarl_manager, TraversalFinder* finder, const Snarl& site,
        const Support& baseline_support, size_t copy_budget, function<void(const Locus&, const Snarl*)> emit_locus) {

    // We need to be an ultrabubble for the traversal finder to work right.
    // TODO: generalize it
    assert(site.type() == ULTRABUBBLE);


#ifdef debug
    cerr << "Site " << site << endl;
#endif


    // Get traversals of this Snarl, with Visits to child Snarls.
    // The 0th is always the reference traversal if we're on a primary path
    vector<SnarlTraversal> here_traversals = finder->find_traversals(site);
    

#ifdef debug
    cerr << "Found " << here_traversals.size() << " traversals" << endl;
#endif

    
    // Make a Locus to hold all our stats for the different traversals
    // available. The 0th will always be the ref traversal if we're on a primary
    // path.
    Locus locus;
    
    vector<Support> supports;
    vector<size_t> traversal_sizes;
    // Calculate the support for all the traversals of this snarl.
    tie(supports, traversal_sizes) = get_traversal_supports_and_sizes(
        augmented, snarl_manager, site, here_traversals);
        
    for (auto& support : supports) {
        // Blit supports over to the locus
        *locus.add_support() = support;
    }
    
    ////////////////////////////////////////////////////////////////////////////

    // look at all the paths for the site and pick the best one    
    function<int(const vector<Support>&, vector<int>)> get_best_allele = [this](
        const vector<Support>& supports, vector<int> skips) {
        int best_allele = -1;
        for(size_t i = 0; i < supports.size(); i++) {
            if(std::find(skips.begin(), skips.end(), i) == skips.end() && (
                   best_allele == -1 || support_val(supports[best_allele]) <= support_val(supports[i]))) {
                best_allele = i;
            }
        }
        return best_allele;
    };

    // Now look at all the paths for the site and pick the best one
    int best_allele = get_best_allele(supports, {});
    
    // We should always have a best allele; we may sometimes have a second best.
    assert(best_allele != -1);
    
#ifdef debug
    cerr << "Choose best allele: " << best_allele << endl;
#endif
    
    // Then recalculate supports assuming we can't count anything shared with that best traversal
    vector<Support> additional_supports;
    tie(additional_supports, std::ignore) = get_traversal_supports_and_sizes(
        augmented, snarl_manager, site, here_traversals, &here_traversals.at(best_allele));
    
    // Then pick the second best one
    int second_best_allele = get_best_allele(additional_supports, {best_allele});

#ifdef debug
    cerr << "Choose second best allele: " << second_best_allele << endl;
#endif

    // Hack for special case where we want to call a multiallelic alt even if the reference
    // has better support than one or both alts
    vector<Support> tertiary_supports;
    int third_best_allele = -1;
    if (second_best_allele != -1) {
        tie(tertiary_supports, std::ignore) = get_traversal_supports_and_sizes(
            augmented, snarl_manager, site, here_traversals, &here_traversals.at(second_best_allele));
        third_best_allele = get_best_allele(tertiary_supports, {best_allele, second_best_allele});
    }    

    // Decide if we're an indel by looking at the traversal sizes
    bool is_indel = traversal_sizes[best_allele] != traversal_sizes[0] ||
        (second_best_allele != -1 && traversal_sizes[0] != traversal_sizes[second_best_allele]);
    bool is_indel_ma_2 = (second_best_allele != -1 && traversal_sizes[0] != traversal_sizes[second_best_allele]);
    bool is_indel_ma_3 = (third_best_allele != -1 && traversal_sizes[0] != traversal_sizes[third_best_allele]);
    
    ////////////////////////////////////////////////////////////////////////////
    
    // Now make a genotype call at this site, up to the allowed copy number
    
    // TODO: Work out how to detect indels when there are nested sites and
    // enable the indel bias multiple again.
    double bias_multiple = 1.0;
    
    // How much support do we have for the top two alleles?
    Support site_support = supports.at(best_allele);
    if(second_best_allele != -1) {
        site_support += supports.at(second_best_allele);
    }
    
    // Pull out the different supports. Some of them may be the same.
    Support best_support = supports.at(best_allele);
    Support second_best_support; // Defaults to 0
    if(second_best_allele != -1) {
        second_best_support = supports.at(second_best_allele);
    }
    Support third_best_support;
    if (third_best_allele != -1) {
        third_best_support = supports.at(third_best_allele);
    }
    
    // As we do the genotype, we also compute the likelihood. Holds
    // likelihood log 10. Starts out at "completely wrong".
    double gen_likelihood = -1 * INFINITY;

    // Minimum allele depth of called alleles
    double min_site_support = 0;
    
    // This is where we'll put the genotype. We only actually add it to the
    // Locus if we are confident enough to actually call.
    Genotype genotype;
    
    // We're going to make some really bad calls at low depth. We can
    // pull them out with a depth filter, but for now just elide them.
    if (support_val(site_support) >= support_val(baseline_support) * min_fraction_for_call * ((double) copy_budget) / 2) {
        // We have enough to emit a call here.
        
        // If best and second best are close enough to be het, we call het.
        // Otherwise, we call hom best.
        
        double bias_limit;
        if (best_allele == 0) {
            // Use ref bias limit
            
            // We decide closeness differently depending on whether best is ref
            // or not. In practice, we use this to slightly penalize homozygous
            // ref calls (by setting max_ref_het_bias higher than max_het_bias)
            // and rather make a less supported alt call instead.  This boost
            // max sensitivity, and because everything is homozygous ref by
            // default in VCF, any downstream filters will effectively reset
            // these calls back to homozygous ref. TODO: This shouldn't apply
            // when off the primary path!
            bias_limit = max_ref_het_bias;
        } else if (is_indel) {
            // This is an indel
            // Use indel bias limit
            bias_limit = max_indel_het_bias;
        } else {
            // Use normal het bias limit
            bias_limit = max_het_bias;
        }
        
#ifdef debug
        cerr << best_allele << ", " << best_support << " and "
            << second_best_allele << ", " << second_best_support << endl;
        
        if (support_val(second_best_support) > 0) {
            cerr << "Bias: (limit " << bias_limit * bias_multiple << "):"
                << support_val(best_support)/support_val(second_best_support) << endl;
        }
        
        cerr << bias_limit * bias_multiple * support_val(second_best_support) << " vs "
            << support_val(best_support) << endl;
            
        cerr << total(second_best_support) << " vs " << min_total_support_for_call << endl;
#endif

        // Call 1/2 : REF-Alt1/Alt2 even if Alt2 has only third best support
        if (copy_budget >= 2 &&
            best_allele == 0 && 
            third_best_allele > 0 &&
            is_indel_ma_3 &&
            max_indel_ma_bias * bias_multiple * support_val(third_best_support) >= support_val(best_support) &&
            total(second_best_support) >= min_total_support_for_call &&
            total(third_best_support) >= min_total_support_for_call) {
            // There's a second best allele and third best allele, and it's not too biased to call,
            // and both alleles exceed the minimum to call them present, and the
            // second-best and third-best alleles have enough support that it won't torpedo the
            // variant.
            
#ifdef debug
            cerr << "Call as second best/third best" << endl;
#endif
            // Say both are present
            genotype.add_allele(second_best_allele);
            genotype.add_allele(third_best_allele);
                        
            // Get minimum support for filter (not assuming it's second_best just to be sure)
            min_site_support = min(total(second_best_support), total(third_best_support));
            
            // Make the call
            *locus.add_genotype() = genotype;
        }
        else if (copy_budget >= 2 &&
            second_best_allele != -1 &&
            bias_limit * bias_multiple * support_val(second_best_support) >= support_val(best_support) &&
            total(best_support) >= min_total_support_for_call &&
            total(second_best_support) >= min_total_support_for_call) {
            // There's a second best allele, and it's not too biased to call,
            // and both alleles exceed the minimum to call them present, and the
            // second-best allele has enough support that it won't torpedo the
            // variant.
            
#ifdef debug
            cerr << "Call as best/second best" << endl;
#endif
            
            // Say both are present
            genotype.add_allele(best_allele);
            genotype.add_allele(second_best_allele);
                        
            // Get minimum support for filter (not assuming it's second_best just to be sure)
            min_site_support = min(total(second_best_support), total(best_support));
            
            // Make the call
            *locus.add_genotype() = genotype;
            
        } else if (copy_budget >= 2 && total(best_support) >= min_total_support_for_call) {
            // The second best allele isn't present or isn't good enough,
            // but the best allele has enough coverage that we can just call
            // two of it.
            
#ifdef debug
            cerr << "Call as best/best" << endl;
#endif
            
            // Say the best is present twice
            genotype.add_allele(best_allele);
            genotype.add_allele(best_allele);
            
            // Get minimum support for filter
            min_site_support = total(best_support);
            
            // Make the call
            *locus.add_genotype() = genotype;

        } else if (copy_budget >= 1 && total(best_support) >= min_total_support_for_call) {
            // We're only supposed to have one copy, and the best allele is good enough to call
            
#ifdef debug
            cerr << "Call as best" << endl;
#endif
            
            // Say the best is present once
            genotype.add_allele(best_allele);
            
            // Get minimum support for filter
            min_site_support = total(best_support);
            
            // Make the call
            *locus.add_genotype() = genotype;
        } else {
            // Either coverage is too low, or we aren't allowed any copies.
            // We can't really call this as anything.
            
#ifdef debug
            cerr << "Do not call" << endl;
#endif
            
            // Don't add the genotype to the locus
        }
    } else {
        // Depth too low. Say we have no idea.
        // TODO: elide variant?
        
        // Don't add the genotype to the locus
    }
    
    // Find the total support for the Locus across all alleles
    Support locus_support;
    for (auto& s : supports) {
        // Sum up all the Supports form all alleles (even the non-best/second-best).
        locus_support += s;
    }
    // Save support
    *locus.mutable_overall_support() = locus_support;

    ////////////////////////////////////////////////////////////////////////////

    // Figure out what child snarls are touched by the paths we have called and
    // how much copy number each should get.
    map<const Snarl*, size_t> child_usage_counts;
    for (size_t i = 0; i < genotype.allele_size(); i++) {
        // For each copy we call as present, find the SnarlTraversal we're
        // asserting
        SnarlTraversal& traversal = here_traversals.at(genotype.allele(i));
        
        for (size_t j = 1; j < traversal.visits_size() - 1; j++) {
            // For each visit to a child snarl
            auto& visit = traversal.visits(j);
            if (visit.node_id() != 0) {
                continue;
            }
        
            // Find the child snarl pointer for the snarl we visit
            const Snarl* child = snarl_manager.manage(visit.snarl());
        
            // Say it's used one more time
            child_usage_counts[child]++;
            
        }
    }

    // Recurse and get traversals for children. We do this for all our children,
    // even the ones called as CN 0, because we need the fully-specified
    // traversals to build our Locus (which needs the alleles we rejected as
    // having no copies).
    unordered_map<const Snarl*, vector<SnarlTraversal>> child_traversals;
    for (const Snarl* child : snarl_manager.children_of(&site)) {
        // Recurse on each child, giving a copy number budget according to the
        // usage count call at this site. This produces fully realized
        // traversals with no Visits to Snarls.
        // Holds ref traversal, best, and optional second best for each child.
        child_traversals[child] = find_best_traversals(augmented, snarl_manager,
            finder, *child, baseline_support, child_usage_counts[child], emit_locus);
    }
    
    for (auto kv : child_traversals) {
        // All children must have at least two traversals (a ref and a best).
        // Off the primary paths, the ref is sort of arbitrary.
        assert(kv.second.size() >= 2);
    }
    
    // Put the best (or ref) traversal for each child in our traversals that
    // visit it (even if that contradicts the calls on the child)
    vector<SnarlTraversal> concrete_traversals;
    for (size_t traversal_number = 0; traversal_number < here_traversals.size(); traversal_number++) {
        // For every abstract traversal of this site, starting with the ref traversal...
        auto& abstract_traversal = here_traversals[traversal_number];
#ifdef debug
        cerr << "Concretizing abstract traversal " << pb2json(abstract_traversal) << endl;
#endif
        
        // Make a "concrete", node-level traversal for every abstract, Snarl-
        // visiting traversal.
        concrete_traversals.emplace_back();
        auto& concrete_traversal = concrete_traversals.back();
        
        for (size_t i = 0; i < abstract_traversal.visits_size(); i++) {
            // Go through all the visits in the abstract traversal
            auto& abstract_visit = abstract_traversal.visits(i);
            
            if (abstract_visit.node_id() != 0) {
                // If they're fully realized, just take them
                *concrete_traversal.add_visits() = abstract_visit;
            } else {
                // If they're visits to children, look up the child
                const Snarl* child = snarl_manager.manage(abstract_visit.snarl());
                
                // Then blit the child's path over. This will be the ref path if
                // we are concrete-izing this snarl's ref traversal, and the
                // best path for the child otherwise. Keep in mind that we may
                // be going through the child backward.
                auto& child_traversal = child_traversals.at(child).at(traversal_number == 0 ? 0 : 1);
                
#ifdef debug
                cerr << "Splicing in child traversal " << pb2json(child_traversal) << endl;
#endif
                
                size_t trav_transfer_start = 0;
                if (i != 0) {
                    // There was a previous visit. It may have been a previous
                    // back-to-back snarl.
                    auto& last_visit = abstract_traversal.visits(i - 1);
                    if (last_visit.node_id() == 0 && to_right_side(last_visit).flip() == to_left_side(abstract_visit)) {
                        // It was indeed a previous back to back site. Don't add the entry node!
#ifdef debug
                        cerr << "Skip entry node for back-to-back sites" << endl;
#endif
                        trav_transfer_start++;
                    }
                }
                for (size_t j = trav_transfer_start; j < child_traversal.visits_size(); j++) {
                    // All the internal visits, in the correct order 
                    *concrete_traversal.add_visits() = abstract_visit.backward() ?
                        reverse(child_traversal.visits(child_traversal.visits_size()- 1 - j)) :
                        child_traversal.visits(j);
                }
            }
        }
#ifdef debug
        cerr << "Finished concrete traversal " << pb2json(concrete_traversals.back()) << endl;
#endif
    }
    
    for (auto& concrete_traversal : concrete_traversals) {
        // Populate the Locus with those traversals by converting to paths
        Path* converted = locus.add_allele();
        
        for (size_t i = 0; i < concrete_traversal.visits_size(); i++) {
            // Convert all the visits to Mappings and stick them in the Locus's Paths
            *converted->add_mapping() = to_mapping(concrete_traversal.visits(i), augmented.graph);
        }
    }
    
    if (locus.genotype_size() > 0) {
        // Emit the locus if we have a call
        emit_locus(locus, &site);
    }
    
    // Build the list of traversals to return as ref, best, second best, with
    // possible repeats.
    vector<SnarlTraversal> to_return{concrete_traversals[0], concrete_traversals[best_allele]};
    if (second_best_allele != -1) {
        to_return.push_back(concrete_traversals[second_best_allele]);
    }
    
    // Return those important traversals
    return to_return;

}

// this was main() in glenn2vcf
void Call2Vcf::call(
    // Augmented graph
    AugmentedGraph& augmented,
    // Should we load a pileup and print out pileup info as comments after
    // variants?
    string pileup_filename) {

    // Toggle support counter
    support_val = use_support_count ? total : support_quality;

    // Set up the graph's paths properly after augmentation modified them.
    augmented.graph.paths.sort_by_mapping_rank();
    augmented.graph.paths.rebuild_mapping_aux();
    
    // Make a list of the specified or autodetected primary reference paths.
    vector<string> primary_path_names = ref_path_names;
    if (primary_path_names.empty()) {
        // Try and guess reference path names for VCF conversion or coverage measurement.
        if (verbose) {
          std:cerr << "Graph has " << augmented.graph.paths.size() << " paths to choose from."
                   << endl;
        }
        if(augmented.graph.paths.size() == 1) {
            // Autodetect the reference path name as the name of the only path
            primary_path_names.push_back((*augmented.graph.paths._paths.begin()).first);
        } else if (augmented.graph.paths.has_path("ref")) {
            // Take any "ref" path.
            primary_path_names.push_back("ref");
        }

        if (verbose && !primary_path_names.empty()) {
            cerr << "Guessed reference path name of " << primary_path_names.front() << endl;
        }
        
    }
    
    // We'll fill this in with a PrimaryPath for every primary reference path
    // that is specified or detected.
    map<string, PrimaryPath> primary_paths;
    for (auto& name : primary_path_names) {
        // Make a PrimaryPath for every primary path we have.
        // Index the primary path and compute the binned supports.   
        primary_paths.emplace(std::piecewise_construct,
            std::forward_as_tuple(name),
            std::forward_as_tuple(augmented, name, ref_bin_size));
    
        auto& primary_path = primary_paths.at(name);
    
        if (verbose) {
            cerr << "Primary path " << name << " average/off-path assumed coverage: "
                << primary_path.get_average_support() << endl;
            cerr << "Mininimum binned average coverage: " << primary_path.get_bin(primary_path.get_min_bin()) << " (bin "
                << (primary_path.get_min_bin() + 1) << " / " << primary_path.get_total_bins() << ")" << endl;
            cerr << "Maxinimum binned average coverage: " << primary_path.get_bin(primary_path.get_max_bin()) << " (bin "
                << (primary_path.get_max_bin() + 1) << " / " << primary_path.get_total_bins() << ")" << endl;
        }
    }
    
    // If applicable, load the pileup.
    // This will hold pileup records by node ID.
    map<int64_t, NodePileup> node_pileups;
    
    function<void(Pileup&)> handle_pileup = [&](Pileup& p) { 
        // Handle each pileup chunk
        for(size_t i = 0; i < p.node_pileups_size(); i++) {
            // Pull out every node pileup
            auto& pileup = p.node_pileups(i);
            // Save the pileup under its node's pointer.
            node_pileups[pileup.node_id()] = pileup;
        }
    };
    if(!pileup_filename.empty()) {
        // We have to load some pileups
        ifstream in;
        in.open(pileup_filename.c_str());
        stream::for_each(in, handle_pileup);
    }
    
    // Parse the translation so we know what original node and offset, if any,
    // each new node came from.
    map<id_t, pos_t> original_positions;
    for (auto& translation : augmented.translations) {
        // TODO: we assume every translation is exactly one old mapping to
        // exactly one new full-node mapping
        auto& new_mapping = translation.to().mapping(0);
        auto& old_mapping = translation.from().mapping(0);
        
        // Store the old source position under the new node ID.
        original_positions[new_mapping.position().node_id()] = make_pos_t(old_mapping.position());
    }
    
    // Make a VCF because we need it in scope later, if we are outputting VCF.
    vcflib::VariantCallFile vcf;
    
    // We also might need to fillin this contig names by path name map
    map<string, string> contig_names_by_path_name;
    
    if (convert_to_vcf) {
        // Do initial setup for VCF output
        
        // Decide on names and lengths for all the primary paths.
        vector<size_t> contig_lengths;
        vector<string> contig_names;
        
        for (size_t i = 0; i < primary_path_names.size(); i++) {
            if (i < contig_name_overrides.size()) {
                // Override this name
                contig_names.push_back(contig_name_overrides.at(i));
            } else {
                // Keep the path name from the graph
                contig_names.push_back(primary_path_names.at(i));
            }
            
            // Allow looking up the assigned contig name later
            contig_names_by_path_name[primary_path_names.at(i)] = contig_names.back();
            
            if (i < length_overrides.size()) {
                // Override this length
                contig_lengths.push_back(length_overrides.at(i));
            } else {
                // Grab the length from the index
                contig_lengths.push_back(primary_paths.at(primary_path_names.at(i)).get_index().sequence.size());
            }
            
            // TODO: is this fall-through-style logic smart, or will we just
            // neglect to warn people that they forgot options by parsing what
            // they said when they provide too few overrides?
        }
    
        // Generate a vcf header. We can't make Variant records without a
        // VariantCallFile, because the variants need to know which of their
        // available info fields or whatever are defined in the file's header,
        // so they know what to output.
        stringstream header_stream;
        write_vcf_header(header_stream, {sample_name}, contig_names, contig_lengths,
            min_mad_for_filter, max_dp_for_filter, max_dp_multiple_for_filter, max_local_dp_multiple_for_filter,
            min_ad_log_likelihood_for_filter);
        
        // Load the headers into a the VCF file object
        string header_string = header_stream.str();
        assert(vcf.openForOutput(header_string));
        
        // Spit out the header
        cout << header_stream.str();
    }
    
    // Find all the top-level sites
    list<const Snarl*> site_queue;
    
    CactusUltrabubbleFinder finder(augmented.graph);
    SnarlManager site_manager = finder.find_snarls();
    
    site_manager.for_each_top_level_snarl_parallel([&](const Snarl* site) {
        // Stick all the sites in this vector.
        #pragma omp critical (sites)
        site_queue.emplace_back(site);
    });
    
    // We're going to run through all the top-level sites and keep just what we
    // can use. If we're converting to VCF it's only stuff on a primary path,
    // and we will break top-level sites to find things on a primary path.
    // Otherwise it's everything.
    vector<const Snarl*> sites;
    
    while(!site_queue.empty()) {
        // Grab the first site
        const Snarl* site = move(site_queue.front());
        site_queue.pop_front();
        
        // If the site is strung on any of the primary paths, find the
        // corresponding PrimaryPath object. Otherwise, leave this null.
        PrimaryPath* primary_path = nullptr;
        {
            auto found = find_path(*site, primary_paths);
            if (found != primary_paths.end()) {
                primary_path = &found->second;
            }
        }
        
        
        if (site->type() == ULTRABUBBLE && primary_path != nullptr) {
            // This site is an ultrabubble on a primary path
        
            // Throw it in the final vector of sites we're going to process.
            sites.push_back(site);
        } else if (site->type() == ULTRABUBBLE && !convert_to_vcf) {
            // This site is an ultrabubble and we can handle things off the
            // primary path.
            
            // Throw it in the final vector of sites we're going to process.
            sites.push_back(site);
            
        } else {
            // The site is not on the primary path or isn't an ultrabubble, but
            // maybe one of its children will meet our requirements.
            
            size_t child_count = site_manager.children_of(site).size();
            
            for(const Snarl* child : site_manager.children_of(site)) {
                // Dump all the children into the queue for separate
                // processing.
                site_queue.emplace_back(child);
            }

            if (verbose) {
                if (child_count) {
                    cerr << "Broke up off-reference site into "
                         << child_count << " children" << endl;
                } else {
                    cerr << "Dropped off-reference site" << endl;
                }
            }
            
        }     
    }

    if (verbose) {
        cerr << "Found " << sites.size() << " sites" << endl;
    }
    
    // Now start looking for traversals of the sites.
    RepresentativeTraversalFinder traversal_finder(augmented, site_manager, max_search_depth, max_search_width,
        max_bubble_paths, [&] (const Snarl& site) -> PathIndex* {
        
        // When the TraversalFinder needs a primary path index for a site, it can look it up with this function.
        auto found = find_path(site, primary_paths);
        if (found != primary_paths.end()) {
            // It's on a path
            return &found->second.get_index();
        } else {
            // It's not on a known primary path, so the TraversalFinder should make its own backbone path
            return nullptr;
        }
    });
    
    // We're going to remember what nodes and edges are covered by sites, so we
    // will know which nodes/edges aren't in any sites and may need generic
    // presence/absence calls.
    set<Node*> covered_nodes;
    set<Edge*> covered_edges;
    
    // When we genotype the sites into Locus objects, we will use this buffer for outputting them.
    vector<Locus> locus_buffer;
    
    // How many sites result in output?
    size_t called_loci = 0;
    
    for(const Snarl* site : sites) {
        // For every site, we're going to make a bunch of Locus objects
        
        // See if the site is on a primary path, so we can use binned support.
        map<string, PrimaryPath>::iterator found_path = find_path(*site, primary_paths);
        
        // We need to figure out how much support a site ought to have.
        // Within its local bin?
        Support baseline_support;
        // On its primary path?
        Support global_baseline_support;
        if (expected_coverage != 0.0) {
            // Use the specified coverage override
            baseline_support.set_forward(expected_coverage / 2);
            baseline_support.set_reverse(expected_coverage / 2);
            global_baseline_support = baseline_support;
        } else if (found_path != primary_paths.end()) {
            // We're on a primary path, so we can find the appropriate bin
        
            // Since the variable part of the site is after the first anchoring node, where does it start?
            // Account for the site possibly being backward on the path.
            size_t variation_start = min(found_path->second.get_index().by_id.at(site->start().node_id()).first
                    + augmented.graph.get_node(site->start().node_id())->sequence().size(),
                found_path->second.get_index().by_id.at(site->end().node_id()).first
                    + augmented.graph.get_node(site->end().node_id())->sequence().size());
            
            // Look in the bins for the primary path to get the support there.
            baseline_support = found_path->second.get_support_at(variation_start);
            
            // And grab the path's overall support
            global_baseline_support = found_path->second.get_average_support();
            
        } else {
            // Just use the primary paths' average support, which may be 0 if there are none.
            // How much support is expected across all the primary paths? May be 0 if there are no primary paths.
            global_baseline_support = PrimaryPath::get_average_support(primary_paths);
            baseline_support = global_baseline_support;
        }
        
        // This function emits the given variant on the given primary path, as
        // VCF. It needs to take the site as an argument because it may be
        // called for children of the site we're working on right now.
        auto emit_variant = [&contig_names_by_path_name, &vcf, &augmented, &original_positions,
            &baseline_support, &global_baseline_support, this](
            const Locus& locus, PrimaryPath& primary_path, const Snarl* site) {
        
            // Note that the locus paths will traverse our site forward, which
            // may make them backward along the primary path.
            bool site_backward = (primary_path.get_index().by_id.at(site->start().node_id()).first >
                primary_path.get_index().by_id.at(site->end().node_id()).first);
        
            // Unpack the genotype back into best and second-best allele
            auto& genotype = locus.genotype(0);
            int best_allele = genotype.allele(0);
            // If we called a single allele, we've lost the second-best allele info. But we won't need it, so we can just say -1.
            int second_best_allele = (genotype.allele_size() >= 2 && genotype.allele(0) != genotype.allele(1)) ?
                genotype.allele(1) :
                -1;
                
            // Populate this with original node IDs, from before augmentation.
            set<id_t> original_nodes;
                
            // Calculate the ID and sequence strings for all the alleles.
            // TODO: we only use some of these
            vector<string> sequences;
            vector<string> id_lists;
            // Also the flags for whether alts are reference (i.e. known)
            vector<bool> is_ref;
            
            for (size_t i = 0; i < locus.allele_size(); i++) {

                // For each allele path in the Locus
                auto& path = locus.allele(i);
#ifdef debug
                cerr << "Extracting allele " << i << ": " << pb2json(path) << endl;
#endif
                // Make a stream for the sequence of the path
                stringstream sequence_stream;
                // And for the description of involved IDs
                stringstream id_stream;
                
                for (size_t j = 0; j < path.mapping_size(); j++) {
                    // For each mapping along the path
                    auto& mapping = path.mapping(j);
                    
                    // Record the sequence
                    string node_sequence = augmented.graph.get_node(mapping.position().node_id())->sequence();
                    if (mapping.position().is_reverse()) {
                        node_sequence = reverse_complement(node_sequence);
                    }
                    sequence_stream << node_sequence;
#ifdef debug
                    cerr << "\tMapping: " << pb2json(mapping) << ", sequence " << node_sequence << endl;
#endif
                    if (j != 0) {
                        // Add a separator
                        id_stream << "_";
                    }
                    // Record the ID
                    id_stream << mapping.position().node_id();
                    
                    if (original_positions.count(mapping.position().node_id())) {
                        // This node is derived from an original graph node. Remember it.
                        original_nodes.insert(id(original_positions[mapping.position().node_id()]));
                    }
                    
                }
                
                // Remember the descriptions of the alleles
                if (site_backward) {
                    sequences.push_back(reverse_complement(sequence_stream.str()));
                } else {
                    sequences.push_back(sequence_stream.str());
                }
#ifdef debug
                cerr << "Recorded allele sequence " << sequences.back() << endl;
#endif
                id_lists.push_back(id_stream.str());
                // And whether they're reference or not
                is_ref.push_back(is_reference(path, augmented));
            }
            
            // Start off declaring the variable part to start at the start of
            // the first anchoring node. We'll clip it back later to just what's
            // after the shared prefix.
            size_t variation_start = min(primary_path.get_index().by_id.at(site->start().node_id()).first,
                primary_path.get_index().by_id.at(site->end().node_id()).first);
        
            // Keep track of the alleles that actually need to go in the VCF:
            // ref, best, and second-best (if any), some of which may overlap.
            // This is the order they will show up in the variant.
            vector<int> used_alleles;
            used_alleles.push_back(0);
            if (best_allele != 0) {
                used_alleles.push_back(best_allele);
            }
            if(second_best_allele != -1 && second_best_allele != 0) {
                used_alleles.push_back(second_best_allele);
            }
            
            // Rewrite the sequences and variation_start to just represent the
            // actually variable part, by dropping any common prefix and common
            // suffix. We just do the whole thing in place, modifying the used
            // entries in sequences.
            
            auto shared_prefix_length = [&](bool backward) {
                size_t shortest_prefix = std::numeric_limits<size_t>::max();
                
                auto here = used_alleles.begin();
                if (here == used_alleles.end()) {
                    // No strings.
                    // Say no prefix is in common...
                    return (size_t) 0;
                }
                auto next = here;
                next++;
                
                if (next == used_alleles.end()) {
                    // Only one string.
                    // Say no prefix is in common...
                    return (size_t) 0;
                }
                
                while (next != used_alleles.end()) {
                    // Consider each allele and the next one after it, as
                    // long as we have both.
                
                    // Figure out the shorter and the longer string
                    string* shorter = &sequences.at(*here);
                    string* longer = &sequences.at(*next);
                    if (shorter->size() > longer->size()) {
                        swap(shorter, longer);
                    }
                
                    // Calculate the match length for this pair
                    size_t match_length;
                    if (backward) {
                        // Find out how far in from the right the first mismatch is.
                        auto mismatch_places = std::mismatch(shorter->rbegin(), shorter->rend(), longer->rbegin());
                        match_length = std::distance(shorter->rbegin(), mismatch_places.first);
                    } else {
                        // Find out how far in from the left the first mismatch is.
                        auto mismatch_places = std::mismatch(shorter->begin(), shorter->end(), longer->begin());
                        match_length = std::distance(shorter->begin(), mismatch_places.first);
                    }
                    
                    // The shared prefix of these strings limits the longest
                    // prefix shared by all strings.
                    shortest_prefix = min(shortest_prefix, match_length);
                
                    here = next;
                    ++next;
                }
                
                // Return the shortest universally shared prefix
                return shortest_prefix;
            };
            // Trim off the shared prefix
            size_t shared_prefix = shared_prefix_length(false);
            for (auto allele : used_alleles) {
                sequences[allele] = sequences[allele].substr(shared_prefix);
            }
            // Add it onto the start coordinate
            variation_start += shared_prefix;
            
            // Then find and trim off the shared suffix
            size_t shared_suffix = shared_prefix_length(true);
            for (auto allele : used_alleles) {
                sequences[allele] = sequences[allele].substr(0, sequences[allele].size() - shared_suffix);
            }
            
            // Make a Variant
            vcflib::Variant variant;
            variant.sequenceName = contig_names_by_path_name.at(primary_path.get_name());
            variant.setVariantCallFile(vcf);
            variant.quality = 0;
            // Position should be 1-based and offset with our offset option.
            variant.position = variation_start + 1 + variant_offset;
            
            // Set the ID based on the IDs of the involved nodes. Note that the best
            // allele may have no nodes (because it's a pure edge)
            variant.id = id_lists.at(best_allele);
            if(second_best_allele != -1 && !id_lists.at(second_best_allele).empty()) {
                // Add the second best allele's nodes in.
                variant.id += "-" + id_lists.at(second_best_allele);
            }
            
            
            if(sequences.at(0).empty() ||
                (best_allele != -1 && sequences.at(best_allele).empty()) ||
                (second_best_allele != -1 && sequences.at(second_best_allele).empty())) {
                
                // Fix up the case where we have an empty allele.
                
                // We need to grab the character before the variable part of the
                // site in the reference.
                assert(variation_start > 0);
                string extra_base = char_to_string(primary_path.get_index().sequence.at(variation_start - 1));
                
                for(auto& seq : sequences) {
                    // Stick it on the front of all the allele sequences
                    seq = extra_base + seq;
                }
                
                // Budge the variant left
                variant.position--;
            }
            
            // Make sure the ref allele is correct
            {
                string real_ref = primary_path.get_index().sequence.substr(
                    variant.position - variant_offset - 1, sequences.front().size());
                string got_ref = sequences.front();
                
                if (real_ref != got_ref) {
                    cerr << "Error: Ref should be " << real_ref << " but is " << got_ref << " at " << variant.position << endl;
                    throw runtime_error("Reference mismatch at site " + pb2json(*site));
                }
            
            }
            
            // Add the ref allele to the variant
            create_ref_allele(variant, sequences.front());
            
            // Add the best allele
            assert(best_allele != -1);
            int best_alt = add_alt_allele(variant, sequences.at(best_allele));
            
            int second_best_alt = (second_best_allele == -1) ? -1 : add_alt_allele(variant, sequences.at(second_best_allele));

            // Say we're going to spit out the genotype for this sample.        
            variant.format.push_back("GT");
            auto& genotype_vector = variant.samples[sample_name]["GT"];

            if (locus.genotype_size() > 0) {
                // We actually made a call. Emit the first genotype, which is the call.
                
                // We need to rewrite the allele numbers to alt numbers, since
                // we aren't keeping all the alleles in the VCF, so we can't use
                // the natural conversion of Genotype to VCF genotype string.
                
                // Emit parts into this stream
                stringstream stream;
                for (size_t i = 0; i < genotype.allele_size(); i++) {
                    // For each allele called as present in the genotype
                    
                    // Convert from allele number to alt number
                    if (genotype.allele(i) == best_allele) {
                        stream << best_alt;
                    } else if (genotype.allele(i) == second_best_allele) {
                        stream << second_best_alt;
                    } else {
                        throw runtime_error("Allele " + to_string(genotype.allele(i)) +
                            " is not best or second-best and has no alt");
                    }
                    
                    if (i + 1 != genotype.allele_size()) {
                        // Write a separator after all but the last one
                        stream << (genotype.is_phased() ? '|' : '/');
                    }
                }
                // Save the finished genotype
                genotype_vector.push_back(stream.str());              
            } else {
                // Say there's no call here
                genotype_vector.push_back("./.");
            }
            
            // Now fill in all the other variant info/format stuff

            if((best_allele != 0 && is_ref.at(best_allele)) || 
                (second_best_allele != 0 && second_best_allele != -1 && is_ref.at(second_best_allele))) {
                // Flag the variant as reference if either of its two best alleles
                // is known but not the primary path. Don't put in a false entry if
                // it isn't known, because vcflib will spit out the flag anyway...
                variant.infoFlags["XREF"] = true;
            }
            
            for (auto id : original_nodes) {
                // Add references to the relevant original nodes
                variant.info["XSEE"].push_back(to_string(id));
            }
            
            for (size_t i = 1; i < variant.alleles.size(); i++) {
                // Claculate the SVLEN for this non-reference allele
                int64_t svlen = (int64_t) variant.alleles.at(i).size() - (int64_t) variant.alleles.at(0).size();
                
                // Add it in
                variant.info["SVLEN"].push_back(to_string(svlen));
            }
            
            // Set up the depth format field
            variant.format.push_back("DP");
            // And expected depth
            variant.format.push_back("XDP");
            // And allelic depth
            variant.format.push_back("AD");
            // And the log likelihood from the assignment of reads among the
            // present alleles
            variant.format.push_back("XADL");
            // And strand bias
            variant.format.push_back("SB");
            // Also the alt allele depth
            variant.format.push_back("XAAD");
            
            // Compute the total support for all the alts that will be appearing
            Support total_support;
            // And total alt allele depth for the alt alleles
            Support alt_support;

            for (int allele : used_alleles) {
                // For all the alleles we are using, look at the support.
                auto& support = locus.support(allele);
                
                // Set up allele-specific stats for the allele
                variant.samples[sample_name]["AD"].push_back(to_string((int64_t)round(total(support))));
                variant.samples[sample_name]["SB"].push_back(to_string((int64_t)round(support.forward())));
                variant.samples[sample_name]["SB"].push_back(to_string((int64_t)round(support.reverse())));
                
                // Sum up into total depth
                total_support += support;
                
                if (allele != 0) {
                    // It's not the primary reference allele
                    alt_support += support;
                }
            }

            // Find the min total support of anything called
            double min_site_support = INFINITY;
            double min_site_quality = INFINITY;
            
            for (size_t i = 0; i < genotype.allele_size(); i++) {
                // Min all the total supports from the non-ref alleles called as present
                min_site_support = min(min_site_support, total(locus.support(genotype.allele(i))));
                min_site_quality = min(min_site_quality, locus.support(genotype.allele(i)).quality());
            }
            
            // Find the binomial bias between the called alleles, if multiple were called.
            double ad_log_likelihood = INFINITY;
            if (second_best_allele != -1) {
                // How many of the less common one do we have?
                size_t successes = round(total(locus.support(second_best_allele)));
                // Out of how many chances
                size_t trials = successes + (size_t) round(total(locus.support(best_allele)));
                
                assert(trials >= successes);
                
                // How weird is that?                
                ad_log_likelihood = binomial_cmf_ln(prob_to_logprob((real_t) 0.5), trials, successes);
                
                assert(!std::isnan(ad_log_likelihood));
                
                variant.samples[sample_name]["XADL"].push_back(to_string(ad_log_likelihood));
            } else {
                // No need to assign reads between two alleles
                variant.samples[sample_name]["XADL"].push_back(".");
            }

            // Set the variant's total depth            
            string depth_string = to_string((int64_t)round(total(total_support)));
            variant.info["DP"].push_back(depth_string); // We only have one sample, so variant depth = sample depth
            
            // And for the sample
            variant.samples[sample_name]["DP"].push_back(depth_string);
            
            // Set the sample's local and global expected depth            
            variant.samples[sample_name]["XDP"].push_back(to_string((int64_t)round(total(baseline_support))));
            variant.samples[sample_name]["XDP"].push_back(to_string((int64_t)round(total(global_baseline_support))));
            
            // And its depth of non-0 alleles
            variant.samples[sample_name]["XAAD"].push_back(to_string((int64_t)round(total(alt_support))));

            // Set the total support quality of the min allele as the variant quality
            variant.quality = min_site_quality;

            // Now do the filters
            variant.filter = "PASS";            
            if (min_site_support < min_mad_for_filter) {
                // Apply Min Allele Depth cutoff across all alleles (even ref)
                variant.filter = "lowad";
            } else if (max_dp_for_filter != 0 && total(total_support) > max_dp_for_filter) {
                // Apply the max depth cutoff
                variant.filter = "highabsdp";
            } else if (max_dp_multiple_for_filter != 0 &&
                total(total_support) > max_dp_multiple_for_filter * total(global_baseline_support)) {
                // Apply the max depth multiple cutoff
                // TODO: Different standard for sites called as haploid
                variant.filter = "highreldp";
            } else if (max_local_dp_multiple_for_filter != 0 &&
                total(total_support) > max_local_dp_multiple_for_filter * total(baseline_support)) {
                // Apply the max local depth multiple cutoff
                // TODO: Different standard for sites called as haoploid
                variant.filter = "highlocaldp";
            } else if (min_ad_log_likelihood_for_filter != 0 &&
                ad_log_likelihood < min_ad_log_likelihood_for_filter) {
                // We have a het, but the assignment of reads between the two branches is just too weird
                variant.filter = "lowxadl";
            }
            
            // Don't bother with trivial calls
            if (write_trivial_calls ||
                (genotype_vector.back() != "./." && genotype_vector.back() != ".|." &&
                 genotype_vector.back() != "0/0" && genotype_vector.back() != "0|0")) {
            
                if(can_write_alleles(variant)) {
                    // No need to check for collisions because we assume sites are correctly found.
                    // Output the created VCF variant.
                    cout << variant << endl;
            
                } else {
                    if (verbose) {
                        cerr << "Variant is too large" << endl;
                    }
                    // TODO: track bases lost again
                }
            }
        };
        
        // Recursively type the site, using that support and an assumption of a diploid sample.
        find_best_traversals(augmented, site_manager, &traversal_finder, *site, baseline_support, 2,
            [&locus_buffer, &emit_variant, &site_manager, &called_loci, &primary_paths, &augmented,
            &covered_nodes, &covered_edges, this](const Locus& locus, const Snarl* site) {
            
            // Now we have the Locus with call information, and the site (either
            // the root snarl we passed in or a child snarl) that the call is
            // for. We need to output the call.
        
            if (convert_to_vcf) {
                // We want to emit VCF
                
                // Look up the path this child site lives on. (TODO: just capture and use the path the parent lives on?)
                auto found_path = find_path(*site, primary_paths);
                if(found_path != primary_paths.end()) {
                    // And this site is on a primary path
                    
                    // Emit the variant for this Locus
                    emit_variant(locus, found_path->second, site);
                }
                // Otherwise discard it as off-path
                // TODO: update bases lost
            } else {
                // Emit the locus itself
                locus_buffer.push_back(locus);
                stream::write_buffered(cout, locus_buffer, locus_buffer_size);
            }
            
            // We called a site
            called_loci++;
            
            // Mark all the nodes and edges in the site as covered
            auto contents = site_manager.deep_contents(site, augmented.graph, true);
            for (auto* node : contents.first) {
                covered_nodes.insert(node);
            }
            for (auto* edge : contents.second) {
                covered_edges.insert(edge);
            }
        });
    }
    
    if (verbose) {
        cerr << "Called " << called_loci << " loci" << endl;
    }
    
    // OK now we have handled all the real sites. But there are still nodes and
    // edges that we might want to call as present or absent.
    
    if (!convert_to_vcf) {
        
        size_t extra_loci = 0;
        
        if (call_other_by_coverage) {
            // We should look at the coverage of things off the primary path and
            // make calls on them.
        
            augmented.graph.for_each_edge([&](Edge* e) {
                // We want to make calls on all the edges that aren't covered yet
                if (covered_edges.count(e)) {
                    // Skip this edge
                    return;
                }
                
                // Make a couple of fake Visits
                Visit from_visit;
                from_visit.set_node_id(e->from());
                from_visit.set_backward(e->from_start());
                Visit to_visit;
                to_visit.set_node_id(e->to());
                to_visit.set_backward(e->to_end());
                
                // Make a Locus for the edge
                Locus locus;
                
                // Give it an allele
                Path* path = locus.add_allele();
                
                // Fill in 
                *path->add_mapping() = to_mapping(from_visit, augmented.graph);
                *path->add_mapping() = to_mapping(to_visit, augmented.graph);
                
                // Set the support
                *locus.add_support() = augmented.edge_supports[e];
                *locus.mutable_overall_support() = augmented.edge_supports[e];
                
                // Decide on the genotype
                Genotype gt;
                
                // TODO: use the coverage bins
                if (support_val(locus.support(0)) > support_val(PrimaryPath::get_average_support(primary_paths)) * 0.25) {
                    // We're closer to 1 copy than 0 copies
                    gt.add_allele(0);
                    
                    if (support_val(locus.support(0)) > support_val(PrimaryPath::get_average_support(primary_paths)) * 0.75) {
                        // We're closer to 2 copies than 1 copy
                        gt.add_allele(0);
                    }
                }
                // Save the genotype with 0, 1, or 2 copies.
                *locus.add_genotype() = gt;
                
                // Send out the locus
                locus_buffer.push_back(locus);
                stream::write_buffered(cout, locus_buffer, locus_buffer_size);
                
                extra_loci++;
                
                // TODO: look at average node coverages and do node loci (in
                // case any nodes have no edges?)
                
            });
        } else {
            // We should just assert the existence of the primary path edges
            // that weren't in something. Everything not asserted will get
            // subsetted out.
            
            for (auto& kv : primary_paths) {
                // For every primary path in the graph
                auto& primary_path = kv.second;
                
                // Remember the end of the previous ndoe
                NodeSide previous_end;
                
                for (auto& offset_and_side : primary_path.get_index()) {
                    // For every NodeSide that happens in the primary path
                    
                    // Grab the side
                    NodeSide here = offset_and_side.second;
                    
                    if (previous_end == NodeSide()) {
                        // Skip the first node and remember its end
                        previous_end = here.flip();
                        continue;
                    }
                    
                    // Find the edge we crossed
                    Edge* crossed = augmented.graph.get_edge(previous_end, here);
                    assert(crossed != nullptr);
                    
                    if (covered_edges.count(crossed)) {
                        // If the edge we crossed is covered by a snarl, don't
                        // emit anything.
                        previous_end = here.flip();
                        continue;
                    }
                    
                    // If the edge we're crossing isn't covered, we should
                    // assert the primary path here.
                    
                    // Make a couple of fake Visits
                    Visit from_visit;
                    from_visit.set_node_id(crossed->from());
                    from_visit.set_backward(crossed->from_start());
                    Visit to_visit;
                    to_visit.set_node_id(crossed->to());
                    to_visit.set_backward(crossed->to_end());
                    
                    // Make a Locus for the edge
                    Locus locus;
                    
                    // Give it an allele
                    Path* path = locus.add_allele();
                    
                    // Fill in 
                    *path->add_mapping() = to_mapping(from_visit, augmented.graph);
                    *path->add_mapping() = to_mapping(to_visit, augmented.graph);
                    
                    // Set the support
                    *locus.add_support() = augmented.get_support(crossed);
                    *locus.mutable_overall_support() = augmented.get_support(crossed);
                    
                    // Decide on the genotype of hom ref.
                    Genotype* gt = locus.add_genotype();
                    gt->add_allele(0);
                    gt->add_allele(0);
                    
                    // Send out the locus
                    locus_buffer.push_back(locus);
                    stream::write_buffered(cout, locus_buffer, locus_buffer_size);
                    
                    extra_loci++;
                    
                    // Make sure to remember the end of the node we just did,
                    // for looking at the next node.
                    previous_end = here.flip();
                    
                }
            }
        }
        
        // Flush the buffer of Locus objects we have to write
        stream::write_buffered(cout, locus_buffer, 0);
        
        if (verbose) {
            cerr << "Called " << extra_loci << " extra loci with copy number estimates" << endl;
        }
        
    }
    
}

bool Call2Vcf::is_reference(const SnarlTraversal& trav, AugmentedGraph& augmented) {
    
    // Keep track of the previous NodeSide
    NodeSide previous;
    
    // We'll call this function with each visit in turn.
    // If it ever returns false, the whole thing is nonreference.
    auto experience_visit = [&](const Visit& visit) {
        // TODO: handle nested sites
        assert(visit.node_id());
        
        if (previous.node != 0) {
            // Consider the edge from the previous visit
            Edge* edge = augmented.graph.get_edge(previous, to_left_side(visit));
            
            if (augmented.edge_calls.at(edge) != CALL_REFERENCE) {
                // Found a novel edge!
                return false;
            }
        }
        
        if (augmented.node_calls.at(augmented.graph.get_node(visit.node_id())) != CALL_REFERENCE) {
            // This node itself is novel
            return false;
        }
        
        // Remember we want an edge from this visit when we look at the next
        // one.
        previous = to_right_side(visit);
        
        // This visit is known.
        return true;
    };
    
    // Experience the entire traversal from start to end
    for (size_t i = 0; i < trav.visits_size(); i++) {
        if (!experience_visit(trav.visits(i))) {
            return false;
        }
    }
    
    // And if we make it through it's a reference traversal.
    return true;
        
}

bool Call2Vcf::is_reference(const Path& path, AugmentedGraph& augmented) {
    
    // The path can't be empty because it's not clear if an empty path should be
    // reference or not.
    assert(path.mapping_size() != 0);
    
    for (size_t i = 0; i < path.mapping_size(); i++) {
        // Check each mapping
        auto& mapping = path.mapping(i);
        
        if (augmented.get_call(augmented.graph.get_node(mapping.position().node_id())) != CALL_REFERENCE) {
            // We use a novel node
            return false;
        }
        
        if (i + 1 < path.mapping_size()) {
            // Also look at the next mapping
            auto& next_mapping = path.mapping(i + 1);
            
            // And see about the edge to it
            Edge* edge = augmented.graph.get_edge(to_right_side(to_visit(mapping)), to_left_side(to_visit(next_mapping)));
            if (augmented.get_call(edge) != CALL_REFERENCE) {
                // We used a novel edge
                return false;
            }
        }
    }
    
    // If we get through everything it's reference.
    return true;
}

}
