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

// We represent support as a pair, but we define math for it.
// We use doubles because we may need fractional math.
// We can't use Support because it's int-based.
typedef pair<double, double> FractionalSupport;

/**
 * Add two FractionalSupport values together, accounting for strand.
 */
FractionalSupport operator+(const FractionalSupport& one, const FractionalSupport& other) {
    return make_pair(one.first + other.first, one.second + other.second);
}

/**
 * Add in a FractionalSupport to another.
 */
FractionalSupport& operator+=(FractionalSupport& one, const FractionalSupport& other) {
    one.first += other.first;
    one.second += other.second;
    return one;
}


/**
 * Scale a support by a factor.
 */
template<typename Scalar>
FractionalSupport operator*(const FractionalSupport& support, const Scalar& scale) {
    return make_pair(support.first * scale, support.second * scale);
}

/**
 * Scale a support by a factor, the other way
 */
template<typename Scalar>
FractionalSupport operator*(const Scalar& scale, const FractionalSupport& support) {
    return make_pair(support.first * scale, support.second * scale);
}

/**
 * Divide a support by a factor.
 */
template<typename Scalar>
FractionalSupport operator/(const FractionalSupport& support, const Scalar& scale) {
    return make_pair(support.first / scale, support.second / scale);
}

/**
 * Upgrade an integral Protobuf Support to a FractionalSupport.
 */
FractionalSupport from_support(const Support& other) {
    return make_pair(other.forward(), other.reverse());
}

/**
 * Downgrade a FractionalSupport to an integral Protobuf suopport.
 */
Support to_support(const FractionalSupport& other) {
    Support to_return;
    to_return.set_forward(other.first);
    to_return.set_reverse(other.second);
    return to_return;
}
    
/**
 * Allow printing a FractionalSupport.
 */
ostream& operator<<(ostream& stream, const FractionalSupport& support) {
    return stream << support.first << "," << support.second;
}

/**
 * Get the total read support in a FractionalSupport.
 */
double total(const FractionalSupport& support) {
    return support.first + support.second;
}

/**
 * Get the strand bias of a FractionalSupport.
 */
double strand_bias(const FractionalSupport& support) {
    return max(support.first, support.second) / (support.first + support.second);
}

/**
 * Get the minimum support of a pair of FractionalSupports, by taking the min in each
 * orientation.
 */
FractionalSupport support_min(const FractionalSupport& a, const FractionalSupport& b) {
    return make_pair(min(a.first, b.first), min(a.second, b.second));
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
 * Write a minimal VCF header for a single-sample file.
 */
void write_vcf_header(ostream& stream, string& sample_name, string& contig_name, size_t contig_size,
                      int min_mad_for_filter) {
    stream << "##fileformat=VCFv4.2" << endl;
    stream << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" << endl;
    stream << "##INFO=<ID=XREF,Number=0,Type=Flag,Description=\"Present in original graph\">" << endl;
    stream << "##INFO=<ID=XSEE,Number=.,Type=String,Description=\"Original graph node:offset cross-references\">" << endl;
    stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << endl;
    stream << "##FILTER=<ID=FAIL,Description=\"Variant does not meet minimum allele read support threshold of " << min_mad_for_filter << "\">" <<endl;
    stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    stream << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << endl;
    stream << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Forward and reverse support for ref and alt alleles.\">" << endl;
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    stream << "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">" << endl;
    stream << "##FORMAT=<ID=AL,Number=.,Type=Float,Description=\"Allelic likelihoods for the ref and alt alleles in the order listed\">" << endl;
    
    if(!contig_name.empty()) {
        // Announce the contig as well.
        stream << "##contig=<ID=" << contig_name << ",length=" << contig_size << ">" << endl;
    }
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << endl;
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
string get_pileup_line(const map<int64_t, NodePileup>& nodePileups,
    const set<pair<int64_t, size_t>>& refCrossreferences,
    const set<pair<int64_t, size_t>>& altCrossreferences) {
    // We'll make a stringstream to write to.
    stringstream out;
    
    out << "#";
    
    for(const auto& xref : refCrossreferences) {
        // For every cross-reference
        if(nodePileups.count(xref.first) && nodePileups.at(xref.first).base_pileup_size() > xref.second) {
            // If we have that base pileup, grab it
            auto basePileup = nodePileups.at(xref.first).base_pileup(xref.second);
            
            out << xref.first << ":" << xref.second << " (ref) " << basePileup.bases() << "\t";
        }
        // Nodes with no pileups (either no pileups were provided or they didn't
        // appear/weren't visited by reads) will not be mentioned on this line
    }
    
    for(const auto& xref : altCrossreferences) {
        // For every cross-reference
        if(nodePileups.count(xref.first) && nodePileups.at(xref.first).base_pileup_size() > xref.second) {
            // If we have that base pileup, grab it
            auto basePileup = nodePileups.at(xref.first).base_pileup(xref.second);
            
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

// this was main() in glenn2vcf
void Call2Vcf::call(
    // Augmented graph
    AugmentedGraph& augmented,
    // Should we load a pileup and print out pileup info as comments after
    // variants?
    string pileupFilename) {
    
    // Set up the graph's paths properly after augmentation modified them.
    augmented.graph.paths.sort_by_mapping_rank();
    augmented.graph.paths.rebuild_mapping_aux();
    
    if(refPathName.empty()) {
        if (verbose) {
          std:cerr << "Graph has " << augmented.graph.paths.size() << " paths to choose from."
                   << endl;
        }
        if(augmented.graph.paths.size() == 1) {
            // Autodetect the reference path name as the name of the only path
            refPathName = (*augmented.graph.paths._paths.begin()).first;
        } else {
            refPathName = "ref";
        }

        if (verbose) {
            cerr << "Guessed reference path name of " << refPathName
                      << endl;
        }
    }
    
    // Follow the reference path and extract indexes we need: index by node ID,
    // index by node start, and the reconstructed path sequence.
    PathIndex index(augmented.graph, refPathName, true);

    // Store support binned along reference path;
    // Last bin extended to include remainder
    refBinSize = min(refBinSize, index.sequence.size());
    vector<FractionalSupport> binnedSupport(max(1, int(index.sequence.size() / refBinSize)),
                                  FractionalSupport(expCoverage / 2, expCoverage /2));
    
    // Crunch the numbers on the reference and its read support. How much read
    // support in total (node length * aligned reads) does the primary path get?
    FractionalSupport primaryPathTotalSupport = FractionalSupport();
    for(auto& pointerAndSupport : augmented.node_supports) {
        if(index.by_id.count(pointerAndSupport.first->id())) {
            // This is a primary path node. Add in the total read bases supporting it
            primaryPathTotalSupport += pointerAndSupport.first->sequence().size() * from_support(pointerAndSupport.second);
            
            // We also update the total for the appropriate bin
            if (expCoverage == 0) {
                int bin = index.by_id[pointerAndSupport.first->id()].first / refBinSize;
                if (bin == binnedSupport.size()) {
                    --bin;
                }
                binnedSupport[bin] = binnedSupport[bin] + 
                    pointerAndSupport.first->sequence().size() * from_support(pointerAndSupport.second);
            }
        }
    }
    // Calculate average support in reads per base
    auto primaryPathAverageSupport = primaryPathTotalSupport / index.sequence.size();
    
    // Average out the support bins too (in place)
    int minBin = -1;
    int maxBin = -1;
    for (int i = 0; i < binnedSupport.size(); ++i) {
        if (expCoverage == 0) {
            binnedSupport[i] = binnedSupport[i] / (
                i < binnedSupport.size() - 1 ? (double)refBinSize :
                (double)(refBinSize + index.sequence.size() % refBinSize));
        }
        if (minBin == -1 || binnedSupport[i] < binnedSupport[minBin]) {
            minBin = i;
        }
        if (maxBin == -1 || binnedSupport[i] > binnedSupport[maxBin]) {
            maxBin = i;
        }
    }

    if (verbose) {
        cerr << "Primary path average coverage: " << primaryPathAverageSupport << endl;
        cerr << "Mininimum binned average coverage: " << binnedSupport[minBin] << " (bin "
                  << (minBin + 1) << " / " << binnedSupport.size() << ")" << endl;
        cerr << "Maxinimum binned average coverage: " << binnedSupport[maxBin] << " (bin "
                  << (maxBin + 1) << " / " << binnedSupport.size() << ")" << endl;
    }
    
    // If applicable, load the pileup.
    // This will hold pileup records by node ID.
    map<int64_t, NodePileup> nodePileups;
    
    function<void(Pileup&)> handlePileup = [&](Pileup& p) { 
        // Handle each pileup chunk
        for(size_t i = 0; i < p.node_pileups_size(); i++) {
            // Pull out every node pileup
            auto& pileup = p.node_pileups(i);
            // Save the pileup under its node's pointer.
            nodePileups[pileup.node_id()] = pileup;
        }
    };
    if(!pileupFilename.empty()) {
        // We have to load some pileups
        ifstream in;
        in.open(pileupFilename.c_str());
        stream::for_each(in, handlePileup);
    }
    
    // Make a VCF because we need it in scope later, if we are outputting VCF.
    vcflib::VariantCallFile vcf;
    
    if (convert_to_vcf) {
        // Generate a vcf header. We can't make Variant records without a
        // VariantCallFile, because the variants need to know which of their
        // available info fields or whatever are defined in the file's header, so
        // they know what to output.
        // Handle length override if specified.
        stringstream headerStream;
        write_vcf_header(headerStream, sampleName, contigName,
                         lengthOverride != -1 ? lengthOverride : (index.sequence.size() + variantOffset),
                         min_mad_for_filter);
        
        // Load the headers into a the VCF file object
        string headerString = headerStream.str();
        assert(vcf.openForOutput(headerString));
        
        // Spit out the header
        cout << headerStream.str();
    }
    
    // Then go through it from the graph's point of view: first over alt nodes
    // backending into the reference (creating things occupying ranges to which
    // we can attribute copy number) and then over reference nodes.

    // We need to track the bases lost.
    size_t basesLost = 0;
    
    
    // Do the new thing where we support multiple alleles

    // Find all the top-level sites
    list<const Snarl*> site_queue;
    
    CactusUltrabubbleFinder finder(augmented.graph, refPathName);
    SnarlManager site_manager = finder.find_snarls();
    
    site_manager.for_each_top_level_snarl_parallel([&](const Snarl* site) {
        if(!index.by_id.count(site->start().node_id()) || !index.by_id.count(site->end().node_id())) {
            // Skip top-level sites that aren't on the reference path at both ends.
            return;
        }
    
        // Stick all the sites in this vector.
        #pragma omp critical (sites)
        site_queue.emplace_back(site);
    });
    
    // We're going to run through all the top-level sites and break up the
    // huge ones and throw out the ones not on the reference, leaving only
    // manageably-sized sites on the ref path. We won't be able to genotype
    // huge translocations or deletions, but we can save those for the
    // proper nested-site-aware genotyper.
    vector<const Snarl*> sites;
    
    while(!site_queue.empty()) {
        // Grab the first site
        const Snarl* site = move(site_queue.front());
        site_queue.pop_front();
        
        if (!index.by_id.count(site->start().node_id())) {
            #pragma omp critical (cerr)
            cerr << "Off-reference left end: " << pb2json(*site) << endl;
        }
        if (!index.by_id.count(site->end().node_id())) {
            #pragma omp critical (cerr)
            cerr << "Off-reference right end: " << pb2json(*site) << endl;
        }
        
        // Sites should not have gotten into the queue without being on the primary path at both ends
        assert(index.by_id.count(site->start().node_id()) && index.by_id.count(site->end().node_id()));
        
        // Where does the variable region start and end on the reference?
        size_t ref_start = index.by_id.at(site->start().node_id()).first +
            augmented.graph.get_node(site->start().node_id())->sequence().size();
        size_t ref_end = index.by_id.at(site->end().node_id()).first;
        
        if(ref_start > ref_end) {
            // Make sure it's the right way around (it will get set straight
            // in the site when we do the bubbling-up).
            swap(ref_start, ref_end);
        }
        
        if(!site_manager.children_of(site).empty() && ref_end > ref_start + max_ref_length) {
            // Site is too big and has children. Split it up and do the
            // children.
            for(const Snarl* child : site_manager.children_of(site)) {
                // Dump all the children into the queue for separate
                // processing.
                
                if(!index.by_id.count(child->start().node_id()) || !index.by_id.count(child->end().node_id())) {
                    // Skip child sites that aren't on the reference path at both ends.
                    continue;
                }
                
                site_queue.emplace_back(child);
            }

            if (verbose) {
                cerr << "Broke up site from " << ref_start << " to " << ref_end << " into "
                          << site_manager.children_of(site).size() << " children" << endl;
            }
            
        } else {
            // With no children, site may still be huge, but it doesn't
            // matter because there's nothing to nest in it, so we can
            // genotype it just fine.
            
            // With children, it's practical to just include the child
            // genotypes in the site genotype.
            
            if(ref_end > ref_start + max_ref_length) {
                // This site is big but we left it anyway.
                if (verbose) {
                    cerr << "Left site from " << ref_start << " to " << ref_end << " with "
                              << site_manager.children_of(site).size() << " children" << endl;
                }
            }
            
            // Make sure start and end are front-ways relative to the ref path.
            if(index.by_id.at(site->start().node_id()).first > index.by_id.at(site->end().node_id()).first) {
                // The site's end happens before its start, flip it around
                site_manager.flip(site);
            }
            
            // Throw it in the final vector of sites we're going to process.
            sites.push_back(site);
        }                
    }

    if (verbose) {
        cerr << "Found " << sites.size() << " sites" << endl;
    }
    
    // Now start looking for traversals of the sites
    RepresentativeTraversalFinder traversal_finder(augmented, site_manager, index, maxDepth, max_bubble_paths);
    
    // When we genotype the sites into Locus objects, we will use this buffer for outputting them.
    vector<Locus> locus_buffer;
    
    for(const Snarl* site : sites) {
        // For every site, we're going to make a variant
        
        // Get the traversals. The ref traversal is the first one. They are all
        // in site orientation, which may oppose primary path orientation.
        auto traversals = traversal_finder.find_traversals(*site);
        
        // Make a Locus to represent this site, with all the Paths of the
        // different alleles.
        Locus locus;
        
        // And keep around a vector of is_reference statuses for all the alleles
        vector<bool> is_ref;
        
        // We turn them all back into NodeTraversal vectors
        vector<pair<vector<NodeTraversal>, bool>> ordered_paths;
        
        for (auto& traversal : traversals) {
            // Convert each traversal to a path
            Path* path = locus.add_allele();
        
            // Start at the start of the site
            *path->add_mapping() = to_mapping(site->start(), augmented.graph);
            for (size_t i = 0; i < traversal.visits_size(); i++) {
                // Then visit each node in turn
                *path->add_mapping() = to_mapping(traversal.visits(i), augmented.graph);
            }
            // Finish up with the site's end
            *path->add_mapping() = to_mapping(site->end(), augmented.graph);
            
            // Save the traversal's reference status
            is_ref.push_back(is_reference(traversal, augmented, index));
            
            // TODO: Also, save the data the old way that powers the VCF output
            // currently
            
            // Make each traversal into a vector of NodeTraversals
            vector<NodeTraversal> traversal_vector;
            // Start at the start of the site
            traversal_vector.push_back(to_node_traversal(site->start(), augmented.graph));
            for (size_t i = 0; i < traversal.visits_size(); i++) {
                // Then visit each node in turn
                traversal_vector.push_back(to_node_traversal(traversal.visits(i), augmented.graph));
            }
            // Finish up with the site's end
            traversal_vector.push_back(to_node_traversal(site->end(), augmented.graph));
            
            if (index.by_id.at(site->start().node_id()).first > index.by_id.at(site->end().node_id()).first) {
                // This site is backward relative to the primary path.
                
                // Turn the traversal around to primary path orientation.
                reverse(traversal_vector.begin(), traversal_vector.end());
                for (auto& trav : traversal_vector) {
                    // Flip every traversal in it
                    trav = trav.reverse();
                }
            }
            
            // Put it in the list, with an annotation saying whether it's a
            // previously-known traversal or not.
            ordered_paths.push_back(make_pair(traversal_vector, is_ref.back()));
        }
        
        // Collect sequences for all the paths
        vector<string> sequences;
        // And the lists of involved IDs that we use for variant IDs
        vector<string> id_lists;
        // Calculate average and min support for all the alts. These both need
        // to be the same type, so they are interchangeable in a reference and
        // we can pick one to use.
        vector<FractionalSupport> min_supports;
        vector<FractionalSupport> average_supports;
        // And the min likelihood along each path
        vector<double> min_likelihoods;
        for(size_t i = 0; i < locus.allele_size(); i++) {
            // Go through all the Paths in the Locus
            const Path& path = locus.allele(i);
            
            // We use this to construct the allele sequence
            stringstream sequence_stream;
            
            // And this to compose the allele's name in terms of node IDs
            stringstream id_stream;
            
            // What's the total support for this path?
            Support total_support;
            
            // And the min
            Support min_support;
            
            // Also, what's the min likelihood
            double min_likelihood = INFINITY;
                            
            for(int64_t i = 1; i + 1 < path.mapping_size(); i++) {
                // For all but the first and last nodes (which are anchors)...
                
                // Grab the mapping
                const Mapping& here = path.mapping(i);
                
                // Look up the node ID we're visiting
                id_t node_id = path.mapping(i).position().node_id();
                // And the Node* in the graph
                Node* node = augmented.graph.get_node(node_id);
                
                // Grab the node sequence in the correct orientation.
                string added_sequence = node->sequence();
                if(path.mapping(i).position().is_reverse()) {
                    // If the node is traversed backward, we need to flip its sequence.
                    added_sequence = reverse_complement(added_sequence);
                }
                
                // Stick the sequence
                sequence_stream << added_sequence;
                
                // Record ID
                id_stream << to_string(node_id);
                if(i + 2 != path.mapping_size()) {
                    // If this is not the last mapping we're going to do, add a separator.
                    id_stream << "_";
                }
                
                // Peek at the next and previous Mappings
                const Mapping& prev = path.mapping(i - 1);
                const Mapping& next = path.mapping(i + 1);
                
                // How much support do we have for visiting this node?
                Support node_support = augmented.node_supports.at(node);
                // Grab the edge we're traversing into the node
                Edge* in_edge = augmented.graph.get_edge(make_pair(to_right_side(to_visit(prev)),
                                                                       to_left_side(to_visit(here))));
                // If there's less support on the in edge than on the node,
                // knock it down. We do this separately in each dimension.
                node_support = support_min(node_support, augmented.edge_supports.at(in_edge));
                
                // Ditto for the edge we're traversing out of the node
                Edge* out_edge = augmented.graph.get_edge(make_pair(to_right_side(to_visit(here)),
                                                                        to_left_side(to_visit(next))));
                node_support = support_min(node_support, augmented.edge_supports.at(out_edge));
                
                
                // Add support in to the total support for the alt. Scale by node length.
                total_support += node->sequence().size() * node_support;
                // Take this as the minimum support if it's the first node, and min it with the min support otherwise.
                min_support = (i == 1 ? node_support : support_min(min_support, node_support));

                // Update minimum likelihood in the alt path
                min_likelihood = min(min_likelihood, augmented.node_likelihoods.at(node));
                
                // TODO: use edge likelihood here too?
                    
            }
            
            if(path.mapping_size() == 2) {
                // We just have the anchoring nodes and the edge between them.
                // Look at that edge specially.
                Edge* edge = augmented.graph.get_edge(make_pair(to_right_side(to_visit(path.mapping(0))),
                                                                    to_left_side(to_visit(path.mapping(1)))));
                
                // Only use the support on the edge
                total_support = augmented.edge_supports.at(edge);
                min_support = total_support;
                
                // And the likelihood on the edge
                min_likelihood = augmented.edge_likelihoods.at(edge);
                
            }
            
            // Calculate the min and average FractionalSupports
            min_supports.push_back(from_support(min_support));
            // The support needs to get divided by bases, unless we're just a
            // single edge empty allele, in which case we're special.
            average_supports.push_back(sequence_stream.str().size() > 0 ? 
                from_support(total_support) / sequence_stream.str().size() : 
                from_support(total_support));
            
            // Fill in the support for this allele in the Locus, by rounding the appropriate FractionalSupport
            *locus.add_support() = to_support(useAverageSupport ? average_supports.back() : min_supports.back());
            
            // Fill in the other vectors
            // TODO: add this info to Locus? Or calculate later when emitting VCF?
            sequences.push_back(sequence_stream.str());
            id_lists.push_back(id_stream.str());
            min_likelihoods.push_back(min_likelihood);
        }
        
        // TODO: complain if multiple copies of the same string exist???
        
        // Decide which support vector we use to actually decide
        vector<FractionalSupport>& supports = useAverageSupport ? average_supports : min_supports;
        
        // Now look at all the paths for the site and pick the top 2.
        int best_allele = -1;
        int second_best_allele = -1;
        for(size_t i = 0; i < locus.allele_size(); i++) {
            if(best_allele == -1 || total(supports[best_allele]) <= total(supports[i])) {
                // We have a new best. Demote the old best.
                second_best_allele = best_allele;
                best_allele = i;
            } else if(second_best_allele == -1 || total(supports[second_best_allele]) <= total(supports[i])) {
                // We're not better than the best, but we can demote the second best.
                second_best_allele = i;
            }
        }
        
        // We should always have a best allele; we may sometimes have a second best.
        assert(best_allele != -1);
        
        if(best_allele == 0 && second_best_allele == -1) {
            // This site isn't variable; don't bother with it.
            continue;
        }
        
        // Since the variable part of the site is after the first anchoring node, where does it start?
        size_t variation_start = index.by_id.at(site->start().node_id()).first
                                 + augmented.graph.get_node(site->start().node_id())->sequence().size();
        
        // Find which coordinate bin the variation start is in, so we can get the typical local support
        int bin = variation_start / refBinSize;
        if (bin == binnedSupport.size()) {
            --bin;
        }
        const FractionalSupport& baseline_support = binnedSupport[bin];

        // Decide if we're an indel. We're an indel if the sequence lengths
        // aren't all equal between the ref and the alleles we're going to call.
        // TODO: is this backwards?
        bool is_indel = sequences.front().size() == sequences[best_allele].size() &&
            (second_best_allele == -1 || sequences.front().size() == sequences[second_best_allele].size());

        // We need to decide what to scale the bias limits by. We scale them up if this is an indel.
        double bias_multiple = is_indel ? indelBiasMultiple : 1.0;
        
        // How much support do we have for the top two alleles?
        FractionalSupport site_support = supports.at(best_allele);
        if(second_best_allele != -1) {
            site_support += supports.at(second_best_allele);
        }
        
        // Pull out the different supports. Some of them may be the same.
        FractionalSupport ref_support = supports.at(0);
        FractionalSupport best_support = supports.at(best_allele);
        FractionalSupport second_best_support = make_pair(0.0, 0.0);
        if(second_best_allele != -1) {
            second_best_support = supports.at(second_best_allele);
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
        if(total(site_support) >= total(baseline_support) * minFractionForCall) {
            // We have enough to emit a call here.
            
            // If best and second best are close enough to be het, we call het.
            // Otherwise, we call hom best.
            
            // We decide closeness differently depending on whether best is ref or not.
            // In practice, we use this to slightly penalize homozygous ref calls
            // (by setting maxRefHetBias higher than maxHetBias) and rather make a less
            // supported alt call instead.  This boost max sensitivity, and because
            // everything is homozygous ref by default in VCF, any downstream filters
            // will effectively reset these calls back to homozygous ref. 
            double bias_limit = (best_allele == 0) ? maxRefHetBias : maxHetBias;
            
#ifdef debug
            cerr << best_allele << ", " << best_support << " and "
                << second_best_allele << ", " << second_best_support << endl;
            
            cerr << bias_limit * bias_multiple * total(second_best_support) << " vs "
                << total(best_support) << endl;
                
            cerr << total(second_best_support) << " vs " << minTotalSupportForCall << endl;
#endif
            
            if(second_best_allele != -1 &&
                bias_limit * bias_multiple * total(second_best_support) >= total(best_support) &&
                total(best_support) >= minTotalSupportForCall &&
                total(second_best_support) >= minTotalSupportForCall) {
                // There's a second best allele, and it's not too biased to
                // call, and both alleles exceed the minimum to call them
                // present.
                
                // Say both are present
                genotype.add_allele(best_allele);
                genotype.add_allele(second_best_allele);
                
                // Compute the likelihood for a best/second best het
                // Quick quality: combine likelihood and depth, using poisson for latter
                // TODO: revize which depth (cur: avg) / likelihood (cur: min) pair to use
                gen_likelihood = ln_to_log10(poisson_prob_ln(total(best_support), 0.5 * total(baseline_support))) +
                    ln_to_log10(poisson_prob_ln(total(second_best_support), 0.5 * total(baseline_support)));
                gen_likelihood += min_likelihoods.at(best_allele) + min_likelihoods.at(second_best_allele);
                
                // Save the likelihood in the Genotype
                genotype.set_likelihood(log10_to_ln(gen_likelihood));
                
                // Get minimum support for filter (not assuming it's second_best just to be sure)
                min_site_support = min(total(second_best_support), total(best_support));
                
                // Make the call
                *locus.add_genotype() = genotype;
                
            } else if(total(best_support) >= minTotalSupportForCall) {
                // The second best allele isn't present or isn't good enough,
                // but the best allele has enough coverage that we can just call
                // two of it.
                
                // Say the best is present twice
                genotype.add_allele(best_allele);
                genotype.add_allele(best_allele);
                
                // Compute the likelihood for hom best allele
                gen_likelihood = ln_to_log10(poisson_prob_ln(total(best_support), total(baseline_support)));
                gen_likelihood += min_likelihoods.at(best_allele);

                // Save the likelihood in the Genotype
                genotype.set_likelihood(log10_to_ln(gen_likelihood));

                // Get minimum support for filter
                min_site_support = total(best_support);
                
                // Make the call
                *locus.add_genotype() = genotype;

            } else {
                // We can't really call this as anything.
                
                // Don't add the genotype to the locus
            }
        } else {
            // Depth too low. Say we have no idea.
            // TODO: elide variant?
            
            // Don't add the genotype to the locus
        }
        
        // Find the total support for the Locus across all alleles
        FractionalSupport locus_support;
        for (auto& s : supports) {
            // Sum up all the FractionalSupports form all alleles (even the non-best/second-best).
            locus_support += s;
        }
        // Round to ints and save
        *locus.mutable_overall_support() = to_support(locus_support);
        
        auto emit_variant = [&](const Locus& locus) {
        
            // Make a Variant
            vcflib::Variant variant;
            variant.sequenceName = contigName;
            variant.setVariantCallFile(vcf);
            variant.quality = 0;
            variant.position = variation_start + 1 + variantOffset;
            
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
                string extra_base = char_to_string(index.sequence.at(variation_start - 1));
                
                for(auto& seq : sequences) {
                    // Stick it on the front of all the allele sequences
                    seq = extra_base + seq;
                }
                
                // Budge the variant left
                variant.position--;
            }
            
            // Add the ref allele to the variant
            create_ref_allele(variant, sequences.front());
            
            // Add the best allele
            assert(best_allele != -1);
            int best_alt = add_alt_allele(variant, sequences.at(best_allele));
            
            int second_best_alt = (second_best_allele == -1) ? -1 : add_alt_allele(variant, sequences.at(second_best_allele));
            
            
            // Say we're going to spit out the genotype for this sample.        
            variant.format.push_back("GT");
            auto& genotype_vector = variant.samples[sampleName]["GT"];

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
            
            // Add depth for the variant and the samples
            FractionalSupport total_support = ref_support;
            if(best_allele != 0) {
                // Add in the depth from the best allele, which is not already counted
                total_support += best_support;
            }
            if(second_best_allele != -1 && second_best_allele != 0) {
                // Add in the depth from the second allele
                total_support += second_best_support;
            }
            string depth_string = to_string((int64_t)round(total(total_support)));
            variant.format.push_back("DP");
            variant.samples[sampleName]["DP"].push_back(depth_string);
            variant.info["DP"].push_back(depth_string); // We only have one sample, so variant depth = sample depth
            
            // Also allelic depths
            variant.format.push_back("AD");
            variant.samples[sampleName]["AD"].push_back(to_string((int64_t)round(total(ref_support))));
            // And strand biases
            variant.format.push_back("SB");
            variant.samples[sampleName]["SB"].push_back(to_string((int64_t)round(ref_support.first)));
            variant.samples[sampleName]["SB"].push_back(to_string((int64_t)round(ref_support.second)));
            // Also allelic likelihoods (from minimum values found on their paths)
            variant.format.push_back("AL");
            variant.samples[sampleName]["AL"].push_back(to_string_ss(min_likelihoods.at(0)));
            if(best_allele != 0) {
                // If our best allele isn't ref, it comes next.
                variant.samples[sampleName]["AD"].push_back(to_string((int64_t)round(total(best_support))));
                variant.samples[sampleName]["SB"].push_back(to_string((int64_t)round(best_support.first)));
                variant.samples[sampleName]["SB"].push_back(to_string((int64_t)round(best_support.second)));
                variant.samples[sampleName]["AL"].push_back(to_string_ss(min_likelihoods.at(best_allele)));
            }
            if(second_best_allele != -1 && second_best_allele != 0) {
                // If our second best allele is real and not ref, it comes next.
                variant.samples[sampleName]["AD"].push_back(to_string((int64_t)round(total(second_best_support))));
                variant.samples[sampleName]["SB"].push_back(to_string((int64_t)round(second_best_support.first)));
                variant.samples[sampleName]["SB"].push_back(to_string((int64_t)round(second_best_support.second)));
                variant.samples[sampleName]["AL"].push_back(to_string_ss(min_likelihoods.at(second_best_allele)));
            }
            
            // And total alt allele depth for the alt alleles
            FractionalSupport alt_support = make_pair(0.0, 0.0);
            if(best_allele != 0) {
                alt_support += best_support;
            }
            if(second_best_allele != -1 && second_best_allele != 0) {
                alt_support += second_best_support;
            }
            variant.format.push_back("XAAD");
            variant.samples[sampleName]["XAAD"].push_back(to_string((int64_t)round(total(alt_support))));

            // Copy over the quality value
            variant.quality = -10. * log10(1. - pow(10, gen_likelihood));

            // Apply Min Allele Depth cutoff and store result in Filter column
            variant.filter = min_site_support >= min_mad_for_filter ? "PASS" : "FAIL";
            
            if(can_write_alleles(variant)) {
                // No need to check for collisions because we assume sites are correctly found.
                
                // Output the created VCF variant.
                cout << variant << endl;
            } else {
                if (verbose) {
                    cerr << "Variant is too large" << endl;
                }
                // TODO: account for the 1 base we added extra if it was a pure
                // insert.
                basesLost += sequences.at(best_allele).size();
            }
            
        };
        
        if (convert_to_vcf) {
            // Emit the variant for this Locus
            emit_variant(locus);
        } else {
            // Emit the locus itself
            locus_buffer.push_back(locus);
            stream::write_buffered(cout, locus_buffer, locus_buffer_size);
        }
    }
    
    if (!convert_to_vcf) {
        // Flush the buffer of Locus objects we have to write
        stream::write_buffered(cout, locus_buffer, locus_buffer_size);
    }
    
    // Announce how much we can't show.
    if (verbose) {
        cerr << "Had to drop " << basesLost << " bp of unrepresentable variation." << endl;
    }
}

bool Call2Vcf::is_reference(const SnarlTraversal& trav, AugmentedGraph& augmented, const PathIndex& primary_path) {
    
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
    
    // Make sure we visit a ref start node
    if (!experience_visit(trav.snarl().start())) {
        return false;
    }
    
    // Then all the internal nodes
    for (size_t i = 0; i < trav.visits_size(); i++) {
        if (!experience_visit(trav.visits(i))) {
            return false;
        }
    }
    
    // And finally the end node
    if (!experience_visit(trav.snarl().end())) {
        return false;
    }
    
    // And if we make it through it's a reference traversal.
    return true;
        
}

}
