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

//#define debug

namespace glenn2vcf {

using namespace vg;

// TODO:
//  - Decide if we need to have sibling alts detect (somehow) and coordinate with each other
//  - Parallelize variant generation
//  - Make variant stamping out some kind of function, don't duplicate the same variant construction code 6 times

// How many bases may we put in an allele in VCF if we expect GATK to be able to
// parse it?
// 0 means no maximum is enforced.
const static int MAX_ALLELE_LENGTH = 0;

// Minimum log likelihood
const static double LOG_ZERO = (double)-1e100;

// convert to string using stringstream (to replace to_string when we want sci. notation)
template <typename T>
std::string to_string_ss(T val) {
    stringstream ss;
    ss << val;
    return ss.str();
}

// Poisson utility functions from Freebayes (when merging into vg
// these should get offloaded into utility.cpp or something
static long double gammaln(
    long double x
    ) {

    long double cofactors[] = { 76.18009173, 
                                -86.50532033,
                                24.01409822,
                                -1.231739516,
                                0.120858003E-2,
                                -0.536382E-5 };    

    long double x1 = x - 1.0;
    long double tmp = x1 + 5.5;
    tmp -= (x1 + 0.5) * log(tmp);
    long double ser = 1.0;
    for (int j=0; j<=5; j++) {
        x1 += 1.0;
        ser += cofactors[j]/x1;
    }
    long double y =  (-1.0 * tmp + log(2.50662827465 * ser));

    return y;
}

static long double factorial(
    int n
    ) {
    if (n < 0) {
        return (long double)0.0;
    }
    else if (n == 0) {
        return (long double)1.0;
    }
    else {
        return exp(gammaln(n + 1.0));
    }
}

long double poissonp(int observed, int expected) {
    return (double) pow((double) expected, (double) observed) * (double) pow(M_E, (double) -expected) / factorial(observed);
}


/**
 * Holds indexes of the reference: position to node, node to position and
 * orientation, and the full reference string.
 */
struct ReferenceIndex {
    // Index from node ID to first position on the reference string and
    // orientation it occurs there.
    std::map<int64_t, std::pair<size_t, bool>> byId;
    
    // Index from start position on the reference to the oriented node that
    // begins there.  Some nodes may be backward (orientation true) at their
    // canonical reference positions. In this case, the last base of the node
    // occurs at the given position.
    std::map<size_t, vg::NodeTraversal> byStart;
    
    // The actual sequence of the reference.
    std::string sequence;
};

/**
 * We need to suppress overlapping variants, but interval trees are hard to
 * write. This accomplishes the collision check with a massive bit vector.
 */
struct IntervalBitfield {
    // Mark every position that's used in a variant
    std::vector<bool> used;
    
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
typedef std::pair<double, double> Support;

/**
 * Add two Support values together, accounting for strand.
 */
Support operator+(const Support& one, const Support& other) {
    return std::make_pair(one.first + other.first, one.second + other.second);
}

/**
 * Add in a Support to another.
 */
Support& operator+=(Support& one, const Support& other) {
    one.first += other.first;
    one.second += other.second;
    return one;
}


/**
 * Scale a support by a factor.
 */
template<typename Scalar>
Support operator*(const Support& support, const Scalar& scale) {
    return std::make_pair(support.first * scale, support.second * scale);
}

/**
 * Scale a support by a factor, the other way
 */
template<typename Scalar>
Support operator*(const Scalar& scale, const Support& support) {
    return std::make_pair(support.first * scale, support.second * scale);
}

/**
 * Divide a support by a factor.
 */
template<typename Scalar>
Support operator/(const Support& support, const Scalar& scale) {
    return std::make_pair(support.first / scale, support.second / scale);
}
    
/**
 * Allow printing a support.
 */
std::ostream& operator<<(std::ostream& stream, const Support& support) {
    return stream << support.first << "," << support.second;
}

/**
 * Get the total read support in a support.
 */
double total(const Support& support) {
    return support.first + support.second;
}

/**
 * Get the strand bias of a support.
 */
double strand_bias(const Support& support) {
    return std::max(support.first, support.second) / (support.first + support.second);
}

/**
 * Get the minimum support of a pair of supports, by taking the min in each
 * orientation.
 */
Support support_min(const Support& a, const Support& b) {
    return std::make_pair(std::min(a.first, b.first), std::min(a.second, b.second));
}

/**
 * Make a letter into a full string because apparently that's too fancy for the
 * standard library.
 */
std::string char_to_string(const char& letter) {
    std::string toReturn;
    toReturn.push_back(letter);
    return toReturn;
}

/**
 * Write a minimal VCF header for a single-sample file.
 */
void write_vcf_header(std::ostream& stream, std::string& sample_name, std::string& contig_name, size_t contig_size,
                      int min_mad_for_filter) {
    stream << "##fileformat=VCFv4.2" << std::endl;
    stream << "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">" << std::endl;
    stream << "##INFO=<ID=XREF,Number=0,Type=Flag,Description=\"Present in original graph\">" << std::endl;
    stream << "##INFO=<ID=XSEE,Number=.,Type=String,Description=\"Original graph node:offset cross-references\">" << std::endl;
    stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">" << std::endl;
    stream << "##FILTER=<ID=FAIL,Description=\"Variant does meat minimum allele read support threshold of " << min_mad_for_filter << "\">" <<endl;
    stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << std::endl;
    stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
    stream << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << std::endl;
    stream << "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Forward and reverse support for ref and alt alleles.\">" << std::endl;
    // We need this field to stratify on for VCF comparison. The info is in SB but vcfeval can't pull it out
    stream << "##FORMAT=<ID=XAAD,Number=1,Type=Integer,Description=\"Alt allele read count.\">" << std::endl;
    stream << "##FORMAT=<ID=AL,Number=.,Type=Float,Description=\"Allelic likelihoods for the ref and alt alleles in the order listed\">" << std::endl;
    
    if(!contig_name.empty()) {
        // Announce the contig as well.
        stream << "##contig=<ID=" << contig_name << ",length=" << contig_size << ">" << std::endl;
    }
    stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << std::endl;
}

/**
 * Create the reference allele for an empty vcflib Variant, since apaprently
 * there's no method for that already. Must be called before any alt alleles are
 * added.
 */
void create_ref_allele(vcflib::Variant& variant, const std::string& allele) {
    // Set the ref allele
    variant.ref = allele;
    
    for(size_t i = 0; i < variant.ref.size(); i++) {
        // Look at all the bases
        if(variant.ref[i] != 'A' && variant.ref[i] != 'C' && variant.ref[i] != 'G' && variant.ref[i] != 'T') {
            // Correct anything bogus (like "X") to N
            variant.ref[i] = 'N';
        }
    }
    
    // Make it 0 in the alleles-by-index list
    variant.alleles.push_back(allele);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();
}

/**
 * Add a new alt allele to a vcflib Variant, since apaprently there's no method
 * for that already.
 *
 * If that allele already exists in the variant, does not add it again.
 *
 * Returns the allele number (0, 1, 2, etc.) corresponding to the given allele
 * string in the given variant. 
 */
int add_alt_allele(vcflib::Variant& variant, const std::string& allele) {
    // Copy the allele so we can throw out bad characters
    std::string fixed(allele);
    
    for(size_t i = 0; i < fixed.size(); i++) {
        // Look at all the bases
        if(fixed[i] != 'A' && fixed[i] != 'C' && fixed[i] != 'G' && fixed[i] != 'T') {
            // Correct anything bogus (like "X") to N
            fixed[i] = 'N';
        }
    }
    
    for(int i = 0; i < variant.alleles.size(); i++) {
        if(variant.alleles[i] == fixed) {
            // Already exists
            return i;
        }
    }

    // Add it as an alt
    variant.alt.push_back(fixed);
    // Make it next in the alleles-by-index list
    variant.alleles.push_back(fixed);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();

    // We added it in at the end
    return variant.alleles.size() - 1;
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
bool mapping_is_perfect_match(const vg::Mapping& mapping) {
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
 * Get the length of a path through nodes, in base pairs.
 */
size_t bp_length(const std::list<vg::NodeTraversal>& path) {
    size_t length = 0;
    for(auto& traversal : path) {
        // Sum up length of each node's sequence
        length += traversal.node->sequence().size();
    }
    return length;
}

/**
 * Get the minimum support of all nodes and edges in path
 */
Support min_support_in_path(vg::VG& graph,
                            const std::map<vg::Node*, Support>& nodeReadSupport,
                            const std::map<vg::Edge*, Support>& edgeReadSupport,
                            const std::list<vg::NodeTraversal>& path) {
    if (path.empty()) {
        return Support();
    }
    auto cur = path.begin();
    auto next = path.begin();
    ++next;
    Support minSupport = nodeReadSupport.count(cur->node) ? nodeReadSupport.at(cur->node) : Support();
    for (; next != path.end(); ++cur, ++next) {
        // check the node support
        Support support = nodeReadSupport.count(next->node) ? nodeReadSupport.at(next->node) : Support();
        minSupport = support_min(minSupport, support);
        
        // check the edge support
        Edge* edge = graph.get_edge(*cur, *next);
        assert(edge != NULL);
        Support edgeSupport = edgeReadSupport.count(edge) ? edgeReadSupport.at(edge) : Support();
        minSupport = support_min(minSupport, edgeSupport);
    }
    return minSupport;
}

/**
 * Do a breadth-first search left from the given node traversal, and return
 * lengths and paths starting at the given node and ending on the indexed
 * reference path. Refuses to visit nodes with no support.
 */
std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> bfs_left(vg::VG& graph,
    vg::NodeTraversal node, const ReferenceIndex& index,
    const std::map<vg::Node*, Support>& nodeReadSupport,
    const std::map<vg::Edge*, Support>& edgeReadSupport,
    int64_t maxDepth = 10, bool stopIfVisited = false) {

    // Holds partial paths we want to return, with their lengths in bp.
    std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> toReturn;
    
    // Do a BFS
    
    // This holds the paths to get to NodeTraversals to visit (all of which will
    // end with the node we're starting with).
    std::list<std::list<vg::NodeTraversal>> toExtend;
    
    // This keeps a set of all the oriented nodes we already got to and don't
    // need to queue again.
    std::set<vg::NodeTraversal> alreadyQueued;
    
    // Start at this node at depth 0
    toExtend.emplace_back(std::list<vg::NodeTraversal> {node});
    // Mark this traversal as already queued
    alreadyQueued.insert(node);
    
    // How many ticks have we spent searching?
    size_t searchTicks = 0;
    // Track how many options we have because size may be O(n).
    size_t stillToExtend = toExtend.size();
    
    while(!toExtend.empty()) {
        // Keep going until we've visited every node up to our max search depth.
        
#ifdef debug
        searchTicks++;
        if(searchTicks % 100 == 0) {
            // Report on how much searching we are doing.
            std::cerr << "Search tick " << searchTicks << ", " << stillToExtend << " options." << std::endl;
        }
#endif
        
        // Dequeue a path to extend.
        // Make sure to move out of the list to avoid a useless copy.
        std::list<vg::NodeTraversal> path(std::move(toExtend.front()));
        toExtend.pop_front();
        stillToExtend--;
        
        // We can't just throw out longer paths, because shorter paths may need
        // to visit a node twice (in opposite orientations) and thus might get
        // rejected later. Or they might overlap with paths on the other side.
        
        // Look up and see if the front node on the path is on our reference
        // path
        if(index.byId.count(path.front().node->id())) {
            // This node is on the reference path. TODO: we don't care if it
            // lands in a place that is itself deleted.
            
            // Say we got to the right place
            toReturn.emplace(bp_length(path), std::move(path));
            
            // Don't bother looking for extensions, we already got there.
        } else if(path.size() <= maxDepth) {
            // We haven't hit the reference path yet, but we also haven't hit
            // the max depth. Extend with all the possible extensions.
            
            // Look left
            vector<vg::NodeTraversal> prevNodes;
            graph.nodes_prev(path.front(), prevNodes);
            
            for(auto prevNode : prevNodes) {
                // For each node we can get to
                vg::Edge* edge = graph.get_edge(prevNode, path.front());
                assert(edge != NULL);
                
                if((!nodeReadSupport.empty() && (!nodeReadSupport.count(prevNode.node) ||
                    total(nodeReadSupport.at(prevNode.node)) == 0)) ||
                   (!edgeReadSupport.empty() && (!edgeReadSupport.count(edge) ||
                    total(edgeReadSupport.at(edge)) == 0))) {
                    
                    // We have no support at all for visiting this node (but we
                    // do have some node read support data)
                    continue;
                }
                
                if(stopIfVisited && alreadyQueued.count(prevNode)) {
                    // We already have a way to get here.
                    continue;
                }
            
                // Make a new path extended left with the node
                std::list<vg::NodeTraversal> extended(path);
                extended.push_front(prevNode);
                toExtend.emplace_back(std::move(extended));
                stillToExtend++;
                
                // Remember we found a way to this node, so we don't try and
                // visit it other ways.
                alreadyQueued.insert(prevNode);
            }
        }
        
    }
    
    return toReturn;
}

/**
 * Flip a NodeTraversal around and return the flipped copy.
 */
vg::NodeTraversal flip(vg::NodeTraversal toFlip) {
    return vg::NodeTraversal(toFlip.node, !toFlip.backward);
}

/**
 * Do a breadth-first search right from the given node traversal, and return
 * lengths and paths starting at the given node and ending on the indexed
 * reference path.
 */
std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> bfs_right(vg::VG& graph,
    vg::NodeTraversal node, const ReferenceIndex& index,
    const std::map<vg::Node*, Support>& nodeReadSupport,
    const std::map<vg::Edge*, Support>& edgeReadSupport,
    int64_t maxDepth = 10, bool stopIfVisited = false) {

    // Look left from the backward version of the node.
    auto toConvert = bfs_left(graph, flip(node), index, nodeReadSupport,
                              edgeReadSupport, maxDepth, stopIfVisited);
    
    // Since we can't modify set records in place, we need to do a copy
    std::set<std::pair<size_t, std::list<vg::NodeTraversal>>> toReturn;
    
    for(auto lengthAndPath : toConvert) {
        // Flip every path to run the other way
        lengthAndPath.second.reverse();
        for(auto& traversal : lengthAndPath.second) {
            // And invert the orientation of every node in the path in place.
            traversal = flip(traversal);
        }
        // Stick it in the new set
        toReturn.emplace(std::move(lengthAndPath));
    }
    
    return toReturn;
}


/**
 * Given a vg graph, an edge or node in the graph, and an index for the
 * reference path, look out from the edge or node in both directions to find a
 * shortest bubble relative to the path, with a consistent orientation. The
 * bubble may not visit the same node twice.
 *
 * Exactly one of edge and node must be null, and one not null.
 *
 * Takes a max depth for the searches producing the paths on each side.
 * 
 * Return the ordered and oriented nodes in the bubble, with the outer nodes
 * being oriented forward along the named path, and with the first node coming
 * before the last node in the reference.  Also return the minimum support
 * found on any edge or node in the bubble (including the reference node endpoints
 * and their edges which aren't stored in the path)
 */
std::pair<Support, std::vector<vg::NodeTraversal> >
find_bubble(vg::VG& graph, vg::Node* node, vg::Edge* edge, const ReferenceIndex& index,
            const std::map<vg::Node*, Support>& nodeReadSupport,
            const std::map<vg::Edge*, Support>& edgeReadSupport, int64_t maxDepth,
            size_t max_bubble_paths) {

    // What are we going to find our left and right path halves based on?
    NodeTraversal left_traversal;
    NodeTraversal right_traversal;

    if(edge != nullptr) {
        // Be edge-based
        
        // Find the nodes at the ends of the edges. Look at them traversed in the
        // edge's local orientation.
        left_traversal = NodeTraversal(graph.get_node(edge->from()), edge->from_start());
        right_traversal = NodeTraversal(graph.get_node(edge->to()), edge->to_end());
        
    } else {
        // Be node-based
        assert(node != nullptr);
        left_traversal = right_traversal = NodeTraversal(node);
    }
    
#ifdef debug
    std::cerr << "Starting from: " << left_traversal << ", " << right_traversal << std::endl;
#endif

    // Find paths on both sides, with nodes on the primary path at the outsides
    // and this edge in the middle. Returns path lengths and paths in pairs in a
    // set.
    auto leftPaths = bfs_left(graph, left_traversal, index, nodeReadSupport,
                              edgeReadSupport, maxDepth);
    auto rightPaths = bfs_right(graph, right_traversal, index, nodeReadSupport,
                                edgeReadSupport, maxDepth);
    
    // Find a combination of two paths which gets us to the reference in a
    // consistent orientation (meaning that when you look at the ending nodes'
    // Mappings in the reference path, the ones with minimal ranks have the same
    // orientations) and which doesn't use the same nodes on both sides.
    // Track support of up to max_bubble_paths combinations, and return the
    // highest
    std::pair<Support, std::vector<NodeTraversal> > bestBubblePath;
    int bubbleCount = 0;
    
    // We need to look in different combinations of lists.
    auto testCombinations = [&](const std::list<std::list<vg::NodeTraversal>>& leftList,
        const std::list<std::list<vg::NodeTraversal>>& rightList) {

        for(auto leftPath : leftList) {
            // Figure out the relative orientation for the leftmost node.
#ifdef debug        
            std::cerr << "Left path: " << std::endl;
            for(auto traversal : leftPath ) {
                std::cerr << "\t" << traversal << std::endl;
            }
#endif    
            // Split out its node pointer and orientation
            auto leftNode = leftPath.front().node;
            auto leftOrientation = leftPath.front().backward;
            
            // Get where it falls in the reference as a position, orientation pair.
            auto leftRefPos = index.byId.at(leftNode->id());
            
            // We have a backward orientation relative to the reference path if we
            // were traversing the anchoring node backwards, xor if it is backwards
            // in the reference path.
            bool leftRelativeOrientation = leftOrientation != leftRefPos.second;
            
            // Make a set of all the nodes in the left path
            std::set<int64_t> leftPathNodes;
            for(auto visit : leftPath) {
                leftPathNodes.insert(visit.node->id());
            }

            // Get the minimum support in the left path
            Support minLeftSupport = min_support_in_path(
                graph, nodeReadSupport, edgeReadSupport, leftPath);
            
            for(auto rightPath : rightList) {
                // Figure out the relative orientation for the rightmost node.
#ifdef debug            
                std::cerr << "Right path: " << std::endl;
                for(auto traversal : rightPath ) {
                    std::cerr << "\t" << traversal << std::endl;
                }
#endif            
                // Split out its node pointer and orientation
                // Remember it's at the end of this path.
                auto rightNode = rightPath.back().node;
                auto rightOrientation = rightPath.back().backward;
                
                // Get where it falls in the reference as a position, orientation pair.
                auto rightRefPos = index.byId.at(rightNode->id());
                
                // We have a backward orientation relative to the reference path if we
                // were traversing the anchoring node backwards, xor if it is backwards
                // in the reference path.
                bool rightRelativeOrientation = rightOrientation != rightRefPos.second;

                // Get the minimum support in the right path
                Support minRightSupport = min_support_in_path(
                    graph, nodeReadSupport, edgeReadSupport, rightPath);
                
                if(leftRelativeOrientation == rightRelativeOrientation &&
                    ((!leftRelativeOrientation && leftRefPos.first < rightRefPos.first) ||
                    (leftRelativeOrientation && leftRefPos.first > rightRefPos.first))) {
                    // We found a pair of paths that get us to and from the
                    // reference without turning around, and that don't go back to
                    // the reference before they leave.

                    // Get the minimum support of combined left and right paths
                    Support minFullSupport = support_min(minLeftSupport, minRightSupport);
                    
                    // Start with the left path
                    std::vector<vg::NodeTraversal> fullPath{leftPath.begin(), leftPath.end()};
                    
                    // We need to detect overlap with the left path
                    bool overlap = false;
                    
                    // If we're starting from an edge, we keep the first node on
                    // the right path. If we're starting from a node, we need to
                    // discard it because it's just another copy of our node
                    // we're starting with.
                    for(auto it = (edge != nullptr ? rightPath.begin() : ++(rightPath.begin())); it != rightPath.end(); ++it) {
                        // For all but the first node on the right path, add that in
                        fullPath.push_back(*it);
                        
                        if(leftPathNodes.count((*it).node->id())) {
                            // We already visited this node on the left side. Try
                            // the next right path instead.
                            overlap = true;
                        }
                    }
                    
                    if(overlap) {
                        // Can't combine this right with this left, as they share
                        // nodes and we can't handle the copy number implications.
                        // Try the next right.
                        // TODO: handle the copy number implications.
                        continue;
                    }
                    
                    if(leftRelativeOrientation) {
                        // Turns out our anchored path is backwards.
                        
                        // Reorder everything the other way
                        std::reverse(fullPath.begin(), fullPath.end());
                        
                        for(auto& traversal : fullPath) {
                            // Flip each traversal
                            traversal = flip(traversal);
                        }
                    }

                    
#ifdef debug        
                    std::cerr << "Merged path:" << std::endl;
                    for(auto traversal : fullPath) {
                        std::cerr << "\t" << traversal << std::endl;
                    }                    
#endif
                    // update our best path by seeing if we've found one with higher min support
                    if (total(minFullSupport) > total(bestBubblePath.first)) {
                        bestBubblePath.first = minFullSupport;
                        bestBubblePath.second = fullPath;
                    }

                    // keep things from getting out of hand
                    if (++bubbleCount >= max_bubble_paths) {
                        return bestBubblePath;
                    }
                }
            }
        }
        
        // Return the best path along with its min support
        // (could be empty)
        return bestBubblePath;
        
    };
    
    // Convert sets to lists, which requires a copy again...
    std::list<std::list<vg::NodeTraversal>> leftConverted;
    for(auto lengthAndPath : leftPaths) {
        leftConverted.emplace_back(std::move(lengthAndPath.second));
    }
    std::list<std::list<vg::NodeTraversal>> rightConverted;
    for(auto lengthAndPath : rightPaths) {
        rightConverted.emplace_back(std::move(lengthAndPath.second));
    }
    
    // Look for a valid combination, or return an empty path if one iesn't
    // found.
    return testCombinations(leftConverted, rightConverted);
    
}

/**
 * Trace out the reference path in the given graph named by the given name.
 * Returns a structure with useful indexes of the reference.
 */
ReferenceIndex trace_reference_path(vg::VG& vg, std::string refPathName, bool verbose) {
    // Make sure the reference path is present
    assert(vg.paths.has_path(refPathName));
    
    // We'll fill this in and then return it.
    ReferenceIndex index;
    
    // We're also going to build the reference sequence string
    std::stringstream refSeqStream;
    
    // What base are we at in the reference
    size_t referenceBase = 0;
    
    // What was the last rank? Ranks must always go up.
    int64_t lastRank = -1;
    
    for(auto mapping : vg.paths.get_path(refPathName)) {
        // All the mappings need to be perfect matches.
        assert(mapping_is_perfect_match(mapping));
    
        if(!index.byId.count(mapping.position().node_id())) {
            // This is the first time we have visited this node in the reference
            // path.
            
            // Add in a mapping.
            index.byId[mapping.position().node_id()] = 
                std::make_pair(referenceBase, mapping.position().is_reverse());
#ifdef debug
            std::cerr << "Node " << mapping.position().node_id() << " rank " << mapping.rank()
                << " starts at base " << referenceBase << " with "
                << vg.get_node(mapping.position().node_id())->sequence() << std::endl;
#endif
            
            // Make sure ranks are monotonically increasing along the path.
            assert(mapping.rank() > lastRank);
            lastRank = mapping.rank();
        } else if (verbose) {
            std::cerr << "Warning: Cycle detected at " << mapping.position().node_id()
                      << " in reference path " << refPathName << endl;
        }
        
        // Find the node's sequence
        std::string sequence = vg.get_node(mapping.position().node_id())->sequence();
        
        while(referenceBase == 0 && sequence.size() > 0 &&
            (sequence[0] != 'A' && sequence[0] != 'T' && sequence[0] != 'C' &&
            sequence[0] != 'G' && sequence[0] != 'N')) {
            
            // If the path leads with invalid characters (like "X"), throw them
            // out when computing reference path positions.
            
            // TODO: this is a hack to deal with the debruijn-brca1-k63 graph,
            // which leads with an X.

            if (verbose) {
                std::cerr << "Warning: dropping invalid leading character "
                          << sequence[0] << " from node " << mapping.position().node_id()
                          << std::endl;
            }
                
            sequence.erase(sequence.begin());
        }
        
        if(mapping.position().is_reverse()) {
            // Put the reverse sequence in the reference path
            refSeqStream << vg::reverse_complement(sequence);
        } else {
            // Put the forward sequence in the reference path
            refSeqStream << sequence;
        }
            
        // Say that this node appears here along the reference in this
        // orientation.
        index.byStart[referenceBase] = vg::NodeTraversal(
            vg.get_node(mapping.position().node_id()), mapping.position().is_reverse()); 
            
        // Whether we found the right place for this node in the reference or
        // not, we still need to advance along the reference path. We assume the
        // whole node (except any leading bogus characters) is included in the
        // path (since it sort of has to be, syntactically, unless it's the
        // first or last node).
        referenceBase += sequence.size();
        
        // TODO: handle leading bogus characters in calls on the first node.
    }
    
    // Create the actual reference sequence we will use
    index.sequence = refSeqStream.str();
    
    // Announce progress.
    if (verbose) {
        std::cerr << "Traced " << referenceBase << " bp reference path " << refPathName << "." << std::endl;
    
        if(index.sequence.size() < 100) {
            std::cerr << "Reference sequence: " << index.sequence << std::endl;
        }
    }
    
    // Give back the indexes we have been making
    return index;
}

/**
 * Given a collection of pileups by original node ID, and a set of original node
 * id:offset cross-references in both ref and alt categories, produce a VCF
 * comment line giving the pileup for each of those positions on those nodes.
 * Includes a trailing newline if nonempty.
 *
 * TODO: VCF comments aren't really a thing.
 */
std::string get_pileup_line(const std::map<int64_t, vg::NodePileup>& nodePileups,
    const std::set<std::pair<int64_t, size_t>>& refCrossreferences,
    const std::set<std::pair<int64_t, size_t>>& altCrossreferences) {
    // We'll make a stringstream to write to.
    std::stringstream out;
    
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
        out << std::endl;
        return out.str();
    } else {
        // Give an empty string.
        return "";
    }
}

    std::map<vg::Node*, Support> nodeReadSupport;
    // And read support for the edges
    std::map<vg::Edge*, Support> edgeReadSupport;
    // This holds all the edges that are deletions, by the pointer to the stored
    // Edge object in the VG graph
    std::set<vg::Edge*> deletionEdges;
    
    // This holds where nodes came from (node and offset) in the original, un-
    // augmented graph. For pieces of original nodes, this is where the piece
    // started. For novel nodes, this is where the piece that thois is an
    // alternative to started.
    std::map<vg::Node*, std::pair<int64_t, size_t>> nodeSources;
    
    // We also need to track what edges and nodes are reference (i.e. already
    // known)
    std::set<vg::Node*> knownNodes;
    std::set<vg::Edge*> knownEdges;

/**
 * Parse tsv into an internal format, where we track status and copy number
 * for nodes and edges.
 */
void parse_tsv(const std::string& tsvFile,
               vg::VG& vg,
               std::map<vg::Node*, Support>& nodeReadSupport,
               std::map<vg::Edge*, Support>& edgeReadSupport,
               std::map<vg::Node*, double>& nodeLikelihood,
               std::map<vg::Edge*, double>& edgeLikelihood,
               std::set<vg::Edge*>& deletionEdges,
               std::map<vg::Node*, std::pair<int64_t, size_t>>& nodeSources,
               std::set<vg::Node*>& knownNodes,
               std::set<vg::Edge*>& knownEdges,
               bool verbose) {
    
    // Open up the TSV-file
    std::stringstream tsvStream(tsvFile);

    // Loop through all the lines
    std::string line;
    size_t lineNumber = 0;
    while(std::getline(tsvStream, line)) {
        // For each line
        
        lineNumber++;
        
        if(line == "") {
            // Skip blank lines
            continue;
        }
        
        // Make a stringstream to read out tokens
        std::stringstream tokens(line);
        
        // Read the kind of line this is ("N"ode or "E"dge)
        std::string lineType;
        tokens >> lineType; 
        
        if(lineType == "N") {
            // This is a line about a node
            
            // Read the node ID
            int64_t nodeId;
            tokens >> nodeId;
            
            if(!vg.has_node(nodeId)) {
                throw std::runtime_error("Line " + std::to_string(lineNumber) + ": Invalid node: " + std::to_string(nodeId));
            }
            
            // Retrieve the node we're talking about 
            vg::Node* nodePointer = vg.get_node(nodeId);
            
            // What kind of call is it? Could be "U"ncalled, or "R"eference
            // (i.e. known in the original graph), which we have special
            // handling for.
            std::string callType;
            tokens >> callType;
                        
            // Read the read support and likelihood
            Support readSupport;
            int other_support = 0;
            double likelihood = 0.;
            tokens >> readSupport.first;
            tokens >> readSupport.second;
            tokens >> other_support;
            tokens >> likelihood;
            
#ifdef debug
            std::cerr << "Line " << std::to_string(lineNumber) << ": Node " << nodeId
                << " has read support " << readSupport.first << "," << readSupport.second << endl;
#endif
            
            // Save it
            nodeReadSupport[nodePointer] = readSupport;
            nodeLikelihood[nodePointer] = likelihood;
            
            if(callType == "R") {
                // Note that this is a reference node
                knownNodes.insert(nodePointer);
            }
            
            // Load the original node ID and offset for this node, if present.
            int64_t originalId;
            size_t originalOffset;
            
            if(tokens >> originalId && tokens >> originalOffset && originalId != 0) {
                nodeSources[nodePointer] = std::make_pair(originalId, originalOffset);
            }
            
        } else if(lineType == "E") {
        
            // Read the edge data
            std::string edgeDescription;
            tokens >> edgeDescription;
            
            // Split on commas. We'd just iterate the regex iterator ourselves,
            // but it seems to only split on the first comma if we do that.
            std::vector<string> parts;
            std::regex comma_re(",");
            std::copy(std::sregex_token_iterator(edgeDescription.begin(), edgeDescription.end(), comma_re, -1), std::sregex_token_iterator(), std::back_inserter(parts));
            
            // We need the four fields to describe an edge.
            assert(parts.size() == 4);
            
            // Parse the from node
            int64_t from = std::stoll(parts[0]);
            // And the from_start flag
            bool fromStart = std::stoi(parts[1]);
            // Make a NodeSide for the from side
            vg::NodeSide fromSide(from, !fromStart);
            // Parse the to node
            int64_t to = std::stoll(parts[2]);
            // And the to_end flag
            bool toEnd = std::stoi(parts[3]);
            // Make a NodeSide for the to side
            vg::NodeSide toSide(to, toEnd);
            
            if(!vg.has_edge(std::make_pair(fromSide, toSide))) {
                // Ensure we really have that edge
                throw std::runtime_error("Line " + std::to_string(lineNumber) + ": Edge " + edgeDescription + " not in graph.");
            }
            
            // Get the edge
            vg::Edge* edgePointer = vg.get_edge(std::make_pair(fromSide, toSide));
            
            // Parse the mode
            std::string mode;
            tokens >> mode;
            
            if(mode == "L" || mode == "R") {
                // This is a deletion edge, or an edge in the primary path that
                // may describe a nonzero-length deletion.
#ifdef debug
                std::cerr << "Line " << std::to_string(lineNumber) << ": Edge "
                    << edgeDescription << " may describe a deletion." << endl;
#endif

                // Say it's a deletion
                deletionEdges.insert(edgePointer);
                
                if(mode == "R") {
                    // The reference edges also get marked as such
                    knownEdges.insert(edgePointer);
#ifdef debug
                    std::cerr << "Line " << std::to_string(lineNumber) << ": Edge "
                        << edgeDescription << " is reference." << endl;
#endif
                }

            }
            
            // Read the read support
            Support readSupport;
            int other_support;
            double likelihood;
            tokens >> readSupport.first;
            tokens >> readSupport.second;
            tokens >> other_support;
            tokens >> likelihood;
            
#ifdef debug
            std::cerr << "Line " << std::to_string(lineNumber) << ": Edge " << edgeDescription
                << " has read support " << readSupport.first << "," << readSupport.second << endl;
#endif
                
            // Save it
            edgeReadSupport[edgePointer] = readSupport;
            edgeLikelihood[edgePointer] = likelihood;
        
        } else {
            // This is not a real kind of line
            throw std::runtime_error("Line " + std::to_string(lineNumber) + ": Unknown line type: " + lineType);
        }
        
    }
    if (verbose) {
        std::cerr << "Loaded " << lineNumber << " lines from tsv buffer" << endl;
    }
}

// this was main() in glenn2vcf
// all that's changed for now is that arguments passed in rather than
// parsed from argc/argv (or passed as files)
int call2vcf(

    // Augmented graph
    vg::VG& vg,
    // "glennfile" as string (relic from old pipeline)
    const std::string& glennFile,
    // Option variables
    // What's the name of the reference path in the graph?
    std::string refPathName,
    // What name should we give the contig in the VCF file?
    std::string contigName,
    // What name should we use for the sample in the VCF file?
    std::string sampleName,
    // How far should we offset positions of variants?
    int64_t variantOffset,
    // How many nodes should we be willing to look at on our path back to the
    // primary path? Keep in mind we need to look at all valid paths (and all
    // combinations thereof) until we find a valid pair.
    int64_t maxDepth,
    // What should the total sequence length reported in the VCF header be?
    int64_t lengthOverride,
    // Should we load a pileup and print out pileup info as comments after
    // variants?
    std::string pileupFilename,
    // What fraction of average coverage should be the minimum to call a variant (or a single copy)?
    // Default to 0 because vg call is still applying depth thresholding
    double minFractionForCall,
    // What fraction of the reads supporting an alt are we willing to discount?
    // At 2, if twice the reads support one allele as the other, we'll call
    // homozygous instead of heterozygous. At infinity, every call will be
    // heterozygous if even one read supports each allele.
    double maxHetBias,
    // Like above, but applied to ref / alt ratio (instead of alt / ref)
    double maxRefHetBias,
    // How much should we multiply the bias limits for indels?
    double indelBiasMultiple,
    // What's the minimum integer number of reads that must support a call? We
    // don't necessarily want to call a SNP as het because we have a single
    // supporting read, even if there are only 10 reads on the site.
    size_t minTotalSupportForCall,
    // Bin size used for counting coverage along the reference path.  The
    // bin coverage is used for computing the probability of an allele
    // of a certain depth
    size_t refBinSize,
    // On some graphs, we can't get the coverage because it's split over
    // parallel paths.  Allow overriding here
    size_t expCoverage,
    // Should we drop variants that would overlap old ones? TODO: we really need
    // a proper system for accounting for usage of graph material.
    bool suppress_overlaps,
    // Should we use average support instead minimum support for our calculations?
    bool useAverageSupport,
    // Should we use the site finder and multiallelic support?
    bool multiallelic_support,
    // What's the max ref length of a site that we genotype as a whole instead
    // of splitting?
    size_t max_ref_length,
    // What's the maximum number of bubble path combinations we can explore
    // while finding one with maximum support?
    size_t max_bubble_paths,
    // what's the minimum minimum allele depth to give a PASS in the filter column
    // (anything below gets FAIL)    
    size_t min_mad_for_filter,
    // print warnings etc. to stderr
    bool verbose) {
    
    vg.paths.sort_by_mapping_rank();
    vg.paths.rebuild_mapping_aux();
    
    if(refPathName.empty()) {
        if (verbose) {
          std:cerr << "Graph has " << vg.paths.size() << " paths to choose from."
                   << std::endl;
        }
        if(vg.paths.size() == 1) {
            // Autodetect the reference path name as the name of the only path
            refPathName = (*vg.paths._paths.begin()).first;
        } else {
            refPathName = "ref";
        }

        if (verbose) {
            std::cerr << "Guessed reference path name of " << refPathName
                      << std::endl;
        }
    }
    
    // Follow the reference path and extract indexes we need: index by node ID,
    // index by node start, and the reconstructed path sequence.
    ReferenceIndex index = trace_reference_path(vg, refPathName, verbose);  
    
    // This holds read support, on each strand, for all the nodes we have read
    // support provided for, by the node pointer in the vg graph.
    std::map<vg::Node*, Support> nodeReadSupport;
    // And read support for the edges
    std::map<vg::Edge*, Support> edgeReadSupport;
    // This maps the likelihood passed from the tsv to the nodes and edges
    // (todo: could save some lookups by lumping with supports)
    std::map<vg::Node*, double> nodeLikelihood;
    std::map<vg::Edge*, double> edgeLikelihood;

    // This holds all the edges that are deletions, by the pointer to the stored
    // Edge object in the VG graph
    std::set<vg::Edge*> deletionEdges;
    
    // This holds where nodes came from (node and offset) in the original, un-
    // augmented graph. For pieces of original nodes, this is where the piece
    // started. For novel nodes, this is where the piece that thois is an
    // alternative to started.
    std::map<vg::Node*, std::pair<int64_t, size_t>> nodeSources;
    
    // We also need to track what edges and nodes are reference (i.e. already
    // known)
    std::set<vg::Node*> knownNodes;
    std::set<vg::Edge*> knownEdges;

    // Parse tsv into an internal format, where we track status and copy number
    // for nodes and edges.
    parse_tsv(glennFile, vg, nodeReadSupport, edgeReadSupport,
              nodeLikelihood, edgeLikelihood, deletionEdges,
              nodeSources, knownNodes, knownEdges, verbose);

    // Store support binned along reference path;
    // Last bin extended to include remainder
    refBinSize = min(refBinSize, index.sequence.size());
    vector<Support> binnedSupport(max(1, int(index.sequence.size() / refBinSize)),
                                  Support(expCoverage / 2, expCoverage /2));
    
    // Crunch the numbers on the reference and its read support. How much read
    // support in total (node length * aligned reads) does the primary path get?
    Support primaryPathTotalSupport = std::make_pair(0.0, 0.0);
    for(auto& pointerAndSupport : nodeReadSupport) {
        if(index.byId.count(pointerAndSupport.first->id())) {
            // This is a primary path node. Add in the total read bases supporting it
            primaryPathTotalSupport += pointerAndSupport.first->sequence().size() * pointerAndSupport.second;
            
            // We also update the total for the appropriate bin
            if (expCoverage == 0) {
                int bin = index.byId[pointerAndSupport.first->id()].first / refBinSize;
                if (bin == binnedSupport.size()) {
                    --bin;
                }
                binnedSupport[bin] = binnedSupport[bin] + 
                    pointerAndSupport.first->sequence().size() * pointerAndSupport.second;
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
        std::cerr << "Primary path average coverage: " << primaryPathAverageSupport << endl;
        std::cerr << "Mininimum binned average coverage: " << binnedSupport[minBin] << " (bin "
                  << (minBin + 1) << " / " << binnedSupport.size() << ")" << endl;
        std::cerr << "Maxinimum binned average coverage: " << binnedSupport[maxBin] << " (bin "
                  << (maxBin + 1) << " / " << binnedSupport.size() << ")" << endl;
    }
    
    // If applicable, load the pileup.
    // This will hold pileup records by node ID.
    std::map<int64_t, vg::NodePileup> nodePileups;
    
    std::function<void(vg::Pileup&)> handlePileup = [&](vg::Pileup& p) { 
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
        std::ifstream in;
        in.open(pileupFilename.c_str());
        stream::for_each(in, handlePileup);
    }
    
    // Generate a vcf header. We can't make Variant records without a
    // VariantCallFile, because the variants need to know which of their
    // available info fields or whatever are defined in the file's header, so
    // they know what to output.
    // Handle length override if specified.
    std::stringstream headerStream;
    write_vcf_header(headerStream, sampleName, contigName,
                     lengthOverride != -1 ? lengthOverride : (index.sequence.size() + variantOffset),
                     min_mad_for_filter);
    
    // Load the headers into a new VCF file object
    vcflib::VariantCallFile vcf;
    std::string headerString = headerStream.str();
    assert(vcf.openForOutput(headerString));
    
    // Spit out the header
    std::cout << headerStream.str();
    
    // Then go through it from the graph's point of view: first over alt nodes
    // backending into the reference (creating things occupying ranges to which
    // we can attribute copy number) and then over reference nodes.

    // We need to track the bases lost.
    size_t basesLost = 0;
    
    if(multiallelic_support) {
        // Do the new thing where we support multiple alleles
    
        // Find all the top-level sites
        std::list<const Snarl*> site_queue;
        
        CactusUltrabubbleFinder finder(vg, refPathName);
        SnarlManager site_manager = finder.find_snarls();
        
        site_manager.for_each_top_level_snarl_parallel([&](const Snarl* site) {
            if(!index.byId.count(site->start().node_id()) || !index.byId.count(site->end().node_id())) {
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
        std::vector<const Snarl*> sites;
        
        while(!site_queue.empty()) {
            // Grab the first site
            const Snarl* site = std::move(site_queue.front());
            site_queue.pop_front();
            
            if (!index.byId.count(site->start().node_id())) {
                #pragma omp critical (cerr)
                cerr << "Off-reference left end: " << pb2json(*site) << endl;
            }
            if (!index.byId.count(site->end().node_id())) {
                #pragma omp critical (cerr)
                cerr << "Off-reference right end: " << pb2json(*site) << endl;
            }
            
            // Sites should not have gotten into the queue without being on the primary path at both ends
            assert(index.byId.count(site->start().node_id()) && index.byId.count(site->end().node_id()));
            
            // Where does the variable region start and end on the reference?
            size_t ref_start = index.byId.at(site->start().node_id()).first +
                vg.get_node(site->start().node_id())->sequence().size();
            size_t ref_end = index.byId.at(site->end().node_id()).first;
            
            if(ref_start > ref_end) {
                // Make sure it's the right way around (it will get set straight
                // in the site when we do the bubbling-up).
                std::swap(ref_start, ref_end);
            }
            
            if(!site_manager.children_of(site).empty() && ref_end > ref_start + max_ref_length) {
                // Site is too big and has children. Split it up and do the
                // children.
                for(const Snarl* child : site_manager.children_of(site)) {
                    // Dump all the children into the queue for separate
                    // processing.
                    
                    if(!index.byId.count(child->start().node_id()) || !index.byId.count(child->end().node_id())) {
                        // Skip child sites that aren't on the reference path at both ends.
                        continue;
                    }
                    
                    site_queue.emplace_back(child);
                }

                if (verbose) {
                    std::cerr << "Broke up site from " << ref_start << " to " << ref_end << " into "
                              << site_manager.children_of(site).size() << " children" << std::endl;
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
                        std::cerr << "Left site from " << ref_start << " to " << ref_end << " with "
                                  << site_manager.children_of(site).size() << " children" << std::endl;
                    }
                }
                
                // Make sure start and end are front-ways relative to the ref path.
                if(index.byId.at(site->start().node_id()).first > index.byId.at(site->end().node_id()).first) {
                    // The site's end happens before its start, flip it around
                    site_manager.flip(site);
                }
                
                // Throw it in the final vector of sites we're going to process.
                sites.push_back(site);
            }                
        }

        if (verbose) {
            std::cerr << "Found " << sites.size() << " sites" << std::endl;
        }
        
        for(const Snarl* site : sites) {
            // For every site, we're going to make a variant
            
            // Get its nodes and edges (including all child sites)
            pair<unordered_set<Node*>, unordered_set<Edge*>> contents = site_manager.deep_contents(site, vg, true);
            
            // Copy its node set
            std::unordered_set<Node*> nodes_left(contents.first);
            
            // Trace the ref path through the site
            std::vector<NodeTraversal> ref_path_for_site;
#ifdef debug
            std::cerr << "Site starts with " << to_node_traversal(site->start(), vg) << " at " << index.byId.at(site->start().node_id()).first
                << " and ends with " << to_node_traversal(site->end(), vg) << " at " << index.byId.at(site->end().node_id()).first << std::endl;
#endif
            
            int64_t ref_node_start = index.byId.at(site->start().node_id()).first;
            while(ref_node_start <= index.byId.at(site->end().node_id()).first) {
            
                // Find the reference node starting here or later.
                auto found = index.byStart.lower_bound(ref_node_start);
                if(found == index.byStart.end()) {
                    throw runtime_error("No ref node found when tracing through site!");
                }
                if((*found).first > index.byId.at(site->end().node_id()).first) {
                    // The next reference node we can find is out of the space
                    // being replaced. We're done.
                    if (verbose) {
                        cerr << "Stopping for out-of-bounds node" << endl;
                    }
                    break;
                }
                
                // Add the traversal to the ref path through the site
                ref_path_for_site.push_back((*found).second);
                
                // Make sure this node is actually in the site
                assert(contents.first.count((*found).second.node));
                
                // Drop it from the set of nodes in the site we haven't visited.
                nodes_left.erase((*found).second.node);
                
                // Next iteration look where this node ends.
                ref_node_start = (*found).first + (*found).second.node->sequence().size();
            }
            
            for(auto node : nodes_left) {
                // Make sure none of the nodes in the site that we didn't visit
                // while tracing along the ref path are on the ref path.
                if(index.byId.count(node->id())) {
                    std::cerr << "Node " << node->id() << " is on ref path at "
                        << index.byId.at(node->id()).first << " but not traced in site "
                        << to_node_traversal(site->start(), vg) << " to " 
                        << to_node_traversal(site->end(), vg) << std::endl;
                    throw runtime_error("Extra ref node found");
                }
            }
            
#ifdef debug
            // Make sure we didn't screw it up
            std::cerr << "Site " << to_node_traversal(site->start(), vg) << " to " << to_node_traversal(site->end(), vg) << ":" << std::endl;
            for(auto& item : ref_path_for_site) {
                std::cerr << "\t" << item << std::endl;
            }
#endif
            assert( ref_path_for_site.front() == to_node_traversal(site->start(), vg) );
            assert( ref_path_for_site.back() == to_node_traversal(site->end(), vg) );

            // check if a path is mostly reference for deciding if we add XREF tag
            auto is_path_known = [&](std::vector<NodeTraversal> path) {
                bool path_is_known = false;
                if(path.size() == 2) {
                    // We just have the anchoring nodes and the edge between them.
                    // Look at that edge specially.
                    vg::Edge* edge = vg.get_edge(make_pair(vg::NodeSide(path[0].node->id(), !path[0].backward),
                                                           vg::NodeSide(path[1].node->id(), path[1].backward)));

                    // The path is known if the edge is known
                    path_is_known = knownEdges.count(edge);
                } else {
                    // How many bases are known (for XREF determination)?
                    assert(path.size() > 2);
                    size_t known_bases = 0;
                    size_t unknown_bases = 0;

                    // don't iterate over referenc anchors at either end
                    for (int i = 1; i < path.size() - 1; ++i) {
                        if (knownNodes.count(path[i].node)) {
                            // This is a reference node.
                            known_bases += path[i].node->sequence().size();
                        } else {
                            unknown_bases += path[i].node->sequence().size();
                        }
                    }
                    // The path is known if more than half of its bases are reference
                    path_is_known = known_bases > 0 && known_bases >= (known_bases + unknown_bases) / 2;
                }
                return path_is_known;
            };
            
            // We need to know all the full-length traversals we're going to
            // consider. We want them in a set so we can find only unique ones.
            // We also store the XREF state here because we want to compute
            // it before the extension is done. 
            std::map<std::vector<NodeTraversal>, bool> site_traversal_set;
            
            // We have this function to extend a partial traversal into a full
            // traversal and add it as a path. The path must already be rooted on
            // the reference in the correct order and orientation.
            auto extend_into_allele = [&](std::vector<NodeTraversal> path) {
                // Splice the ref path through the site and the bubble's path
                // through the site together.
                vector<NodeTraversal> extended_path;

                // Check if the path is xref before adding a bunch of reference stuff to it
                bool path_known = is_path_known(path);
                
                for(auto& traversal : path) {
                    // Make sure the site actually has the nodes we're visiting.
                    assert(contents.first.count(traversal.node));
#ifdef debug
                    if(index.byId.count(traversal.node->id())) {
                        std::cerr << "Path member " << traversal << " lives on ref at "
                        << index.byId.at(traversal.node->id()).first << std::endl;
                    } else {
                        std::cerr << "Path member " << traversal << " does not live on ref" << std::endl;
                    }
#endif
                }
                
                size_t ref_path_index = 0;
                size_t bubble_path_index = 0;
                
                while(ref_path_for_site.at(ref_path_index) != path.at(bubble_path_index)) {
                    // Collect NodeTraversals from the ref path until we hit the one
                    // at which the bubble path starts.
                    extended_path.push_back(ref_path_for_site[ref_path_index++]);
                }
                while(bubble_path_index < path.size()) {
                    // Then take the whole bubble path
                    extended_path.push_back(path[bubble_path_index++]);
                }
                while(ref_path_index < ref_path_for_site.size() && ref_path_for_site.at(ref_path_index) != path.back()) {
                    // Then skip ahead to the matching point in the ref path
                    ref_path_index++;
                }
                if(ref_path_index == ref_path_for_site.size()) {
                    // We ran out of ref path before finding the place to leave the alt.
                    // This must be a backtracking loop or something; start over from the beginning.
                    ref_path_index = 0;
                    
                    while(ref_path_index < ref_path_for_site.size() && ref_path_for_site.at(ref_path_index) != path.back()) {
                        // Then skip ahead to the matching point in the ref path
                        ref_path_index++;
                    }
                    
                    if(ref_path_index == ref_path_for_site.size()) {
                        // Still couldn't find it!
                        std::stringstream err;
                        err << "Couldn't find " << path.back() << " in ref path of site " << to_node_traversal(site->start(), vg) << " to " << to_node_traversal(site->end(), vg) << endl;
                        throw std::runtime_error(err.str());
                    }
                }
                // Skip the matching NodeTraversal
                ref_path_index++;
                while(ref_path_index < ref_path_for_site.size()) {
                    // Then take the entier rest of the ref path
                    extended_path.push_back(ref_path_for_site[ref_path_index++]);
                }
                
                // Now add it to the set
                site_traversal_set.insert(make_pair(extended_path, path_known));
            
            };

            for(Node* node : contents.first) {
                // Find the bubble for each node
                
                if(total(nodeReadSupport.at(node)) == 0) {
                    // Don't bother with unsupported nodes
                    continue;
                }
                
                if(index.byId.count(node->id())) {
                    // Don't try to pathfind to the reference for reference nodes.
                    continue;
                }
                
                std::pair<Support, std::vector<NodeTraversal>> sup_path = find_bubble(
                    vg, node, nullptr, index, nodeReadSupport,
                    edgeReadSupport, maxDepth, max_bubble_paths);

                std::vector<NodeTraversal>& path = sup_path.second;
                
                if(path.empty()) {
                    // We couldn't find a path back to the primary path. Discard
                    // this material.
                    if (verbose) {
                        cerr << "Warning: No path found for node " << node->id() << endl;
                    }
                    basesLost += node->sequence().size();
                    // TODO: what if it's already in another bubble/the node is deleted?
                    continue;
                }
                
                // Extend it out into an allele
                extend_into_allele(path);
                
            }
            
            for(Edge* edge : contents.second) {
                // Go through all the edges
                
                if(!index.byId.count(edge->from()) || !index.byId.count(edge->to())) {
                    // Edge doesn't touch reference at both ends. Don't use it
                    // because for some reason it makes performance worse
                    // overall.
                    continue;
                }
                
                // Find a path based around this edge
                std::pair<Support, std::vector<NodeTraversal>> sup_path = find_bubble(
                    vg, nullptr, edge, index, nodeReadSupport,
                    edgeReadSupport, maxDepth, max_bubble_paths);
                std::vector<NodeTraversal>& path = sup_path.second;
                
#ifdef debug
                std::cerr << "Edge " << edge->from() << " to " << edge->to() << " yields:" << std::endl;
                for(auto& traversal : path) {
                    std::cerr << "\t" << traversal << std::endl;
                }
#endif
                
                if(path.empty()) {
                    // We couldn't find a path back to the primary path. Discard
                    // this material.
                    if (verbose) {
                        cerr << "Warning: No path found for edge " << edge->from() << "," << edge->to() << endl;
                    }
                    // TODO: bases lost
                    // TODO: what if it's already in another bubble/the node is deleted?
                    continue;
                }
                
                // Extend it out into an allele
                extend_into_allele(path);
            }
            
            // Throw the ref allele out of the set
            if(site_traversal_set.count(ref_path_for_site)) {
                site_traversal_set.erase(ref_path_for_site);
            }
            
            // Make it the first in the ordered alleles
            std::vector<std::pair<std::vector<NodeTraversal>, bool>> ordered_paths {make_pair(ref_path_for_site, true)};
            // Then add all the rest
            std::copy(site_traversal_set.begin(), site_traversal_set.end(), std::back_inserter(ordered_paths));
            
            // Collect sequences for all the paths
            std::vector<std::string> sequences;
            // And the lists of involved IDs that we use for variant IDs
            std::vector<string> id_lists;
            // Calculate average and min support for all the alts.
            std::vector<Support> min_supports;
            std::vector<Support> average_supports;
            // And the min likelihood along each path
            std::vector<double> min_likelihoods;
            // Is the path as a whole known or novel?
            std::vector<bool> is_known;
            for(auto& path_xref : ordered_paths) {
                // Go through all the paths
                auto& path = path_xref.first;
                bool path_is_known = path_xref.second;
                
                // We use this to construct the allele sequence
                std::stringstream sequence_stream;
                
                // And this to compose the allele's name in terms of node IDs
                std::stringstream id_stream;
                
                // What's the total support for this path?
                Support total_support = std::make_pair(0.0, 0.0);
                
                // And the min
                Support min_support = std::make_pair(INFINITY, INFINITY);
                
                // Also, what's the min likelihood
                double min_likelihood = INFINITY;
                                
                for(int64_t i = 1; i < path.size() - 1; i++) {
                    // For all but the first and last nodes (which are anchors),
                    // grab their sequences in the correct orientation.
                    
                    std::string added_sequence = path[i].node->sequence();
                
                    if(path[i].backward) {
                        // If the node is traversed backward, we need to flip its sequence.
                        added_sequence = reverse_complement(added_sequence);
                    }
                    
                    // Stick the sequence
                    sequence_stream << added_sequence;
                    
                    // Record ID
                    id_stream << std::to_string(path[i].node->id());
                    if(i != path.size() - 2) {
                        // Add a separator (-2 since the last thing is path is an
                        // anchoring reference node)
                        id_stream << "_";
                    }
                    
                    // How much support do we have for visiting this node?
                    Support node_support = nodeReadSupport.at(path[i].node);
                    // Grab the edge we're traversing into the node
                    vg::Edge* in_edge = vg.get_edge(make_pair(vg::NodeSide(path[i-1].node->id(), !path[i-1].backward),
                                                              vg::NodeSide(path[i].node->id(), path[i].backward)));
                    // If there's less support on the in edge than on the node,
                    // knock it down. We do this separately in each dimension.
                    node_support = support_min(node_support, edgeReadSupport.at(in_edge));
                    
                    // Ditto for the edge we're traversing out of the node
                    vg::Edge* out_edge = vg.get_edge(make_pair(vg::NodeSide(path[i].node->id(), !path[i].backward),
                                                               vg::NodeSide(path[i+1].node->id(), path[i+1].backward)));
                    node_support = support_min(node_support, edgeReadSupport.at(out_edge));
                    
                    
                    // Add support in to the total support for the alt. Scale by node length.
                    total_support += path[i].node->sequence().size() * node_support;
                    min_support = support_min(min_support, node_support);

                    // Update minimum likelihood in the alt path
                    min_likelihood = std::min(min_likelihood, nodeLikelihood.at(path[i].node));
                    
                    // TODO: use edge likelihood here too?
                        
                }
                
                if(path.size() == 2) {
                    // We just have the anchoring nodes and the edge between them.
                    // Look at that edge specially.
                    vg::Edge* edge = vg.get_edge(make_pair(vg::NodeSide(path[0].node->id(), !path[0].backward),
                                                           vg::NodeSide(path[1].node->id(), path[1].backward)));
                    
                    // Only use the support on the edge
                    total_support = edgeReadSupport.at(edge);
                    min_support = total_support;
                    
                    // And the likelihood on the edge
                    min_likelihood = edgeLikelihood.at(edge);
                    
                }
                
                // Fill in the vectors
                sequences.push_back(sequence_stream.str());
                id_lists.push_back(id_stream.str());
                min_supports.push_back(min_support);
                // The support needs to get divided by bases, unless we're just a
                // single edge empty allele, in which case we're special.
                average_supports.push_back(sequences.back().size() > 0 ? total_support / sequences.back().size() : total_support);
                min_likelihoods.push_back(min_likelihood);
                is_known.push_back(path_is_known);
            }
            
            // TODO: complain if multiple copies of the same string exist???
            
            // Decide which support vector we use to actually decide
            std::vector<Support>& supports = useAverageSupport ? average_supports : min_supports;
            
            // Now look at all the paths for the site and pick the top 2.
            int best_allele = -1;
            int second_best_allele = -1;
            
            for(size_t i = 0; i < ordered_paths.size(); i++) {
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
            size_t variation_start = index.byId.at(site->start().node_id()).first
                                     + vg.get_node(site->start().node_id())->sequence().size();
            
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
                std::string extra_base = char_to_string(index.sequence.at(variation_start - 1));
                
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
            auto& genotype = variant.samples[sampleName]["GT"];

            // find which coordinate bin we're in, so we can get the typical local support
            int bin = variation_start / refBinSize;
            if (bin == binnedSupport.size()) {
                --bin;
            }
            const Support& baseline_support = binnedSupport[bin];

            // Decide if we're an indel. We're an indel if the sequence lengths
            // aren't all equal between the ref and the alleles we're going to call.
            bool is_indel = sequences.front().size() == sequences[best_allele].size() &&
                (second_best_allele == -1 || sequences.front().size() == sequences[second_best_allele].size());

            // We need to decide what to scale the bias limits by. We scale them up if this is an indel.
            double bias_multiple = is_indel ? indelBiasMultiple : 1.0;
            
            // How much support do we have for the top two alleles?
            Support site_support = supports.at(best_allele);
            if(second_best_allele != -1) {
                site_support += supports.at(second_best_allele);
            }
            
            // Pull out the different supports. Some of them may be the same.
            Support ref_support = supports.at(0);
            Support best_support = supports.at(best_allele);
            Support second_best_support = std::make_pair(0.0, 0.0);
            if(second_best_allele != -1) {
                second_best_support = supports.at(second_best_allele);
            }
            
            // As we do the genotype, we also compute the likelihood. Holds
            // likelihood log 10. Starts out at "completely wrong".
            double gen_likelihood = -1 * INFINITY;

            // Minimum allele depth of called alleles
            double min_site_support = 0;
            
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
                std::cerr << best_allele << ", " << best_support << " and "
                    << second_best_allele << ", " << second_best_support << std::endl;
                
                std::cerr << bias_limit * bias_multiple * total(second_best_support) << " vs "
                    << total(best_support) << std::endl;
                    
                std::cerr << total(second_best_support) << " vs " << minTotalSupportForCall << std::endl;
#endif
                
                if(second_best_allele != -1 &&
                    bias_limit * bias_multiple * total(second_best_support) >= total(best_support) &&
                    total(best_support) >= minTotalSupportForCall &&
                    total(second_best_support) >= minTotalSupportForCall) {
                    // There's a second best allele, and it's not too biased to
                    // call, and both alleles exceed the minimum to call them
                    // present.
                    genotype.push_back(std::to_string(best_alt) + "/" + std::to_string(second_best_alt));
                    
                    // Compute the likelihood for a best/second best het
                    // Quick quality: combine likelihood and depth, using poisson for latter
                    // TODO: revize which depth (cur: avg) / likelihood (cur: min) pair to use
                    gen_likelihood = log10(poissonp(total(best_support), 0.5 * total(baseline_support))) +
                        log10(poissonp(total(second_best_support), 0.5 * total(baseline_support)));
                    gen_likelihood += min_likelihoods.at(best_allele) + min_likelihoods.at(second_best_allele);
                    // Get minimum support for filter (not assuming it's second_best just to be sure)
                    min_site_support = std::min(total(second_best_support), total(best_support));
                    
                } else if(total(best_support) >= minTotalSupportForCall) {
                    // The second best allele isn't present or isn't good enough,
                    // but the best allele has enough coverage that we can just call
                    // two of it.
                    genotype.push_back(std::to_string(best_alt) + "/" + std::to_string(best_alt));
                    
                    // Compute the likelihood for hom best allele
                    gen_likelihood = log10(poissonp(total(best_support), total(baseline_support)));
                    gen_likelihood += min_likelihoods.at(best_allele);

                    // Get minimum support for filter
                    min_site_support = total(best_support);

                
                } else {
                    // We can't really call this as anything.
                    genotype.push_back("./.");
                }
            } else {
                // Depth too low. Say we have no idea.
                // TODO: elide variant?
                genotype.push_back("./.");
            }
            
            // Now fill in all the other variant info/format stuff

            if((best_allele != 0 && is_known.at(best_allele)) || 
                (second_best_allele != 0 && second_best_allele != -1 && is_known.at(second_best_allele))) {
                // Flag the variant as reference if either of its two best alleles
                // is known but not the primary path. Don't put in a false entry if
                // it isn't known, because vcflib will spit out the flag anyway...
                variant.infoFlags["XREF"] = true;
            }
            
            // Add depth for the variant and the samples
            Support total_support = ref_support;
            if(best_allele != 0) {
                // Add in the depth from the best allele, which is not already counted
                total_support += best_support;
            }
            if(second_best_allele != -1 && second_best_allele != 0) {
                // Add in the depth from the second allele
                total_support += second_best_support;
            }
            std::string depth_string = std::to_string((int64_t)round(total(total_support)));
            variant.format.push_back("DP");
            variant.samples[sampleName]["DP"].push_back(depth_string);
            variant.info["DP"].push_back(depth_string); // We only have one sample, so variant depth = sample depth
            
            // Also allelic depths
            variant.format.push_back("AD");
            variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(ref_support))));
            // And strand biases
            variant.format.push_back("SB");
            variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(ref_support.first)));
            variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(ref_support.second)));
            // Also allelic likelihoods (from minimum values found on their paths)
            variant.format.push_back("AL");
            variant.samples[sampleName]["AL"].push_back(to_string_ss(min_likelihoods.at(0)));
            if(best_allele != 0) {
                // If our best allele isn't ref, it comes next.
                variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(best_support))));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(best_support.first)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(best_support.second)));
                variant.samples[sampleName]["AL"].push_back(to_string_ss(min_likelihoods.at(best_allele)));
            }
            if(second_best_allele != -1 && second_best_allele != 0) {
                // If our second best allele is real and not ref, it comes next.
                variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(second_best_support))));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(second_best_support.first)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(second_best_support.second)));
                variant.samples[sampleName]["AL"].push_back(to_string_ss(min_likelihoods.at(second_best_allele)));
            }
            
            // And total alt allele depth for the alt alleles
            Support alt_support = std::make_pair(0.0, 0.0);
            if(best_allele != 0) {
                alt_support += best_support;
            }
            if(second_best_allele != -1 && second_best_allele != 0) {
                alt_support += second_best_support;
            }
            variant.format.push_back("XAAD");
            variant.samples[sampleName]["XAAD"].push_back(std::to_string((int64_t)round(total(alt_support))));

            // Copy over the quality value
            variant.quality = -10. * log10(1. - pow(10, gen_likelihood));

            // Apply Min Allele Depth cutoff and store result in Filter column
            variant.filter = min_site_support >= min_mad_for_filter ? "PASS" : "FAIL";
            
            if(can_write_alleles(variant)) {
                // No need to check for collisions because we assume sites are correctly found.
                
                // Output the created VCF variant.
                std::cout << variant << std::endl;
            } else {
                if (verbose) {
                    std::cerr << "Variant is too large" << std::endl;
                }
                // TODO: account for the 1 base we added extra if it was a pure
                // insert.
                basesLost += sequences.at(best_allele).size();
            }
        }
        
    } else {
        // Do the old thing with no support for multiallelic sites
        
        // We need to track reference regions taken up by other variants. Holds
        // a map from start position to past_end position. Guaranteed not to have
        // overlaps if used.
        // Only make it nonzero size if we're going to use it.
        IntervalBitfield occupied_regions(suppress_overlaps ? index.sequence.size() : 0);
    
        vg.for_each_node([&](vg::Node* node) {
            // Look at every node in the graph and spit out variants for the ones
            // that are non-reference, but for which we can push copy number to the
            // reference path, greedily.
        
            // Ensure this node is nonreference
            if(index.byId.count(node->id())) {
                // Skip reference nodes
                return;
            }
            
            if(total(nodeReadSupport.at(node)) > 0) {
                // We have copy number on this node.
                
                // Find a path to the primary reference from here
                std::pair<Support, std::vector<NodeTraversal>> sup_path = find_bubble(
                    vg, node, nullptr, index, nodeReadSupport,
                    edgeReadSupport, maxDepth, max_bubble_paths);
                std::vector<NodeTraversal>& path = sup_path.second;
                
                if(path.empty()) {
                    // We couldn't find a path back to the primary path. Discard
                    // this material.
                    basesLost += node->sequence().size();
                    return;
                }
                
                // Turn it into a substitution/insertion
                
                // The position we have stored for this start node is the first
                // position along the reference at which it occurs. Our bubble
                // goes forward in the reference, so we must come out of the
                // opposite end of the node from the one we have stored.
                auto referenceIntervalStart = index.byId.at(path.front().node->id()).first +
                    path.front().node->sequence().size();
                
                // The position we have stored for the end node is the first
                // position it occurs in the reference, and we know we go into
                // it in a reference-concordant direction, so we must have our
                // past-the-end position right there.
                auto referenceIntervalPastEnd = index.byId.at(path.back().node->id()).first;
                
                // We'll fill in this stream with all the node sequences we visit on
                // the path, except for the first and last.
                std::stringstream altStream;
                
                
                if(referenceIntervalPastEnd - referenceIntervalStart == 0) {
                    // If this is an insert, make sure we have the 1 base before it.
                    
                    // TODO: we should handle an insert at the very beginning
                    assert(referenceIntervalStart > 0);
                    
                    // Budge left and add that character to the alt as well
                    referenceIntervalStart--;
                    altStream << index.sequence[referenceIntervalStart];
                }
                
                // Variants should be reference if most of their bases are
                // reference, and novel otherwise. A single SNP base on a known
                // insert should not make it a novel insert.
                size_t knownAltBases = 0;
                
                // How many alt bases are there overall, both known and novel?
                size_t altBases = 0;
                
                // We also need a list of all the alt node IDs for naming the
                // variant.
                std::stringstream idStream;
                
                // And collections of all the involved node pointers, on both the
                // ref and alt sides
                std::set<vg::Node*> refInvolvedNodes;
                std::set<vg::Node*> altInvolvedNodes;
                
                // And we want to know how much read support in total ther alt has
                // (support * node length)
                Support altReadSupportTotal = std::make_pair(0.0, 0.0);
                Support altReadSupportMin = sup_path.first;
                
                // And we keep track of the ref and alt nodes with the lowest likelihoods
                // (for insertion, we use the bypass edge for the ref likelihood)
                std::pair<vg::Node*, double> altMinLikelihood(nullptr, LOG_ZERO);
                std::pair<vg::Node*, double> refMinLikelihood(nullptr, LOG_ZERO);

                for(int64_t i = 1; i < path.size() - 1; i++) {
                    // For all but the first and last nodes, grab their sequences in
                    // the correct orientation.
                    
                    std::string addedSequence = path[i].node->sequence();

                    if(path[i].backward) {
                        // If the node is traversed backward, we need to flip its sequence.
                        addedSequence = vg::reverse_complement(addedSequence);
                    }
                    
                    // Stick the sequence
                    altStream << addedSequence;
                    
                    // Record ID
                    idStream << std::to_string(path[i].node->id());
                    if(i != path.size() - 2) {
                        // Add a separator (-2 since the last thing is path is an
                        // anchoring reference node)
                        idStream << "_";
                    }
                    
                    // Record involvement
                    altInvolvedNodes.insert(path[i].node);

                    Support nodeSupport;
                    if(nodeReadSupport.count(path[i].node)) {
                        // We have read support for this node. 
                        nodeSupport = nodeReadSupport.at(path[i].node);
                    }
                    vg::Edge* pathEdge = vg.get_edge(make_pair(vg::NodeSide(path[i-1].node->id(), !path[i-1].backward),
                                                               vg::NodeSide(path[i].node->id(), path[i].backward)));
                    if (edgeReadSupport.count(pathEdge)) {
                        // We have read support for edge into node, take minimum between it and the node.
                        nodeSupport = support_min(nodeSupport, edgeReadSupport[pathEdge]);
                    } else {
                        throw runtime_error("No support on edge!");
                    }
                    
                    //Add support in to the total support for the alt.
                    altReadSupportTotal += path[i].node->sequence().size() * nodeSupport;

                    // Update minimum likelihood in the alt path
                    if(nodeLikelihood.count(path[i].node)) {
                        double likelihood = nodeLikelihood.at(path[i].node);
                        if (altMinLikelihood.first == nullptr || likelihood < altMinLikelihood.second) {
                            altMinLikelihood = make_pair(path[i].node, likelihood);
                        }
                    }
                        
                    if(knownNodes.count(path[i].node)) {
                        // This is a reference node.
                        knownAltBases += path[i].node->sequence().size();
                    }
                    // We always need to add in the length of the node to the total
                    // length
                    altBases += path[i].node->sequence().size();
                }


                // Find the primary path nodes that are being skipped over/bypassed
                
                // First collect all the IDs of nodes we aren't skipping because
                // they're in the alt. Don't count those.
                std::set<int64_t> altIds;
                for(auto& visit : path) {
                    altIds.insert(visit.node->id());
                }

                // Holds total primary path base readings observed (read support * node length).
                Support refReadSupportTotal = std::make_pair(0.0, 0.0);
                Support refReadSupportMin = std::make_pair(std::numeric_limits<double>::max(),
                                                           std::numeric_limits<double>::max());
                
                // And total bases of material we looked at on the primary path and
                // not this alt
                size_t refBases = 0;
                int64_t refNodeStart = referenceIntervalStart;
                
                while(refNodeStart < referenceIntervalPastEnd) {
                
                    // Find the reference node starting here or later. Remember that
                    // a variant anchored at its left base to a reference position
                    // may have no node starting right where it starts.
                    auto found = index.byStart.lower_bound(refNodeStart);
                    if(found == index.byStart.end()) {
                        // No reference nodes here! That's a bit weird. But stop the
                        // loop.
                        break;
                    }
                    if((*found).first >= referenceIntervalPastEnd) {
                        // The next reference node we can find is out of the space
                        // being replaced. We're done.
                        break;
                    }
                    
                    // Pull out the reference node we located
                    auto* refNode = (*found).second.node;
                    
                    // Record involvement
                    refInvolvedNodes.insert(refNode);
                
                    // Next iteration look where this node ends.
                    refNodeStart = (*found).first + refNode->sequence().size();
                
                    if(altIds.count(refNode->id())) {
                        // This node is also involved in the alt we did take, so
                        // skip it. TODO: work out how to deal with shared nodes.
#ifdef debug
                        std::cerr << "Node " << refNode->id() << " also used in alt" << std::endl;
#endif
                        continue;
                    }
                    
                    // Say we saw these bases, which may or may not have been called present
                    refBases += refNode->sequence().size();
#ifdef debug
                    std::cerr << "Node " << refNode->id() << " has " << nodeReadSupport.at(refNode) << " copies" << std::endl;
#endif
                    
                    // Count the bases we see not deleted
                    refReadSupportTotal += refNode->sequence().size() * nodeReadSupport.at(refNode);
                    refReadSupportMin = support_min(refReadSupportMin, nodeReadSupport.at(refNode));

                    // Update minimum likelihood in the ref path
                    if(nodeLikelihood.count(refNode)) {
                        double likelihood = nodeLikelihood.at(refNode);
                        if (refMinLikelihood.first == nullptr || likelihood < refMinLikelihood.second) {
                            refMinLikelihood = make_pair(refNode, likelihood);
                        }
                    }

                }
                
#ifdef debug
                std::cerr << idStream.str() << " ref alternative: " << refReadSupportTotal << "/" << refBases
                    << " from " << referenceIntervalStart << " to " << referenceIntervalPastEnd << std::endl;
#endif

                // We divide the read support of stuff passed over by the total
                // bases of stuff passed over to get the average read support for
                // the primary path allele.
                Support refReadSupportAverage = refBases == 0 ? std::make_pair(0.0, 0.0) : refReadSupportTotal / refBases;
                
                // And similarly for the alt
                Support altReadSupportAverage = altBases == 0 ? std::make_pair(0.0, 0.0) : altReadSupportTotal / altBases;

                if(refBases == 0) {
                    // There's no reference node; we're a pure insert. Like with
                    // deletions, we should look at the edge that bypasses us to see
                    // if there's any support for it.
                    
                    // We eant an edge from the end of the node before us to the
                    // start of the node after us.
                    std::pair<vg::NodeSide, vg::NodeSide> edgeWanted = std::make_pair(
                        vg::NodeSide(path.front().node->id(), true),
                        vg::NodeSide(path.back().node->id()));
                    
                    if(vg.has_edge(edgeWanted)) {
                        // We found it!
                        vg::Edge* bypass = vg.get_edge(edgeWanted);
                        
                        // Any reads supporting the edge bypassing the insert are
                        // really ref support reads, and should count as supporting
                        // the whole ref allele.
                        refReadSupportTotal = edgeReadSupport.count(bypass) ? edgeReadSupport.at(bypass) : std::make_pair(0.0, 0.0); 
                        refReadSupportAverage = refReadSupportTotal;
                        refReadSupportMin = refReadSupportTotal;

                        // set minimum likelihood in the ref path using the edge
                        if(edgeLikelihood.count(bypass)) {
                            double likelihood = edgeLikelihood.at(bypass);
                            assert(refMinLikelihood.first == nullptr);
                            refMinLikelihood = make_pair(nullptr, likelihood);
                        }
                    }
                }
                // Otherwise if there's no edge or no support for that edge, the ref support should stay 0.
                
                // Make the variant and emit it.
                std::string refAllele = index.sequence.substr(
                    referenceIntervalStart, referenceIntervalPastEnd - referenceIntervalStart);
                std::string altAllele = altStream.str();
                
                // Make a Variant
                vcflib::Variant variant;
                variant.sequenceName = contigName;
                variant.setVariantCallFile(vcf);
                variant.quality = 0;
                variant.position = referenceIntervalStart + 1 + variantOffset;
                variant.id = idStream.str();
                
                if(knownAltBases > 0 && knownAltBases >= altBases / 2) {
                    // Flag the variant as reference. Don't put in a false entry if
                    // it isn't, because vcflib will spit out the flag anyway...
                    variant.infoFlags["XREF"] = true;
                }
                
                // We need to deduplicate the corss-references, because multiple
                // involved nodes may cross-reference the same original node and
                // offset, and because the same original node and offset can be
                // referenced by both ref and alt paths.
                std::set<std::pair<int64_t, size_t>> refCrossreferences;
                std::set<std::pair<int64_t, size_t>> altCrossreferences;
                std::set<std::pair<int64_t, size_t>> crossreferences;
                
                for(auto* node : refInvolvedNodes) {
                    // Every involved node gets its original node:offset recorded as
                    // an XSEE cross-reference.
                    
                    if(nodeSources.count(node)) {
                        // It has a source. Find it
                        auto& source = nodeSources.at(node);
                        // Then add it to be referenced.
                        refCrossreferences.insert(source);
                        crossreferences.insert(source);
                    }
                }
                for(auto* node : altInvolvedNodes) {
                    // Every involved node gets its original node:offset recorded as
                    // an XSEE cross-reference.
                    
                    if(nodeSources.count(node)) {
                        // It has a source. Find it
                        auto& source = nodeSources.at(node);
                        // Then add it to be referenced.
                        altCrossreferences.insert(source);
                        crossreferences.insert(source);
                    }
                }
                for(auto& crossreference : crossreferences) {
                    variant.info["XSEE"].push_back(std::to_string(crossreference.first) + ":" +
                        std::to_string(crossreference.second));
                }
                
                // Initialize the ref allele
                create_ref_allele(variant, refAllele);
                
                // Add the alt allele
                int altNumber = add_alt_allele(variant, altAllele);
                
                // Say we're going to spit out the genotype for this sample.        
                variant.format.push_back("GT");
                auto& genotype = variant.samples[sampleName]["GT"];

                // find which bin we're in
                int bin = referenceIntervalStart / refBinSize;
                if (bin == binnedSupport.size()) {
                    --bin;
                }
                const Support& baselineSupport = binnedSupport[bin];

                Support refSupport = useAverageSupport ? refReadSupportAverage : refReadSupportMin;
                Support altSupport = useAverageSupport ? altSupport : altReadSupportMin;
                
                // We need to decide what to scale the bias limits by. We scale them up if this is an indel.
                double biasMultiple = (altAllele.size() == refAllele.size()) ? 1.0 : indelBiasMultiple;
                
                // We're going to make some really bad calls at low depth. We can
                // pull them out with a depth filter, but for now just elide them.
                if(total(refSupport + altSupport) >= total(baselineSupport) * minFractionForCall) {
                    if(total(refSupport) > maxRefHetBias * biasMultiple * total(altSupport) &&
                        total(refReadSupportTotal) >= minTotalSupportForCall) {
                        // Biased enough towards ref, and ref has enough total reads.
                        // Say it's hom ref
                        genotype.push_back("0/0");
                    } else if(total(altSupport) > maxHetBias * biasMultiple * total(refSupport)
                        && total(altReadSupportTotal) >= minTotalSupportForCall) {
                        // Say it's hom alt
                        genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
                    } else if(total(refReadSupportTotal) >= minTotalSupportForCall &&
                        total(altReadSupportTotal) >= minTotalSupportForCall) {
                        // Say it's het
                        genotype.push_back("0/" + std::to_string(altNumber));
                    } else {
                        // We can't really call this as anything.
                        genotype.push_back("./.");
                    }
                } else {
                    // Depth too low. Say we have no idea.
                    // TODO: elide variant?
                    genotype.push_back("./.");
                }
                // TODO: use legit thresholds here.
                
                // Add depth for the variant and the samples
                std::string depthString = std::to_string((int64_t)round(total(refSupport + altSupport)));
                variant.format.push_back("DP");
                variant.samples[sampleName]["DP"].push_back(depthString);
                variant.info["DP"].push_back(depthString); // We only have one sample, so variant depth = sample depth
                
                // Also allelic depths
                variant.format.push_back("AD");
                variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(refSupport))));
                variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(altSupport))));
                
                // Also strand biases
                variant.format.push_back("SB");
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(refSupport.first)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(refSupport.second)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(altSupport.first)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(altSupport.second)));

                // And total alt allele depth
                variant.format.push_back("XAAD");
                variant.samples[sampleName]["XAAD"].push_back(std::to_string((int64_t)round(total(altSupport))));

                // Also allelic likelihoods (from minimum values found on their paths)
                variant.format.push_back("AL");
                variant.samples[sampleName]["AL"].push_back(to_string_ss(refMinLikelihood.second));
                variant.samples[sampleName]["AL"].push_back(to_string_ss(altMinLikelihood.second));

                // Quick quality: combine likelihood and depth, using poisson for latter
                // todo: revize which depth (cur: avg) / likelihood (cur: min) pair to use
                double genLikelihood;
                double min_site_support = 0;
                if (genotype.back() == "0/0") {
                    genLikelihood = log10(poissonp(total(refSupport), total(baselineSupport)));
                    genLikelihood += refMinLikelihood.second;
                    min_site_support =  total(refSupport);
                } else if (genotype.back() == "1/1") {
                    genLikelihood = log10(poissonp(total(altSupport), total(baselineSupport)));
                    genLikelihood += altMinLikelihood.second;
                    min_site_support = total(altSupport);
                } else {
                    genLikelihood = log10(poissonp(total(refSupport), 0.5 * total(baselineSupport))) +
                        log10(poissonp(total(altSupport), 0.5 * total(baselineSupport)));
                    genLikelihood += refMinLikelihood.second + altMinLikelihood.second;
                    min_site_support = std::min(total(refSupport), total(altSupport));
                }
                variant.quality = -10. * log10(1. - pow(10, genLikelihood));

                // Apply Min Allele Depth cutoff and store result in Filter column
                variant.filter = min_site_support >= min_mad_for_filter ? "PASS" : "FAIL";
                
                
#ifdef debug
                std::cerr << "Found variant " << refAllele << " -> " << altAllele
                    << " caused by nodes " <<  variant.id
                    << " at 1-based reference position " << variant.position
                    << std::endl;
#endif

                if(can_write_alleles(variant)) {
                    if(!suppress_overlaps || !occupied_regions.collides(referenceIntervalStart, referenceIntervalPastEnd)) {
                        // Variant doesn't intersect with something we care about
                        
                        // Output the created VCF variant.
                        std::cout << variant << std::endl;
                        
                        // Output the pileup line, which will be nonempty if we have pileups
                        std::cout << get_pileup_line(nodePileups, refCrossreferences, altCrossreferences);
                    
                        if(suppress_overlaps) {
                            // Mark region as occupied
                            occupied_regions.add(referenceIntervalStart, referenceIntervalPastEnd);
                        }
                    
                    } else {
                        if (verbose) {
                            std::cerr << "Variant collides with already-emitted variant" << std::endl;
                        }
                        basesLost += altAllele.size();
                    }
                } else {
                    if (verbose) {
                        std::cerr << "Variant is too large" << std::endl;
                    }
                    // TODO: account for the 1 base we added extra if it was a pure
                    // insert.
                    basesLost += altAllele.size();
                }
                
                
            }
            
            
        });
        
        for(vg::Edge* deletion : deletionEdges) {
            // Make deletion variants for each deletion edge
            
            // Make a string naming the edge
            std::string edgeName = std::to_string(deletion->from()) +
                (deletion->from_start() ? "L" : "R") + "->" +
                std::to_string(deletion->to()) + (deletion->to_end() ? "R" : "L");
            
            if(!index.byId.count(deletion->from()) || !index.byId.count(deletion->to())) {
                // This deletion edge does not cover a reference interval.
                // TODO: take into account its presence when pushing copy number.
#ifdef debug
                std::cerr << "Deletion edge " << edgeName << " does not cover a reference interval. Skipping!" << endl;
#endif
                continue;
            }
            
            // Where are we from and to in the reference (leftmost position and
            // relative orientation)
            auto& fromPlacement = index.byId.at(deletion->from());
            auto& toPlacement = index.byId.at(deletion->to());
            
#ifdef debug
            std::cerr << "Node " << deletion->from() << " is at ref position " 
                << fromPlacement.first << " orientation " << fromPlacement.second << std::endl;
            std::cerr << "Node " << deletion->to() << " is at ref position "
                << toPlacement.first << " orientation " << toPlacement.second << std::endl;
#endif
            
            // Are we attached to the reference-relative left or right of our from
            // base?
            bool fromFirst = fromPlacement.second != deletion->from_start();
            
            // And our to base?
            bool toLast = toPlacement.second != deletion->to_end();
            
            // What base should the from end really be on? This is the non-deleted
            // base outside the deletion on the from end.
            int64_t fromBase = fromPlacement.first + (fromFirst ? 0 : vg.get_node(deletion->from())->sequence().size() - 1);
            
            // And the to end?
            int64_t toBase = toPlacement.first + (toLast ? vg.get_node(deletion->to())->sequence().size() - 1 : 0);

            if(toBase <= fromBase) {
                // Our edge ought to be running backward.
                if(!(fromFirst && toLast)) {
                    // We're not a proper deletion edge in the backwards spelling
                    // Discard the edge
                    if (verbose) {
                        std::cerr << "Improper deletion edge " << edgeName << std::endl;
                    }
                    basesLost += toBase - fromBase;
                    continue;
                } else {
                    // Just invert the from and to bases.
                    std::swap(fromBase, toBase);
#ifdef debug
                    std::cerr << "Inverted deletion edge " << edgeName << std::endl;
#endif
                }
            } else if(fromFirst || toLast) {
                // We aren't a proper deletion edge in the forward spelling either.
                if (verbose) {
                    std::cerr << "Improper deletion edge " << edgeName << std::endl;
                }
                basesLost += fromBase - toBase;
                continue;
            }
            
            
            if(toBase <= fromBase + 1) {
                // No bases were actually deleted. Maybe this is just a normal reference edge.
                continue;
            }
            
#ifdef debug
            std::cerr << "Deletion " << edgeName << " bookended by " << fromBase << " and " << toBase << std::endl;
#endif
            
            // What original node:offset places do we care about?
            std::set<std::pair<int64_t, size_t>> crossreferences;
            
            // Guess the copy number of the deletion.
            // Holds total base copies (node length * read support) observed as not deleted.
            Support refReadSupportTotal = std::make_pair(0.0, 0.0);
            Support refReadSupportMin = std::make_pair(std::numeric_limits<double>::max(),
                                                       std::numeric_limits<double>::max());
            std::pair<vg::Node*, double> refMinLikelihood(NULL, LOG_ZERO);
            
            int64_t deletedNodeStart = fromBase + 1;
            while(deletedNodeStart != toBase) {
#ifdef debug
                std::cerr << "Next deleted node starts at " << deletedNodeStart << std::endl;
#endif
            
                // Find the deleted node starting here in the reference
                auto* deletedNode = index.byStart.at(deletedNodeStart).node;
                // We know the next reference node should start just after this one.
                // Even if it previously existed in the reference.
                deletedNodeStart += deletedNode->sequence().size();
                
                // Count the read observations we see not deleted
                refReadSupportTotal += deletedNode->sequence().size() * nodeReadSupport.at(deletedNode);
                refReadSupportMin = support_min(refReadSupportMin, nodeReadSupport.at(deletedNode));

                // Update minimum node likelihood
                double likelihood = nodeLikelihood.at(deletedNode);
                if (refMinLikelihood.first == NULL || likelihood < refMinLikelihood.second) {
                    refMinLikelihood = make_pair(deletedNode, likelihood);
                }
                
                if(nodeSources.count(deletedNode)) {
                    // Add the beginning of this node as a see also. The deletion is
                    // stored at the beginning of the first node, and the other
                    // locations will help us get an idea of how much support the
                    // deleted stuff has.
                    crossreferences.insert(nodeSources.at(deletedNode));
                }
            }
            
            // We divide the total read support by the total bases deleted to get
            // the average read support for the reference.
            Support refReadSupportAverage = refReadSupportTotal / (toBase - fromBase - 1);
            
            // Get the support for the edge itself?
            Support altReadSupportTotal = edgeReadSupport.count(deletion) ? edgeReadSupport[deletion] : std::make_pair(0.0, 0.0);
            double altMinLikelihood = edgeLikelihood.count(deletion) ? edgeLikelihood[deletion] : LOG_ZERO;
            // No sense averaging the deletion edge read support because there are no bases.
            
            // What copy number do we call for the deletion?
            int64_t copyNumberCall = 2;

            // Now we know fromBase is the last non-deleted base and toBase is the
            // first non-deleted base. We'll make an alt replacing the first non-
            // deleted base plus the deletion with just the first non-deleted base.
            // Rename everything to the same names we were using before.
            size_t referenceIntervalStart = fromBase;
            size_t referenceIntervalPastEnd = toBase;
     
            // find which bin we're in
            int bin = referenceIntervalStart / refBinSize;
            if (bin == binnedSupport.size()) {
                --bin;
            }
            const Support& baselineSupport = binnedSupport[bin];

            Support refSupport = useAverageSupport ? refReadSupportAverage : refReadSupportMin;
            // altSupport comes from single edge, so will always just be altReadSupportTotal

            // These are always indels, so no need to decide on a variable bias multiple

            // We're going to make some really bad calls at low depth. We can
            // pull them out with a depth filter, but for now just elide them.
            if(total(refSupport + altReadSupportTotal) >= total(baselineSupport) * minFractionForCall) {
                if(total(refSupport) > maxRefHetBias * indelBiasMultiple * total(altReadSupportTotal) &&
                    total(refReadSupportTotal) >= minTotalSupportForCall) {
                    // Say it's hom ref
                    copyNumberCall = 0;
                } else if(total(altReadSupportTotal) > maxHetBias * indelBiasMultiple * total(refSupport) &&
                    total(altReadSupportTotal) >= minTotalSupportForCall) {
                    // Say it's hom alt
                    copyNumberCall = 2;
                } else if(total(refReadSupportTotal) >= minTotalSupportForCall &&
                    total(altReadSupportTotal) >= minTotalSupportForCall) {
                    // Say it's het
                    copyNumberCall = 1;
                } else {
                    // We're not biased enough towards either homozygote, but we
                    // don't have enough support for each allele to call het.
                    // TODO: we just don't call
                    copyNumberCall = 0;
                }
            } else {
                // Depth too low. Don't call.
                copyNumberCall = 0;
            }
            // TODO: use legit thresholds here.
            
            
            if(copyNumberCall == 0) {
                // Actually don't call a deletion
                continue;
            }
                   
            // Make the variant and emit it.
            std::string refAllele = index.sequence.substr(
                referenceIntervalStart, referenceIntervalPastEnd - referenceIntervalStart);
            std::string altAllele = index.sequence.substr(referenceIntervalStart, 1);
            
            // Make a Variant
            vcflib::Variant variant;
            variant.sequenceName = contigName;
            variant.setVariantCallFile(vcf);
            variant.quality = 0;
            variant.position = referenceIntervalStart + 1 + variantOffset;
            variant.id = edgeName;
            
            if(knownEdges.count(deletion)) {
                // Mark it as reference if it is a reference edge. Apparently vcflib
                // does not check the flag value when serializing, so don't put in a
                // flase entry if it's not a reference edge.œ
                variant.infoFlags["XREF"] = true;
#ifdef debug
                std::cerr << edgeName << " is a reference deletion" << std::endl;
#endif
            }
            
            for(auto& crossreference : crossreferences) {
                // Add in all the deduplicated cross-references
                variant.info["XSEE"].push_back(std::to_string(crossreference.first) + ":" +
                    std::to_string(crossreference.second));
            }
            
            // Initialize the ref allele
            create_ref_allele(variant, refAllele);
            
            // Add the alt allele
            int altNumber = add_alt_allele(variant, altAllele);
            
            // Say we're going to spit out the genotype for this sample.        
            variant.format.push_back("GT");
            auto& genotype = variant.samples[sampleName]["GT"];
            
            if(copyNumberCall == 1) {
                // We're allele alt and ref heterozygous.
                genotype.push_back("0/" + std::to_string(altNumber));
            } else if(copyNumberCall == 2) {
                // We're alt homozygous, other overlapping variants notwithstanding.
                genotype.push_back(std::to_string(altNumber) + "/" + std::to_string(altNumber));
            } else {
                // We're something weird
                throw std::runtime_error("Invalid copy number for deletion: " + std::to_string(copyNumberCall));
            }
            
            // Add depth for the variant and the samples
                std::string depthString = std::to_string((int64_t)round(total(refSupport + altReadSupportTotal)));
                variant.format.push_back("DP");
                variant.samples[sampleName]["DP"].push_back(depthString);
                variant.info["DP"].push_back(depthString); // We only have one sample, so variant depth = sample depth
                
                // Also allelic depths
                variant.format.push_back("AD");
                variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(refSupport))));
                variant.samples[sampleName]["AD"].push_back(std::to_string((int64_t)round(total(altReadSupportTotal))));
                
                // Also strand biases
                variant.format.push_back("SB");
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(refSupport.first)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(refSupport.second)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(altReadSupportTotal.first)));
                variant.samples[sampleName]["SB"].push_back(std::to_string((int64_t)round(altReadSupportTotal.second)));
                
                // And total alt allele depth
                variant.format.push_back("XAAD");
                variant.samples[sampleName]["XAAD"].push_back(std::to_string((int64_t)round(total(altReadSupportTotal))));

                // Also allelic likelihoods (from minimum values found on their paths)
                variant.format.push_back("AL");
                variant.samples[sampleName]["AL"].push_back(to_string_ss(refMinLikelihood.second));
                variant.samples[sampleName]["AL"].push_back(to_string_ss(altMinLikelihood));

                // Quick quality: combine likelihood and depth, using poisson for latter
                // todo: revize which depth (cur: avg) / likelihood (cur: min) pair to use
                double genLikelihood;
                double min_site_support = 0;
                if (genotype.back() == "0/0") {
                    genLikelihood = log10(poissonp(total(refSupport), total(baselineSupport)));
                    genLikelihood += refMinLikelihood.second;
                    min_site_support = total(refSupport);
                } else if (genotype.back() == "1/1") {
                    genLikelihood = log10(poissonp(total(altReadSupportTotal), total(baselineSupport)));
                    genLikelihood += altMinLikelihood;
                    min_site_support = total(altReadSupportTotal);
                } else {
                    genLikelihood = log10(poissonp(total(refSupport), 0.5 * total(baselineSupport))) +
                        log10(poissonp(total(altReadSupportTotal), 0.5 * total(baselineSupport)));
                    genLikelihood += refMinLikelihood.second + altMinLikelihood;
                    min_site_support = std::min(total(refSupport), total(altReadSupportTotal));
                }
                variant.quality = -10. * log10(1. - pow(10, genLikelihood));

                // Apply Min Allele Depth cutoff and store result in Filter column
                variant.filter = min_site_support >= min_mad_for_filter ? "PASS" : "FAIL";
                
#ifdef debug
            std::cerr << "Found variant " << refAllele << " -> " << altAllele
                << " caused by edge " <<  variant.id
                << " at 1-based reference position " << variant.position
                << std::endl;
#endif

            if(can_write_alleles(variant)) {
                if(!suppress_overlaps || !occupied_regions.collides(referenceIntervalStart, referenceIntervalPastEnd)) {
                    // Variant doesn't intersect with something we care about
                    
                    // Output the created VCF variant.
                    std::cout << variant << std::endl;
                    
                    // Output the pileup line, which will be nonempty if we have pileups
                    // We only have ref crossreferences here. TODO: make the ref and alt
                    // labels make sense for deletions/re-design the way labeling works.
                    std::cout << get_pileup_line(nodePileups, crossreferences, std::set<std::pair<int64_t, size_t>>());
                
                    if(suppress_overlaps) {
                        // Mark region as occupied
                        occupied_regions.add(referenceIntervalStart, referenceIntervalPastEnd);
                    }
                
                } else {
                    if (verbose) {
                        std::cerr << "Variant collides with already-emitted variant" << std::endl;
                    }
                    basesLost += altAllele.size();
                }
            } else {
                if (verbose) {
                    std::cerr << "Variant is too large" << std::endl;
                }
                // TODO: Drop the anchoring base that doesn't really belong to the
                // deletion, when we can be consistent with inserts.
                basesLost += altAllele.size();
            }
            
        }
        
    }
    
    // Announce how much we can't show.
    if (verbose) {
        std::cerr << "Had to drop " << basesLost << " bp of unrepresentable variation." << std::endl;
    }
    
    return 0;
}

}
