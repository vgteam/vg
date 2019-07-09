#ifndef VG_NODETRAVERSAL_HPP_INCLUDED
#define VG_NODETRAVERSAL_HPP_INCLUDED

#include <vg/vg.pb.h>
#include "hash_map.hpp"


namespace vg {

using namespace std;

/// Represents a node traversed in a certain orientation. The default orientation
/// is start to end, but if `backward` is set, represents the node being
/// traversed end to start. A list of these can serve as an edit-free version of
/// a path, especially if supplemented with a length and an initial node offset.
/// A path node has a left and a right side, which are the start and end of the
/// node if it is forward, or the end and start of the node if it is backward.
class NodeTraversal {
public:
    /// What Node is being traversed?
    Node* node;
    /// In what orientation is it being traversed?
    bool backward;

    
    /// Make a NodeTraversal that traverses the given Node in the given orientation.
    /// We don't want Node*s to turn into NodeTraversals when we aren't expecting it, so this is explicit.
    explicit inline NodeTraversal(Node* node, bool backward = false): node(node), backward(backward) {
        // Nothing to do
    }

    /// Create a NodeTraversal of no node.
    inline NodeTraversal(): NodeTraversal(nullptr) {
        // Nothing to do
    }

    /// Equality operator.
    inline bool operator==(const NodeTraversal& other) const {
        return node == other.node && backward == other.backward;
    }

    /// Inequality operator.
    inline bool operator!=(const NodeTraversal& other) const {
        return node != other.node || backward != other.backward;
    }

    /// Comparison operator for sorting in sets and maps.
    /// Make sure to sort by node ID and not pointer value, because people will expect that.
    inline bool operator<(const NodeTraversal& other) const {
        if(node == nullptr && other.node != nullptr) {
            // We might have a null node when they don't.
            return true;
        }
        if(other.node == nullptr && node != nullptr) {
            // They might have a null node when we don't.
            return false;
        }
        if(other.node == nullptr && node == nullptr) {
            // We bith might have null nodes.
            return backward < other.backward;
        }
        // Now we know none of the nodes are null. Sort by actual ID.
        return node->id() < other.node->id() || (node == other.node && backward < other.backward);
    }

    /// Reverse complement the node traversal, returning a traversal of the same node in the opposite direction.
    inline NodeTraversal reverse(void) const {
        return NodeTraversal(node, !backward);
    }
};

/// Print the given NodeTraversal.
inline ostream& operator<<(ostream& out, const NodeTraversal& nodetraversal) {
    return out << (nodetraversal.node ? nodetraversal.node->id() : (int64_t) 0)  << " " << (nodetraversal.backward ? "rev" : "fwd");
}

/// hash function for NodeTraversals
template<>
struct wang_hash<NodeTraversal> {
    size_t operator()(const NodeTraversal& x) const {
        return wang_hash<pair<vg::Node*, bool>>()(make_pair(x.node, x.backward));
    }
};

}   // namespace vg

namespace std {
    /// hash function for NodeTraversals
    template<>
    struct hash<vg::NodeTraversal> {
        size_t operator()(const vg::NodeTraversal& trav) const {
            return hash<pair<vg::Node*, bool > >()(make_pair(trav.node, trav.backward));
        }
    };
}


#endif
