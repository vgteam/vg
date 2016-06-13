#ifndef VG_NODETRAVERSAL_HPP
#define VG_NODETRAVERSAL_HPP

#include "vg.pb.h"

namespace vg {

using namespace std;

// Represents a node traversed in a certain orientation. The default orientation
// is start to end, but if `backward` is set, represents the node being
// traversed end to start. A list of these can serve as an edit-free version of
// a path, especially if supplemented with a length and an initial node offset.
// A path node has a left and a right side, which are the start and end of the
// node if it is forward, or the end and start of the node if it is backward.
class NodeTraversal {
public:
    Node* node;
    bool backward;

    inline NodeTraversal(Node* node, bool backward = false): node(node), backward(backward) {
        // Nothing to do
    }

    inline NodeTraversal(): NodeTraversal(nullptr) {
        // Nothing to do
    }

    inline bool operator==(const NodeTraversal& other) const {
        return node == other.node && backward == other.backward;
    }

    inline bool operator!=(const NodeTraversal& other) const {
        return node != other.node || backward != other.backward;
    }

    inline bool operator<(const NodeTraversal& other) const {
        return node < other.node || (node == other.node && backward < other.backward);
    }

    // reverse complement the node traversal
    inline NodeTraversal reverse(void) const {
        return NodeTraversal(node, !backward);
    }

};

}

#endif
