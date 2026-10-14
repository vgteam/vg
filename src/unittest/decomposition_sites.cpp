/// \file decomposition_sites.cpp
///
/// Unit tests for the decomposition_sites functions, which must give the same
/// sites on a SnarlManagerDecomposition and on a SnarlDistanceIndex built
/// from the same decomposition.

#include "catch.hpp"
#include "../decomposition_sites.hpp"
#include "../handle.hpp"
#include "../vg.hpp"
#include "../integrated_snarl_finder.hpp"
#include "support/decomposition_bench.hpp"
#include "support/decomposition_comparison.hpp"
#include "support/randomness.hpp"

#include <map>
#include <string>
#include <vector>

namespace vg::unittest {

/**
 * @brief Require that the bounds of a site's or chain's traversals agree with
 * their flips.
 *
 * @post For each realizable traversal that starts and ends at bounds, the
 * flipped traversal's start bound is the traversal's end bound read the
 * other way, and its end bound is the traversal's start bound read the other
 * way.
 */
static void require_flipped_bounds(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                   const net_handle_t& net) {
    decomposition.for_each_traversal(net, [&decomposition, &graph](const net_handle_t& traversal) {
        if (decomposition.starts_at_tip(traversal) || decomposition.ends_at_tip(traversal)) {
            return;
        }
        SiteBounds bounds = site_bounds(decomposition, graph, traversal);
        SiteBounds flipped = site_bounds(decomposition, graph, decomposition.flip(traversal));
        REQUIRE(flipped.start == graph.flip(bounds.end));
        REQUIRE(flipped.end == graph.flip(bounds.start));
    });
}

/**
 * @brief Require that a site and every site under it fit together.
 *
 * @post The site, its child chains and the sites under it satisfy
 * require_flipped_bounds().
 * @post Each child site has the site as its parent, and is one deeper.
 */
static void require_site_subtree(const SnarlDecomposition& decomposition, const HandleGraph& graph,
                                 const net_handle_t& site) {
    TraversalLabeler labeler(decomposition, graph);
    require_flipped_bounds(decomposition, graph, site);
    for (const net_handle_t& chain : child_chains(decomposition, site)) {
        require_flipped_bounds(decomposition, graph, chain);
        for (const net_handle_t& child : child_sites(decomposition, chain)) {
            REQUIRE(labeler.bounds_label(parent_site(decomposition, child)) == labeler.bounds_label(site));
            REQUIRE(depth(decomposition, child) == depth(decomposition, site) + 1);
            require_site_subtree(decomposition, graph, child);
        }
    }
}

/**
 * @brief Require that the sites of a decomposition fit together.
 *
 * @post Each top-level site has depth 0, has a parent for which is_root() is
 * true, and satisfies require_site_subtree().
 */
static void require_site_tree(const SnarlDecomposition& decomposition, const HandleGraph& graph) {
    for (const net_handle_t& site : top_level_sites(decomposition)) {
        REQUIRE(depth(decomposition, site) == 0);
        REQUIRE(decomposition.is_root(parent_site(decomposition, site)));
        require_site_subtree(decomposition, graph, site);
    }
}

/// @brief A summary of one site's place in the tree of sites, in labels that
/// compare between decompositions.
struct SiteFacts {
    size_t depth;
    bool is_leaf;
    /// @brief The bounds label of the parent site, or "root".
    std::string parent;
    /// @brief For each child chain, its bounds label and the bounds labels of
    /// its sites in order.
    std::vector<std::pair<std::string, std::vector<std::string>>> child_chains;

    bool operator==(const SiteFacts& other) const = default;
};

/// @brief Sites of one decomposition, by bounds label.
using SiteFactsByLabel = std::map<std::string, SiteFacts, std::less<>>;

/// @brief Find the facts about a site and every site under it.
static void collect_site_facts(const SnarlDecomposition& decomposition, const TraversalLabeler& labeler,
                               const net_handle_t& site, SiteFactsByLabel& found) {
    net_handle_t parent = parent_site(decomposition, site);
    SiteFacts facts {depth(decomposition, site), is_leaf(decomposition, site),
                     decomposition.is_root(parent) ? "root" : labeler.bounds_label(parent), {}};
    for (const net_handle_t& chain : child_chains(decomposition, site)) {
        std::vector<std::string> sites;
        for (const net_handle_t& child : child_sites(decomposition, chain)) {
            sites.push_back(labeler.bounds_label(child));
            collect_site_facts(decomposition, labeler, child, found);
        }
        facts.child_chains.emplace_back(labeler.bounds_label(chain), sites);
    }
    std::ranges::sort(facts.child_chains);
    found[labeler.bounds_label(site)] = facts;
}

/// @brief Find the facts about every site of a decomposition.
static SiteFactsByLabel all_site_facts(const SnarlDecomposition& decomposition, const HandleGraph& graph) {
    TraversalLabeler labeler(decomposition, graph);
    SiteFactsByLabel found;
    for (const net_handle_t& site : top_level_sites(decomposition)) {
        collect_site_facts(decomposition, labeler, site, found);
    }
    return found;
}

/**
 * @brief Require that the decomposition_sites functions say the same things
 * on the bench's two decompositions.
 *
 * @post Both decompositions' sites fit together, as require_site_tree()
 * checks.
 * @post Both decompositions have the same sites, matched by bounds label,
 * and each has the same depth, the same parent, whether it is a leaf, and
 * the same child chains holding the same sites in the same order.
 */
static void require_same_sites(const DecompositionBench& bench) {
    require_site_tree(bench.manager_view(), bench.graph());
    require_site_tree(bench.distance_index(), bench.graph());
    REQUIRE(all_site_facts(bench.manager_view(), bench.graph()) ==
            all_site_facts(bench.distance_index(), bench.graph()));
}

TEST_CASE("decomposition_sites gives the same sites on a SnarlManagerDecomposition and a distance index",
          "[decomposition_sites]") {

    SECTION("Three nested snarls give the same sites") {
        VG graph;
        Node* n1 = graph.create_node("GCA");
        Node* n2 = graph.create_node("T");
        Node* n3 = graph.create_node("G");
        Node* n4 = graph.create_node("CTGA");
        Node* n5 = graph.create_node("GCA");
        Node* n6 = graph.create_node("T");
        Node* n7 = graph.create_node("G");
        Node* n8 = graph.create_node("CTGA");
        graph.create_edge(n1, n2);
        graph.create_edge(n1, n8);
        graph.create_edge(n2, n3);
        graph.create_edge(n2, n6);
        graph.create_edge(n3, n4);
        graph.create_edge(n3, n5);
        graph.create_edge(n4, n5);
        graph.create_edge(n5, n7);
        graph.create_edge(n6, n7);
        graph.create_edge(n7, n8);
        IntegratedSnarlFinder finder(graph);
        DecompositionBench bench(graph, finder);
        require_same_sites(bench);

        const SnarlManagerDecomposition& view = bench.manager_view();
        TraversalLabeler labeler(view, graph);
        vector<net_handle_t> top_level = top_level_sites(view);
        REQUIRE(top_level.size() == 1);
        REQUIRE(labeler.bounds_label(top_level[0]) == "1+/8+");
        vector<net_handle_t> chains = child_chains(view, top_level[0]);
        REQUIRE(chains.size() == 1);
        vector<net_handle_t> middle = child_sites(view, chains[0]);
        REQUIRE(middle.size() == 1);
        REQUIRE(labeler.bounds_label(middle[0]) == "2+/7+");
        REQUIRE(depth(view, middle[0]) == 1);
        REQUIRE(!is_leaf(view, middle[0]));
        vector<net_handle_t> inner = child_sites(view, child_chains(view, middle[0])[0]);
        REQUIRE(inner.size() == 1);
        REQUIRE(labeler.bounds_label(inner[0]) == "3+/5+");
        REQUIRE(depth(view, inner[0]) == 2);
        REQUIRE(is_leaf(view, inner[0]));
    }

    SECTION("Random graphs give the same sites") {
        std::default_random_engine generator(test_seed_source());
        for (double chain_flip_probability : {0.0, 0.5}) {
            for (size_t repeat = 0; repeat < 100; repeat++) {
                bdsg::HashGraph graph = random_flipped_graph(generator);
                DecompositionBench bench(graph, fuzzed_decomposition(graph, chain_flip_probability, generator));
                require_same_sites(bench);
            }
        }
    }
}

TEST_CASE("site_bounds reads in where a traversal starts and out where it ends", "[decomposition_sites]") {
    VG graph;
    Node* n1 = graph.create_node("GCA");
    Node* n2 = graph.create_node("T");
    Node* n3 = graph.create_node("G");
    Node* n4 = graph.create_node("CTGA");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n4);
    // Turn around on the end of node 2, so the snarl from node 1 to node 4
    // can be entered and left at its start.
    graph.create_edge(n2, n2, false, true);
    IntegratedSnarlFinder finder(graph);
    DecompositionBench bench(graph, finder);
    const SnarlManagerDecomposition& view = bench.manager_view();
    vector<net_handle_t> sites = top_level_sites(view);
    REQUIRE(sites.size() == 1);
    handle_t start = graph.get_handle(n1->id());
    handle_t end = graph.get_handle(n4->id());
    REQUIRE(site_bounds(view, graph, sites[0]) == SiteBounds {start, end});
    REQUIRE(site_bounds(view, graph, view.flip(sites[0])) == SiteBounds {graph.flip(end), graph.flip(start)});
    net_handle_t start_to_start = view.get_parent_traversal(view.get_bound(sites[0], false, true),
                                                            view.get_bound(sites[0], false, false));
    REQUIRE(view.starts_at(start_to_start) == SnarlDecomposition::START);
    REQUIRE(view.ends_at(start_to_start) == SnarlDecomposition::START);
    REQUIRE(site_bounds(view, graph, start_to_start) == SiteBounds {start, graph.flip(start)});
    require_flipped_bounds(view, graph, sites[0]);
}

TEST_CASE("decomposition_sites follows the manager's snarl tree", "[decomposition_sites]") {
    VG graph;
    Node* n1 = graph.create_node("GCA");
    Node* n2 = graph.create_node("T");
    Node* n3 = graph.create_node("G");
    Node* n4 = graph.create_node("CTGA");
    Node* n5 = graph.create_node("GCA");
    Node* n6 = graph.create_node("T");
    graph.create_edge(n1, n2);
    graph.create_edge(n1, n6);
    graph.create_edge(n2, n3);
    graph.create_edge(n2, n4);
    graph.create_edge(n3, n5);
    graph.create_edge(n4, n5);
    graph.create_edge(n5, n6);
    IntegratedSnarlFinder finder(graph);
    SnarlManager manager = finder.find_snarls();
    SnarlManagerDecomposition view(manager, graph);
    manager.for_each_snarl_preorder([&manager, &view, &graph](const Snarl* snarl) {
        if (manager.shallow_contents(snarl, graph, false).first.empty()) {
            // The view hides snarls with no children.
            return;
        }
        net_handle_t site = view.get_snarl_net(snarl);
        REQUIRE(site_bounds(view, graph, site) ==
                SiteBounds {graph.get_handle(snarl->start().node_id(), snarl->start().backward()),
                            graph.get_handle(snarl->end().node_id(), snarl->end().backward())});
        net_handle_t parent = parent_site(view, site);
        if (manager.parent_of(snarl) == nullptr) {
            REQUIRE(view.is_root(parent));
        } else {
            REQUIRE(view.get_snarl(parent) == manager.parent_of(snarl));
        }
        REQUIRE(is_leaf(view, site) == manager.all_children_trivial(snarl, graph));
    });
}

}
