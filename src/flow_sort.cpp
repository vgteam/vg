#include "vg.hpp"
#include <vg/io/stream.hpp>
#include "aligner.hpp"
#include <vg/vg.pb.h>
#include "flow_sort.hpp"
#include <raptor2/raptor2.h>

namespace vg {

using namespace std;
FlowSort::FlowSort(VG& vg):vg(vg) {
    
}

std::unique_ptr< list<NodeTraversal> >  FlowSort::max_flow_sort(const string& ref_name, bool isGrooming) 
{
    std::unique_ptr <list<NodeTraversal> > sorted_nodes(new list<NodeTraversal>());
    if (vg.size() <= 1) return sorted_nodes;
    // Topologically sort, which orders and orients all the nodes.
    
    
    vg.paths.sort_by_mapping_rank();
    
    flow_sort_nodes(*sorted_nodes, ref_name, isGrooming);
    list<NodeTraversal>::reverse_iterator n = sorted_nodes->rbegin();
    int i = 0;
    for ( ; i < vg.graph.node_size() && n != sorted_nodes->rend();
          ++i, ++n) {
        // Put the nodes in the order we got
        vg.swap_nodes(vg.graph.mutable_node(i), (*n).node);
    }
    vg.rebuild_indexes();
    return sorted_nodes;
}

void FlowSort::fast_linear_sort(const string& ref_name, bool isGrooming)
{
    if (vg.size() <= 1) return;
    // Topologically sort, which orders and orients all the nodes.
    list<NodeTraversal> sorted_nodes;
    vg.paths.sort_by_mapping_rank();
    
    //get weighted graph
    WeightedGraph w_graph;
    w_graph.construct(*this,ref_name, isGrooming);

    //all nodes size
    id_t nodes_size =vg.node_by_id.size();

    //create set of sinks and sources
    std::set<id_t> sources;
    //index - degree of nodes
    std::vector<std::set<id_t>> nodes_degree;
    int cur_node_degree;
    //find sources
    for (auto const &entry : vg.node_by_id)
    {
        if(!w_graph.edges_in_nodes.count(entry.first) || 
                w_graph.edges_in_nodes[entry.first].size() == 0)
        {
            sources.insert(entry.first);
            continue;
        }
        cur_node_degree = get_node_degree(w_graph, entry.first);
        if (cur_node_degree < 0)
            continue;

        if (cur_node_degree + 1 > nodes_degree.size())
            nodes_degree.resize(cur_node_degree + 1);
        nodes_degree[cur_node_degree].insert(entry.first);
        //if(edges_out_nodes.find(entry.first) != edges_out_nodes.end() )
        //   sinks.insert(entry.first);
    }

    id_t next = -1;
    //std::vector<id_t> sorted;


    int remaining_nodes_size = nodes_size;
    set<id_t>::iterator i_last_el;

    while (remaining_nodes_size > 0)
    {   //Get next vertex to delete
        if (next == -1)
        {
            if (sources.size() > 0) //|| sinks.size() > 0)
            {
                i_last_el = sources.end();
                --i_last_el;
                next = *i_last_el;
                sources.erase(i_last_el);
            }
            else
                next = find_max_node(nodes_degree);
        }
        else
        {
            i_last_el = sources.find(next);
            if(i_last_el != sources.end())
                sources.erase(i_last_el);
        }

        //add next node
        //sorted.push_back(next);
        NodeTraversal node = NodeTraversal(vg.node_by_id[next], false);
        sorted_nodes.push_back(node);

        //remove edges related with node
        next = get_next_node_recalc_degrees(w_graph, nodes_degree,sources,
                                            next);
        remaining_nodes_size--;
    }

    //output
    list<NodeTraversal>::iterator n = sorted_nodes.begin();
    int i = 0;
    for ( ; i < vg.graph.node_size() && n != sorted_nodes.end();
          ++i, ++n) {
        // Put the nodes in the order we got
        vg.swap_nodes(vg.graph.mutable_node(i), (*n).node);
    }
    vg.rebuild_indexes();
}

void FlowSort::flow_sort_nodes(list<NodeTraversal>& sorted_nodes, 
        const string& ref_name, bool isGrooming) 
{

    //assign weight to edges
    //edge weight is determined as number of paths, that go through the edge
    WeightedGraph w_graph;
    w_graph.construct(*this, ref_name, isGrooming);
    list<mapping_t> ref_path(vg.paths.get_path(ref_name).begin(),
                             vg.paths.get_path(ref_name).end());
    ref_path.reverse();
    Growth growth;

    for(auto const &mapping : ref_path) 
    {
        growth.backbone.insert(mapping.node_id());
        growth.ref_path.push_back(mapping.node_id());
    }
    for (auto const &entry : vg.node_by_id) 
    {
        growth.nodes.insert(entry.first);
    }

    set<id_t> unsorted_nodes(growth.nodes.begin(), growth.nodes.end());
    
    find_in_out_web(sorted_nodes, growth, w_graph, unsorted_nodes, -1, 
            false, 0);
   
    if (sorted_nodes.size() > vg.graph.node_size()) 
    {
#ifdef debug
        cerr << "Failed to sort graph " << endl;
        cerr << "sorted: " << sorted_nodes.size() << " in graph: " << endl;
#endif
        return;
    }
    while (sorted_nodes.size() < vg.graph.node_size()) 
    {
#ifdef debug
        cerr << "additional sorting for missing nodes" << endl;
        cerr << "unsorted: " << unsorted_nodes.size() << endl;
        cerr << "sorted: " << sorted_nodes.size() << " in graph: " <<
                graph.node_size() << endl;
        cerr << "unsorted: ";
        
        for (auto const &uns : unsorted_nodes) 
        {
            cerr << uns << " ";
        }
        
        cerr << endl;
#endif
        
        list<NodeTraversal> sorted_nodes_new;
    
        set<id_t> unsorted_nodes_new(growth.nodes.begin(), growth.nodes.end());
        Growth growth_new;
        growth_new.nodes.insert(growth.nodes.begin(), growth.nodes.end());
        for (auto const &node : sorted_nodes) 
        {
            growth_new.backbone.insert(node.node->id());
            growth_new.ref_path.push_back(node.node->id());
        }
      
        WeightedGraph wgraph_new;
        wgraph_new.construct(*this, ref_name);
        
        find_in_out_web(sorted_nodes_new, growth_new, wgraph_new, 
                unsorted_nodes_new, -1, false, 0);

        sorted_nodes = sorted_nodes_new;
        if (sorted_nodes.size() != vg.graph.node_size() && 
                unsorted_nodes.size() == unsorted_nodes_new.size())
        {
#ifdef debug          
            cerr << "Failed to insert missing nodes"<< endl;
#endif
            break;
        }
        
        unsorted_nodes = unsorted_nodes_new;
    }
    
    if (unsorted_nodes.size() != 0) 
    {
        for (auto const &id : unsorted_nodes) 
        {
            NodeTraversal node = NodeTraversal(vg.node_by_id[id], false);
            sorted_nodes.push_back(node);
        }
    }
}

/* return degree of the node in the WeightedGraph*/
int FlowSort::get_node_degree(FlowSort::WeightedGraph& wg, id_t node_id)
{

    std::vector<Edge*> in_edges = wg.edges_in_nodes[node_id];
    int degree = 0;
    for(auto const &edge : in_edges)
        degree -= wg.edge_weight[edge];

    std::vector<Edge*> out_edges = wg.edges_out_nodes[node_id];
    for(auto const &edge : out_edges)
        degree += wg.edge_weight[edge];

    return degree;
}

/* Weight is assigned to edges as number of paths, that go through hat edge.
   Path goes through the edge, if both adjacent nodes of the edge are mapped to that path.*/
void 
FlowSort::WeightedGraph::construct(FlowSort& fs, const string& ref_name, bool isGrooming)
{
    int ref_weight = fs.vg.paths._paths.size();
    if (ref_weight < FlowSort::DEFAULT_PATH_WEIGHT) {
        ref_weight = FlowSort::DEFAULT_PATH_WEIGHT;
    }

    fs.vg.flip_doubly_reversed_edges();
    //"bad" edges
    map<id_t, set<Edge*>> minus_start;//vertex - edges
    map<id_t, set<Edge*>> minus_end;//vertex - edges
    set<id_t> nodes;
    id_t start_ref_node = 0;
    for (auto &edge : fs.vg.edge_index)
    {
        id_t from = edge.first->from();
        id_t to = edge.first->to();
        if (edge.first->from_start() || edge.first->to_end())
        {
                minus_start[from].insert(edge.first);
                minus_end[to].insert(edge.first);
        }
        else
        {
            edges_out_nodes[edge.first->from()].push_back(edge.first);
            edges_in_nodes[edge.first->to()].push_back(edge.first);
        }
        //Node* nodeTo = node_by_id[edge.first->to()];
        nodes.insert(edge.first->from());
        nodes.insert(edge.first->to());

        //assign weight to the minimum number of paths of the adjacent nodes
        auto from_node_mapping = fs.vg.paths.get_node_mapping_by_path_name(from);
//        NodeMapping to_node_mapping = paths.get_node_mapping(to);
        int weight = 1;

        for (auto const &path_mapping : from_node_mapping) {
            string path_name = path_mapping.first;
            if (fs.vg.paths.are_consecutive_nodes_in_path(from, to, path_name))
            {
                if (path_name == ref_name)
                {
                    weight += ref_weight;
                    if(start_ref_node == 0)
                        start_ref_node = edge.first->from();
                }
                weight++;
            }
        }

        edge_weight[edge.first] = weight;
#ifdef debug        
        cerr << from << "->" << to << " " << weight << endl;
#endif
    }
    if (isGrooming)
    {
        // get connected components
        vector<set<id_t>> ccs = fs.get_cc_in_wg(edges_in_nodes, edges_out_nodes, nodes, start_ref_node);
        // grooming
        id_t main_cc = 0;
        if(ccs.size() > 1)
        {
            for(id_t j = 0; j < ccs.size(); j++)
            {
                if(j != main_cc)
                    fs.groom_components(edges_in_nodes, edges_out_nodes, ccs[j], ccs[main_cc], minus_start, minus_end);

            }
        }
        fs.vg.rebuild_edge_indexes();
    }
//    return WeightedGraph(edges_out_nodes, edges_in_nodes, edge_weight);
}


void FlowSort::groom_components(EdgeMapping& edges_in, EdgeMapping& edges_out, set<id_t>& isolated_nodes, set<id_t>& main_nodes,
                          map<id_t, set<Edge*>> &minus_start, map<id_t, set<Edge*>> &minus_end)
{
    vector<Edge*> from_minus_edges;
    vector<Edge*> to_minus_edges;
    auto nodes_it = isolated_nodes.begin();
    //find all "bad" edeges
    while (nodes_it != isolated_nodes.end())
    {
        if(minus_start.find(*nodes_it) != minus_start.end())
        {
            for(auto& e: minus_start[*nodes_it])
            {
                //find edge between main component and isolated component
                if(main_nodes.find(e->from()) != main_nodes.end() || main_nodes.find(e->to()) != main_nodes.end() )
                {
                    if(e->from_start())
                        from_minus_edges.push_back(e);
                    if(e->to_end())
                        to_minus_edges.push_back(e);
                    //break;
                }
            }
        }
        id_t lol;
        if(minus_end.find(*nodes_it) != minus_end.end())
        {
            lol = *nodes_it;
            for(auto& e: minus_end[*nodes_it])
            {
                //find edge between main component and isolated component
                if(main_nodes.find(e->from()) != main_nodes.end() || main_nodes.find(e->to()) != main_nodes.end() )
                {
                    if(e->from_start())
                        from_minus_edges.push_back(e);
                    if(e->to_end())
                        to_minus_edges.push_back(e);
                    //break;
                }
            }
        }
        nodes_it++;
    }
    //reverse "bad" edges
    for(auto& e: from_minus_edges)
    {
        if(main_nodes.find(e->from()) != main_nodes.end())
        {
            from_simple_reverse_orientation(e);
            update_in_out_edges(edges_in,edges_out, e);
            continue;
        }
        if(main_nodes.find(e->to()) != main_nodes.end())
        {
            from_simple_reverse(e);
            update_in_out_edges(edges_in,edges_out, e);
        }
    }
    for(auto& e: to_minus_edges)
    {
        if(main_nodes.find(e->to()) != main_nodes.end())
        {
            to_simple_reverse_orientation(e);
            update_in_out_edges(edges_in,edges_out, e);
            continue;
        }
        if(main_nodes.find(e->from()) != main_nodes.end())
        {
            to_simple_reverse(e);
            update_in_out_edges(edges_in,edges_out, e);
        }
    }

    nodes_it = isolated_nodes.begin();
    Edge* internal_edge;
    vector<Edge*> edges_to_flip;
   
    //reverse all internal edges
    while (nodes_it != isolated_nodes.end())
    {    
        //reverse-complement node sequence
        Node* node = vg.get_node(*nodes_it);
        string* sequence = node->mutable_sequence();
        if (sequence) {
            reverse(sequence->begin(), sequence->end());
            for (string::size_type i = 0; i < sequence->size(); ++i) {
                char& nucleotide = sequence->at(i);
                if (COMPLEMENTARY_NUCLEOTIDES.count(nucleotide)) {
                    nucleotide = COMPLEMENTARY_NUCLEOTIDES.at(nucleotide);
                }
            }
        }
        
        for(auto& e:edges_in[*nodes_it])
        {                
            // if edge is internal
            if(isolated_nodes.find(e->from()) != isolated_nodes.end())
            {
                internal_edge = e;
                edges_to_flip.push_back(e);
            }
        }
        nodes_it++;
    }
    for(auto& e: edges_to_flip)
    {
        internal_edge = e;
        erase_in_out_edges(edges_in, edges_out, internal_edge);
        reverse_edge(internal_edge);
        update_in_out_edges(edges_in,edges_out, internal_edge);
    }
    //insert isolated set to main set
    main_nodes.insert(isolated_nodes.begin(), isolated_nodes.end());
}

void FlowSort::update_in_out_edges(EdgeMapping& edges_in, EdgeMapping& edges_out, Edge* e)
{
    edges_in[e->to()].push_back(e);
    edges_out[e->from()].push_back(e);
}

void FlowSort::erase_in_out_edges(EdgeMapping& edges_in, EdgeMapping& edges_out, Edge* e)
{
    int i = 0;
    while(*(edges_in[e->to()].begin() + i) != e)
        i++;
    edges_in[e->to()].erase(edges_in[e->to()].begin() + i);
    i = 0;
    while(*(edges_out[e->from()].begin() + i) != e)
        i++;
    edges_out[e->from()].erase(edges_out[e->from()].begin() + i);
}

void FlowSort::reverse_from_start_to_end_edge(Edge* &e)
{
    e->set_from_start(false);
    e->set_to_end(false);
    reverse_edge(e);
}


void FlowSort::reverse_edge(Edge* &e)
{
    id_t tmp_vrtx = e->to();
    e->set_to(e->from());
    e->set_from(tmp_vrtx);
}


// a(from_start ==true) -> b        =>        not a (from_start == false)  -> b
id_t FlowSort::from_simple_reverse(Edge* &e)
{
    e->set_from_start(false);
    return e->from();
}

// b(from_start ==true) -> a        =>        not a (from_start == false)  -> b
id_t FlowSort::from_simple_reverse_orientation(Edge* &e)
{
    e->set_from_start(false);
    reverse_edge(e);
    return e->from();
}

// a -> b (to_end ==true)       =>        a -> not b(to_end ==false)
id_t FlowSort::to_simple_reverse(Edge* &e)
{
    e->set_to_end(false);
    return e->to();
}

// b -> a (to_end ==true)       =>        a -> not b(to_end ==false)
id_t FlowSort::to_simple_reverse_orientation(Edge* &e)
{
    reverse_edge(e);
    e->set_to_end(false);
    return e->to();
}

vector<set<id_t>> FlowSort::get_cc_in_wg(EdgeMapping& edges_in,EdgeMapping& edges_out,
                                   const set<id_t>& all_nodes, id_t start_ref_node)
{
    set<id_t> nodes(all_nodes.begin(), all_nodes.end());
    vector<set<id_t>> result;
    bool main_cc = true;
    id_t s = start_ref_node;
    nodes.erase(nodes.find(s));
    while(nodes.size() > 0)
    {
        set<id_t> visited;
        if(!main_cc)
        {
            s = *nodes.begin();
            nodes.erase(nodes.begin());
        }
        main_cc = false;
        std::stack<id_t> q({ s });
        while (!q.empty())
        {
            s = q.top();
            q.pop();
            if (visited.find(s) != visited.end())
                continue;
            nodes.erase(s);
            visited.insert(s);
            for(const auto& e: edges_in[s])
            {
                id_t from = e->from();
                if (visited.find(from ) == visited.end())
                    q.push(from);
            }
            for(const auto& e: edges_out[s])
            {
                id_t to = e->to();
                if (visited.find(to) == visited.end())
                    q.push(to);
            }
        }
        result.push_back(visited);
    }
    return result;
}

/* Iterate all edges adjacent to node, recalc degrees of related nodes.
 * If node has no incoming edges we add it to the sources and return as next node.
 * */
id_t FlowSort::get_next_node_recalc_degrees(WeightedGraph& wg, std::vector<std::set<id_t>>& degrees,std::set<id_t> &sources,
                                     id_t node)
{
    id_t result = -1;
    //recalc degrees
    id_t cur_node;
    int new_degree;
    int cur_degree;
    //remove from degrees (if node is not source)
    int node_degree = get_node_degree(wg, node);
    if(node_degree >= 0 && (node_degree < degrees.size()))
    {
        set<id_t>::iterator i_el = degrees[node_degree].find(node);
        if(i_el != degrees[node_degree].end())
            degrees[node_degree].erase(i_el);
    }


    //decrease number of in_edges from current node
    for(const auto& edge: wg.edges_in_nodes[node])
    {
        cur_node = edge->from();
        cur_degree = get_node_degree(wg, cur_node);
        new_degree = cur_degree - wg.edge_weight[edge];

        //decrease out_edges from current node
        if(cur_degree >= 0)
            degrees[cur_degree].erase(cur_node);
        //add node to new_degree set
        if(new_degree >= 0)
        {
            if (new_degree + 1> degrees.size())
                degrees.resize(new_degree + 1);
            degrees[new_degree].insert(cur_node);
        }


        //remove current edge from wg
        std::vector<Edge*>& related_edges = wg.edges_out_nodes[cur_node];
        related_edges.erase(std::remove(related_edges.begin(), related_edges.end(), edge), related_edges.end());
    }
    for(const auto& edge: wg.edges_out_nodes[node])
    {
        cur_node = edge->to();
        //check if node is new source
        std::vector<Edge*>& related_edges = wg.edges_in_nodes[cur_node];
        //add related node to sources
        if(related_edges.size() == 1)
        {
            sources.insert(cur_node);
            result = cur_node;
        }

        cur_degree = get_node_degree(wg, cur_node);
        new_degree = cur_degree + wg.edge_weight[edge];

        //decrease in_edges from current node
        if(cur_degree >= 0)
            degrees[cur_degree].erase(cur_node);
        //add node to new_degree set (if not source)
        if(new_degree >= 0 && related_edges.size() != 1)
        {
            if (new_degree + 1> degrees.size())
                degrees.resize(new_degree + 1);
            degrees[new_degree].insert(cur_node);
        }
        //remove current edge from wg
        related_edges.erase(std::remove(related_edges.begin(), related_edges.end(), edge), related_edges.end());
    }

    return result;
}

/*Find node with max degree*/
id_t FlowSort::find_max_node(std::vector<std::set<id_t>> nodes_degree)
{
    int start_size = nodes_degree.size()-1;
    int last_ind = start_size;
    while(nodes_degree[last_ind].size() == 0)
    {
        last_ind--;
    }
    if(last_ind != start_size)
        nodes_degree.resize(last_ind + 1);
    set<id_t>::iterator i_last_el = nodes_degree[last_ind].end();
    return *(--i_last_el);
}

/*  Method finds min cut in a given set of nodes, then removes min cut edges,
    finds in- and -out growth from the reference path and calls itself on them recursively.
*/
void FlowSort::find_in_out_web(list<NodeTraversal>& sorted_nodes,
                            Growth& io_growth,
                            WeightedGraph& w_graph,
                            set<id_t>& unsorted_nodes, 
                            id_t start, bool in_out, int count) 
{
#ifdef debug
    cerr << "enter recursion: " << count << endl;
#endif
    
    set<id_t>& backbone = io_growth.backbone;
    set<id_t>& nodes = io_growth.nodes;
    list<id_t>& ref_path = io_growth.ref_path;

    //for efficiency if size of the backbone == size of the nodes
    //we just add all backbone to sorted nodes and quit
    if (nodes.size() == backbone.size()) 
    {
        if (in_out) {
            ref_path.reverse();
        }
        for (auto const &id: ref_path) 
        {
            if (!unsorted_nodes.count(id) || id == start) 
            {
                continue;
            }
            NodeTraversal node (vg.node_by_id[id], false);
            sorted_nodes.push_back(node);
            unsorted_nodes.erase(id);
//            cerr << "erasing " << id << endl;
        }
#ifdef debug        
        cerr << "backbone simple from " << start << ": ";
        for (auto const &id : ref_path) 
        {
            cerr << id << " ";
        }
        cerr << endl;
        cerr << "leaving recursion: " << count << endl;
#endif
        return;
    }

    EdgeMapping& edges_out_nodes = w_graph.edges_out_nodes;
    EdgeMapping& edges_in_nodes = w_graph.edges_in_nodes;
    map<Edge*, int>& edge_weight = w_graph.edge_weight;

    //add source and sink nodes
    id_t source = vg.graph.node_size() + 1;
    id_t sink = vg.graph.node_size() + 2;
    id_t graph_size = sink + 1;

    set<Edge*> out_joins;
    set<Edge*> in_joins;

    map<id_t, map<id_t, int>> weighted_sub_graph;

    //determine max outgoing weight
    int max_weight = 0;
    for(auto const &node : nodes) 
    {
        int current_weight = 0;
        for (auto const &edge : edges_out_nodes[node]) 
        {
            weighted_sub_graph[edge->from()][edge->to()] = 0;
            if (!nodes.count(edge->to())) 
            {
                continue;
            }
            if (backbone.count(node) && backbone.count(edge->to())) 
            {             
                continue;
            }
            current_weight += edge_weight[edge];
            if (backbone.count(node) && !backbone.count(edge->to())) 
            {
                out_joins.insert(edge);
            }

            if (!backbone.count(node) && backbone.count(edge->to())) 
            {
                in_joins.insert(edge);
                continue;
            }
            weighted_sub_graph[edge->from()][edge->to()] = edge_weight[edge];
        }

        if (current_weight > max_weight) 
        {
            max_weight = current_weight;
        }
    }
    //set weight for source edges
    for (auto const &edge : out_joins) 
    {
        weighted_sub_graph[source][edge->from()] = max_weight + 1;
    }
    //set weight for sink edges
    for (auto const &edge : in_joins) 
    {
        weighted_sub_graph[edge->from()][sink] = edge_weight[edge];
    }

    //find min-cut edges
    vector<pair<id_t,id_t>> cut = min_cut(weighted_sub_graph, nodes, source,
                                                    sink, edges_out_nodes, in_joins);
    //remove min cut edges
    for(auto const &edge : cut) 
    {
        id_t to = edge.second;
        if (to == sink) 
        {
            for (auto const &in_join : in_joins) 
            {
                if (in_join->from() == edge.first) 
                {
                    to = in_join->to();
                }
            }
        }
//        cerr << " removing min-cut edge " <<  edge.first << "->" << to << endl;  
        remove_edge(edges_out_nodes, edge.first, to, false);
        remove_edge(edges_in_nodes, to, edge.first, true);
    }

    set<id_t> visited;
    
    visited.insert(backbone.begin(), backbone.end());

    for (auto const &current_id: ref_path) 
    {
        list<NodeTraversal> sort_out;
        NodeTraversal node = NodeTraversal(vg.node_by_id[current_id], false);
        //out growth
        process_in_out_growth(edges_out_nodes, current_id, io_growth,
                    w_graph, visited, sort_out, false, unsorted_nodes, 
                in_out, count);
        
        if (sort_out.size() != 0) 
        {
#ifdef debug            
            cerr << "out growth " << current_id << ": ";
            for (auto const &n : sort_out) 
            {
                cerr << n.node->id() << " ";
            }
            cerr << endl;
#endif
            sorted_nodes.insert(sorted_nodes.end(), sort_out.begin(), sort_out.end());
        }
        //add backbone node to the result
        if (unsorted_nodes.count(current_id) || current_id == start) 
        {
            sorted_nodes.push_back(node);
            if (current_id != start) 
            {
                unsorted_nodes.erase(current_id);
            }
        }
    }
#ifdef debug    
    cerr << "backbone backward from " << start << ": ";
    for (auto const &id : ref_path) 
    {
        cerr << id << " ";
    }
    cerr << endl;
    
    cerr << "after backward from " << start << ": ";
    for (auto const &id : sorted_nodes)
    {
        cerr << id.node->id() << " ";
    }
    cerr << endl;
#endif    
    
    list<NodeTraversal> new_sorted;
    sorted_nodes.reverse();
    for (auto const &node : sorted_nodes) 
    {
        id_t current_id = node.node->id();     
        list<NodeTraversal> sort_in;
        //in growth
        process_in_out_growth(edges_in_nodes, current_id, io_growth,
                    w_graph, visited, sort_in, true, unsorted_nodes, 
                    in_out, count);
        if (sort_in.size() != 0) 
        {
#ifdef debug            
            cerr << "in growth "<< current_id << ": ";
            for (auto const &n : sort_in) 
            {
                cerr << n.node->id() << " ";
            }
            cerr << endl;
#endif            
            new_sorted.insert(new_sorted.end(), sort_in.begin(), sort_in.end());
        }
        if (current_id != start) 
        {
            new_sorted.push_back(node);
            unsorted_nodes.erase(current_id);
        }
    }
    if (!in_out) 
    {
        new_sorted.reverse();
    }
#ifdef debug    
    cerr << "backbone forward from " << start << ": ";
    for (auto const &id : sorted_nodes) 
    {
        cerr << id.node->id() << " ";
    }
    cerr << endl;
#endif
    sorted_nodes = new_sorted;
#ifdef debug    
    cerr << "growth sorted from " << start << ": ";
    for (auto const &id : sorted_nodes) 
    {
        cerr << id.node->id() << " ";
    }
 
    cerr << endl;
    cerr << "leaving recursion: " << count << endl;
#endif
}
/*
        Determines the presence of a in- out- growth, finds its backbone and calls min cut algorithm.
     */
void FlowSort::process_in_out_growth(EdgeMapping& nodes_to_edges, id_t current_id,
        Growth& io_growth,
        WeightedGraph& w_graph,
        set<id_t>& visited,
        list<NodeTraversal>& sorted_nodes,
        bool reverse,
        set<id_t>& unsorted_nodes,
        bool in_out, int count) 
{

    if (!nodes_to_edges.count(current_id) || nodes_to_edges[current_id].size() == 0)
        return;
    
    set<id_t>& backbone = io_growth.backbone;
    set<id_t>& nodes = io_growth.nodes;
    
    Growth growth_new;
    id_t start_node = current_id;
    
    while (true) 
    {
        if (growth_new.backbone.count(start_node) || (start_node != current_id && 
                visited.count(start_node))) 
        {
            break;
        }
        growth_new.backbone.insert(start_node);
        growth_new.ref_path.push_back(start_node);
        int weight = 0;
        Edge* next_edge;
        //take edges with maximum weight to the new reference path
        for (auto const &out_edge : nodes_to_edges[start_node]) 
        {
            if (w_graph.edge_weight[out_edge] > weight) 
            {
                id_t tmp = reverse ? out_edge->from() : out_edge->to();
                if (!nodes.count(tmp) || backbone.count(tmp) ||
                        growth_new.backbone.count(tmp) || visited.count(tmp)) 
                {
                    continue;
                }
                start_node = tmp;
                weight = w_graph.edge_weight[out_edge];
                next_edge = out_edge;

            }
        }
    }

    mark_dfs(nodes_to_edges, current_id, growth_new.nodes, visited, reverse, nodes, backbone);
    if (growth_new.nodes.size() == 1 && growth_new.nodes.count(current_id)) 
    {
        return;
    }
    if (!reverse)
    {
        growth_new.ref_path.reverse();
    }
    find_in_out_web(sorted_nodes, growth_new, w_graph, unsorted_nodes, 
            current_id, reverse , count+1);
   
}

void FlowSort::remove_edge(EdgeMapping& nodes_to_edges, id_t node, id_t to, bool reverse)
{
    if (nodes_to_edges.count(node))
    {
        auto& edges = nodes_to_edges[node];
        auto it = begin(edges);
        while (it != end(edges)) 
        {
            id_t next = reverse ? (*it)->from() : (*it)->to();
            if (next == to) 
            {
                edges.erase(it);
            } else 
            {
                ++it;
            }
        }
    }
}

/*
    Adds the nodes in dfs discovery order to sorted_nodes
 */
void FlowSort::mark_dfs(EdgeMapping& edges_to_nodes, id_t s,
            set<id_t>& new_nodes, set<id_t>& visited, bool reverse,
        set<id_t>& nodes,
        set<id_t>& backbone) 
{
    visited.insert(s);
    new_nodes.insert(s);
    std::stack<id_t> q({ s });
    while(!q.empty()) 
    {
        s = q.top();
        q.pop();
        if (edges_to_nodes.count(s)) 
        {
            for (auto const &edge: edges_to_nodes[s]) 
            {
                id_t next_node = reverse ? edge->from() : edge->to();
                if (!visited.count(next_node) && nodes.count(next_node)
                        && !backbone.count(next_node)) 
                {
                    visited.insert(next_node);
                    new_nodes.insert(next_node);
                    q.push(next_node);
                }
            }
        }
    }

}

/* Returns true if there is a path from source 's' to sink 't' in
   residual graph. Also fills parent[] to store the path */
bool FlowSort::bfs(set<id_t>& nodes, map<id_t, map<id_t, int>>& edge_weight, id_t s,
            id_t t, map<id_t, id_t>& parent) 
{
    // Create a visited array and mark all vertices as not visited
    set<id_t> visited;

    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    std::queue<id_t> q({ s });
    visited.insert(s);
    parent[s] = -1;

    // Standard BFS Loop
    while (!q.empty())
    {
        id_t u = q.front();
        q.pop();
        for (auto const &weight : edge_weight[u])
        {
            if (weight.second <= 0)
            {
                continue;
            }
            id_t v = weight.first;
            if (!visited.count(v)) 
            {
                q.push(v);
                parent[v] = u;
                visited.insert(v);
                  if (v == t) 
                {
                    return true;
                }
            }
        }
    }
    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited.count(t) != 0);

}


void FlowSort::dfs(set<id_t>& nodes, id_t s, set<id_t>& visited, 
        map<id_t, map<id_t, int>>& edge_weight) 
{
    visited.insert(s);
    std::stack<id_t> q({ s });
    while (!q.empty()) {
        s = q.top();
        q.pop();
        for (auto const &weight : edge_weight[s]) 
        {
            id_t to = weight.first;
            if (weight.second && !visited.count(to)) 
            {
                visited.insert(to);
                q.push(to);
            }
        }
    }
}

// Returns the minimum s-t cut
vector<pair<id_t,id_t>> FlowSort::min_cut(map<id_t, map<id_t, int>>& graph_weight,
                set<id_t>& nodes,
                id_t s, id_t t,
                EdgeMapping& edges_out_nodes,
                set<Edge*>& in_joins) 
{
    id_t u, v;
    nodes.insert(s);
    nodes.insert(t);
    map<id_t, id_t> parent_map;

    // Augment the flow while there is path from source to sink
    while (bfs(nodes, graph_weight, s, t, parent_map)) 
    {
        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        int path_flow = numeric_limits<int>::max();
        for (v = t; v != s; v = parent_map[v]) 
        {
            u = parent_map[v];
            if (!graph_weight.count(u) || !graph_weight[u].count(v)) 
            {
                graph_weight[u][v] = 0;
            }
            path_flow = min(path_flow, graph_weight[u][v]);
        }
        // update residual capacities of the edges and reverse edges
        // along the path
        for (v = t; v != s; v = parent_map[v]) 
        {
            u = parent_map[v];
            if (!graph_weight.count(u) || !graph_weight[u].count(v)) 
            {
                graph_weight[u][v] = 0;
            }
             if (!graph_weight.count(v) || !graph_weight[v].count(u)) 
            {
                graph_weight[v][u] = 0;
            }
            graph_weight[u][v] -= path_flow;
            graph_weight[v][u] += path_flow;
        }
    }

    // Flow is maximum now, find vertices reachable from s

    set<id_t> visited_set;
    dfs(nodes, s, visited_set, graph_weight);
    vector<pair<id_t,id_t>> min_cut;

    nodes.erase(s);
    nodes.erase(t);

    for(auto const node_id : nodes) 
    {
        for (auto const &edge : edges_out_nodes[node_id])
        {
            id_t to = edge->to();
            if (!nodes.count(to)) 
            {
                continue;
            }

            if (in_joins.count(edge)) 
            {
                 to = t;
            }
            if (visited_set.count(node_id) && !visited_set.count(to))
            {
                min_cut.push_back(pair<id_t, id_t>(node_id, to));
            }
        }
    }
    return min_cut;
}

} // end namespace
