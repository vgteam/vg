#include "0_compare_gbwt_to_graph.hpp"
namespace vg {
namespace algorithms{

void compare_gbwt_to_graph(gbwt::GBWT& gbwt, HandleGraph& graph)
{
    // cerr << "Does graph contain all the nodes in gbwt?" << endl;
    // gbwt::node_type node = gbwt.firstNode();

    cerr << "Does gbwt contain all the nodes in graph?" << endl;
    graph.for_each_handle([&](const handle_t handle)
    {
        // (graph.get_id(handle), graph.get_is_reverse(handle))
        gbwt::node_type node = handle_to_gbwt(graph, handle); 

        bool contains = gbwt.contains(node);
        if (!contains)
        {
            cerr << "gbwt does not contain node: " << graph.get_id(handle) << endl;
        }
        // else
        // {
        //     cerr << "gbwt contains node:" << graph.get_id(handle) << endl;
        // }
    });
    // for (handle : )
}



}
}