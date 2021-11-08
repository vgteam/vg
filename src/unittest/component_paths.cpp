#include "catch.hpp"
#include "../utility.hpp"
#include "../algorithms/component_paths.hpp"

#include <bdsg/hash_graph.hpp>

namespace vg {
namespace unittest {

using namespace std;

set<set<path_handle_t>> normalize(vector<unordered_set<path_handle_t>> result) {
    set<set<path_handle_t>> return_val;
    for (auto& comp_set : result) {
        return_val.emplace(comp_set.begin(), comp_set.end());
    }
    return return_val;
}

TEST_CASE("Parallel component paths produces correct results", "[comppathset]") {
    
    
    int thread_count_pre = get_thread_count();
    
    SECTION("On a graph with several components") {
        
        bdsg::HashGraph graph;
        
        handle_t prev;
        path_handle_t p1 = graph.create_path_handle("1");
        path_handle_t p2 = graph.create_path_handle("2");
        path_handle_t p3 = graph.create_path_handle("3");
        path_handle_t p4 = graph.create_path_handle("4");
        path_handle_t p5 = graph.create_path_handle("5");
        path_handle_t p6 = graph.create_path_handle("6");
        path_handle_t p7 = graph.create_path_handle("7");
        path_handle_t p8 = graph.create_path_handle("8");
        path_handle_t p9 = graph.create_path_handle("9");
        path_handle_t p10 = graph.create_path_handle("10");
        path_handle_t p11 = graph.create_path_handle("11");
        path_handle_t p12 = graph.create_path_handle("12");
        path_handle_t p13 = graph.create_path_handle("13");
        for (int i = 0; i < 256; ++i) {
            handle_t h = graph.create_handle("A");
            if (i != 0 &&
                i != 64 &&
                i != 128 &&
                i != 192 &&
                i != 208 &&
                i != 224 &&
                i != 240) {
                graph.create_edge(prev, h);
            }
            
            // first component: separated by a spacer
            if (i < 32) {
                graph.append_step(p1, h);
            }
            if (i >= 48 && i < 64) {
                graph.append_step(p2, h);
            }
            // second component: overlapping
            if (i >= 64 && i < 100) {
                graph.append_step(p3, h);
            }
            if (i >= 80 && i < 128) {
                graph.append_step(p4, h);
            }
            // third component: large path with small overlapping paths
            if (i >= 128 && i < 192) {
                graph.append_step(p5, h);
            }
            if (i >= 128 && i < 140) {
                graph.append_step(p6, h);
            }
            if (i >= 150 && i < 165) {
                graph.append_step(p7, h);
            }
            if (i >= 160 && i < 180) {
                graph.append_step(p8, h);
            }
            if (i >= 185 && i < 192) {
                graph.append_step(p9, h);
            }
            // several smaller components with one path each
            if (i >= 192 && i < 208) {
                graph.append_step(p10, h);
            }
            if (i >= 208 && i < 224) {
                graph.append_step(p11, h);
            }
            if (i >= 224 && i < 240) {
                graph.append_step(p12, h);
            }
            if (i >= 240) {
                graph.append_step(p13, h);
            }
            prev = h;
        }
        
        auto serial_result = normalize(algorithms::component_paths(graph));
        
        int repetitions = 10;
        for (int num_threads : {1, 2, 4, 8, 16}) {
            omp_set_num_threads(num_threads);
            for (int i = 0; i < repetitions; ++i) {
                auto parallel_result = normalize(algorithms::component_paths_parallel(graph));
                REQUIRE(parallel_result == serial_result);
            }
        }
    }
    
    SECTION("On a graph with many small components") {

        bdsg::HashGraph graph;

        for (int i = 0; i < 4096; ++i) {
            handle_t h = graph.create_handle("A");
            path_handle_t p = graph.create_path_handle(to_string(i));
            graph.append_step(p, h);
        }

        auto serial_result = normalize(algorithms::component_paths(graph));

        int repetitions = 10;
        for (int num_threads : {1, 2, 4, 8, 16}) {
            omp_set_num_threads(num_threads);
            for (int i = 0; i < repetitions; ++i) {
                auto parallel_result = normalize(algorithms::component_paths_parallel(graph));
                REQUIRE(parallel_result == serial_result);
            }
        }
    }
    
    omp_set_num_threads(thread_count_pre);
}

// TODO: test that will require second int vector
        
}
}
