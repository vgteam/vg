//
//  suffix_tree.cpp
//  
//
//  Created by Jordan Eizenga on 1/25/17.
//
//

#include <stdio.h>
#include <iostream>
#include <unordered_map>

#include "msa_converter.hpp"
#include "catch.hpp"

using namespace std;

namespace vg {
    namespace unittest {
        
        TEST_CASE("MSAConverter can build from each kind of input", "[msa]") {
            
            
            SECTION("MSAConverter can build from FASTA input") {
                
                string input = ">seq1 description 1\nAAGTGATAGAGATAAA\nGTGATAGAGATA\n>seq2 description2\n-AAGTATAGACATAAA\nGTTAGAGATA--\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                VG graph = msa_converter.make_graph();
            }
            
            SECTION("MSAConverter can build from Clustal input") {
                
                string input = "CLUSTAL O(1.2.4) multiple sequence alignment\n\nGI262359905      TGTACCTTGATTTCGTATTCTGAGAGGCTGCTGCTTAGCGGTAGCCCCT-TGGTTTCCGT\nGI528476558      ---------------------------GTGGAAGTGTTTGCTACCAAGTTTATTTGCAGT\nref              ---------------------------GTGGAAGTGTTTGCTACCAAGTTTATTTGCAGT\n                **    *    * ** *   * *  ** * **\n\nGI262359905      GGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTGCGCGGCGTGAGCT\nGI528476558      GTTAACAGCACAACATTTACAAAACGTATTTTGTACAATCAAGTCTTCACTGCCCT----\nref              GTTAACAGCACAACATTTACAAAACGTATTTTGTACAATCAAGTCTTCACTGCCCT----\n                *  *** * * * *      *  **  **    ** **    * ** * * **                \n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "clustal");
                
                VG graph = msa_converter.make_graph();
            }
        }
        
        TEST_CASE("MSAConverter produces correct graphs from MSAs", "[msa]") {
            
            SECTION("MSAConverter produces one node from a completely matching alignment") {
                
                string input = ">seq1\nAAA\n>seq2\nAAA\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                VG graph = msa_converter.make_graph();
                
                Graph& g = graph.graph;
                
                REQUIRE(g.node_size() == 1);
                REQUIRE(g.edge_size() == 0);
                REQUIRE(g.node(0).sequence() == "AAA");
            }
            
            SECTION("MSAConverter respects the max node length") {
                
                string input = ">seq1\nAAA\n>seq2\nAAA\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                VG graph = msa_converter.make_graph(true, 1);
                
                Graph& g = graph.graph;
                
                REQUIRE(g.node_size() == 3);
                REQUIRE(g.edge_size() == 2);
            }
            
            SECTION("MSAConverter splits columns that mismatch into multiple nodes") {
                
                string input = ">seq1\nATG\n>seq2\nACG\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                VG graph = msa_converter.make_graph();
                
                Graph& g = graph.graph;
                
                unordered_map<string, int64_t> nodes;
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    nodes[n.sequence()] = n.id();
                }
                
                REQUIRE(g.node_size() == 4);
                REQUIRE(g.edge_size() == 4);
                
                REQUIRE(nodes.count("A"));
                REQUIRE(nodes.count("C"));
                REQUIRE(nodes.count("G"));
                REQUIRE(nodes.count("T"));
                
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["C"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["T"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["T"], true), NodeSide(nodes["G"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["C"], true), NodeSide(nodes["G"], false)));
            }
            
            SECTION("MSAConverter adds edges over gap") {
                
                string input = ">seq1\nA-G\n>seq2\nACG\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                VG graph = msa_converter.make_graph();
                
                Graph& g = graph.graph;
                
                unordered_map<string, int64_t> nodes;
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    nodes[n.sequence()] = n.id();
                }
                
                REQUIRE(g.node_size() == 3);
                REQUIRE(g.edge_size() == 3);
                
                REQUIRE(nodes.count("A"));
                REQUIRE(nodes.count("C"));
                REQUIRE(nodes.count("G"));
                
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["C"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["G"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["C"], true), NodeSide(nodes["G"], false)));
            }
            
            SECTION("MSAConverter handles overlapping gaps") {
                
                string input = ">seq1\nAA--GTT\n>seq2\nAAACGTT\n>seq3\nAAA--TT\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                VG graph = msa_converter.make_graph();
                
                Graph& g = graph.graph;
                
                unordered_map<string, int64_t> nodes;
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    nodes[n.sequence()] = n.id();
                }
                
                REQUIRE(g.node_size() == 5);
                REQUIRE(g.edge_size() == 6);
                
                REQUIRE(nodes.count("AA"));
                REQUIRE(nodes.count("A"));
                REQUIRE(nodes.count("C"));
                REQUIRE(nodes.count("G"));
                REQUIRE(nodes.count("TT"));
                
                REQUIRE(graph.has_edge(NodeSide(nodes["AA"], true), NodeSide(nodes["A"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["AA"], true), NodeSide(nodes["G"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["C"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["TT"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["C"], true), NodeSide(nodes["G"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["G"], true), NodeSide(nodes["TT"], false)));
            }
            
            SECTION("MSAConverter handles nested gaps") {
                
                string input = ">seq1\nAAACGTT\n>seq2\nAA---TT\n>seq3\nAAA-GTT\n";
                
                istringstream strm(input);
                
                MSAConverter msa_converter;
                msa_converter.load_alignments(strm, "fasta");
                
                VG graph = msa_converter.make_graph();
                
                Graph& g = graph.graph;
                
                unordered_map<string, int64_t> nodes;
                for (size_t i = 0; i < g.node_size(); i++) {
                    const Node& n = g.node(i);
                    nodes[n.sequence()] = n.id();
                }
                
                REQUIRE(g.node_size() == 5);
                REQUIRE(g.edge_size() == 6);
                
                REQUIRE(nodes.count("AA"));
                REQUIRE(nodes.count("A"));
                REQUIRE(nodes.count("C"));
                REQUIRE(nodes.count("G"));
                REQUIRE(nodes.count("TT"));
                
                REQUIRE(graph.has_edge(NodeSide(nodes["AA"], true), NodeSide(nodes["A"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["AA"], true), NodeSide(nodes["TT"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["C"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["A"], true), NodeSide(nodes["G"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["C"], true), NodeSide(nodes["G"], false)));
                REQUIRE(graph.has_edge(NodeSide(nodes["G"], true), NodeSide(nodes["TT"], false)));
            }
        }
    }
}
