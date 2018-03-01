/// \file multipath_alignment.cpp
///  
/// unit tests for multipath alignment construction and utility functions
///

#include <stdio.h>
#include <iostream>
#include "json2pb.h"
#include "vg.pb.h"
#include "vg.hpp"
#include "multipath_alignment.hpp"
#include "catch.hpp"
#include "utility.hpp"

namespace vg {
    namespace unittest {
        
        TEST_CASE( "Multipath alignments correctly identify source subpaths",
                  "[alignment][multipath][mapping]" ) {
            
            SECTION( "Multipath alignment can identify source subpath in a linear subpath structure") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("GCATCTGA");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                Mapping* mapping1 = subpath0->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                Mapping* mapping2 = subpath1->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(4);
                
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(multipath_aln.start_size() == 1);
                REQUIRE(multipath_aln.start(0) == 0);
            }
            
            SECTION( "Multipath alignment can identify source subpath in a forked subpath structure") {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("GCATCTGA");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(3);
                subpath2->add_next(3);
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                
                Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                
                
                identify_start_subpaths(multipath_aln);
                
                REQUIRE(multipath_aln.start_size() == 2);
                bool found_1 = false;
                bool found_2 = false;
                
                
                REQUIRE(multipath_aln.start_size() == 2);
                
                for (int i = 0; i < multipath_aln.start_size(); i++) {
                    if (multipath_aln.start(i) == 0) {
                        found_1 = true;
                    }
                    else if (multipath_aln.start(i) == 2) {
                        found_2 = true;
                    }
                    else {
                        REQUIRE(false);
                    }
                }
                REQUIRE( found_1 );
                REQUIRE( found_2 );
            }
        }
        
        TEST_CASE( "Multipath alignments correctly identifies optimal alignment within subpath DAG",
                  "[alignment][multipath][mapping]" ) {
            
            SECTION( "Multipath alignment can identify optimal alignment between two disjoint paths" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("T");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                
                // set scores
                subpath0->set_score(1);
                subpath1->set_score(0);
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(2);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(3);
                
                // get optimal alignment
                identify_start_subpaths(multipath_aln);
                Alignment aln;
                optimal_alignment(multipath_aln, aln);
                
                REQUIRE(aln.path().mapping_size() == multipath_aln.subpath(0).path().mapping_size());
                for (int i = 0; i < aln.path().mapping_size(); i++) {
                    REQUIRE(aln.path().mapping(i).position().node_id() == multipath_aln.subpath(0).path().mapping(i).position().node_id());
                }
                REQUIRE(aln.score() == multipath_aln.subpath(0).score());
                
            }
            
            SECTION( "Multipath alignment can identify optimal alignment among paths that intersect" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                
                string read = string("GCAGCTGA");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                Subpath* subpath4 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(2);
                subpath1->add_next(2);
                subpath2->add_next(3);
                subpath2->add_next(4);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(1);
                subpath3->set_score(0);
                subpath4->set_score(4);
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                
                Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(4);
                
                Mapping* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(5);
                
                // get optimal alignment
                identify_start_subpaths(multipath_aln);
                Alignment aln;
                optimal_alignment(multipath_aln, aln);
                
                // follows correct path
                REQUIRE(aln.path().mapping_size() == 3);
                REQUIRE(aln.path().mapping(0).position().node_id() == 1);
                REQUIRE(aln.path().mapping(1).position().node_id() == 3);
                REQUIRE(aln.path().mapping(2).position().node_id() == 5);
                
                // has correct score
                REQUIRE(aln.score() == 8);
                
            }
            
            SECTION( "Multipath alignment correctly merge Mappings while finding optimal alignment" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("GTTGA");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                
                string read = string("GCAGTGACTGA");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                Subpath* subpath4 = multipath_aln.add_subpath();
                Subpath* subpath5 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(2);
                subpath1->add_next(2);
                subpath2->add_next(3);
                subpath3->add_next(4);
                subpath3->add_next(5);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(3);
                subpath3->set_score(2);
                subpath4->set_score(0);
                subpath5->set_score(4);
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                Edit* edit2 = mapping2->add_edit();
                edit2->set_from_length(2);
                edit2->set_to_length(2);
                
                Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(3);
                mapping3->mutable_position()->set_offset(2);
                Edit* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                
                Mapping* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(4);
                
                Mapping* mapping5 = subpath5->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(5);
                
                // get optimal alignment
                identify_start_subpaths(multipath_aln);
                Alignment aln;
                optimal_alignment(multipath_aln, aln);
                
                // follows correct path
                REQUIRE(aln.path().mapping_size() == 3);
                REQUIRE(aln.path().mapping(0).position().node_id() == 1);
                REQUIRE(aln.path().mapping(1).position().node_id() == 3);
                REQUIRE(aln.path().mapping(2).position().node_id() == 5);
                
                // has correct ranks
                REQUIRE(aln.path().mapping(0).rank() == 1);
                REQUIRE(aln.path().mapping(1).rank() == 2);
                REQUIRE(aln.path().mapping(2).rank() == 3);
                
                // has correct score
                REQUIRE(aln.score() == 12);
                
                // has correct edits
                REQUIRE(aln.path().mapping(1).position().offset() == 0);
                REQUIRE(aln.path().mapping(1).edit(0).from_length() == 5);
                REQUIRE(aln.path().mapping(1).edit(0).to_length() == 5);
            }
        }
        
        TEST_CASE("Multipath alignment correctly identifies suboptimal alignments") {
        
            SECTION( "Multipath alignment can identify optimal and suboptimal alignment between two disjoint paths" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("G");
                Node* n4 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n2);
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n4);
                graph.create_edge(n3, n4);
                
                string read = string("T");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                
                // set scores
                subpath0->set_score(1);
                subpath1->set_score(0);
                
                // designate mappings
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(2);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(3);
                
                // get top 10 alignments
                auto top10 = optimal_alignments(multipath_aln, 10);
                
                // There should be two alignments
                REQUIRE(top10.size() == 2);
                
                // The best alignment should be the same one we got for this example before as optimal
                REQUIRE(top10[0].path().mapping_size() == multipath_aln.subpath(0).path().mapping_size());
                for (int i = 0; i < top10[0].path().mapping_size(); i++) {
                    REQUIRE(top10[0].path().mapping(i).position().node_id() == multipath_aln.subpath(0).path().mapping(i).position().node_id());
                }
                REQUIRE(top10[0].score() == multipath_aln.subpath(0).score());
                
                // The second best alignment should be that 0-score subpath
                REQUIRE(top10[1].path().mapping_size() == multipath_aln.subpath(1).path().mapping_size());
                for (int i = 0; i < top10[1].path().mapping_size(); i++) {
                    REQUIRE(top10[1].path().mapping(i).position().node_id() == multipath_aln.subpath(1).path().mapping(i).position().node_id());
                }
                REQUIRE(top10[1].score() == multipath_aln.subpath(1).score());
                
            }
            
            SECTION( "Multipath alignment correctly merge Mappings while finding optimal and suboptimal alignments" ) {
                
                VG graph;
                
                Node* n1 = graph.create_node("GCA");
                Node* n2 = graph.create_node("T");
                Node* n3 = graph.create_node("GTTGA");
                Node* n4 = graph.create_node("A");
                Node* n5 = graph.create_node("CTGA");
                
                graph.create_edge(n1, n3);
                graph.create_edge(n2, n3);
                graph.create_edge(n3, n4);
                graph.create_edge(n3, n5);
                
                string read = string("GCAGTGACTGA");
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                Subpath* subpath4 = multipath_aln.add_subpath();
                Subpath* subpath5 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(2);
                subpath1->add_next(2);
                subpath2->add_next(3);
                subpath3->add_next(4);
                subpath3->add_next(5);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(3);
                subpath3->set_score(2);
                subpath4->set_score(0);
                subpath5->set_score(4);
                
                // designate mappings
                // TODO: Note that these edits aren't quite realistic
                Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
                mapping0->mutable_position()->set_node_id(1);
                
                Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
                mapping1->mutable_position()->set_node_id(2);
                
                Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
                mapping2->mutable_position()->set_node_id(3);
                Edit* edit2 = mapping2->add_edit();
                edit2->set_from_length(2);
                edit2->set_to_length(2);
                
                Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
                mapping3->mutable_position()->set_node_id(3);
                mapping3->mutable_position()->set_offset(2);
                Edit* edit3 = mapping3->add_edit();
                edit3->set_from_length(3);
                edit3->set_to_length(3);
                
                Mapping* mapping4 = subpath4->mutable_path()->add_mapping();
                mapping4->mutable_position()->set_node_id(4);
                
                Mapping* mapping5 = subpath5->mutable_path()->add_mapping();
                mapping5->mutable_position()->set_node_id(5);
                
                // get top 10 alignments
                identify_start_subpaths(multipath_aln);
                vector<Alignment> top10 = optimal_alignments(multipath_aln, 10);
                
                SECTION("Top alignment is correct") {
                    // Exists
                    REQUIRE(top10.size() > 0);
                    
                    // follows correct path
                    REQUIRE(top10[0].path().mapping_size() == 3);
                    REQUIRE(top10[0].path().mapping(0).position().node_id() == 1);
                    REQUIRE(top10[0].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[0].path().mapping(2).position().node_id() == 5);
                    
                    // has correct score
                    REQUIRE(top10[0].score() == 12);
                    
                    // has correct edits
                    REQUIRE(top10[0].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[0].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[0].path().mapping(1).edit(0).to_length() == 5);
                }
                
                SECTION("Secondary alignment is correct") {
                    
                    // Exists
                    REQUIRE(top10.size() > 1);
                    
                    // follows correct path
                    REQUIRE(top10[1].path().mapping_size() == 3);
                    REQUIRE(top10[1].path().mapping(0).position().node_id() == 2);
                    REQUIRE(top10[1].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[1].path().mapping(2).position().node_id() == 5);
                    
                    // has correct score
                    REQUIRE(top10[1].score() == 9);
                    
                    // has correct edits
                    REQUIRE(top10[1].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[1].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[1].path().mapping(1).edit(0).to_length() == 5);
                
                }
                
                SECTION("Tertiary alignment is correct") {
                    
                    // Exists
                    REQUIRE(top10.size() > 2);
                    
                    // follows correct path
                    REQUIRE(top10[2].path().mapping_size() == 3);
                    REQUIRE(top10[2].path().mapping(0).position().node_id() == 1);
                    REQUIRE(top10[2].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[2].path().mapping(2).position().node_id() == 4);
                    
                    // has correct score
                    REQUIRE(top10[2].score() == 8);
                    
                    // has correct edits
                    REQUIRE(top10[2].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[2].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[2].path().mapping(1).edit(0).to_length() == 5);
                
                }
                
                SECTION("Quaternary alignment is correct") {
                    
                    // Exists
                    REQUIRE(top10.size() > 3);
                    
                    // follows correct path
                    REQUIRE(top10[3].path().mapping_size() == 3);
                    REQUIRE(top10[3].path().mapping(0).position().node_id() == 2);
                    REQUIRE(top10[3].path().mapping(1).position().node_id() == 3);
                    REQUIRE(top10[3].path().mapping(2).position().node_id() == 4);
                    
                    // has correct score
                    REQUIRE(top10[3].score() == 5);
                    
                    // has correct edits
                    REQUIRE(top10[3].path().mapping(1).position().offset() == 0);
                    REQUIRE(top10[3].path().mapping(1).edit(0).from_length() == 5);
                    REQUIRE(top10[3].path().mapping(1).edit(0).to_length() == 5);
                
                }
                
                SECTION("Quinary alignment does not exist") {
                    REQUIRE(top10.size() < 5);
                }
            }
        
        }
        
        TEST_CASE( "Reverse complementing multipath alignments works correctly",
                  "[alignment][multipath][mapping]" ) {
            
            VG graph;
            
            Node* n1 = graph.create_node("GCA");
            Node* n2 = graph.create_node("T");
            Node* n3 = graph.create_node("GCC");
            Node* n4 = graph.create_node("A");
            Node* n5 = graph.create_node("CTGA");
            
            graph.create_edge(n1, n3);
            graph.create_edge(n2, n3);
            graph.create_edge(n3, n4);
            graph.create_edge(n3, n5);
            
            string read = string("CACCCTGA");
            MultipathAlignment multipath_aln;
            multipath_aln.set_sequence(read);
            
            // add subpaths
            Subpath* subpath0 = multipath_aln.add_subpath();
            Subpath* subpath1 = multipath_aln.add_subpath();
            Subpath* subpath2 = multipath_aln.add_subpath();
            Subpath* subpath3 = multipath_aln.add_subpath();
            Subpath* subpath4 = multipath_aln.add_subpath();
            Subpath* subpath5 = multipath_aln.add_subpath();
            
            // set edges between subpaths
            subpath0->add_next(2);
            subpath0->add_next(5);
            subpath1->add_next(2);
            subpath2->add_next(3);
            subpath2->add_next(4);
            
            // set scores
            subpath0->set_score(3);
            subpath1->set_score(0);
            subpath2->set_score(1);
            subpath3->set_score(0);
            subpath4->set_score(4);
            subpath5->set_score(0);
            
            // designate mappings
            Mapping* mapping0 = subpath0->mutable_path()->add_mapping();
            mapping0->mutable_position()->set_node_id(1);
            mapping0->mutable_position()->set_offset(1);
            Edit* edit00 = mapping0->add_edit();
            edit00->set_from_length(2);
            edit00->set_to_length(2);
            
            Mapping* mapping1 = subpath1->mutable_path()->add_mapping();
            mapping1->mutable_position()->set_node_id(2);
            mapping1->mutable_position()->set_offset(1);
            Edit* edit10 = mapping1->add_edit();
            edit10->set_from_length(0);
            edit10->set_to_length(2);
            edit10->set_sequence("CA");
            
            Mapping* mapping2 = subpath2->mutable_path()->add_mapping();
            mapping2->mutable_position()->set_node_id(3);
            Edit* edit20 = mapping2->add_edit();
            edit20->set_from_length(1);
            edit20->set_to_length(0);
            Edit* edit21 = mapping2->add_edit();
            edit21->set_from_length(2);
            edit21->set_to_length(2);
            
            Mapping* mapping3 = subpath3->mutable_path()->add_mapping();
            mapping3->mutable_position()->set_node_id(4);
            Edit* edit30 = mapping3->add_edit();
            edit30->set_from_length(0);
            edit30->set_to_length(4);
            edit30->set_sequence("CTGA");
            
            Mapping* mapping4 = subpath4->mutable_path()->add_mapping();
            mapping4->mutable_position()->set_node_id(5);
            Edit* edit40 = mapping4->add_edit();
            edit40->set_from_length(4);
            edit40->set_to_length(4);
            
            Mapping* mapping5 = subpath5->mutable_path()->add_mapping();
            mapping5->mutable_position()->set_node_id(3);
            mapping5->mutable_position()->set_offset(3);
            Edit* edit50 = mapping5->add_edit();
            edit50->set_from_length(0);
            edit50->set_to_length(6);
            edit50->set_sequence("CCCTGA");
            
            // identify starts
            multipath_aln.add_start(0);
            multipath_aln.add_start(1);
            
            MultipathAlignment rc_multipath_aln;
            
            auto node_length = [&graph](int64_t node_id) { return graph.get_node(node_id)->sequence().length(); };
            rev_comp_multipath_alignment(multipath_aln, node_length, rc_multipath_aln);
            
            // all subpaths preserved
            REQUIRE(rc_multipath_aln.subpath_size() == multipath_aln.subpath_size());
            
            // subpaths are in reversed order
            for (int i = 0; i < multipath_aln.subpath_size(); i++) {
                REQUIRE(rc_multipath_aln.subpath(i).score() == multipath_aln.subpath(multipath_aln.subpath_size() - i - 1).score());
                REQUIRE(rc_multipath_aln.subpath(i).path().mapping(0).position().node_id() == multipath_aln.subpath(multipath_aln.subpath_size() - i - 1).path().mapping(0).position().node_id());
            }
            
            // positions are on reverse
            for (int i = 0; i < rc_multipath_aln.subpath_size(); i++) {
                for (int j = 0; j < rc_multipath_aln.subpath(i).path().mapping_size(); j++) {
                    REQUIRE(rc_multipath_aln.subpath(i).path().mapping(j).position().is_reverse());
                }
            }
            
            // offsets set correctly
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(1).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).position().offset() == 1);
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).position().offset() == 0);
            REQUIRE(rc_multipath_aln.subpath(5).path().mapping(0).position().offset() == 0);
            
            // edits are correct
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).edit(0).to_length() == 6);
            REQUIRE(rc_multipath_aln.subpath(0).path().mapping(0).edit(0).sequence() == "TCAGGG");
            
            REQUIRE(rc_multipath_aln.subpath(1).path().mapping(0).edit(0).from_length() == 4);
            REQUIRE(rc_multipath_aln.subpath(1).path().mapping(0).edit(0).to_length() == 4);
            
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).edit(0).to_length() == 4);
            REQUIRE(rc_multipath_aln.subpath(2).path().mapping(0).edit(0).sequence() == "TCAG");
            
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(0).from_length() == 2);
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(0).to_length() == 2);
            
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(1).from_length() == 1);
            REQUIRE(rc_multipath_aln.subpath(3).path().mapping(0).edit(1).to_length() == 0);
            
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).edit(0).from_length() == 0);
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).edit(0).to_length() == 2);
            REQUIRE(rc_multipath_aln.subpath(4).path().mapping(0).edit(0).sequence() == "TG");
            
            REQUIRE(rc_multipath_aln.subpath(5).path().mapping(0).edit(0).from_length() == 2);
            REQUIRE(rc_multipath_aln.subpath(5).path().mapping(0).edit(0).to_length() == 2);
            
            // starts are switched correctly
            for (int i : {0, 1, 2}) {
                bool found = false;
                for (int64_t j = 0; j < rc_multipath_aln.start_size(); j++) {
                    found = found || rc_multipath_aln.start(j) == i;
                }
                REQUIRE(found);
            }
            
            // now check that the in place reverse complement works correctly
            // (it should match the already verified reverse complement)
            rev_comp_multipath_alignment_in_place(&multipath_aln, node_length);
            
            REQUIRE(rc_multipath_aln.subpath_size() == multipath_aln.subpath_size());
            
            for (int64_t i = 0; i < multipath_aln.subpath_size(); i++) {
                const Subpath& subpath_1 = multipath_aln.subpath(i);
                const Subpath& subpath_2 = rc_multipath_aln.subpath(i);
                REQUIRE(subpath_1.score() == subpath_2.score());
                const Path& path_1 = subpath_1.path();
                const Path& path_2 = subpath_2.path();
                REQUIRE(path_1.mapping_size() == path_2.mapping_size());
                for (int64_t j = 0; j < path_1.mapping_size(); j++) {
                    const Mapping& mapping_1 = path_1.mapping(j);
                    const Mapping& mapping_2 = path_2.mapping(j);
                    REQUIRE(mapping_1.edit_size() == mapping_2.edit_size());
                    for (int64_t k = 0; k < mapping_1.edit_size(); k++) {
                        const Edit& edit_1 = mapping_1.edit(k);
                        const Edit& edit_2 = mapping_2.edit(k);
                        REQUIRE(edit_1.from_length() == edit_2.from_length());
                        REQUIRE(edit_1.to_length() == edit_2.to_length());
                        REQUIRE(edit_1.sequence() == edit_2.sequence());
                    }
                }
                REQUIRE(subpath_1.next_size() == subpath_2.next_size());
                // these might be out of order
                for (int64_t j = 0; j < subpath_1.next_size(); j++) {
                    bool found = false;
                    for (int64_t k = 0; k < subpath_2.next_size(); k++) {
                        found = found || subpath_1.next(j) == subpath_2.next(k);
                    }
                    REQUIRE(found);
                }
            }
            
            // starts are switched correctly
            REQUIRE(rc_multipath_aln.start_size() == multipath_aln.start_size());
            for (int64_t i = 0; i < rc_multipath_aln.start_size(); i++) {
                bool found = false;
                for (int64_t j = 0; j < multipath_aln.start_size(); j++) {
                    found = found || rc_multipath_aln.start(i) == multipath_aln.start(j);
                }
                REQUIRE(found);
            }
        }
        
        TEST_CASE( "Algorithm returns correct connected components for MultipathAlignments",
                  "[alignment][multipath]" ){
            
            SECTION("Works for a single component MultipathAlignment") {
                
                
                string read = "CACCCTGA";
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(2);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(4);
                subpath2->set_score(7);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 1);
                
                bool found_1 = false;
                for (size_t i = 0; i < comps.size(); i++) {
                    bool this_is_1 = true;
                    for (int j : {0, 1, 2}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        this_is_1 = this_is_1 && found_this_item;
                    }
                    found_1 = this_is_1 || found_1;
                }
                
                REQUIRE(found_1);
            }
            
            SECTION("Works for a two component MultipathAlignment") {
                
                string read = "CACCCTGA";
                MultipathAlignment multipath_aln;
                multipath_aln.set_sequence(read);
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(4);
                subpath2->set_score(7);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 2);
                
                bool found_1 = false;
                bool found_2 = false;
                for (size_t i = 0; i < comps.size(); i++) {
                    bool is_this_comp;
                    is_this_comp = true;
                    for (int j : {0, 1}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_1 = is_this_comp || found_1;
                    
                    is_this_comp = true;
                    for (int j : {2}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_2 = is_this_comp || found_2;
                }
                
                REQUIRE(found_1);
                REQUIRE(found_2);
            }
            
            SECTION("Works for a multi-component MultipathAlignment with complicated edge structure") {
                
                MultipathAlignment multipath_aln;
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                Subpath* subpath4 = multipath_aln.add_subpath();
                Subpath* subpath5 = multipath_aln.add_subpath();
                Subpath* subpath6 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(5);
                subpath2->add_next(4);
                subpath3->add_next(5);
                
                // set scores
                subpath0->set_score(3);
                subpath1->set_score(0);
                subpath2->set_score(1);
                subpath3->set_score(0);
                subpath4->set_score(4);
                subpath5->set_score(0);
                subpath5->set_score(9);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 3);
                
                bool found_1 = false;
                bool found_2 = false;
                bool found_3 = false;
                for (size_t i = 0; i < comps.size(); i++) {
                    bool is_this_comp;
                    is_this_comp = true;
                    for (int j : {0, 1, 3, 5}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_1 = is_this_comp || found_1;
                    
                    is_this_comp = true;
                    for (int j : {2, 4}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_2 = is_this_comp || found_2;
                    
                    is_this_comp = true;
                    for (int j : {6}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    found_3 = is_this_comp || found_3;
                }
                
                REQUIRE(found_1);
                REQUIRE(found_2);
                REQUIRE(found_3);
            }
        }
        
        TEST_CASE("We can extract subgraphs of a MultipathAlignment",
                  "[alignment][multipath]" ) {
            
            SECTION("Works correctly when splitting by connected components") {
                
                MultipathAlignment multipath_aln;
                
                // add subpaths
                Subpath* subpath0 = multipath_aln.add_subpath();
                Subpath* subpath1 = multipath_aln.add_subpath();
                Subpath* subpath2 = multipath_aln.add_subpath();
                Subpath* subpath3 = multipath_aln.add_subpath();
                Subpath* subpath4 = multipath_aln.add_subpath();
                Subpath* subpath5 = multipath_aln.add_subpath();
                Subpath* subpath6 = multipath_aln.add_subpath();
                
                // set edges between subpaths
                subpath0->add_next(1);
                subpath1->add_next(5);
                subpath2->add_next(4);
                subpath3->add_next(5);
                
                // set scores (hijacking them here to label subpaths with their original index)
                subpath0->set_score(0);
                subpath1->set_score(1);
                subpath2->set_score(2);
                subpath3->set_score(3);
                subpath4->set_score(4);
                subpath5->set_score(5);
                subpath5->set_score(6);
                
                auto comps = connected_components(multipath_aln);
                
                REQUIRE(comps.size() == 3);
                
                for (size_t i = 0; i < comps.size(); i++) {
                    
                    MultipathAlignment sub;
                    extract_sub_multipath_alignment(multipath_aln, comps[i], sub);
                                        
                    bool is_this_comp;
                    bool finds_all_subpaths;
                    
                    is_this_comp = true;
                    for (int j : {0, 1, 3, 5}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    
                    if (is_this_comp) {
                        
                        REQUIRE(sub.subpath_size() == 4);
                        
                        finds_all_subpaths = true;
                        for (int j : {0, 1, 3, 5}) {
                            bool found_this_item = false;
                            for (size_t k = 0; k < sub.subpath_size(); k++) {
                                if (sub.subpath(k).score() == j) {
                                    REQUIRE(sub.subpath(k).next_size() == multipath_aln.subpath(j).next_size());
                                    for (size_t l = 0; l < sub.subpath(k).next_size(); l++) {
                                        bool found_this_next = false;
                                        for (size_t m = 0; m < multipath_aln.subpath(j).next_size(); m++) {
                                            if (sub.subpath(sub.subpath(k).next(l)).score() == multipath_aln.subpath(j).next(m)) {
                                                found_this_next = true;
                                            }
                                        }
                                        found_this_item = found_this_item && found_this_next;
                                    }
                                }
                            }
                            finds_all_subpaths = finds_all_subpaths && found_this_item;
                        }
                    }
                    
                    is_this_comp = true;
                    for (int j : {2, 4}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    
                    if (is_this_comp) {
                        
                        REQUIRE(sub.subpath_size() == 2);
                        
                        finds_all_subpaths = true;
                        for (int j : {2, 4}) {
                            bool found_this_item = false;
                            for (size_t k = 0; k < sub.subpath_size(); k++) {
                                if (sub.subpath(k).score() == j) {
                                    REQUIRE(sub.subpath(k).next_size() == multipath_aln.subpath(j).next_size());
                                    for (size_t l = 0; l < sub.subpath(k).next_size(); l++) {
                                        bool found_this_next = false;
                                        for (size_t m = 0; m < multipath_aln.subpath(j).next_size(); m++) {
                                            if (sub.subpath(sub.subpath(k).next(l)).score() == multipath_aln.subpath(j).next(m)) {
                                                found_this_next = true;
                                            }
                                        }
                                        found_this_item = found_this_item && found_this_next;
                                    }
                                }
                            }
                            finds_all_subpaths = finds_all_subpaths && found_this_item;
                        }
                    }
                    
                    is_this_comp = true;
                    for (int j : {6}) {
                        bool found_this_item = false;
                        for (size_t k = 0; k < comps[i].size(); k++) {
                            if (comps[i][k] == j) {
                                found_this_item = true;
                                break;
                            }
                        }
                        is_this_comp = is_this_comp && found_this_item;
                    }
                    
                    if (is_this_comp) {
                        
                        REQUIRE(sub.subpath_size() == 1);
                        
                        finds_all_subpaths = true;
                        for (int j : {6}) {
                            bool found_this_item = false;
                            for (size_t k = 0; k < sub.subpath_size(); k++) {
                                if (sub.subpath(k).score() == j) {
                                    REQUIRE(sub.subpath(k).next_size() == multipath_aln.subpath(j).next_size());
                                    for (size_t l = 0; l < sub.subpath(k).next_size(); l++) {
                                        bool found_this_next = false;
                                        for (size_t m = 0; m < multipath_aln.subpath(j).next_size(); m++) {
                                            if (sub.subpath(sub.subpath(k).next(l)).score() == multipath_aln.subpath(j).next(m)) {
                                                found_this_next = true;
                                            }
                                        }
                                        found_this_item = found_this_item && found_this_next;
                                    }
                                }
                            }
                            finds_all_subpaths = finds_all_subpaths && found_this_item;
                        }
                    }
                }
            }
        }
    }
}



