/// \file annotation.cpp
///  
/// Unit tests for the annotations system on Alignments and MultipathAlignments
///

#include <iostream>
#include <string>
#include "../json2pb.h"
#include "../vg.pb.h"
#include "../annotation.hpp"
#include "catch.hpp"


namespace vg {
namespace unittest {
using namespace std;

TEST_CASE("Annotations can be populated", "[alignment][annotation]") {
    
    Alignment aln;
   
    SECTION("Referencing an annotation in the mutable map can create it") {
   
        string key_name = "cats";
        string key2_name = "dogs";
        
        // It's undefined behavior to have a completely empty Value in the annotations Struct.
        // So put in a NullValue-holding Value.
        (*aln.mutable_annotation()->mutable_fields())[key_name].set_null_value(google::protobuf::NullValue::NULL_VALUE);
        
        (*aln.mutable_annotation()->mutable_fields())[key2_name].set_number_value(1000);
        
        REQUIRE(aln.annotation().fields().at(key_name).kind_case() == google::protobuf::Value::KindCase::kNullValue);
        REQUIRE(aln.annotation().fields().at(key2_name).number_value() == 1000);
        
        set_annotation(&aln, "exploding", true);
        
        REQUIRE(get_annotation<double>(aln, "dogs") == 1000.0);
        REQUIRE(get_annotation<bool>(aln, "exploding") == true);
       
        string serialized = pb2json(aln);
        
        // Serialize and deserialize
        Alignment aln2;
        json2pb(aln2, serialized.c_str(), serialized.size());
       
        REQUIRE(aln2.annotation().fields().count(key_name));
        REQUIRE(aln2.annotation().fields().count(key2_name));
    }
    
}
   
}
}
        
