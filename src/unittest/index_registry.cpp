/// \file indexed_vg.cpp
///  
/// unit tests for the vg-file-backed handle graph implementation

#include <iostream>
#include "../index_registry.hpp"
#include "catch.hpp"

namespace vg {
namespace unittest {

// expose the private functions to testing
class TestIndexRegistry : public IndexRegistry {
public:
    using IndexRegistry::dependency_order;
    using IndexRegistry::make_plan;
    using IndexRegistry::get_index;
};

using namespace std;
    
TEST_CASE("IndexRegistry works on a dummy dependency graph", "[indexregistry]") {
    
    TestIndexRegistry registry;
    
    // name the indexes
    registry.register_index("FASTA");
    registry.register_index("VCF");
    registry.register_index("GFA");
    registry.register_index("VG");
    registry.register_index("Pruned VG");
    registry.register_index("XG");
    registry.register_index("GCSA+LCP");
    
    // make some dummy recipes that don't actually do anything
    registry.register_recipe("VG", {"FASTA", "VCF"},
                             [&] (const vector<const IndexFile*>& inputs) {
        vector<string> filename(1, "vg-file");
        return filename;
    });
    registry.register_recipe("VG", {"GFA"},
                             [&] (const vector<const IndexFile*>& inputs) {
        vector<string> filename(1, "vg-file");
        return filename;
    });
    registry.register_recipe("XG", {"GFA"},
                             [&] (const vector<const IndexFile*>& inputs) {
        vector<string> filename(1, "xg-file");
        return filename;
    });
    registry.register_recipe("XG", {"VG"},
                             [&] (const vector<const IndexFile*>& inputs) {
        vector<string> filename(1, "xg-file");
        return filename;
    });
    registry.register_recipe("Pruned VG", {"VG"},
                             [&] (const vector<const IndexFile*>& inputs) {
        vector<string> filenames{"pruned-vg-file"};
        return filenames;
    });
    registry.register_recipe("GCSA+LCP", {"Pruned VG"},
                             [&] (const vector<const IndexFile*>& inputs) {
        vector<string> filenames{"gcsa-file", "lcp-file"};
        return filenames;
    });
    
    
    SECTION("Impossible and possible plans can be identified") {
        
        bool caught = false;
        try {
            registry.make_plan({"XG"});
        }
        catch (InsufficientInputException ex) {
            caught = true;
        }
        REQUIRE(caught);
        registry.provide("VCF", "vcf-name");
        caught = false;
        try {
            registry.make_plan({"XG"});
        }
        catch (InsufficientInputException ex) {
            caught = true;
        }
        REQUIRE(caught);
        
        registry.provide("FASTA", "fasta-name");
        // we now should have sufficient input to make this
        auto plan = registry.make_plan({"XG"});
        
        REQUIRE(plan.size() == 2);
        REQUIRE(plan[0].first == "VG");
        REQUIRE(plan[0].second == 0);
        REQUIRE(plan[1].first == "XG");
        REQUIRE(plan[1].second == 1);
    }
    
    SECTION("Plans can select preferred recipes") {
        
        registry.provide("VCF", "vcf-name");
        registry.provide("FASTA", "fasta-name");
        registry.provide("GFA", "gfa-name");
        
        auto plan = registry.make_plan({"XG"});
        REQUIRE(plan.size() == 1);
        REQUIRE(plan[0].first == "XG");
        REQUIRE(plan[0].second == 0);
    }
    
    SECTION("Plans can be made for multiple indexes") {
        
        registry.provide("VCF", "vcf-name");
        registry.provide("FASTA", "fasta-name");
        
        auto plan = registry.make_plan({"XG", "GCSA+LCP"});
        REQUIRE(plan.size() == 4);
        
        // TODO: this is ugly, i should have just used a map
        vector<int> which_item(plan.size(), -1);
        for (int i = 0; i < plan.size(); ++i) {
            if (plan[i].first == "XG") {
                which_item[0] = i;
            }
            else if (plan[i].first == "GCSA+LCP") {
                which_item[1] = i;
            }
            else if (plan[i].first == "Pruned VG") {
                which_item[2] = i;
            }
            else if (plan[i].first == "VG") {
                which_item[3] = i;
            }
        }
        
        // did we find them all?
        for (int i = 0; i < which_item.size(); ++i) {
            REQUIRE(which_item[i] != -1);
        }
        // are they in a feasible order?
        REQUIRE(which_item[0] > which_item[3]);
        REQUIRE(which_item[2] > which_item[3]);
        REQUIRE(which_item[1] > which_item[2]);
        // are they the recipes we expect?
        REQUIRE(plan[which_item[0]].second == 1);
        REQUIRE(plan[which_item[1]].second == 0);
        REQUIRE(plan[which_item[2]].second == 0);
        REQUIRE(plan[which_item[3]].second == 0);
    }
    
    SECTION("Midpoints of a pipeline can be provided directly") {
        
        registry.provide("VG", "vg-name");
        
        auto plan = registry.make_plan({"XG"});
        REQUIRE(plan.size() == 1);
        REQUIRE(plan[1].first == "XG");
        REQUIRE(plan[1].second == 1);
    }
}

}
}
