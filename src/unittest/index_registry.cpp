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
    
TEST_CASE("IndexRegistry can make plans on a dummy recipe graph", "[indexregistry]") {
    
    TestIndexRegistry registry;
    
    // name the indexes
    registry.register_index({"FASTA"}, "fasta");
    registry.register_index({"VCF"}, "vcf");
    registry.register_index({"GFA"}, "gfa");
    registry.register_index({"VG"}, "vg");
    registry.register_index({"Pruned VG"}, "pruned.vg");
    registry.register_index({"XG"}, "xg");
    registry.register_index({"GCSA", "LCP"}, "gcsa_lcp");
    registry.register_index({"Trivial Snarls"}, "snarls");
    registry.register_index({"Distance"}, "dist");
    
    // make some dummy recipes that don't actually do anything
    registry.register_recipe({"VG"}, {{"FASTA"}, {"VCF"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filename(1, "vg-file");
        return filename;
    });
    registry.register_recipe({"VG"}, {{"GFA"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filename(1, "vg-file");
        return filename;
    });
    registry.register_recipe({"XG"}, {{"GFA"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filename(1, "xg-file");
        return filename;
    });
    registry.register_recipe({"XG"}, {{"VG"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filename(1, "xg-file");
        return filename;
    });
    registry.register_recipe({"Pruned VG"}, {{"VG"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filenames{"pruned-vg-file"};
        return filenames;
    });
    registry.register_recipe({"GCSA", "LCP"}, {{"Pruned VG"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filenames{"gcsa-file", "lcp-file"};
        return filenames;
    });
    registry.register_recipe({"Trivial Snarls"}, {{"XG"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filenames(1, "snarls-file");
        return filenames;
    });
    registry.register_recipe({"Trivial Snarls"}, {{"VG"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filenames(1, "snarls-file");
        return filenames;
    });
    registry.register_recipe({"Distance"}, {{"XG"}, {"Trivial Snarls"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filenames(1, "dist-file");
        return filenames;
    });
    registry.register_recipe({"Distance"}, {{"VG"}, {"Trivial Snarls"}},
                             [&] (const vector<const IndexFile*>& inputs,
                                  const IndexingPlan* plan,
                                  const IndexName& constructing) {
        vector<string> filenames(1, "dist-file");
        return filenames;
    });
    
    
    SECTION("Impossible and possible plans can be identified") {
        
        bool caught = false;
        try {
            registry.make_plan({{"XG"}});
        }
        catch (InsufficientInputException ex) {
            caught = true;
        }
        REQUIRE(caught);
        registry.provide({"VCF"}, "vcf-name");
        caught = false;
        try {
            registry.make_plan({{"XG"}});
        }
        catch (InsufficientInputException ex) {
            caught = true;
        }
        REQUIRE(caught);
        
        registry.provide({"FASTA"}, "fasta-name");
        // we now should have sufficient input to make this
        auto plan = registry.make_plan({{"XG"}});
        
        REQUIRE(plan.steps.size() == 2);
        REQUIRE(plan.steps[0].first.size() == 1);
        REQUIRE(plan.steps[0].first.count("VG"));
        REQUIRE(plan.steps[0].second == 0);
        REQUIRE(plan.steps[1].first.size() == 1);
        REQUIRE(plan.steps[1].first.count("XG"));
        REQUIRE(plan.steps[1].second == 1);
    }
    
    SECTION("Plans can select preferred recipes") {
        
        registry.provide({"VCF"}, "vcf-name");
        registry.provide({"FASTA"}, "fasta-name");
        registry.provide({"GFA"}, "gfa-name");
        
        auto plan = registry.make_plan({{"XG"}});
        REQUIRE(plan.steps.size() == 1);
        REQUIRE(plan.steps[0].first.size() == 1);
        REQUIRE(plan.steps[0].first.count("XG"));
        REQUIRE(plan.steps[0].second == 0);
    }
    
    SECTION("Plans can be made for multiple indexes") {
        
        registry.provide({"VCF"}, "vcf-name");
        registry.provide({"FASTA"}, "fasta-name");
        
        auto plan = registry.make_plan({{"XG"}, {"GCSA", "LCP"}});
        REQUIRE(plan.steps.size() == 4);
        
        // Work out when everything is made 
        map<IndexName, size_t> made_at_step;
        for (size_t i = 0; i < plan.steps.size(); ++i) {
            made_at_step[plan.steps[i].first] = i;
        }
        
        // Make sure we make what we want
        REQUIRE(made_at_step.count({"VG"}));
        REQUIRE(made_at_step.count({"XG"}));
        REQUIRE(made_at_step.count({"Pruned VG"}));
        REQUIRE(made_at_step.count({"GCSA", "LCP"}));
        
        // are they in a feasible order?
        REQUIRE(made_at_step.at({"XG"}) > made_at_step.at({"VG"}));
        REQUIRE(made_at_step.at({"Pruned VG"}) > made_at_step.at({"VG"}));
        REQUIRE(made_at_step.at({"GCSA", "LCP"}) > made_at_step.at({"Pruned VG"}));
        
        // are they the recipes we expect?
        REQUIRE(plan.steps[made_at_step.at({"VG"})].second == 0);
        REQUIRE(plan.steps[made_at_step.at({"XG"})].second == 1);
        REQUIRE(plan.steps[made_at_step.at({"Pruned VG"})].second == 0);
        REQUIRE(plan.steps[made_at_step.at({"GCSA", "LCP"})].second == 0);
        
    }
    
    SECTION("Midpoints of a pipeline can be provided directly") {
        
        registry.provide({"VG"}, "vg-name");
        
        auto plan = registry.make_plan({{"XG"}});
        REQUIRE(plan.steps.size() == 1);
        REQUIRE(plan.steps[0].first.count("XG"));
        REQUIRE(plan.steps[0].second == 1);
    }
    
    SECTION("Impossible plans with some inputs available can be identified") {
        
        registry.provide({"Trivial Snarls"}, "snarls-name");
        
        bool caught = false;
        try {
            auto plan = registry.make_plan({{"Distance"}});
        }
        catch (InsufficientInputException ex) {
            caught = true;
        }
        REQUIRE(caught);
    }
    
    // Now add a fake simplification
    registry.register_joint_recipe({{"VG"}, {"XG"}}, {{"FASTA"}, {"VCF"}},
                                   [&] (const vector<const IndexFile*>& inputs,
                                        const IndexingPlan* plan) {
    
        vector<string> vg_filename(1, "vg-file");
        vector<string> xg_filename(1, "xg-file");
        return vector<vector<string>>{vg_filename, xg_filename};
    });
    
    SECTION("Plans are simplified when appropriate") {
        
        registry.provide({"VCF"}, "vcf-name");
        registry.provide({"FASTA"}, "fasta-name");
        
        auto plan = registry.make_plan({{"XG"}, {"GCSA", "LCP"}});
        REQUIRE(plan.steps.size() == 4);
        
        // Work out when everything is made 
        map<IndexName, size_t> made_at_step;
        for (size_t i = 0; i < plan.steps.size(); ++i) {
            made_at_step[plan.steps[i].first] = i;
        }
        
        // Make sure we make what we want
        REQUIRE(made_at_step.count({"VG"}));
        REQUIRE(made_at_step.count({"XG"}));
        REQUIRE(made_at_step.count({"Pruned VG"}));
        REQUIRE(made_at_step.count({"GCSA", "LCP"}));
        
        // are they in a feasible order?
        REQUIRE(made_at_step.at({"XG"}) > made_at_step.at({"VG"}));
        REQUIRE(made_at_step.at({"Pruned VG"}) > made_at_step.at({"VG"}));
        REQUIRE(made_at_step.at({"GCSA", "LCP"}) > made_at_step.at({"Pruned VG"}));
        
        // are they the recipes we expect?
        REQUIRE(plan.steps[made_at_step.at({"VG"})].second == 2);
        REQUIRE(plan.steps[made_at_step.at({"XG"})].second == 2);
        REQUIRE(plan.steps[made_at_step.at({"Pruned VG"})].second == 0);
        REQUIRE(plan.steps[made_at_step.at({"GCSA", "LCP"})].second == 0);
    }
}

}
}
