#ifndef VG_INDEX_REGISTRY_HPP_INCLUDED
#define VG_INDEX_REGISTRY_HPP_INCLUDED

#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <functional>
#include <string>
#include <memory>
#include <stdexcept>

namespace vg {

using namespace std;

// Forward declarations
class IndexFile;
class IndexRecipe;
class InsufficientInputException;
class IndexRegistry;
class IndexingPlan;

/**
 * Names an index. Each index may have multiple components.
 */
using IndexName = set<string>;
// We use a sorted set so that we can know which item in the name corresponds to which file

/**
 * Names a recipe in the collection of registered recipes.
 */
using RecipeName = pair<IndexName, size_t>;

/**
 * Is a recipe to create the files (returned by name) associated with some
 * index, from a series of input indexes, given the plan it is being
 * generated for and the index being generated.
 */
using RecipeFunc = function<vector<string>(const vector<const IndexFile*>&,const IndexingPlan*,const IndexName&)>;

/**
 * Is a recipe to create the files (returned by name) associated with some
 * indexes, from a series of input indexes, given the plan they are being
 * generated for.
 */
using JointRecipeFunc = function<vector<vector<string>>(const vector<const IndexFile*>&,const IndexingPlan*)>;

/**
 * A struct namespace for global handling of parameters used by
 * the IndexRegistry
 */
struct IndexingParameters {
    // enums for categorical options
    enum MutableGraphImplementation {HashGraph, PackedGraph, ODGI, VG};
    
    // the actual parameters
    
    // the format that "VG" indexes will be created in [HashGraph]
    static MutableGraphImplementation mut_graph_impl;
    // the maximum node length for graphs that are created from VCFs [32]
    static int max_node_size;
    // during pruning, remove nodes with degree higher than this [128] // TODO: is this a good default?
    static int pruning_max_node_degree;
    // during pruning, identify complex regions using walks of this length [24]
    static int pruning_walk_length;
    // during pruning, remove edges if a walk contains this many branching edges [3]
    static int pruning_max_edge_count;
    // during pruning, remove any isolated components with at less than this total seq length [33]
    static int pruning_min_component_size;
    // length of the k-mers indexed in GCSA2 before any doubling steps [16]
    static int gcsa_initial_kmer_length;
    // number of k-mer length doubling steps in GCSA2 [4]
    static int gcsa_doubling_steps;
    // number of gbwt nodes inserted at a time in dynamic gbwt [100M]
    static int gbwt_insert_batch_size;
    // the sampling interval in the GBWT suffix array [1024]
    static int gbwt_sampling_interval;
    // should the haplotype-transcript GBWT be made bidirectional [false]
    static bool bidirectional_haplo_tx_gbwt;
    // whether indexing algorithms will log progress (if available) [false]
    static bool verbose;
};

/**
 * A struct namespace for standard inputs
 */
struct VGIndexes {
    /// A complete index registry for VG mapping utilities
    static IndexRegistry get_vg_index_registry();
    /// A list of the identifiers of the default indexes to run vg map
    static vector<IndexName> get_default_map_indexes();
    /// A list of the identifiers of the default indexes to run vg mpmap
    static vector<IndexName> get_default_mpmap_indexes();
    /// A list of the identifiers of the default indexes to run vg giraffe
    static vector<IndexName> get_default_giraffe_indexes();
};

/**
 * A plan for producing indexes, which knows what should be saved and what should be ephemeral.
 * Wants to be nested inside IndexRegistry, but you can't forward-declare a nested class.
 */
class IndexingPlan {

public:
    // The IndexRegistry is responsible for setting us up.
    friend class IndexRegistry;

    /// Returns true if the given index is to be intermediate under the given
    /// plan, and false if it is to be preserved.
    bool is_intermediate(const IndexName& identifier) const;
    
    /// Get the prefix with which to save the given index's files, accounting
    /// for whether it is intermediate or not.
    string prefix(const IndexName& identifier) const;
    
    /// Get the suffix with which to save the given index's files.
    string suffix(const IndexName& identifier) const;
    
    /// The steps to be invoked in the plan. May be empty before the plan is
    /// actually planned.
    vector<RecipeName> steps;
    
private:
    /// The indexes to create as outputs.
    set<IndexName> targets;
    
    /// The registry that the plan is using.
    /// The registry must not move while the plan is in use.
    /// Can't be const because we need to get_work_dir() on it, which may
    /// create the work directory.
    IndexRegistry* registry;
};

/**
 * An object that can record methods to produce indexes and design
 * workflows to create a set of desired indexes
 */
class IndexRegistry {
public:

    // IndexingPlan can't be a child class, but it needs to be able to
    // get_index, so it has to be a friend.
    friend class IndexingPlan;

    /// Constructor
    IndexRegistry() = default;
    
    /// Destructor to clean up temp files.
    ~IndexRegistry();
    
    // Because we own temporary files and unique pointers, we should not be copied.
    IndexRegistry(const IndexRegistry& other) = delete;
    IndexRegistry& operator=(const IndexRegistry& other) = delete;
    
    // And we need to be moved carefully
    IndexRegistry(IndexRegistry&& other);
    IndexRegistry& operator=(IndexRegistry&& other);
    
    
    /// Prefix for all saved outputs
    void set_prefix(const string& prefix);
    
    /// Should intermediate files be saved to the output directory
    /// or the temp directory?
    void set_intermediate_file_keeping(bool keep_intermediates);
    
    /// Register an index containing the given identifier
    void register_index(const IndexName& identifier, const string& suffix);
    
    /// Register a recipe to produce an index using other indexes
    /// or input files. Also takes a for output as input
    RecipeName register_recipe(const IndexName& identifier,
                               const vector<IndexName>& input_identifiers,
                               const RecipeFunc& exec);
                        
    /// Register a recipe to produce multiple indexes.
    /// Individual index recipes must still be registered; this recipe will be
    /// used to simplify the plan when the individual recipes with the same
    /// input set are all called.
    ///
    /// Joint recipes may also be used to generate single indexes if they have
    /// higher priority than other recipes.
    ///
    /// All output identifiers must be distinct, and there must be at least
    /// one. Only one multi-index recipe can be applied to simplify the
    /// production of any given index in a plan.
    void register_joint_recipe(const vector<IndexName>& identifiers,
                               const vector<IndexName>& input_identifiers,
                               const JointRecipeFunc& exec);
    
    /// Indicate a serialized file that contains some identified index
    void provide(const IndexName& identifier, const string& filename);
    
    /// Indicate a list of serialized files that contains some identified index
    void provide(const IndexName& identifier, const vector<string>& filenames);
    
    /// Get a list of all indexes that have already been completed or provided
    vector<IndexName> completed_indexes() const;
    
    /// Create and execute a plan to make the indicated indexes using provided inputs
    /// If provided inputs cannot create the desired indexes, throws a
    /// InsufficientInputException.
    void make_indexes(const vector<IndexName>& identifiers);
    
    /// Returns the recipe graph in dot format
    string to_dot() const;
    
    /// Returns the recipe graph in dot format with a plan highlighted
    string to_dot(const vector<IndexName>& targets) const;
    
protected:
    
    /// get a topological ordering of all registered indexes in the dependency DAG
    vector<IndexName> dependency_order() const;
    
    /// generate a plan to create the indexes
    IndexingPlan make_plan(const vector<IndexName>& end_products) const;
    
    /// access index file
    IndexFile* get_index(const IndexName& identifier);
    
    /// access const index file
    const IndexFile* get_index(const IndexName& identifier) const;
    
    /// Function to get and/or initialize the temporary directory in which indexes will live
    string get_work_dir();
    
    /// The storage struct for named indexes. Ordered so it is easier to key on index names.
    map<IndexName, unique_ptr<IndexFile>> registry;
    
    unordered_set<string> registered_suffixes;
    
    /// Record that, when the given input indexes are available, this
    /// collection of recipes is efficient to run together.
    vector<pair<vector<IndexName>, vector<RecipeName>>> simplifications;
    
    /// Temporary directory in which indexes will live
    string work_dir;
    
    /// filepath that will prefix all saved output
    string output_prefix = "index";
    
    /// should intermediate files end up in the scratch or the output directory?
    bool keep_intermediates = false;
};

/**
 * An object that generically represents a serializable index or input file
 */
class IndexFile {
public:
    
    /// Create a new IndexFile with a unique identifier
    IndexFile(const IndexName& identifier, const string& suffix);
        
    /// Get the globally unique identifier for this index
    const IndexName& get_identifier() const;
    
    /// Returns the suffix to be used for this index
    const string& get_suffix() const;
    
    /// Get the filename(s) that contain this index
    const vector<string>& get_filenames() const;
    
    /// Get the list of recipes that this index can be built with
    const vector<IndexRecipe>& get_recipes() const;
    
    /// Identify a serialized file that already contains this index
    void provide(const vector<string>& filenames);
    
    /// Describe a recipe with a lower preference order than all existing recipes
    /// for creating this index, if there are any (i.e., recipes must be added in
    /// preference order). Recipes should return the filepath(s) to their output.
    ///
    /// Returns the name assigned to the recipe.
    RecipeName add_recipe(const vector<const IndexFile*>& inputs,
                          const RecipeFunc& exec);
    
    /// Returns true if the index has already been built or provided
    bool is_finished() const;
    
    /// Build the index using the recipe with the provided priority.
    /// Expose the plan so that the recipe knows where it is supposed to go.
    void execute_recipe(size_t recipe_priority, const IndexingPlan* plan);
    
    /// Returns true if the index was provided through provide method
    bool was_provided_directly() const;
    
private:
    
    // the global identifier for the
    IndexName identifier;
    
    // the suffix it adds to output files
    string suffix;
    
    // the filename(s) associated with the index
    vector<string> filenames;
    
    // the priority-ordered recipes to make this index file
    vector<IndexRecipe> recipes;
    
    // keep track of whether the index was provided directly
    bool provided_directly = false;
};

/**
 * struct that indicates a method to produce and serialize an index
 */
struct IndexRecipe {
    IndexRecipe(const vector<const IndexFile*>& inputs,
                const RecipeFunc& exec);
    // execute the recipe and return the filename(s) of the indexes created
    vector<string> execute(const IndexingPlan* plan, const IndexName& constructing);
    vector<const IndexFile*> inputs;
    RecipeFunc exec;
};


/**
 * Exception that is thrown to indicate the input data is insufficient
 * to create some index(es)
 */
class InsufficientInputException : public runtime_error {
public:
    InsufficientInputException() = delete;
    InsufficientInputException(const IndexName& target,
                               const IndexRegistry& registry);
    const char* what() const throw ();
private:
    IndexName target;
    vector<IndexName> inputs;
};

}

#endif
