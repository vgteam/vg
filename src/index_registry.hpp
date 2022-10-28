#ifndef VG_INDEX_REGISTRY_HPP_INCLUDED
#define VG_INDEX_REGISTRY_HPP_INCLUDED

#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <string>
#include <memory>
#include <stdexcept>
#include <limits>

namespace vg {

using namespace std;

// Forward declarations
class IndexFile;
class IndexRecipe;
class InsufficientInputException;
class IndexRegistry;
class IndexingPlan;
class AliasGraph;

/**
 * A unique identifier for an Index
 */
using IndexName = string;

/**
 * A group of indexes that can be made simultaneously
 */
using IndexGroup = set<IndexName>;

/**
 * Names a recipe in the collection of registered recipes.
 */
using RecipeName = pair<IndexGroup, size_t>;

/**
 * Is a recipe to create the files (returned by name) associated with some
 * index, from a series of input indexes, given the plan it is being
 * generated for and the index being generated.
 */
using RecipeFunc = function<vector<vector<string>>(const vector<const IndexFile*>&,
                                                   const IndexingPlan*,
                                                   AliasGraph&,
                                                   const IndexGroup&)>;

/**
 * Is a recipe to create the files (returned by name) associated with some
 * indexes, from a series of input indexes, given the plan they are being
 * generated for.
 */
using JointRecipeFunc = function<vector<vector<vector<string>>>(const vector<const IndexFile*>&,const IndexingPlan*)>;

/**
 * A struct namespace for global handling of parameters used by
 * the IndexRegistry
 */
struct IndexingParameters {
    // enums for categorical options
    enum MutableGraphImplementation {HashGraph, PackedGraph, ODGI, VG};
    enum Verbosity {None = 0, Basic = 1, Debug = 2};
    
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
    // the feature from column 3 of the GTF/GFF that should be added ["gene"]
    static string gff_feature_name;
    // transcript tag in GTF/GFF ["transcript_id"]
    static string gff_transcript_tag;
    // if true, minimizer index uses bounded syncmers, otherwise uses minimizers [false]
    static bool use_bounded_syncmers;
    // length of k-mer used in minimizer index [29]
    static int minimizer_k;
    // length of window if using minimizers [11]
    static int minimizer_w;
    // length of internal s-mer if using bounded syncmers [18]
    static int minimizer_s;
    // the number of paths that will make up the path cover GBWT [16]
    static int path_cover_depth;
    // the number of haplotypes to downsample to in giraffe's GBWT [64]
    static int giraffe_gbwt_downsample;
    // sample subpaths of this length (in nodes) [4]
    static int downsample_context_length;
    // augment the existing GBWT instead of downsampling it if the number of haplotypes is < this * giraffe_gbwt_downsample [3]
    static int downsample_threshold;
    // actually use this fraction of the maximum memory to give slosh for bad estmates [0.75]
    static double max_memory_proportion;
    // aim to have X timese as many chunks as threads [2]
    static double thread_chunk_inflation_factor;
    // whether indexing algorithms will log progress (if available) [Basic]
    static Verbosity verbosity;
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
    
    IndexingPlan() = default;
    ~IndexingPlan() = default;
    
    /// Get the suffix with which to save the given index's files.
    string output_filepath(const IndexName& identifier) const;
    
    /// Get the suffix with which to save the given index's files.
    string output_filepath(const IndexName& identifier, size_t chunk, size_t num_chunks) const;
    
    /// Ge the steps of the plan
    const vector<RecipeName>& get_steps() const;
    
    /// Returns true if the given index is to be intermediate under the given
    /// plan, and false if it is to be preserved.
    bool is_intermediate(const IndexName& identifier) const;
    
    /// TODO: is this where this function wants to live?
    int64_t target_memory_usage() const;
    
    /// Returns the recipes in the plan that depend on this index, including the one in which
    /// it was created (if any)
    set<RecipeName> dependents(const IndexName& identifier) const;
    
protected:
    
    /// The steps to be invoked in the plan. May be empty before the plan is
    /// actually planned.
    vector<RecipeName> steps;
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
 * workflows to create a set of desired indexes.
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
    
    /// Get the current prefix for saving output files.
    string get_prefix() const;
    
    /// Should intermediate files be saved to the output directory
    /// or the temp directory?
    void set_intermediate_file_keeping(bool keep_intermediates);
    
    /// Register an index containing the given identifier
    void register_index(const IndexName& identifier, const string& suffix);
    
    /// Register a recipe to produce an index using other indexes
    /// or input files. Recipes registered earlier will have higher priority.
    RecipeName register_recipe(const vector<IndexName>& identifiers,
                               const vector<IndexName>& input_identifiers,
                               const RecipeFunc& exec);
                        
//    /// Register a recipe to produce multiple indexes.
//    /// Individual index recipes must still be registered; this recipe will be
//    /// used to simplify the plan when the individual recipes with the same
//    /// input set are all called.
//    ///
//    /// Joint recipes may also be used to generate single indexes if they have
//    /// higher priority than other recipes.
//    ///
//    /// All output identifiers must be distinct, and there must be at least
//    /// one. Only one multi-index recipe can be applied to simplify the
//    /// production of any given index in a plan.
//    void register_joint_recipe(const vector<IndexName>& identifiers,
//                               const vector<IndexName>& input_identifiers,
//                               const JointRecipeFunc& exec);
    
    /// Indicate a serialized file that contains some identified index
    void provide(const IndexName& identifier, const string& filename);
    
    /// Indicate a list of serialized files that contains some identified index
    void provide(const IndexName& identifier, const vector<string>& filenames);
    
    /// Return true if the given index is available and can be require()'d, and
    /// false otherwise.
    bool available(const IndexName& identifier) const;
    
    /// Get the filename(s) associated with the given index. Aborts if the
    /// index is not a known type, or if it is not provided or made.
    vector<string> require(const IndexName& identifier) const;
    
    /// Set the maximum memory that indexing should try to consume (note: this is
    /// not strictly adhered to due to difficulties in estimating memory use)
    void set_target_memory_usage(int64_t bytes);
    
    /// Get the maximum memory we will try to consume
    int64_t get_target_memory_usage() const;
    
    /// Get the amount of free memory
    static int64_t get_system_memory();
    
    /// Get a list of all indexes that have already been completed or provided
    vector<IndexName> completed_indexes() const;
    
    /// Create and execute a plan to make the indicated indexes using provided inputs
    /// If provided inputs cannot create the desired indexes, throws a
    /// InsufficientInputException.
    /// When completed, all requested index files will be available via require().
    void make_indexes(const vector<IndexName>& identifiers);
    
    /// Returns the recipe graph in dot format
    string to_dot() const;
    
    /// Returns the recipe graph in dot format with a plan highlighted
    string to_dot(const vector<IndexName>& targets) const;
    
    /// Determine if a VCF file is phased or not
    static bool vcf_is_phased(const string& filepath);
    
    /// Determine if a GFA has haplotypes as W-lines
    static bool gfa_has_haplotypes(const string& filepath);
    
    /// Discard any provided or constructed indexes
    void reset();
    
protected:
    
    /// get a topological ordering of all registered indexes in the dependency DAG
    vector<IndexGroup> dependency_order() const;
    
    /// generate a plan to create the indexes
    IndexingPlan make_plan(const IndexGroup& end_products) const;
    
    /// use a recipe identifier to get the recipe
    const IndexRecipe& get_recipe(const RecipeName& recipe_name) const;
    
    /// Build the index using the recipe with the provided priority.
    /// Expose the plan so that the recipe knows where it is supposed to go.
    vector<vector<string>> execute_recipe(const RecipeName& recipe_name, const IndexingPlan* plan,
                                          AliasGraph& alias_graph);
    
    /// access index file
    IndexFile* get_index(const IndexName& identifier);
    
    /// access const index file
    const IndexFile* get_index(const IndexName& identifier) const;
    
    bool all_finished(const vector<const IndexFile*>& inputs) const;
    
    bool all_finished(const IndexGroup& inputs) const;
    
    /// Function to get and/or initialize the temporary directory in which indexes will live
    string get_work_dir();
    
    /// The storage struct for named indexes. Ordered so it is easier to key on index names.
    map<IndexName, unique_ptr<IndexFile>> index_registry;
    
    /// All of the suffixes that have been registered by indexes
    unordered_set<string> registered_suffixes;
    
    /// The storage struct for recipes, which may make index
    map<IndexGroup, vector<IndexRecipe>> recipe_registry;
    
    /// Record that, when the given input indexes are available, this
    /// collection of recipes is efficient to run together.
    vector<pair<vector<IndexName>, vector<RecipeName>>> simplifications;
    
    /// Temporary directory in which indexes will live
    string work_dir;
    
    /// filepath that will prefix all saved output
    string output_prefix = "index";
    
    /// should intermediate files end up in the scratch or the output directory?
    bool keep_intermediates = false;
    
    /// the max memory we will *attempt* to use
    int64_t target_memory_usage = numeric_limits<int64_t>::max();
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
    
    /// Identify a serialized file that already contains this index
    void provide(const vector<string>& filenames);
    
    /// Assign constructed filenames to this index
    void assign_constructed(const vector<string>& filenames);
    
    /// Returns true if the index has already been built or provided
    bool is_finished() const;
    
    /// Returns true if the index was provided through provide method
    bool was_provided_directly() const;
    
    /// Discard any constructed or provided indexes
    void reset();
    
private:
    
    // the global identifier for the
    IndexName identifier;
    
    // the suffix it adds to output files
    const string suffix;
    
    // the filename(s) associated with the index
    vector<string> filenames;
    
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
    vector<vector<string>> execute(const IndexingPlan* plan, AliasGraph& alias_graph,
                                   const IndexGroup& constructing) const;
    IndexGroup input_group() const;
    vector<const IndexFile*> inputs;
    RecipeFunc exec;
};

/**
 * Class to keep track of which indexes are aliasing other indexes
 */
class AliasGraph {
public:
    AliasGraph() = default;
    ~AliasGraph() = default;
    
    /// Record that one index is aliasing another
    void register_alias(const IndexName& aliasor, const IndexFile* aliasee);
    
    /// Return a list of all indexes that are being aliased by non-intermediate
    /// indexes. If the aliasee is non-intermediate itself, it ill be listed among the
    /// aliases too.
    vector<pair<IndexName, vector<IndexName>>> non_intermediate_aliases(const IndexingPlan* plan,
                                                                        bool keep_all) const;
    
private:
    
    // graph aliasees to their aliasors
    unordered_map<IndexName, vector<IndexName>> graph;
    
};


/**
 * Exception that is thrown to indicate the input data is insufficient
 * to create some index(es)
 */
class InsufficientInputException : public runtime_error {
public:
    InsufficientInputException() = delete;
    InsufficientInputException(const IndexName& target,
                               const IndexRegistry& registry) noexcept;
    const char* what() const noexcept;
private:
    IndexName target;
    vector<IndexName> inputs;
};


/**
 * An exception that indicates that we must rewind the plan to re-create some indexes
 */
class RewindPlanException : public std::exception {
public:
    
    RewindPlanException() = delete;
    RewindPlanException(const string& msg, const IndexGroup& rewind_to) noexcept;
    ~RewindPlanException() noexcept = default;
    
    const char* what() const noexcept;
    const IndexGroup& get_indexes() const noexcept;
    
private:
    
    const string msg;
    IndexGroup indexes;
    
};

}

#endif
