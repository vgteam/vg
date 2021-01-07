#ifndef VG_INDEX_REGISTRY_HPP_INCLUDED
#define VG_INDEX_REGISTRY_HPP_INCLUDED

#include <vector>
#include <unordered_map>
#include <functional>
#include <string>
#include <memory>

namespace vg {

using namespace std;

// Forward declarations
class IndexFile;
class IndexRecipe;
class InsufficientInputException;

/**
 * An object that can record methods to produce indexes and design
 * workflows to create a set of desired indexes
 */
class IndexRegistry {
public:
    /// Constructor
    IndexRegistry() = default;
    
    /// Register an index with the given identifier
    void register_index(const string& identifier);
    /// Register a recipe to produce an index using other indexes
    /// or input files
    void register_recipe(const string& identifier,
                         const vector<string>& input_identifiers,
                         const function<vector<string>(const vector<const IndexFile*>&)>& exec);
    /// Indicate a serialized file that contains some identified index
    void provide(const string& identifier, const string& filename);
    /// Indicate a list of serialized files that contains some identified index
    void provide(const string& identifier, const vector<string>& filenames);
    /// Get a list of all indexes that have already been completed or provided
    vector<string> completed_indexes() const;
    /// Create a plan to make the indicated indexes using provided inputs
    /// and execute the plan.
    /// If provided inputs cannot create the desired indexes, throws a
    /// InsufficientInputException.
    void make_indexes(const vector<string>& identifiers);
    
protected:
    
    /// get a topological ordering of all registered indexes in the dependency DAG
    vector<string> dependency_order() const;
    
    /// generate a plan to create the indexes
    vector<pair<string, size_t>> make_plan(const vector<string>& end_products) const;
    
    /// access index file
    IndexFile* get_index(const string& identifier);
    
    /// access const index file
    const IndexFile* get_index(const string& identifier) const;
    
    /// the storage struct for named indexes
    unordered_map<string, unique_ptr<IndexFile>> registry;
    
};

/**
 * An object that generically represents a serializable index or input file
 */
class IndexFile {
public:
    
    /// Create a new IndexFile with a unique identifier
    IndexFile(const string& identifier);
        
    /// Get the globally unique identifier for this index
    const string& get_identifier() const;
    
    /// Get the filename(s) that contain this index
    const vector<string>& get_filenames() const;
    
    /// Get the list of recipes that this index can be built with
    const vector<IndexRecipe>& get_recipes() const;
    
    /// Identify a serialized file that already contains this index
    void provide(const vector<string>& filenames);
    
    /// Describe a recipe with a lower preference order than all existing recipes
    /// for creating this index, if there are any (i.e., recipes must be added in
    /// preference order). Recipes should return the filepath(s) to their output.
    void add_recipe(const vector<const IndexFile*>& inputs,
                    const function<vector<string>(const vector<const IndexFile*>&)>& exec);
    
    /// Returns true if the index has already been built or provided
    bool is_finished() const;
    
    /// Build the index using the recipe with the provided priority
    void execute_recipe(size_t recipe_priority);
    
private:
    
    // the global identifier for the
    string identifier;
    
    // the filename(s) associated with the index
    vector<string> filenames;
    
    // the priority-ordered recipes to make this index file
    vector<IndexRecipe> recipes;
};

/**
 * struct that indicates a method to produce and serialize an index
 */
struct IndexRecipe {
    IndexRecipe(const vector<const IndexFile*>& inputs,
                const function<vector<string>(const vector<const IndexFile*>&)>& exec);
    // execute the recipe and return the filename(s) of the indexes created
    vector<string> execute();
    vector<const IndexFile*> inputs;
    function<vector<string>(const vector<const IndexFile*>&)> exec;
};


/**
 * Exception that is thrown to indicate the input data is insufficient
 * to create some index(es)
 */
class InsufficientInputException : public runtime_error {
public:
    InsufficientInputException() = delete;
    InsufficientInputException(const string& target,
                               const IndexRegistry& registry);
    const char* what() const throw ();
private:
    string target;
    vector<string> inputs;
};

}

#endif
