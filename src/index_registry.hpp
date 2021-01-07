#ifndef VG_INDEX_REGISTRY_HPP_INCLUDED
#define VG_INDEX_REGISTRY_HPP_INCLUDED

#include <vector>
#include <unordered_map>
#include <functional>
#include <string>

namespace vg {

using namespace std;

class IndexFile;
class IndexRecipe;
class InsufficientInputException;

class IndexRegistry {
public:
    IndexRegistry() = default;
    void register_index(const string& identifier);
    void register_recipe(const string& identifier,
                         const vector<string>& input_identifiers,
                         const function<vector<string>(const vector<const IndexFile*>&)>& exec);
    void provide(const string& identifier, const string& filename);
    void provide(const string& identifier, const vector<string>& filenames);
    void make_indexes(const vector<string>& identifiers);
    
    vector<string> provided_indexes() const;
    
protected:
    
    /// get a topological ordering of all registered indexes in the dependency DAG
    vector<string> dependency_order() const;
    
    /// generate a plan
    vector<pair<string, size_t>> make_plan(const vector<string>& end_products) const;
    
    IndexFile* get_index(const string& identifier);
    
    const IndexFile* get_index(const string& identifier) const;
    
    unordered_map<string, unique_ptr<IndexFile>> registry;
    
};

class IndexFile {
public:
    
    /**
     * Create a new IndexFile with a unique identifier
     */
    IndexFile(const string& identifier);
        
    /**
     * Get the identifier.
     */
    const string& get_identifier() const;
    
    const vector<string>& get_filenames() const;
    
    const vector<IndexRecipe>& get_recipes() const;
    
    void provide(const vector<string>& filenames);
    
    /// Describe a recipe with a lower preference order than all existing recipes
    /// for creating this index, if there are any (i.e., recipes must be added in
    /// preference order). Recipes should return the filepath(s) to their output.
    void add_recipe(const vector<const IndexFile*>& inputs,
                    const function<vector<string>(const vector<const IndexFile*>&)>& exec);
    
    ///
    bool is_finished() const;
    
    void execute_recipe(size_t recipe_priority);
    
private:
    
    // the global identifier for the
    string identifier;
    
    // the filename(s) associated with the index
    vector<string> filenames;
    
    // the priority-ordered recipes to make this index file
    vector<IndexRecipe> recipes;
};

struct IndexRecipe {
    IndexRecipe(const vector<const IndexFile*>& inputs,
                const function<vector<string>(const vector<const IndexFile*>&)>& exec);
    // execute the recipe and return the filename(s) of the indexes created
    vector<string> execute();
    vector<const IndexFile*> inputs;
    function<vector<string>(const vector<const IndexFile*>&)> exec;
};


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
