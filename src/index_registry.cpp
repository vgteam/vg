// index_registry.cpp: index registry system implementation

#include "index_registry.hpp"

#include <unordered_set>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>

#include "hash_map.hpp"

//#define debug_index_registry

namespace vg {

void IndexRegistry::make_indexes(const vector<string>& identifiers) {
    auto plan = make_plan(identifiers);
    for (auto& step : plan) {
        get_index(step.first)->execute_recipe(step.second);
    }
}

void IndexRegistry::register_index(const string& identifier) {
    // Add this index to the registry
    if (registry.count(identifier)) {
        cerr << "error:[IndexRegistry] index registry contains a duplicated identifier: " << identifier << endl;
        exit(1);
    }
    registry[identifier] = unique_ptr<IndexFile>(new IndexFile(identifier));
}


void IndexRegistry::provide(const string& identifier, const string& filename) {
    provide(identifier, vector<string>(1, filename));
}

void IndexRegistry::provide(const string& identifier, const vector<string>& filenames) {
    get_index(identifier)->provide(filenames);
}

vector<string> IndexRegistry::provided_indexes() const {
    vector<string> indexes;
    for (const auto& index : registry) {
        if (index.second->is_finished()) {
            indexes.push_back(index.first);
        }
    }
    return indexes;
}

void IndexRegistry::register_recipe(const string& identifier,
                                    const vector<string>& input_identifiers,
                                    const function<vector<string>(const vector<const IndexFile*>&)>& exec) {
    vector<const IndexFile*> inputs;
    for (const auto& input_identifier : input_identifiers) {
        inputs.push_back(get_index(input_identifier));
    }
    get_index(identifier)->add_recipe(inputs, exec);
}

IndexFile* IndexRegistry::get_index(const string& identifier) {
    return registry.at(identifier).get();
}

const IndexFile* IndexRegistry::get_index(const string& identifier) const {
    return registry.at(identifier).get();
}

vector<string> IndexRegistry::dependency_order() const {
    
#ifdef debug_index_registry
    cerr << "finding topological order in dependency graph" << endl;
#endif
    
    // assign each index file an index in a vector (arbitrarily)
    unordered_map<string, size_t> graph_idx;
    vector<string> graph_label;
    for (const auto& idx_file : registry) {
        graph_idx[idx_file.first] = graph_label.size();
        graph_label.push_back(idx_file.first);
    }
    
    // build the dependency graph
    vector<vector<size_t>> dependency_graph(graph_label.size());
    for (size_t i = 0; i < dependency_graph.size(); ++i) {
        auto index = get_index(graph_label[i]);
        for (const auto& recipe : index->get_recipes()) {
            for (auto input : recipe.inputs) {
                dependency_graph[graph_idx[input->get_identifier()]].push_back(i);
            }
        }
    }
    
    // deduplicate any edges
    for (auto& adj : dependency_graph) {
        sort(adj.begin(), adj.end());
        adj.resize(unique(adj.begin(), adj.end()) - adj.begin());
    }
    
    // kahn's algorithm to determine a topological order
    vector<size_t> in_degree(dependency_graph.size(), 0);
    for (auto& adj : dependency_graph) {
        for (size_t i : adj) {
            ++in_degree[i];
        }
    }
    
    vector<size_t> stack;
    for (size_t i = 0; i < dependency_graph.size(); ++i) {
        if (in_degree[i] == 0) {
            stack.push_back(i);
        }
    }
    
    vector<size_t> order;
    while (!stack.empty()) {
        size_t i = stack.back();
        stack.pop_back();
        order.push_back(i);
        for (size_t j : dependency_graph[i]) {
            --in_degree[j];
            if (in_degree[j] == 0) {
                stack.push_back(j);
            }
        }
    }
    
    if (order.size() != dependency_graph.size()) {
        cerr << "error:[IndexFile] index dependency graph is not a DAG" << endl;
        exit(1);
    }
    
    // convert to return format
    vector<string> ordered_identifiers(order.size());
    for (size_t i = 0; i < order.size(); ++i) {
        ordered_identifiers[i] = graph_label[order[i]];
    }
    
#ifdef debug_index_registry
    for (const auto& identifier : ordered_identifiers) {
        cerr << "\t" << identifier << endl;
    }
#endif
    
    return ordered_identifiers;
}

vector<pair<string, size_t>> IndexRegistry::make_plan(const vector<string>& end_products) const {
    
#ifdef debug_index_registry
    cerr << "generating plan for indexes:" << endl;
    for (const auto& product : end_products) {
        cerr << "\t" << product << endl;
    }
#endif
    
    // get the dependency ordering of the indexes
    vector<string> identifier_order = dependency_order();
    unordered_map<string, size_t> dep_order_of_identifier;
    for (size_t i = 0; i < identifier_order.size(); ++i) {
        dep_order_of_identifier[identifier_order[i]] = i;
    }
    
    // TODO: I'm sure there's a more elegant implementation of this algorithm
    unordered_set<pair<string, size_t>> plan_elements;
    for (const auto& product : end_products) {
        
        // records of (identifier, lowest level requester, ordinal index of recipe selected)
        vector<tuple<size_t, size_t, size_t>> plan_path;
        
        // map dependency priority to lowest level priority that requested this and
        // the number of requesters
        map<size_t, pair<size_t, size_t>, greater<size_t>> queue;
        queue[dep_order_of_identifier[product]] = pair<size_t, size_t>(identifier_order.size(), 1);
        
        while (!queue.empty()) {
            
            // get the latest file in the dependency order
            // that we have left to build
            auto it = queue.begin();
            plan_path.emplace_back(it->first, it->second.first, 0);
            
#ifdef debug_index_registry
            cerr << "dequeue " << identifier_order[it->first] << " requested from " << (it->second.first == identifier_order.size() ? string("PLAN TARGET") : identifier_order[it->second.first]) << " and " << (it->second.second - 1) << " other indexes" << endl;
#endif
            queue.erase(it);
            
            if (get_index(identifier_order[get<0>(plan_path.back())])->is_finished()) {
                // this index has been provided, we don't need to use a recipe
#ifdef debug_index_registry
                cerr << "file has been provided as input" << endl;
#endif
                continue;
            }
            else if (!get_index(identifier_order[get<0>(plan_path.back())])->get_recipes().empty()) {
                
                // this index can be created by a recipe
                const auto& recipe = get_index(identifier_order[get<0>(plan_path.back())])->get_recipes().front();
#ifdef debug_index_registry
                cerr << "index can be made by a recipe requiring";
                for (auto input : recipe.inputs) {
                    cerr << " " << input->get_identifier();
                }
                cerr << endl;
#endif
                
                for (auto input : recipe.inputs) {
                    size_t dep_order = dep_order_of_identifier[input->get_identifier()];
                    auto f = queue.find(dep_order);
                    if (f == queue.end()) {
                        // no lower-level index has requested this one yet
                        queue[dep_order] = pair<size_t, size_t>(get<0>(plan_path.back()), 1);
                    }
                    else {
                        // record that one more index is requesting this one
                        f->second.second++;
                    }
                    
                }
            }
            else {
                // we've reached a file that needs to be provided but we don't have it,
                // so now we backtrack until hitting something that has remaining
                // lower priority recipes
#ifdef debug_index_registry
                cerr << "file cannot be made from existing inputs" << endl;
#endif
                while (!plan_path.empty() &&
                       get<2>(plan_path.back()) == get_index(identifier_order[get<0>(plan_path.back())])->get_recipes().size()) {
                    // there are no remaining recipes to build the last index in the plan
                    
                    // remove items off the plan path until we get to the index that requested
                    // this one
                    size_t requester = get<1>(plan_path.back());
#ifdef debug_index_registry
                    cerr << "pruning path to previous requester: " << identifier_order[requester] << endl;
#endif
                    while (!plan_path.empty() && get<0>(plan_path.back()) != requester) {
                        
                        auto index = get_index(identifier_order[get<0>(plan_path.back())]);
                        
                        if (!index->is_finished() && get<2>(plan_path.back()) < index->get_recipes().size()) {
                            // this index was using a recipe, we need to update its dependencies
                            // that are currently in the queue
                            const auto& recipe = index->get_recipes().at(get<2>(plan_path.back()));
                            for (auto input : recipe.inputs) {
                                size_t input_dep_order = dep_order_of_identifier[input->get_identifier()];
                                auto q = queue.find(input_dep_order);
                                if (q != queue.end()) {
                                    // there is now one fewer index requesting this index as input
                                    --q->second.second;
                                    if (q->second.second == 0) {
                                        // this is the only index that's requesting this queued index,
                                        // so we can remove it from the queue
                                        queue.erase(q);
                                    }
                                }
                            }
                        }
                        
                        plan_path.pop_back();
                    }
                    
                    if (!plan_path.empty()) {
                        // the requester should now use its next highest priority recipe
                        ++get<2>(plan_path.back());
                    }
                }
                
                if (!plan_path.empty()) {
                    const auto& recipe = get_index(identifier_order[get<0>(plan_path.back())])->get_recipes()[get<2>(plan_path.back())];
                    
#ifdef debug_index_registry
                    cerr << "advancing to recipe " << get<2>(plan_path.back()) << " for index " << identifier_order[get<0>(plan_path.back())] << ", which requires";
                    for (auto input : recipe.inputs) {
                        cerr << " " << input->get_identifier();
                    }
                    cerr << endl;
#endif
                    for (auto input : recipe.inputs) {
                        size_t dep_order = dep_order_of_identifier[input->get_identifier()];
                        auto f = queue.find(dep_order);
                        if (f == queue.end()) {
                            // no lower-level index has requested this one yet
                            queue[dep_order] = pair<size_t, size_t>(get<0>(plan_path.back()), 1);
                        }
                        else {
                            // record that one more index is requesting this one
                            f->second.second++;
                        }
                        
                    }
                }
            }
            
        }
        
#ifdef debug_index_registry
        cerr << "final plan path for index " << product << ":" << endl;
        for (auto path_elem : plan_path) {
            cerr << "\t" << identifier_order[get<0>(path_elem)] << ", from " << (get<1>(path_elem) == identifier_order.size() ? "PLAN START" : identifier_order[get<1>(path_elem)]) << ", recipe " << get<2>(path_elem) << endl;
        }
#endif
        
        if (plan_path.empty()) {
            // we don't have enough of the inputs to create this index
            throw InsufficientInputException(product, *this);
        }
        
        // record the elements of this plan
        for (size_t i = 0; i < plan_path.size(); ++i) {
            plan_elements.emplace(identifier_order[get<0>(plan_path[i])], get<2>(plan_path[i]));
        }
    }
    
    // convert the aggregated plan elements into a forward ordered plan
    vector<pair<string, size_t>> plan(plan_elements.begin(), plan_elements.end());
    sort(plan.begin(), plan.end(), [&](const pair<string, size_t>& a, const pair<string, size_t>& b) {
        return dep_order_of_identifier[a.first] < dep_order_of_identifier[b.first];
    });
#ifdef debug_index_registry
    cerr << "full plan including provided files:" << endl;
    for (auto plan_elem : plan) {
        cerr << "\t" << plan_elem.first << " " << plan_elem.second << endl;
    }
#endif
    
    // and remove the input data from the plan
    plan.resize(remove_if(plan.begin(), plan.end(), [&](const pair<string, size_t>& recipe_choice) {
        return get_index(recipe_choice.first)->is_finished();
    }) - plan.begin());
    
    return plan;
}

IndexFile::IndexFile(const string& identifier) : identifier(identifier) {
    // nothing more to do
}

bool IndexFile::is_finished() const {
    return !filenames.empty();
}

const string& IndexFile::get_identifier() const {
    return identifier;
}

const vector<string>& IndexFile::get_filenames() const {
    return filenames;
}

const vector<IndexRecipe>& IndexFile::get_recipes() const {
    return recipes;
}

void IndexFile::provide(const vector<string>& filenames) {
    this->filenames = filenames;
}

void IndexFile::execute_recipe(size_t recipe_priority) {
    assert(recipe_priority < recipes.size());
    auto& recipe = recipes[recipe_priority];
    for (auto input : recipe.inputs) {
        assert(input->is_finished());
    }
    filenames = recipe.execute();
}

void IndexFile::add_recipe(const vector<const IndexFile*>& inputs,
                           const function<vector<string>(const vector<const IndexFile*>&)>& exec) {
    recipes.emplace_back(inputs, exec);
}

IndexRecipe::IndexRecipe(const vector<const IndexFile*>& inputs,
                         const function<vector<string>(const vector<const IndexFile*>&)>& exec) :
    exec(exec), inputs(inputs)
{
    // nothing more to do
}

vector<string> IndexRecipe::execute() {
    return exec(inputs);
}

InsufficientInputException::InsufficientInputException(const string& target,
                                                       const IndexRegistry& registry) :
    runtime_error("Insufficient input to create " + target), target(target), inputs(registry.provided_indexes())
{
    // nothing else to do
}

const char* InsufficientInputException::what() const throw () {
    stringstream ss;
    ss << "Inputs" << endl;
    for (const auto& input : inputs) {
        ss << "\t" << input << endl;
    }
    ss << "are insufficient to create target index " << target << endl;
    return ss.str().c_str();
}

}

