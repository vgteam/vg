#ifndef VG_UTILITY_HPP_INCLUDED
#define VG_UTILITY_HPP_INCLUDED

#include <string>
#include <vector>
#include <sstream>
#include <omp.h>
#include <signal.h>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <unordered_set>
#include <random>
#include <type_traits>
#include <regex>
#include <signal.h>
#include <unistd.h>
#include <vg/vg.pb.h>
#include "types.hpp"
#include "sha1.hpp"
#include "Variant.h"

namespace vg {

using namespace std;

char reverse_complement(const char& c);
string reverse_complement(const string& seq);
void reverse_complement_in_place(string& seq);
/// Return True if the given string is entirely Ns of either case, and false
/// otherwise.
bool is_all_n(const string& seq);
/// Return the number of Ns as a fraction of the total sequence length
/// (or 0 if the sequence is empty)
double get_fraction_of_ns(const string& seq);
/// Return the number of threads that OMP will produce for a parallel section.
/// TODO: Assumes that this is the same for every parallel section.
int get_thread_count(void);
/// Decide on and apply a sensible OMP thread count. Pay attention to
/// OMP_NUM_THREADS if set, the "hardware concurrency", and container limit
/// information that may be available in /proc.
void choose_good_thread_count();
string wrap_text(const string& str, size_t width);
bool is_number(const string& s);

// split a string on any character found in the string of delimiters (delims)
// if max_cuts specified, only split at the first <max_cuts> delimiter occurrences
std::vector<std::string>& split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems,
                                       size_t max_cuts = numeric_limits<size_t>::max());
std::vector<std::string> split_delims(const std::string &s, const std::string& delims,
                                      size_t max_cuts = numeric_limits<size_t>::max());

/// Check if a string starts with another string
bool starts_with(const std::string& value, const std::string& prefix);

const std::string sha1sum(const std::string& data);
const std::string sha1head(const std::string& data, size_t head);

/// Return true if a character is an uppercase A, C, G, or T, and false otherwise.
bool isATGC(const char& b);
bool allATGC(const string& s);
bool allATGCN(const string& s);
string nonATGCNtoN(const string& s);
/// Convert known IUPAC ambiguity codes (which we don't support) to N (which we
/// do), while leaving any other garbage to trigger validation checks later.
string allAmbiguousToN(const string& s);
// Convert ASCII-encoded DNA to upper case
string toUppercase(const string& s);
void toUppercaseInPlace(string& s);

// write a fasta sqeuence
void write_fasta_sequence(const std::string& name, const std::string& sequence, ostream& os, size_t width=80);



template<typename T, typename V>
set<T> map_keys_to_set(const map<T, V>& m) {
    set<T> r;
    for (auto p : m) r.insert(p.first);
    return r;
}

// pairwise maximum
template<typename T>
vector<T> pmax(const std::vector<T>& a, const std::vector<T>& b) {
    std::vector<T> c;
    assert(a.size() == b.size());
    c.reserve(a.size());
    std::transform(a.begin(), a.end(), b.begin(),
                   std::back_inserter(c),
                   [](T a, T b) { return std::max<T>(a, b); });
    return c;
}

// maximum of all vectors
template<typename T>
vector<T> vpmax(const std::vector<std::vector<T>>& vv) {
    std::vector<T> c;
    if (vv.empty()) return c;
    c = vv.front();
    typename std::vector<std::vector<T> >::const_iterator v = vv.begin();
    ++v; // skip the first element
    for ( ; v != vv.end(); ++v) {
        c = pmax(c, *v);
    }
    return c;
}

/**
 * Compute the sum of the values in a collection. Values must be default-
 * constructable (like numbers are).
 */
template<typename Collection>
typename Collection::value_type sum(const Collection& collection) {

    // Set up an alias
    using Item = typename Collection::value_type;

    // Make a new zero-valued item to hold the sum
    auto total = Item();
    for(auto& to_sum : collection) {
        total += to_sum;
    }

    return total;

}


/**
 * Temporary files and directories. Create with create() or create_directory()
 * and remove with remove(). All temporary files and directories will be
 * deleted when the program exits normally or with std::exit(). The files will
 * be created in a directory determined from environment variables, though this
 * can be overridden with set_dir().
 * The interface is thread-safe.
 *
 * The temporary directory will be propagated to submodules (gbwt, gcsa2, xg).
 */
namespace temp_file {

    /// Create a temporary file starting with the given base name
    string create(const string& base);

    /// Create a temporary file
    string create();
    
    /// Create a temporary directory
    string create_directory();

    /// Remove a temporary file or directory. File or directory must have been
    /// created by create() or create_directory() and not any other means.
    void remove(const string& filename);

    /// Set a directory for placing temporary files and directories in,
    /// overriding system defaults and environment variables.
    void set_dir(const string& new_temp_dir);

    /// Reset the temporary directory back to the default based on environment
    /// variables and system defaults.
    void set_system_dir();

    /// Get the current location for temporary files and directories.
    string get_dir();

} // namespace temp_file

// Code to detect if a variant lacks an ID and give it a unique but repeatable
// one.
string get_or_make_variant_id(const vcflib::Variant& variant);
string make_variant_id(const vcflib::Variant& variant);

// TODO: move these to genotypekit on a VCF emitter?

/**
 * Create the reference allele for an empty vcflib Variant, since apaprently
 * there's no method for that already. Must be called before any alt alleles are
 * added.
 */
void create_ref_allele(vcflib::Variant& variant, const std::string& allele);

/**
 * Add a new alt allele to a vcflib Variant, since apaprently there's no method
 * for that already.
 *
 * If that allele already exists in the variant, does not add it again.
 *
 * Retuerns the allele number (0, 1, 2, etc.) corresponding to the given allele
 * string in the given variant. 
 */
int add_alt_allele(vcflib::Variant& variant, const std::string& allele);

/**
 * We have a transforming map function that we can chain.
 */ 
template <template <class T, class A = std::allocator<T>> class Container, typename Input, typename Output>
Container<Output> map_over(const Container<Input>& in, const std::function<Output(const Input&)>& lambda) {
    Container<Output> to_return;
    for (const Input& item : in) {
        to_return.push_back(lambda(item));
    }
    return to_return;
}

/**
 * We have a wrapper of that to turn a container reference into a container of pointers.
 */
template <template <class T, class A = std::allocator<T>> class Container, typename Item>
Container<const Item*> pointerfy(const Container<Item>& in) {
    return map_over<Container, Item, const Item*>(in, [](const Item& item) -> const Item* {
        return &item;
    });
}

// Simple little tree
template<typename T>
struct TreeNode {
    T v;
    vector<TreeNode<T>*> children;
    TreeNode<T>* parent;
    TreeNode() : parent(0) {}
    ~TreeNode() { for (auto c : children) { delete c; } }
    void for_each_preorder(function<void(TreeNode<T>*)> lambda) {
        lambda(this);
        for (auto c : children) {
            c->for_each_preorder(lambda);
        }
    }
    void for_each_postorder(function<void(TreeNode<T>*)> lambda) {
        for (auto c : children) {
            c->for_each_postorder(lambda);
        }
        lambda(this);
    }
};
    

    
template<typename T>
struct Tree {
    typedef TreeNode<T> Node;
    Node* root;
    Tree(Node* r = 0) : root(r) { }
    ~Tree() { delete root; }
    void for_each_preorder(function<void(Node*)> lambda) {
        if (root) root->for_each_preorder(lambda);
    }
    void for_each_postorder(function<void(Node*)> lambda) {
       if (root) root->for_each_postorder(lambda);
    }

};

// vector containing positive integer values in [begin, end)
vector<size_t> range_vector(size_t begin, size_t end);
    
// vector containing positive integer values in [0, end)
inline vector<size_t> range_vector(size_t end) {
    return range_vector(0, end);
}

struct IncrementIter {
public:
    IncrementIter(size_t number) : current(number) {
        
    }
    
    inline IncrementIter& operator=(const IncrementIter& other) {
        current = other.current;
        return *this;
    }
    
    inline bool operator==(const IncrementIter& other) const {
        return current == other.current;
    }
    
    inline bool operator!=(const IncrementIter& other) const {
        return current != other.current;
    }
    
    inline IncrementIter operator++() {
        current++;
        return *this;
    }
    
    inline IncrementIter operator++( int ) {
        IncrementIter temp = *this;
        current++;
        return temp;
    }
    
    inline size_t operator*(){
        return current;
    }
    
private:
    size_t current;
};
    
size_t integer_power(size_t x, size_t power);

/// Computes base^exponent in log(exponent) time
size_t integer_power(uint64_t base, uint64_t exponent);
/// Computes base^exponent mod modulus in log(exponent) time without requiring more
/// than 64 bits to represent exponentiated number
size_t modular_exponent(uint64_t base, uint64_t exponent, uint64_t modulus);

/// Returns a uniformly random DNA sequence of the given length
string random_sequence(size_t length);

/// Returns a uniformly random DNA sequence sequence deterministically from a seed
string pseudo_random_sequence(size_t length, uint64_t seed);

/// Escape "%" to "%25"
string percent_url_encode(const string& seq);
string replace_in_string(string subject, const string& search, const string& replace);

/// AN RNG that can skip initialization and any hashing of the seed until it is
/// needed. Not thread safe, not even a little bit.
class LazyRNG {
public:
    /// Make a new LazyRNG. Seed-generating closure will not be used after
    /// object is destroyed.
    LazyRNG(const std::function<string(void)>& get_seed);
    
    /// Get a random number, computing the seed and initilaizing the RNG if
    /// that has not yet happened.
    minstd_rand::result_type operator()();
private:
    /// Closure used to generate the seed. Makes sure to copy and not store a
    /// reference to a temporary closure.
    std::function<string(void)> get_seed;
    /// Backing RNG, or empty.
    unique_ptr<minstd_rand> rng;
};

/// Flip a coin with 50% probability against the given RNG.
bool deterministic_flip(LazyRNG& rng);

/// Given a pair of random access iterators defining a range, deterministically
/// shuffle the contents of the range based on the given RNG. Allows one RNG
/// from deterministic_start() to be used for multiple shuffles.
template<class RandomIt>
void deterministic_shuffle(RandomIt begin, RandomIt end, LazyRNG& rng) {
    // Perform Knuth shuffle algorithm using RNG
    int64_t width = end - begin;
    for (int64_t i = 1; i < width; i++) {
        std::swap(*(begin + (rng() % (i + 1))), *(begin + i));
    }
}

/// Return true if a is larger than b, or else equal to b and wins a coin flip.
template<typename Number>
bool deterministic_beats(const Number& a, const Number& b, LazyRNG& rng) {
    return (a > b || (a == b && deterministic_flip(rng)));
}

// For seed generation, we aren't just using our normal hash functions, because
// we want it to be more semantic and not dependent on memory addresses.

/// Make seeds for Alignments based on their sequences.
inline string make_shuffle_seed(const Alignment& aln) {
    return aln.sequence();
}

/// Make seeds for pointers to things we can make seeds for.
template<typename T>
inline string make_shuffle_seed(const T* ptr) {
    return make_shuffle_seed(*ptr);
}

/// Make seeds for pairs of things we can make seeds for.
template<typename T1, typename T2>
inline string make_shuffle_seed(const pair<T1, T2>& p) {
    return make_shuffle_seed(p.first) + make_shuffle_seed(p.second);
}

/// Do a deterministic shuffle with automatic seed determination.
template<class RandomIt>
void deterministic_shuffle(RandomIt begin, RandomIt end) {
    LazyRNG rng([&]() {
        return make_shuffle_seed(*begin);
    });

    deterministic_shuffle(begin, end, rng);
}

/**
 * Sort the items between the two given random-access iterators, as with
 * std::sort. Deterministically shuffle the ties, if any, at the top end.
 */
template<class RandomIt, class Compare>
void sort_shuffling_ties(RandomIt begin, RandomIt end, Compare comp, LazyRNG& rng) {
    
    // Sort everything
    std::stable_sort(begin, end, comp);
    
    // Comparison returns true if first argument must come before second, and
    // false otherwise. So the ties will be a run where the top thing doesn't
    // necessarily come before each other thing (i.e. comparison returns
    // false).
    
    // Count the ties at the top
    RandomIt ties_end = begin;
    while (ties_end != end && !comp(*begin, *ties_end)) {
        // We haven't hit the end of the list, and the top thing isn't strictly better than this thing.
        // So mark it as a tie and advance.
        ++ties_end;
    }
    
    if (begin != ties_end) {
        // Shuffle the ties.
        deterministic_shuffle(begin, ties_end, rng);
    }

}

/**
 * Sort the items between the two given random-access iterators, as with
 * std::sort. Deterministically shuffle the ties, if any, at the top end, using
 * automatic seed determination as defined by a make_shuffle_seed() overload
 * for the collection's item type.
 */
template<class RandomIt, class Compare>
void sort_shuffling_ties(RandomIt begin, RandomIt end, Compare comp) {

    LazyRNG rng([&]() {
        return make_shuffle_seed(*begin);
    });
    
    // Make the seed and start the RNG using the pre-defined seed making
    // approaches
    sort_shuffling_ties(begin, end, comp, rng);

}

/// Compose the translations from two graph operations, both of which involved oriented transformations.
unordered_map<id_t, pair<id_t, bool>> overlay_node_translations(const unordered_map<id_t, pair<id_t, bool>>& over,
                                                                const unordered_map<id_t, pair<id_t, bool>>& under);

/// Compose the translations from two graph operations, the first of which involved oriented transformations.
unordered_map<id_t, pair<id_t, bool>> overlay_node_translations(const unordered_map<id_t, id_t>& over,
                                                                const unordered_map<id_t, pair<id_t, bool>>& under);

/// Compose the translations from two graph operations, the second of which involved oriented transformations.
unordered_map<id_t, pair<id_t, bool>> overlay_node_translations(const unordered_map<id_t, pair<id_t, bool>>& over,
                                                                const unordered_map<id_t, id_t>& under);

/// Compose the translations from two graph operations, neither of which involved oriented transformations.
unordered_map<id_t, id_t> overlay_node_translations(const unordered_map<id_t, id_t>& over,
                                                    const unordered_map<id_t, id_t>& under);
    

/// Return true if there's a command line argument (i.e. input file name) waiting to be processed. 
bool have_input_file(int& optind, int argc, char** argv);

/// Get a callback with an istream& to an open file if a file name argument is
/// present after the parsed options, or print an error message and exit if one
/// is not. Handles "-" as a filename as indicating standard input. The reference
/// passed is guaranteed to be valid only until the callback returns. Bumps up
/// optind to the next argument if a filename is found.
///
/// Warning: If you're reading a HandleGraph via VPKG::load_one (as is the pattern in vg)
///          it is best to use get_input_file_name() below instead, and run load_one on that.
///          This allows better GFA support because it allows memmapping the file directly
void get_input_file(int& optind, int argc, char** argv, function<void(istream&)> callback);

/// Parse out the name of an input file (i.e. the next positional argument), or
/// throw an error. File name must be nonempty, but may be "-" or may not exist.
string get_input_file_name(int& optind, int argc, char** argv, bool test_open = true);

/// Parse out the name of an output file (i.e. the next positional argument), or
/// throw an error. File name must be nonempty.
string get_output_file_name(int& optind, int argc, char** argv);

/// Get a callback with an istream& to an open file. Handles "-" as a filename as
/// indicating standard input. The reference passed is guaranteed to be valid
/// only until the callback returns.
void get_input_file(const string& file_name, function<void(istream&)> callback);

/// Split off the extension from a filename and return both parts. 
pair<string, string> split_ext(const string& filename);

/// Determine if a file exists.
/// Only works for files readable by the current user.
bool file_exists(const string& filename);

/// Parse a command-line argument string. Exits with an error if the string
/// does not contain exactly an item fo the appropriate type.
template<typename Result>
Result parse(const string& arg);

/// Parse a command-line argument C string. Exits with an error if the string
/// does not contain exactly an item fo the appropriate type.
template<typename Result>
Result parse(const char* arg);

/// Parse the appropriate type from the string to the destination value.
/// Return true if parsing is successful and false (or throw something) otherwise.
template<typename Result>
bool parse(const string& arg, Result& dest);

// Do one generic implementation for signed integers that fit in a long long.
// Cram the constraint into the type of the output parameter.
template<typename Result>
bool parse(const string& arg, typename enable_if<sizeof(Result) <= sizeof(long long) &&
    is_integral<Result>::value &&
    is_signed<Result>::value, Result>::type& dest) {
    
    // This will hold the next character after the number parsed
    size_t after;
    long long buffer = std::stoll(arg, &after);
    if (buffer > numeric_limits<Result>::max() || buffer < numeric_limits<Result>::min()) {
        // Out of range
        return false;
    }
    dest = (Result) buffer;
    return(after == arg.size());    
}

// Do another generic implementation for unsigned integers
template<typename Result>
bool parse(const string& arg, typename enable_if<sizeof(Result) <= sizeof(unsigned long long) &&
    is_integral<Result>::value &&
    !is_signed<Result>::value, Result>::type& dest) {
    
    // This will hold the next character after the number parsed
    size_t after;
    unsigned long long buffer = std::stoull(arg, &after);
    if (buffer > numeric_limits<Result>::max() || buffer < numeric_limits<Result>::min()) {
        // Out of range
        return false;
    }
    dest = (Result) buffer;
    return(after == arg.size());    
}              

// We also have an implementation for doubles (defined in the cpp)
template<>
bool parse(const string& arg, double& dest);

// And one for regular expressions
template<>
bool parse(const string& arg, std::regex& dest);

// Implement the first version in terms of the second, for any type
template<typename Result>
Result parse(const string& arg) {
    Result to_return;
    bool success;
    try {
        success = parse<Result>(arg, to_return);
    } catch(exception& e) {
        success = false;
    }
    if (success) {
        // Parsing worked
        return to_return;
    } else {
        // Parsing failed
        cerr << "error: could not parse " << typeid(to_return).name() << " from argument \"" << arg << "\"" << endl;
        exit(1);
    }
}

// Implement the C string version in terms of that
template<typename Result>
Result parse(const char* arg) {
    return parse<Result>(string(arg));
}
 
}

#endif
