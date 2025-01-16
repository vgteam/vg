#include "utility.hpp"
#include "statistics.hpp"

#include <set>
#include <mutex>
#include <dirent.h>
#include <thread>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cctype>
// We don't define _GNU_SOURCE to get the cpuset functions since we will
// already have it for libstdc++ on the platforms where we need them
#include <sched.h>


// For setting the temporary directory in submodules.
#include <gcsa/utils.h>
#include <gbwt/utils.h>
#include <xg.hpp>

namespace vg {

static const char complement[256] = {'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 16
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 24
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 32
                                     'N', 'N', 'N', '$', '#', 'N', 'N', 'N', // 40 GCSA stop/start characters
                                     'N', 'N', 'N', 'N', 'N', '-', 'N', 'N', // 48
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 56
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 64
                                     'N', 'T', 'V', 'G', 'H', 'N', 'N', 'C', // 72
                                     'D', 'N', 'N', 'M', 'N', 'K', 'N', 'N', // 80
                                     'N', 'Q', 'Y', 'W', 'A', 'A', 'B', 'S', // 88
                                     'N', 'R', 'N', 'N', 'N', 'N', 'N', 'N', // 96
                                     'N', 't', 'v', 'g', 'h', 'N', 'N', 'c', // 104
                                     'd', 'N', 'N', 'm', 'N', 'k', 'n', 'N', // 112
                                     'N', 'q', 'y', 'w', 'a', 'a', 'b', 's', // 120
                                     'N', 'r', 'N', 'N', 'N', 'N', 'N', 'N', // 128
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 136
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 144
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 152
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 160
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 168
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 176
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 184
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 192
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 200
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 208
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 216
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 224
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 232
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 240
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 248
                                     'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};// 256

char reverse_complement(const char& c) {
    return complement[c];
}

string reverse_complement(const string& seq) {
    string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        c = complement[c];
    }
    return rc;
}
    
void reverse_complement_in_place(string& seq) {
    size_t swap_size = seq.size() / 2;
    for (size_t i = 0, j = seq.size() - 1; i < swap_size; i++, j--) {
        char tmp = seq[i];
        seq[i] = complement[seq[j]];
        seq[j] = complement[tmp];
    }
    
    if (seq.size() % 2) {
        seq[swap_size] = complement[seq[swap_size]];
    }
}

bool is_all_n(const string& seq) {
    for (auto& c : seq) {
        if (c != 'N' && c != 'n') {
            return false;
        }
    }
    return true;
}

double get_fraction_of_ns(const string& seq) {
    double n_frac = 0.;
    if (!seq.empty()) {
        size_t n_count = 0;
        for (char c : seq) {
            if (c == 'n' || c == 'N') {
                ++n_count;
            }
        }
        n_frac = (double)n_count/(double) seq.length();
    }
    return n_frac;
}

int get_thread_count(void) {
    int thread_count = 1;
#pragma omp parallel
    {
#pragma omp master
        thread_count = omp_get_num_threads();
    }
    return thread_count;
}

void choose_good_thread_count() {
    // If we leave this at 0, we won't apply any thread count and leave whatever OMP defaults to.
    int count = 0;

    if (count == 0) {
        // First priority: OMP_NUM_THREADS
        const char* value = getenv("OMP_NUM_THREADS");
        if (value && *value != '\0') {
            // Read the value. Throws if it isn't a legit number.
            count = std::stoi(value);
        }
    }

    if (count == 0) {
        // Next priority: /sys/fs/cgroup/cpu/cpu.cfs_quota_us over /sys/fs/cgroup/cpu/cpu.cfs_period_us
        ifstream quota_file("/sys/fs/cgroup/cpu/cpu.cfs_quota_us");
        ifstream period_file("/sys/fs/cgroup/cpu/cpu.cfs_period_us");

        if (quota_file && period_file) {
            // Read the period and quota
            int64_t quota;
            quota_file >> quota;
            int64_t period;
            period_file >> period;

            if (quota >= 0 && period != 0) {
                // Compute how many threads we may use.
                // May come out to 0, in which case it is ignored.
                count = (int) ceil(quota / (double) period);
            }
        }
    }

#if !defined(__APPLE__) && defined(_GNU_SOURCE)
    if (count == 0) {
        // Next priority: CPU affinity mask (used by Slurm)
        cpu_set_t mask;
        if (sched_getaffinity(getpid(), sizeof(cpu_set_t), &mask)) {
            // TODO: If you have >1024 bits in your mask, glibc can't deal and you will get EINVAL.
            // We're supposed to then try increasingly large dynamically-allocated CPU flag sets until we find one that works.
            auto problem = errno;
            std::cerr << "warning[vg]: Cannot determine CPU count from affinity mask: " << strerror(problem) << std::endl;
        } else {
            // We're also supposed to intersect this mask with the actual
            // existing processors, in case somebody flags on way more
            // processors than actually exist. But Linux doesn't seem to do
            // that by default, so we don't worry about it.
            count = CPU_COUNT(&mask);
        }
    }
#endif

    if (count == 0) {
        // Next priority: SLURM_JOB_CPUS_PER_NODE
        const char* value = getenv("SLURM_JOB_CPUS_PER_NODE");
        if (value && *value != '\0') {
            // Read the value. Throws if it isn't a legit number.
            count = std::stoi(value);
        }
    }

    if (count == 0) {
        // Next priority: hardware concurrency as reported by the STL.
        // This may itself be 0 if ungettable.
        count = std::thread::hardware_concurrency();
    }
    
    if (count != 0) {
        omp_set_num_threads(count);
    }
}

std::vector<std::string> &split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems, size_t max_cuts) {
    size_t start = string::npos;
    size_t cuts = 0;
    for (size_t i = 0; i < s.size(); ++i) {
        if (delims.find(s[i]) != string::npos && cuts < max_cuts) {
            if (start != string::npos && i > start) {
                elems.push_back(s.substr(start, i - start));
            }
            start = string::npos;
            ++cuts;
        } else if (start == string::npos) {
            start = i;
        }
    }
    if (start != string::npos && start < s.size()) {
        elems.push_back(s.substr(start, s.size() - start));
    }
    return elems;
}

std::vector<std::string> split_delims(const std::string &s, const std::string& delims, size_t max_cuts) {
    std::vector<std::string> elems;
    return split_delims(s, delims, elems, max_cuts);
}

bool starts_with(const std::string& value, const std::string& prefix) {
#if __cplusplus > 201703L
    // C++20 provides this
    return value.starts_with(prefix);
#else
    // Before then, C++ is terrible and we have to do it ourselves.
    if (value.size() < prefix.size()) {
        return false;
    }
    return std::equal(prefix.begin(), prefix.end(), value.begin());
#endif
}

const std::string sha1sum(const std::string& data) {
    SHA1 checksum;
    checksum.update(data);
    return checksum.final();
}

const std::string sha1head(const std::string& data, size_t head) {
    return sha1sum(data).substr(0, head);
}

string wrap_text(const string& str, size_t width) {
    stringstream w;
    size_t j = 0;
    for (auto c : str) {
        if (j++ > 50) {
            if (c == ' ') {
                w << "\n";
                j = 0;
            } else {
                w << c;
            }
        } else {
            w << c;
        }
    }
    return w.str();
}

bool is_number(const std::string& s) {
    return !s.empty()
        && std::find_if(s.begin(), s.end(),
                        [](char c) { return !std::isdigit(c); }) == s.end();
}

bool isATGC(const char& b) {
    return (b == 'A' || b == 'T' || b == 'G' || b == 'C');
}

bool allATGC(const string& s) {
    for (string::const_iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
}

bool allATGCN(const string& s) {
    for (string::const_iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C' && b != 'N') {
            return false;
        }
    }
    return true;
}

string nonATGCNtoN(const string& s) {
    auto n = s;
    for (string::iterator c = n.begin(); c != n.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C' && b != 'N') {
            *c = 'N';
        }
    }
    return n;
}

string allAmbiguousToN(const string& s) {
    auto n = s;
    for (string::iterator c = n.begin(); c != n.end(); ++c) {
        char b = *c;
        if (b == 'M' || b == 'R' || b == 'W' || b == 'S' || b == 'Y' ||
            b == 'K' || b == 'V' || b == 'H' || b == 'D' || b == 'B') {
            // Replace known IUPAC ambiguity codes with N.
            // Leave weird things like '-' or '*' which should be errors.
            *c = 'N';
        }
    }
    return n;
}

string toUppercase(const string& s) {
    auto n = s;
    for (string::iterator c = n.begin(); c != n.end(); ++c) {
        if (*c >= 'a' && *c <= 'z') {
            *c -= 'a' - 'A';
        }
    }
    return n;
}

void toUppercaseInPlace(string& s) {
    for (int i = 0; i < s.size(); ++i) {
        s[i] = toupper(s[i]);
    }
}

void write_fasta_sequence(const std::string& name, const std::string& sequence, ostream& os, size_t width) {
    os << ">" << name << "\n";
    for (size_t written = 0; written < sequence.length(); written += width) {
        os << sequence.substr(written, min(width, sequence.length() - written)) << "\n";
    }
    os << flush;
}

namespace temp_file {

// We use this to make the API thread-safe
static recursive_mutex monitor;

static string temp_dir;

/// Because the names are in a static object, we can delete them when
/// std::exit() is called.
struct Handler {
    /// What temp files have we handed out?
    unordered_set<string> filenames;
    
    /// What temp directories have we handed out
    unordered_set<string> dirnames;
    
    /// Place where all temporary files for this vg invocation live.
    string parent_directory;
    
    /// Make sure the parent directory for all our temporary files and
    /// directories exists.
    void ensure_parent_directory() {
        if (parent_directory.empty()) {
            // Make a parent directory for our temp files
            string tmpdirname = get_dir() + "/vg-XXXXXX";
            auto got = mkdtemp(&tmpdirname[0]);
            if (got != nullptr) {
                // Save the directory we got
                parent_directory = got;
            } else {
                cerr << "[vg utility.cpp]: couldn't create temp directory: " << tmpdirname << endl;
                exit(1);
            }
        }
    }
    
    /// Delete a directory and all files in it.
    static void remove_directory(const string& name) {
        // Open it up to get the files
        auto directory = opendir(name.c_str());
        
        dirent* dp;
        while ((dp = readdir(directory)) != nullptr) {
            // For every item still in it
            
            if (strcmp(dp->d_name, ".") == 0 || strcmp(dp->d_name, "..") == 0) {
                // This is a special here or up entry, so skip it.
                continue;
            }
            
            // Compute the full path
            string path = name + "/" + dp->d_name;
            
            struct stat dp_stat;
            stat(path.c_str(), &dp_stat);
            if (S_ISDIR(dp_stat.st_mode) && !S_ISLNK(dp_stat.st_mode)) {
                // It's a directory and may have stuff in it.
                // It isn't a link to some other random place (unless the
                // current user is tampering with their own temp files to
                // delete their own stuff).
                // Clear it out.
                remove_directory(path);
            } else {
                // Normal file or symlink.
                // Delete just it.
                std::remove(path.c_str());
            }
        }
        closedir(directory);
        
        // Delete the directory itself
        std::remove(name.c_str());
    }

    ~Handler() {
        // No need to lock in static destructor
        for (auto& filename : filenames) {
            std::remove(filename.c_str());
        }
        for (auto& dirname : dirnames) {
            remove_directory(dirname);
        }
        if (!parent_directory.empty()) {
            // There may be extraneous files in the directory still (like .fai files)
            remove_directory(parent_directory);
        }
    }
};
// make a static instance so that its destructor is called from std::exit()
static Handler handler;

string create(const string& base) {
    lock_guard<recursive_mutex> lock(monitor);

    handler.ensure_parent_directory(); 

    string tmpname = handler.parent_directory + "/" + base + "XXXXXX";
    // hack to use mkstemp to get us a safe temporary file name
    int fd = mkstemp(&tmpname[0]);
    if(fd != -1) {
        // we don't leave it open; we are assumed to open it again externally
        close(fd);
    } else {
        cerr << "[vg utility.cpp]: couldn't create temp file on base "
             << base << " : " << tmpname << endl;
        exit(1);
    }
    handler.filenames.insert(tmpname);
    return tmpname;
}

string create() {
    // No need to lock as we call this thing that locks
    return create("vg-");
}

string create_directory() {
    lock_guard<recursive_mutex> lock(monitor);
    
    handler.ensure_parent_directory(); 

    string tmpname = handler.parent_directory + "/dir-XXXXXX";
    auto got = mkdtemp(&tmpname[0]);
    if (got != nullptr) {
        // Save the directory we got
        handler.dirnames.insert(got);
        return got;
    } else {
        cerr << "[vg utility.cpp]: couldn't create temp directory: " << tmpname << endl;
        exit(1);
    }
}

void remove(const string& filename) {
    lock_guard<recursive_mutex> lock(monitor);
    
    if (handler.dirnames.count(filename)) {
        handler.remove_directory(filename);
        handler.dirnames.erase(filename);
    } else if (handler.filenames.count(filename)) {
        std::remove(filename.c_str());
        handler.filenames.erase(filename);
    } else {
        // Probably e.g. a .fai file. Just remove it and fail if it's a
        // nonempty directory.
        std::remove(filename.c_str());
    }
}

void forget() {
    lock_guard<recursive_mutex> lock(monitor);
    
    // forget in our submodules as well
    xg::temp_file::forget();
    gbwt::TempFile::forget();
    gcsa::TempFile::forget();
    
    handler.filenames.clear();
    handler.dirnames.clear();
    handler.parent_directory.clear();
}

void set_system_dir() {
    const char* system_temp_dir = nullptr;
    for(const char* var_name : {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"}) {
        if (system_temp_dir == nullptr) {
            system_temp_dir = getenv(var_name);
        }
    }
    set_dir(system_temp_dir == nullptr ? "/tmp" : system_temp_dir);
}

void set_dir(const string& new_temp_dir) {
    lock_guard<recursive_mutex> lock(monitor);
    temp_dir = new_temp_dir;

    // Several submodules use their own temporary directories.
    gcsa::TempFile::setDirectory(temp_dir);
    gbwt::TempFile::setDirectory(temp_dir);
    xg::temp_file::set_dir(temp_dir);
}

string get_dir() {
    lock_guard<recursive_mutex> lock(monitor);
    return temp_dir;
}

} // namespace temp_file

string get_or_make_variant_id(const vcflib::Variant& variant) {

     if(!variant.id.empty() && variant.id != ".") {
        // We assume all the actually filled in ID fields in a VCF are unique.
        return variant.id;
    } else {
        // Synthesize a name for the variant

        return make_variant_id(variant);

    }
}

string make_variant_id(const vcflib::Variant& variant) {
    // Synthesize a name for the variant

    // Let's just hash
    SHA1 hasher;

    // Turn the variant into a string, leaving out the actual calls and any
    // assigned ID. Note that this keeps the modified 0-based position.
    std::stringstream variant_stringer;
    variant_stringer << variant.sequenceName << '\n';
    variant_stringer << variant.position << '\n';
    // Case insensitive so we can use this function to get ids consistent with
    // vg construct (which converts to uppercase before assigning ids)
    variant_stringer << toUppercase(variant.ref) << '\n';
    for (auto& alt : variant.alt) {
      variant_stringer << toUppercase(alt) << '\n';
    }
    hasher.update(variant_stringer.str());

    // Name the variant with the hex hash. Will be unique unless two
    // identical variants are in the file.
    return hasher.final();
}
    
vector<size_t> range_vector(size_t begin, size_t end) {
    size_t len = end - begin;
    vector<size_t> range(len, begin);
    for (size_t i = 1; i < len; i++) {
        range[i] = begin + i;
    }
    return range;
}

std::vector<size_t> stack_permutations(const std::vector<size_t>& bottom, const std::vector<size_t>& top) {
    std::vector<size_t> result;
    result.reserve(top.size());
    for (auto& index : top) {
        result.push_back(bottom[index]);
    }
    return result;
}

bool have_input_file(int& optind, int argc, char** argv) {

    if (optind >= argc) {
        // Out of arguments
        return false;
    }
    
    if (argv[optind][0] == '\0') {
        // File name is empty
        return false;
    }
    
    // Otherwise we found one
    return true;
}

void get_input_file(int& optind, int argc, char** argv, function<void(istream&)> callback) {
    
    // Just combine the two operations below in the way they're supposed to be used together    
    get_input_file(get_input_file_name(optind, argc, argv), callback);

}

string get_input_file_name(int& optind, int argc, char** argv, bool test_open) {

    if (optind >= argc) {
        // Complain that the user didn't specify a filename
        cerr << "error:[get_input_file_name] specify input filename, or \"-\" for standard input" << endl;
        exit(1);
    }
    
    string file_name(argv[optind++]);
    
    if (file_name.empty()) {
        cerr << "error:[get_input_file_name] specify a non-empty input filename" << endl;
        exit(1);
    }

    if (test_open && file_name != "-") {
        ifstream file_stream(file_name);
        if (!file_stream) {
            cerr << "error:[get_input_file_name] unable to open input file: " << file_name << endl;
        }
    }
    
    return file_name;
    
}

string get_output_file_name(int& optind, int argc, char** argv) {

    if (optind >= argc) {
        // Complain that the user didn't specify a filename
        cerr << "error:[get_output_file_name] specify output filename" << endl;
        exit(1);
    }
    
    string file_name(argv[optind++]);
    
    if (file_name.empty()) {
        cerr << "error:[get_output_file_name] specify a non-empty output filename" << endl;
        exit(1);
    }
    
    return file_name;
    
}

void get_input_file(const string& file_name, function<void(istream&)> callback) {

    if (file_name == "-") {
        // Just use standard input
        callback(std::cin);
    } else {
        // Open a file
        ifstream in;
        in.open(file_name.c_str());
        if (!in.is_open()) {
            // The user gave us a bad filename
            cerr << "error:[get_input_file] could not open file \"" << file_name << "\"" << endl;
            exit(1);
        }
        callback(in);
    }
    
}

pair<string, string> split_ext(const string& filename) {
    pair<string, string> parts;

    size_t dot = filename.rfind('.');
    
    if (dot == string::npos) {
        // Put it all in the first part
        parts.first = filename;
    } else {
        // Split on either side of the dot.
        parts.first = filename.substr(0, dot);
        parts.second = filename.substr(dot + 1);
    }
    
    return parts;
}

string file_base_name(const string& filename) {
    size_t slash = filename.rfind('/');
    if (slash == string::npos) {
        return split_ext(filename).first;
    } else {
        return split_ext(filename.substr(slash + 1)).first;
    }
}

bool file_exists(const string& filename) {
    // TODO: use C++17 features to actually poll existence.
    // For now we see if we can open it.
    ifstream in(filename);
    return in.is_open();
}

    
void create_ref_allele(vcflib::Variant& variant, const std::string& allele) {
    // Set the ref allele
    variant.ref = allele;
    
    for(size_t i = 0; i < variant.ref.size(); i++) {
        // Look at all the bases
        if(variant.ref[i] != 'A' && variant.ref[i] != 'C' && variant.ref[i] != 'G' && variant.ref[i] != 'T') {
            // Correct anything bogus (like "X") to N
            variant.ref[i] = 'N';
        }
    }
    
    // Make it 0 in the alleles-by-index list
    variant.alleles.push_back(allele);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();
}

int add_alt_allele(vcflib::Variant& variant, const std::string& allele) {
    // Copy the allele so we can throw out bad characters
    std::string fixed(allele);
    
    for(size_t i = 0; i < fixed.size(); i++) {
        // Look at all the bases
        if(fixed[i] != 'A' && fixed[i] != 'C' && fixed[i] != 'G' && fixed[i] != 'T') {
            // Correct anything bogus (like "X") to N
            fixed[i] = 'N';
        }
    }
    
    for(int i = 0; i < variant.alleles.size(); i++) {
        if(variant.alleles[i] == fixed) {
            // Already exists
            return i;
        }
    }

    // Add it as an alt
    variant.alt.push_back(fixed);
    // Make it next in the alleles-by-index list
    variant.alleles.push_back(fixed);
    // Build the reciprocal index-by-allele mapping
    variant.updateAlleleIndexes();

    // We added it in at the end
    return variant.alleles.size() - 1;
}

// exponentiation by squaring as implemented in
// https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
size_t integer_power(uint64_t base, uint64_t exponent) {
    size_t result = 1;
    while (exponent) {
        if (exponent & 1) {
            result *= base;
        }
        exponent >>= 1;
        base *= base;
    }
    return result;
}

// modular exponent by squaring as described in
// https://en.wikipedia.org/wiki/Modular_exponentiation`
size_t modular_exponent(uint64_t base, uint64_t exponent, uint64_t modulus) {
    if (modulus == 1) {
        return 0;
    }
    size_t result = 1;
    base = base % modulus;
    while (exponent) {
        if (exponent & 1) {
            result = (result * base) % modulus;
        }
        exponent >>= 1;
        base = (base * base) % modulus;
    }
    return result;
}
    
default_random_engine random_sequence_gen(102);
    
string random_sequence(size_t length) {
    static const string alphabet = "ACGT";
    uniform_int_distribution<char> distr(0, 3);
    
    string seq(length, '\0');
    for (size_t i = 0; i < length; i++) {
        seq[i] = alphabet[distr(random_sequence_gen)];
    }
    return seq;
}

string pseudo_random_sequence(size_t length, uint64_t seed) {
    static const string alphabet = "ACGT";
    mt19937_64 gen(1357908642ull * seed + 80085ull);
    vg::uniform_int_distribution<char> distr(0, 3);
    
    string seq(length, '\0');
    for (size_t i = 0; i < length; i++) {
        seq[i] = alphabet[distr(gen)];
    }
    return seq;
}

string replace_in_string(string subject,
                         const string& search,
                         const string& replace) {
    size_t pos = 0;
    while ((pos = subject.find(search, pos)) != string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}

string percent_url_encode(const string& seq) {
    return replace_in_string(seq, "%", "%25");
}

LazyRNG::LazyRNG(const std::function<string(void)>& get_seed) : get_seed(get_seed) {
    // Nothing to do
}

minstd_rand::result_type LazyRNG::operator()() {
    if (!rng) {
        // Make sure the RNG is initialized
        string seed = get_seed();
        
        // Turn the string into a 32-bit number.
        uint32_t seedNumber = 0;
        for (uint8_t byte : seed) {
            // Sum up with primes and overflow.
            // TODO: this is a bit of a bad hash function but it should be good enough.
            seedNumber = seedNumber * 13 + byte;
        }
        
        rng = make_unique<minstd_rand>(seedNumber);
    }
    return (*rng)();
}

bool deterministic_flip(LazyRNG& rng) {
    return rng() % 2;
}

unordered_map<id_t, pair<id_t, bool>> overlay_node_translations(const unordered_map<id_t, pair<id_t, bool>>& over,
                                                                const unordered_map<id_t, pair<id_t, bool>>& under) {
    
    unordered_map<id_t, pair<id_t, bool>> overlaid;
    overlaid.reserve(over.size());
    for (const pair<id_t, pair<id_t, bool>>& node_trans : over) {
        const pair<id_t, bool>& trans_thru = under.at(node_trans.second.first);
        overlaid[node_trans.first] = make_pair(trans_thru.first, trans_thru.second != node_trans.second.second);
    }
    return overlaid;
}

unordered_map<id_t, pair<id_t, bool>> overlay_node_translations(const unordered_map<id_t, id_t>& over,
                                                                const unordered_map<id_t, pair<id_t, bool>>& under) {
    unordered_map<id_t, pair<id_t, bool>> overlaid;
    overlaid.reserve(over.size());
    for (const pair<id_t, id_t>& node_trans : over) {
        overlaid[node_trans.first] = under.at(node_trans.second);
    }
    return overlaid;
}

unordered_map<id_t, pair<id_t, bool>> overlay_node_translations(const unordered_map<id_t, pair<id_t, bool>>& over,
                                                                const unordered_map<id_t, id_t>& under) {
    unordered_map<id_t, pair<id_t, bool>> overlaid;
    overlaid.reserve(over.size());
    for (const pair<id_t, pair<id_t, bool>>& node_trans : over) {
        overlaid[node_trans.first] = make_pair(under.at(node_trans.second.first), node_trans.second.second);
    }
    return overlaid;
}

unordered_map<id_t, id_t> overlay_node_translations(const unordered_map<id_t, id_t>& over,
                                                    const unordered_map<id_t, id_t>& under) {
    unordered_map<id_t, id_t> overlaid;
    overlaid.reserve(over.size());
    for (const pair<id_t, id_t>& node_trans : over) {
        overlaid[node_trans.first] = under.at(node_trans.second);
    }
    return overlaid;
}

template<>
bool parse(const string& arg, double& dest) {
    size_t after;
    dest = std::stod(arg, &after);
    return(after == arg.size());
}

template<>
bool parse(const string& arg, float& dest) {
    size_t after;
    dest = std::stof(arg, &after);
    return(after == arg.size());
}


template<>
bool parse(const string& arg, std::regex& dest) {
    // This throsw std::regex_error if it can't parse.
    // That contains a kind of useless error code that we can't turn itno a string without switching on all the values.
    dest = std::regex(arg);
    return true;
}

template<>
bool parse(const string& arg, pos_t& dest) {
    const char* cursor = arg.c_str();
    // The STL cheats and breaks constness here.
    char* next = nullptr;
    // Read the node ID
    get_id(dest) = std::strtoull(cursor, &next, 10);
    if (next == nullptr || *next == '\0' || id(dest) == 0) {
        // Out of parts or didn't get an ID
        return false;
    }
    if (*next == '+') {
        get_is_rev(dest) = false;
    } else if (*next == '-') {
        get_is_rev(dest) = true;
    } else {
        // No separator character
        return false;
    }
    cursor = next;
    ++cursor;
    if (*cursor == '\0') {
        // No offset
        return false;
    }
    // Parse the offset
    get_offset(dest) = std::strtoull(cursor, &next, 10);
    if (next == nullptr || *next != '\0') {
        // We didn't consume the rest of the string
        return false;
    }
    return true;
}

}
