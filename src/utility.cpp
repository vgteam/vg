#include "utility.hpp"

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

int get_thread_count(void) {
    int thread_count = 1;
#pragma omp parallel
    {
#pragma omp master
        thread_count = omp_get_num_threads();
    }
    return thread_count;
}

std::vector<std::string> &split_delims(const std::string &s, const std::string& delims, std::vector<std::string> &elems) {
    char* tok;
    char cchars [s.size()+1];
    char* cstr = &cchars[0];
    strcpy(cstr, s.c_str());
    tok = strtok(cstr, delims.c_str());
    while (tok != NULL) {
        elems.push_back(tok);
        tok = strtok(NULL, delims.c_str());
    }
    return elems;
}
std::vector<std::string> split_delims(const std::string &s, const std::string& delims) {
    std::vector<std::string> elems;
    return split_delims(s, delims, elems);
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

bool allATGC(const string& s) {
    for (string::const_iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
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

string tmpfilename(const string& base) {
    string tmpname = base + "XXXXXXXX";
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
    return tmpname;
}

string find_temp_dir() {
    // We need to find the system temp directory.
    const char* system_temp_dir = nullptr;
    
    for(const char* var_name : {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"}) {
        // Try all these env vars in order
        if (system_temp_dir == nullptr) {
            system_temp_dir = getenv(var_name);
        }
    }
    if (system_temp_dir == nullptr) {
        // Then if none were set default to /tmp
        system_temp_dir = "/tmp";
    }
    
    return std::string(system_temp_dir);
}

string tmpfilename() {
    // Make a temp file in the system temp directory.
    return tmpfilename(find_temp_dir() + "/vg");
}

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
    variant_stringer << variant.ref << '\n';
    for (auto& alt : variant.alt) {
        variant_stringer << alt << '\n';
    }
    hasher.update(variant_stringer.str());

    // Name the variant with the hex hash. Will be unique unless two
    // identical variants are in the file.
    return hasher.final();
}

double median(std::vector<int> &v) {
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    int vn = v[n];
    if (v.size()%2 == 1) {
        return vn;
    } else {
        std::nth_element(v.begin(), v.begin()+n-1, v.end());
        return 0.5*(vn+v[n-1]);
    }
}
    
vector<size_t> range_vector(size_t begin, size_t end) {
    size_t len = end - begin;
    vector<size_t> range(len, begin);
    for (size_t i = 1; i < len; i++) {
        range[i] = begin + i;
    }
    return range;
}

UnionFind::UnionFind(size_t size) {
    uf_nodes.reserve(size);
    for (size_t i = 0; i < size; i++) {
        uf_nodes.emplace_back(i);
    }
}
    
UnionFind::~UnionFind() {
    // nothing to do
}

size_t UnionFind::size() {
    return uf_nodes.size();
}

size_t UnionFind::find_group(size_t i) {
    vector<size_t> path;
    // traverse tree upwards
    while (uf_nodes[i].head != i) {
        path.push_back(i);
        i = uf_nodes[i].head;
    }
    // compress path
    unordered_set<size_t>& head_children = uf_nodes[i].children;
    for (size_t p = 1; p < path.size(); p++) {
        size_t j = path[p - 1];
        uf_nodes[j].head = i;
        uf_nodes[path[p]].children.erase(j);
        head_children.insert(j);
    }
    // note: don't need to compress path for the final index since it
    // already points to the head
    return i;
}

void UnionFind::union_groups(size_t i, size_t j) {
    size_t head_i = find_group(i);
    size_t head_j = find_group(j);
    //cerr << "union " << i << ", " << j << " in groups " << head_i << ", " << j << endl;
    if (head_i == head_j) {
        // the indices are already in the same group
        return;
    }
    else {
        // use rank as a pivot to determine which group to make the head
        UFNode& node_i = uf_nodes[head_i];
        UFNode& node_j = uf_nodes[head_j];
        if (node_i.rank > node_j.rank) {
            node_j.head = head_i;
            node_i.children.insert(head_j);
            node_i.size += node_j.size;
            //cerr << head_j << " head to " << head_i << ", " << head_i << " size to " << node_i.size << endl;
        }
        else {
            node_i.head = head_j;
            node_j.children.insert(head_i);
            node_j.size += node_i.size;
            //cerr << head_i << " head to " << head_j << ", " << head_j << " size to " << node_j.size << endl;
            
            if (node_j.rank == node_i.rank) {
                node_j.rank++;
                //cerr << head_j << " rank to " << node_j.rank << endl;
            }
        }
    }
}

size_t UnionFind::group_size(size_t i) {
    return uf_nodes[find_group(i)].size;
}

vector<size_t> UnionFind::group(size_t i) {
    vector<size_t> to_return;
    // go to head of group
    vector<size_t> stack{find_group(i)};
    // traverse tree downwards to find all indices in group
    while (!stack.empty()) {
        size_t curr = stack.back();
        stack.pop_back();
        to_return.push_back(curr);
        unordered_set<size_t>& children = uf_nodes[curr].children;
        for (size_t child : children) {
            stack.push_back(child);
        }
    }
    return to_return;
}

vector<vector<size_t>> UnionFind::all_groups() {
    vector<vector<size_t>> to_return(uf_nodes.size());
    for (size_t i = 0; i < uf_nodes.size(); i++) {
        to_return[find_group(i)].push_back(i);
    }
    auto new_end = std::remove_if(to_return.begin(), to_return.end(),
                                  [](const vector<size_t>& grp) { return grp.empty(); });
    to_return.resize(new_end - to_return.begin());
    return to_return;
}
    
string UnionFind::current_state() {
    stringstream strm;
    for (size_t i = 0; i < uf_nodes.size(); i++) {
        strm << "Node " << i << ": " << endl;
        strm << "\tHead: " << uf_nodes[i].head << endl;
        strm << "\tRank: " << uf_nodes[i].rank << endl;
        if (uf_nodes[i].head == i) {
            strm << "\tSize: " << uf_nodes[i].size << endl;
        }
        if (!uf_nodes[i].children.empty()) {
            strm << "\tChildren:" << endl;
            for (size_t child : uf_nodes[i].children) {
                strm << "\t\t" << child << endl;
            }
        }
    }
    return strm.str();
}

void get_input_file(int& optind, int argc, char** argv, function<void(istream&)> callback) {
    
    // Just combine the two operations below in the way they're supposed to be used together    
    get_input_file(get_input_file_name(optind, argc, argv), callback);

}

string get_input_file_name(int& optind, int argc, char** argv) {

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

double phi(double x1, double x2) {
    return (std::erf(x2/std::sqrt(2)) - std::erf(x1/std::sqrt(2)))/2;
}

// Modified from qnorm function in R source:
// https://svn.r-project.org/R/trunk/src/nmath/qnorm.c
double normal_inverse_cdf(double p) {
    assert(0.0 < p && p < 1.0);
    double q, r, val;
    
    q = p - 0.5;
    
    /*-- use AS 241 --- */
    /* double ppnd16_(double *p, long *ifault)*/
    /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
     
     Produces the normal deviate Z corresponding to a given lower
     tail area of P; Z is accurate to about 1 part in 10**16.
     
     (original fortran code used PARAMETER(..) for the coefficients
     and provided hash codes for checking them...)
     */
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
        val =
        q * (((((((r * 2509.0809287301226727 +
                   33430.575583588128105) * r + 67265.770927008700853) * r +
                 45921.953931549871457) * r + 13731.693765509461125) * r +
               1971.5909503065514427) * r + 133.14166789178437745) * r +
             3.387132872796366608)
        / (((((((r * 5226.495278852854561 +
                 28729.085735721942674) * r + 39307.89580009271061) * r +
               21213.794301586595867) * r + 5394.1960214247511077) * r +
             687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */
        
        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = 1.0 - p;
        else
            r = p;
        
        r = sqrt(- log(r));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
        
        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                         .0227238449892691845833) * r + .24178072517745061177) *
                       r + 1.27045825245236838258) * r +
                      3.64784832476320460504) * r + 5.7694972214606914055) *
                    r + 4.6303378461565452959) * r +
                   1.42343711074968357734)
            / (((((((r *
                     1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                   .14810397642748007459) * r + .68976733498510000455) *
                 r + 1.6763848301838038494) * r +
                2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                         2.71155556874348757815e-5) * r +
                        .0012426609473880784386) * r + .026532189526576123093) *
                      r + .29656057182850489123) * r +
                     1.7848265399172913358) * r + 5.4637849111641143699) *
                   r + 6.6579046435011037772)
            / (((((((r *
                     2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                    r + 1.8463183175100546818e-5) * r +
                   7.868691311456132591e-4) * r + .0148753612908506148525)
                 * r + .13692988092273580531) * r +
                .59983220655588793769) * r + 1.);
        }
        
        if(q < 0.0)
            val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return val;
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

// https://stackoverflow.com/a/19039500/238609
double slope(const std::vector<double>& x, const std::vector<double>& y) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto a    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    return a;
}

//https://stats.stackexchange.com/a/7459/14524
// returns alpha parameter of zipf distribution
double fit_zipf(const vector<double>& y) {
    // assume input is log-scaled
    // fit a log-log model
    assert(y.size());
    vector<double> ly(y.size());
    for (int i = 0; i < ly.size(); ++i) {
        //cerr << y[i] << " ";
        ly[i] = log(y[i]);
    }
    //cerr << endl;
    vector<double> lx(y.size());
    for (int i = 1; i <= lx.size(); ++i) {
        lx[i-1] = log(i);
    }
    return -slope(lx, ly);
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
}
