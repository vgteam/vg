#include "gcsa_helper.hpp"

#include <vg/io/vpkg.hpp>

namespace vg {

//------------------------------------------------------------------------------

void load_gcsa(gcsa::GCSA& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading GCSA from " << filename << std::endl;
    }
    std::unique_ptr<gcsa::GCSA> loaded = vg::io::VPKG::load_one<gcsa::GCSA>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_gcsa()] cannot load GCSA " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    index = std::move(*loaded);
}

void load_lcp(gcsa::LCPArray& lcp, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Loading LCP from " << filename << std::endl;
    }
    std::unique_ptr<gcsa::LCPArray> loaded = vg::io::VPKG::load_one<gcsa::LCPArray>(filename);
    if (loaded.get() == nullptr) {
        std::cerr << "error: [load_lcp()] cannot load LCP " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
    lcp = std::move(*loaded);
}

void save_gcsa(const gcsa::GCSA& index, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving GCSA to " << filename << std::endl;
    }
    if (!sdsl::store_to_file(index, filename)) {
        std::cerr << "error: [save_gcsa()] cannot write GCSA to " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

void save_lcp(const gcsa::LCPArray& lcp, const std::string& filename, bool show_progress) {
    if (show_progress) {
        std::cerr << "Saving LCP to " << filename << std::endl;
    }
    if (!sdsl::store_to_file(lcp, filename)) {
        std::cerr << "error: [save_gcsa()] cannot write LCP to " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

//------------------------------------------------------------------------------

} // namespace vg
