/** \file describe_main.cpp
 *
 * Defines the "vg describe" subcommand, which attempts to identify and describe a file based on readily available information.
 *
 * TODO: Try identifying GAM, GAF, GFA, snarls...
 * TODO: Can we get any useful information from graphs?
 * TODO: Could we do this better with vg::io?
 */

#include "subcommand.hpp"

#include "../recombinator.hpp"
#include "../xg.hpp"
#include "../zip_code.hpp"

#include <gbwt/gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/gbz.h>
#include <gbwtgraph/minimizer.h>
#include <gcsa/files.h>

#include <bdsg/hash_graph.hpp>
#include <bdsg/packed_graph.hpp>
#include <bdsg/snarl_distance_index.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

#include <getopt.h>
#include <arpa/inet.h>

using namespace vg;

//----------------------------------------------------------------------------

// Functions that attempt to describe an already identified file.
// Exceptions have been set on the input stream, and the caller will handle them.

void describe_gbwt(std::ifstream& in, const std::string& index_type, std::ostream& out);
void describe_r_index(std::ifstream& in, const std::string& index_type, std::ostream& out);
void describe_gbwtgraph(std::ifstream& in, const std::string& index_type, std::ostream& out);
void describe_gbz(std::ifstream& in, const std::string& index_type, std::ostream& out);
void describe_minimizer_index(std::ifstream& in, const std::string& index_type, std::ostream& out);

void describe_gcsa(std::ifstream& in, const std::string& index_type, std::ostream& out);
void describe_lcp(std::ifstream& in, const std::string& index_type, std::ostream& out);

void describe_haplotypes(std::ifstream& in, const std::string& index_type, std::ostream& out);

//----------------------------------------------------------------------------

template<typename T>
void load_value(std::ifstream& in, T& value) {
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
}

template<typename T>
T load_and_validate_header(std::ifstream& in) {
    T header;
    load_value(in, header);
    try {
        header.check();
    } catch (const sdsl::simple_sds::InvalidData& e) {
        throw std::runtime_error("Cannot validate header: " + std::string(e.what()));
    }
    return header;
}

//----------------------------------------------------------------------------

void help_describe(char** argv) {
    std::cerr << "usage: " << argv[0] << " " << argv[1] << " <file1> [<file2> ...]" << std::endl;
    std::cerr << std::endl;

    std::cerr << "Identify and describe files based on readily available information." << std::endl;
    std::cerr << std::endl;

    std::cerr << "Options:" << std::endl;
    std::cerr << "  -h, --help   print this help message and exit" << std::endl;
    std::cerr << std::endl;
}

void parse_args(int argc, char** argv) {
    if (argc < 3) {
        help_describe(argv);
        std::exit(EXIT_FAILURE);
    }

    static struct option long_options[] = {
        { "help", no_argument, nullptr, 'h' },
        { 0, 0, 0, 0 }
    };
    optind = 2; // force optind past command positional argument
    while (true) {
        int option_index = 0;
        int c = getopt_long(argc, argv, "h?", long_options, &option_index);
        if (c == -1) { break; } // End of options.
        switch (c)
        {
        case 'h':
        case '?':
            help_describe(argv);
            std::exit(EXIT_FAILURE);
        default:
            std::abort();
        }
    }
}

int main_describe(int argc, char** argv) {
    parse_args(argc, argv);

    // Libbdsg structures are kind of inconvenient and the files have the magic numbers in big-endian form.
    bdsg::HashGraph hg;
    std::uint32_t hg_magic = htonl(hg.get_magic_number());
    bdsg::PackedGraph pg;
    std::uint32_t pg_magic = htonl(pg.get_magic_number());
    xg::XG xg;
    std::uint32_t xg_magic = htonl(xg.get_magic_number());
    bdsg::SnarlDistanceIndex dist;
    std::uint32_t dist_magic = htonl(dist.get_magic_number());


    // Map from 32-bit magic numbers to (index type, description function) pairs.
    const std::unordered_map<std::uint32_t, std::pair<std::string, void(*)(std::ifstream&, const std::string&, std::ostream&)>> magic_map = {
        { gbwt::GBWTHeader::TAG, { "GBWT", describe_gbwt } },
        { gbwt::FastLocate::Header::TAG, { "R-index", describe_r_index } },
        { gbwtgraph::GBWTGraph::Header::TAG, { "GBWTGraph", describe_gbwtgraph } },
        { gbwtgraph::GBZ::Header::TAG, { "GBZ", describe_gbz } },
        { gbwtgraph::MinimizerHeader::TAG, { "MinimizerIndex", describe_minimizer_index } },
        { gcsa::GCSAHeader::TAG, { "GCSA", describe_gcsa } },
        { gcsa::LCPHeader::TAG, { "LCP", describe_lcp } },
        { Haplotypes::Header::MAGIC_NUMBER, { "Haplotypes", describe_haplotypes } },
        { hg_magic, { "HashGraph", nullptr } },
        { pg_magic, { "PackedGraph", nullptr } },
        { xg_magic, { "XG", nullptr } },
        { dist_magic, { "SnarlDistanceIndex", nullptr } },
        { ZipCodeCollection::get_magic_number(), { "ZipCodeCollection", nullptr } }
    };

    for (int i = 2; i < argc; ++i) {
        std::string filename = argv[i];
        std::ifstream in(filename, std::ios_base::binary);
        if (!in) {
            std::cerr << "error: [vg describe] cannot open file " << filename << std::endl;
            continue;
        }
        in.exceptions(std::ifstream::badbit | std::ifstream::failbit);

        try {
            // First try identifying the file by a magic number.
            std::uint32_t magic;
            load_value(in, magic);
            in.seekg(0); // rewind
            auto iter = magic_map.find(magic);
            if (iter != magic_map.end()) {
                std::cout << "File " << filename << " is " << iter->second.first << std::endl;
                std::cout << std::endl;
                if (iter->second.second != nullptr) {
                    // Now that we have identified the file, unexpected EOF is also an error.
                    in.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
                    iter->second.second(in, iter->second.first, std::cout);
                }
                continue;
            }

            // TODO: Other ways to identify files?

            std::cout << "File " << filename << " could not be identified" << std::endl;
            std::cout << std::endl;
        } catch (const std::runtime_error& e) {
            std::cerr << "error: [vg describe] failed to read file " << filename << ": " << e.what() << std::endl;
            continue;
        }
    }

    return 0;
}

static vg::subcommand::Subcommand vg_chains("describe", "identify and describe files", vg::subcommand::WIDGET, main_describe);

//----------------------------------------------------------------------------

void describe_gbwt(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    gbwt::GBWTHeader header = load_and_validate_header<gbwt::GBWTHeader>(in);
    bool simple_sds = header.get(gbwt::GBWTHeader::FLAG_SIMPLE_SDS);
    out << "  Version " << header.version << " in " << (simple_sds ? "Simple-SDS" : "SDSL") << " format" << std::endl;
    out << "  VG node ids in [" << ((header.offset + 1) / 2) << ", " << (header.alphabet_size / 2) << ")" << std::endl;
    out << "  GBWT node ids in [" << (header.offset + 1) << ", " << header.alphabet_size << ")" << std::endl;
    bool bidirectional = header.get(gbwt::GBWTHeader::FLAG_BIDIRECTIONAL);
    if (bidirectional) {
        out << "  Bidirectional index with " << (header.sequences / 2) << " paths of total length " << (header.size / 2) << std::endl;
    } else {
        out << "  Unidirectional index with " << header.sequences << " paths of total length " << header.size << std::endl;
    }
    bool metadata = header.get(gbwt::GBWTHeader::FLAG_METADATA);
    out << "  " << (metadata ? "Contains" : "Does not contain") << " path metadata" << std::endl;
    out << std::endl;

    out << index_type << " tags:" << std::endl;
    gbwt::Tags tags;
    if (simple_sds) {
        tags.simple_sds_load(in);
    } else {
        tags.load(in);
    }
    for (const auto& tag : tags.tags) {
        out << "  " << tag.first << " = " << tag.second << std::endl;
    }
    out << std::endl;
}

void describe_r_index(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    gbwt::FastLocate::Header header = load_and_validate_header<gbwt::FastLocate::Header>(in);
    out << "  Version " << header.version << std::endl;
    out << "  Maximum sequence length: " << header.max_length << std::endl;
    out << std::endl;
}

//----------------------------------------------------------------------------

void describe_gbwtgraph(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    gbwtgraph::GBWTGraph::Header header = load_and_validate_header<gbwtgraph::GBWTGraph::Header>(in);
    out << "  Version " << header.version << " in " << (header.get(gbwtgraph::GBWTGraph::Header::FLAG_SIMPLE_SDS) ? "Simple-SDS" : "SDSL") << " format" << std::endl;
    out << "  " << header.nodes << " nodes" << std::endl;
    out << "  " << (header.get(gbwtgraph::GBWTGraph::Header::FLAG_TRANSLATION) ? "Contains" : "Does not contain") << " segment to node translation" << std::endl;
    out << std::endl;
}

void describe_gbz(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    gbwtgraph::GBZ::Header header = load_and_validate_header<gbwtgraph::GBZ::Header>(in);
    out << "  Version " << header.version << std::endl;
    out << std::endl;

    out << index_type << " tags:" << std::endl;
    gbwt::Tags tags;
    tags.simple_sds_load(in);
    for (const auto& tag : tags.tags) {
        out << "  " << tag.first << " = " << tag.second << std::endl;
    }
    out << std::endl;

    describe_gbwt(in, "GBWT", out);
}

void describe_minimizer_index(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    gbwtgraph::MinimizerHeader header = load_and_validate_header<gbwtgraph::MinimizerHeader>(in);
    out << "  Version " << header.version << std::endl;
    if (header.flags & gbwtgraph::MinimizerHeader::FLAG_SYNCMERS) {
        out << "  Syncmers with k = " << header.k << ", s = " << header.w_or_s << std::endl;
    } else {
        size_t iterations = header.get_int(gbwtgraph::MinimizerHeader::FLAG_WEIGHT_MASK, gbwtgraph::MinimizerHeader::FLAG_WEIGHT_OFFSET);
        if (iterations > 0) {
            out << "  Weighted minimizers with k = " << header.k << ", w = " << header.w_or_s << ", iterations = " << iterations << std::endl;
        } else {
            out << "  Minimizers with k = " << header.k << ", w = " << header.w_or_s << std::endl;
        }
    }
    out << "  " << header.keys << " keys (" << header.unique << " unique) with a total of " << header.values << " occurrences" << std::endl;
    size_t key_bits = header.get_int(gbwtgraph::MinimizerHeader::FLAG_KEY_MASK, gbwtgraph::MinimizerHeader::FLAG_KEY_OFFSET);
    size_t payload_words = header.get_int(gbwtgraph::MinimizerHeader::FLAG_PAYLOAD_MASK, gbwtgraph::MinimizerHeader::FLAG_PAYLOAD_OFFSET);
    out << "  " << key_bits << "-bit keys and " << payload_words << "-word payloads" << std::endl;
    out << std::endl;

    out << index_type << " tags:" << std::endl;
    gbwt::Tags tags;
    tags.load(in);
    for (const auto& tag : tags.tags) {
        out << "  " << tag.first << " = " << tag.second << std::endl;
    }
    out << std::endl;
}

//----------------------------------------------------------------------------

void describe_gcsa(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    gcsa::GCSAHeader header = load_and_validate_header<gcsa::GCSAHeader>(in);
    out << "  Version " << header.version << std::endl;
    out << "  " << header.path_nodes << " path nodes of max length " << header.order << std::endl;
    out << "  " << header.edges << " edges" << std::endl;
    out << std::endl;
}

void describe_lcp(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    gcsa::LCPHeader header = load_and_validate_header<gcsa::LCPHeader>(in);
    out << "  Version " << header.version << std::endl;
    out << "  Size " << header.size << " with branching factor " << header.branching << std::endl;
    out << std::endl;
}

//----------------------------------------------------------------------------

void describe_haplotypes(std::ifstream& in, const std::string& index_type, std::ostream& out) {
    out << index_type << " header:" << std::endl;
    Haplotypes::Header header;
    load_value(in, header);
    // The header does not have a check() function.
    out << "  Version " << header.version << std::endl;
    out << "  " << header.construction_jobs << " GBWT construction jobs" << std::endl;
    out << "  " << header.top_level_chains << " top-level chains with " << header.total_subchains << " subchains" << std::endl;
    out << "  " << header.total_kmers << " kmers of length " << header.k << std::endl;
    out << std::endl;
}

//----------------------------------------------------------------------------
