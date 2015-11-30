#ifndef DECON_HPP
#define DECON_HPP
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include <queue>
#include "Variant.h"
#include "index.hpp"
#include "path.hpp"
#include "vg.hpp"
#include "vg.pb.h"
#include "Fasta.h"
#include "xg.hpp"
#include "position.hpp"


namespace vg{
    using namespace std;
    class Deconstructor{
        public:

            Deconstructor();
            ~Deconstructor();
            void clear();
            void set_xg(string xg);
            void enumerate_path_names_in_index();
            void set_reference(string ref_file);
            void set_index(string index_file);
            void set_graph(VG* v);

            /**
             * Project a path onto another path,
             * much like a projection of p2 onto p1 in space.
             */
             vector<Mapping> relative_mapping(Path& p1, Path& p2);

            // bool surject_alignment(Alignment& source,
            //                         set<string> path_names,
            //                         Alignment& surjection,
            //                         string& path_name,
            //                         int64_t& path_pos,
            //                         int window);

            /**
             * Convenience function.
             * Project a mapping with name 'p2' onto p1.
             * Return the projection as a mapping (i.e. an alignment of p2 on p1).
             */
            // Path relative_mapping(Path& p1, string p2);
            void enumerate_graph_paths();
            /**
            *
            */
            vector<vcflib::Variant> get_variants(string region_name);
            /**
             * Build a vcf record from two paths, with the
             * second path argument taken as the reference.
             */
            vcflib::Variant path_to_variants(Path variant, Path ref);
            vcflib::Variant pathname_to_variants(string variant, Path ref);
            /**
             * Turn a vector of variants into a proper VCF.
             */
            void write_variants(string filename, vector<vcflib::Variant> variants);
        private:
            // TODO Should be able to handle XG or VG indices
            string index_file;
            string reference;
            string xg_file;
            VG* vgraph;
            vector<string> ref_paths;
            map<string, int64_t> inter_ref_and_index;

    };
}
#endif
