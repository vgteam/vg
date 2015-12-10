#ifndef DECON_HPP
#define DECON_HPP
#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <unordered_map>
#include <queue>
#include <stack>
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
            void set_index(string ref_index);
            void set_graph(VG* v);
            void b_call(string pathname);
            void indel_caller(string pathname);

            /**
             * Project a path onto another path,
             * much like a projection of p2 onto p1 in space.
             */
             //list<Mapping> relative_mapping(Path& p1, Path& p2);

            void get_variants_using_edges_from_file(string pathfile);
            void get_variants_using_edges(string pathname);
           /**
             * Build a vcf record from a mapping.
             */
            vcflib::Variant mapping_to_simple_variant(Mapping m, int64_t alt_id);


            Mapping node_id_to_mapping(int64_t node_id);
            /**
             * Turn a vector of variants into a proper VCF.
             */
            void write_variants(string filename, vector<vcflib::Variant> variants);

        private:
            // TODO Should be able to handle XG or VG indices
            string ref_index;
            string reference;
            string xg_file;
            VG* vgraph;
            vector<string> ref_paths;
            map<string, int64_t> inter_ref_and_index;
            int in_degree(Node n);
            bool on_ref(Node n, int64_t path);
            Node get_anchor_node(Node current, int64_t path);
            int inDegree(Node n);
            bool beenVisited(Node n, map<int64_t, int> node_to_level);
            Node get_alleles(Node n, int64_t path_id, map<int64_t, int>& node_to_level);
            pair<Node, Mapping> map_between_nodes(Node a, Node b);
            bool on_ref(int64_t node_id, int64_t path);

    };
}
#endif
