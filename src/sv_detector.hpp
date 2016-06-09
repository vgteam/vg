#include <vector>
#include <map>
#include <unordered_map>
#include "vg.hpp"
#include "path.hpp"
#include "vg.pb.h"
#include "filter.hpp"
#include "stream.hpp"
#include <regex>

using namespace std;
using namespace vg;

struct SV_Position{
    pair<Position, Position> start_and_end;
    string start_sequence_source;
    long start_sequence_position;
};

struct Breakpoint{
    vg::id_t node_id;
    int offset;
    bool approximate;
    int bound = 0;
    string source_sequence = "";
    bool novel = false;
};

struct Evidence{
    vector<Alignment> reads;
    vector<pair< string, int> > ref_and_alt_counts;

};


struct StructuralVariant{
    pair <SV_Position, SV_Position> positions;
    pair <Breakpoint, Breakpoint> breakpoints;
    bool novel;
    string path;
    string to_vcf(){
        return "";
    }
};

struct Deletion : StructuralVariant {
    string sequence = "";
};

struct Insertion{

};

struct Inversion{

};

struct Translocation{

};

struct WigglePos{
    
};


struct SingleNucleotideVariant{
    vector<string> ref;
    vector<string> alt;
    int position;
    string to_vcf(){
        return "";
    }
};

class SV_DETECTOR{
    
    public:
        SV_DETECTOR();
        SV_DETECTOR(VG* g);
        ~SV_DETECTOR();

        void cache_paths_of_graph();
        /** Returns a vector of StructuralVariant objects, which can be turned into VCFs **/
        vector<StructuralVariant> gam_to_known_sv(string gamfile);

        /** Returns a list of alt_path names that the alignment matches to **/
        vector<string> alignment_to_known_sv(Alignment aln);

        /** Looks for signatures based on Filter and returns StructuralVariant objects **/
        vector<StructuralVariant> alignment_to_putative_sv(Alignment aln);

        //vector<Deletion> alignment_to_deletion(Alignment& aln);

    private:
        Filter read_filter;
        vg::VG* my_graph;
        map<string, list<Mapping> > name_to_paths;
    
        map<string, pair<int, int> > known_to_ref_alt_count;
        map<vg::id_t, string> node_id_to_path;
        map<vg::id_t, string> node_to_variant;

        // Allow <wiggle> basepairs of variance at the tips of SVs
        /*
        int wiggle = 10;
        map<int64_t, StructuralVariant> known_node_vg::id_to_sv ;
        map<int64_t, SingleNucleotideVariant> known_node_vg::id_to_snp;
        map<int64_t, StructuralVariant> unknown_node_vg::id_to_sv;
        map<Position, Deletion> pos_to_del;

        map<WigglePos, vector<Position> > w_to_p;

        unordered_map<Deletion, pair<int, int> > del_to_refCount_altCount;

        void gam_sr(Alignment& aln);

        void alignment_split_read(Alignment& aln);
        void alignment_read_pair(Alignment& aln);
        void alignment_known(Alignment& aln);
        void calculate_avg_depth(vg::VG graph); 
        */
};
