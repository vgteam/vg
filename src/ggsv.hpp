#include <vector>
#include <map>
#include <unordered_map>
#include "vg.hpp"
#include "filter.hpp"


struct StructuralVariant{
    string type;
    id_t start_id
    int start;
    id_t end_id
    int end;
    string path;
    string to_vcf(){
        return "";
    }
};

struct SV_Position{
    pair<Position, Position> start_and_end;
    string start_sequence_source;
    long start_sequence_position;
};

struct Breakpoint{
    id_t node_id;
    int offset;
    bool approximate;
    int bound = 0;
    string source_sequence = "";
    bool novel = false;
};

struct Evidence{
    vector<Alignment> reads;

};

struct Deletion{
    SV_Position position;
    pair<Breakpoint, Breakpoint> breakpoints;
    bool novel = false;
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

class GGSV{
    
    public:
        GGSV();
        GGSV(string indexfile);
        ~GGSV();
        void set_index(indexfile);
        vector<StructuralVariant> alignment_to_known_sv(Alignment aln);
        vector<StructuralVariant> alignment_to_putative_sv(Alignment aln);
        vector<Deletion> alignment_to_deletion(Alignment& aln);

    private:
        Filter read_filter;
        // Allow <wiggle> basepairs of variance at the tips of SVs
        int wiggle = 10;
        map<int64_t, StructuralVariant> known_node_id_to_sv ;
        map<int64_t, SingleNucleotideVariant> known_node_id_to_snp;
        map<int64_t, StructuralVariant> unknown_node_id_to_sv;
        map<Position, Deletion> pos_to_del;

        map<WigglePos, vector<Position> > w_to_p;

        unordered_map<Deletion, pair<int, int> > del_to_refCount_altCount;

        void gam_sr(Alignment& aln);

        void alignment_split_read(Alignment& aln);
        void alignment_read_pair(Alignment& aln);
        void alignment_known(Alignment& aln);
        void calculate_avg_depth(vg::VG graph);
}
