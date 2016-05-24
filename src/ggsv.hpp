#include <vector>
#include <map>
#include <unordered_map>

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

struct Deletion{

};

struct Insertion{

};

struct Inversion{

};

struct Translocation{

};

struct Breakpoint{
    id_t node_id;
    int offset;
    bool approximate;
    int bound = 0;
    string source_sequence = "";
    bool novel = true;
};

struct SingleNucleotideVariant{
    string ref;
    string alt;
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

    private:
        map<int64_t, StructuralVariant> known_node_id_to_sv ;
        map<int64_t, SingleNucleotideVariant> known_node_id_to_snp;
        map<int64_t, StructuralVariant> unknown_node_id_to_sv;

        unordered_map<StructuralVariant, pair<int, int> > known_sv_to_refCount_altCount;
}
