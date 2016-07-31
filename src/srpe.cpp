#include "srpe.hpp"

using namespace std;
namespace vg{

    string SRPE::bucket(Position p){
        int offset = p.offset();
        bucket pos = offset  - (offset % my_bucket_size);
        vg::id_t node_id = p.node_id();
        stringstream sst;
        sst << node_id << "_" << offset;
        return sst.str();
        
    }
    void SRPE::apply(Alignment& aln, vector<Alignment>& suspect){

    }

    void SRPE::apply(Alignment& aln, vector<string>& suspect){
        

    }

    void SRPE::pre_depth(Alignment& aln){

    }

    void SRPE::call_depth_variants(void){

    }

    void SRPE::call(void){

    }

    string SRPE::pos_to_string(Positoon& p){
        string stream sst;
        sst << p.node_id() << "_" << p.offset();
        return sst.str();

    }

}

/*      
        private:
        Filter my_filter;
        map<Position, int> depth_at;
        int avg_depth;
// int depth_window;
// Consider variants within my_pos_window_size of each other on either
// side (WITHIN THE SAME NODE) the same variant.
int my_pos_window_size = 20;


*/
