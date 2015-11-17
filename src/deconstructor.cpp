#include "deconstructor.hpp"

using namespace std;

namespace vg {
    Deconstructor::Deconstructor(){
        this.index = nullptr;
    }

    Deconstructor::Deconstructor(Index i) index(i){

    }

    Deconstructor::~Deconstructor(){
        clear()
    }

    void clear(){
        
    }

    /**
     * Set the deconstructor's index, be it a rocksdb-backed disk index
     * or an XG + gcsa index
     */
    void Deconstructor::set_index(Index index, bool isXG){
        if (isXG){
        }
        else{
        }
    }

    /**
     * This one is still a bit unsettled in implementation, but the gist
     * is to:
     * 1. Open the index
     * 2. 
     */
    Path Deconstructor::relative_mapping(Path& p1, Path& p2){
        VG graph;
        Path ret;
        
        this.index.open_read_only(db_name)
    }

    Path Deconstructor::relative_mapping(Path& p1, string p2){
        VG graph;
        Path ret;
    }

}
