#include "deconstructor.hpp"

using namespace std;


namespace vg {
	Deconstructor::Deconstructor(){

	}

	Deconstructor::Deconstructor(xg::XG x) : my_xg(x){

	}

	Deconstructor::Deconstructor(VG v) {
		my_vg = v;;
  	
		init();
	}

	Deconstructor::~Deconstructor(){
	}

	void Deconstructor::init(){
	}

/**
*   BFS through a superbubble and fill out the corresponding SuperBubble
*   struct.
*/
SuperBubble Deconstructor::report_superbubble(int64_t start, int64_t end){
	
	SuperBubble sb;
	return sb;
}

/**
* Implementation of method SuperBubble in Brankovic 2015.
* Takes a graph (vg), topologically sorts it, and marches
* backward through the ordered nodes looking for candidates superbubbles.
* Makes a SuperBubble struct if one is found, and appends it to a vector which
* is returned. If none are found, the returned vector is empty.
* Modifies global variables!
* TODO need to call unroll/unfold/dagify first
* TODO should really use more pointers.
*/
vector<SuperBubble> Deconstructor::get_all_superbubbles(){
	vector<SuperBubble> ret;;



	return ret;
}


vector<int64_t> Deconstructor::nt_to_ids(deque<NodeTraversal>& nt){
	vector<int64_t> ret = vector<int64_t>(nt.size(), 0);
	int64_t count = 0;
	for (auto n: nt){
		ret[(n.node->id() - 1)] = count;
		count++;
	}

	return ret;
}

}
