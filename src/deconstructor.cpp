#include "deconstructor.hpp"

using namespace std;


namespace vg {
	Deconstructor::Deconstructor(){

	}

	Deconstructor::Deconstructor(xg::XG x) : my_xg(x){

	}

	Deconstructor::Deconstructor(VG v) : my_vg(v){
		init();
	}

	Deconstructor::~Deconstructor(){
	}

	void Deconstructor::init(){
		my_vg.topological_sort(nodes_in_topo_order);
		ord_ID = nt_to_ids(nodes_in_topo_order);
		previous_entrances = vector<int64_t> (ord_ID.size(), -1);
		alt_entrances = vector<int64_t> (ord_ID.size(), -1);
	}

	void Deconstructor::report(SuperBubble sb){
		my_super_bubbles.push_back(sb);
		cerr << sb.start_node << "\t" << sb.end_node << endl;
	}

	/**
	Use a basic BFS to find nodes on path
	deconstruct(path p){

}
*/

/**
* Another function taken from Brankovic. Takes a potential start and
* end index for a superbubble.
* Returns: a SuperBubble struct if one is found. An empty one if
* no superbubbles are found.
* Modifies global variables!
* Should probably have a matching report function which actually returns a bloody superbubble
* and tosses it in the global namespace.
*/
SuperBubble Deconstructor::report_superbubble(int64_t start, int64_t end){
	// This seems iffy, check it.
	if ((start == -1 || end == -1) || ord_ID[start] >= ord_ID[end]){
		delete_tail(candidates);
	}
	int64_t s = previous_entrances[end];
	int64_t valid;
	while (ord_ID[s] >= ord_ID[start]){
		valid = validate_superbubble(s, end);
		//TODO check this
		if ((s == valid || valid == alt_entrances[s]) || valid == -1){
			break;
		}
		alt_entrances[s] = valid;
		s = valid;
	}
	delete_tail(candidates);
	if (valid == s){
		//TODO MAke a superbubble obj and toss it up!
		SuperBubble sb;
		report(sb);
		//March back from end to s in candidates
		while (tail(candidates) != s){
			if(is_exit(tail(candidates))){
				report_superbubble(next(s), tail(candidates));
			}
			else{
				delete_tail(candidates);
			}
		}
	}


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
	int prev_entr = -1;

	for (int ind = 0; ind < nodes_in_topo_order.size(); ind++){
		int64_t vert_id = nodes_in_topo_order[ind].node->id();
		alt_entrances[vert_id] = -1;
		previous_entrances[vert_id] = prev_entr;
		if (is_exit(vert_id)){
			insert_exit(vert_id);
		}
		if (is_entrance(vert_id)){
			insert_entrance(vert_id);
			prev_entr = vert_id;
		}
	}

	while(!candidates.empty()){
		if (is_entrance( tail(candidates) )){
			delete_tail( candidates );
		}
		else{
			// Could use pointers, but we use indices instead.
			// Perhaps wiser to use iterators?? TODO
			report_superbubble(head(candidates), tail(candidates));
		}
	}


	return ret;
}

/**
* Take in the indices of potential superbubble entrances/exits
* and return the index of the actual entrance OR -1 if none found
* and superbubble isn't so super.
* Modifies global variables!
*/
int Deconstructor::validate_superbubble(int start, int end){
	//NodeTraversal s = nodes_in_topo_order[start];
	//NodeTraversal e = nodes_in_topo_order[end];
	int outchild = rangemax(ord_ID, start, end-1);
	int outparent = rangemin(ord_ID, start+1, end);
	if (outchild != end){
		return -1;
	}
	if (outparent == start){
		return start;
	}
	else if (is_entrance(outparent)){
		return outparent;
	}
	else{
		return previous_entrances[outparent];
	}
}

vector<int64_t> Deconstructor::emit_nodes_on_path_through_superbubble(Path p, vector<SuperBubble> subs){

}

vector<int64_t> Deconstructor::emit_nodes_off_path_through_superbubble(SuperBubble sb){

}

vector<vcflib::Variant> Deconstructor::sb_to_variants(SuperBubble sb){

}

/** Checks whether a vertex V satisifes Lemma 3,
* that v has one parent vertex.
*/
bool Deconstructor::is_exit(int64_t i){
	Node v = vertex(i);
	if (v.backward){
		vector<pair<id_t, bool>>& edges_end = my_vg.edges_end(v.node);
		return edges_end.size() <= 1;
	}
	else{
		vector<pair<id_t, bool>>& edges_start = my_vg.edges_start(v.node);
		return edges_start.size() <= 1;
	}
}

bool Deconstructor::is_entrance(int64_t i){

	Node v = vertex(i);
	if (!v.backward){
		vector<pair<id_t, bool>>& edges_end = my_vg.edges_end(v.node);
		return edges_end.size() <= 1;
	}
	else{
		vector<pair<id_t, bool>>& edges_start = my_vg.edges_start(v.node);
		return edges_start.size() <= 1;
	}
}

int64_t Deconstructor::rangemin(deque<int64_t> a, int i, int j){
	int64_t min = INT64_MAX;
	for (int ind = i; ind < j; ind++){
		if (a[ind] < min){
			min = a[ind];
		}
	}
	return min;
}

int next(int n){
	return n + 1;
}

int64_t Deconstructor::rangemax(deque<int64_t> a, int i, int j){
	int64_t max = -1;
	for (int ind = i; ind < j; ind++){
		if (max < a[ind]){
			max = a[ind];
		}
	}
	return max;
}

void Deconstructor::insert_exit(int64_t n){
	candidates.push_back(n);
}
void Deconstructor::insert_entrance( int64_t n){
	candidates.push_back(n);
}

int Deconstructor::head(deque<int>& cands){
	return cands.front();
}

int Deconstructor::tail(deque<int64_t>& cands){
	return cands.back();
}

void Deconstructor::delete_tail(deque<int64_t>& cands){
	cands.pop_back();
}

Node Deconstructor::vertex(int64_t i){
	return *(my_vg.get_node(i));
}

vector<int64_t> Deconstructor::nt_to_ids(deque<NodeTraversal> nt){
	deque<int64_t> ret = vector<int64_t>(nt.size(), 0);
	int count = 0;
	for (auto n: nt){
		ret[n.node->id()] = count;
		count++;
	}
	return ret;
}

}
