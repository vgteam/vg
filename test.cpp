#include <iostream>
#include <fstream>
#include "pb2json.h"
#include "vg.pb.h"

using namespace std;
using namespace google::protobuf;
using namespace vg;

int main(int argc,char *argv[])
{
    Node n;
    n.set_sequence("AAAACACGAAGGTGCTCACTGAAACATGGGAACCAAAGTTTCCCACAACATAAGGAGCAGAGTGAAACTGCAGAGGTACAATCATGGAATTCCAGAAAATG");
    n.set_name("14:20000000-20000100");
    n.set_id(1);
    n.add_prev(2);
    n.add_prev(3);
    for (int i = 0; i < n.prev_size(); ++i) {
        cout << n.prev(i) << endl;
    }

    Graph g;
    Node* m = g.add_nodes();
    *m = n;
    m = g.add_nodes();
    *m = n;
    m = g.add_nodes();
    *m = n;

    ofstream of("test.vg");
    n.SerializeToOstream(&of);
	char *json2 = pb2json(n);
	cout<<json2<<endl;
	free(json2);

	char *json3 = pb2json(g);
	cout<<json3<<endl;
	free(json3);
    return 0;
}
