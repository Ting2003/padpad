// ----------------------------------------------------------------//
// Filename : circuit.h
// Author : Xiao Zigang <zxiao2@illinois.edu>
//
// declaration of Circuit class
// used to construct the circuit network
// ----------------------------------------------------------------//
// - Ting Yu - Tue Feb 8 5:45 pm 2011
//   * added the ostream<< func in .h file
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#ifndef __CIRCUIT_H__
#define __CIRCUIT_H__

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <tr1/unordered_map>
#include <map>
#include <queue>
#include <list>
#include <cmath>
#include "cholmod.h"
#include "global.h"
#include "node.h"
#include "net.h"
#include "vec.h"
#include "triplet.h"
#include "pad.h"
#include "block.h"
using namespace std;
using namespace std::tr1;

typedef vector<double> DoubleVector;
typedef vector<Net *> NetPtrVector;
typedef list<Net *> NetPtrList;
typedef vector<Node *> NodePtrVector;
typedef NetPtrVector NetList;

// functor of translating Node * to void *
namespace std{ 
	namespace tr1{
		template<> struct hash< Node * >{
			size_t operator()( const Node * x ) const {
				return hash< const char* >()( (char*) x );
			}
		};
	}
}

class Circuit{
public:
	Circuit(string name="");
	~Circuit();
	void check_sys() const;
	friend class Block;
	// can be written as inline to speed up
	Node * get_node(string name);
	Net * get_net(string name);
	string get_name() const;

	static size_t get_total_num_layer();

	// add a node into nodelist
	bool add_node(Node * nd);

	// add a net into netset
	bool add_net(Net * net);

	bool has_node(string name) const;
	bool has_net(string name) const;

	// sort nodes according to predefined order
	void sort_nodes();

	// solve for node voltage
	void solve();
	
	void set_blocklist(Node * nd);

	static void set_parameters(double, double, double, size_t, int);
	static void get_parameters(double&, double&, double&, size_t&, int&);

	friend ostream & operator << (ostream & os, const Circuit & ckt);
	friend class Parser;

	////// new functions for pad /////
	double locate_maxIRdrop();
	double locate_special_maxIRdrop();
	void mark_special_nodes();
	void build_pad_set();
	////// new member for pad //////
	
	double max_IRdrop;
	vector<Pad*> pad_set;
	vector<Node*> origin_pad_set;
	vector<Node*> special_nodes;
	// mapping from name to Node object pointer
	unordered_map<string, Node*> map_node_pt;

	////// new functions for pad /////
	void assign_distance(Node *nds, Node *nd, double dist);
	void print_pad_map();
	void clear_flags();
	double update_pad_pos(double ref_drop_value, size_t i);
	double update_pad_pos_all(vector<double> ref_drop_vec);
	void round_data(double &data);
	Node * pad_projection(Pad *pad, Node *nd);
	void project_pads();
	bool has_node_pt(string pt_name) const;
	Node * get_node_pt(string pt_name);
	void build_map_node_pt();
	void relocate_pads();
	void relocate_pads_graph();
	void restore_pad_set(vector<Node*>&pad_set_old);
	void assign_pad_set(vector<Node*>&pad_set_old);
	void rebuild_voltage_nets();
	void print_pad_set();
	void extract_pads(int pad_number);
	void print_matlab();
	void clear_pad_control_nodes();
	void update_pad_control_nodes(vector<double> & ref_drop_vec, size_t iter);
	void extract_min_max_pads(vector<double> ref_drop_vec);
	void extract_min_max_pads_new(vector<double> ref_drop_vec);

	void build_graph();
	Pad *find_nbr_pad(Pad *pad);
	double get_distance(Node *na, Node *nb);
	void graph_move_pads(vector<double> ref_drop_vec);
	int locate_max_drop_pad(vector<double> vec);
	double calc_avg_ref_drop(vector<double> &ref_drop_vec);
	double calc_avg_ref(vector<double> ref_drop_vec);
	double locate_ref(size_t i);
	void dynamic_update_violate_ref(vector<double> & ref_drop_vec);
	bool print_flag(Node *nd);
	void move_violate_pads(vector<double> ref_drop_vec);
	void modify_newxy();
	void resolve_direct();
	void resolve_queue(vector<Node *> origin_pad_set);
	void solve_queue(vector<Node *> pad_set_old);
	void initialize_queue(vector<Node *> pad_set_old, queue <Node*> &q);
	double update_single_iter(vector<Node *> pad_set_old);
	double update_value(Node *nd);
	size_t update_queue(queue<Node *>&q, Node *nd);
	void solve_GS();
	void print_all_control_nodes();
	//////// end functions for pad ////////

	// C style output
	void print();
	cholmod_common c, *cm;
	size_t peak_mem;
	size_t CK_mem;
private:
	// member functions
	void solve_LU();
	void solve_LU_core();

	bool solve_IT();
	void solve_block_LU();

	bool solve_pcg();
	//bool solve_block_pcg();

public:	
	// initialize things before solve_iteration
	void solve_init();
private:
	// updates nodes value in each iteration
	double solve_iteration();
	void block_init();
	void update_block_geometry();

	// methods of stamping the matrix
	void stamp_by_set(Matrix & A, double* b);
	void stamp_resistor(Matrix & A, Net * net);
	void stamp_current(double* b, Net * net);
	void stamp_VDD(Matrix & A, double* b, Net * net);
	
	void make_A_symmetric(double *bp);
	void make_A_symmetric_block();

	void stamp_block_matrix();
	void stamp_boundary_matrix();
	void stamp_boundary_net(Net * net);
	void stamp_block_resistor(Net *net, Matrix * A);
	void stamp_block_current(Net * net);
	void stamp_block_VDD(Net * net, Matrix * A);

	void update_block_rhs(Block & block, int dir);

	//  ******* method for PCG method  ********
	// solve circuit with preconditioned pcg method
	void copy_node_voltages_block(bool from=true);

	// after solving, copy node voltage from replist to nodes
	void get_voltages_from_LU_sol(double* x);
	void get_voltages_from_block_LU_sol();
	void get_vol_mergelist();

	Vec compute_precondition(const Matrix & ML, const Matrix & D, 
			const Matrix & MU, Vec &r);
	void init_precondition(const Matrix &A, Matrix &ML, Matrix &D,
			Matrix &MU);
	// searching all the blocks for global residue
	Vec get_block_res(const Vec& b, const Vec& xk1);
	// searching all the blocks for global zk1
	Vec compute_block_precondition( Vec &r);

	void set_len_per_block();
	void find_block_size ();
	void block_boundary_insert_net(Net * net);
	void find_block_base();

	void partition_circuit();
	double modify_voltage(Block & block, double* x_old);

	void node_voltage_init();
	void solve_one_block(size_t block_id);

	void select_omega();

	void set_type(CIRCUIT_TYPE type){circuit_type = type;};

	void get_samples();

	bool check_diverge() const;

	void merge_along_dir(Node *, DIRECTION dir);
	Node * merge_along_dir_one_pass(Node *, DIRECTION dir, bool remove);
	void merge_node(Node * node);

	// ************** member variables *******************
	NodePtrVector nodelist;		// a set of nodes
	NodePtrVector replist;		// a set of representative nodes
	NodePtrVector mergelist;	// nodes for merging
	NetList net_set[NUM_NET_TYPE];// should be the same as size of NET_TYPE
	// defines the net direction in layers
	static vector<LAYER_DIR> layer_dir;
	vector<int> layers;
	
	// mapping from name to Node object pointer
	unordered_map<string, Node*> map_node;

	// mapping from Net pointer to their index in netlist
	unordered_map<Net*, size_t> net_id;

	// circuit name
	string name;

	// blocks
	BlockInfo block_info;
	size_t x_min, y_min, x_max, y_max;

	// control variables
	static double EPSILON;
	static double OMEGA;
	static double OVERLAP_RATIO;
	static size_t MAX_BLOCK_NODES;
	static int MODE; // 0 = IT, 1 = LU

	CIRCUIT_TYPE circuit_type;

	NodePtrVector sample;

	double VDD;
};

inline size_t Circuit::get_total_num_layer(){return layer_dir.size();}

// adds a node into nodelist
inline bool Circuit::add_node(Node * node){
	nodelist.push_back(node);
	map_node[node->name] = node;
	return true;
}

// adds a net into netset
inline bool Circuit::add_net(Net * net){
	if( net->type == RESISTOR )
		net_id[net] = net_set[net->type].size();
	net_set[net->type].push_back(net);
	return true;
}

// fina a node by name
inline bool Circuit::has_node(string name) const{
	if( map_node.find(name) != map_node.end() ) return true;
	return false;
}

// get a node by name
inline Node * Circuit::get_node(string name){
	unordered_map<string, Node*>::const_iterator it = map_node.find(name);
	if( it != map_node.end() ) return it->second;
	else return NULL;
}

// fina a node by pt
inline bool Circuit::has_node_pt(string pt_name) const{
	if( map_node_pt.find(pt_name) != map_node_pt.end() ) return true;
	return false;
}

// get a node by pt
inline Node * Circuit::get_node_pt(string pt_name){
	unordered_map<string, Node*>::const_iterator it = map_node_pt.find(pt_name);
	if( it != map_node_pt.end() ) return it->second;
	else return NULL;
}

inline void Circuit::merge_node(Node * node){
	for(DIRECTION dir = WEST; dir <= NORTH; dir=DIRECTION(dir+1)){
		// test whether this line has been processed
		if( node->end[dir] != node ) continue;

		// probe for one step, if the line is only one step, don't do it.
		Node * next = node->get_nbr_node(dir);
		if( next == NULL || !next->is_mergeable() ) continue;
		merge_along_dir(node, dir);
	}
}

/*
// find a net by name
inline bool Circuit::has_net(string name) const{
	if( map_net.find(name) != map_net.end() ) return true;
	return false;
}


// get a net by name
inline Net * Circuit::get_net(string name){return map_net[name];}
*/

bool compare_node_ptr(const Node *a, const Node *b);
bool compare_pads(const pair<Node*, double> a, const 
	pair<Node*, double> b);
bool pad_equal(Pad *pa, Pad *pb);
ostream & operator << (ostream & os, const NodePtrVector & nodelist);
ostream & operator << (ostream & os, const NetList & nets);
//ostream & operator << (ostream & os, const vector<Block > & block_info);
#endif
