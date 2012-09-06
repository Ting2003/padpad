// ----------------------------------------------------------------//
// Filename : circuit.cpp
// Author : Zigang Xiao <zxiao2@illinois.edu>
//          Ting Yu <tingyu1@illinois.edu>
//
// implementation file of circuit.h
// ----------------------------------------------------------------//
// - Zigang Xiao - Sun Jan 30 18:35:56 CST 2011
//   * Add UMFPACK support
// - Zigang Xiao - Tue Jan 25 17:19:21 CST 2011
//   * added framework of PCG
// - Zigang Xiao - Tue Jan 18 21:46:13 CST 2011
//   * added solve() and related function
// - Zigang Xiao - Sun Jan 16 16:04:03 CST 2011
//   * added this log

#include <fstream>
#include <string>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <utility>
#include <cassert>
#include <vector>
#include "cholmod.h"
#include "umfpack.h"
#include "circuit.h"
#include "util.h"
#include "algebra.h"
#include "node.h"
#include "pad.h"
using namespace std;

double Circuit::EPSILON = 1e-5;
size_t Circuit::MAX_BLOCK_NODES =5500;
double Circuit::OMEGA = 1.0;//1.2;
double Circuit::OVERLAP_RATIO = 0.2;
int    Circuit::MODE = 0;
const int MAX_ITERATION = 1000;//1000000;//1000;
const int SAMPLE_INTERVAL = 5;
const size_t SAMPLE_NUM_NODE = 10;
const double MERGE_RATIO = 0.3;
//////////////////////////////////////////////////////////////////////////
// Constructor and utility functions goes here

vector<LAYER_DIR> Circuit::layer_dir(MAX_LAYER);

// constructor of Circuit class, name is optional
Circuit::Circuit(string _name): 
	max_IRdrop(0), name(_name),
	x_min(INFTY),y_min(INFTY),x_max(0),y_max(0),
	circuit_type(UNKNOWN),VDD(0.0) {
	// add ground node
	Node * gnd = new Node(string("0"), Point(-1,-1,-1));
	gnd->rep = gnd;
	this->add_node(gnd);

	for(int i=0;i<MAX_LAYER;i++)
		layer_dir[i]=NA;
	peak_mem = 0;
	CK_mem = 0;
}

// Trick: do not release memory to increase runtime
Circuit::~Circuit(){
	pad_set.clear();
	special_nodes.clear();
	map_node_pt.clear();
	for(size_t i=0;i<nodelist.size();i++) delete nodelist[i];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		for(it=ns.begin();it!=ns.end();++it)
			delete *it;
	}
}

void Circuit::check_sys() const{
	clog<<"**** CHECKING SYSTEM ENVIRONMENT ****"<<endl;
	clog<<"* int size     = "<< sizeof(int)<<endl;
	clog<<"* long size    = "<< sizeof(long)<<endl;
	clog<<"* size_t size  = "<< sizeof(size_t)<<endl;
	clog<<"* UF_long size = "<< sizeof(UF_long)<<endl;
	clog<<"* Max nodelist = "<<(size_t)nodelist.max_size()<<endl;
	clog<<"****            END              ****"<<endl<<endl;
}

// functor to be used in STL sort function
// order: y > x > z > flag 
// input: two node a, b
// return true if a < b, false o/w
// note that ground note are put to last
bool compare_node_ptr(const Node * a, const Node * b){
	if( a->is_ground() ) return false;
	if (b->is_ground() ) return true;

	if( a->pt.y == b->pt.y ){
		if( a->pt.x == b->pt.x ){
			if( a->pt.z == b->pt.z ){
				return a->isX();
			}
			else{
				return (a->pt.z > b->pt.z);// top down
			}
		}
		else
			return ( a->pt.x < b->pt.x );
	}
	else
		return (a->pt.y < b->pt.y);
}

bool compare_values(double a, double b){
	return (a<b);
}

// sort the nodes according to their coordinate 
void Circuit::sort_nodes(){
	sort(nodelist.begin(), nodelist.end(), compare_node_ptr);
	// update node id mapping, 
	// NOTE: ground node will be the last
}

string Circuit::get_name() const{return this->name;}

ostream & operator << (ostream & os, const NodePtrVector & nodelist){
	for(size_t i=0;i<nodelist.size();i++)
		os<<*nodelist[i]<<endl;
	return os;
}

ostream & operator << (ostream & os, const NetList & nets){
	NetList::const_iterator it;
	for(it=nets.begin();it!=nets.end();++it)
		if( (*it) != NULL ) os<<**it<<endl;
	return os;
}

ostream & operator << (ostream & os, const Circuit & ckt){
	os<<"Circuit ["<<ckt.name<<"] info:"<<endl;

	os<<"==== Nodes ===="<<endl;
	os<<ckt.nodelist;

	os<<"==== Reps  ===="<<endl;
	os<<ckt.replist;

	os<<"==== Nets  ===="<<endl;
	os<<ckt.net_set[RESISTOR];

	return os;
}

void Circuit::print(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%s  %.5e\n", nodelist[i]->name.c_str(), 
				nodelist[i]->value);
	}
}

///////////////////////////////////////////////////////////////////////////////
// Computation Functions

// initialization before solving the circuit
// 1. sort the nodes
// 2. set node representatives
// 3. find node in which block, update count
// 4. get representative lists
void Circuit::solve_init(){
	sort_nodes();

	size_t size = nodelist.size() - 1;
	Node * p = NULL;
	for(size_t i=0, nr=0;i<size;i++){
		p=nodelist[i];

		// test if it can be merged
		if( p->is_mergeable() ){
			mergelist.push_back(p);
			continue;
		}

		Net * net = p->nbr[TOP];
		merge_node(p);

		// find the VDD value
		if( p->isX() ) VDD = p->get_value();

		// test short circuit
		if( !p->isX() && // X must be representative 
		    net != NULL &&
		    fzero(net->value) ){
			// TODO: ensure ab[1] is not p itself
			assert( net->ab[1] != p );
			p->rep = net->ab[1]->rep;
		} // else the representative is itself

		// push the representatives into list
		if( p->rep == p ) {
			replist.push_back(p);
			p->rid = nr++;
		}
	}// end of for i

	//size_t n_merge = mergelist.size();
	//size_t n_nodes = nodelist.size();
	//size_t n_reps  = replist.size();
	//double ratio = n_merge / (double) (n_merge + n_reps);
	/*clog<<"mergeable  "<<n_merge<<endl;
	clog<<"replist    "<<n_reps <<endl;
	clog<<"nodelist   "<<n_nodes<<endl;
	clog<<"ratio =    "<<ratio  <<endl;*/

	net_id.clear();
	build_pad_set();
}

void Circuit::mark_special_nodes(){
	special_nodes.clear();
	int type = CURRENT;
	Net *net;
	Node *nd;
	// push back all current nodes
	for(size_t i=0;i<net_set[type].size();i++){
		net = net_set[type][i];
		nd = net->ab[0];
		if(nd->is_ground())
			nd = net->ab[1];
		//if(nd->pt.x%10==0 &&
			//nd->pt.y%10==0){// i%4000==0){
			special_nodes.push_back(nd);
			//clog<<"special_nodes: "<<*nd<<endl;
		//}
	}
	//clog<<"total special nodes: "<<special_nodes.size()<<endl;
	/*for(size_t i=0;i<nodelist.size()-1;i++){
		if(nodelist[i]->name == "n0_0_0" ||
		   nodelist[i]->name =="n0_1_2")//"n0_150_100")
			special_nodes.push_back(nodelist[i]);
	}*/
}

// partition the circuit to X_BLOCKS * Y_BLOCKS blocks
// according to the node size. Assuming they are distributed
// uniformly at random
void Circuit::partition_circuit(){
	size_t num_nodes = replist.size();
	size_t num_blocks =  num_nodes / MAX_BLOCK_NODES;
	if( num_nodes % MAX_BLOCK_NODES > 0 ) ++num_blocks;
	size_t len_x = x_max-x_min;
	size_t len_y = y_max-y_min;
	size_t X_BLOCKS, Y_BLOCKS;

	// Extreme case: only one node in x/y axis
	if( num_blocks == 1 ){
		X_BLOCKS = Y_BLOCKS = 1;
	}
	else if( len_x == 0 ){
		X_BLOCKS = 1;
		Y_BLOCKS = num_blocks;
	}
	else if (len_y == 0 ){
		Y_BLOCKS = 1;
		X_BLOCKS = num_blocks;
	}
	else{// compute length of x and y according to their ratio
		double ratio = double(len_x) / double(len_y);
		X_BLOCKS = sqrt(num_blocks/ratio);
		Y_BLOCKS = X_BLOCKS*ratio;
		if( X_BLOCKS * Y_BLOCKS < num_blocks ) Y_BLOCKS+=1; 
		num_blocks = X_BLOCKS * Y_BLOCKS;
	}
	//X_BLOCKS =2;
	//Y_BLOCKS =1;
	clog<<"num_blocks: "<<X_BLOCKS<<" / "<<Y_BLOCKS <<endl;
	block_info.X_BLOCKS = X_BLOCKS;
	block_info.Y_BLOCKS = Y_BLOCKS;
	block_info.resize(X_BLOCKS * Y_BLOCKS);
}

// build up block info
// 1. Find block divide point
// 2. Set each node in replist into block
// 3. Compute block size
// 4. Insert boundary netlist into map
void Circuit::block_init(){
	block_info.set_len_per_block(x_min, x_max, y_min, y_max, OVERLAP_RATIO);
	block_info.update_block_geometry();
	find_block_size();
	copy_node_voltages_block();

	//stamp_boundary_matrix();
	stamp_block_matrix();
}

void Circuit::block_boundary_insert_net(Net * net){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
	if(nd[0]->is_ground() || nd[1]->is_ground()) return;

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};

	// idea: for each block where node_a is in
	// find whether node_b is also in the block
	// if yes, then this net is not a boundary net
	// for this block
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		// for each block k in ls[j]
		for(size_t k=0;k<p->size();k++){
			// find whether block_id in another list
			size_t block_id = (*p)[k];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				Block & blk = block_info[block_id];
				blk.boundary_netlist.push_back(net);
			}
		}// end of for k
	}// end of for j
}

// stamp the nets by sets, block version
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_block_matrix(){
	size_t num_blocks = block_info.size();
	Matrix A[num_blocks];
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				Net * net = *it;
				if( net == NULL ) continue;
				if(net->ab[0]->is_ground() || 
				   net->ab[1]->is_ground()) 
					continue;
				assert( fzero(net->value) == false );
				stamp_block_resistor(*it, A);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it)
				stamp_block_current((*it));
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_block_VDD((*it), A);
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
	make_A_symmetric_block();
	clock_t t1, t2;
	t1 = clock();
	// after stamping, convert A to column compressed form
	for(size_t i=0;i<num_blocks;i++){
		if(block_info[i].count>0){
			A[i].set_row(block_info[i].count);		
			block_info[i].CK_decomp(A[i], cm, peak_mem, CK_mem);
		}
	}
	t2 = clock();
	clog<<"decomp time for CK is: "<<1.0*(t2-t1) / CLOCKS_PER_SEC<<endl;
	//clog<<"peak memory for cholmod: "<<peak_mem / 1e9<<" e+06"<<endl;
	//clog<<"CK_mem is: "<<CK_mem / 1e6 <<" G"<<endl;
}

// 1. mark rep nodes into corresponding blocks
// 2. find block size of replist nodes
// 3. allocate rhs size of each block
// 4. find local block index for each node
void Circuit::find_block_size(){
	// start assign replist node into blocks
	size_t n = replist.size();
	Node * p = NULL;
	
	// find which blocks the node belongs to,
	// and its index in the block
	for(size_t i=0;i<n;i++){
		p=replist[i];
		set_blocklist(p);
	}
	// for each block, allocate resource
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		block.allocate_resource(cm);
		//for(int k=0;k<4;k++) block.A_nbr[k].set_row(block.count);
	}
}

void Circuit::solve(){
	// solve_init();
	/*if( MODE == 0 )
		solve_IT();
	else*/
		solve_LU();
}

// solve Circuit
bool Circuit::solve_IT(){
	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	//solve_init();

	if( replist.size() <= 2*MAX_BLOCK_NODES ){
		clog<<"Replist is small, use direct LU instead."<<endl;
		solve_LU_core();
		return true;
	}
	
	//select_omega();
	partition_circuit();

	// Only one block, use direct LU instead
	if( block_info.size() == 1 ){
		clog<<"Block size = 1, use direct LU instead."<<endl;
		solve_LU_core();
		return true;
	}

	// cm declared in circuit class
	//cholmod_common c, *cm;
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;
	cm->final_ll = true;
	//cholmod_print_common("first_cm",cm);
	block_init();
	//cholmod_print_common("stamp_cm",cm);

	clog<<"e="<<EPSILON
	    <<"\to="<<OMEGA
	    <<"\tr="<<OVERLAP_RATIO
	    <<"\tb="<<MAX_BLOCK_NODES
	    <<"\tmode="<<MODE<<endl;

	int iter = 0;	
	double diff=0;
	bool successful = false;
	
	clock_t t1, t2;
	t1 = clock();
	while( iter < MAX_ITERATION ){
		diff = solve_iteration();
		iter++;
		clog<<"iter, diff: "<<iter<<" "<<diff<<endl;
		if( diff < EPSILON ){
			successful = true;
			break;
		}
	}
	t2 = clock();
	clog<<"solving iteration use: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	clog<<"# iter: "<<iter<<endl;
	get_voltages_from_block_LU_sol();
	get_vol_mergelist();
	//clog<<"before free. "<<endl;
	for(size_t i=0;i<block_info.size();i++){
		if(block_info[i].count > 0)
			block_info[i].free_block_cholmod(cm);
	}
	//clog<<"after free. "<<endl;
	
	cholmod_finish(cm);
	return successful;
}

// TODO: add comment
void Circuit::node_voltage_init(){
	for(size_t i=0;i<block_info.size();i++){
		Block & block = block_info[i];
		for(size_t j=0;j<block.count;j++){
			block.xp[i] = VDD;
			//block.x[i] = VDD;
			block.nodes[i]->value = VDD;
		}
	}
}

// One iteration during solving the circuit, for any block B:
// 1. update the righthand-side of the matrix of B
// 2. solve the matrix
// 3. update node voltages
// 4. track the maximum error of solution
double Circuit::solve_iteration(){	
	double diff = .0, max_diff = .0;
	for(size_t i=0;i<block_info.size();i++){
		Block &block = block_info[i];
		if( block.count == 0 ) continue;

		block.update_rhs();
		// backup the old voltage value
		
		double *x_old;
		x_old = new double [block.count];
		for(size_t k=0; k<block.count;k++){
			x_old[k] = block.xp[k];
		}
		//cout<<"Matrix A for block: "<<block.bid<<endl;
		block.solve_CK(cm);
		block.xp = static_cast<double *>(block.x_ck->x); 
		//cout<<"block_index: "<<block.bid<<endl;
		//for(size_t j=0;j<block_info[i].count;j++)
			//cout<<j<<" "<<block_info[i].bp[j]<<" "<<block_info[i].xp[j]<<endl;

		// modify node voltage with OMEGA and old voltage value
		diff = modify_voltage(block, x_old);
		delete [] x_old;

		//diff = distance_inf( block.x, x_old );
		if( max_diff < diff ) max_diff = diff;
	}
	return max_diff;
}

double Circuit::modify_voltage(Block & block, double * x_old){
	double max_diff = 0.0;
	OMEGA = 1.0;
	for(size_t i=0;i<block.count;i++){
		block.xp[i] = (1-OMEGA) * x_old[i] + OMEGA * block.xp[i];
		block.nodes[i]->value = block.xp[i];
		double diff = fabs(x_old[i] - block.xp[i]);
		if( diff > max_diff ) max_diff = diff;
	}
	return max_diff;
}

// stamp the matrix and solve
void Circuit::solve_LU_core(){
	size_t n = replist.size();	// replist dosn't contain ground node
	if( n == 0 ) return;		// No node
	//Vec b(n), x(n);
	cholmod_common c, *cm;
	cholmod_dense *b, *x;
	cm = &c;
	cholmod_start(cm);
	cm->print = 5;
	b = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	x = cholmod_zeros(n, 1, CHOLMOD_REAL, cm);
	double *bp, *xp;
	bp = static_cast<double *> (b->x);
	Matrix A;
	stamp_by_set(A, bp);
	make_A_symmetric(bp);
	//A.merge();

	A.set_row(replist.size());
	Algebra::solve_CK(A, x, b, cm, peak_mem, CK_mem);
	
	xp = static_cast<double *> (x->x);
	// Vec b contains result, copy it back to nodelist
	get_voltages_from_LU_sol(xp);
	get_vol_mergelist();
	cholmod_free_dense(&x, cm);
	cholmod_free_dense(&b, cm);
	cholmod_finish(&c);
}

// solve the node voltages using direct LU
void Circuit::solve_LU(){
	//solve_init();
	solve_LU_core();
}

// given vector x that obtained from LU, set the value to the corresponding
// node in nodelist
void Circuit::get_voltages_from_LU_sol(double * x){
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		size_t id = node->rep->rid;	// get rep's id in Vec
		double v = x[id];		// get its rep's value
		node->value = v;
	}
}

// compute value of mergelist nodes
void Circuit::get_vol_mergelist(){
	DIRECTION p, q;
	for(size_t i=0;i<mergelist.size();i++){
		Node * node = mergelist[i];
		// check direction
		if( node->nbr[WEST] != NULL ){
			p = WEST;
			q = EAST;
		}
		else{
			p = SOUTH;
			q = NORTH;
		}
		// assign vol value to node
		// left end node and its value
		double r1 = node->eqvr[p];
		double v1 = node->end[p]->value;
		//clog<<" left node: "<<r1<<" / "<<v1<<endl;
		//clog<<" left end "<<node->end[p]->name<<endl;
		// right end node and its value
		double r2 = node->eqvr[q];
		double v2 = node->end[q]->value;
		//clog<<"right node: "<<r2<<" / "<<v2<<endl;
		//clog<<"right end "<<node->end[q]->name<<endl;
		// value for node
		if(v1 > v2){
			node->value = v2 + (v1 - v2) * r2 / (r1 + r2);
		}
		else{
			node->value = v1 + (v2 - v1)  * r1 / (r1 + r2);
		}
		//clog<<" node "<<*node<<endl;
	}
}

// copy solution of block into circuit
void Circuit::get_voltages_from_block_LU_sol(){

	size_t block_id;
	for(size_t i=0;i<nodelist.size()-1;i++){
		Node * node = nodelist[i];
		if( node->is_mergeable() ) continue;
		block_id = node->rep->blocklist[0];
		Block &block = block_info[block_id];
		size_t id = node->rep->id_in_block[0];
		//Vec &p = block.x;
		//double v = p[id];		// get its rep's value
		double v = block.xp[id];
		node->value = v;
	}
}

// solve the circuit using PCG method
/*bool Circuit::solve_pcg(){
	solve_init();
	// remember to exclude ground node
	size_t n = nodelist.size()-1;

	// xk is for k-th iteration, and xk1 is for (k+1)-th iteration
	Vec b(n);
	Vec pk(n), pk1(n);
	Vec rk(n), rk1(n);
	Vec zk(n), zk1(n);
	Vec xk(n), xk1(n);
	double alpha, beta;
	bool successful = false;
	Matrix A, ML, MU, D;

	this->stamp_by_set(A, b);
	
	// initialize first iteration	
	init_precondition(A, ML, D, MU);
	copy_node_voltages(xk1, true);
	rk1 = b - A * xk1; 
	zk1 = compute_precondition(ML, D, MU, rk1);// solve M*z0=r0
	pk1 = zk1;
	int k=0; 
	while(k++ < MAX_ITERATION){
		// new values becomes old values
		zk = zk1;
		pk = pk1;
		xk = xk1;
		double diff=1.0; 

		alpha = (rk * zk) / (pk * A * pk);
		xk1 = xk + alpha * pk;
		rk1 = rk - alpha * (A * pk);
		diff = find_diff(rk1);
		if( diff < EPSILON ) {
			successful = true;
			break;
		}
		zk1 = compute_precondition(ML, D, MU, rk1);
		beta = (rk1 * zk1) / (rk * zk);
		pk1 = zk1 + beta * pk; 
	}
	if( successful ) copy_node_voltages(xk1, false);
	return successful;
}

// M1= (D/OMEGA+L); M2= (D/OMEGA+L'); D:diagonal
void Circuit::init_precondition(const Matrix &A, 
		Matrix &ML, Matrix & D, Matrix &MU){
	A.diagonal_split(ML, D, MU);
	ML *= 1/OMEGA;
	ML += D;

	MU *= 1/OMEGA;
	MU += D;
}

// solve M*z=r; 
// compute a preconditioner for matrix A = ML + D + MU
Vec Circuit::compute_precondition(const Matrix & ML, const Matrix & D, 
			const Matrix & MU, Vec &r){
	size_t n = nodelist.size()-1;
	Vec z1(n), z2(n), z(n);
		
	// solve (D/OMEGA+L)*z=r
	// D is the diagonal of M, L is the lower triangular
	//z1 = Algebra::solve(ML, r);	
	//Algebra::solve(ML, r, z1);	

	// solve (OMEGA/(2-OMEGA)*D^-1)*z=r: D is the diagonal of M,
	z2 = (2-OMEGA)/OMEGA * (D * z1);
	
	// solve (D/OMEGA+L')*z=r: D is the diagonal of M, L' is the upper 
	//z  = Algebra::solve(MU, r);
	//Algebra::solve(MU, r, z);
	return z;
}
*/

// 1. copy node voltages from the circuit to a Vec
//    from = true then copy circuit to x
//    else copy from x to circuit
// 2. map block voltage into global_index
void Circuit::copy_node_voltages_block(bool from){
	size_t id;
	if( from == true ){
		for(size_t i=0;i<replist.size();i++){
			Node *node = replist[i];
			const vector<size_t> &block_id = node->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->id_in_block[j];
				block.xp[id] = replist[i]->value;
				//block.x[id] = replist[i]->value;
				block.nodes[id] = replist[i];
			}
		}
	}
	else{
		for(size_t i=0;i<nodelist.size()-1;i++){
			Node * node = nodelist[i];
			const vector<size_t> &block_id = 
				node->rep->get_block_id();
			for(size_t j=0;j<block_id.size();j++){
				Block &block = block_info[block_id[j]];
				id = node->rep->id_in_block[j];
				node->value = block.xp[id];
				//Vec &p = block.x;
				//node->value = p[id];
			}
		}
	}
}

// stamp the net in each set, 
// *NOTE* at the same time insert the net into boundary netlist
void Circuit::stamp_by_set(Matrix & A, double* b){
	for(int type=0;type<NUM_NET_TYPE;type++){
		NetList & ns = net_set[type];
		NetList::iterator it;
		switch(type){
		case RESISTOR:
			for(it=ns.begin();it!=ns.end();++it){
				if( (*it) == NULL ) continue;
				assert( fzero((*it)->value) == false );
				stamp_resistor(A, (*it));
				//block_boundary_insert_net(ns[i]);
			}
			break;
		case CURRENT:
			for(it=ns.begin();it!=ns.end();++it)
				stamp_current(b, (*it));
			break;
		case VOLTAGE:
			for(it=ns.begin();it!=ns.end();++it){
				if( fzero((*it)->value)  && 
				    !(*it)->ab[0]->is_ground() &&
				    !(*it)->ab[1]->is_ground() )
					continue; // it's a 0v via
				stamp_VDD(A, b, (*it));
			}
			break;
		default:
			report_exit("Unknwon net type\n");
			break;
		}
	}
}

void Circuit::make_A_symmetric(double *b){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	Node *p, *q;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
		// node a points to X node
		if((*it)->ab[0]->isX()){
			p = (*it)->ab[0]; q = (*it)->ab[1];
		}
		else if((*it)->ab[1]->isX()){
			p = (*it)->ab[1]; q = (*it)->ab[0];
		}
		else continue;
		size_t id = q->rep->rid;
		double G = 1.0 / (*it)->value;
		b[id] += p->value * G;
	}
}

void Circuit::make_A_symmetric_block(){
	int type = RESISTOR;
	NetList & ns = net_set[type];
	NetList::iterator it;
	vector<size_t> *p, *q;
	Node *nk, *nl;

	for(it=ns.begin();it!=ns.end();it++){
		if( (*it) == NULL ) continue;
			assert( fzero((*it)->value) == false );
	
		Net *net = *it;
		Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};
		if(!nd[0]->isX() && !nd[1]->isX()) continue;
		double G;
		G = 1./net->value;
		vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
		vector<size_t >::const_iterator it_1;
		for(size_t j=0;j<2;j++){
			p = ls[j];
			q = ls[1-j];
			nk = nd[j];
		       	nl = nd[1-j];
			for(size_t i=0;i<p->size();i++){
				// find whether block_id in another list
				size_t block_id = (*p)[i];
				it_1 = find( (*q).begin(), (*q).end(), block_id);
				// 2 nodes in the same block
				if(it_1!=(*q).end() && !nk->isX()){
					Block &block = block_info[block_id];
					size_t k1 = nk->id_in_block[i];
					block.bp[k1] += G *(nl->value);
				}
			}
		}
	}
}
void Circuit::stamp_resistor(Matrix & A, Net * net){
	double G;
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	size_t k = nk->rid;
	size_t l = nl->rid;
	G = 1./net->value;
	if( !nk->isX() ) {
		A.push_back(k,k, G);
		if(!nl->isX() && l < k) // store lower triangular
			A.push_back(k,l,-G);
	}

	if( !nl->isX() ) {
		A.push_back(l,l, G);
		if(!nk->isX() && k < l) // store ower triangular
			A.push_back(l,k,-G);
	}
}

// stamp a current source
void Circuit::stamp_current(double * b, Net * net){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && !nk->isX() ) { 
		size_t k = nk->rid;
		b[k] += -net->value;
	}
	if( !nl->is_ground() && !nl->isX() ) {
		size_t l = nl->rid;
		b[l] +=  net->value;
	}
}

// stamp a voltage source
void Circuit::stamp_VDD(Matrix & A, double * b, Net * net){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	size_t id = X->rep->rid;
	A.push_back(id, id, 1.0);
	Net * south = X->rep->nbr[SOUTH];
	if( south != NULL &&
	    south->type == CURRENT ){
		// this node connects to a VDD and a current
		assert( feqn(1.0, b[id]) ); // the current should be stamped
		b[id] = net->value;	    // modify it
	}
	else
		b[id] += net->value;
}

// =========== stamp block version of matrix =======

void Circuit::stamp_block_resistor(Net * net, Matrix * A){
	Node * nd[] = {net->ab[0]->rep, net->ab[1]->rep};

	double G;	
	G = 1./net->value;

	vector<size_t> * ls[] = {&nd[0]->blocklist, &nd[1]->blocklist};
	vector<size_t>::const_iterator it;
	for(size_t j=0;j<2;j++){
		vector<size_t> *p = ls[j], *q = ls[1-j];
		Node *nk = nd[j], *nl = nd[1-j];
		for(size_t i=0;i<p->size();i++){
			// find whether block_id in another list
			size_t block_id = (*p)[i];
			it = find( (*q).begin(), (*q).end(), block_id);

			// block_id not in another list, net is boundary
			if( it == (*q).end() ){
				Block & blk = block_info[block_id];
				blk.boundary_netlist.push_back(net);
				if( !nk->isX() ) {
					// stamp value into block_ids
					size_t k1 = nk->id_in_block[i];
					Matrix &pk = A[block_id];	
					pk.push_back(k1,k1, G);
				}
			}
			// else 2 nodes belongs to the same block
			// stamp resistor
			else if( !nk->isX() ) {
				size_t k1 = nk->id_in_block[i];
				Matrix &pk = A[block_id];

				size_t j1 = it - (*q).begin();
				size_t l1 = nl->id_in_block[j1];

				pk.push_back(k1,k1, G);
				if(!nl->isX() && l1 < k1) // only store the lower triangular part
					pk.push_back(k1,l1,-G);
			}
		}// end of for k
	}// end of for j	
}

void Circuit::stamp_block_current(Net * net){
	Node * nk = net->ab[0]->rep;
	Node * nl = net->ab[1]->rep;
	if( !nk->is_ground() && !nk->isX() ) { 
		for(size_t i=0;i<nk->blocklist.size();i++){
			size_t block_idk = nk->blocklist[i];
			Block &block_k = block_info[block_idk];
			//Vec & pk = block_k.b;
			size_t k = nk->id_in_block[i];
			block_k.bp[k] += -net->value;
			//pk[k] += -net->value;
		}
	}
	if( !nl->is_ground() && !nl->isX() ) {
		for(size_t i=0;i<nl->blocklist.size();i++){
			size_t block_idl = nl->blocklist[i];
			Block & block_l = block_info[block_idl];
			//Vec & pl = block_l.b;
			size_t l = nl->id_in_block[i];
			block_l.bp[l] += net->value;
			//pl[l] +=  net->value;
		}
	}
}

void Circuit::stamp_block_VDD(Net * net, Matrix * A){
	// find the non-ground node
	Node * X = net->ab[0];
	if( X->is_ground() ) X = net->ab[1];
	for(size_t i=0;i<X->rep->id_in_block.size();i++){
		size_t block_id = X->rep->blocklist[i];
		Block &block = block_info[block_id];	
		Matrix & p = A[block_id];
		//Vec & q = block.b;
		size_t id =X->rep->id_in_block[i];
		p.push_back(id, id, 1.0);
		Net * south = X->rep->nbr[SOUTH];
		if( south != NULL &&
	   	 south->type == CURRENT ){
		// this node connects to a VDD and a current
		// the current should be stamped
			//assert( feqn(1.0, q[id]) ); 
			assert( feqn(1.0, block.bp[id]) );
			block.bp[id] = net->value;
			//q[id] = net->value;	    // modify it
		}
		else{
			block.bp[id] += net->value;
			//q[id] += net->value;
		}
	}
}

// set the block_id in a node
// according to the following order
// where 0 is the original block the node should be 
// 8 1 2
// 7 0 3
// 6 5 4
//
void Circuit::set_blocklist(Node * nd){
	const double len_per_block_x = block_info.len_per_block_x;
	const double len_per_block_y = block_info.len_per_block_y;
	const double len_ovr_x = block_info.len_ovr_x;
	const double len_ovr_y = block_info.len_ovr_y;
	const size_t X_BLOCKS = block_info.X_BLOCKS;
	const size_t Y_BLOCKS = block_info.Y_BLOCKS;
	const long x = nd->pt.x - x_min;	// point relative location
	const long y = nd->pt.y - y_min;
	size_t bx0 = x / len_per_block_x;
	size_t by0 = y / len_per_block_y;
	long bx, by;		// block index
	double lx, ly, ux, uy;	// block bounding box

	const long dx[]={0, 0, 1, 1,  1,  0, -1, -1, -1};
	const long dy[]={0, 1, 1, 0, -1, -1, -1,  0, 1};

	// test whether the node is in one of the nine blocks
	for(int i=0;i<9;i++){
		bx = bx0 + dx[i];
		by = by0 + dy[i];

		// check the block index is valid
		if( bx < 0 || bx >= (long)X_BLOCKS ) continue;
		if( by < 0 || by >= (long)Y_BLOCKS ) continue;

		// compute block coordinate
		lx = bx * len_per_block_x - len_ovr_x;
		ly = by * len_per_block_y - len_ovr_y;
		ux = (bx+1) * len_per_block_x + len_ovr_x;
		uy = (by+1) * len_per_block_y + len_ovr_y;

		// check if the point is in the block
		if( !(x>=lx && x<=ux && y>=ly && y<=uy) ) continue;	

		size_t id = by * X_BLOCKS + bx;
		assert( id<X_BLOCKS* Y_BLOCKS );

		/*
		vector<size_t >::const_iterator it;
		it = find(nd->blocklist.begin(),nd->blocklist.end(), id);
		if(it!=nd->blocklist.end()){
			printf("id=%ld, i=%d\n", id,i);
			printf("xbase,ybase=%lf %lf\n", lx_base, ly_base);
			printf("lx,ly=%lf %lf\n", lx, ly);
			printf("ux,uy=%lf %lf\n", ux, uy);
			printf("x,y=%ld %ld\n", x, y);
			printf("bx,by=%ld %ld\n", bx, by);
			printf("xc,yc=%ld %ld\n", x_center, y_center);
			continue;
		}
		*/
		nd->blocklist.push_back(id);
		Block & block = block_info[id];
		nd->id_in_block.push_back(block.count++);
	}
}

void Circuit::get_parameters(
		double & epsilon,
		double & omega,
		double & overlap_ratio,
		size_t & max_block_nodes,
		int & mode){
	epsilon		= EPSILON;
	omega		= OMEGA;
	overlap_ratio	= OVERLAP_RATIO; 
	max_block_nodes	= MAX_BLOCK_NODES;
	mode		= MODE;
}

// default values of these parameters are at the begining of this file
void Circuit::set_parameters(
		double epsilon, 
		double omega, 
		double overlap_ratio,
		size_t max_block_nodes,
		int mode){
	EPSILON		= epsilon;
	OMEGA		= omega;
	OVERLAP_RATIO 	= overlap_ratio;
	MAX_BLOCK_NODES	= max_block_nodes;
	MODE		= mode;
}

// choose an appropriate omega for the circuit s.t.
// - node size (use replist)
// - type (c4 or wb)
// - number of layers
void Circuit::select_omega(){
	double omega=OMEGA;
	size_t num_nodes = replist.size();
	size_t num_layers = layers.size();
	if( num_nodes < 0.05e6 )
		omega=1.0;
	else if (num_nodes < 0.2e6 )
		omega = 1.1;
	else if (num_nodes < 0.3e6 )
		omega = 1.2;
	else if (num_nodes < 0.5e6 )
		omega = 1.3;
	else if (num_nodes < 1.2e6 )
		omega = 1.4;
	else
		omega = 1.5;

	if( circuit_type == WB && num_layers >= 8 ) omega += 0.2;

	if( circuit_type == C4 ) omega += 0.1;

	if( name == "GND" && num_nodes < 1.2e6) omega -= 0.1;

	if( omega >= 1.6 ) omega = 1.6;
	if( omega <= 1.0 ) omega = 1.0;

	OMEGA = omega;
}

// Randomly choose a number of sample nodes to monitor
void Circuit::get_samples(){
	size_t num_nodes = replist.size();
	srand(time(NULL));
	while(sample.size()<SAMPLE_NUM_NODE){
		int id = rand() % num_nodes;
		sample.push_back(replist[id]);
	}
}

bool Circuit::check_diverge() const{
	for(size_t i=0;i<SAMPLE_NUM_NODE;i++){
		double x = sample[i]->value;
		if(VDD > 0){
			if( x < 0.0 || x > VDD ) return true;
		}
		else
			if(x<0.0) return true;
	}
	return false;
}

Node * Circuit::merge_along_dir_one_pass(Node * start, DIRECTION dir, bool remove){
	double sum = 0.0;
	DIRECTION ops = get_opposite_dir(dir);
	Node * p = start;

	// traverse along the direction, sum the resistor value and set the node end
	while(1){
		p = p->get_nbr_node(dir);
		p->end[ops] = start;
		Net * net = p->nbr[ops];
		sum += net->value;
		p->eqvr[ops] = sum;
		if( remove ) {
			size_t id = net_id[net];
			net_set[RESISTOR][id] = NULL;
			delete net;
		}
		if( !p->is_mergeable() ) break;
	}

	return p;	// return end point
}

// merge a line along direction
void Circuit::merge_along_dir(Node * node, DIRECTION dir){
	// two pass traversal
	DIRECTION ops = get_opposite_dir(dir);
	node->end[dir] = merge_along_dir_one_pass(node, dir, false);
	Node * other = node->end[dir];
	other->end[ops] = node;
	merge_along_dir_one_pass(other, ops, true);
	//assert( ret == node );

	// add a new net between `node' and its end
	Net * net = new Net(RESISTOR, node->eqvr[dir], node, other);
	node->nbr[dir] = other->nbr[ops] = net;
	this->add_net(net);
}

double Circuit::locate_maxIRdrop(){
	max_IRdrop = 0;
	for(size_t i=0;i<nodelist.size()-1;i++){
		double IR_drop = VDD - nodelist[i]->value;		
		if(IR_drop > max_IRdrop)
			max_IRdrop = IR_drop;
	}
	return max_IRdrop;
}

double Circuit::locate_special_maxIRdrop(){
	double max_IRdrop = 0;
	Node *nd;
	for(size_t j=0;j<special_nodes.size(); 
			j++){
		nd = special_nodes[j];
		double IR_drop = VDD - nd->value;	
		if(IR_drop > max_IRdrop)
			max_IRdrop = IR_drop;
	}
	return max_IRdrop;
}

void Circuit::build_pad_set(){
	pad_set.resize(0);
	//static Pad pad_a;
	for(size_t i=0;i<nodelist.size()-1;i++){
		if(nodelist[i]->isX()){
			Pad *pad_ptr = new Pad();
			pad_ptr->node = nodelist[i];
			pad_set.push_back(pad_ptr);
		}
	}
	//for(size_t j=0;j<pad_set.size();j++)
		//clog<<"pad: "<<*pad_set[j]->node<<endl;
}

void Circuit::extract_pads(int pad_number){
	vector<Node*> pair_first;
	vector<double> pair_second;
	pair<Node*, double> pair_nd;
	//int pad_number = 5;
	double distance = 0;
	map<Node*, double>::iterator it;

	clear_pad_control_nodes();
	for(size_t i=0;i<special_nodes.size();i++){
		int count = 0;
		pair_first.clear();
		pair_second.clear();
		Node *nd = special_nodes[i];
		// search for closest pads
		for(size_t j=0;j<pad_set.size();j++){
			Node *ptr = pad_set[j]->node;
			distance = get_distance(ptr, nd);

			if(count < pad_number){
				pair_first.push_back(ptr);
				pair_second.push_back(distance);
				count++;
			}else{// substitute the pad node
				double max_dist = 0;
				size_t max_index = 0;
				for(size_t k=0;k<pair_second.size();k++){
					if(pair_second[k]>max_dist){
						max_dist = pair_second[k];
						max_index = k;
					}
				}
				if(distance < max_dist){ 
					pair_first[max_index] = ptr;
					pair_second[max_index] = distance;
				}
			}
		}
		// then map these distance into pads
		for(size_t j=0;j<pair_first.size();j++){
			Node *ptr = pair_first[j];
			//if(nd->name == "n0_0_0")
			//clog<<"ptr: "<<*ptr<<endl;
			for(size_t k=0;k<pad_set.size();k++){
				if(pad_set[k]->node->name == ptr->name){
					// control nodes
					pair_nd.first = nd;
					// distance
					pair_nd.second = pair_second[j];
					pad_set[k]->control_nodes.insert(pair_nd);
					break;
				}
			}
		}
	}
	//print_pad_map();	
	pair_first.clear();
	pair_second.clear();
	//print_pad_map();
}

void Circuit::extract_min_max_pads_new(vector<double> ref_drop_vec){
	Node *nd;
	Pad *pad;
	size_t max_index = 0;
	double max = 0;

	vector<Node*> min_pads;
	vector<Node*> max_pads;
	vector<bool >temp_flag;

	map<Node *, double>::iterator it;
	temp_flag.resize(pad_set.size());
	for(size_t i=0;i<temp_flag.size();i++)
		temp_flag[i] = false;
	size_t id_minpad;
	double drop = 0;
	double avg_ref = calc_avg_ref(ref_drop_vec);
	double avg_drop = VDD - avg_ref;
	
	for(size_t i=0;i<pad_set.size();i++){
		pad = pad_set[i];
		nd = pad->node;
		drop = VDD - ref_drop_vec[i];
		if(drop>max){
			max = drop;
			max_index = i;
		}
		if(drop < 0.7*avg_drop){
			min_pads.push_back(nd);	
		}	
	}
	double max_id;
	do{
		double max_temp = 0;
		max_id = -1;
		for(size_t j=0;j<ref_drop_vec.size();j++){
			if(VDD - ref_drop_vec[j] < max*0.9)
				continue;
			if(temp_flag[j] ==  false){
				if(max_id == -1)
					max_id = j;
				if(VDD - ref_drop_vec[j] > max_temp){
					max_temp = VDD - ref_drop_vec[j];
					max_id = j;
				}
			}
		}
		if(max_id == -1) break;
		temp_flag[max_id] = true;
		if(max_temp >= max*0.9)	
			max_pads.push_back(pad_set[max_id]->node);
	}while(max_id != -1);

	temp_flag.clear();
	Node *new_pad;
	Pad * pad_ptr;

	// set nd into the weighted center
	// start to map min pads into max pads locations
	for(size_t j=0;j<min_pads.size();j=j+2){
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				min_pads[j]->name){
				id_minpad = k;
			}
		}
	
		size_t i = j % max_pads.size();
		//size_t i = locate_max_pad(max_pads, iter);
		Node * nd = max_pads[i];
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				nd->name){	
				pad_ptr = pad_set[k];
				double ref_drop_value = ref_drop_vec[k];

				new_pad = pad_projection(pad_ptr, min_pads[j]);
			if(min_pads[j]->name == "n0_135_104" ||
			min_pads[j]->name == "n0_30_74"||
			min_pads[j]->name == "n0_81_29" ||
			min_pads[j]->name =="n0_186_47"||
			min_pads[j]->name == "n0_255_59"||
			min_pads[j]->name == "n0_67_148")
			clog<<"old pad / new pad: "<<*min_pads[j]<<" "<<*new_pad<<endl;

				size_t m = id_minpad;
				pad_set[m]->node->disableX();
				pad_set[m]->node->value = 0;
				pad_set[m]->node = new_pad;
				pad_set[m]->visit_flag = true;
				// already taken care of
				pad_set[m]->control_nodes.clear();
				break;
			}
		}

	}
	// next step is to insert the min pads into max pads area
	min_pads.clear();
	max_pads.clear();
	
}

void Circuit::extract_min_max_pads(vector<double> ref_drop_vec){
	double min = VDD-ref_drop_vec[0];
	double max = 0;	
	size_t min_index=0;
	size_t max_index=0;
	vector<Node*> min_pads;
	vector<Node*> max_pads;
	size_t id_minpad = 0;
	/*vector<bool >temp_flag;
	temp_flag.resize(pad_set.size());
	for(size_t i=0;i<temp_flag.size();i++)
		temp_flag[i] = false;*/

	for(size_t j=0;j<ref_drop_vec.size();j++){
		if(VDD-ref_drop_vec[j]>max){
			max = VDD - ref_drop_vec[j];
			max_index = j;
		}
		if(VDD-ref_drop_vec[j]<min){
			min = VDD - ref_drop_vec[j];
			min_index=j;
		}
	}
	for(size_t j=0;j<ref_drop_vec.size();j++){
		if(VDD - ref_drop_vec[j]<= max*0.5){
			min_pads.push_back(pad_set[j]->node);
		}
		if(VDD-ref_drop_vec[j]>=max*0.8){
			max_pads.push_back(pad_set[j]->node);
		}
	}
	/*double max_id;
	do{
		double max_temp = 0;
		max_id = -1;
		for(size_t j=0;j<ref_drop_vec.size();j++){
			if(VDD - ref_drop_vec[j] < max*0.9)
				continue;
			if(temp_flag[j] ==  false){
				if(max_id == -1)
					max_id = j;
				if(VDD - ref_drop_vec[j] > max_temp){
					max_temp = VDD - ref_drop_vec[j];
					max_id = j;
				}
			}
		}
		if(max_id == -1) break;
		temp_flag[max_id] = true;
		if(max_temp >= max*0.8)	
			max_pads.push_back(pad_set[max_id]->node);
	}while(max_id != -1);

	temp_flag.clear();*/

	Node *new_pad;
	Pad * pad_ptr;

	// set nd into the weighted center
	// start to map min pads into max pads locations
	for(size_t j=0;j<min_pads.size();j++){
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				min_pads[j]->name){
				id_minpad = k;
			}
		}
		size_t i = j % max_pads.size();
		Node * nd = max_pads[i];
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				nd->name){	
				pad_ptr = pad_set[k];
				double ref_drop_value = ref_drop_vec[k];
				new_pad = pad_projection(pad_ptr, min_pads[j]);
				size_t m = id_minpad;		
				pad_set[m]->node->disableX();
				pad_set[m]->node->value = 0;
				pad_set[m]->node = new_pad;
				//pad_set[m]->visit_flag = true;
				// already taken care of
				pad_set[m]->control_nodes.clear();
				break;
			}
		}

	}
	//for(size_t j=0;j<min_pads.size();j++){
		//clog<<"min_pads: "<<*min_pads[j]<<endl;
	//}
	//for(size_t j=0;j<max_pads.size();j++){
		//clog<<"max_pads: "<<*max_pads[j]<<endl;
	//}

	min_pads.clear();
	max_pads.clear();
}

// tune 50% nodes with the small IR drops
void Circuit::update_pad_control_nodes(vector<double> & ref_drop_value, size_t iter){
	ref_drop_value.resize(pad_set.size());
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		double middle_value = locate_ref(i);
		ref_drop_value[i] = middle_value;
		//clog<<"middle value: "<<middle_value<<endl;
	}
}

/*void Circuit::update_queue(Node *nds, queue<Node*> &q, Node *nd, int &count, double &dist){
	Net * net; Node *nbr;
	Node *na, *nb;
	for(int i=0;i<6;i++){
		//if(count >=pad_number) break;
		net = nd->nbr[i];
		if(net==NULL) continue;
		// find the other node, which is nbr
		na = net->ab[0];
		nb = net->ab[1];
		if(nd->name == na->name)
			nbr = nb;
		else	nbr = na;
		if(!nbr->is_ground()&& !nbr->region_flag){
			//cout<<"nbr: "<<*nbr<<endl;
			dist += net->value;
			q.push(nbr);
			nbr->region_flag = true;
			if(nbr->isX()== true){
				count ++;
				assign_distance(nds, nbr, dist);
				
				//cout<<"add a new pad: "<<*nbr<<" "<<count<<endl;
			}
		}
	}
}*/

void Circuit::assign_distance(Node *nds, Node *nd, double dist){
	Node *ndt;
	pair<Node*, double> pad_pair;

	
	for(size_t i=0;i<pad_set.size();i++){
		ndt = pad_set[i]->node;
		if(nd->name == ndt->name){
			cout<<"pad: "<<*ndt<<endl;
			cout<<"source node, dist: "<<*nds<<" "<<dist<<endl;
			pad_pair.first = nds;
			pad_pair.second = dist*fabs(VDD-nds->value);
			pad_set[i]->control_nodes.insert(pad_pair);
			break;
		}
	}
}


void Circuit::print_pad_map(){
	Node *nd;
	Pad *pad;
	pair<Node*, double> pad_pair;
	map<Node*, double>::iterator it;
	for(size_t i=0;i<pad_set.size();i++){
		nd = pad_set[i]->node;
		pad = pad_set[i];
		if(nd->name == "n0_477_17"){// ||
		   //nd->name == "n0_255_59" ||
		   //nd->name == "n0_30_74"){
		for(it = pad->control_nodes.begin();it != pad->control_nodes.end();it++){
			printf("%ld %ld  %.5e\n", it->first->pt.y+1, it->first->pt.x+1, it->first->value);
		}
		}
	}
}

void Circuit::clear_flags(){
	Pad *pad;
	for(size_t i=0;i<pad_set.size();i++){
		pad = pad_set[i];
		pad->visit_flag = false;
		pad->fix_flag = false;
		pad->violate_flag = false;
		//nd->region_flag = false;
		//nd->fix_flag = false;
		//nd->visit_flag = false;
	}
}

void Circuit::clear_pad_control_nodes(){
	Node *nd;
	for(size_t i=0;i<pad_set.size();i++){
		pad_set[i]->control_nodes.clear();
	}
}

/*void Circuit::relocate_pads(){
	vector<Node*> pad_set_old;
	double dist = 0;
	double new_dist = 0;
	//for(size_t i=0;i<6;i++){
	pad_set_old.resize(pad_set.size());
	assign_pad_set(pad_set_old);
	// store a original copy of the pad set
	
	// build up the map for nodes in PAD layer
	build_map_node_pt();

	mark_special_nodes();

	vector<double> ref_drop_vec;
	//print_pad_set();
	for(size_t i=0;i<1;i++){
		int pad_number = 1;
		origin_pad_set.resize(pad_set.size());
		assign_pad_set(origin_pad_set);
		
		// find control nodes for each pad
		extract_pads(pad_number);
		// find the tune spot for control nodes	
		update_pad_control_nodes(ref_drop_vec, i);
		// find new point for all pads	
		dist = update_pad_pos_all(ref_drop_vec);
		// move the low 10% pads into high 10% 
		// pad area 
		if(i==0)
			extract_min_max_pads(ref_drop_vec);
		// update the old pad set value
		assign_pad_set(pad_set_old);
		// actual move pads into the new spots
		project_pads();

		rebuild_voltage_nets();
		solve_LU_core();
		double max_IR = locate_maxIRdrop();	
		//double max_IRS = locate_special_maxIRdrop();
		clog<<"i, max_IR is: "<<i<<" "<<max_IR<<endl;
		//clog<<"max_IRS is: "<<max_IRS<<endl<<endl;
	}
	ref_drop_vec.clear();
	map_node_pt.clear();
	print_pad_set();
}*/

void Circuit::relocate_pads_graph(){
	vector<Node*> pad_set_old;
	double dist = 0;
	double new_dist = 0;
	pad_set_old.resize(pad_set.size());
	assign_pad_set(pad_set_old);
	// store a original copy of the pad set
	
	// build up the map for nodes in PAD layer
	build_map_node_pt();
	mark_special_nodes();
	
	vector<double> ref_drop_vec;
	//print_pad_set();
	for(size_t i=0;i<12;i++){
		int pad_number = 1;
		origin_pad_set.resize(pad_set.size());
		assign_pad_set(origin_pad_set);
		// build pad connection graph
		build_graph();
		// find control nodes for each pad
		extract_pads(pad_number);
		// find the tune spot for control nodes	
		update_pad_control_nodes(ref_drop_vec, i);
		//print_all_control_nodes();	
		if(i>=6)
			dynamic_update_violate_ref(ref_drop_vec);
		// find new point for all pads	
		dist = update_pad_pos_all(ref_drop_vec);		
		// move the low 10% pads into high 10% 
		// pad area 
		if(i==0)
			extract_min_max_pads(ref_drop_vec);
		// update the old pad set value
		assign_pad_set(pad_set_old);
		
		move_violate_pads(ref_drop_vec);	
		
		// actual move pads into the new spots
		//project_pads();

		// move pads according to graph contraints
		graph_move_pads(ref_drop_vec);	
		
		clear_flags();
		// actual move pads into the new spots
		// project_pads();
		
		resolve_direct();
		//resolve_queue(origin_pad_set);
		//solve_GS();
		//clog<<"max_IRS is: "<<max_IRS<<endl<<endl;
	}
	ref_drop_vec.clear();
	map_node_pt.clear();
	origin_pad_set.clear();
	pad_set_old.clear();
	//print_pad_set();
}

double Circuit::calc_avg_ref_drop(vector<double> &ref_drop_vec){
	Node *pad;
	Pad *pad_ptr;
	double max_drop, min_drop;
	double sum_max = 0;
	double sum_min = 0;
	double sum_diff = 0;
	size_t count = 0;
	double ref_drop_value = 0;

	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		count ++;	
		ref_drop_value = ref_drop_vec[i];
		map<Node *, double>::iterator it;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		max_drop = 0;
		min_drop = -1;
		
		for(it = pad_ptr->control_nodes.begin();
		    it != pad_ptr->control_nodes.end();
		    it++){
			if(it->second > ref_drop_value)
				continue;
			  if(it->second > max_drop)
				max_drop = it->second;
			  if(min_drop == -1)
				min_drop = it->second;
			  else if(it->second < min_drop)
				min_drop = it->second;
		}
		pad_ptr->data = max_drop - min_drop;
		sum_diff += pad_ptr->data;
	}
	double avg_drop = sum_diff / count;
	return avg_drop;
}

double Circuit::calc_avg_ref(vector<double> ref_drop_vec){
	Node *pad;
	Pad *pad_ptr;
	double sum_ref = 0;
	size_t count = 0;
	double ref_drop_value = 0;

	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		count ++;	
		ref_drop_value = ref_drop_vec[i];
		sum_ref += ref_drop_value;
	}
	double avg_ref = sum_ref / count;
	return avg_ref;
}

double Circuit::update_pad_pos_all(vector<double> ref_drop_vec){
	double total_dist = 0;
	double dist = 0;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		double ref_drop_value = ref_drop_vec[i];

		dist = update_pad_pos(ref_drop_value, i);
		total_dist += dist;
	}
	return total_dist;
}

void Circuit::modify_newxy(){
	Pad *pad;
	Node *nd;
	double delta_x;
	double delta_y;
	for(size_t i=0;i<pad_set.size();i++){
		pad = pad_set[i];
		nd = pad->node;
		if(pad->newx == nd->pt.x &&
		   pad->newy == nd->pt.y) 
			continue;
		delta_x = pad->newx - nd->pt.x;
		delta_y = pad->newy - nd->pt.y;
		delta_x /=2;
		delta_y /=2;
		pad->newx += delta_x;
		pad->newy += delta_y;
		
		round_data(pad->newx);
		round_data(pad->newy);
	}
}
// decide pad's new pos with the weights
// need to be tuned
double Circuit::update_pad_pos(double ref_drop_value, size_t i){
	double total_dist=0;
	Node *pad;
	Pad *pad_ptr;
	Node *nd;
	double weight = 0;
	double distance = 0;
	double pad_newx;
	double pad_newy;
	map<Node *, double>::iterator it;

		double sum_weight = 0;
		double weighted_x =0;
		double weighted_y =0;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		for(it = pad_ptr->control_nodes.begin();
		    it != pad_ptr->control_nodes.end();
		    it++){
			if(it->second > ref_drop_value)
				continue;
		 	  
			if((pad_set[i]->node->name == "n0_30_74" ||
			    pad_set[i]->node->name == "n0_135_104"||
			    pad_set[i]->node->name == "n0_255_59")){
				//clog<<"data: "<<pad_set[i]->data<<endl;
			  //cout<<"control node: "<<*it->first<<" "<<it->second<<endl;
			 printf("%ld %ld  %.5e\n", it->first->pt.y+1, it->first->pt.x+1, it->first->value);
			
			}
			nd = it->first;
			weight = 1.0/it->second;
			weighted_x += weight * nd->pt.x;
			weighted_y += weight * nd->pt.y;
			sum_weight += weight; 	
		}
		
		if(sum_weight !=0){
			pad_newx = weighted_x / sum_weight;
			pad_newy = weighted_y / sum_weight;

			// pad_newx = (pad_newx - pad->pt.x)/2+pad->pt.x;
			//pad_newy = (pad_newy - pad->pt.y)/2+pad->pt.y;
			round_data(pad_newx);
			round_data(pad_newy);

			/*if(pad->name == "n0_67_159")
				clog<<"band pad, new: "<<*pad_ptr->node<<" "<<pad_newx<<" "<<pad_newy<<endl;*/

			if((pad_ptr->node->pt.x > 300 || pad_ptr->node->pt.y > 150) && (pad_newx <= 300 && pad_newy <= 150)){
				//clog<<"band pad, new: "<<*pad_ptr->node<<" "<<pad_newx<<" "<<pad_newy<<endl;
				pad_ptr->control_nodes.clear();
				pad_ptr->visit_flag = true;
				
			}else{
				

				pad_ptr->newx = pad_newx;
				pad_ptr->newy = pad_newy;
			}
		}else{
			pad_ptr->newx = pad->pt.x;
			pad_ptr->newy = pad->pt.y;
		}
 
		double dist = sqrt(weighted_x*weighted_x 			 + weighted_y*weighted_y);
		//total_dist += temp;
	//}

	//clog<<"dist: "<<total_dist<<endl<<endl;
	return dist;
}

void Circuit::project_pads(){
	Node *pad;
	Node *new_pad;
	Pad *pad_ptr;

	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		// search for nearest node to (pad_newx,
		// pad_newy)
		
		new_pad = pad_projection(pad_ptr, pad);
		if(pad->name == "n0_135_104" ||
			pad->name == "n0_255_159"||
			pad->name == "n0_30_74")
			clog<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;
		// update pad information
		pad_ptr->node = new_pad;
		pad_ptr->control_nodes.clear();	
	}
}

// round double, depends on whether the frac >=0.5
void Circuit::round_data(double &data){
	double fractpart, intpart;
	fractpart = modf(data, &intpart);
	if(fractpart >= 0.5)
		data = ceil(data);
	else
		data = floor(data);
}

// expand from (x,y) to nearest node in grid
// fits for non-uniform grid
Node * Circuit::pad_projection(Pad *pad, Node *nd){
	//Node *nd;
	queue<Point> q;
	Point pt;
	Point pt_cur;
	stringstream sstream;
	string pt_name;
	Node *nd_new=NULL;
	//int gap - 10;
	double dx[4] = {30, 0, -30, 0};
	double dy[4] = {0, 30, 0, -30};

	//nd = pad->node;
	pt.z = nd->get_layer();
	pt.x = pad->newx;
	pt.y = pad->newy;
	if(pad->newx == nd->pt.x && 
		pad->newy == nd->pt.y){
		nd_new = nd;	
		return nd_new;
	}

	sstream<<pt.z<<"_"<<pt.x<<"_"<<pt.y; 
	pt_name = sstream.str();
	// first see if this node is on grid
	// and if it is occupied by pad or not
	//clog<<"orig pad: "<<*nd<<endl;
	if(has_node_pt(pt_name)){
		nd_new = get_node_pt(pt_name);
		//clog<<"has new node: "<<*nd_new<<" "<<nd_new->isX()<<endl;
		// if this node is not occupied by pad
		if(!nd_new->isX()){
			nd->disableX();
			nd->value = 0;
			nd_new->enableX();
			nd_new->value = VDD;
			return nd_new;
		}
	}
	bool return_flag = false;
	// else start to search for node
	q.push(pt);
	// if not, expand it to neighboring area
	while(!q.empty()&& return_flag == false){
		pt_cur = q.front();
		Point pt_nbr = pt_cur;
		//expand_pad_pos(q, pt_cur);	
		for(size_t i=0;i<4;i++){
			pt_nbr.x = pt_cur.x + dx[i];
			pt_nbr.y = pt_cur.y + dy[i];
			stringstream sstream;
			string pt_name;
			sstream << pt_nbr.z<<"_"<<
				pt_nbr.x<<"_"<<
				pt_nbr.y;
			pt_name = sstream.str();
			if(has_node_pt(pt_name)){
				nd_new = get_node_pt(pt_name);

				//cout<<"new name: "<<*nd_new<<" "<<nd_new->isX()<<endl;
				if(!nd_new->isX()){
					nd->disableX();
					nd->value = 0;
					nd_new->enableX();
					nd_new->value = VDD;
					return_flag = true;
					break;
				}
			}
			q.push(pt_nbr);
		}
		q.pop();
	}
	while(!q.empty()){
		q.pop();
	}
	if(return_flag == true)
		return nd_new;
	clog<<"no point for new pad. return. "<<endl;
	return NULL;
}

void Circuit::build_map_node_pt(){
	if(pad_set.size()==0)
		clog<<"no pad on grid. ERROR"<<endl;
	// ref layer
	int ref_layer = pad_set[0]->node->get_layer();

	Node *nd;
	pair<string, Node*> pt_pair;
	for(size_t i=0;i<nodelist.size()-1;i++){
		nd = nodelist[i];
		if(nd->get_layer()!=ref_layer)
			continue;
		stringstream sstream;
		sstream<<ref_layer<<"_"<<nd->pt.x<<
			"_"<<nd->pt.y;
		pt_pair.first = sstream.str();
		//cout<<"string: "<<pt_pair.first<<endl;
		pt_pair.second = nd;
		map_node_pt.insert(pt_pair);
	}
}

void Circuit::restore_pad_set(vector<Node*>&pad_set_old){
	Node *nd_old=NULL;
	Node *nd_new=NULL;
	for(size_t i=0;i<pad_set_old.size();i++){
		// have to use find
		if(pad_set[i]->node !=
				pad_set_old[i]){
			nd_old = pad_set_old[i];
			nd_new = pad_set[i]->node;
			//clog<<"nd_old / nd_new: "<<
			//*nd_old<<" "<<*nd_new<<endl;
			nd_new->disableX();
			nd_new->value = 0;
			nd_old->enableX();
			nd_old->value = VDD;		
		}
		pad_set[i]->node = nd_old;	
	}
}

void Circuit::assign_pad_set(vector<Node*>&pad_set_old){
	//clog<<"assign pad set."<<endl;
	for(size_t i=0;i<pad_set_old.size();i++){
		pad_set_old[i] = pad_set[i]->node;
		//clog<<"pad: "<<i<<" "<<*pad_set_old[i]<<endl;	
	}
}

// modify voltage net set with rm_node and add_node
void Circuit::rebuild_voltage_nets(){
	int type = VOLTAGE;
	size_t index_rm_net = 0;
	Net *net=NULL;
	Net *add_net=NULL;
	Node *nd_ori=NULL;
	Node *nd_new=NULL;
	Node *rm_node=NULL;
	Node *add_node=NULL;
	vector<Net*> rm_net;
	// delete all origin pad set
	// and build nets of new pad set
	for(size_t i=0;i<origin_pad_set.size();i++){
		rm_node = origin_pad_set[i];
		//rm_node = pad_set_old[i];
		add_node = pad_set[i]->node;
		//clog<<"nd_old, nd_new: "<<*rm_node<<" "
			//<<*add_node<<endl;

		for(size_t i=0;i<net_set[type].size();i++){
			net = net_set[type][i];
			if(net->ab[0]->name == rm_node->name ||
					net->ab[1]->name == rm_node->name){
				index_rm_net = i;
				rm_net.push_back(net);
				break;
			}
		}
		add_net = new Net(VOLTAGE, VDD, add_node, 
				nodelist[nodelist.size()-1]);
		net_set[type][index_rm_net] = add_net;
	}
	for(size_t i=0;i<rm_net.size();i++){
		delete rm_net[i];
	}
	origin_pad_set.clear();
}

void Circuit::print_pad_set(){
	for(size_t i=0;i<pad_set.size();i++){
		clog<<"pad: "<<*pad_set[i]->node<<endl;
		//printf("%d %d\n", pad_set[i]->node->pt.x+1, pad_set[i]->node->pt.y+1);
	}		
}

void Circuit::print_matlab(){
	// uncomment this if want to output to a file
	//freopen("output.txt","w",stdout);

	// don't output ground node
	size_t count = 0;
	for(size_t i=0;i<nodelist.size()-1;i++){
		printf("%ld %ld  %.5e\n", nodelist[i]->pt.y+1, nodelist[i]->pt.x+1, 
				VDD-nodelist[i]->value);
		if(nodelist[i]->value==VDD)
			count += 1;
	}
	clog<<"VDD # is: "<<count<<endl;
}

// build graph for pad nodes
void Circuit::build_graph(){
	Pad *pad;
	Pad *pad_nbr;
	Node *nd;
	bool flag_pad = false;
	bool flag_nbr = false;
	// clear content
	for(size_t i=0;i<pad_set.size();i++){
		pad_set[i]->nbrs.clear();
	}
	// find nbr pad nodes
	for(size_t i=0;i<pad_set.size();i++){
		//clog<<"pad: "<<*pad->node<<endl;
		flag_pad = false;
		flag_nbr = false;
		pad = pad_set[i];
		pad_nbr = find_nbr_pad(pad);
		for(size_t j=0;j<pad_nbr->nbrs.size();j++){
			if(pad_nbr->nbrs[j]->node->name== pad->node->name)
				flag_pad = true;
				break;
		}
		for(size_t j=0;j<pad->nbrs.size();j++){
			if(pad->nbrs[j]->node->name== pad_nbr->node->name)
				flag_nbr = true;
				break;
		}
		if(flag_pad == false){
			pad->nbrs.push_back(pad_nbr);
		}
		if(flag_nbr == false)
			pad_nbr->nbrs.push_back(pad);
	}
	/*for(size_t i=0;i<pad_set.size();i++){
		Pad *pad = pad_set[i];
		clog<<"pad: "<<*pad->node<<endl;
		for(size_t j=0;j<pad->nbrs.size();j++){
			clog<<"nbr: "<<*pad->nbrs[j]->node<<endl;
		}
	}*/
}

// use Euclidiean distance to locate nearest nbr pad
Pad * Circuit::find_nbr_pad(Pad *pad){
	Pad * nbr;
	double distance=-1;
	double min_dist=0;
	bool flag = false;
	size_t min_index=0;
	for(size_t i=0;i<pad_set.size();i++){
		nbr = pad_set[i];
		if(nbr->node->name == pad->node->name)
			continue;
		distance = get_distance(nbr->node, pad->node);
		if(flag == false){
			flag = true;
			min_dist = distance;
			min_index = i;
		}else{
			if(distance < min_dist){
				min_dist = distance;	
				min_index = i;
			}	
		}
	}
	return pad_set[min_index];
}

double Circuit::get_distance(Node *na, Node *nb){
	double distance = 0;
	double delta_x = 0;
	double delta_y = 0;
	delta_x=(na->pt.x-nb->pt.x);
	delta_y=(na->pt.y-nb->pt.y);
	delta_x *= delta_x;
	delta_y *= delta_y;
	distance = sqrt(delta_x + delta_y);
	return distance;
}

void Circuit::graph_move_pads(vector<double> ref_drop_vec){
	Node *new_pad;
	int id=0;
	do{
		id = locate_max_drop_pad(ref_drop_vec);
		if(id==-1) break;
		Pad *pad_ptr = pad_set[id];
		Pad *pad_nbr = NULL;
		Node *pad = pad_ptr->node;
		//clog<<endl<<"pad: "<<*pad<<endl;
		new_pad = pad_projection(pad_ptr, pad);

		//bool flag = print_flag(pad);
		//if(flag == true || new_pad->name == "n0_86_30" || new_pad->name =="n0_477_17" || pad->name == "n0_477_17")
			//clog<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;

		//clog<<"new_pad: "<<*new_pad<<endl;
		pad_ptr->visit_flag = true;
		for(size_t i=0;i<pad_ptr->nbrs.size();i++){
			pad_nbr = pad_ptr->nbrs[i];
			if(pad_nbr->fix_flag == false){
				pad_nbr->fix_flag = true;
				//clog<<"pad_nbr: "<<*pad_nbr->node<<endl;
			}
		}

		pad_ptr->node = new_pad;
		pad_ptr->control_nodes.clear();
	}while(id != -1);
}

// locate id that has minimum value and not visited or fixed 
int Circuit::locate_max_drop_pad(vector<double> vec){
	int min_id = -1;
	double min_ref = 0;
	bool flag = false;
	for(size_t i=0;i<vec.size();i++){
		if(pad_set[i]->visit_flag == true ||
			pad_set[i]->fix_flag ==  true)
			continue;
		//clog<<"i, vec: "<<i<<" "<<vec[i]<<endl;
		if(flag == false){
			flag = true;
			min_ref = vec[i];
			min_id = i;
		}
		else if(vec[i] < min_ref){
			min_ref = vec[i];
			min_id = i;
		}
	}
	/*double max_id=-1;
	double max_ref = 0;
	for(size_t i=0;i<vec.size();i++){
		if(pad_set[i]->visit_flag == true ||
			pad_set[i]->fix_flag ==  true)
			continue;
		//clog<<"i, vec: "<<i<<" "<<vec[i]<<endl;
		if(vec[i] > max_ref){
			max_ref = vec[i];
			max_id = i;
		}
	}*/
	//clog<<"min_ref, min_pad: "<<min_ref<<" "<<*pad_set[min_id]->node<<endl;
	return min_id;
	//return max_id;
}

// locate the tune spot for the control nodes.
double Circuit::locate_ref(size_t i){
	Pad *pad_ptr;
	Node *pad;
	map<Node*, double>::iterator it;
	Node *nd;
	double weight = 0;
	//vector<double> drop_vec;
	pad_ptr = pad_set[i];
	pad = pad_ptr->node;
	pad_ptr->drop_vec.clear();
	for(it = pad_ptr->control_nodes.begin();
			it != pad_ptr->control_nodes.end();
			it++){
		nd = it->first;
		weight = nd->value;
		if(weight <0)
			weight *=10;

		pad_ptr->control_nodes[nd] = weight;
		pad_ptr->drop_vec.push_back(nd->value); 
	}
	sort(pad_ptr->drop_vec.begin(), pad_ptr->drop_vec.end(),
			compare_values);
	pad_ptr->ratio = 2;
	size_t id = pad_ptr->drop_vec.size() / pad_ptr->ratio;
	double middle_value = pad_ptr->drop_vec[id];
	//drop_vec.clear();
	return middle_value;
}

void Circuit::dynamic_update_violate_ref(vector<double> & ref_drop_vec){
	//for(size_t j=0;j<2;j++){
	double avg_drop = calc_avg_ref_drop(ref_drop_vec);
	Pad *pad_ptr;
	Node *pad;
	//cout<<"j: "<<j<<endl;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;

		if(pad_ptr->data >= 2*avg_drop){
			pad_ptr->violate_flag = true;
			double ratio_new = pad_ptr->ratio * 2;
			size_t id = pad_ptr->drop_vec.size() / ratio_new;
			pad_ptr->ratio = ratio_new;
			double middle_value = pad_ptr->drop_vec[id];
			ref_drop_vec[i] = middle_value;
		}
	}
	
	extract_min_max_pads_new(ref_drop_vec);	
}

bool Circuit::print_flag(Node *pad){
	bool flag = false;
	if(pad->name == "n0_135_104" ||
	   pad->name == "n0_255_59"||
	   pad->name == "n0_30_74")
	   /*pad->name == "n0_67_148" ||
	   pad->name == "n0_114_171" ||
	   pad->name == "n0_162_159" ||
	   pad->name == "n0_216_171" ||
	   pad->name == "n0_268_177")*/
		flag = true;
	return flag;
}

void Circuit::move_violate_pads(vector<double> ref_drop_vec){
	Pad *pad_ptr;
	Node * pad;
	Node * new_pad;
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;

		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		// if violate, move this pad
		if(pad_ptr->violate_flag == true){
			new_pad = pad_projection(pad_ptr, pad);
			if(pad->name == "n0_135_104" ||
			pad->name == "n0_255_159"||
			pad->name == "n0_30_74")
			clog<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;

			//clog<<"old pad / new pad: "<<*pad<<" "<<*new_pad<<endl;
			pad_ptr->node = new_pad;
			pad_ptr->control_nodes.clear();
			pad_ptr->visit_flag = true;
		}
	}
}

void Circuit::resolve_direct(){
	clock_t t1, t2;
	t1 = clock();
	rebuild_voltage_nets();
	solve();
	//solve_LU_core();
	double max_IR = locate_maxIRdrop();	
	//double max_IRS = locate_special_maxIRdrop();
	clog<<"max_IR by cholmod is: "<<max_IR<<endl;
	t2 = clock();
		clog<<"single solve by cholmod is: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
}

void Circuit::resolve_queue(vector<Node *> pad_set_old){
	clock_t t1, t2;
	t1 = clock();	
	//rebuild_voltage_nets();
	solve_queue(pad_set_old);
	//solve_LU_core();
	double max_IR = locate_maxIRdrop();	
	//double max_IRS = locate_special_maxIRdrop();
	clog<<"max_IR by queue is: "<<max_IR<<endl;
	t2 = clock();
		clog<<"single solve by queue is: "<<1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
}

void Circuit::solve_queue(vector<Node *> pad_set_old){
	double max_diff = 1;
	double eps0 = 1e-10;
	size_t iter = 0;	
	 while(max_diff > eps0){
	//for(size_t i=0;i<5;i++){
		iter++;
		for(size_t j=0;j<nodelist.size()-1;j++)
			nodelist[j]->flag_visited = false;

		max_diff = update_single_iter(pad_set_old);
		clog<<"iter, max_diff: "<<iter<<" "<<max_diff<<endl;
	}
}

double Circuit::update_single_iter(vector<Node *> pad_set_old){
	queue <Node *> q;
	double diff = 1;
	double eps1 = 1e-8;
	Node *nd;
	double max_diff = 0;
	size_t cur_count = 0;
	size_t count_total = 0;
	// inserted varied pads into queue
	initialize_queue(pad_set_old, q);
	// start to expand and update nbr node values
	while(!q.empty() && abs(diff) > eps1){
	//for(size_t j=0;j<120;j++){
		nd = q.front();
		//cur_count++;
		nd->flag_visited = true;
		if(nd->isX() == false){
			// update node value
			diff = update_value(nd);
			if(abs(diff) > max_diff)
				max_diff = abs(diff);
			//if(abs(diff) < eps1)
				//clog<<"find bound node, diff: "<<*nd<<" "<<abs(diff)<<" "<<cur_count<<endl;
		}
		size_t count = update_queue(q, nd);
		count_total += count;
		q.pop();
	}
	//clog<<"total: "<<count_total<<" in the queue. "<<endl;
	while(!q.empty()){
		q.pop();
	}
	return max_diff;
}

double Circuit::update_value(Node *nd){
	if(nd->isX()==true) {
		return 0;
	}
	double h = 0.06;
	double omega = 1;//2-h;
	double V_old=0;
	double V_temp = 0;
	double G = 0;
	Net *net;
	Node *nbr, *na, *nb;
	double sum = 0;
	double current = 0;
	net = NULL;
	nbr = NULL; na = NULL; nb = NULL;

	V_old = nd->value;

	// update nd->value
	for(int i=0;i<6;i++){
		net = nd->nbr[i];
		if(net ==NULL) continue;
		G = 1.0/net->value;
		na = net->ab[0]; nb = net->ab[1];
		if(nd->name == na->name) nbr = nb;
		else	nbr = na;
		if(!nbr->is_ground()){
			sum += G;
			V_temp += G*nbr->value;
		}
	}
	if(nd->nbr[BOTTOM]== NULL) current = 0;
	else	current = -nd->nbr[BOTTOM]->value;
	V_temp += current;
	V_temp /=sum;
	//if(iter ==1 && nd->name==rm_pad->name)
		nd->value  = V_temp;
	//else
		//nd->value = (1-omega)*nd->value + omega*V_temp;
 	
	double V_improve = fabs(nd->value - V_old);
	//cout<<"V_diff: "<<V_improve<<endl;

	return V_improve;
}

// add new nbr nodes into queue
size_t Circuit::update_queue(queue<Node *>&q, Node *nd){
	Net * net; Node *nbr;
	Node *na, *nb;
	size_t count = 0;
	for(int i=0;i<6;i++){
		net = nd->nbr[i];
		if(net==NULL) continue;
		// find the other node, which is nbr
		na = net->ab[0];
		nb = net->ab[1];
		if(nd->name == na->name)
			nbr = nb;
		else	nbr = na;
		if(!nbr->is_ground()&& nbr->flag_visited == false){
			q.push(nbr);
			count ++;
			//cout<<"push q: "<<*nbr<<endl;
			nbr->flag_visited = true;
		}
	}
	return count;
}

void Circuit::initialize_queue(vector<Node *> pad_set_old, queue <Node*> &q){
	Pad *B;
	Node *na;
	Node *nb;
	bool flag = false;
	/*for(size_t i=0;i<nodelist.size()-1;i++){
		if(nodelist[i]->value == 0)
			q.push(nodelist[i]);
	}*/
	// insert moved old pads into queue
	for(size_t i=0;i<pad_set_old.size();i++){
		na = pad_set_old[i];
		flag = false;
		for(size_t j=0;j<pad_set.size();j++){
			B = pad_set[j];
			nb = B->node;
			if(nb->name == na->name){
				flag = true;
				break;
			}
		}
		if(flag == true)
			continue;
		q.push(na);	
		
	}
	for(size_t i=0;i<nodelist.size()-1;i++){
		if(nodelist[i]->name == "n0_48_400" ||
		   nodelist[i]->name == "n0_82_400" ||
		   nodelist[i]->name == "n0_122_400")
			//clog<<"max node: "<<*nodelist[i]<<" "<<nodelist[i]->isX()<<endl;
			q.push(nodelist[i]);
	}

	// insert moved new pads into queue
	for(size_t i=0;i<pad_set.size();i++){
		B = pad_set[i];
		nb = B->node;
		flag = false;
		for(size_t j=0;j<pad_set_old.size();j++){
			na = pad_set_old[j];
			if(nb->name == na->name){
				flag = true;
				break;
			}
		}
		if(flag == true)
			continue;
		q.push(nb);	
		
	}
}

void Circuit::solve_GS(){
	double max_diff = 1;
	int iter = 0;
	double omega=1.3;// 2-0.06;
	Node *max_nd;
	//clog<<"nodelist.size: "<<nodelist.size()-1<<endl;
	while(max_diff >1e-8){// && iter < 500){
		max_diff = 0;
		for(size_t i=0;i<nodelist.size()-1;i++){
			Node *nd = nodelist[i];
			
			if(nd->isX()==true) continue;
			double V_old=0;
			double V_temp = 0;
			double G = 0;
			Net *net;
			Node *nbr, *na, *nb;
			double sum = 0;
			net = NULL;
			nbr = NULL; na = NULL; nb = NULL;

			V_old = nd->value;

			// update nd->value
			for(int i=0;i<6;i++){
				net = nd->nbr[i];
				if(net ==NULL) continue;
				G = 1.0/net->value;
				na = net->ab[0]; nb = net->ab[1];
				if(nd->name == na->name) nbr = nb;
				else	nbr = na;
				if(!nbr->is_ground()){
					sum += G;
					V_temp += G*nbr->value;
				}
			}

			V_temp += -nd->nbr[BOTTOM]->value;//rhs[nd->rid];
			V_temp /=sum;
				nd->value = (1-omega)*nd->value + omega * V_temp;

			double diff = fabs(nd->value - V_old);
			if(diff > max_diff) {
				max_diff = diff;
				max_nd = nd;
			}
		}
		iter++;
		clog<<"iter, diff: "<<iter<<" "<<max_diff<<endl;
	}
}

void Circuit::print_all_control_nodes(){
}
