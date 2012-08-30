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
double Circuit::OMEGA = 1.2;
double Circuit::OVERLAP_RATIO = 0;//0.2;
int    Circuit::MODE = 0;
const int MAX_ITERATION = 1;//1000;//1000000;//1000;
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
	clog<<"total special nodes: "<<special_nodes.size()<<endl;
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
	X_BLOCKS =2;
	Y_BLOCKS =1;
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
	//if( MODE == 0 )
	//	solve_IT();
	//else
		solve_LU();
}

// solve Circuit
bool Circuit::solve_IT(){
	// did not find any `X' node
	if( circuit_type == UNKNOWN )
		circuit_type = C4;
	solve_init();

	/*if( replist.size() <= 2*MAX_BLOCK_NODES ){
		clog<<"Replist is small, use direct LU instead."<<endl;
		solve_LU_core();
		return true;
	}
	*/
	select_omega();
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
	solve_init();
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
	mark_special_nodes();
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

// expand the regions for all the special nodes
void Circuit::expand_region(){
	Node *nd;
	for(size_t i=0;i<special_nodes.size();i++){
		nd = special_nodes[i];
		cout<<endl<<"special node: "<<*nd<<endl;
		// mark nodes with region flag
		clock_t t1, t2;
		t1 = clock();
		expand_region_of_a_node(nd);
		t2 = clock();
		//cout<<i<<" expand region. "<<
			//1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
		// find the weighted shortest_path for 
		// all nodes in the region of this node
		/*t1 = clock();
		find_shortest_paths(nd);
		t2 = clock();
		//cout<<i<<" find shortest path cost: "<<
			//1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
		//print_distance(nd);
		t1 = clock();
		map_min_dist_to_pad(nd);*/
		t2 = clock();
		//cout<<i<<" map mpad cost: "<<
			//1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
		// clear region, traversal, visit flag
		t1 = clock();
		clear_flags();
		t2 = clock();
		//cout<<i<<" clear flags cost: "<<
			//1.0*(t2-t1)/CLOCKS_PER_SEC<<endl;
	}

	//print_pad_map();
}

// expand the region for each node, covering 10 pads
void Circuit::expand_region_of_a_node(Node *nds){
	// stop when reaching pad_number
	int pad_number = 5;//pad_set.size();//2; 	
	int count = 0;
	double dist = 0;
	queue<Node*> q;
	q.push(nds);
	nds->region_flag = true;
	if(nds->isX())
		count ++;
	//cout<<"nds: "<<*nds<<endl;
	while(count < pad_number || !q.empty()){
		Node * nd = q.front();
		// add neighboring nodes into queue
		update_queue(nds, q, nd, count, dist);
		q.pop();
		if(count == pad_number){
			while(!q.empty()){
				q.pop();
			}
			break;
		}
	}
}

void Circuit::extract_pads(int pad_number){
	vector<Node*> pair_first;
	vector<double> pair_second;
	pair<Node*, double> pair_nd;
	//int pad_number = 5;
	double distance = 0;
	double delta_x = 0;
	double delta_y = 0;
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
			delta_x=(ptr->pt.x-nd->pt.x);
			delta_y=(ptr->pt.y-nd->pt.y);
			delta_x *= delta_x;
			delta_y *= delta_y;
			distance = sqrt(delta_x + delta_y);

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

void Circuit::extract_min_max_pads(vector<double> ref_drop_vec){
	double min = VDD-ref_drop_vec[0];
	double max = 0;	
	size_t min_index=0;
	size_t max_index=0;
	vector<Node*> min_pads;
	vector<Node*> max_pads;
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
		if(VDD - ref_drop_vec[j]<max*0.2){
			min_pads.push_back(pad_set[j]->node);
		}
	//}
	//for(size_t j=0;j<ref_drop_vec.size();j++){
		if(VDD-ref_drop_vec[j]>max*0.9){
			max_pads.push_back(pad_set[j]->node);
		}
	}

	Node *new_pad;
	Pad * pad_ptr;

	// set nd into the weighted center
	// start to map min pads into max pads locations
	for(size_t j=0;j<min_pads.size();j++){
		size_t i = j % max_pads.size();
		Node * nd = max_pads[i];
		//clog<<"min_pads, max_pads: "<<*min_pads[j]<<" "<<*max_pads[i]<<endl;
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				nd->name){	
				pad_ptr = pad_set[k];
				double ref_drop_value = ref_drop_vec[k];
				new_pad = pad_projection(pad_ptr, min_pads[j]);
				pad_ptr->node->enableX();
				pad_ptr->node->value = VDD;	
				break;
			}
		}
		for(size_t k=0;k<pad_set.size();k++){
			if(pad_set[k]->node->name == 
				min_pads[j]->name){
				//clog<<"min_pad, new: "<<*pad_set[k]->node<<" "<<*new_pad<<endl;
				pad_set[k]->node->disableX();
				pad_set[k]->node->value = 0;
				pad_set[k]->node = new_pad;
				// already taken care of
				pad_set[k]->control_nodes.clear();
				break;
			}
		}

	}
	min_pads.clear();
	max_pads.clear();
}

// tune 50% nodes with the small IR drops
void Circuit::update_pad_control_nodes(vector<double> & ref_drop_value, size_t iter){
	Pad *pad_ptr;
	Node *pad;
	map<Node*, double>::iterator it;
	Node *nd;
	double weight=0;
	vector<double> drop_vec;
	ref_drop_value.resize(pad_set.size());
	for(size_t i=0;i<pad_set.size();i++){
		if(pad_set[i]->control_nodes.size()==0)
			continue;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		drop_vec.clear();
		for(it = pad_ptr->control_nodes.begin();
		    it != pad_ptr->control_nodes.end();
		    it++){
			nd = it->first;
			weight = nd->value;
			if(weight <0)
				weight *=10;

			pad_ptr->control_nodes[nd] = weight;
			drop_vec.push_back(nd->value); 
		}
				// sort control_nodes in ascending order
		sort(drop_vec.begin(), drop_vec.end(),
			compare_values);
		//if(pad_set[i]->node->name=="n0_350_50")
			//for(size_t j=0;j<drop_vec.size();j++)
				//clog<<"drop_vec: "<<drop_vec[j]<<endl;

		double middle_value = drop_vec[drop_vec.size()/2];
		//if((VDD-middle_value) / max_IRdrop <=0.5)
			ref_drop_value[i] = middle_value;
		//else
			//ref_drop_value[i] = drop_vec[drop_vec.size()-1];
		
		//if(pad_set[i]->node->name=="n0_350_50")
			//clog<<"middle value: "<<middle_value<<endl; 
	}
	drop_vec.clear();
}

void Circuit::update_queue(Node *nds, queue<Node*> &q, Node *nd, int &count, double &dist){
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
}

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

// find weighted paths from source node nds to all nodes
// in the region
void Circuit::find_shortest_paths(Node *nds){
	Node *nds_old;
	Node *nds_new = NULL;
	nds_old = nds;
	nds_old->distance = 0;
	nds_old->visit_flag = true;
	nds_old->traverse_flag = true;
	// stores all front end nodes
	vector<Node*> front_nodes;
	do{
		nds_new = update_distance(nds_old, 
			front_nodes);	
		nds_old = nds_new;
		//cout<<endl;
	}while(nds_new !=NULL);
}

Node* Circuit::update_distance(Node *nd, vector<Node *> &front_nodes){
	Net * net; Node *nbr;
	Node *na, *nb;
	for(int i=0;i<6;i++){
		net = nd->nbr[i];
		if(net==NULL) continue;
		// find neighboring nodes
		na = net->ab[0];
		nb = net->ab[1];
		if(nd->name == na->name)
			nbr = nb;
		else	nbr = na;
		if(nbr->visit_flag == true) continue;
		if(nbr->region_flag == false) continue;
		if(!nbr->is_ground()&& !nbr->visit_flag)
{
			// assign initial weight
			double distance = nd->distance 
				+ net->value;
			if(nbr->distance == -1 || 
			   distance < nbr->distance){
				//cout<<"nbr, distance, calc: "<<*nbr<<" "<<nbr->distance<<" "<<distance<<endl;
				nbr->distance = distance;
			}

			// if node is not traversed yet	
			if(!nbr->traverse_flag){
				front_nodes.push_back(nbr);
			}
		}
		nbr->traverse_flag = true;
	}
	// select the least distance node from 
	// front_nodes
	Node *nds_new = min_dist_front_node(front_nodes);
	if(nds_new !=NULL){
		//cout<<"new min dist node: "<<*nds_new<<endl;
		nds_new->visit_flag = true;
	}
	return nds_new;
}

Node* Circuit::min_dist_front_node(vector<Node*> front_nodes){
	Node *nd;
	Node *min_dist_nd;
	min_dist_nd = NULL;
	double min_dist = -1.0;
	for(size_t i=0;i<front_nodes.size();i++){
		nd = front_nodes[i];
		if(nd->visit_flag == true) continue;
		//cout<<"front_node: "<<*nd<<" "<<nd->distance<<endl;
		if(min_dist == -1.0){
			min_dist = nd->distance;
			min_dist_nd = nd;
		}			
		else{
			if(nd->distance < min_dist){
				min_dist = nd->distance;
				min_dist_nd = nd;
			}
		}
	}
	return min_dist_nd;
}

void Circuit::print_distance(Node *nd){
	cout<<endl<<"special node: "<<*nd<<endl;
	for(size_t j=0;j<nodelist.size()-1;j++){
		if(nodelist[j]->region_flag==true)
			cout<<"nd, dist: "<<nodelist[j]->name<<" "<<nodelist[j]->distance<<endl;
	}	
}

// map the shortest distance from node to pad
void Circuit::map_min_dist_to_pad(Node *nds){
	Pad *nd;
	Node *ndt;
	pair<Node*, double> pad_pair;
	//cout<<"nds: "<<*nds<<endl;
	for(size_t i=0;i<pad_set.size();i++){
		nd = pad_set[i];
		ndt = nd->node;
		if(!ndt->region_flag) continue;
		//cout<<"pad belongs to region: "<<*ndt<<endl;
		pad_pair.first = nds;
		pad_pair.second = ndt->distance;
		nd->control_nodes.insert(pad_pair);
	}	
}

void Circuit::print_pad_map(){
	Pad *nd;
	pair<Node*, double> pad_pair;
	map<Node*, double>::iterator it;
	for(size_t i=0;i<pad_set.size();i++){
		nd = pad_set[i];
		if((nd->node->name != "n0_15_20" && 
			nd->node->name != "n0_11_78"))// &&
			//if(nd->node->name !="n0_52_60")
			//if(nd->node->name !="n0_74_91") 
			//if(nd->node->name !="n0_77_26")
			continue;
		//cout<<endl<<"pad_node: "<<*nd->node<<endl;
		//cout<<"v"<<i<<" "<<nd->node->name<<" "<<"0 "<<VDD<<endl;
		clog<<"control_nodes.size()"<<nd->control_nodes.size()<<endl;
		for(it = nd->control_nodes.begin();it != nd->control_nodes.end();it++){
			//cout<<"control node: "<<*it->first<<" "<<it->second<<endl;
			printf("%ld %ld  %.5e\n", it->first->pt.y+1, it->first->pt.x+1, it->first->value);
			
		}
	}
}

void Circuit::clear_flags(){
	Node *nd;
	for(size_t i=0;i<nodelist.size()-1;i++){
		nd = nodelist[i];
		nd->region_flag = false;
		nd->traverse_flag = false;
		nd->visit_flag = false;
		nd->distance = -1;	
	}
}

void Circuit::clear_pad_control_nodes(){
	Node *nd;
	for(size_t i=0;i<pad_set.size();i++){
		pad_set[i]->control_nodes.clear();
	}
}

void Circuit::relocate_pads(){
	vector<Node*> pad_set_old;
	double dist = 0;
	double new_dist = 0;
	//for(size_t i=0;i<6;i++){
	pad_set_old.resize(pad_set.size());
	assign_pad_set(pad_set_old);
	// store a original copy of the pad set
	
	// build up the map for nodes in PAD layer
	build_map_node_pt();

	vector<double> ref_drop_vec;
	//print_pad_set();
	for(size_t i=0;i<6;i++){
		int pad_number = 1;
		origin_pad_set.resize(pad_set.size());
		assign_pad_set(origin_pad_set);
		//clog<<"i, pad set size: "<<i<<" "<<pad_set.size()<<endl;
		extract_pads(pad_number);
		update_pad_control_nodes(ref_drop_vec, i);
			
		dist = update_pad_pos_all(ref_drop_vec);
		//if(i==0)
			extract_min_max_pads(ref_drop_vec);
		assign_pad_set(pad_set_old);
		project_pads();

		//print_pad_set();
		//dist = update_pad_pos();
		//clog<<"finish first update pad pos. "<<endl;
		//project_pads();
		//clog<<"finish first project pads. "<<endl;

		rebuild_voltage_nets();
		solve_LU_core();

		double max_IR = locate_maxIRdrop();	
		double max_IRS = locate_special_maxIRdrop();
		clog<<"max_IR is: "<<max_IR<<endl;
		clog<<"max_IRS is: "<<max_IRS<<endl<<endl;
	}
	ref_drop_vec.clear();
	map_node_pt.clear();
	print_pad_set();
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
	//double ref_drop_value=0;

	//for(size_t i=0;i<pad_set.size();i++){
		//if(pad_set[i]->control_nodes.size()==0)
			//continue;

		//ref_drop_value = ref_drop_vec[i];

		double sum_weight = 0;
		double weighted_x =0;
		double weighted_y =0;
		size_t count = 0;
		pad_ptr = pad_set[i];
		pad = pad_ptr->node;
		//cout<<endl<<"i, pad_node: "<<i<<" "<<*pad<<endl;
		for(it = pad_ptr->control_nodes.begin();
		    it != pad_ptr->control_nodes.end();
		    it++){
			if(it->second > ref_drop_value)
				continue;
			if((pad->name == "n0_15_20" || 
			pad->name == "n0_11_78")){
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

			pad_ptr->newx = pad_newx;
			pad_ptr->newy = pad_newy;
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
		if(pad->name == "n0_15_20" ||
			//pad->name == "n0_50_150"||
			pad->name == "n0_11_78")
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

	nd = pad->node;
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
					nd_new->enableX();
					nd_new->value = VDD;
					return_flag = true;
					//clog<<"break. "<<endl;
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
