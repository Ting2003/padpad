#include <iomanip>
#include "global.h"
#include "node.h"
#include "pad.h"
using namespace std;

Pad::Pad(): node(NULL){
	control_nodes.clear();
	nbrs.clear();
	newx = 0;
	newy = 0;
	visit_flag = false;
	fix_flag = false;
}

Pad::~Pad(){
	control_nodes.clear();
	nbrs.clear();
}
