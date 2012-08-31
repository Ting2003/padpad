#ifndef __PAD_H__
#define __PAD_H__

#include <string>
#include <algorithm>
#include "global.h"
#include "point.h"
#include "node.h"
#include <vector>
#include <iostream>
#include <map>
//#include "block.h"
using namespace std;

class Pad{
public:
	Pad();
	~Pad();
	
	Node * node;
	// record the controlled nodes and weighting
	map<Node*, double> control_nodes;
	vector<Pad*> nbrs;
	double newx;
	double newy;
	bool visit_flag;
	bool fix_flag;
};
#endif
