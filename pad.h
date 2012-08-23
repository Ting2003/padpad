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
	
	Node * node;
	// record the controlled nodes and weighting
	map<Node*, double> control_nodes;
};
#endif
