#include <iomanip>
#include "global.h"
#include "node.h"
#include "pad.h"
using namespace std;

Pad::Pad(): node(NULL){
	control_nodes.clear();
}