#ifndef _CIRCULAR_QUEUE_H
#define _CIRCULAR_QUEUE_H
#include <vector>
#include <iostream>
#include <fstream>
#include "node.h"
class CircularQueue {
public:
    int queueHead_, queueTail_ ;
    int maxSize_ ;
    vector<Node*> queue_ ;

    CircularQueue(int maxSize); 
    ~CircularQueue();
    void insert(Node* node);       
    Node* extractFront(); 
    bool isEmpty();
    void reset();
    friend ostream & operator <<(ostream &os,
const CircularQueue &q);  
} ;

#endif

