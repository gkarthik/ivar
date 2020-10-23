#include <iostream> 
using namespace std;

#ifndef parse_gff
#define parse_gff

// Structure to represent an interval 
class Interval 
{   public:
    Interval(int val1, int val2): low(std::min(val1, val2)), high(std::max(val1, val2)) {}  // constructor
    int low, high; 
}; 
// Structure to represent a node in Interval Search Tree 
class ITNode 
{ 
    /*
    public:
    ITNode(Interval *values): data(value), left(nullptr), right(nullptr) {}  // constructor
    int max; 
    // Getters - access member functions
    Interval getData()const;
    ITNode getLeft()const;
    ITNode getRight()const;

    // Setters - access member functions
    void setLeft(ITNode *node);
    void setRight(ITNode *node);
    */
    public:
    ITNode(Interval value): data(new Interval(value)), left(nullptr), right(nullptr), max(value.high) {}  // constructor
    int max; 
    Interval *data;  // pointer to node's interval data object
    ITNode *left, *right; // pointer to node's left & right child node objects
}; 

/////////////////////////////////////////////////////////////////////////////////////////
// IntervalTree class
class IntervalTree{
private:
        ITNode *_root;
        void insert(ITNode *root, Interval data);
        bool overlapSearch(ITNode *root, Interval data);
        void inOrder(ITNode * root);
        
public:

        IntervalTree();  // constructor     



        void insert(Interval data){ insert(_root, data);}

        bool overlapSearch(Interval data){ return overlapSearch(_root, data);} 

        void inOrder() {inOrder(_root);}

};

#endif
