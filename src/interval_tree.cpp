#include <iostream> 
using namespace std;

// Structure to represent an interval 
class Interval 
{   public:

    int low, high; 
}; 

// Structure to represent a node in Interval Search Tree 
class ITNode 
{ 
    public:

    ITNode(Interval value): data(value), left(NULLptr), right(nullptr) {}  // constructor

    Interval getData()const;
    ITNode getLeft()const;
    ITNode getRight()const;

    private:

    Interval *data;  // pointer to node's interval data object
    int max; 
    ITNode *left, *right; // pointer to node's left & right child node objects
}; 

// Simple access functions to retrieve node's interval data
Interval ITNode::getData()const{
    return data;
}

// Simple access functions to retrieve node's left child
Interval ITNode::getLeft()const{
    return left;
}

// Simple access functions to retrieve node's right child
Interval ITNode::getRight()const{
    return right;
}

/////////////////////////////////////////////////////////////////////////////////////////
// IntervalTree class
class IntervalTree{
private:
        ITNode *_root;
        void insert(ITNode *treeNode, int data);
        
public:

        IntervalTree();  // constructor     
        ~IntervalTree();     // destractor

        void insert(Interval data){ insert(_root, data);}       

        void delete() {delete(_root);}

        bool isBalanced(){return isBalanced(_root);        }

};

// A utility function to insert a new Interval Search Tree Node 
// This is similar to BST Insert.  Here the low value of interval 
// is used tomaintain BST property 
void IntervalTree::insert(ITNode *root, Interval data)
{ 
    // Base case: Tree is empty, new node becomes root 
    if (root == NULL) {
        root = new ITNode(data);           
        _root = root;
    }
    else {
    // Get low value of interval at root 
    int l = root->data->low; 
  
    // If root's low value is smaller, then new interval goes to 
    // left subtree 
    if (data.low < l) 
        root->left = insert(root->left, data); 
  
    // Else, new node goes to right subtree. 
    else
        root->right = insert(root->right, data); 
  
    // Update the max value of this ancestor if needed 
    if (root->max < data.high) 
        root->max = data.high; 
    }
} 
  
// A utility function to check if given two intervals overlap 
bool IntervalTree::doOVerlap(Interval i1, Interval i2) 
{ 
    if (i1.low <= i2.high && i2.low <= i1.high) 
        return true; 
    return false; 
} 
  
// The main function that searches a given interval i in a given 
// Interval Tree. 
bool IntervalTree::*overlapSearch(ITNode *root, Interval i) 
{ 
    // Base Case, tree is empty 
    if (root == NULL) return false; 
  
    // If given interval overlaps with root 
    if (doOVerlap(*(root->i), i)) 
        return true; 
  
    // If left child of root is present and max of left child is 
    // greater than or equal to given interval, then i may 
    // overlap with an interval in left subtree 
    if (root->left != NULL && root->left->max >= i.low) 
        return overlapSearch(root->left, i); 
  
    // Else interval can only overlap with right subtree 
    return overlapSearch(root->right, i);
}