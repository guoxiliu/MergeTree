#ifndef MERGETREE_H
#define MERGETREE_H

#include "Utils.h"

using namespace std;

// define the tree node in merge tree algorithm 
struct node{
  node(int id):idx(id),numChildren(0),parent(nullptr){}
  int idx;
  int numChildren;
  node* parent;
  vector<node*> children;
};

// define graph node used in constructJoin and construcSplit
struct vNode{
  vNode(){
  }
  node* jNode;
  node* sNode;
};

/**
 * Merge Tree Class.
 * The merge tree is created by combining join tree and split tree, which 
 * captures the evolution of the superlevel or sublevel sets.
 */ 
class MergeTree{
  public:
    MergeTree(vtkUnstructuredGrid *p);
    int build();  // Wrap function for compute JT, ST and CT
    static vector<vtkIdType> argsort(vector<vtkIdType>, vtkUnstructuredGrid*, bool=true);

    vector<vtkIdType> MaximaQuery();// return maxima points in the simplicial 
    vtkIdType ComponentMaximumQuery(vtkIdType&, float&); // return nodeId within the superlevel component  that has maximum scalar

  protected:
    vtkUnstructuredGrid* usgrid;  // Unstructed grid
    vector<vtkIdType> SetMin;
    vector<vtkIdType> SetMax;
    //vector<SuperArc> arcs;  // Save the point set?

  private:
    vector<vtkIdType> argsort();  // Sort the vertex ids based on the scalar values
    vtkSmartPointer<vtkIdList> getConnectedVertices(vtkSmartPointer<vtkUnstructuredGrid>, int);
    float* getScalar();

    void constructJoin(vector<vtkIdType>&);   // Construct the join tree.
    void constructSplit(vector<vtkIdType>&);  // Construct the split tree.
    void mergeJoinSplit(vector<node*>&, vector<node*>&);  // Merge the split and join tree.
  
  private:
    vector<node*> joinTree;   // Represent the join tree
    vector<node*> splitTree;  // Represent the split tree
    vector<vNode*> graph;
    
  public:
    vector<node*> mergeTree;
};


#endif
