#ifndef MERGETREE_H
#define MERGETREE_H

#include "Utils.h"

using namespace std;

// define the tree node in merge tree algorithm 
struct node{
  node(long int id, vtkIdType vtkId):idx(id),vtkIdx(vtkId),numChildren(0),parent(nullptr){}
  long int idx;
  vtkIdType vtkIdx;
  int numChildren;
  node* parent;
  vector<node*> children;
};

// define graph node used in constructJoin and construcSplit
/* struct vNode{
  vNode(){
  }
  node* jNode;
  node* sNode;
}; */

/**
 * Merge Tree Class.
 * The merge tree is created by combining join tree and split tree, which 
 * captures the evolution of the superlevel or sublevel sets.
 */ 
class MergeTree{
  public:
    MergeTree(vtkImageData *p);
    int build();  // Wrap function for compute JT, ST and CT
    int build(vector<vtkIdType>&);

    vector<vtkIdType> MaximaQuery(const set<pair<vtkIdType, vtkIdType>> &);   // return all local maxima in the simplicial complex
    vtkIdType ComponentMaximumQuery(vtkIdType&, float&);  // return vertexId within the superlevel component that has maximum scalar function value
  
  protected:
    vtkImageData* sgrid;  // Unstructed grid
    vector<vtkIdType> SetMin;
    vector<vtkIdType> SetMax;

  private:
    void constructJoin(vector<vtkIdType>&);   // Construct the join tree.
    void constructSplit(vector<vtkIdType>&);  // Construct the split tree.
    void mergeJoinSplit(vector<node*>&, vector<node*>&);  // Merge the split and join tree.
    vtkSmartPointer<vtkIdList> getConnectedVertices(vtkIdType, const vector<vtkIdType>&);
    void getEdgeList(vector<set<vtkIdType>>&)
  
  private:
    vector<node*> joinTree;   // Represent the join tree
    vector<node*> splitTree;  // Represent the split tree
    vector<set<vtkIdType>> edgeList;
    //vector<vNode*> graph;
    
  public:
    vector<node*> mergeTree;
};


#endif
