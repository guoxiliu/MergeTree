#ifndef MERGETREE_H
#define MERGETREE_H

#include "Utils.h"

using namespace std;


// struct superArc{
//   superArc *parent;
//   vector<superArc*> children;
//   list<vtkIdType> vertexList;

//   superArc():parent(nullptr){}
// };

struct node{
  vtkIdType vtkIdx;
  node *parent;
  vector<node *> children;
  node(vtkIdType id):vtkIdx(id), parent(nullptr), children(vector<node *>()){}
};


/**
 * Merge Tree Class.
 * The merge tree is created by combining join tree and split tree, which 
 * captures the evolution of the superlevel or sublevel sets.
 */ 
class MergeTree{
  public:
    MergeTree(vtkImageData*);
    MergeTree(vtkImageData*, vector<vtkIdType>);
    int build();  // Wrap function for compute JT, ST and CT
    int build(vector<vtkIdType>&);

    vector<vtkIdType> MaximaQuery(const set<pair<vtkIdType, vtkIdType>> &);   // return all local maxima in the simplicial complex
    vtkIdType ComponentMaximumQuery(vtkIdType&, float&);  // return vertexId within the superlevel component that has maximum scalar function value
  
  protected:
    vtkImageData* sgrid;  // Unstructed grid

  private:
    vector<vtkIdType> vertexList;
    void constructJoin(vector<vtkIdType>&);   // Construct the join tree.
    void constructSplit(vector<vtkIdType>&);  // Construct the split tree.
    void mergeJoinSplit(vector<node*>&, vector<node*>&);  // Merge the split and join tree.
    vector<vtkIdType> getLowerLinks(vtkIdType);
    vector<vtkIdType> getUpperLinks(vtkIdType);
  
    vector<node*> joinTree;   // Represent the join tree
    vector<node*> splitTree;  // Represent the split tree
    vector<node*> mergeTree;
    
   
};


#endif
