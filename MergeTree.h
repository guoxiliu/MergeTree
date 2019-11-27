#ifndef MERGETREE_H
#define MERGETREE_H

#include <map>
#include <vector>
#include <queue>
#include <set>
#include <array>
#include <algorithm>
#include <numeric>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

// TODO: Super Arc class.
/* class SuperArc{
  public: 
    SuperArc(){ vertexIds = vector<vtkIdType>(); }
    inline vector<vtkIdType> getVertexIds(){
      return vertexIds;
    }
    inline void addVertex(vtkIdType id){
      vertexIds.push_back(id);
    }
  
  private:
    vector<vtkIdType> vertexIds; 
}; */

//define the tree node in merge tree algorithm 
struct node{
  node(int id):idx(id),numChildren(0),parent(nullptr){}
  int idx;
  int numChildren;
  node* parent;
};

// define graph node used in constructJoin and construcSplit
struct vNode{
  vNode(){}
  node* jNode;
  node* sNode;
}
/**
 * Merge Tree Class.
 * The merge tree is created by combining join tree and split tree, which 
 * captures the evolution of the superlevel or sublevel sets.
 */ 
class MergeTree{
  public:
    // struct Comp{
    //   bool operator()(int lhs, int rhs) const{
    //     return scalarValue[lhs] < scalarValue[rhs];
    //   }
    // };
    MergeTree(vtkUnstructuredGrid *p);
    int build();    // Wrap function for compute JT, ST and CT
    void output();

  protected:
    vector<int> Set;    // Used for Union-Find algorithm
    vtkUnstructuredGrid* usgrid;    // Unstructed grid
    vector<int> SetMin;
    vector<int> SetMax
    //vector<SuperArc> arcs;    // Save the point set?
    // vector<double> scalarValue; // store the scalar function value
    // map<int,vector<int>,Comp> neighbors; //store points id & its neighbors, naturally ordered

  private:
    vector<vtkIdType> argsort(); // Sort the vertex ids based on the scalar values
    int findSet(vtkIdType); // Find the group id of a given vertex id.
    void unionSet(vtkIdType, vtkIdType);  // Do union of two groups.
    int findSetMax(vtkIdType); // get the maximum element of a set
    int findSetMin(vtkIdType); // get the minimum element of a set

    vtkSmartPointer<vtkIdList> getConnectedVertices(vtkSmartPointer<vtkUnstructuredGrid>, int);

    void constructJoin(vector<vtkIdType>&); // Construct the join tree.
    void constructSplit(vector<vtkIdType>&); // Construct the split tree.
    void mergeJoinSplit(vector<node*>&, vector<node*>&); // Merge the split and join tree.
  
  private:
    vector<node*> joinTree;//represent the join tree
    vector<node*> splitTree;//represent the split tree
    vector<vNode*> graph;
    
  public:
    vector<node*> mergeTree;
};

#endif
