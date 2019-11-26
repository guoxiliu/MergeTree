#ifndef MERGETREE_H
#define MERGETREE_H

#include <map>
#include <vector>
#include <numeric>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

// TODO: Super Arc class.
class SuperArc{
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
};


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
    vector<int> group;    // Used for Union-Find algorithm
    vtkUnstructuredGrid* usgrid;    // Unstructed grid
    vector<SuperArc> arcs;    // Save the point set?
    // vector<double> scalarValue; // store the scalar function value
    // map<int,vector<int>,Comp> neighbors; //store points id & its neighbors, naturally ordered

  private:
    vector<vtkIdType> argsort(); // Sort the vertex ids based on the scalar values
    int findGroup(vtkIdType); // Find the group id of a given vertex id.
    void unionGroup(vtkIdType, vtkIdType);  // Do union of two groups.
    void constructJoin(); // Construct the join tree.
    void constructSplit(); // Construct the split tree.
    void mergeJoinSplit(); // Merge the split and join tree.
};

#endif
