#ifndef MERGETREE_H
#define MERGETREE_H

#include <vector>
#include <map>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

// TODO: Super Arc class.
class SuperArc{

};


/**
 * Merge Tree Class.
 * The merge tree is created by combining join tree and split tree, which 
 * captures the evolution of the superlevel or sublevel sets.
 */ 
class MergeTree{

public:
  struct Comp{
    bool operator(int lhs, int rhs) const{
      return scalar[lhs] < scalar[rhs];
    }
  };
public:
  MergeTree(vtkSmartPointer<vtkUnstructuredGrid>);
  void computeMergeTree();//wrap function for compute JT, ST and CT
  void output();

protected:
  void setupData(); // initialize data from vtkUnstructedGrid
  //void orderVertices(); // order vertices by scalarValue

public:

private:
  void computeJoinTree();
  void computeSplitTree();

  int findGroup(int); // Find the group id of a given vertex id.
  void unionGroup(int, int);  // Do union of two groups.

protected:
  vector<int> group;    // Used for Union-Find algorithm
  vtkSmartPointer<vtkUnstructuredGrid> usgrid;    // Unstructed grid
  vector<SuperArc> arcs;    // Save the point set?
  vector<double> scalarValue; // store the scalar function value
  map<int,double,Comp> neighbors; //store points id & its scalar value, naturally ordered
};

#endif
