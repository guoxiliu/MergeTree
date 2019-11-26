#ifndef MERGETREE_H
#define MERGETREE_H

#include <vector>
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
  MergeTree(vtkSmartPointer<vtkUnstructuredGrid> p){
    usgrid = p;
    group = vector<int>(p->GetNumberOfPoints(), -1);
  }

  // TODO: Implement the topological queries.

  // TODO: Maxima query

  // TODO: Component maximum query


protected:
  vector<int> group;    // Used for Union-Find algorithm
  vtkSmartPointer<vtkUnstructuredGrid> usgrid;    // Unstructed grid
  vector<SuperArc> arcs;    // Save the point set?

private:
  // Find the group id of a given vertex id.
  int findGroup(int i){
    if(group[i] == -1)
      return i;
    return findGroup(group[i]);
  }

  // Do union of two groups.
  void unionGroup(int x, int y){
    int xset = findGroup(x);
    int yset = findGroup(y);
    if(xset != yset){
      group[xset] = yset;
    }
  }
};

#endif
