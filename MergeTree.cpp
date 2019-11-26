#include "MergeTree.h"


MergeTree::MergeTree(vtkSmartPointer<vtkUnstructuredGrid> p):usgrid(p),group(vector<int>(p->GetNumberOfPoints(), -1)){}

int MergeTree::findGroup(int i){
    if(group[i] == -1)
      return i;
    group[i] = findGroup(group[i]);
    return group[i];
}

void MergeTree::unionGroup(int x, int y){
    int xset = findGroup(x);
    int yset = findGroup(y);
    if(xset != yset){
      group[xset] = yset;
    }
}

void MergeTree::setupData(){

}

