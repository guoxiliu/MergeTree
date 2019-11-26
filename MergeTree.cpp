#include "MergeTree.h"

// Constructor.
MergeTree::MergeTree(vtkUnstructuredGrid *p){
  usgrid = p;
  group = vector<int>(p->GetNumberOfPoints(), -1);
}

// Find the group id of a given vertex id.
int MergeTree::findGroup(vtkIdType i){
    if(group[i] == -1)
      return i;
    group[i] = findGroup(group[i]);
    return group[i];
}

// Do union of two groups.
void MergeTree::unionGroup(vtkIdType x, vtkIdType y){
    int xset = findGroup(x);
    int yset = findGroup(y);
    if(xset != yset){
      group[xset] = yset;
    }
}

// Sort the scalar values while keeping track of the indices.
vector<vtkIdType> MergeTree::argsort(){
  vector<vtkIdType> indices(usgrid->GetNumberOfPoints());
  iota(indices.begin(), indices.end(), 0);
  vtkDataArray *scalarfield = usgrid->GetPointData()->GetArray(0);
  switch(scalarfield->GetDataType()){
    case VTK_FLOAT:
    {
      float *scalarData = (float *)scalarfield->GetVoidPointer(0);
      sort(indices.begin(), indices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
      break;
    }
    case VTK_DOUBLE:
    {
      double *scalarData = (double *)scalarfield->GetVoidPointer(0);
      sort(indices.begin(), indices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
      break;
    }
    default:
    {
      cout << "Type of scalarfield: " << scalarfield->GetDataType() << ", " << scalarfield->GetDataTypeAsString() << endl;
      break;
    }
  }

  return indices;
}

// Build the merge tree.
int MergeTree::build(){
  vector<vtkIdType> sortedIndices = argsort();
  constructJoin();
  constructSplit();
  mergeJoinSplit();
  return 0;
}

// Construct the join tree.
void MergeTree::constructJoin(){

}

// Construct the split tree.
void MergeTree::constructSplit(){

}

// Merge the split and join tree.
void MergeTree::mergeJoinSplit(){

}