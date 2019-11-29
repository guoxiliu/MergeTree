#include "MergeTree.h"


// Constructor.
MergeTree::MergeTree(vtkUnstructuredGrid *p){
  usgrid = p;
  Set = vector<int>(p->GetNumberOfPoints());
  SetMax = vector<int>(p->GetNumberOfPoints());
  SetMin = vector<int>(p->GetNumberOfPoints());
  graph = vector<vNode*>(p->GetNumberOfPoints(), new vNode());
}

// Find the Set id of a given vertex id.
int MergeTree::findSet(vtkIdType i){
    if(Set[i] == i)
      return i;
    Set[i] = findSet(Set[i]);
    return Set[i];
}

int MergeTree::findSetMax(vtkIdType i){
  if(SetMax[i] == i)
      return i;
    SetMax[i] = findSet(SetMax[i]);
    return SetMax[i];
}

int MergeTree::findSetMin(vtkIdType i){
  if(SetMin[i] == i)
      return i;
    SetMin[i] = findSet(SetMin[i]);
    return SetMin[i];
}

// Do union of two Sets.
void MergeTree::unionSet(vtkIdType x, vtkIdType y){
    int xset = findSet(x);
    int yset = findSet(y);
    if(xset != yset){
      Set[xset] = yset;  
    }

    // make the root with higher scalar value
    xset = findSetMax(x);
    yset = findSetMax(y);
    if(xset != yset){
      if(xset < yset){
        SetMax[xset] = yset;
      }else{
        SetMax[yset] = xset;
      } 
    }

    // make the root with lower scalar value
    xset = findSetMin(x);
    yset = findSetMin(y);
    if(xset != yset){
      if(xset < yset){
        SetMax[yset] = xset;
      }else{
        SetMax[xset] = yset;
      } 
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
  constructJoin(sortedIndices);
  constructSplit(sortedIndices);
  mergeJoinSplit(joinTree, splitTree);
  return 0;
}

// Construct the join tree.
void MergeTree::constructJoin(vector<vtkIdType>& sortedIndices){

  // record [vtkId] = pos in sortedindices
  vector<int> sortedIds(sortedIndices.size());
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  }
  iota(Set.begin(), Set.end(), 0);
  iota(SetMax.begin(), SetMax.end(), 0);
  iota(SetMin.begin(), SetMin.end(), 0);
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    node *ai = new node(i);
    joinTree.push_back(ai);
    graph[i]->jNode = ai;
    // get the neighbors of ai
    vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(usgrid, sortedIndices[i]);
    for(vtkIdType adj = 0; adj < connectedVertices->GetNumberOfIds();++adj){
      // index of adj in sortedIndices
      vtkIdType j = sortedIds[connectedVertices->GetId(adj)];
      if(j < i && findSet(i) != findSet(j)){
          int k = findSetMax(j);
          graph[k]->jNode->parent = ai;
          ai->numChildren += 1;
          unionSet(i,j);
      }
    }
  }
}

// Construct the split tree.
void MergeTree::constructSplit(vector<vtkIdType>& sortedIndices){
  // record [vtkId] = pos in sortedindices
  vector<int> sortedIds(sortedIndices.size());
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  }

  iota(Set.begin(), Set.end(), 0);
  iota(SetMax.begin(), SetMax.end(), 0);
  iota(SetMin.begin(), SetMin.end(), 0);
  for (int i = sortedIndices.size()-1; i >= 0; --i){
    node *bi = new node(i);
    splitTree.push_back(bi);
    graph[i]->sNode = bi;
    // get the neighbors of bi
    vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(usgrid, sortedIndices[i]);
    for(vtkIdType adj = 0; adj < connectedVertices->GetNumberOfIds(); ++adj){
      // index of adj in sortedIndices
      vtkIdType j = sortedIndices[connectedVertices->GetId(adj)];
      if(j > i && findSet(i) != findSet(j)){
        int k = findSetMin(j);
        graph[k]->sNode->parent = bi;
        bi->numChildren += 1;
        unionSet(i,j);
      }
    }
  }
  // reverse the order of split nodes from 0 to n-1
  reverse(splitTree.begin(), splitTree.end());
}

// Merge the split and join tree.
void MergeTree::mergeJoinSplit(vector<node*>& joinTree, vector<node*>& splitTree){
  queue<int> Q;
  for(int i = 0; i < usgrid->GetNumberOfPoints(); ++i){
    node *ci = new node(i);
    mergeTree.push_back(ci);
    if(joinTree[i]->numChildren + splitTree[i]->numChildren == 1){
      Q.push(i);
    }
  }

  while (Q.size() > 1){
    int i = Q.back();
    Q.pop();
    int k;
    if(joinTree[i]->numChildren == 0){
      k = joinTree[i]->parent->idx;
      mergeTree[i]->parent = mergeTree[k]; 
    }else{
      k = splitTree[i]->parent->idx;
      mergeTree[i]->parent = mergeTree[k];
    }
    // delete ai from join Tree
    if(joinTree[i]->parent){
      joinTree[i]->parent->numChildren -= 1;
    }
    joinTree[i] = nullptr;
    
    //delete bi from split tree
    if(splitTree[i]->parent){
      splitTree[i]->parent->numChildren -= 1;
    }
    splitTree[i] = nullptr;

    if(joinTree[k]->numChildren + splitTree[k]->numChildren == 1){
      Q.push(k);
    }
    cout<< "MergeTree is built" << endl;
  }
}


vtkSmartPointer<vtkIdList> MergeTree::getConnectedVertices(vtkSmartPointer<vtkUnstructuredGrid> usgrid, int id){
  vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();

  // get all cells that vertex 'id' is a part of 
  vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
  usgrid->GetPointCells(id,cellIdList);

  set<vtkIdType> neighbors;
  for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds();++i){
    // get the points of each cell
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    usgrid->GetCellPoints(cellIdList->GetId(i),pointIdList);

    for(vtkIdType j = 0; j <pointIdList->GetNumberOfIds(); ++j){
      neighbors.insert(pointIdList->GetId(j));   
    }
  }
  neighbors.erase(id);
  for(auto it = neighbors.begin(); it != neighbors.end(); ++it){
    connectedVertices->InsertNextId(*it);
  }
  return connectedVertices;
}