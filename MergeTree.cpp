#include "MergeTree.h"

// Constructor.
MergeTree::MergeTree(vtkUnstructuredGrid *p){
  usgrid = p;
  //Set = vector<int>(p->GetNumberOfPoints());
  SetMax = vector<vtkIdType>(p->GetNumberOfPoints());
  SetMin = vector<vtkIdType>(p->GetNumberOfPoints());
  //graph = vector<vNode*>(p->GetNumberOfPoints()， new vNode());
}

// Sort the scalar values while keeping track of the indices.
vector<vtkIdType> MergeTree::argsort(){
  vector<vtkIdType> indices(usgrid->GetNumberOfPoints());
  iota(indices.begin(), indices.end(), 0);
  return argsort(indices, usgrid);
}

// Sort the scalar values while keeping track of the indices.
vector<vtkIdType> MergeTree::argsort(vector<vtkIdType> vertexSet,vtkUnstructuredGrid* usgrid, bool increasing){
  vtkDataArray *scalarfield = usgrid->GetPointData()->GetArray(0);
  switch(scalarfield->GetDataType()){
    case VTK_FLOAT:
    {
      float *scalarData = (float *)scalarfield->GetVoidPointer(0);
      //float *scalarData = (float*) getScalar();
      if(increasing){
        sort(vertexSet.begin(), vertexSet.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
      }else{
        sort(vertexSet.begin(), vertexSet.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] > scalarData[i2];});
      }
      break;
    }
    case VTK_DOUBLE:
    {
      double *scalarData = (double *)scalarfield->GetVoidPointer(0);
      //double* scalarData = (double*) getScalar();
      if(increasing){
        sort(vertexSet.begin(), vertexSet.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
      }else{
        sort(vertexSet.begin(), vertexSet.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] > scalarData[i2];});
      }
      break;
    }
    default:
    {
      cout << "Type of scalarfield: " << scalarfield->GetDataType() << ", " << scalarfield->GetDataTypeAsString() << endl;
      break;
    }
  }

  return vertexSet;
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
  //iota(Set.begin(), Set.end(), 0);
  iota(SetMax.begin(), SetMax.end(), 0);
  //iota(SetMin.begin(), SetMin.end(), 0);
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    node *ai = new node(i);
    joinTree.push_back(ai);
    graph.push_back(new vNode());
    graph[i]->jNode = ai;
    // get the neighbors of ai
    vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(usgrid, sortedIndices[i]);
    for(vtkIdType adj = 0; adj < connectedVertices->GetNumberOfIds();++adj){
      // index of adj in sortedIndices
      vtkIdType j = sortedIds[connectedVertices->GetId(adj)];
      if(j < i && findSet(SetMax, i) != findSet(SetMax, j)){
        int k = findSet(SetMax, j);
        graph[k]->jNode->parent = ai;
        ai->children.push_back(graph[k]->jNode);
        ai->numChildren += 1;
        unionSet(SetMax, i, j);
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

  // iota(Set.begin(), Set.end(), 0);
  // iota(SetMax.begin(), SetMax.end(), 0);
  iota(SetMin.begin(), SetMin.end(), 0);
  for (int i = sortedIndices.size()-1; i >= 0; --i){
    node *bi = new node(i);
    splitTree.push_back(bi);
    graph[i]->sNode = bi;
    // get the neighbors of bi
    vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(usgrid, sortedIndices[i]);
    for(vtkIdType adj = 0; adj < connectedVertices->GetNumberOfIds(); ++adj){
      // index of adj in sortedIndices
      vtkIdType j = sortedIds[connectedVertices->GetId(adj)];
      if(j > i && findSet(SetMin, i) != findSet(SetMin, j)){
        int k = findSet(SetMin, j);
        graph[k]->sNode->parent = bi;
        bi->children.push_back(graph[k]->sNode);
        bi->numChildren += 1;
        unionSet(SetMin, i,j);
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
    int i = Q.front();
    Q.pop();
    int k;
    if(joinTree[i]->numChildren == 0){
      k = joinTree[i]->parent->idx;
      mergeTree[i]->parent = mergeTree[k]; 
      mergeTree[k]->numChildren += 1;
      mergeTree[k]->children.push_back(mergeTree[i]);
      // delete ai from join Tree
      //if(joinTree[i]->parent){
      auto it = find(joinTree[i]->parent->children.begin(), joinTree[i]->parent->children.end(),joinTree[i]);
      joinTree[i]->parent->children.erase(it);
      joinTree[i]->parent->numChildren -= 1;
      //}
      joinTree[i] = nullptr;

      // connect bi's parent with bi's children

      //if bi is the root of split tree, just make it's children's parent to null
      if(!splitTree[i]->parent){
        for(auto child : splitTree[i]->children){
          child->parent = nullptr;
        }
      }else{
        it = find(splitTree[i]->parent->children.begin(), splitTree[i]->parent->children.end(), splitTree[i]);
        splitTree[i]->parent->children.erase(it);
        splitTree[i]->parent->numChildren -= 1;

        for(auto child : splitTree[i]->children){
          child->parent = splitTree[i]->parent;
          splitTree[i]->parent->children.push_back(child);
          splitTree[i]->parent->numChildren += 1;
        }
      }
      splitTree[i] = nullptr;
    }else{
      k = splitTree[i]->parent->idx;
      mergeTree[i]->parent = mergeTree[k];
      mergeTree[k]->numChildren += 1;
      mergeTree[k]->children.push_back(mergeTree[i]);
      //delete bi from split tree
      //if(splitTree[i]->parent){
      auto it = find(splitTree[i]->parent->children.begin(), splitTree[i]->parent->children.end(), splitTree[i]);
      splitTree[i]->parent->children.erase(it);
      splitTree[i]->parent->numChildren -= 1;
      //}
      splitTree[i] = nullptr;

      // connect ai's parent with ai's children
      // if ai is the root of split tree, just make it's children's parent to null
       if(!joinTree[i]->parent){
        for(auto child : joinTree[i]->children){
          child->parent = nullptr;
        }
      }else{
        it = find(joinTree[i]->parent->children.begin(), joinTree[i]->parent->children.end(),joinTree[i]);
        joinTree[i]->parent->children.erase(it);
        joinTree[i]->parent->numChildren -= 1;
      
        for(auto child : joinTree[i]->children){
          child->parent = joinTree[i]->parent;
          joinTree[i]->parent->children.push_back(child);
          joinTree[i]->parent->numChildren += 1;
        }
      }
      joinTree[i] = nullptr;
    }
     
    if(joinTree[k]->numChildren + splitTree[k]->numChildren == 1){
      Q.push(k);
    }
  }
  cout<< "MergeTree is built" << endl;
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

// get the nodeId (id in .vtu) of local maximum
vector<vtkIdType> MergeTree::MaximaQuery(){
  vector<vtkIdType> sortedIndices = argsort();
  vector<int> sortedIds(sortedIndices.size());
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  }
  vector<vtkIdType> maxima;
  //iterate mergeTree to find local maximum
  for(auto node:mergeTree){
    if(node->numChildren == 0 && node->parent->idx < node->idx){
      maxima.push_back(sortedIds[node->idx]);
    }
  }
  return maxima;
}

// return nodeId within the superlevel component  that has maximum scalar
vtkIdType MergeTree::ComponentMaximumQuery(vtkIdType& v,float& level){
  vector<vtkIdType> sortedIndices = argsort();
  vector<int> sortedIds(sortedIndices.size());
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  }
  
  int idx = sortedIds[v];
  int compMax = idx;
  float* scalarData = (float*) getScalar();
  if(scalarData[v] < level){
    // return bottom of superlevel  component
    for(auto id:sortedIndices){
      if(scalarData[id] >= level)
        return id;
    }
    return v;
  }else{
    queue<node*> nodes;
    nodes.push(mergeTree[idx]);
    while(nodes.size()){
      node* n = nodes.front();
      compMax = n->idx > compMax? n->idx: compMax;
      nodes.pop();
      if(n->parent && n->parent->idx >idx){
        nodes.push(n->parent);
      }
      for(auto child : n->children){
        if(child->idx > idx){
          nodes.push(n->parent);
        }
      }
    }
  }
  return sortedIndices[compMax];
}

void* MergeTree::getScalar(){
  vtkDataArray *scalarfield = usgrid->GetPointData()->GetArray(0);
  void *scalarData = scalarfield->GetVoidPointer(0);
  return scalarData;
}