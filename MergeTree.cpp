#include "MergeTree.h"

// Constructor.
MergeTree::MergeTree(vtkUnstructuredGrid *p){
  usgrid = p;
  //Set = vector<int>(p->GetNumberOfPoints());
  //graph = vector<vNode*>(p->GetNumberOfPoints()， new vNode());
}

// Build the merge tree.
int MergeTree::build(){
  vector<vtkIdType> indices(usgrid->GetNumberOfPoints());
  iota(indices.begin(), indices.end(), 0);
  vector<vtkIdType> sortedIndices = argsort(indices, usgrid);  
  constructJoin(sortedIndices);
  constructSplit(sortedIndices);
  mergeJoinSplit(joinTree, splitTree);
  return 0;
}
int MergeTree::build(vector<vtkIdType>& points){
  vector<vtkIdType> sortedIndices = argsort(points, usgrid);
  constructJoin(sortedIndices);
  constructSplit(sortedIndices);
  mergeJoinSplit(joinTree, splitTree);
  return 0;
}

// Construct the join tree.
void MergeTree::constructJoin(vector<vtkIdType>& sortedIndices){
  // record [vtkId] = pos in sortedindices
  map<vtkIdType, int> sortedIds;
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  }
  SetMax = vector<vtkIdType>(sortedIndices.size());
  //iota(Set.begin(), Set.end(), 0);
  iota(SetMax.begin(), SetMax.end(), 0);
  //iota(SetMin.begin(), SetMin.end(), 0);
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    node *ai = new node(i, sortedIndices[i]);
    joinTree.push_back(ai);
    graph.push_back(new vNode());
    graph[i]->jNode = ai;
    // get the neighbors of ai
    vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(sortedIndices[i], sortedIndices);
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
  map<vtkIdType, int> sortedIds;
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  }

  // iota(Set.begin(), Set.end(), 0);
  // iota(SetMax.begin(), SetMax.end(), 0);
  SetMin = vector<vtkIdType>(sortedIndices.size());
  iota(SetMin.begin(), SetMin.end(), 0);
  for (int i = sortedIndices.size()-1; i >= 0; --i){
    node *bi = new node(i, sortedIndices[i]);
    splitTree.push_back(bi);
    graph[i]->sNode = bi;
    // get the neighbors of bi
    vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(sortedIndices[i], sortedIndices);
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
  for(int i = 0; i < joinTree.size(); ++i){
    node *ci = new node(i, joinTree[i]->vtkIdx);
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

  printf("Merge tree is built!\n");
}

// Get the neighbors of the given vertex.
vtkSmartPointer<vtkIdList> MergeTree::getConnectedVertices(vtkIdType vertexId, const vector<vtkIdType> &vertexList){
  vtkSmartPointer<vtkIdList> connectedVertices = vtkSmartPointer<vtkIdList>::New();

  // get all cells that vertex 'id' is a part of 
  vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
  usgrid->GetPointCells(vertexId,cellIdList);


  set<vtkIdType> neighbors;
  for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds();++i){
    // get the points of each cell
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    usgrid->GetCellPoints(cellIdList->GetId(i),pointIdList);

    for(vtkIdType j = 0; j <pointIdList->GetNumberOfIds(); ++j){
      neighbors.insert(pointIdList->GetId(j));   
    }
  }
  neighbors.erase(vertexId);

  // add the neighbors only in the given vertex set
  set<vtkIdType> vertexSet(vertexList.begin(), vertexList.end());
  for(auto it = neighbors.begin(); it != neighbors.end(); ++it){
    if(vertexSet.find(*it) != vertexSet.end())
      connectedVertices->InsertNextId(*it);
  }

  return connectedVertices;
}

// Return all local maxima in the simplicial complex.
vector<vtkIdType> MergeTree::MaximaQuery(const set<pair<vtkIdType, vtkIdType>> &bridgeSet){
  // collect the higher end vertices from the bridge set
  set<vtkIdType> highEndVertices;
  for(auto it = bridgeSet.begin(); it != bridgeSet.end(); it++){
    highEndVertices.insert(it->second);   // the higher end vertex has the smaller scalar value
  }

  map<vtkIdType,int> sortedIds;
  for(unsigned int i = 0; i < mergeTree.size(); ++i){
    sortedIds[mergeTree[i]->vtkIdx] = mergeTree[i]->idx;
  }
  vector<vtkIdType> maxima;
  //iterate mergeTree to find local maximum
  for(auto node:mergeTree){
    if(node->numChildren == 0 && node->parent->idx < node->idx){
      if(highEndVertices.find(node->vtkIdx) == highEndVertices.end())
        maxima.push_back(node->vtkIdx);
    }
  }
  return maxima;
}

// Return the vertex within the superlevel component that has maximum scalar function value.
vtkIdType MergeTree::ComponentMaximumQuery(vtkIdType& v, float& level){
  // vector<vtkIdType> sortedIndices = argsort();
  map<vtkIdType,int> sortedIds;
  for(unsigned int i = 0; i < mergeTree.size(); ++i){
    //sortedIds[sortedIndices[i]] = i;
    sortedIds[mergeTree[i]->vtkIdx] = mergeTree[i]->idx;
  }
  
  int idx = sortedIds[v];
  int compMax = idx;
  float* scalarData = (float*) getScalar(usgrid);
  if(scalarData[v] < level){
    // return bottom of superlevel  component
    for(auto node:mergeTree){
      if(scalarData[node->vtkIdx] >= level)
        return node->vtkIdx;
    }
    return v;
  }else{
    queue<node*> nodes;
    nodes.push(mergeTree[idx]);
    set<int> visitedVertices;
    while(nodes.size()){
      node* n = nodes.front();
      compMax = n->idx > compMax? n->idx: compMax;
      nodes.pop();
      if(n->parent && n->parent->idx >idx && visitedVertices.find(n->parent->vtkIdx) == visitedVertices.end()){
        nodes.push(n->parent);
        visitedVertices.insert(n->parent->vtkIdx);
      }
      for(auto child : n->children){
        if(child->idx > idx && visitedVertices.find(child->vtkIdx) == visitedVertices.end()){
          nodes.push(n->parent);
          visitedVertices.insert(child->vtkIdx);
        }
      }
    }
  }
  return mergeTree[compMax]->vtkIdx;
}

