#include "MergeTree.h"

// Constructor.
MergeTree::MergeTree(vtkImageData *p){
  sgrid = p;
  //Set = vector<int>(p->GetNumberOfPoints());
  //graph = vector<vNode*>(p->GetNumberOfPoints()， new vNode());
  //getEdgeList(edgeList);
}

// Build the merge tree.
int MergeTree::build(){
  vector<vtkIdType> indices(sgrid->GetNumberOfPoints());
  iota(indices.begin(), indices.end(), 0);
  vector<vtkIdType> sortedIndices = argsort(indices, sgrid);  
  
  constructJoin(sortedIndices);
 
  constructSplit(sortedIndices);
  
  mergeJoinSplit(joinTree, splitTree);
  return 0;
}
int MergeTree::build(vector<vtkIdType>& points){
  vector<vtkIdType> sortedIndices = argsort(points, sgrid);
  constructJoin(sortedIndices);
  constructSplit(sortedIndices);
  mergeJoinSplit(joinTree, splitTree);
  return 0;
}

// Construct the join tree.
void MergeTree::constructJoin(vector<vtkIdType>& sortedIndices){
  // record [vtkId] = pos in sortedindices
  /* map<vtkIdType, int> sortedIds;
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  } */
  //printf("finished sortedIds!\n");
  //float* scalarData = (float*) getScalar(sgrid);
  //printf("Got the scalarData!\n");
  SetMax = vector<vtkIdType>(sortedIndices.size());
  //iota(Set.begin(), Set.end(), 0);
  iota(SetMax.begin(), SetMax.end(), 0);
  printf("finished iota setMax!\n");
  //iota(SetMin.begin(), SetMin.end(), 0);
  for(auto i = 0; i < sortedIndices.size(); ++i){
    node *ai = new node(i, sortedIndices[i]);
    joinTree.push_back(ai);
    //graph.push_back(new vNode());
    //graph[i]->jNode = ai;
    // get the neighbors of ai
    //vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(sortedIndices[i], sortedIndices);
    //printf("Got connected vertices!\n");
    set<vtkIdType> neighbors = edgeList[ai->vtkIdx];
    if(i %10000 == 0)
      printf("i = %d\n", i );
    //for(vtkIdType adj = 0; adj < connectedVertices->GetNumberOfIds();++adj){
      // value of adj 
    for(vtkIdType vtkIdx : neighbors){
      //vtkIdType vtkIdx = connectedVertices->GetId(adj);
      auto j = distance(sortedIndices.begin(), find(sortedIndices.begin(), sortedIndices.end(), vtkIdx));
      //printf("Got j!\n");
      if(j < i && findSet(SetMax, i) != findSet(SetMax, j)){
        auto k = findSet(SetMax, j);
        //graph[k]->jNode->parent = ai;
        joinTree[k]->parent = ai;
        //ai->children.push_back(graph[k]->jNode);
        ai->children.push_back(joinTree[k]);
        ai->numChildren += 1;
        unionSet(SetMax, i, j);
      }
    }
  }
  printf("Join tree is built!\n");
}

// Construct the split tree.
void MergeTree::constructSplit(vector<vtkIdType>& sortedIndices){
  // record [vtkId] = pos in sortedindices
  /* map<vtkIdType, int> sortedIds;
  for(unsigned int i = 0; i < sortedIndices.size(); ++i){
    sortedIds[sortedIndices[i]] = i;
  }
 */
  // iota(Set.begin(), Set.end(), 0);
  // iota(SetMax.begin(), SetMax.end(), 0);
  SetMin = vector<vtkIdType>(sortedIndices.size());
  iota(SetMin.begin(), SetMin.end(), 0);
  for (auto i = sortedIndices.size()-1; i >= 0; --i){
    node *bi = new node(i, sortedIndices[i]);
    splitTree.push_back(bi);
    //graph[i]->sNode = bi;
    // get the neighbors of bi
    //vtkSmartPointer<vtkIdList> connectedVertices = getConnectedVertices(sortedIndices[i], sortedIndices);
    set<vtkIdType> neighbors = edgeList[bi->vtkIdx];
    //for(vtkIdType adj = 0; adj < connectedVertices->GetNumberOfIds(); ++adj){
      // index of adj in sortedIndices
    for(vtkIdType vtkIdx : neighbors){
      //vtkIdType vtkIdx = connectedVertices->GetId(adj);
      auto j = distance(sortedIndices.begin(), find(sortedIndices.begin(), sortedIndices.end(), vtkIdx));
      
      if(j > i && findSet(SetMin, i) != findSet(SetMin, j)){
        auto k = findSet(SetMin, j);
        //graph[k]->sNode->parent = bi;
        splitTree[k]->parent = bi;
        //bi->children.push_back(graph[k]->sNode);
        bi->children.push_back(splitTree[k]);
        bi->numChildren += 1;
        unionSet(SetMin, i,j);
      }
    }
  }
  // reverse the order of split nodes from 0 to n-1
  reverse(splitTree.begin(), splitTree.end());
  printf("Split tree is built!\n");
}

// Merge the split and join tree.
void MergeTree::mergeJoinSplit(vector<node*>& joinTree, vector<node*>& splitTree){
  
  queue<long int> Q;
  for(auto i = 0; i < joinTree.size(); ++i){
    node *ci = new node(i, joinTree[i]->vtkIdx);
    mergeTree.push_back(ci);
    if(joinTree[i]->numChildren + splitTree[i]->numChildren == 1){
      Q.push(i);
    }
  }

  while (Q.size() > 1){
    auto i = Q.front();
    Q.pop();
    long int k;
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
  sgrid->GetPointCells(vertexId,cellIdList);


  set<vtkIdType> neighbors;
  for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds();++i){
    // get the points of each cell
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    sgrid->GetCellPoints(cellIdList->GetId(i),pointIdList);

    for(vtkIdType j = 0; j <pointIdList->GetNumberOfIds(); ++j){
      neighbors.insert(pointIdList->GetId(j));   
    }
  }
  neighbors.erase(vertexId);

  // add the neighbors only in the given vertex set
  //set<vtkIdType> vertexSet(vertexList.begin(), vertexList.end());
  for(auto it = neighbors.begin(); it != neighbors.end(); ++it){
    if(find(vertexList.begin(), vertexList.end(), *it) != vertexList.end())
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

  /* map<vtkIdType,int> sortedIds;
  for(long int i = 0; i < mergeTree.size(); ++i){
    sortedIds[mergeTree[i]->vtkIdx] = mergeTree[i]->idx;
  } */
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
  /* map<vtkIdType,int> sortedIds;
  
  for(int i = 0; i < treeSize; ++i){
    //sortedIds[sortedIndices[i]] = i;
    sortedIds[mergeTree[i]->vtkIdx] = mergeTree[i]->idx;
  } */
  //auto idx = sortedIds[v];

  long int idx;
  auto treeSize = mergeTree.size();
  for(auto i = 0; i < treeSize; ++i){
    if(mergeTree[i]->vtkIdx == v){
      idx = i;
    }
  }
  queue<node*> nodes;
  set<vtkIdType> visitedVertices;
  float* scalarData = (float*) getScalar(sgrid);

  auto compMax = scalarData[v] < level ? treeSize:idx;

  nodes.push(mergeTree[idx]);
  while(nodes.size()){
    node* n = nodes.front();
    nodes.pop();
    visitedVertices.insert(n->vtkIdx);

    if(scalarData[v] < level){
    // return bottom of superlevel  component
      if(scalarData[n->vtkIdx] >= level && n->idx < compMax)
        compMax = n->vtkIdx;  
    }else{
      compMax = n->idx > compMax? n->idx: compMax;
    }
    
    if(n->parent && n->parent->idx > idx && visitedVertices.find(n->parent->vtkIdx) == visitedVertices.end())
        nodes.push(n->parent);
    for(auto child : n->children){
      if(child->idx > idx && visitedVertices.find(child->vtkIdx) == visitedVertices.end())
        nodes.push(child);
    }
  }
  return compMax == treeSize ? v : mergeTree[compMax]->vtkIdx;
}

/* void MergeTree::getEdgeList(vector<set<vtkIdType>> &edgeList){
  //float *scalars = (float *)getScalar(sgrid);
  edgeList = vector<set<vtkIdType>>(sgrid->GetNumberOfPoints());
  int cellNum = sgrid->GetNumberOfCells();
  for(int i = 0; i < cellNum; i++){
    // get the points of each cell
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    sgrid->GetCellPoints(i, pointIdList);
    int cellSize = pointIdList->GetNumberOfIds();
    for(vtkIdType j = 0; j < cellSize; ++j){
      for(vtkIdType k = 0; k < cellSize ; ++k){
        if(k != j)
        edgeList[pointIdList->GetId(j)].insert(pointIdList->GetId(k));
      }
    }
  }
} */