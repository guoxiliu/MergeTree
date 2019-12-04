#include "MergeTree.h"

// Constructor.
MergeTree::MergeTree(vtkImageData *p){
  sgrid = p;
  vertexList = vector<vtkIdType>(sgrid->GetNumberOfPoints());
  iota(vertexList.begin(), vertexList.end(), 0);
}

MergeTree::MergeTree(vtkImageData *p, vector<vtkIdType> idlist){
  sgrid = p;
  vertexList = idlist;
}

// Build the merge tree.
int MergeTree::build(){
  vector<vtkIdType> sortedIndices = argsort(vertexList, sgrid);
  constructJoin(sortedIndices);
  constructSplit(sortedIndices);
  // mergeJoinSplit(joinTree, splitTree);
  return 0;
}

// Construct the join tree.
void MergeTree::constructJoin(vector<vtkIdType>& sortedIndices){
  int regionSize = sortedIndices.size();
  vector<vtkIdType> component(regionSize, -1);

  joinTree = vector<node*>(regionSize, nullptr);
  for(int i = 0; i < regionSize; i++){
    node *newnode = new node(sortedIndices[i]);
    joinTree.at(sortedIndices[i]) = newnode;

    vector<vtkIdType> lowerLinks = getLowerLinks(sortedIndices[i]);
    for(vtkIdType &vj : lowerLinks){
      // find the set of vi and vj
      // the scalar value of j should be lower
      vtkIdType iset = findSet(component, sortedIndices[i]);
      vtkIdType jset = findSet(component, vj);

      if(iset != jset){
        joinTree[jset]->parent = joinTree[iset];
        joinTree[iset]->children.push_back(joinTree[jset]);
        unionSet(component, iset, jset);
      }
    }
  }
}

// Construct the split tree.
void MergeTree::constructSplit(vector<vtkIdType>& sortedIndices){
  int regionSize = sortedIndices.size();
  vector<vtkIdType> component(regionSize, -1);

  splitTree = vector<node*>(regionSize, nullptr);
  for(int i = regionSize-1; i >= 0; i--){
    node *newnode = new node(sortedIndices[i]);
    splitTree.at(sortedIndices[i]) = newnode;

    vector<vtkIdType> upperLinks = getUpperLinks(sortedIndices[i]);
    for(vtkIdType &vj : upperLinks){
      // find the set of vi and vj
      // the scalar value of j should be greater
      vtkIdType iset = findSet(component, sortedIndices[i]);
      vtkIdType jset = findSet(component, vj);

      if(iset != jset){
        splitTree[jset]->parent = splitTree[iset];
        splitTree[iset]->children.push_back(splitTree[jset]);
        unionSet(component, jset, iset);
      }
    }
  }
}

/*
// Merge the split and join tree.
void MergeTree::mergeJoinSplit(vector<node*>& joinTree, vector<node*>& splitTree){
  
  queue<int> Q;
  for(unsigned int i = 0; i < joinTree.size(); ++i){
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
*/

// Get the lower links of the given vertex.
vector<vtkIdType> MergeTree::getLowerLinks(vtkIdType vertexId){
  float *scalars = (float*)getScalar(sgrid);

  // get all cells that vertex 'id' is a part of 
  vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
  sgrid->GetPointCells(vertexId,cellIdList);

  set<vtkIdType> neighbors;
  for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds();++i){
    // get the points of each cell
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    sgrid->GetCellPoints(cellIdList->GetId(i), pointIdList);

    for(vtkIdType j = 0; j <pointIdList->GetNumberOfIds(); ++j){
      if(scalars[pointIdList->GetId(j)] < scalars[vertexId])
        neighbors.insert(pointIdList->GetId(j));
    }
  }

  vector<vtkIdType> lowerLinks(neighbors.begin(), neighbors.end());

  return lowerLinks;
}

// Get the upper links of the given vertex.
vector<vtkIdType> MergeTree::getUpperLinks(vtkIdType vertexId){
  float *scalars = (float*)getScalar(sgrid);

  // get all cells that vertex 'id' is a part of 
  vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
  sgrid->GetPointCells(vertexId,cellIdList);

  set<vtkIdType> neighbors;
  for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds();++i){
    // get the points of each cell
    vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
    sgrid->GetCellPoints(cellIdList->GetId(i), pointIdList);

    for(vtkIdType j = 0; j <pointIdList->GetNumberOfIds(); ++j){
      if(scalars[pointIdList->GetId(j)] > scalars[vertexId])
        neighbors.insert(pointIdList->GetId(j));
    }
  }

  vector<vtkIdType> lowerLinks(neighbors.begin(), neighbors.end());

  return lowerLinks;
}

/*
// Return all local maxima in the simplicial complex.
vector<vtkIdType> MergeTree::MaximaQuery(const set<pair<vtkIdType, vtkIdType>> &bridgeSet){
  // collect the higher end vertices from the bridge set
  set<vtkIdType> lowEndVertices;
  for(auto it = bridgeSet.begin(); it != bridgeSet.end(); it++){
    lowEndVertices.insert(it->first);   // the higher end vertex has the smaller scalar value
  }

  map<vtkIdType,int> sortedIds;
  for(unsigned int i = 0; i < mergeTree.size(); ++i){
    sortedIds[mergeTree[i]->vtkIdx] = mergeTree[i]->idx;
  }
  vector<vtkIdType> maxima;
  //iterate mergeTree to find local maximum
  for(auto node:mergeTree){
    if(node->numChildren == 0 && node->parent->idx < node->idx){
      if(lowEndVertices.find(node->vtkIdx) == lowEndVertices.end())
        maxima.push_back(node->vtkIdx);
    }
  }
  return maxima;
}

// Return the vertex within the superlevel component that has maximum scalar function value.
vtkIdType MergeTree::ComponentMaximumQuery(vtkIdType& v, float& level){
  // vector<vtkIdType> sortedIndices = argsort();
  map<vtkIdType,int> sortedIds;
  int treeSize = mergeTree.size();
  for(int i = 0; i < treeSize; ++i){
    //sortedIds[sortedIndices[i]] = i;
    sortedIds[mergeTree[i]->vtkIdx] = mergeTree[i]->idx;
  }
  int idx = sortedIds[v];
  queue<node*> nodes;
  set<vtkIdType> visitedVertices;
  float* scalarData = (float*) getScalar(sgrid);

  int compMax = scalarData[v] < level ? treeSize:idx;

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
*/
