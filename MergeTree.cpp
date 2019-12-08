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
  auto start = chrono::high_resolution_clock::now();
  vector<vtkIdType> sortedIndices = argsort(vertexList, sgrid);
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << "Sort cost: " << duration.count() << " microseconds" <<endl;
  
  start = chrono::high_resolution_clock::now();
  constructJoin(sortedIndices);
  stop = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << "Join tree cost: " << duration.count() << " microseconds" <<endl;
  
  start = chrono::high_resolution_clock::now();
  constructSplit(sortedIndices);
  stop = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << "Split tree cost: " << duration.count() << " microseconds" <<endl;
  
  start = chrono::high_resolution_clock::now();
  mergeJoinSplit(joinTree, splitTree);
  stop = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << "Merge Join Split cost: " << duration.count() << " microseconds" <<endl;
  
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
  printf("Join tree built!\n");
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
    if(i%1000 == 0){
        for(auto id: upperLinks){
            printf("i = %d, vj = %d\n", i, (int)(id));
        }
    }
    for(vtkIdType &vj : upperLinks){
      // find the set of vi and vj
      // the scalar value of j should be greater
      vtkIdType iset = findSet(component, sortedIndices[i]);
      vtkIdType jset = findSet(component, vj);

      if(iset != jset){
        splitTree[jset]->parent = splitTree[iset];
        splitTree[iset]->children.push_back(splitTree[jset]);
        unionSet(component, iset, jset);
      }
    }
  }
  printf("Split tree built!\n");
}


// Merge the split and join tree.
void MergeTree::mergeJoinSplit(vector<node*>& joinTree, vector<node*>& splitTree){
  
  queue<int> Q;
  mergeTree = vector<node*>(joinTree.size(), nullptr);
  for(unsigned int i = 0; i < joinTree.size(); ++i){
    node *ci = new node(i);
    mergeTree.at(i) = ci;
    if(joinTree[i]->children.size() + splitTree[i]->children.size() == 1){
      Q.push(i);
    }
  }

  while (Q.size() > 1){
    int i = Q.front();
    Q.pop();
    vtkIdType k;
    if(joinTree[i]->children.size() == 0){
      k = joinTree[i]->parent->vtkIdx;
      mergeTree[i]->parent = mergeTree[k]; 
      mergeTree[k]->children.push_back(mergeTree[i]);
      // delete ai from join Tree
      auto it = find(joinTree[i]->parent->children.begin(), joinTree[i]->parent->children.end(),joinTree[i]);
      joinTree[i]->parent->children.erase(it);
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

        for(auto child : splitTree[i]->children){
          child->parent = splitTree[i]->parent;
          splitTree[i]->parent->children.push_back(child);
        }
      }
      splitTree[i] = nullptr;
    }else{
      k = splitTree[i]->parent->vtkIdx;
      mergeTree[i]->parent = mergeTree[k];
      mergeTree[k]->children.push_back(mergeTree[i]);
      //delete bi from split tree
      auto it = find(splitTree[i]->parent->children.begin(), splitTree[i]->parent->children.end(), splitTree[i]);
      splitTree[i]->parent->children.erase(it);
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
        
        for(auto child : joinTree[i]->children){
          child->parent = joinTree[i]->parent;
          joinTree[i]->parent->children.push_back(child);
        }
      }
      joinTree[i] = nullptr;
    }
    if(joinTree[k]->children.size() + splitTree[k]->children.size() == 1){
      Q.push(k);
    }
  }

  printf("Merge tree built!\n");
}


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
      vtkIdType vtkIdx = pointIdList->GetId(j);
      if((scalars[vtkIdx] < scalars[vertexId])|| (scalars[vtkIdx] == scalars[vertexId] && vtkIdx < vertexId))
        neighbors.insert(vtkIdx);
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
	printf("# of points in the cell is %d\n", (int) pointIdList->GetNumberOfIds());
    for(vtkIdType j = 0; j <pointIdList->GetNumberOfIds(); ++j){
      vtkIdType vtkIdx = pointIdList->GetId(j);
      if((scalars[vtkIdx] > scalars[vertexId]) || (scalars[vtkIdx] == scalars[vertexId] && vtkIdx > vertexId))
        neighbors.insert(vtkIdx);
    }
  }

  vector<vtkIdType> lowerLinks(neighbors.begin(), neighbors.end());

  return lowerLinks;
}


// Return all local maxima in the simplicial complex.
vector<vtkIdType> MergeTree::MaximaQuery(const set<pair<vtkIdType, vtkIdType>> &bridgeSet){
  // collect the higher end vertices from the bridge set
  set<vtkIdType> lowEndVertices;
  for(auto it = bridgeSet.begin(); it != bridgeSet.end(); it++){
    lowEndVertices.insert(it->first);   // the higher end vertex has the smaller scalar value
  }
  float* scalarData = (float*) getScalar(sgrid);

  vector<vtkIdType> maxima;
  //iterate mergeTree to find local maximum
  for(auto node:mergeTree){
    if(node->children.size() == 0 && scalarData[node->parent->vtkIdx] < scalarData[node->vtkIdx]){
      if(lowEndVertices.find(node->vtkIdx) == lowEndVertices.end())
        maxima.push_back(node->vtkIdx);
    }
  }
  return maxima;
}

// Return the vertex within the superlevel component that has maximum scalar function value.
vtkIdType MergeTree::ComponentMaximumQuery(vtkIdType& v, float& level){
  queue<node*> nodes;
  set<vtkIdType> visitedVertices;
  float* scalarData = (float*) getScalar(sgrid);
  if(scalarData[v] < level)
    return v;
  vtkIdType compMax = v;

  nodes.push(mergeTree[v]);
  while(nodes.size()){
    node* n = nodes.front();
    nodes.pop();
    visitedVertices.insert(n->vtkIdx);

    compMax = scalarData[n->vtkIdx] > scalarData[compMax]? n->vtkIdx: compMax;
    
    if(n->parent && scalarData[n->parent->vtkIdx] > level && visitedVertices.find(n->parent->vtkIdx) == visitedVertices.end())
        nodes.push(n->parent);
    for(auto child : n->children){
      if(scalarData[child->vtkIdx] > level && visitedVertices.find(child->vtkIdx) == visitedVertices.end())
        nodes.push(child);
    }
  }
  return compMax;
}

