#include "MergeTree.h"

/**
 * Constructor.
 */ 
MergeTree::MergeTree(vtkImageData *p){
  sgrid = p;
  vertexList = vector<vtkIdType>(sgrid->GetNumberOfPoints());
  iota(vertexList.begin(), vertexList.end(), 0);
  sgrid->GetDimensions(dimension);
}

MergeTree::MergeTree(vtkImageData *p, vector<vtkIdType> idlist){
  sgrid = p;
  vertexList = idlist;
  sgrid->GetDimensions(dimension);
}

/**
 *  A wrapper function to build the merge tree.
 */ 
int MergeTree::build(){
  // auto start = chrono::high_resolution_clock::now();
  vector<size_t> sortedIndices = indexSort(vertexList, sgrid);
  // auto stop = chrono::high_resolution_clock::now();
  // auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  // printf("Index sort cost: %lld\n", duration.count());

  // start = chrono::high_resolution_clock::now();
  constructJoin(sortedIndices);
  // stop = chrono::high_resolution_clock::now();
  // duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  // printf("Join tree cost: %lld\n", duration.count());

  // start = chrono::high_resolution_clock::now();
  constructSplit(sortedIndices);
  // stop = chrono::high_resolution_clock::now();
  // duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  // printf("Split tree cost: %lld\n", duration.count());

  // start = chrono::high_resolution_clock::now();
  mergeJoinSplit();
  // stop = chrono::high_resolution_clock::now();
  // duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  // printf("Merge tree cost: %lld\n", duration.count());
  
  return 0;
}

/**
 * Construct the join tree.
 */ 
void MergeTree::constructJoin(vector<size_t>& sortedIndices){
  float *scalars = (float *)getScalar(sgrid);
  int regionSize = sortedIndices.size();
  vector<vtkIdType> component(regionSize, -1);

  joinTree = vector<node*>(regionSize, nullptr);
  for(int i = 0; i < regionSize; i++){
    size_t idx = sortedIndices[i];
    node *newnode = new node(idx);
    joinTree[idx] = newnode;

    vector<vtkIdType> neighbors = getConnectedVertices(vertexList[idx], sgrid, dimension);
    for(vtkIdType &vj : neighbors){
      // see if the vertex is in the range
      if(vj < vertexList.front() || vj > vertexList.back()) 
        continue;
      if((scalars[vj] < scalars[vertexList[idx]]) || (scalars[vj] == scalars[vertexList[idx]] && vj < vertexList[idx])){
        // find the set of vi and vj
        // the scalar value of j should be lower
        vtkIdType iset = findSet(component, idx);
        vtkIdType jset = findSet(component, vj-vertexList.front());

        if(iset != jset){
          joinTree[jset]->parent = joinTree[iset];
          joinTree[iset]->children.push_back(joinTree[jset]);
          unionSet(component, iset, jset);
        }
      }
    }
  }
  // printf("Join tree built!\n");
}

/**
 * Construct the split tree.
 */ 
void MergeTree::constructSplit(vector<size_t>& sortedIndices){
  float *scalars = (float *)getScalar(sgrid);
  int regionSize = sortedIndices.size();
  vector<vtkIdType> component(regionSize, -1);

  splitTree = vector<node*>(regionSize, nullptr);
  for(int i = regionSize-1; i >= 0; i--){
    size_t idx = sortedIndices[i];
    node *newnode = new node(idx);
    splitTree[idx] = newnode;

    vector<vtkIdType> neighbors = getConnectedVertices(vertexList[idx], sgrid, dimension);
    for(vtkIdType &vj : neighbors){
      // find the set of vi and vj
      // the scalar value of j should be greater
      if (vj < vertexList.front() || vj > vertexList.back())
        continue;
      if((scalars[vj] > scalars[vertexList[idx]]) || (scalars[vj] == scalars[vertexList[idx]] && vj > vertexList[idx])){
          vtkIdType iset = findSet(component, idx);
          vtkIdType jset = findSet(component, vj-vertexList.front());

          if(iset != jset){
            splitTree[jset]->parent = splitTree[iset];
            splitTree[iset]->children.push_back(splitTree[jset]);
            unionSet(component, iset, jset);
          }
      }
    }
  }
  // printf("Split tree created!\n");
}


/**
 * Merge the split and join tree.
 */ 
void MergeTree::mergeJoinSplit(){
  
  queue<int> leavesQueue;
  mergeTree = vector<node*>(joinTree.size(), nullptr);

  // construct a queue of leaves
  for(unsigned int i = 0; i < joinTree.size(); ++i){
    node *ci = new node(i);
    mergeTree.at(i) = ci;
    if(joinTree[i]->children.size() + splitTree[i]->children.size() == 1){
      leavesQueue.push(i);
    }
  }

  while (!leavesQueue.empty()){
    int i = leavesQueue.front();
    leavesQueue.pop();
    vtkIdType k;
    // if ai is the lower leaf, i.e., from the join tree
    if(joinTree[i]->children.size() == 0){
      // add (ai, bi) to the merge tree
      k = joinTree[i]->parent->vtkIdx;
      mergeTree[i]->parent = mergeTree[k]; 
      mergeTree[k]->children.push_back(mergeTree[i]);

      // delete ai from join tree
      auto it = find(joinTree[i]->parent->children.begin(), joinTree[i]->parent->children.end(),joinTree[i]);
      joinTree[i]->parent->children.erase(it);
      joinTree[i] = nullptr;

      // delete ai from split tree
      // connect bi's parent with bi's children
      // if ai is the root of split tree, just make its children's parents to null
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

    // if vi is the upper leaf, i.e., from the split tree
    }else{
      // add (ai, bi) to the merge tree
      k = splitTree[i]->parent->vtkIdx;
      mergeTree[i]->parent = mergeTree[k];
      mergeTree[k]->children.push_back(mergeTree[i]);

      //delete ai from split tree
      auto it = find(splitTree[i]->parent->children.begin(), splitTree[i]->parent->children.end(), splitTree[i]);
      splitTree[i]->parent->children.erase(it);
      splitTree[i] = nullptr;

      // delete ai from the join tree
      // connect bi's parent with bi's children
      // if bi is the root of join tree, just make its children's parents to null
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
    // if bi is a leaf, then enqueue
    if(joinTree[k]->children.size() + splitTree[k]->children.size() == 1){
      leavesQueue.push(k);
    }
  }

  // printf("Merge tree built!\n");
}


/**
 * Return all local maxima in the simplicial complex.
 */ 
vector<vtkIdType> MergeTree::MaximaQuery(const set<pair<vtkIdType, vtkIdType>> &bridgeSet){
  // collect the higher end vertices from the bridge set
  set<vtkIdType> lowEndVertices;
  float* scalarData = (float*) getScalar(sgrid);
  for(auto it = bridgeSet.begin(); it != bridgeSet.end(); it++){
    lowEndVertices.insert(it->first);   // the lower end vertex has the smaller scalar value
    if(scalarData[it->first] == scalarData[it->second])
      lowEndVertices.insert(it->second);
  }

  vector<vtkIdType> maxima;
  //iterate mergeTree to find local maximum
  for(auto node:mergeTree){
    if(node->children.size() == 0 && scalarData[vertexList[node->parent->vtkIdx]] < scalarData[vertexList[node->vtkIdx]]){
      if(lowEndVertices.find(vertexList[node->vtkIdx]) == lowEndVertices.end())
        maxima.push_back(vertexList[node->vtkIdx]);
    }
  }
  return maxima;
}

/**
 * Return the vertex within the superlevel component that has maximum scalar function value.
 */ 
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
