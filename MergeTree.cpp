#include "MergeTree.h"

// Constructor.
MergeTree::MergeTree(vtkImageData *p){
  sgrid = p;
  //Set = vector<int>(p->GetNumberOfPoints());
  //graph = vector<vNode*>(p->GetNumberOfPoints()， new vNode());
}

// Build the merge tree.
int MergeTree::build(){
  vector<vtkIdType> indices(sgrid->GetNumberOfPoints());
  iota(indices.begin(), indices.end(), 0);
  return this->build(indices);
}
int MergeTree::build(vector<vtkIdType>& points){
  vector<vtkIdType> sortedIndices = argsort(points, sgrid);
  constructJoin(sortedIndices);
  // constructSplit(sortedIndices);
  // mergeJoinSplit(joinTree, splitTree);
  return 0;
}

// Construct the join tree.
void MergeTree::constructJoin(vector<vtkIdType>& sortedIndices){
  float *scalars = (float *)getScalar(sgrid);
  int regionSize = sortedIndices.size();
  vector<vtkIdType> component(regionSize, -1);
  printf("Component size: %d\n", regionSize);
  // if(lowerLinks.empty()){
  //   printf("OH YES!!!!!\n");
  // }else{
  //   fprintf(stderr, "FUCK! This is WRONG!\n");
  // }

  for(int i = 0; i < regionSize; ++i){
    if(i % 10000 == 0) printf("i: %d\n", i); 
    // node *ai = new node(sortedIndices[i]);
    // joinTree.push_back(ai);
    // get the neighbors of ai
    vector<vtkIdType> lowerLinks = getLowerLinks(sortedIndices[i]);
    if(lowerLinks.empty()){
      superArc *newarc = new superArc();
      component[sortedIndices[i]] = i;
      newarc->vertexList.push_back(sortedIndices[i]);
      vertexArcMap[sortedIndices[i]] = newarc;
      continue;
    }
    for(vtkIdType adj = 0; adj < lowerLinks.size(); ++adj){
      // find the set of i and j
      int iset = findSet(component, sortedIndices[i]);
      int jset = findSet(component, lowerLinks[adj]);
      
      if(jset == -1){
        if(iset != -1)
          fprintf(stderr, "THIS DOES NOT MAKE ANY SENSE!\n");
      // they are in the same component
      }else if(iset == -1){
        component[sortedIndices[i]] = jset;
        superArc *arc = vertexArcMap[lowerLinks[adj]];
        for(auto iter = arc->vertexList.begin(); iter != arc->vertexList.end(); iter++){
          if(scalars[sortedIndices[i]] > *iter){
            arc->vertexList.insert(iter, sortedIndices[i]);
            break;
          }
        }
        vertexArcMap[sortedIndices[i]] = arc;
      }else if(iset != jset){
        superArc *arci = vertexArcMap[sortedIndices[i]];
        superArc *arcj = vertexArcMap[lowerLinks[adj]];
        unionSet(component, iset, jset);
        superArc *newarc = new superArc();
        newarc->vertexList.push_back(sortedIndices[i]);
        arci->parent = newarc;
        arcj->parent = newarc;
        vertexArcMap[sortedIndices[i]] = newarc;
        newarc->children.push_back(arci);
        newarc->children.push_back(arcj);
      }
    }
  }
}

// Construct the split tree.
void MergeTree::constructSplit(vector<vtkIdType>& sortedIndices){
  // record [vtkId] = pos in sortedindices
  // map<vtkIdType, int> sortedIds;
  // for(unsigned int i = 0; i < sortedIndices.size(); ++i){
  //   sortedIds[sortedIndices[i]] = i;
  // }

  // set<vtkIdType> vertexSet(sortedIndices.begin(), sortedIndices.end());
  // SetMin = vector<vtkIdType>(sortedIndices.size());
  // iota(SetMin.begin(), SetMin.end(), 0);
  // for (int i = sortedIndices.size()-1; i >= 0; --i){
  //   node *bi = new node(i, sortedIndices[i]);
  //   splitTree.push_back(bi);
  //   graph[i]->sNode = bi;
  //   // get the neighbors of bi
  //   vector<vtkIdType> lowerLinks = getLowerLinks(sortedIndices[i]);
  //   for(vtkIdType adj = 0; adj < lowerLinks.size(); ++adj){
  //     if(vertexSet.find(lowerLinks[adj]) == vertexSet.end())
  //       continue;
  //     // index of adj in sortedIndices
  //     if(findSet(SetMin, sortedIndices[i]) != findSet(SetMin, lowerLinks[adj])){
  //       int k = findSet(SetMin, lowerLinks[adj]);
  //       graph[k]->sNode->parent = bi;
  //       bi->children.push_back(graph[k]->sNode);
  //       bi->numChildren += 1;
  //       unionSet(SetMin, i,j);
  //     }
  //   }
  // }
  // // reverse the order of split nodes from 0 to n-1
  // reverse(splitTree.begin(), splitTree.end());
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

// Get the neighbors of the given vertex.
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
