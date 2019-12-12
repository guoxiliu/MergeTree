
#include "Utils.h"

/**
 * Get the void pointer of the scalar data.
 */ 
void* getScalar(vtkImageData* sgrid){
  vtkDataArray *scalarfield = sgrid->GetPointData()->GetArray(0);
  void *scalarData = scalarfield->GetVoidPointer(0);
  return scalarData;
}

/**
 * Sort the scalar values while keeping track of the indices.
 */  
vector<size_t> indexSort(const vector<vtkIdType>& vertexList, vtkImageData* sgrid, bool increasing){
  float *scalarData = (float*)getScalar(sgrid);
  vector<size_t> idx(vertexList.size());
  iota(idx.begin(), idx.end(), 0);

  if(increasing){
    stable_sort(idx.begin(), idx.end(), [scalarData, vertexList](vtkIdType i1, vtkIdType i2) {return scalarData[vertexList[i1]] < scalarData[vertexList[i2]];});
  }else{
    stable_sort(idx.begin(), idx.end(), [scalarData, vertexList](vtkIdType i1, vtkIdType i2) {return scalarData[vertexList[i1]] > scalarData[vertexList[i2]];});
  }
  return idx;
}

/**
 * Sort the scalar values while keeping track of the indices.
 */ 
vector<vtkIdType> argsort(const vector<vtkIdType>& vertexList, vtkImageData* sgrid, bool increasing){
  float *scalarData = (float*)getScalar(sgrid);
  vector<vtkIdType> sortedVertices(vertexList.begin(), vertexList.end());
  if(increasing){
    stable_sort(sortedVertices.begin(), sortedVertices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
  }else{
    stable_sort(sortedVertices.begin(), sortedVertices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] > scalarData[i2];});
  }
  return sortedVertices;
}

/**
 * Find the set id of a given vertex id.
 */
vtkIdType findSet(vector<vtkIdType> &group, vtkIdType i){
  if(group[i] == -1)
    return i;
  group[i] = findSet(group, group[i]);
  return group[i];
}

/**
 * Do union of two sets.
 */ 
void unionSet(vector<vtkIdType> &group, vtkIdType i, vtkIdType j){
  vtkIdType iset = findSet(group, i);
  vtkIdType jset = findSet(group, j);
  if(iset != jset)
    group[jset] = iset;
}

/**
 * Get connected vertices with a given vertex id.
 */ 
vector<vtkIdType> getConnectedVertices(vtkIdType id, const vtkImageData *sgrid, int dim[3]){
  vector<vtkIdType> connectedVertices;
  int ids[3];
  ids[0] = id % dim[0];
  ids[2] = id / (dim[0]*dim[1]);
  ids[1] = (id % (dim[0]*dim[1]) / dim[0]);
  for(int i = 0; i < 3; i++){
    ids[i] -= 1;
    if(ids[i] >= 0 && ids[i] < dim[i]){
      connectedVertices.push_back((ids[0]+ids[1]*dim[0]+ids[2]*dim[0]*dim[1]));
    }
    ids[i] += 2;
    if(ids[i] >= 0 && ids[i] < dim[i]){
      connectedVertices.push_back((ids[0]+ids[1]*dim[0]+ids[2]*dim[0]*dim[1]));
    }
    ids[i] -= 1;
  }

  return connectedVertices;
}


/**
 * Decompose the domain according to the number of threads.
 * Also create the global bridge set at the same time.
 */ 
void decompose(int numThreads, vtkImageData *sgrid, vector<vector<vtkIdType>> &regions, set<pair<vtkIdType, vtkIdType>> &gBridgeSet){

  // initialize regions
  vtkIdType totalVertices = sgrid->GetNumberOfPoints();
  vtkIdType regionPoints = totalVertices / numThreads;
  regions = vector<vector<vtkIdType>>(numThreads);

  vtkIdType startId = 0;
  for(int i = 0; i < numThreads-1; i++){
    regions[i] = vector<vtkIdType>(regionPoints);
    iota(regions[i].begin(), regions[i].end(), startId);
    startId += regionPoints;
  }

  regions[numThreads-1] = vector<vtkIdType>(regionPoints + totalVertices%numThreads);
  iota(regions[numThreads-1].begin(), regions[numThreads-1].end(), startId);

  // Create the global bridge set
  float *scalars = (float *)getScalar(sgrid);
  gBridgeSet = set<pair<vtkIdType, vtkIdType>>();
  vtkIdType totalCells = sgrid->GetNumberOfCells();

  // loop all the cells in the grid
  for(vtkIdType i = 0; i < totalCells; i++){
  vtkCell *cell = sgrid->GetCell(i);
    for(int e = 0; e < cell->GetNumberOfEdges(); e++){
      vtkCell *cellEdge = cell->GetEdge(e);
      vtkIdList *pointIdList = cellEdge->GetPointIds();
      if(pointIdList->GetId(0)/regionPoints != pointIdList->GetId(1)/regionPoints){
        pair<vtkIdType, vtkIdType> edge(pointIdList->GetId(0), pointIdList->GetId(1));
        if(scalars[edge.first] > scalars[edge.second]){
          swap(edge.first, edge.second);
        }else if(scalars[edge.first] == scalars[edge.second] && edge.first > edge.second){
          swap(edge.first, edge.second);
        }
        gBridgeSet.insert(edge);
      }
    }
  }
}

/**
 * Get the local bridge set.
 */
set<pair<vtkIdType, vtkIdType>> getLocalBridgeSet(const set<pair<vtkIdType, vtkIdType>> &globalBridgeSet, const vector<vtkIdType> &vertexList){
  set<pair<vtkIdType, vtkIdType>> localBridgeSet;
  for(auto iter = globalBridgeSet.begin(); iter != globalBridgeSet.end(); iter++){
    if(iter->first >= vertexList.front() && iter->first <= vertexList.back()){
    // if the first vertex is in the region
      localBridgeSet.insert(*iter);
    }else if(iter->second >= vertexList.front() && iter->second <= vertexList.back()){
    // if the second vertex is in the region
      localBridgeSet.insert(*iter);
    }
  }
  return localBridgeSet;
}

/**
 * Get the reduced bridge set.
 */ 
set<pair<vtkIdType, vtkIdType>> getReducedBridgeSet(const set<pair<vtkIdType, vtkIdType>> &bridgeSet, const vector<vtkIdType> &vertexList, vtkImageData *sgrid){
  // initialize
  int dimension[3];
  sgrid->GetDimensions(dimension);
  int regionSize = vertexList.size();
  set<pair<vtkIdType, vtkIdType>> reducedBS;
  float *scalars = (float *)getScalar(sgrid);

  vector<vtkIdType> component(regionSize, -1);
  vector<vtkIdType> sortedVertices = argsort(vertexList, sgrid, false);

  // loop the vertex ids in decreasing order
  for(int i = 0; i < regionSize; i++){
    // build the upper link set of the current vertex id
    vector<vtkIdType> neighbors = getConnectedVertices(sortedVertices[i], sgrid, dimension);
    for(vtkIdType &vj : neighbors){
      // upper links
      if((scalars[vj] >= scalars[sortedVertices[i]]) || (scalars[vj] == scalars[sortedVertices[i]] && vj > sortedVertices[i])){
        // connect inside region
        pair<vtkIdType, vtkIdType> edge(sortedVertices[i], vj);
        if(bridgeSet.find(edge) == bridgeSet.end()){
          unionSet(component, sortedVertices[i], vj);
        }
      }
    }
    // connect between regions
    for(vtkIdType &vj : neighbors){
      if((scalars[vj] >= scalars[sortedVertices[i]]) || (scalars[vj] == scalars[sortedVertices[i]] && vj > sortedVertices[i])){
        if(findSet(component, sortedVertices[i]) != findSet(component, vj)){
          unionSet(component, sortedVertices[i], vj);
          reducedBS.insert(pair<vtkIdType, vtkIdType>(vertexList[i], vj));
        }
      }
    }
  }

  return reducedBS;
}
