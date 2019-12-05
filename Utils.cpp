
#include "Utils.h"

/**
 * Get the void pointer of the scalar data.
 */ 
void* getScalar(vtkImageData* sgrid) {
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
  group[jset] = iset;
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

  // loop all the cells in the grid
  // gBridgeSet = set<pair<vtkIdType, vtkIdType>>();
  // vtkIdType totalCells = sgrid->GetNumberOfCells();
  // for(vtkIdType i = 0; i < totalCells; i++){
  //   vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
  //   sgrid->GetCellPoints(i, cellPointIds);
  //   vtkIdType cellSize = cellPointIds->GetNumberOfIds();
  //   vector<int> regionIds(cellSize);
  //   for(vtkIdType j = 0; j < cellSize; j++){
  //     regionIds[j] = vertexRegions[cellPointIds->GetId(j)];
  //   }

  //   // see if all the vertices of the cell are in the same region
  //   for(vtkIdType j = 0; j < cellSize-1; j++){
  //     for(vtkIdType k = j+1; k < cellSize; k++){
  //       if(regionIds[j] != regionIds[k]){
  //         pair<vtkIdType, vtkIdType> edge(cellPointIds->GetId(j), cellPointIds->GetId(k));
  //         if(scalars[edge.first] > scalars[edge.second]) 
  //           swap(edge.first, edge.second);
  //         gBridgeSet.insert(edge);
  //       }
  //     }
  //   }
  // }
}

/**
 * Get the local bridge set.
 */
set<pair<vtkIdType, vtkIdType>> getLocalBridgeSet(const set<pair<vtkIdType, vtkIdType>> &globalBridgeSet, const vector<vtkIdType> &vertexList){
  // create the local bridge set first
  set<vtkIdType> vertexSet(vertexList.begin(), vertexList.end());
  set<pair<vtkIdType, vtkIdType>> localBridgeSet;
  for(auto iter = globalBridgeSet.begin(); iter != globalBridgeSet.end(); iter++){
    if(vertexSet.find(iter->first) != vertexSet.end() || vertexSet.find(iter->second) != vertexSet.end()){
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
  set<pair<vtkIdType, vtkIdType>> reducedBS;
  int regionSize = vertexList.size();
  vector<vtkIdType> component(sgrid->GetNumberOfPoints(), -1);
  float *scalars = (float *)getScalar(sgrid);

  vector<vtkIdType> sortedVertices = argsort(vertexList, sgrid, false);
  // loop the vetex ids in decreasing order
  for(int i = 0; i < regionSize; i++){
   // build the upper link set of the current vertex id
   vector<vtkIdType> upperLinks;
   // first get the cells incident to the vertex
   vtkSmartPointer<vtkIdList> pointCellIds = vtkSmartPointer<vtkIdList>::New();
   sgrid->GetPointCells(sortedVertices[i], pointCellIds);
   vtkIdType pointCellSize = pointCellIds->GetNumberOfIds();
   for(vtkIdType j = 0; j < pointCellSize; j++){
     // for each cell get its vertices
     vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
     sgrid->GetCellPoints(pointCellIds->GetId(j), cellPointIds);
     vtkIdType cellSize = cellPointIds->GetNumberOfIds();
     for(vtkIdType k = 0; k < cellSize; k++){
       // see if the scalar value is greater than the current one
       if(scalars[cellPointIds->GetId(k)] > scalars[sortedVertices[i]]){
         upperLinks.push_back(cellPointIds->GetId(k));
       }
     }
   }
   
   // connect inside region
   for(vtkIdType &j:upperLinks){
     pair<vtkIdType, vtkIdType> edge(sortedVertices[i], j);
     if(bridgeSet.find(edge) == bridgeSet.end()){
       unionSet(component, sortedVertices[i], j);
     }
   }

   // connect between regions
   for(vtkIdType &j:upperLinks){
     if(findSet(component, sortedVertices[i]) != findSet(component, j)){
       unionSet(component, sortedVertices[i], j);
       reducedBS.insert(pair<vtkIdType, vtkIdType>(vertexList[i], j));
     }
   }
  }

  // printf("Size of the reduced bridge set: %zu\n", reducedBS.size());
  return reducedBS;
}
