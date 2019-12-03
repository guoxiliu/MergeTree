
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
vector<vtkIdType> argsort(const vector<vtkIdType>& vertexSet, vtkImageData* sgrid, bool increasing){
  vector<vtkIdType> sortedVertices(vertexSet.begin(), vertexSet.end());
  vtkDataArray *scalarfield = sgrid->GetPointData()->GetArray(0);
  switch(scalarfield->GetDataType()){
    case VTK_FLOAT:
    {
      float *scalarData = (float *)scalarfield->GetVoidPointer(0);
      //float *scalarData = (float*) getScalar();
      if(increasing){
        sort(sortedVertices.begin(), sortedVertices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
      }else{
        sort(sortedVertices.begin(), sortedVertices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] > scalarData[i2];});
      }
      break;
    }
    case VTK_DOUBLE:
    {
      double *scalarData = (double *)scalarfield->GetVoidPointer(0);
      //double* scalarData = (double*) getScalar();
      if(increasing){
        sort(sortedVertices.begin(), sortedVertices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
      }else{
        sort(sortedVertices.begin(), sortedVertices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] > scalarData[i2];});
      }
      break;
    }
    default:
    {
      cout << "Type of scalarfield: " << scalarfield->GetDataType() << ", " << scalarfield->GetDataTypeAsString() << endl;
      break;
    }
  }
  return sortedVertices;
}

/**
 * Find the set id of a given vertex id.
 */
int findSet(vector<vtkIdType> &group, vtkIdType i){
  if(group[i] == i)
    return i;
  group[i] = findSet(group, group[i]);
  return group[i];
}

/**
 * Do union of two sets.
 */ 
void unionSet(vector<vtkIdType> &group, vtkIdType i, vtkIdType j){
  
  int iset = findSet(group, i);
  int jset = findSet(group, j);
  
  if(i > j){ // union for SetMax
    if(iset < jset){
      group[iset] = jset;
    }else{
      group[jset] = iset;
    }
  }else{ // union for SetMin
    if(iset < jset){
      group[jset] = iset;
    }else{
      group[iset] = jset;
    }
  } 
}

/**
 * Recursive bisection.
 */
void bisect(int count, int idx, vtkImageData *sgrid, vector<vector<vtkIdType>> &regions){

  if(count <= 0) return;

  vector<vector<vtkIdType>> newRegions;

  // loop through every region
  for(unsigned int i = 0; i < regions.size(); i++){
    // find the minimum and maximum idx coordinates
    double minVal = DBL_MAX, maxVal = DBL_MIN;
    vtkIdType totalPoints = regions[i].size();

    // find the minimum and maximum coordinate value of current region.
    for(vtkIdType j = 0; j < totalPoints; j++){
      double xyz[3];
      sgrid->GetPoint(regions[i][j], xyz);
      if(xyz[idx] < minVal) minVal = xyz[idx];
      if(xyz[idx] > maxVal) maxVal = xyz[idx];
    }

    // bisect the current region
    double center = (minVal + maxVal) / 2.0;
    vector<vtkIdType> first, second;
    for(vtkIdType j = 0; j < totalPoints; j++){
      double xyz[3];
      sgrid->GetPoint(regions[i][j], xyz);
      if(xyz[idx] < center) first.push_back(regions[i][j]);
      else second.push_back(regions[i][j]);
    }
    newRegions.push_back(first);
    newRegions.push_back(second);
  }

  // do bisection for the next dimension
  regions = newRegions;
  bisect(count/2, (idx+1)%3, sgrid, regions);
}

/**
 * Decompose the domain according to the number of threads.
 * Also create the global bridge set at the same time.
 */ 
void decompose(int numThreads, vtkImageData *sgrid, vector<vector<vtkIdType>> &regions, set<pair<vtkIdType, vtkIdType>> &gBridgeSet){

  // initialize regions
  vtkIdType totalVertices = sgrid->GetNumberOfPoints();
  regions = vector<vector<vtkIdType>>(1, vector<vtkIdType>(totalVertices));
  iota(regions[0].begin(), regions[0].end(), 0);

  // partition the vertices into different regions
  bisect(numThreads/2, 0, sgrid, regions);


  // create the vertex-region map
  unordered_map<vtkIdType, int> vertexRegions;
  for(unsigned int i = 0; i < regions.size(); i++){
    for(unsigned int j = 0; j < regions[i].size(); j++){
      vertexRegions[regions[i][j]] = i;
    }
  }

  float *scalars = (float *)getScalar(sgrid);

  // loop all the cells in the grid
  gBridgeSet = set<pair<vtkIdType, vtkIdType>>();
  vtkIdType totalCells = sgrid->GetNumberOfCells();
  for(vtkIdType i = 0; i < totalCells; i++){
    vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
    sgrid->GetCellPoints(i, cellPointIds);
    vtkIdType cellSize = cellPointIds->GetNumberOfIds();
    vector<int> regionIds(cellSize);
    for(vtkIdType j = 0; j < cellSize; j++){
      regionIds[j] = vertexRegions[cellPointIds->GetId(j)];
    }

    // see if all the vertices of the cell are in the same region
    for(vtkIdType j = 0; j < cellSize-1; j++){
      for(vtkIdType k = j+1; k < cellSize; k++){
        if(regionIds[j] != regionIds[k]){
          pair<vtkIdType, vtkIdType> edge(cellPointIds->GetId(j), cellPointIds->GetId(k));
          if(scalars[edge.first] < scalars[edge.second]) 
            swap(edge.first, edge.second);
          gBridgeSet.insert(edge);
        }
      }
    }
  }
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
  vector<vtkIdType> component(sgrid->GetNumberOfPoints());
  iota(component.begin(), component.end(), 0);
  float *scalars = (float *)getScalar(sgrid);

  vector<vtkIdType> sortedVertices = argsort(vertexList, sgrid, false);
  // loop the vetex ids in decreasing order
  for(int i = 0; i < regionSize; i++){
   // build the super link set of the current vertex id
   vector<vtkIdType> superLinks;
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
       if(scalars[cellPointIds->GetId(k)] < scalars[sortedVertices[i]]){
         superLinks.push_back(cellPointIds->GetId(k));
       }
     }
   }
   
   // connect inside region
   for(vtkIdType &j:superLinks){
     pair<vtkIdType, vtkIdType> edge(sortedVertices[i], j);
     if(bridgeSet.find(edge) == bridgeSet.end()){
       unionSet(component, sortedVertices[i], j);
     }
   }

   // connect between regions
   for(vtkIdType &j:superLinks){
     if(findSet(component, sortedVertices[i]) != findSet(component, j)){
       unionSet(component, sortedVertices[i], j);
       reducedBS.insert(pair<vtkIdType, vtkIdType>(vertexList[i], j));
     }
   }
  }

  // printf("Size of the reduced bridge set: %zu\n", reducedBS.size());
  return reducedBS;
}
