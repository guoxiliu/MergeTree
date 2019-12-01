
#include "Utils.h"

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
void unionSet(vector<int> &group, vtkIdType i, vtkIdType j){
  
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
void bisect(int count, int idx, vtkUnstructuredGrid *usgrid, vector<vector<vtkIdType>> &regions){

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
      usgrid->GetPoint(regions[i][j], xyz);
      if(xyz[idx] < minVal) minVal = xyz[idx];
      if(xyz[idx] > maxVal) maxVal = xyz[idx];
    }

    // bisect the current region
    double center = (minVal + maxVal) / 2.0;
    vector<vtkIdType> first, second;
    for(vtkIdType j = 0; j < totalPoints; j++){
      double xyz[3];
      usgrid->GetPoint(regions[i][j], xyz);
      if(xyz[idx] < center) first.push_back(regions[i][j]);
      else second.push_back(regions[i][j]);
    }
    newRegions.push_back(first);
    newRegions.push_back(second);
  }

  // do bisection for the next dimension
  regions = newRegions;
  bisect(count/2, (idx+1)%3, usgrid, regions);
}

/**
 * Decompose the domain according to the number of threads.
 * Also create the global bridge set at the same time.
 */ 
void decompose(int numThreads, vtkUnstructuredGrid *usgrid, vector<vector<vtkIdType>> &regions, set<vtkIdType> &gBridgeSet){

  // initialize regions
  vtkIdType totalVertices = usgrid->GetNumberOfPoints();
  regions = vector<vector<vtkIdType>>(1, vector<vtkIdType>(totalVertices));
  iota(regions[0].begin(), regions[0].end(), 0);

  // partition the vertices into different regions
  bisect(numThreads/2, 0, usgrid, regions);


  // create the vertex-region map
  unordered_map<vtkIdType, int> vertexRegions;
  for(unsigned int i = 0; i < regions.size(); i++){
    for(unsigned int j = 0; j < regions[i].size(); j++){
      vertexRegions[regions[i][j]] = i;
    }
  }

  // loop all the cells in the grid
  gBridgeSet = set<vtkIdType>();
  vtkIdType totalCells = usgrid->GetNumberOfCells();
  for(vtkIdType i = 0; i < totalCells; i++){
    vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
    usgrid->GetCellPoints(i, cellPointIds);
    vtkIdType cellSize = cellPointIds->GetNumberOfIds();
    vector<int> regionIds(cellSize);
    for(vtkIdType j = 0; j < cellSize; j++){
      regionIds[j] = vertexRegions[cellPointIds->GetId(j)];
    }

    // see if all the vertices of the cell are in the same region
    bool inBridge = false;
    for(vtkIdType j = 1; j < cellSize; j++){
      if(regionIds[0] != regionIds[j])
        inBridge = true;
    }
    if(inBridge) gBridgeSet.insert(i);
  }
}

/**
 * Get the reduced bridge set.
 */ 
set<pair<vtkIdType, vtkIdType>> getReducedBrdigeSet(set<pair<vtkIdType, vtkIdType>> bridgeSet, vector<vtkIdType> vertexSet, vtkUnstructuredGrid *usgrid){
  // initialize
  set<pair<vtkIdType, vtkIdType>> reducedBS;
  int regionSize = vertexSet.size();
  vector<vtkIdType> component(vertexSet.begin(), vertexSet.end());

  // loop the vetex ids in decreasing order
  for(int i = regionSize-1; i >= 0; i--){
    // build the super link set of the current vertex id
    vector<vtkIdType> superLinks;
    // first get the cells incident to the vertex
    vtkSmartPointer<vtkIdList> pointCellIds = vtkSmartPointer<vtkIdList>::New();
    usgrid->GetPointCells(vertexSet[i], pointCellIds);
    vtkIdType pointCellSize = pointCellIds->GetNumberOfIds();
    for(vtkIdType j = 0; j < pointCellSize; j++){
      // for each cell get its vertices
      vtkSmartPointer<vtkIdList> cellPointIds = vtkSmartPointer<vtkIdList>::New();
      usgrid->GetCellPoints(pointCellIds->GetId(j), cellPointIds);
      vtkIdType cellSize = cellPointIds->GetNumberOfIds();
      for(vtkIdType k = 0; k < cellSize; k++){
        // see if the vertex id is greater than the current one
        if(pointCellIds->GetId(j) > vertexSet[i]){
          superLinks.push_back(pointCellIds->GetId(j));
        }
      }
    }
    
    // connect inside region
    for(vtkIdType &j:superLinks){
      pair<vtkIdType, vtkIdType> edge(vertexSet[i], j);
      if(bridgeSet.find(edge) == bridgeSet.end()){
        unionSet(component, vertexSet[i], j);
      }
    }

    // connect between regions
    for(vtkIdType &j:superLinks){
      if(findSet(component, vertexSet[i]) != findSet(component, j)){
        unionSet(component, vertexSet[i], j);
        reducedBS.insert(pair<vtkIdType, vtkIdType>(vertexSet[i], j));
      }
    }
  }
  return reducedBS;
}
