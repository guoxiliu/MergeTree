#ifndef UTILS_H
#define UTILS_H

#include <set>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <float.h>
#include <stdio.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

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
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    usgrid->GetCellPoints(i, cellIds);
    vtkIdType cellSize = cellIds->GetNumberOfIds();
    int regionIds[cellSize];
    for(vtkIdType j = 0; j < cellSize; j++){
      regionIds[j] = vertexRegions[cellIds->GetId(j)];
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




#endif