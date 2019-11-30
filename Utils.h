#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <numeric>
#include <float.h>
#include <stdio.h>
#include <vtkPointData.h>
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
    int totalPoints = regions[i].size();

    // find the minimum and maximum coordinate value of current region.
    for(int j = 0; j < totalPoints; j++){
      double xyz[3];
      usgrid->GetPoint(regions[i][j], xyz);
      if(xyz[idx] < minVal) minVal = xyz[idx];
      if(xyz[idx] > maxVal) maxVal = xyz[idx];
    }

    // bisect the current region
    double center = (minVal + maxVal) / 2.0;
    vector<vtkIdType> first, second;
    for(int j = 0; j < totalPoints; j++){
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
void decompose(int numThreads, vtkUnstructuredGrid *usgrid, vector<vector<vtkIdType>> &regions, vector<vtkIdType> gBridgeSet){

  // initialize regions
  int totalVertices = usgrid->GetNumberOfPoints();
  regions = vector<vector<vtkIdType>>(1, vector<vtkIdType>(totalVertices));
  iota(regions[0].begin(), regions[0].end(), 0);

  // partition the vertices into different regions
  bisect(numThreads/2, 0, usgrid, regions);


  // test if the partition works well
  printf("Total regions: %lu\n", regions.size());
  for(unsigned int i = 0; i < regions.size(); i++){
    printf("Region %u: %lu\n", i, regions[i].size());
  }
}


#endif