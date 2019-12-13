#ifndef UTILS_H
#define UTILS_H

#include <set>
#include <map>
#include <list>
#include <queue>
#include <vector>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <omp.h>
#include <float.h>
#include <stdio.h>
#include <vtkCell.h>
#include <vtkIdList.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <chrono>

using namespace std;

void* getScalar(vtkImageData *);
vtkIdType findSet(vector<vtkIdType> &, vtkIdType);
void unionSet(vector<vtkIdType> &, vtkIdType, vtkIdType);

vector<size_t> indexSort(const vector<vtkIdType> &, vtkImageData *, bool=true);
vector<vtkIdType> argsort(const vector<vtkIdType> &, vtkImageData *, bool=true);
vector<vtkIdType> getConnectedVertices(vtkIdType, const vtkImageData *, int[3]);
void decompose(int, vtkImageData *, vector<vector<vtkIdType>> &, set<pair<vtkIdType, vtkIdType>> &);
set<pair<vtkIdType, vtkIdType>> getLocalBridgeSet(const set<pair<vtkIdType, vtkIdType>> &, const vector<vtkIdType> &);
set<pair<vtkIdType, vtkIdType>> getReducedBridgeSet(const set<pair<vtkIdType, vtkIdType>> &, const vector<vtkIdType> &, vtkImageData *);

#endif
