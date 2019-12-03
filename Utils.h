#ifndef UTILS_H
#define UTILS_H

#include <set>
#include <map>
#include <queue>
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

int findSet(vector<vtkIdType> &, vtkIdType);
void unionSet(vector<vtkIdType> &, vtkIdType, vtkIdType);
void bisect(int, int, vtkUnstructuredGrid *, vector<vector<vtkIdType>> &);
void decompose(int, vtkUnstructuredGrid *, vector<vector<vtkIdType>> &, set<pair<vtkIdType, vtkIdType>> &);
set<pair<vtkIdType, vtkIdType>> getReducedBridgeSet(const set<pair<vtkIdType, vtkIdType>> &, vector<vtkIdType>, vtkUnstructuredGrid *);

#endif
