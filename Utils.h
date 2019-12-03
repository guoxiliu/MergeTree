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

void* getScalar(vtkUnstructuredGrid*);
int findSet(vector<vtkIdType> &, vtkIdType);
void unionSet(vector<vtkIdType> &, vtkIdType, vtkIdType);
vector<vtkIdType> argsort(const vector<vtkIdType>&, vtkUnstructuredGrid*, bool=true);
void bisect(int, int, vtkUnstructuredGrid *, vector<vector<vtkIdType>> &);
void decompose(int, vtkUnstructuredGrid *, vector<vector<vtkIdType>> &, set<pair<vtkIdType, vtkIdType>> &);
set<pair<vtkIdType, vtkIdType>> getLocalBridgeSet(const set<pair<vtkIdType, vtkIdType>> &, const vector<vtkIdType> &);
set<pair<vtkIdType, vtkIdType>> getReducedBridgeSet(const set<pair<vtkIdType, vtkIdType>> &, const vector<vtkIdType> &, vtkUnstructuredGrid *);

#endif
