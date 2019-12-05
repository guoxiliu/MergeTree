#ifndef UTILS_H
#define UTILS_H

#include <set>
#include <map>
#include <list>
#include <queue>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <float.h>
#include <stdio.h>
#include <vtkIdList.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <chrono>

using namespace std;

void* getScalar(vtkImageData*);
vtkIdType findSet(vector<vtkIdType> &, vtkIdType);
void unionSet(vector<vtkIdType> &, vtkIdType, vtkIdType);
vector<vtkIdType> argsort(const vector<vtkIdType>&, vtkImageData*, bool=true);
void bisect(int, int, vtkImageData *, vector<vector<vtkIdType>> &);
void decompose(int, vtkImageData *, vector<vector<vtkIdType>> &, set<pair<vtkIdType, vtkIdType>> &);
set<pair<vtkIdType, vtkIdType>> getLocalBridgeSet(const set<pair<vtkIdType, vtkIdType>> &, const vector<vtkIdType> &);
set<pair<vtkIdType, vtkIdType>> getReducedBridgeSet(const set<pair<vtkIdType, vtkIdType>> &, const vector<vtkIdType> &, vtkImageData *);

#endif
