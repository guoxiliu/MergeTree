#include <vector>
#include <numeric>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

using namespace std;

// TODO: Super Arc class.
class SuperArc{
  public: 
    SuperArc(){ vertexIds = vector<vtkIdType>(); }
    inline vector<vtkIdType> getVertexIds(){
      return vertexIds;
    }
    inline void addVertex(vtkIdType id){
      vertexIds.push_back(id);
    }
  private:
    vector<vtkIdType> vertexIds; 
};


/**
 * Merge Tree Class.
 * The merge tree is created by combining join tree and split tree, which 
 * captures the evolution of the superlevel or sublevel sets.
 */ 
class MergeTree{
public:
  MergeTree(vtkUnstructuredGrid *p){
    usgrid = p;
    group = vector<int>(p->GetNumberOfPoints(), -1);
  }

  // Build the merge tree.
  int build(){
    vector<vtkIdType> sortedIndices = argsort();
    constructJoin();
    constructSplit();
    mergeJoinSplit();
    return 0;
  }

  // TODO: Implement the topological queries.

  // TODO: Maxima query

  // TODO: Component maximum query


protected:
  vector<int> group;    // Used for Union-Find algorithm
  vtkUnstructuredGrid *usgrid;    // Unstructed grid
  vector<SuperArc> arcs;    // Save the point set?

private:
  // Sort the scalar values while keeping track of the indices.
  vector<vtkIdType> argsort(){
    vector<vtkIdType> indices(usgrid->GetNumberOfPoints());
    iota(indices.begin(), indices.end(), 0);
    vtkDataArray *scalarfield = usgrid->GetPointData()->GetArray(0);
    switch(scalarfield->GetDataType()){
      case VTK_FLOAT:
      {
        float *scalarData = (float *)scalarfield->GetVoidPointer(0);
        sort(indices.begin(), indices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
        break;
      }
      case VTK_DOUBLE:
      {
        double *scalarData = (double *)scalarfield->GetVoidPointer(0);
        sort(indices.begin(), indices.end(), [scalarData](vtkIdType i1, vtkIdType i2) {return scalarData[i1] < scalarData[i2];});
        break;
      }
      default:
      {
        cout << "Type of scalarfield: " << scalarfield->GetDataType() << ", " << scalarfield->GetDataTypeAsString() << endl;
        break;
      }
    }

    return indices;
  }

  // Find the group id of a given vertex id.
  int findGroup(vtkIdType i){
    if(group[i] == -1)
      return i;
    group[i] = findGroup(group[i]);
    return group[i];
  }

  // Do union of two groups.
  void unionGroup(vtkIdType x, vtkIdType y){
    int xset = findGroup(x);
    int yset = findGroup(y);
    if(xset != yset){
      group[xset] = yset;
    }
  }

  // Construct the join tree.
  void constructJoin(){

  }

  // Construct the split tree.
  void constructSplit(){

  }

  // Merge the split and join tree.
  void mergeJoinSplit(){

  }
};
