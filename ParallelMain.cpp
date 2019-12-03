#include "MergeTree.h"
#include <omp.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

using namespace std;

const int threadNum = 8;

int main ( int argc, char *argv[] )
{
  // parse command line arguments
  if(argc < 2){
    fprintf(stderr, "Usage: %s Filename(.vtu)\n", argv[0]);
    return 1;
  }

  string filename = argv[1];
  string extension = filename.substr(filename.length() - 3);

  if(extension != "vtu"){
    fprintf(stderr, "The file extension should be .vtu!\n");
    return 2;
  }

  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkIdType cellNum = reader->GetNumberOfCells();
  vtkIdType pointNum = reader->GetNumberOfPoints();
  printf("There are %lld cells in the unstructed grid.\n", cellNum);
  printf("There are %lld points in the unstructed grid.\n", pointNum);

  // Partition the dataset 
  vtkUnstructuredGrid *usgrid = reader->GetOutput();
  usgrid->BuildLinks();
  vector<vector<vtkIdType>> regions;
  set<pair<vtkIdType, vtkIdType>> globalBridgeSet;
  decompose(threadNum, usgrid, regions, globalBridgeSet);
  
  // Test the domain decomposition and global bridge set
  // for(unsigned i = 0; i < regions[0].size(); i++){
  //   printf("[%lld], ", regions[0][i]);
  // }
  // printf("\n");
  // printf("Size of global bridge set: %zu\n", globalBridgeSet.size());
  // for(auto iter = globalBridgeSet.begin(); iter != globalBridgeSet.end(); iter++){
  //   printf("<%lld, %lld>\n", (*iter).first, (*iter).second);
  // }

  // Test argsort function
  // vector<vtkIdType> sortedIndices = MergeTree::argsort(regions[0], usgrid, false);
  // float *scalars = (float *)(usgrid->GetPointData()->GetArray(0)->GetVoidPointer(0));
  // for (unsigned int i = 0; i < sortedIndices.size(); i++) {
  //   printf("id: %lld, %.3f\n", sortedIndices[i], scalars[sortedIndices[i]]);
  // }


  
  // OpenMP routine
  omp_set_num_threads(threadNum);
  #pragma omp parallel
  {
    unsigned int tid = omp_get_thread_num();

    if(tid < regions.size()){
      // printf("Thread id = %d\n", tid);

      // TODO: Construct the local merge tree with the vertex set
      
      // TODO: Construct the reduced bridge set
      set<pair<vtkIdType, vtkIdType>> reducedBridgeSet = getReducedBridgeSet(globalBridgeSet, regions[tid], usgrid);
      
    }

    // if(tid == 0){
    //   nthreads = omp_get_num_threads();
    //   printf("Number of threads = %d\n", nthreads);
    // }
  }

  return 0;
}
