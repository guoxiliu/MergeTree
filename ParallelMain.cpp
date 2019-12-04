#include "MergeTree.h"
#include <omp.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>

using namespace std;

const int threadNum = 8;

int main ( int argc, char *argv[] )
{
  // parse command line arguments
  if(argc < 2){
    fprintf(stderr, "Usage: %s Filename(.vti)\n", argv[0]);
    return 1;
  }

  string filename = argv[1];
  string extension = filename.substr(filename.length() - 3);

  if(extension != "vti"){
    fprintf(stderr, "The file extension should be .vti!\n");
    return 2;
  }

  vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkIdType cellNum = reader->GetNumberOfCells();
  vtkIdType pointNum = reader->GetNumberOfPoints();
  printf("There are %lld cells in the data.\n", cellNum);
  printf("There are %lld points in the data.\n", pointNum);

  // Partition the dataset 
  vtkImageData *sgrid = reader->GetOutput();
  vector<vector<vtkIdType>> regions;
  set<pair<vtkIdType, vtkIdType>> globalBridgeSet;

  // Test merge tree
  MergeTree testTree(sgrid);
  testTree.build();
  // decompose(threadNum, sgrid, regions, globalBridgeSet);
  
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
  // vector<vtkIdType> sortedIndices = argsort(regions[0], sgrid, false);
  // float *scalars = (float *)(sgrid->GetPointData()->GetArray(0)->GetVoidPointer(0));
  // for (unsigned int i = 0; i < sortedIndices.size(); i++) {
  //   printf("id: %lld, %.3f\n", sortedIndices[i], scalars[sortedIndices[i]]);
  // }
  
  // vector<vtkIdType> allVertices(sgrid->GetNumberOfPoints());
  // iota(allVertices.begin(), allVertices.end(), 0);
  // set<pair<vtkIdType, vtkIdType>> reducedGlobalBS = getReducedBridgeSet(globalBridgeSet, allVertices, sgrid);
  // Test the reduced global bridge set
  // printf("Size of reduced global bridge set: %zu\n", reducedGlobalBS.size());
  // for(auto iter = globalBridgeSet.begin(); iter != globalBridgeSet.end(); iter++){
  //   printf("<%lld, %lld>\n", (*iter).first, (*iter).second);
  // }
  


  
  // // OpenMP routine
  // vector<vtkIdType> maxima;   // use for maxima query
  // omp_set_num_threads(threadNum);
  // #pragma omp parallel
  // {
  //   unsigned int tid = omp_get_thread_num();

  //   if(tid < regions.size()){
  //     // printf("Thread id = %d\n", tid);

  //     // Construct the local merge tree with the vertex set
  //     MergeTree localMergeTree(sgrid);
  //     localMergeTree.build(regions[tid]);
      
  //     // Construct the reduced bridge set
  //     //set<pair<vtkIdType, vtkIdType>> localBridgeSet = getLocalBridgeSet(reducedGlobalBS, regions[tid]);

  //     // Perform queries
  //     //vector<vtkIdType> regionMaxima = localMergeTree.MaximaQuery(localBridgeSet);
  //     //maxima.insert(maxima.end(), regionMaxima.begin(), regionMaxima.end());
  //   }

  //   // if(tid == 0){
  //   //   nthreads = omp_get_num_threads();
  //   //   printf("Number of threads = %d\n", nthreads);
  //   // }
  // }

  // printf("The size of the maxima is %zu\n", maxima.size());
  // for (unsigned int i = 0; i < maxima.size(); i++) {
  //   printf("maxima[%u]: %lld\n", i, maxima[i]);
  // }

  return 0;
}
