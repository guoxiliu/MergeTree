#include "MergeTree.h"
#include <omp.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <chrono>

using namespace std;

const int threadNum = 16;
const bool runMode = 0;

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

  if(runMode){
    // Test merge tree
    MergeTree testTree(sgrid);
    auto start = chrono::high_resolution_clock::now();
    testTree.build();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    printf("Build merge tree cost: %lld\n", duration.count());
  }else{

    vector<vector<vtkIdType>> regions;
    set<pair<vtkIdType, vtkIdType>> globalBridgeSet;
    auto start = chrono::high_resolution_clock::now();
    decompose(threadNum, sgrid, regions, globalBridgeSet);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    printf("Decomposition cost: %lld\n", duration.count());
    
    // Test the domain decomposition and global bridge set
    // for(size_t i = 0; i < regions.size(); i++){
    //   printf("region %lu: <%lld, %lld> \n", i, regions[i].front(), regions[i].back());
    // }
    // printf("\n");
    printf("Size of global bridge set: %zu\n", globalBridgeSet.size());
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
    //set<pair<vtkIdType, vtkIdType>> reducedGlobalBS = getReducedBridgeSet(globalBridgeSet, allVertices, sgrid);
    // Test the reduced global bridge set
    // printf("Size of reduced global bridge set: %zu\n", reducedGlobalBS.size());
    // for(auto iter = globalBridgeSet.begin(); iter != globalBridgeSet.end(); iter++){
    //   printf("<%lld, %lld>\n", (*iter).first, (*iter).second);
    // }
    

    // OpenMP routine
    start = chrono::high_resolution_clock::now();
    vector<vtkIdType> maxima;   // use for maxima query
    omp_set_num_threads(threadNum);
    #pragma omp parallel
    {
      unsigned int tid = omp_get_thread_num();

      if(tid < regions.size()){
        // printf("Thread id = %d\n", tid);

        // Construct the local merge tree with the vertex set
        MergeTree localMergeTree(sgrid, regions[tid]);
        localMergeTree.build();
        
        // Construct the reduced bridge set

        // Perform queries
        //vector<vtkIdType> regionMaxima = localMergeTree.MaximaQuery(localBridgeSet);
        //maxima.insert(maxima.end(), regionMaxima.begin(), regionMaxima.end());
      }

      // if(tid == 0){
      //   nthreads = omp_get_num_threads();
      //   printf("Number of threads = %d\n", nthreads);
      // }
    }

    stop = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    printf("Build tree cost: %lld\n", duration.count());

    // printf("The size of the maxima is %zu\n", maxima.size());
    // for (unsigned int i = 0; i < maxima.size(); i++) {
    //   printf("maxima[%u]: %lld\n", i, maxima[i]);
    // }

  }

  return 0;
}
