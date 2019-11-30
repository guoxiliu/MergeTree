#include "Utils.h"
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
    return EXIT_FAILURE;
  }

  string filename = argv[1];
  string extension = filename.substr(filename.length() - 3);

  if(extension != "vtu"){
    fprintf(stderr, "The file extension should be .vtu!\n");
    return -1;
  }

  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkIdType cellNum = reader->GetNumberOfCells();
  vtkIdType pointNum = reader->GetNumberOfPoints();
  printf("There are %lld cells in the unstructed grid.\n", cellNum);
  printf("There are %lld points in the unstructed grid.\n", pointNum);

  // Partition the dataset 
  vector<vector<vtkIdType>> regions;
  set<vtkIdType> globalBridgeSet;
  decompose(threadNum, reader->GetOutput(), regions, globalBridgeSet);
  printf("Size of global bridge set: %lu\n", globalBridgeSet.size());

  // OpenMP test
  // omp_set_num_threads(threadNum);
  // #pragma omp parallel
  // {
  //   int nthreads, tid;
  //   tid = omp_get_thread_num();
  //   printf("Thread id = %d\n", tid);
  //   if(tid == 0){
  //     nthreads = omp_get_num_threads();
  //     printf("Number of threads = %d\n", nthreads);
  //   }
  // }

  return EXIT_SUCCESS;
}
