#include "MergeTree.h"
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataReader.h>

using namespace std;

int main ( int argc, char *argv[] )
{
  //parse command line arguments
  if(argc < 2){
    cerr << "Usage: " << argv[0] << " Filename(.vti)" << endl;
    return EXIT_FAILURE;
  }


  string filename = argv[1];
  string extension = filename.substr(filename.length() - 3);

  if(extension != "vti"){
    cout << "The file extension should be .vti!\n";
    return -1;
  }

  vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkIdType cellNum = reader->GetNumberOfCells();
  vtkIdType pointNum = reader->GetNumberOfPoints();
  cout << "There are " << cellNum << " cells in the triangulation.\n";
  cout << "There are " << pointNum << " points in the triangulation.\n";

  // Create the merge tree here.
  MergeTree testTree(reader->GetOutput());
  auto start = chrono::high_resolution_clock::now();
  testTree.build();
  auto stop = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << "Build merge Tree cost: " << duration.count() << " microseconds" <<endl;
  // Test the queries here.
  set<pair<vtkIdType, vtkIdType>> emptyBridgeSet;
  start = chrono::high_resolution_clock::now();
  vector<vtkIdType> maxima = testTree.MaximaQuery(emptyBridgeSet);
  stop = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << "MaximaQuery cost: " << duration.count() << " microseconds" <<endl;

  printf("The size of the maxima is %zu\n", maxima.size());
  /* for (unsigned int i = 0; i < maxima.size(); i++) {
    printf("maxima[%u]: %lld\n", i, maxima[i]);
  } */
  start = chrono::high_resolution_clock::now();
  vtkIdType v = 0;
  float *scalar = (float*) getScalar(reader->GetOutput());
  float level = scalar[v];
  vtkIdType  CompMaxima = testTree.ComponentMaximumQuery(v,level);
  stop = chrono::high_resolution_clock::now();
  duration = chrono::duration_cast<chrono::microseconds>(stop - start);
  cout << "ComponentMaximaQuery cost: " << duration.count() << " microseconds" <<endl;

 printf("the component maxima is %d\n", (int)CompMaxima);	
  return EXIT_SUCCESS;
}
