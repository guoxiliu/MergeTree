#include "MergeTree.h"
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLUnstructuredGridReader.h>

using namespace std;

int main ( int argc, char *argv[] )
{
  //parse command line arguments
  if(argc < 2){
    cerr << "Usage: " << argv[0] << " Filename(.vtu)" << endl;
    return EXIT_FAILURE;
  }


  string filename = argv[1];
  string extension = filename.substr(filename.length() - 3);

  if(extension != "vtu"){
    cout << "The file extension should be .vtu!\n";
    return -1;
  }

  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkIdType cellNum = reader->GetNumberOfCells();
  vtkIdType pointNum = reader->GetNumberOfPoints();
  cout << "There are " << cellNum << " cells in the triangulation.\n";
  cout << "There are " << pointNum << " points in the triangulation.\n";

  // TODO: Create the merge tree here.
  MergeTree testTree(reader->GetOutput());
  testTree.build();

  // TODO: Test the queries here.


  return EXIT_SUCCESS;
}
