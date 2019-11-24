#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkNamedColors.h>

int main ( int argc, char *argv[] )
{
  //parse command line arguments
  if(argc < 2)
  {
    std::cerr << "Usage: " << argv[0]
              << " Filename(.vtu)" << std::endl;
    return EXIT_FAILURE;
  }

  std::string filename = argv[1];

  //read all the data from the file
  vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
    vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkIdType cellNum = reader->GetNumberOfCells();
  vtkIdType pointNum = reader->GetNumberOfPoints();
  std::cout << "There are " << cellNum << " cells in the triangulation.\n";
  std::cout << "There are " << pointNum << " points in the triangulation.\n";

  return EXIT_SUCCESS;
}