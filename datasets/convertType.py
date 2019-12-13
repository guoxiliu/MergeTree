
import os
import vtk
import numpy as np
from vtk.util import numpy_support

# The function is modified based on https://code.ornl.gov/rwp/javelin/blob/master/javelin/io.py.
def numpy_to_vti(array, origin, spacing, filename):
    """
    This function writes a vti file from a numpy array.

    Parameters
    ---------------

    array: `ndarray`
        input numpy array
    origin: `array like object`
        the origin of the array
    spacing: `array like object`
        the step in eac  h dimension
    filename: str
        output filename (.vti)

    Returns
    ---------------
    None
    """

    if array.ndim != 3:
        raise ValueError("Only works with 3 dimensional arrays")

    vtkArray = numpy_support.numpy_to_vtk(num_array=array.ravel(), deep=True,
                            array_type=vtk.VTK_FLOAT)

    imageData = vtk.vtkImageData()
    imageData.SetOrigin(origin)
    imageData.SetSpacing(spacing)
    imageData.SetDimensions(array.shape)
    imageData.GetPointData().SetScalars(vtkArray)

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(imageData)
    writer.Write()


# Read the raw file 
dir_name = os.path.dirname(os.path.abspath(__file__))
base_filename = "silicium_98x34x34_uint8"
filename_suffix = ".raw"

inputfile = os.path.join(dir_name, base_filename + filename_suffix)
print("Reading from " + inputfile)

# Set arguments for the function
scalars = np.fromfile(inputfile, dtype='uint8')
scalars = np.reshape(scalars, (98, 34, 34))

origin = (0, 0, 0)
spacing = (1, 1, 1)

head_filename = base_filename.split('_', 1)[0]
outfile = os.path.join(dir_name, head_filename + '.vti')

numpy_to_vti(scalars, origin, spacing, outfile)

# toydata = np.random.randint(20, size=(20, 20, 20), dtype='uint8')
# print(toydata.ravel())
# outfilename = "toy.vti"
# numpy_to_vti(toydata, origin, spacing, outfilename)

print("Finished!")