#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

// #include <filesystem>
// namespace fs = std::filesystem;

#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>

#include <vtkPLYWriter.h>
#include <vtkXMLPolyDataWriter.h>

int main ( int argc, char *argv[] )
{
    // Ensure a filename was specified
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " directory" << endl;
        return EXIT_FAILURE;
    }

    // // Get the filename from the command line
    // std::string inputFilename = argv[1];
    
    for(auto& p: fs::recursive_directory_iterator(argv[1]))
    {
        if (fs::path(p).extension() == ".vtk")
        {
            // Get all data from the file
            vtkSmartPointer<vtkGenericDataObjectReader> reader =
                vtkSmartPointer<vtkGenericDataObjectReader>::New();
            reader->SetFileName(p.path().c_str());
            reader->Update();

            std::cout << "VTK file name: " << p.path().c_str() << endl;
            std::cout << "Num of scalar fields: " 
                      << reader->GetNumberOfScalarsInFile() << endl;

            if(reader->IsFilePolyData())
            {
                fs::path outputFileName = p.path();
                outputFileName.replace_extension(".vtp");
    
                vtkSmartPointer<vtkXMLPolyDataWriter> writer =
                    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
                writer->SetFileName(outputFileName.c_str());
                writer->SetInputConnection(reader->GetOutputPort());
                writer->Update();
            }
            else
            {
                std::cerr << "VTK file is not poly data" << endl;
                return EXIT_FAILURE;
            }
        }
    }
    
    return EXIT_SUCCESS;
}
