cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(vtk2vtp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(VTK
  COMPONENTS 
    vtkCommonCore
    vtkCommonDataModel
    vtkIOLegacy
    vtkIOPLY
    vtkIOXML)
if (NOT VTK_FOUND)
  message("Skipping vtk2vtp: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()

message (STATUS "VTK_VERSION: ${VTK_VERSION}")

add_executable(vtk2vtp MACOSX_BUNDLE vtk2vtp.cxx)

target_link_libraries(vtk2vtp PRIVATE ${VTK_LIBRARIES} stdc++fs)

vtk_module_autoinit(TARGETS vtk2vtp MODULES ${VTK_LIBRARIES})
